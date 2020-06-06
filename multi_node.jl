# Contains a multi-node implementation of the Shallow Water Equation Solver. It
# performs a domain-decomposition with the help of the built-in Distributed
# package of julia using a task-based approach to multiprocessing.
# Start it with julia -p NUM_PROCESSORS or with the help of the machine list

using Distributed
using ArgParse
@everywhere using Printf

# The cd are necessary since it seems julia has problems with finding relative
# paths on remote workers
@everywhere cd("src/")
@everywhere include("CartesianBlock.jl")
@everywhere include("multi_node_setup.jl")
@everywhere include("single_node_setup.jl")
@everywhere include("writer.jl")
@everywhere include("godunov.jl")
@everywhere include("boundary.jl")

@everywhere cd("../scenario/")
@everywhere include("radial_dam_break.jl")

@everywhere cd("../")


# Some stub constants for the simulation, has to be set by command-line options
# later on
const offset_x = 0.;
const offset_y = 0.;
const time_end = 15.;
const cfl_number = 0.4;

# Given the the number of processors available, try to distribute the the number
# of blocks in x and y direction
# Best results achieved for quadratic numbers
function compute_number_of_block_rows(number_of_processors::Int)
    number_of_rows::Int = Int(floor(sqrt(number_of_processors)))
    while number_of_processors % number_of_rows != 0
        number_of_rows -= 1
    end
    return number_of_rows
end

# Transforms the 2d index to a block in the domain to its corresponding CPU
@everywhere @inline function processor_2d_to_id(i, j, number_of_blocks_x)
    # Addition of plus one necessary since CPU with id 1 is the master
    return (j-1) * number_of_blocks_x + i + 1
end

function main()
    # Parsing the arguments
    arg_parse_settings = ArgParseSettings()
    @add_arg_table! arg_parse_settings begin
        "-x"
            help = "number of cells in x"
            arg_type = Int
            default = 400
        "-y"
            help = "number of cells in y"
            arg_type = Int
            default = 400
        "-o"
            help = "output_file_name"
            arg_type = String
            default = "multi_test"
        "-c"
            help = "number of checkpoints"
            arg_type = Int
            default = 20
        "--no-io"
            help = "Do not output to netcdf"
            action = :store_true
            default = false
    end

    parsed_args = parse_args(arg_parse_settings)
    num_cells_x = parsed_args["x"]
    num_cells_y = parsed_args["y"]
    output_name = parsed_args["o"]
    num_checkpoints = parsed_args["c"]
    no_io::Bool = parsed_args["no-io"]


    println()
    println("Welcome to the distributed SWE equation solver")
    println("This solver is distributed using Julia's builtin primitives")
    println("----------------------------------------------------------")
    if nprocs() == 1
        println("ERROR: You have to start this program in Julia'a multiprocessing mode")
        println("ERROR: do julia -p [NUM_CORES] multi_node.jl [ARGS]")
        exit(1)
    end
    number_of_processors = nworkers()
    number_of_blocks_y::Int = compute_number_of_block_rows(number_of_processors)
    number_of_blocks_x::Int = number_of_processors / number_of_blocks_y
    # This must not necessarily equal the number of CPUs
    number_of_blocks_total = number_of_blocks_x * number_of_blocks_y
    println("--> $(number_of_processors) processors provided")
    println("--> Using $number_of_blocks_x blocks in x direction")
    println("--> Using $number_of_blocks_y blocks in y direction")
    println("--> Using $(number_of_blocks_total) blocks (=CPUs) in total")
    println()

    # Instantiate main block
    # Instantiate decomposed domain on each process of the distributed task
    # Then have a function that transfers the halo cells from one block to another
    simulation_multi_node = SWE_Simulation_Distributed(
        Array{
            Future, 2
        }(undef, number_of_blocks_x, number_of_blocks_y),
        # Contains information on the time integration
        Time_Mesh(
            time_end,
            num_checkpoints,
            range(0, time_end; length=num_checkpoints),
        ),
        # Set the initial time to zero
        0.0,
        Block_Mesh(
            range(offset_x, offset_x+domain_size_x; length=number_of_blocks_x+1),
            range(offset_y, offset_y+domain_size_y; length=number_of_blocks_y+1),
            number_of_blocks_x,
            number_of_blocks_y,
            number_of_blocks_total,
        ),
        # Create an empty array for the handles to the nc datasets, each block
        # will have its own netcdf file
        Array{Future, 2}(undef, number_of_blocks_x, number_of_blocks_y)
    )

    # Calculate the cell widths. This are valid for all blocks since we are
    # using an equidistant cartesian grid
    cell_width_x = domain_size_x / num_cells_x
    cell_width_y = domain_size_y /num_cells_y

    # Create the array of channels. These channels are buffered with one element
    # and send the domain boundary layer exchange. We need two per edge, since
    # data has to be sent in both directions
    #
    # Since data is only exchanged at the interior edges, the size of this array
    # equals the number of blocks in the given direction (direction parallel to
    # the edge normal) minus 1 times the number of blocks in the direction
    # perpendicular to the edge normal
    channels_left = Array{RemoteChannel, 2}(
        undef,
        number_of_blocks_x,
        number_of_blocks_y,
    )
    channels_top = Array{RemoteChannel, 2}(
        undef,
        number_of_blocks_x,
        number_of_blocks_y,
    )
    channels_right = Array{RemoteChannel, 2}(
        undef,
        number_of_blocks_x,
        number_of_blocks_y,
    )
    channels_bottom = Array{RemoteChannel, 2}(
        undef,
        number_of_blocks_x,
        number_of_blocks_y,
    )
    # Instantiate the RemoteChannels with the underlying Channel structure being
    # located on the respective core
    for i in 1:number_of_blocks_x
        for j in 1:number_of_blocks_y
            # Instantiate each channel with a buffer of one
            channels_left[i, j] = RemoteChannel(() -> Channel{SWE_Copy_Fields}(1), processor_2d_to_id(i, j, number_of_blocks_x))
            channels_top[i, j] = RemoteChannel(() -> Channel{SWE_Copy_Fields}(1), processor_2d_to_id(i, j, number_of_blocks_x))
            channels_right[i, j] = RemoteChannel(() -> Channel{SWE_Copy_Fields}(1), processor_2d_to_id(i, j, number_of_blocks_x))
            channels_bottom[i, j] = RemoteChannel(() -> Channel{SWE_Copy_Fields}(1), processor_2d_to_id(i, j, number_of_blocks_x))
        end
    end


    # Determine how many cells we use per block. We approach it that we try to
    # assign the biggest number without remainder to the n-1 blocks in a
    # row/column and then assign the rest of the cells to the leftover block

    # Create all blocks on the remote processes and add them to the list of
    # references to these remote datastructures
    for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
        for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
            # The number of cells for the current considered block
            number_of_cells_this_block_x::Int = 0
            number_of_cells_this_block_y::Int = 0
            # Calculate the local number of cells in x for this block
            if i < simulation_multi_node.block_mesh.number_of_blocks_x
                number_of_cells_this_block_x = floor(num_cells_x /
                    simulation_multi_node.block_mesh.number_of_blocks_x)
            else
                number_of_cells_this_block_x = num_cells_x -
                    (simulation_multi_node.block_mesh.number_of_blocks_x - 1) *
                    (floor(num_cells_x / simulation_multi_node.block_mesh.number_of_blocks_x))
            end
            # Calculate the local number of cells in y for this block
            if j < simulation_multi_node.block_mesh.number_of_blocks_y
                number_of_cells_this_block_y = floor(num_cells_y /
                    simulation_multi_node.block_mesh.number_of_blocks_y)
            else
                number_of_cells_this_block_y = num_cells_y -
                    (simulation_multi_node.block_mesh.number_of_blocks_y - 1) *
                    (floor(num_cells_y / simulation_multi_node.block_mesh.number_of_blocks_y))
            end

            println("CPU $(processor_2d_to_id(i, j, number_of_blocks_x) - 1) ",
                "spans x=[$(simulation_multi_node.block_mesh.block_mesh_x[i]),",
                " $(simulation_multi_node.block_mesh.block_mesh_x[i+1])] ",
                "and y=[$(simulation_multi_node.block_mesh.block_mesh_y[j]), ",
                "$(simulation_multi_node.block_mesh.block_mesh_y[j+1])] ")
            println(
                "with $(number_of_cells_this_block_x) cells in x ",
                "and $(number_of_cells_this_block_y) in y"
            )

            simulation_multi_node.block_references[i, j] =
                @spawnat processor_2d_to_id(i, j, number_of_blocks_x) SWE_Simulation(
                    Simulation_Block(
                        Layout(
                            # Calculate the offset by the position of the block_mesh
                            simulation_multi_node.block_mesh.block_mesh_x[i],
                            simulation_multi_node.block_mesh.block_mesh_y[j],
                            # Set the size of this simulation block (not including the
                            # halo cells)
                            simulation_multi_node.block_mesh.block_mesh_x[i+1] -
                                simulation_multi_node.block_mesh.block_mesh_x[i],
                            simulation_multi_node.block_mesh.block_mesh_y[j+1] -
                                simulation_multi_node.block_mesh.block_mesh_y[j],
                            # Set the number of cells for this block
                            number_of_cells_this_block_x,
                            number_of_cells_this_block_y,
                            # Set the cell widths
                            cell_width_x,
                            cell_width_y,
                            # Create the mesh for within the block
                            range(
                                simulation_multi_node.block_mesh.block_mesh_x[i],
                                simulation_multi_node.block_mesh.block_mesh_x[i+1];
                                length=number_of_cells_this_block_x+1
                            ),
                            range(
                                simulation_multi_node.block_mesh.block_mesh_y[j],
                                simulation_multi_node.block_mesh.block_mesh_y[j+1];
                                length=number_of_cells_this_block_y+1
                            ),
                        ),
                        # These arrays save the field quantities (h, hu, hv). We have to
                        # add one cell in each direction. This will be used as a halo
                        # cell either for imprinting the Boundary Condition or for
                        # exchanging with the neighbor domain
                        SWE_Fields(
                            zeros(number_of_cells_this_block_x + 2, number_of_cells_this_block_y + 2),
                            zeros(number_of_cells_this_block_x + 2, number_of_cells_this_block_y + 2),
                            zeros(number_of_cells_this_block_x + 2, number_of_cells_this_block_y + 2),
                        ),
                        # The bathymetry data
                        zeros(number_of_cells_this_block_x + 2, number_of_cells_this_block_y + 2),
                        # Assign the boundaries. The outermost blocks obviously have
                        # true BC at their outer edges. All other edges have the
                        # connected BC with the block next to it
                        #
                        # The connection is realized by RemoteChannels. Convention: The
                        # block is sending on its channel and receiving on the other's
                        Boundary_Collection(
                            Boundary(
                                i > 1 ? CONNECT : radial_dam_break_get_boundary_type(), 
                                i > 1 ? channels_left[i-1, j] : RemoteChannel(),
                                i > 1 ? channels_right[i, j] : RemoteChannel(),
                            ),
                            Boundary(
                                j < number_of_blocks_y ? CONNECT : radial_dam_break_get_boundary_type(),
                                j < number_of_blocks_y ? channels_top[i, j+1] : RemoteChannel(),
                                j < number_of_blocks_y ? channels_bottom[i, j] : RemoteChannel(),
                            ),
                            Boundary(
                                i < number_of_blocks_x ? CONNECT : radial_dam_break_get_boundary_type(),
                                i < number_of_blocks_x ? channels_right[i+1, j] : RemoteChannel(),
                                i < number_of_blocks_x ? channels_left[i, j] : RemoteChannel(),
                            ),
                            Boundary(
                                j > 1 ? CONNECT : radial_dam_break_get_boundary_type(),
                                j > 1 ? channels_bottom[i, j-1] : RemoteChannel(),
                                j > 1 ? channels_top[i, j] : RemoteChannel(),
                            ),
                        ),
                    ),
                    # Contains information on the time integration
                    Time_Mesh(
                        time_end,
                        num_checkpoints,
                        range(0, time_end; length=num_checkpoints),
                    ),
                    # Set initial time to zero, but this quantity won't be mutated.
                    # Just for compatibility with the single-node implementation
                    0.0
                )
        end
    end
    println()

    # Instantiate the container for all the selected simulation settings
    simulation_settings = SWE_Simulation_Settings(
        output_name,
        cfl_number,
    )

    # Imprint the bathymetry on all blocks
    for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
        for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
            @spawnat processor_2d_to_id(i, j, number_of_blocks_x) radial_dam_break_imprint_bathymetry!(
                fetch(simulation_multi_node.block_references[i, j])
            )
        end
    end

    # Imprint the initial condition on all blocks
    for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
        for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
            @spawnat processor_2d_to_id(i, j, number_of_blocks_x) radial_dam_break_imprint_initial_condition!(
                fetch(simulation_multi_node.block_references[i, j])
            )
        end
    end

    # Create the netCDF file for each block
    if !no_io
        for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
            for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                simulation_multi_node.nc_data_sets[i, j] =
                    @spawnat processor_2d_to_id(i, j, number_of_blocks_x) create_output_file(
                    "$(simulation_settings.output_file_name)_$(processor_2d_to_id(i, j, number_of_blocks_x))",
                    fetch(simulation_multi_node.block_references[i, j])
                )
            end
        end
    end

    # Write the initial state to the cdf file for each block
    if !no_io
        @sync for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
            for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                @spawnat processor_2d_to_id(i, j, number_of_blocks_x) write_fields!(
                    fetch(simulation_multi_node.nc_data_sets[i, j]),
                    fetch(simulation_multi_node.block_references[i, j]),
                    1
                )
            end
        end
    end

    # The futures to the computed wave speeds of every block
    max_wave_speed_references = Array{Future, 2}(
        undef,
        simulation_multi_node.block_mesh.number_of_blocks_x,
        simulation_multi_node.block_mesh.number_of_blocks_y,
    )

    # Iterate over all checkpoints
    @time for i_checkpoint in 2:num_checkpoints

        # Integrate until next checkpoint is reached
        while simulation_multi_node.time <
                simulation_multi_node.time_mesh.time_nodes[i_checkpoint]
            println("Simulating at time $(simulation_multi_node.time)",
                "/$(simulation_multi_node.time_mesh.time_end)")            

            # Clear the flux fields on each node
            fluxes_references = Array{Future, 2}(
                undef,
                simulation_multi_node.block_mesh.number_of_blocks_x,
                simulation_multi_node.block_mesh.number_of_blocks_y,
            )
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    fluxes_references[i, j] = @spawnat processor_2d_to_id(i, j, number_of_blocks_x) SWE_Fields(
                        zeros(
                            fetch(simulation_multi_node.block_references[i, j]).current.layout.num_interior_cells_x + 2,
                            fetch(simulation_multi_node.block_references[i, j]).current.layout.num_interior_cells_y + 2,
                        ),
                        zeros(
                            fetch(simulation_multi_node.block_references[i, j]).current.layout.num_interior_cells_x + 2,
                            fetch(simulation_multi_node.block_references[i, j]).current.layout.num_interior_cells_y + 2,
                        ),
                        zeros(
                            fetch(simulation_multi_node.block_references[i, j]).current.layout.num_interior_cells_x + 2,
                            fetch(simulation_multi_node.block_references[i, j]).current.layout.num_interior_cells_y + 2,
                        ),
                    )
                end
            end


            # The maximum wave speed is relevant for the CFL condition
            max_wave_speed = 0.0

            # (1) set values in ghost layer (either by Boundary Condition of
            # copied layer from neighboring domain)
            # Has to be synced otherwise a
            # fast processor would look when taking from the channel
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    @spawnat processor_2d_to_id(i, j, number_of_blocks_x) queue_copy_layer(
                        fetch(simulation_multi_node.block_references[i, j])
                    )
                end
            end
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    @spawnat processor_2d_to_id(i, j, number_of_blocks_x) update_boundaries!(
                        fetch(simulation_multi_node.block_references[i, j])
                    )
                end
            end

            # (2) Compute numerical fluxes In our distributed case we hold the
            # future to the wave_speed calculation and can then have a barrier
            # to wait for each flux compuation to finish and find the overall
            # max wave speed Spawn the task to compute the numerical fluxes on
            # each block by its corresponding processor
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    max_wave_speed_references[i, j] = @spawnat processor_2d_to_id(i, j, number_of_blocks_x) calculate_numerical_fluxes!(
                        fetch(fluxes_references[i, j]),
                        fetch(simulation_multi_node.block_references[i, j]),
                    )
                end
            end
            # Collect the maximum wave speed over all blocks
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    max_wave_speed = max(
                        max_wave_speed,
                        fetch(max_wave_speed_references[i, j])
                    )
                end
            end

            # (3) Calculate the new time step (maximum allows time_step due to
            # the wave speeds)
            time_step = compute_max_time_step(
                cell_width_x,
                cell_width_y,
                max_wave_speed,
                simulation_settings.clf_number,
            )

            # (4) Update the cell values Godunov style first order
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    @spawnat processor_2d_to_id(i, j, number_of_blocks_x) update_cells!(
                        fetch(simulation_multi_node.block_references[i, j]),
                        fetch(fluxes_references[i, j]),
                        time_step,
                    )
                end
            end

            # (5) Refresh to current simulation time
            simulation_multi_node.time += time_step
        end

        # Save the fields
        if !no_io
            println("-> Saving Fields")
            for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
                for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                    @spawnat processor_2d_to_id(i, j, number_of_blocks_x) write_fields!(
                        fetch(simulation_multi_node.nc_data_sets[i, j]),
                        fetch(simulation_multi_node.block_references[i, j]),
                        i_checkpoint,
                    )
                end
            end
        end
    end

    # Close the connection to the dataset handle
    if !no_io
        for i in 1:simulation_multi_node.block_mesh.number_of_blocks_x
            for j in 1:simulation_multi_node.block_mesh.number_of_blocks_y
                @spawnat processor_2d_to_id(i, j, number_of_blocks_x) close_output_file(
                    fetch(simulation_multi_node.nc_data_sets[i, j])
                )
            end
        end
    end
end

# Start the main function
@time main()