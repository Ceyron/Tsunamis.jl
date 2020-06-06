#=
   Contains a single-node sequential implementation of an SWE solver for the simulation of tsunamis
=#

#import Plots
using Printf
using Distributed
using ArgParse

include("src/CartesianBlock.jl")
include("src/single_node_setup.jl")
include("src/writer.jl")
include("src/godunov.jl")
include("scenario/radial_dam_break.jl")
include("src/boundary.jl")


# Some stub constants for the simulation, has to be set by command-line options
# later on
const offset_x = 0.;
const offset_y = 0.;
const time_end = 15.;
const cfl_number = 0.4;

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
            default = "test"
        "-c"
            help = "number of checkpoints"
            arg_type = Int
            default = 20
        "--no-io"
            help = "Do not write an output file"
            action = :store_true
    end

    parsed_args = parse_args(arg_parse_settings)
    num_cells_x = parsed_args["x"]
    num_cells_y = parsed_args["y"]
    output_name = parsed_args["o"]
    num_checkpoints = parsed_args["c"]
    no_io = parsed_args["no-io"]

    println()
    println("Welcome to the SWE solver using Julia")
    println("-------------------------------------")

    # Get the size of the simulation domain
    domain_size_x, domain_size_y = radial_dam_break_get_domain_size()

    # Instantiate main block
    simulation_single_node = SWE_Simulation(
        Simulation_Block(
            Layout(
                # No offset since it is single block in a sequential implementation
                offset_x,
                offset_y,
                # Set the size of the simulation domain (halo cells do not count
                # towards this)
                domain_size_x,
                domain_size_y,
                # Set number of interior cells
                num_cells_x,
                num_cells_y,
                # Calculate and set the cell widths
                domain_size_x / num_cells_x,
                domain_size_y /num_cells_y,
                # Create the mesh, for N interior cells we have N+1 edges in each direction
                range(offset_x, offset_x+domain_size_x; length=num_cells_x+1),
                range(offset_y, offset_y+domain_size_y; length=num_cells_y+1),
            ),
            SWE_Fields(
                # Account for the additional halo layers of the boundary
                zeros(num_cells_x + 2, num_cells_y + 2),
                zeros(num_cells_x + 2, num_cells_y + 2),
                zeros(num_cells_x + 2, num_cells_y + 2),
            ),
            # The Bathymetry data
            zeros(num_cells_x + 2, num_cells_y + 2),
            # No particular distribution occurs in the single_node case,
            # therefore set all boundary conditions to periodic
            Boundary_Collection(
                Boundary(
                    radial_dam_break_get_boundary_type(),
                    RemoteChannel(),
                    RemoteChannel(),
                ),
                Boundary(
                    radial_dam_break_get_boundary_type(),
                    RemoteChannel(),
                    RemoteChannel(),
                ),
                Boundary(
                    radial_dam_break_get_boundary_type(),
                    RemoteChannel(),
                    RemoteChannel(),
                ),
                Boundary(
                    radial_dam_break_get_boundary_type(),
                    RemoteChannel(),
                    RemoteChannel(),
                ),
            ),
        ),
        # Contains information on the time integration
        Time_Mesh(
            time_end,
            num_checkpoints,
            range(0, time_end; length=num_checkpoints),
        ),
        # Set initial time to zero
        0.0
    )

    # Instantiate the container for all the selected simulation settings
    simulation_settings = SWE_Simulation_Settings(
        output_name,
        cfl_number,
    )

    # Imprint the bathymetry
    radial_dam_break_imprint_bathymetry!(simulation_single_node)

    # Imprint the initial condition
    radial_dam_break_imprint_initial_condition!(simulation_single_node)

    # Create the netCDF file
    if !no_io
        nc_data_set = create_output_file(
            simulation_settings.output_file_name,
            simulation_single_node
        )
    end

    # Write the initial state to the cdf file
    if !no_io
        write_fields!(nc_data_set, simulation_single_node, 1)
    end
    # Instantiate the flux field struct to be used over the iterations, they
    # contain the accumulated (~= cummulative) fluxes summed up from the fluxes
    # over each edge For loop convenience in the update routines we also
    # calculate flux summations for the halo cells (bottom and left) even though they are not
    # updated
    fluxes = SWE_Fields(
        Array{Float64, 2}(undef, num_cells_x + 2, num_cells_y + 2),
        Array{Float64, 2}(undef, num_cells_x + 2, num_cells_y + 2),
        Array{Float64, 2}(undef, num_cells_x + 2, num_cells_y + 2),
    )

    # Iterate over all checkpoints
    @time for i_checkpoint in 2:num_checkpoints

        # Integrate until next checkpoint is reached
        while simulation_single_node.time <
                simulation_single_node.time_mesh.time_nodes[i_checkpoint]
            println("Simulating at time $(simulation_single_node.time)",
                "/$(simulation_single_node.time_mesh.time_end)")            

            # Clear the flux fields
            fluxes.h .= 0.0
            fluxes.hu .= 0.0
            fluxes.hv .= 0.0


            # The maximum wave speed is relevant for the CFL condition
            max_wave_speed = 0.0

            # (1) set values in ghost layer
            update_boundaries!(simulation_single_node)

            # (2) Compute numerical fluxes
            max_wave_speed = calculate_numerical_fluxes!(
                fluxes,
                simulation_single_node
            )

            # (3) Calculate the new time step (maximum allows time_step due to
            # the wave speeds)
            time_step = compute_max_time_step(
                simulation_single_node.current.layout.cell_width_x,
                simulation_single_node.current.layout.cell_width_y,
                max_wave_speed,
                simulation_settings.clf_number,
            )

            # (4) Update the cell values Godunov style first order
            update_cells!(simulation_single_node, fluxes, time_step)

            # (5) Refresh to current simulation time
            simulation_single_node.time += time_step
        end

        # Save the fields
        if !no_io
            println("-> Saving Fields")
            write_fields!(nc_data_set, simulation_single_node, i_checkpoint)
        end
    end
    # Close the connection to the dataset handle
    if !no_io
        close_output_file(nc_data_set)
    end
end


# Start the main function
@time main()