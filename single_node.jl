#=
   Contains a single-node sequential implementation of an SWE solver for the simulation of tsunamis
=#

#import Plots
using Printf

# Use Plotly backend in an electron chromium window (nice visuals :D)
#Plots.plotlyjs()

include("src/CartesianBlock.jl")
include("src/single_node_setup.jl")
include("src/writer.jl")


# Some stub constants for the simulation, has to be set by command-line options
# later on
const offset_x = 0.;
const offset_y = 0.;
const size_x = 5000.;
const size_y = 5000.;
const num_cells_x = 50;
const num_cells_y = 50;
const time_end = 15.;
const num_checkpoints = 20;
const output_name = "single_run"

function main()
    println()
    println("Welcome to the SWE solver using Julia")
    println()

    # Instantiate main block
    simulation_single_node = SWE_Simulation(
        Simulation_Block(
            Layout(
                # No offset since it is single block in a sequential implementation
                offset_x,
                offset_y,
                # Set the size of the simulation domain (halo cells do not count
                # towards this)
                size_x,
                size_y,
                # Set number of interior cells
                num_cells_x,
                num_cells_y,
                # Calculate and set the cell widths
                size_x / num_cells_x,
                size_y /num_cells_y,
                # Create the mesh, for N interior cells we have N+1 edges in each direction
                range(offset_x, offset_x+size_x; length=num_cells_x),
                range(offset_y, offset_y+size_y; length=num_cells_y),
            ),
            SWE_Fields(
                # Account for the additional halo layers of the boundary
                zeros(num_cells_x + 2, num_cells_y + 2),
                zeros(num_cells_x + 2, num_cells_y + 2),
                zeros(num_cells_x + 2, num_cells_y + 2),
            ),
            # Set initial time to zero
            0.0,
            # The Bathymetry data
            zeros(num_cells_x + 2, num_cells_y + 2),
            # No particular distribution occurs in the single_node case,
            # therefore set all boundary conditions to periodic
            Boundary_Collection(
                Boundary(
                    periodic,
                    "none"
                ),
                Boundary(
                    periodic,
                    "none"
                ),
                Boundary(
                    periodic,
                    "none"
                ),
                Boundary(
                    periodic,
                    "none"
                ),
            ),
        ),
        # Contains information on the time integration
        Time_Mesh(
            time_end,
            num_checkpoints,
            range(0, time_end; length=num_checkpoints),
        ),
    )

    # Instantiate the container for all the selected simulation settings
    simulation_settings = SWE_Simulation_Settings(
        output_name,
    )

    # Create the netCDF file
    nc_data_set = create_output_file(simulation_settings, simulation_single_node)

    # Imprint the initial condition
    simulation_single_node.current.fields.h[22:28, 22:28] .= 3.0

    # Write the initial state to the cdf file
    write_fields!(nc_data_set, simulation_single_node, 1)

    # Print the head
    println("Welcome to SWE solver on julia")
    println("time|todo")

    # Instantiate the flux field struct to be used over the iterations
    fluxes = Fluxes(
        Flux_Field(
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
        ),
        Flux_Field(
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
        ),
        Flux_Field(
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
            zeros(num_cells_x, num_cells_y),
        ),
        # The max wave speed standard value
        0.0,
    )

    # Iterate over all checkpoints
    for i_checkpoint in 2:num_checkpoints


        # Integrate until next checkpoint is reached
        while simulation_single_node.current.time <
                simulation_single_node.time_mesh.time_nodes[i_checkpoint]
            # (1) set values in ghost layer

            # (2) Compute numerical fluxes

            # (3) Calculate the new time step (maximum allows time_step due to
            # the wave speeds)
            time_step = 0.5

            # (4) Update the cell values Gudonov style first order

            # (5) Refresh to current simulation time
            simulation_single_node.current.time += time_step
            
        end

        # Save the fields
        write_fields!(nc_data_set, simulation_single_node, i_checkpoint)
    end

    # Close the connection to the dataset handle
    close_output_file(nc_data_set)




end


main()