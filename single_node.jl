#=
   Contains a single-node sequential implementation of an SWE solver for the simulation of tsunamis
=#

import Plots
Plots.plotlyjs()

include("src/CartesianBlock.jl")

# The main simulation data structure
#
# Contains the current struct which is mutated over the course of the simulation
# as well as the required history of the saved structs over the course of the
# simulation
mutable struct SWE_Simulation
    current::Simulation_Block
    checkpoints::Array{SWE_Fields, 1}
end


const size_x = 50;
const size_y = 50;
const time_end = 15;
const num_checkpoints = 20;

function main()
    println("Hello World");

    # Instantiate main block
    simulation_single_node = SWE_Simulation(
        Simulation_Block(
            SWE_Fields(
                # Account for the additional halo layers of the boundary
                zeros(size_x+2, size_y+2),
                zeros(size_x+2, size_y+2),
                zeros(size_x+2, size_y+2),
            ),
            zeros(size_x+2, size_y+2),
            # No particular distribution occurs in the single_node case,
            # therefore set all boundary conditions to periodic
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
        Array{SWE_Fields, 1}(undef, num_checkpoints),
    )

    # Imprint the initial condition
    simulation_single_node.current.fields.h[22:28, 22:28] .= 3.0

    display(Plots.surface(1:52, 1:52, simulation_single_node.current.fields.h))
    readline()

end


main()