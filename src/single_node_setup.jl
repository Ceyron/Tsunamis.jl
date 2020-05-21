include("CartesianBlock.jl")

struct Time_Mesh
    time_end::Float64
    num_checkpoints::Int64
    time_nodes::Array{Float64, 1}
end

# The main simulation data structure
#
# Contains the current struct which is mutated over the course of the simulation
# as well the time information
mutable struct SWE_Simulation
    current::Simulation_Block
    time_mesh::Time_Mesh
end

struct SWE_Simulation_Settings
    output_file_name::AbstractString
end