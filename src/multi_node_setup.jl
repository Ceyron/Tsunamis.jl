@everywhere using Distributed
@everywhere import NCDatasets

@everywhere include("CartesianBlock.jl")
@everywhere include("single_node_setup.jl")


struct Time_Mesh
    time_end::Float64
    num_checkpoints::Int64
    time_nodes::Array{Float64, 1}
end

struct Block_Mesh
    block_mesh_x::Array{Float64, 1}
    block_mesh_y::Array{Float64, 1}
    number_of_blocks_x::Int
    number_of_blocks_y::Int
    number_of_blocks_total::Int
end

# The main simulation data structure
#
# Contains the current struct which is mutated over the course of the simulation
# as well the time information
#
# In the distributed case we safe references (=Futures) to the Simulation Data
# Structures on the other cores based on their block position
mutable struct SWE_Simulation_Distributed
    block_references::Array{Future, 2}
    time_mesh::Time_Mesh
    time::Float64
    block_mesh::Block_Mesh
    nc_data_sets::Array{Future, 2}
end

struct SWE_Simulation_Settings
    output_file_name::AbstractString
    clf_number::Float64
end