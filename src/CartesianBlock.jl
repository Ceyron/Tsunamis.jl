"""
Implements the main data structure and all routines of the main simulation loop
"""


mutable struct SWE_Fields
# Contains the 3 conserved field variables (that depend on 2D-space and time) as
# well as the bathymetry (this also helps to represent the version of the SWE
# that can work without source terms -> transferring the source terms into a
# separate variable in the hyperbolic equation)
#
# The simulation domain is 2d with halo cells around to support for boundaries and decomposed synchronization. 
    h::Array{Float64, 2}
    hu::Array{Float64, 2}
    hv::Array{Float64, 2}
end

@enum Boundary_Type begin
    periodic
    connected
end

struct Boundary
    type::Boundary_Type
    connector::String
end

# Contains all the 4 edges of a Cartesian Block
struct Boundary_Collection
    left::Boundary
    top::Boundary
    right::Boundary
    bottom::Boundary
end

# Defines the structure of a simulation Block
#
#
#
#
#
#
#  ^y
#  |
#  --> x
struct Layout
    offset_x::Float64
    offset_y::Float64
    size_x::Float64
    size_y::Float64
    num_interior_cells_x::UInt64
    num_interior_cells_y::UInt64
    cell_width_x::Float64
    cell_width_y::Float64
end

mutable struct Simulation_Block
    layout::Layout
    fields::SWE_Fields
    time::Float64
    bathymetry::Array{Float64, 2}
    boundaries::Boundary_Collection
end

# Contains 4 arrays of the same size as the simulation domain that contain the
# numerical fluxes at the respective edges
mutable struct Flux_Field
    left::Array{Float64, 2}
    top::Array{Float64, 2}
    right::Array{Float64, 2}
    bottom::Array{Float64, 2}
end

# Contains all the fluxes for a simulation block
mutable struct Fluxes
    h::Flux_Field
    hu::Flux_Field
    hv::Flux_Field
    max_wave_speed::Float64
end