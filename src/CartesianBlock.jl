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
    # Boundary cell mimics interior cell next to it and reflects the perpendicular momentum
    WALL
    # Boundary cell mimics interior cell next to it in all conserved quantities
    OUTFLOW
    # Boundary cell values are copied from the interior values of the adjacent block
    CONNECT
    # The Boundary is not updated at all
    PASSIVE
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
    mesh_x::Array{Float64, 1}
    mesh_y::Array{Float64, 1}
end

mutable struct Simulation_Block
    layout::Layout
    fields::SWE_Fields
    time::Float64
    bathymetry::Array{Float64, 2}
    boundaries::Boundary_Collection
end