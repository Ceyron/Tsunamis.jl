"""
Implements the main data structure and all routines of the main simulation loop
"""


# Contains the 3 conserved field variables (that depend on 2D-space and time) as
# well as the bathymetry (this also helps to represent the version of the SWE
# that can work without source terms -> transferring the source terms into a
# separate variable in the hyperbolic equation)
mutable struct SWE_Fields
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

mutable struct Simulation_Block
    fields::SWE_Fields
    bathymetry::Array{Float64, 2}
    left::Boundary
    top::Boundary
    right::Boundary
    bottom::Boundary
end