# Comprises routines to update the boundary conditions

include("CartesianBlock.jl")
include("single_node_setup.jl")

function update_boundaries!(
    simulation_data::SWE_Simulation
)
    # Iterate over the boundary edges and set them according to the type

    # LEFT
    if simulation_data.current.boundaries.left.type == WALL
        # Height (copied)
        simulation_data.current.fields.h[1, 2:end-1] =
            simulation_data.current.fields.h[1, 2:end-1]
        # Momentum in x (mirrored)
        simulation_data.current.fields.hu[1, 2:end-1] = -
            simulation_data.current.fields.hu[1, 2:end-1]
        # Momentum in y (copied)
        simulation_data.current.fields.hv[1, 2:end-1] =
            simulation_data.current.fields.hv[1, 2:end-1]
    elseif simulation_data.current.boundaries.left.type == OUTFLOW
        # Height (copied)
        simulation_data.current.fields.h[1, 2:end-1] =
            simulation_data.current.fields.h[1, 2:end-1]
        # Momentum in x (copied)
        simulation_data.current.fields.hu[1, 2:end-1] =
            simulation_data.current.fields.hu[1, 2:end-1]
        # Momentum in y (copied)
        simulation_data.current.fields.hv[1, 2:end-1] =
            simulation_data.current.fields.hv[1, 2:end-1]
    elseif simulation_data.current.boundaries.left.type == PASSIVE
        # Nothing happens
    end

    # TOP
    if simulation_data.current.boundaries.top.type == WALL
        # Height (copied)
        simulation_data.current.fields.h[2:end-1, end] =
            simulation_data.current.fields.h[2:end-1, end-1]
        # Momentum x (copied)
        simulation_data.current.fields.hu[2:end-1, end] =
            simulation_data.current.fields.hu[2:end-1, end-1]
        # Momentum y (mirrored)
        simulation_data.current.fields.hv[2:end-1, end] = -
            simulation_data.current.fields.hv[2:end-1, end]
    elseif simulation_data.current.boundaries.top.type == OUTFLOW
        # Height (copied)
        simulation_data.current.fields.h[2:end-1, end] =
            simulation_data.current.fields.h[2:end-1, end-1]
        # Momentum x (copied)
        simulation_data.current.fields.hu[2:end-1, end] =
            simulation_data.current.fields.hu[2:end-1, end-1]
        # Momentum y (copied)
        simulation_data.current.fields.hv[2:end-1, end] =
            simulation_data.current.fields.hv[2:end-1, end]
    elseif simulation_data.current.boundaries.top.type == PASSIVE
        # Nothing happens
    end

    # RIGHT
    if simulation_data.current.boundaries.right.type == WALL
        # Height (copied)
        simulation_data.current.fields.h[end, 2:end-1] =
            simulation_data.current.fields.h[end-1, 2:end-1]
        # Momentum x (mirrored)
        simulation_data.current.fields.hu[end, 2:end-1] = -
            simulation_data.current.fields.hu[end-1, 2:end-1]
        # Momentum y (copied)
        simulation_data.current.fields.hv[end, 2:end-1] =
            simulation_data.current.fields.hv[end-1, 2:end-1]
    elseif simulation_data.current.boundaries.right.type == OUTFLOW
        # Height (copied)
        simulation_data.current.fields.h[end, 2:end-1] =
            simulation_data.current.fields.h[end-1, 2:end-1]
        # Momentum x (copied)
        simulation_data.current.fields.hu[end, 2:end-1] =
            simulation_data.current.fields.hu[end-1, 2:end-1]
        # Momentum y (copied)
        simulation_data.current.fields.hv[end, 2:end-1] =
            simulation_data.current.fields.hv[end-1, 2:end-1]
    elseif simulation_data.current.boundaries.right.type == PASSIVE
        # Nothing happens
    end

    # BOTTOM
    if simulation_data.current.boundaries.bottom.type == WALL
        # Height (copied)
        simulation_data.current.fields.h[2:end-1, 1] =
            simulation_data.current.fields.h[2:end-1, 2]
        # Momentum x (copied)
        simulation_data.current.fields.hu[2:end-1, 1] =
            simulation_data.current.fields.hu[2:end-1, 2]
        # Momentum y (mirrored)
        simulation_data.current.fields.hv[2:end-1, 1] = -
            simulation_data.current.fields.hv[2:end-1, 2]
    elseif simulation_data.current.boundaries.bottom.type == OUTFLOW
        # Height (copied)
        simulation_data.current.fields.h[2:end-1, 1] =
            simulation_data.current.fields.h[2:end-1, 2]
        # Momentum x (copied)
        simulation_data.current.fields.hu[2:end-1, 1] =
            simulation_data.current.fields.hu[2:end-1, 2]
        # Momentum y (copied)
        simulation_data.current.fields.hv[2:end-1, 1] =
            simulation_data.current.fields.hv[2:end-1, 2]
    elseif simulation_data.current.boundaries.bottom.type == PASSIVE
        # Nothing happens
    end

    # Choose the corner boundary cells to achieve a zero Riemann solution

    # BOTTOM LEFT
    simulation_data.current.fields.h[1, 1] =
        simulation_data.current.fields.h[2, 2]
    simulation_data.current.fields.hu[1, 1] =
        simulation_data.current.fields.h[2, 2]
    simulation_data.current.fields.hv[1, 1] =
        simulation_data.current.fields.h[2, 2]
    
    # TOP LEFT
    simulation_data.current.fields.h[1, end] =
        simulation_data.current.fields.h[2, end-1]
    simulation_data.current.fields.hu[1, end] =
        simulation_data.current.fields.hu[2, end-1]
    simulation_data.current.fields.hv[1, end] =
        simulation_data.current.fields.hv[2, end-1]
    
    # TOP RIGHT
    simulation_data.current.fields.h[end, end] =
        simulation_data.current.fields.h[end-1, end-1]
    simulation_data.current.fields.hu[end, end] =
        simulation_data.current.fields.hu[end-1, end-1]
    simulation_data.current.fields.hv[end, end] =
        simulation_data.current.fields.hv[end-1, end-1]
    
    # BOTTOM RIGHT
    simulation_data.current.fields.h[end, 1] =
        simulation_data.current.fields.h[end-1, 2]
    simulation_data.current.fields.hu[end, 1] =
        simulation_data.current.fields.hu[end-1, 2]
    simulation_data.current.fields.hv[end, 1] =
        simulation_data.current.fields.hv[end-1, 2]
end