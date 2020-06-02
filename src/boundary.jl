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
        @. @views simulation_data.current.fields.h[1, 2:end-1] =
            simulation_data.current.fields.h[2, 2:end-1]
        # Momentum in x (mirrored)
        @. @views simulation_data.current.fields.hu[1, 2:end-1] .=
            - simulation_data.current.fields.hu[2, 2:end-1]
        # Momentum in y (copied)
        @. @views simulation_data.current.fields.hv[1, 2:end-1] =
            simulation_data.current.fields.hv[2, 2:end-1]
    elseif simulation_data.current.boundaries.left.type == OUTFLOW
        # Height (copied)
        @. @views simulation_data.current.fields.h[1, 2:end-1] =
            simulation_data.current.fields.h[2, 2:end-1]
        # Momentum in x (copied)
        @. @views simulation_data.current.fields.hu[1, 2:end-1] =
            simulation_data.current.fields.hu[2, 2:end-1]
        # Momentum in y (copied)
        @. @views simulation_data.current.fields.hv[1, 2:end-1] =
            simulation_data.current.fields.hv[2, 2:end-1]
    elseif simulation_data.current.boundaries.left.type == CONNECT
        recv_left::SWE_Copy_Fields = take!(
            simulation_data.current.boundaries.left.connector_recv
        )
        # Height (from left domain)
        @. @views simulation_data.current.fields.h[1, 2:end-1] =
            recv_left.h
        # Momentum in x (from left domain)
        @. @views simulation_data.current.fields.hu[1, 2:end-1] =
            recv_left.hu
        # Momentum in y (from left domain)
        @. @views simulation_data.current.fields.hv[1, 2:end-1] =
            recv_left.hv
    elseif simulation_data.current.boundaries.left.type == PASSIVE
        # Nothing happens
    end

    # TOP
    if simulation_data.current.boundaries.top.type == WALL
        # Height (copied)
        @. @views simulation_data.current.fields.h[2:end-1, end] =
            simulation_data.current.fields.h[2:end-1, end-1]
        # Momentum x (copied)
        @. @views simulation_data.current.fields.hu[2:end-1, end] =
            simulation_data.current.fields.hu[2:end-1, end-1]
        # Momentum y (mirrored)
        @. @views simulation_data.current.fields.hv[2:end-1, end] .=
            - simulation_data.current.fields.hv[2:end-1, end-1]
    elseif simulation_data.current.boundaries.top.type == OUTFLOW
        # Height (copied)
        @. @views simulation_data.current.fields.h[2:end-1, end] =
            simulation_data.current.fields.h[2:end-1, end-1]
        # Momentum x (copied)
        @. @views simulation_data.current.fields.hu[2:end-1, end] =
            simulation_data.current.fields.hu[2:end-1, end-1]
        # Momentum y (copied)
        @. @views simulation_data.current.fields.hv[2:end-1, end] =
            simulation_data.current.fields.hv[2:end-1, end-1]
    elseif simulation_data.current.boundaries.top.type == CONNECT
        recv_top::SWE_Copy_Fields = take!(
            simulation_data.current.boundaries.top.connector_recv
        )
        # Height (from left domain)
        @. @views simulation_data.current.fields.h[2:end-1, end] =
            recv_top.h
        # Momentum in x (from left domain)
        @. @views simulation_data.current.fields.hu[2:end-1, end] =
            recv_top.hu
        # Momentum in y (from left domain)
        @. @views simulation_data.current.fields.hv[2:end-1, end] =
            recv_top.hv
    elseif simulation_data.current.boundaries.top.type == PASSIVE
        # Nothing happens
    end

    # RIGHT
    if simulation_data.current.boundaries.right.type == WALL
        # Height (copied)
        @. @views simulation_data.current.fields.h[end, 2:end-1] =
            simulation_data.current.fields.h[end-1, 2:end-1]
        # Momentum x (mirrored)
        @. @views simulation_data.current.fields.hu[end, 2:end-1] .=
            - simulation_data.current.fields.hu[end-1, 2:end-1]
        # Momentum y (copied)
        @. @views simulation_data.current.fields.hv[end, 2:end-1] =
            simulation_data.current.fields.hv[end-1, 2:end-1]
    elseif simulation_data.current.boundaries.right.type == OUTFLOW
        # Height (copied)
        @. @views simulation_data.current.fields.h[end, 2:end-1] =
            simulation_data.current.fields.h[end-1, 2:end-1]
        # Momentum x (copied)
        @. @views simulation_data.current.fields.hu[end, 2:end-1] =
            simulation_data.current.fields.hu[end-1, 2:end-1]
        # Momentum y (copied)
        @. @views simulation_data.current.fields.hv[end, 2:end-1] =
            simulation_data.current.fields.hv[end-1, 2:end-1]
    elseif simulation_data.current.boundaries.right.type == CONNECT
        recv_right::SWE_Copy_Fields = take!(
            simulation_data.current.boundaries.right.connector_recv
        )
        # Height (from left domain)
        @. @views simulation_data.current.fields.h[end, 2:end-1] =
            recv_right.h
        # Momentum in x (from left domain)
        @. @views simulation_data.current.fields.hu[end, 2:end-1] =
            recv_right.hu
        # Momentum in y (from left domain)
        @. @views simulation_data.current.fields.hv[end, 2:end-1] =
            recv_right.hv
    elseif simulation_data.current.boundaries.right.type == PASSIVE
        # Nothing happens
    end

    # BOTTOM
    if simulation_data.current.boundaries.bottom.type == WALL
        # Height (copied)
        @. @views simulation_data.current.fields.h[2:end-1, 1] =
            simulation_data.current.fields.h[2:end-1, 2]
        # Momentum x (copied)
        @. @views simulation_data.current.fields.hu[2:end-1, 1] =
            simulation_data.current.fields.hu[2:end-1, 2]
        # Momentum y (mirrored)
        @. @views simulation_data.current.fields.hv[2:end-1, 1] .= 
            - simulation_data.current.fields.hv[2:end-1, 2]
    elseif simulation_data.current.boundaries.bottom.type == OUTFLOW
        # Height (copied)
        @. @views simulation_data.current.fields.h[2:end-1, 1] =
            simulation_data.current.fields.h[2:end-1, 2]
        # Momentum x (copied)
        @. @views simulation_data.current.fields.hu[2:end-1, 1] =
            simulation_data.current.fields.hu[2:end-1, 2]
        # Momentum y (copied)
        @. @views simulation_data.current.fields.hv[2:end-1, 1] =
            simulation_data.current.fields.hv[2:end-1, 2]
    elseif simulation_data.current.boundaries.bottom.type == CONNECT
        recv_bottom::SWE_Copy_Fields = take!(
            simulation_data.current.boundaries.bottom.connector_recv
        )
        # Height (from left domain)
        @. @views simulation_data.current.fields.h[2:end-1, 1] =
            recv_bottom.h
        # Momentum in x (from left domain)
        @. @views simulation_data.current.fields.hu[2:end-1, 1] =
            recv_bottom.hu
        # Momentum in y (from left domain)
        @. @views simulation_data.current.fields.hv[2:end-1, 1] =
            recv_bottom.hv
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

# Enqueues the copied layer in the channel connecting it to the neighboring
# domain
function queue_copy_layer(
    simulation_data::SWE_Simulation
)
    # LEFT
    if simulation_data.current.boundaries.left.type == CONNECT
        to_send_left = SWE_Copy_Fields(
            simulation_data.current.fields.h[2, 2:end-1],
            simulation_data.current.fields.hu[2, 2:end-1],
            simulation_data.current.fields.hv[2, 2:end-1],
        )
        put!(simulation_data.current.boundaries.left.connector_send, to_send_left)
    end

    # TOP
    if simulation_data.current.boundaries.top.type == CONNECT
        to_send_top = SWE_Copy_Fields(
            simulation_data.current.fields.h[2:end-1, end-1],
            simulation_data.current.fields.hu[2:end-1, end-1],
            simulation_data.current.fields.hv[2:end-1, end-1],
        )
        put!(simulation_data.current.boundaries.top.connector_send, to_send_top)
    end

    # LEFT
    if simulation_data.current.boundaries.right.type == CONNECT
        to_send_right = SWE_Copy_Fields(
            simulation_data.current.fields.h[end-1, 2:end-1],
            simulation_data.current.fields.hu[end-1, 2:end-1],
            simulation_data.current.fields.hv[end-1, 2:end-1],
        )
        put!(simulation_data.current.boundaries.right.connector_send, to_send_right)
    end

    # LEFT
    if simulation_data.current.boundaries.bottom.type == CONNECT
        to_send_bottom = SWE_Copy_Fields(
            simulation_data.current.fields.h[2:end-1, 2],
            simulation_data.current.fields.hu[2:end-1, 2],
            simulation_data.current.fields.hv[2:end-1, 2],
        )
        put!(simulation_data.current.boundaries.bottom.connector_send, to_send_bottom)
    end
end