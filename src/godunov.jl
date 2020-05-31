# Contains routines to update the simulation block in a godunov style fashion by
# a first order explicit time integration
include("single_node_setup.jl")
include("CartesianBlock.jl")
include("hlle.jl")

function calculate_numerical_fluxes!(
    numerical_fluxes::SWE_Fields,
    simulation_data::SWE_Simulation
    )
    # Initialize and collect over all over cells
    max_wave_speed = 0.0
    
    # Iterate over each and calculate the flux in each direction perpendicular
    # to the edge, then add this to the summed fluxes for each cell
    #
    # One could understand this iteration as if we are iterating over the
    # interior cells and bottom row and left column on halo cells and considering
    # the right and top edge of these cells
    for j in 1:simulation_data.current.layout.num_interior_cells_y+1
        for i in 1:simulation_data.current.layout.num_interior_cells_x+1
            # The quantities to be calculated for each vertical edge
            flux_h_left = 0.0
            flux_h_right = 0.0
            flux_hu_left = 0.0
            flux_hu_right = 0.0

            # The quantities to be calculated for each horizontal edge
            flux_h_top = 0.0
            flux_h_bottom = 0.0
            flux_hv_top = 0.0
            flux_hv_bottom = 0.0

            # The maximum wave speed at the edge currently considered by the Riemann solver
            max_edge_speed = 0.0

            ####
            # Calculate right edge
            ####

            # Solve Riemann Problem
            @inbounds flux_h_left, flux_h_right,
            flux_hu_left, flux_hu_right,
            max_edge_speed = solve_riemann_hlle(
                simulation_data.current.fields.h[i, j],
                simulation_data.current.fields.h[i+1, j],
                simulation_data.current.fields.hu[i, j],
                simulation_data.current.fields.hu[i+1, j],
                simulation_data.current.bathymetry[i, j],
                simulation_data.current.bathymetry[i+1, j],
            )

            # Update the accumulated fluxes of the cells adjacent to the right
            # edge (i.e. the current aiding iterating cell and the one right of
            # it)
            @inbounds numerical_fluxes.h[i, j] +=
                flux_h_left / simulation_data.current.layout.cell_width_x
            @inbounds numerical_fluxes.h[i+1, j] +=
                flux_h_right / simulation_data.current.layout.cell_width_x
            @inbounds numerical_fluxes.hu[i, j] +=
                flux_hu_left / simulation_data.current.layout.cell_width_x
            @inbounds numerical_fluxes.hu[i+1, j] +=
                flux_hu_right / simulation_data.current.layout.cell_width_x

            # Update wave speed
            max_wave_speed = max(max_wave_speed, max_edge_speed, )


            ####
            # Calculate top edge
            ####

            # Solve Riemann Problem
            @inbounds flux_h_bottom, flux_h_top,
            flux_hv_bottom, flux_hv_top,
            max_edge_speed = solve_riemann_hlle(
                simulation_data.current.fields.h[i, j],
                simulation_data.current.fields.h[i, j+1],
                simulation_data.current.fields.hv[i, j],
                simulation_data.current.fields.hv[i, j+1],
                simulation_data.current.bathymetry[i, j],
                simulation_data.current.bathymetry[i, j+1],
            )

            # Update the accumulated fluxes of the cells adjacent to the top
            # edge (i.e. the current aiding iterating cell and the one on top of
            # it)
            # TODO: Don't we have to divide by cell_width_x because we need the length of the edge?!
            @inbounds numerical_fluxes.h[i, j] +=
                flux_h_bottom / simulation_data.current.layout.cell_width_y
            @inbounds numerical_fluxes.h[i, j+1] +=
                flux_h_top / simulation_data.current.layout.cell_width_y
            @inbounds numerical_fluxes.hv[i, j] +=
                flux_hv_bottom / simulation_data.current.layout.cell_width_y
            @inbounds numerical_fluxes.hv[i, j+1] +=
                flux_hv_top / simulation_data.current.layout.cell_width_y

            # Update wave speed
            max_wave_speed = max(max_wave_speed, max_edge_speed, )
        end
    end

    return max_wave_speed
end

function compute_max_time_step(
    simulation_data::SWE_Simulation,
    max_wave_speed::Float64,
    cfl_condition::Float64,
)
    if max_wave_speed < 0.00000
        # This shouldn't be the case
        println("WARNING: Max wave speed is too low")
    end
    # The CFL condition ensures that the numerical solution is stable
    max_time_step = min(
        simulation_data.current.layout.cell_width_x / max_wave_speed,
        simulation_data.current.layout.cell_width_y / max_wave_speed,
    )

    max_time_step *= cfl_condition

    return max_time_step
end

# Used the accumulated fluxes from all surrounding edges of a cell to update the
# value. Requires a purposefully chosen time_step to be numerically stable
function update_cells!(
    simulation_data::SWE_Simulation,
    numerical_fluxes::SWE_Fields,
    time_step::Float64,
    )

    # Loop over all interior cells
    for j in 2:simulation_data.current.layout.num_interior_cells_y+1
        for i in 2:simulation_data.current.layout.num_interior_cells_x+1
            # Perform the explicit integration
            @inbounds simulation_data.current.fields.h[i, j] -=
                time_step * numerical_fluxes.h[i, j]
            @inbounds simulation_data.current.fields.hu[i, j] -=
                time_step * numerical_fluxes.hu[i, j]
            @inbounds simulation_data.current.fields.hv[i, j] -=
                time_step * numerical_fluxes.hv[i, j]
        end
    end
end
