# A scenario with no bathymetry involved and peak in the center of the domain

include("../src/single_node_setup.jl")

const domain_size_x = 1000.0
const domain_size_y = 1000.0
const center_radius = 100.0
const position_center_x = domain_size_x / 2.0
const position_center_y = domain_size_y / 2.0

const height_water_center = 15.0
const height_water_else = 10.0

function radial_dam_break_get_domain_size()
    return domain_size_x, domain_size_y
end

function radial_dam_break_get_boundary_type()
    return OUTFLOW
end

function radial_dam_break_imprint_initial_condition!(
    simulation_data::SWE_Simulation,
)
    # Set everywhere to the standard quantity (including the boundary)
    simulation_data.current.fields.h .= height_water_else

    # Check for which cells lie in the central region that has elevated water
    # height
    for i in 2:simulation_data.current.layout.num_interior_cells_x+1
        for j in 2:simulation_data.current.layout.num_interior_cells_y+1
            # Calculate the center position of the current interior cell
            position_x = simulation_data.current.layout.mesh_x[i-1] +
                0.5 * simulation_data.current.layout.cell_width_x
            position_y = simulation_data.current.layout.mesh_y[j-1] +
                0.5 * simulation_data.current.layout.cell_width_y


            in_center_radius = sqrt(
                (position_x - position_center_x)^2 +
                (position_y - position_center_y)^2
            ) < center_radius

            if in_center_radius
                simulation_data.current.fields.h[i, j] = height_water_center
            end
        end
    end
end

function radial_dam_break_imprint_bathymetry!(
    simulation_data::SWE_Simulation
)
    simulation_data.current.bathymetry .= 0
end