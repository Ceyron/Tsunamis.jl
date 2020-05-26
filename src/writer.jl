#=
Implements routines to write the results of the shallow water simulation to a netcdf file that is readable by Paraview
=#

import NCDatasets
import NetCDF

include("single_node_setup.jl")

# Creates the nc dataset to which we can write
function create_output_file(
    simulation_settings::SWE_Simulation_Settings,
    simulation_data::SWE_Simulation,
    )

    # Interpolate the file name
    nc_file_name = "$(simulation_settings.output_file_name).nc"

    # If the file already exists, ask the user to overwrite it
    if isfile(nc_file_name)
        println("WARNING: The file $nc_file_name already exists")
        println("WARNING: Press enter to overwrite it")
        readline()
        rm(nc_file_name)
    end

    x_in_ds = (deepcopy(simulation_data.current.layout.mesh_x) .-
        simulation_data.current.layout.cell_width_x ./ 2.)[1:end]
    y_in_ds = (deepcopy(simulation_data.current.layout.mesh_y) .-
        simulation_data.current.layout.cell_width_y ./ 2.)[1:end]
    time_in_ds = deepcopy(simulation_data.time_mesh.time_nodes)

    # Attribute Dictionaries for dimensions
    x_attributes = Dict(
        "longname" => "x coordinate",
        "units" => "meter"
    )
    y_attributes = Dict(
        "longname" => "y coordinate",
        "units" => "meter"
    )
    time_attributes = Dict(
        "longname" => "Time",
        "units" => "second"
    )
    # Saved variables attributes
    h_attributes = Dict(
        "longname" => "water height above bathymetry",
        "units" => "m"
    )
    hu_attributes = Dict(
        "longname" => "x momentum",
        "units" => "m*m/s"
    )
    hv_attributes = Dict(
        "longname" => "y momentum",
        "units" => "m*m/s"
    )
    b_attributes = Dict(
        "longname" => "bathymetry",
        "units" => "m"
    )
    NetCDF.nccreate(
        nc_file_name,
        "h",
        "x",
        x_in_ds,
        x_attributes,
        "y",
        y_in_ds,
        y_attributes,
        "time",
        time_in_ds,
        time_attributes,
        atts=h_attributes
    )
    NetCDF.nccreate(
        nc_file_name,
        "hu",
        "x",
        x_in_ds,
        x_attributes,
        "y",
        y_in_ds,
        y_attributes,
        "time",
        time_in_ds,
        time_attributes,
        atts=hu_attributes
    )
    NetCDF.nccreate(
        nc_file_name,
        "hv",
        "x",
        x_in_ds,
        x_attributes,
        "y",
        y_in_ds,
        y_attributes,
        "time",
        time_in_ds,
        time_attributes,
        atts=hv_attributes
    )
    NetCDF.nccreate(
        nc_file_name,
        "b",
        "x",
        x_in_ds,
        x_attributes,
        "y",
        y_in_ds,
        y_attributes,
        atts=b_attributes
    )

    # Now use the other library to read in the dataset
    nc_data_set = NCDatasets.Dataset(nc_file_name, "a")
    @views nc_data_set["b"][:, :] =
        simulation_data.current.bathymetry[2:end-1, 2:end-1]

    return nc_data_set
end

# Write the current fields in the netcdf output
function write_fields!(
    data_set_handle,
    simulation_data::SWE_Simulation,
    time_index::Int
    )

    # Define the slicer to extract the interior cells
    slicer_x = 2:simulation_data.current.layout.num_interior_cells_x+1
    slicer_y = 2:simulation_data.current.layout.num_interior_cells_y+1

    # Write the fields to the open netcdf handle
    @views data_set_handle["h"][:, :, time_index] =
        simulation_data.current.fields.h[slicer_x, slicer_y]
    @views data_set_handle["hu"][:, :, time_index] =
        simulation_data.current.fields.hu[slicer_x, slicer_y]
    @views data_set_handle["hv"][:, :, time_index] =
        simulation_data.current.fields.hv[slicer_x, slicer_y]
end


# Properly closes the handle
function close_output_file(
    data_set_handle
)
    NCDatasets.close(data_set_handle)
end