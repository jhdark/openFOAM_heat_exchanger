import meshio


def convert_med_to_xdmf(
    medfilename,
    cell_file="mesh_domains.xdmf",
    facet_file="mesh_boundaries.xdmf",
    cell_type="tetra",
    facet_type="triangle",
):
    """_summary_
    Args:
        medfilename (_type_): _description_
        cell_file (str, optional): _description_. Defaults to "mesh_domains.xdmf".
        facet_file (str, optional): _description_. Defaults to "mesh_boundaries.xdmf".
        cell_type (str, optional): _description_. Defaults to "tetra".
        facet_type (str, optional): _description_. Defaults to "triangle".
    Returns:
        dict, dict: the correspondance dict, the cell types
    """
    msh = meshio.read(medfilename)

    correspondance_dict = msh.cell_tags

    cell_data_types = msh.cell_data_dict["cell_tags"].keys()

    for mesh_block in msh.cells:
        if mesh_block.type == cell_type:
            meshio.write_points_cells(
                cell_file,
                msh.points,
                [mesh_block],
                cell_data={"f": [-1 * msh.cell_data_dict["cell_tags"][cell_type]]},
            )
        elif mesh_block.type == facet_type:
            meshio.write_points_cells(
                facet_file,
                msh.points,
                [mesh_block],
                cell_data={"f": [-1 * msh.cell_data_dict["cell_tags"][facet_type]]},
            )

    return correspondance_dict, cell_data_types


if __name__ == "__main__":

    correspondance_dict, cell_data_types = convert_med_to_xdmf(
        medfilename="Mesh_heat_exchanger.med",
        cell_file="mesh_domains.xdmf",
        facet_file="mesh_boundaries.xdmf",
        cell_type="tetra",
        facet_type="triangle",
    )

    print(correspondance_dict)