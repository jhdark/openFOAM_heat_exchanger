from festim_script import read_openfoam_data
from dolfinx.io import VTXWriter
from mpi4py import MPI


def export_openfoam_data(T_cold, T_hot, u_cold, u_hot):
    """
    Export OpenFOAM data to VTX files.
    """

    print("Exporting OpenFOAM data")
    writer_T_cold = VTXWriter(
        MPI.COMM_WORLD,
        "openfoam_data/temperature_cold.bp",
        T_cold,
        "BP5",
    )
    writer_T_hot = VTXWriter(
        MPI.COMM_WORLD,
        "openfoam_data/temperature_hot.bp",
        T_hot,
        "BP5",
    )
    writer_u_cold = VTXWriter(
        MPI.COMM_WORLD,
        "openfoam_data/u_cold.bp",
        u_cold,
        "BP5",
    )
    writer_u_hot = VTXWriter(
        MPI.COMM_WORLD,
        "openfoam_data/u_hot.bp",
        u_hot,
        "BP5",
    )

    writer_T_cold.write(t=0)
    writer_T_hot.write(t=0)
    writer_u_cold.write(t=0)
    writer_u_hot.write(t=0)


if __name__ == "__main__":
    # read openfoam data
    T_cold, T_hot, u_cold, u_hot = read_openfoam_data()

    export_openfoam_data(T_cold, T_hot, u_cold, u_hot)
