import festim as F
from foam2dolfinx import OpenFOAMReader
from dolfinx import fem
from festim.helpers import nmm_interpolate
from dolfinx.io import VTXWriter
from mpi4py import MPI
import adios4dolfinx
from dolfinx import log


def read_openfoam_data():
    """
    Read OpenFOAM data from a file and return the temperature and velocity fields.
    """
    print("Reading OpenFOAM data...")
    openfoam_reader = OpenFOAMReader(
        filename="openfoam_data/heatexchanger.foam", cell_type=10
    )
    T_cold = openfoam_reader.create_dolfinx_function(t=700, name="T", subdomain="cold")
    T_hot = openfoam_reader.create_dolfinx_function(t=700, name="T", subdomain="hot")
    u_cold = openfoam_reader.create_dolfinx_function(t=700, name="U", subdomain="cold")
    u_hot = openfoam_reader.create_dolfinx_function(t=700, name="U", subdomain="hot")

    return T_cold, T_hot, u_cold, u_hot


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


def build_festim_model(T_cold, T_hot, u_cold, u_hot):
    print("Building FESTIM model...")
    id_hot = 6
    id_cold = 7
    if_hot_inlet = 9
    if_cold_inlet = 8

    my_model = F.HydrogenTransportProblem()

    my_mesh = F.MeshFromXDMF(
        volume_file="mesh_files/mesh_domains.xdmf",
        facet_file="mesh_files/mesh_boundaries.xdmf",
    )
    my_model.mesh = my_mesh
    volume_meshtags = my_mesh.define_volume_meshtags()

    H = F.Species(name="H", mobile=True)
    my_model.species = [H]

    mat_1 = F.Material(name="mat_1", D_0=1e-06, E_D=0)

    vol_hot = F.VolumeSubdomain(id=id_hot, material=mat_1)
    vol_cold = F.VolumeSubdomain(id=id_cold, material=mat_1)
    hot_inlet = F.SurfaceSubdomain(id=if_hot_inlet)
    cold_inlet = F.SurfaceSubdomain(id=if_cold_inlet)

    my_model.subdomains = [vol_hot, vol_cold, hot_inlet, cold_inlet]

    # interpolate temperature fields onto one function
    print("Interpolating temperature field")
    V = fem.functionspace(my_mesh.mesh, ("P", 1))
    my_temperature_field = fem.Function(V)
    nmm_interpolate(
        my_temperature_field,
        T_cold,
        cells=volume_meshtags.find(id_cold),
        padding=1e-4,
    )
    nmm_interpolate(
        my_temperature_field,
        T_hot,
        cells=volume_meshtags.find(id_hot),
        padding=1e-4,
    )
    my_model.temperature = my_temperature_field

    hot_inlet_conc = F.FixedConcentrationBC(subdomain=hot_inlet, value=0, species=H)
    cold_inlet_conc = F.FixedConcentrationBC(subdomain=cold_inlet, value=1e5, species=H)
    my_model.boundary_conditions = [hot_inlet_conc, cold_inlet_conc]

    hot_adv = F.AdvectionTerm(velocity=u_hot, subdomain=vol_hot, species=H)
    cold_adv = F.AdvectionTerm(velocity=u_cold, subdomain=vol_cold, species=H)
    my_model.advection_terms = [hot_adv, cold_adv]

    dt = F.Stepsize(
        initial_value=1, growth_factor=1.1, cutback_factor=0.9, target_nb_iterations=5
    )
    my_model.settings = F.Settings(
        atol=1e-08,
        rtol=1e-10,
        max_iterations=30,
        transient=True,
        final_time=1500,
        stepsize=dt,
    )

    mobile = F.VTXSpeciesExport(filename="results/mobile_concentration.bp", field=H)
    my_model.exports = [mobile]

    return my_model


if __name__ == "__main__":
    # read openfoam data
    T_cold, T_hot, u_cold, u_hot = read_openfoam_data()

    # export_openfoam_data(T_cold, T_hot, u_cold, u_hot)

    # build FESTIM model
    my_model = build_festim_model(T_cold, T_hot, u_cold, u_hot)

    # log.set_log_level(log.LogLevel.INFO)

    my_model.initialise()
    my_model.run()
