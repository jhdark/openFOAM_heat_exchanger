import festim as F
from foam2dolfinx import OpenFOAMReader
import numpy as np
from dolfinx import fem
from festim.helpers import nmm_interpolate
import basix
import ufl
from dolfinx.io import VTXWriter
from mpi4py import MPI
import adios4dolfinx

# import temperature and velocity data from openfoam
openfoam_reader = OpenFOAMReader("openfoam_data/heatexchanger.foam", cell_type=10)

print("Reading OpenFOAM data...")
T_cold = openfoam_reader.create_dolfinx_function(t=700, name="T", subdomain="cold")
T_hot = openfoam_reader.create_dolfinx_function(t=700, name="T", subdomain="hot")
u_cold = openfoam_reader.create_dolfinx_function(t=700, name="U", subdomain="cold")
u_hot = openfoam_reader.create_dolfinx_function(t=700, name="U", subdomain="hot")

# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/temperature_hot.bp",
#     T_hot,
#     "BP5",
# )
# writer.write(t=0)
# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/temperature_cold.bp",
#     T_cold,
#     "BP5",
# )
# writer.write(t=0)
# quit()

# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/u_cold.bp",
#     u_cold,
#     "BP5",
# )
# writer.write(t=0)
# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/u_hot.bp",
#     u_hot,
#     "BP5",
# )
# writer.write(t=0)

print("Building FESTIM model...")
# FESTIM model
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

# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/mesh.bp",
#     my_mesh.mesh,
#     "BP5",
# )
# writer.write(t=0)
# quit()

H = F.Species(name="H", mobile=True)
my_model.species = [H]

mat_1 = F.Material(name="mat_1", D_0=1e-06, E_D=0)

vol_hot = F.VolumeSubdomain(id=id_hot, material=mat_1)
vol_cold = F.VolumeSubdomain(id=id_cold, material=mat_1)
hot_inlet = F.SurfaceSubdomain(id=if_hot_inlet)
cold_inlet = F.SurfaceSubdomain(id=if_cold_inlet)

my_model.subdomains = [vol_hot, vol_cold, hot_inlet, cold_inlet]


# interpolate temperature fields onto one function
print("Interpolating temperature field...")
T_ele = basix.ufl.element(
    basix.ElementFamily.P,
    my_mesh.mesh.basix_cell(),
    1,
    basix.LagrangeVariant.equispaced,
)
V = fem.functionspace(my_mesh.mesh, T_ele)
my_temperature_field = fem.Function(V)


mesh_T = adios4dolfinx.read_mesh("results/temperature_field_cp.bp", MPI.COMM_WORLD)
V_T = fem.functionspace(mesh_T, T_ele)
T = fem.Function(V_T)
T.name = "T"
adios4dolfinx.read_function("results/temperature_field_cp.bp", T, time=0.0)
nmm_interpolate(my_temperature_field, T)
quit()

nmm_interpolate(
    my_temperature_field, T_hot, cells=volume_meshtags.find(id_hot), padding=1e-4
)
# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/T_with_hot.bp",
#     my_temperature_field,
#     "BP5",
# )
# writer.write(t=0)
nmm_interpolate(
    my_temperature_field, T_cold, cells=volume_meshtags.find(id_cold), padding=1e-4
)
# writer = VTXWriter(
#     MPI.COMM_WORLD,
#     "results/T_with_hot_and_cold.bp",
#     my_temperature_field,
#     "BP5",
# )
# writer.write(t=0)

# # convert to Kelvin
# my_temperature_field.x.array[:] += 273.15

my_model.temperature = my_temperature_field
print("Temperature field interpolated.")

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
# my_model.settings = F.Settings(
#     atol=1e-08, rtol=1e-10, max_iterations=30, transient=False
# )

mobile = F.VTXSpeciesExport(filename="results/mobile_concentration.bp", field=H)
my_model.exports = [mobile]


from dolfinx import log

log.set_log_level(log.LogLevel.INFO)

my_model.initialise()
my_model.run()
