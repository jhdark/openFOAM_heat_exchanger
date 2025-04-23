import festim as F
from foam2dolfinx import OpenFOAMReader
from dolfinx import fem
from festim.helpers import nmm_interpolate
import basix
from dolfinx.io import VTXWriter
from mpi4py import MPI
import adios4dolfinx


# import temperature and velocity data from openfoam
openfoam_reader = OpenFOAMReader("openfoam_data/heatexchanger.foam", cell_type=10)

print("Reading OpenFOAM data...")
T_cold = openfoam_reader.create_dolfinx_function(t=700, name="T", subdomain="cold")
T_hot = openfoam_reader.create_dolfinx_function(t=700, name="T", subdomain="hot")

my_mesh = F.MeshFromXDMF(
    volume_file="mesh_files/mesh_domains.xdmf",
    facet_file="mesh_files/mesh_boundaries.xdmf",
)
volume_meshtags = my_mesh.define_volume_meshtags()
id_hot = 6
id_cold = 7

# interpolate temperature fields onto one function
print("Interpolating temperature fields...")
T_ele = basix.ufl.element(
    basix.ElementFamily.P,
    my_mesh.mesh.basix_cell(),
    5,
    basix.LagrangeVariant.equispaced,
)
V = fem.functionspace(my_mesh.mesh, T_ele)
my_temperature_field = fem.Function(V)

nmm_interpolate(my_temperature_field, T_hot, cells=volume_meshtags.find(id_hot))
nmm_interpolate(my_temperature_field, T_cold, cells=volume_meshtags.find(id_cold))
writer = VTXWriter(
    MPI.COMM_WORLD,
    "results/temperature_field.bp",
    my_temperature_field,
    "BP5",
)
writer.write(t=0)

print("Writing field to file...")

my_temperature_field.name = "T"
adios4dolfinx.write_mesh("results/temperature_field_cp.bp", mesh=my_mesh.mesh)
adios4dolfinx.write_function(
    "results/temperature_field_cp.bp", my_temperature_field, time=0.0
)
