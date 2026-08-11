import cadquery as cq
import pyvista as pv
from microgen import Phase, Rve, Tpms, meshPeriodic
from microgen.shape.surface_functions import gyroid

rve = Rve(dim=1.0)
density_to_loads = [0.6]
# for density in density_to_loads:
#     geometry = Tpms(
#         surface_function=gyroid,
#         density=density,
#         resolution=30,
#     )
#     shape = geometry.generate(type_part="sheet")
#     cq.exporters.export(shape, f"gyroid{density * 100}.step")
#     phases_cut = [Phase(shape)]
#     meshPeriodic(
#         mesh_file=f"gyroid{density * 100}.step",
#         rve=rve,
#         listPhases=phases_cut,
#         order=1,
#         size=0.45,
#         output_file=f"simuEF/cellules/gyroid{density * 100}.vtk",
#     )

mesh = pv.read("simuEF/cellules/RhombicDodecahedron40.vtk")
mesh.plot(show_edges=False, color="white", show_scalar_bar=False, cpos="xy")
