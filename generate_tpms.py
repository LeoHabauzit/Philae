from microgen import Tpms
from microgen.shape.surface_functions import gyroid
import cadquery as cq
from microgen import Phase, Rve, meshPeriodic
import pyvista as pv

rve = Rve(dim=1.0)
density_to_loads = [0.8, 0.9]
for density in density_to_loads:
    geometry = Tpms(
        surface_function=gyroid,
        density=density,
        resolution=30,
    )
    shape = geometry.generate(type_part="sheet")
    cq.exporters.export(shape, f"gyroid{density * 100}.step")
    phases_cut = [Phase(shape)]
    meshPeriodic(
        mesh_file=f"gyroid{density * 100}.step",
        rve=rve,
        listPhases=phases_cut,
        order=1,
        size=0.45,
        output_file=f"simuEF/cellules/gyroid{density * 100}.vtk",
    )

    # mesh = pv.read(f"gyroid{density}.vtk")
    # mesh.plot(show_edges=True)
