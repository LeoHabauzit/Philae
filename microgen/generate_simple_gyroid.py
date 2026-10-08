import cadquery as cq
import numpy as np
import pyvista as pv
from cadquery.vis import show
from microgen import Tpms
from microgen.shape import strut_lattice
from microgen.shape.surface_functions import gyroid

box = cq.Workplane("XY").box(5, 5, 10)
# def gradeing_function(x, y, z):
#     return gyroid(x, y, z) + 0.1 * (x**2 + y**2 + z**2)
geometry = Tpms(
    surface_function=gyroid,
    cell_size=1.0,
    density=0.25,
    # thickness=0.1,
    repeat_cell=(1, 1, 1),
    resolution=30,
)

geometry.sheet.plot(color="white")
print(geometry.offset)
sheet = geometry.sheet
# print(sheet.offset)
plotter = pv.Plotter()

plotter.add_mesh(
    sheet,
    color="white",
    smooth_shading=False,
)
plotter.add_measurement_widget()
plotter.show()
