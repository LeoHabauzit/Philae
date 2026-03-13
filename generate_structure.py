##vesrion de microgen a installé en local avec un pip install -e .
from microgen import (
    BodyCenteredCubic,
    Cubic,
    Cuboctahedron,
    Diamond,
    FaceCenteredCubic,
    Octahedron,
    OctetTruss,
    RhombicCuboctahedron,
    RhombicDodecahedron,
    TruncatedCube,
    TruncatedCuboctahedron,
    TruncatedOctahedron,
)
from microgen import rve
import cadquery as cq
from microgen import Phase, Rve, meshPeriodic
import os
import pyvista as pv

rve = Rve(dim=1.0)
lattice_shapes_name = [
    "BodyCenteredCubic",
    "Cubic",
    "Cuboctahedron",
    "Diamond",
    "FaceCenteredCubic",
    "TruncatedCuboctahedron",
    "Octahedron",
    "OctetTruss",
    "TruncatedCube",
    "RhombicCuboctahedron",
    "RhombicDodecahedron",
    "TruncatedOctahedron",
]


shape = RhombicCuboctahedron
shape_name = shape.__name__


lattice_shape = shape(density=0.4)
shapeStep = lattice_shape.generate()
cq.exporters.export(shapeStep, f"{shape_name}.step")

print("strut_radius=", lattice_shape.strut_radius)
phases_cut = [Phase(shapeStep)]
meshPeriodic(
    mesh_file=f"{shape_name}.step",
    rve=rve,
    listPhases=phases_cut,
    order=1,
    size=0.05,
    output_file=f"{shape_name}.vtk",
)

mesh = pv.read(f"{shape_name}.vtk")
mesh.plot(show_edges=True)
