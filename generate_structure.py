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
    BodyCenteredCubic,
    # Cubic,
    # Cuboctahedron,
    # Diamond,
    # FaceCenteredCubic,
    TruncatedCuboctahedron,
    Octahedron,
    OctetTruss,
    TruncatedCube,
    # RhombicCuboctahedron,
    # RhombicDodecahedron,
    TruncatedOctahedron,
]


# shape = FaceCenteredCubic
# shape_name = shape.__name__


def generate_structure(shape_):
    lattice_shape = shape(density=0.4)
    shapeStep = lattice_shape.generate()
    cq.exporters.export(shapeStep, f"{shape_}.step")

    print("strut_radius=", lattice_shape.strut_radius)
    phases_cut = [Phase(shapeStep)]
    meshPeriodic(
        mesh_file=f"{shape_}.step",
        rve=rve,
        listPhases=phases_cut,
        order=1,
        size=0.05,
        output_file=f"simuEF/cellules/{shape}40.vtk",
    )

    # mesh = pv.read(f"simuEF/cellules/{shape}40.vtk")
    # mesh.plot(show_edges=True)

    os.remove(f"{shape_}.step")


for shape in lattice_shapes_name:
    shape_name = shape.__name__
    print(shape_name)
    generate_structure(shape_name)
