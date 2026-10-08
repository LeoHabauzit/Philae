"""
This script generates a gyroid-based TPMS (Triply Periodic Minimal Surface) structure
using *cadquery*, exports the geometry to a STEP file, performs periodic
meshing, and optionally remeshes the result for finite element method (FEM) simulations.

Steps performed:
1. Generate a TPMS geometry using the gyroid surface function.
2. Wrap the geometry into a microgen Phase object.
3. Export the geometry as a STEP file.
4. Import the STEP file and create a mesh with periodic constraints and export as VTK.
5. Optionally remesh the mesh while preserving periodicity.

Generated files:
- gyroid.step: CAD representation of the TPMS geometry.
- gyroid.vtk: Initial meshed structure.
- remeshed_gyroid_mesh.vtk (optional): Quality-improved periodic mesh for FEM.

Dependencies:
- microgen (with the [cad] extra: cadquery-ocp-novtk via microgen.cad)
- pyvista
"""

from pathlib import Path

import pyvista as pv
from microgen import Phase, Rve, Tpms, mesh_periodic

# from microgen.phase import from_cad
from microgen.remesh import remesh_keeping_boundaries_for_fem
from microgen.shape.surface_functions import gyroid

# 1. Generate a TPMS geometry using the gyroid surface function.
geometry = Tpms(
    surface_function=gyroid,
    density=0.25,
    resolution=40,
)

# 2. Wrap the geometry into a microgen Phase object.
shape = geometry.generate_cad()
phases = []
phases.append(Phase.from_cad(shape))
rve = Rve(dim=1.0)

# # 3. Export the geometry as a STEP file
step_file = "gyroid.step"

# # 4. Import the STEP file and create a mesh with periodic constraints and export as VTK.

vtk_file = "gyroid.vtk"

mesh_periodic(
    mesh_file=step_file,
    rve=rve,
    list_phases=phases,
    order=1,
    size=0.025,
    output_file=vtk_file,
    tol=1e-8,
)

# # 5. Optionally remesh the mesh while preserving periodicity.

initial_gyroid = pv.UnstructuredGrid(vtk_file)
max_element_edge_length = 0.02
remeshed_gyroid = remesh_keeping_boundaries_for_fem(
    initial_gyroid,
    periodic=True,
    hmax=max_element_edge_length,
)
remeshed_vtk_file = "remeshed_gyroid_mesh.vtk"
remeshed_gyroid.save(remeshed_vtk_file)
