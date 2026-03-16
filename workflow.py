from pathlib import Path
import numpy as np
from tools_homogeneisation import *

import os
import fedoo as fd

# from fedoo.core.boundary_conditions import ListBC, BoundaryCondition

# import Pat

from simuEF.tools_fea import *

np.float_ = np.float64
# cell = "RhombicDodecahedron40"
# meshfile = f"cellules/{cell}.vtk"

material_law = "SMAUT"

props = read_props("simuEF/params_sma_init.txt")

typesim_to_loads = define_typesim_to_loads(0.05)
cell = "Cuboctahedron40"
for typesim in typesim_to_loads.keys():
    load = typesim_to_loads.get(typesim)
    cell_fea(props, material_law, typesim, load, cell)
    process_data_fea(typesim, cell)
    erase_fea_file(typesim)
# run_homogeneisation(cell=cell)


# cell = "Cuboctahedron40"
# density_to_load = [30, 40, 50, 60]
# np.float_ = np.float64
# for density in density_to_load:
#     print(f"gyroid{density}")
#     cell = f"gyroid{density}"
#     for typesim in typesim_to_loads.keys():
#         load = typesim_to_loads.get(typesim)
#         cell_fea(props, material_law, typesim, load, cell)
#         process_data_fea(typesim, cell)
#         erase_fea_file(typesim)
#     run_homogeneisation(cell=cell)
