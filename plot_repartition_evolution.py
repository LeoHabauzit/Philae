import os

os.environ["OMP_NUM_THREADS"] = "1"  # avant tout import de simcoon/fedoo
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

try:
    from scipy.sparse.linalg._dsolve.linsolve import useUmfpack as _scipy_uu

    if not hasattr(_scipy_uu, "u"):
        _scipy_uu.u = True  # active umfpack, évite le crash dans base.py
except Exception:
    pass

from pathlib import Path

import fedoo as fd
import numpy as np

# from fedoo.core.boundary_conditions import ListBC, BoundaryCondition
# import Pat
from simuEF.tools_fea import *
from tools_homogeneisation import *

cell = "RhombicCuboctahedron40"
typesim = "tension"
material_law = "SMAUT"
props = read_props("simuEF/params_sma_init.txt")
typesim_to_loads = define_typesim_to_loads(0.05)
load = typesim_to_loads.get("tension")
# cell_fea(props, material_law, typesim, load, cell)
# process_data_fea(typesim, cell)
# erase_fea_file(typesim)


process_element_repartition(typesim)
