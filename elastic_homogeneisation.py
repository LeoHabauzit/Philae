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

run_linear_homogenization(
    cell="PeriodicGyroid25", young_modulus=67538, poisson_ratio=0.42
)
