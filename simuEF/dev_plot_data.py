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
import sys
from pathlib import Path
import fedoo as fd
from tools_fea import read_props

sys.path.append(str(Path(__file__).resolve().parent.parent))
from lstm.tools_database import read_data
import numpy as np

read_data("simuEF/train_fea.csv")
