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
import time

sys.path.append(str(Path(__file__).resolve().parent.parent))
from tools_homogeneisation import *
from tools_database import generate_data_csv, read_data

# import torch
read_data("lstm/dataset/test_fea.csv", i=0)
# read_data("lstm/dataset/test_dataset.csv", i=1)
# read_data("lstm/dataset/predictions.csv", i=1)
