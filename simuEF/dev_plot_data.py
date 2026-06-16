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
from lstm.tools_database import read_data, read_multiple_data
import numpy as np

# read_data("simuEF/train_fea_50_cases.csv", i=0)
# read_data("simuEF/csv_files/compare_fea_sym2.csv", i=0)
files = {
    "simuEF/csv_files/ref.csv": {
        "color": "black",
        "linestyle": "-",
        "label": "Référence",
    },
    "simuEF/csv_files/sym1.csv": {
        "color": "red",
        "linestyle": "--",
        "label": "sym1",
    },
    # "simuEF/csv_files/compare_fea_sym2.csv": {
    #     "color": "blue",
    #     "linestyle": ":",
    #     "label": "sym2",
    # },
}
read_multiple_data(files)
