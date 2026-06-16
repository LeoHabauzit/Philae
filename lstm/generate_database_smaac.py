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

# import matplotlib as mpl
t1 = time.time()
cell = "RhombicDodecahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)

# props = read_props("simuEF/params_sma_init.txt")
generate_data_csv(
    finalprops, 3000, csv_file="lstm/dataset/validation_dataset.csv", stress_lim=150
)
# read_data("lstm/train_dataset.csv", i=1)
t2 = time.time()
print(f"Time taken: {t2 - t1:.2f} seconds")
