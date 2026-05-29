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
generate_data_csv(finalprops, 1000)
# read_data("lstm/train_dataset.csv", i=1)
t2 = time.time()
print(f"Time taken: {t2 - t1:.2f} seconds")
