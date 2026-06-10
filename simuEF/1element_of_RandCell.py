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

import fedoo as fd
import numpy as np
import matplotlib.pyplot as plt
from tools_fea import plot_data_6D, read_props

element = 10000
dataset = fd.read_data("simuEF/fea_RandCell.fdz")
n_iter = dataset.n_iter
print("n_iter =", n_iter)
time = np.linspace(0, 1, n_iter + 1)
stress_array = np.zeros((6, n_iter + 1))
strain_array = np.zeros((6, n_iter + 1))
mean_strain_array = np.zeros((6, n_iter + 1))
mean_stress_array = np.zeros((6, n_iter + 1))
all_components = ["XX", "YY", "ZZ", "XY", "XZ", "YZ"]
time = np.linspace(0, 1, n_iter + 1)
for k in range(dataset.n_iter):
    dataset.load(k)
    for i, component in enumerate(all_components, start=0):
        data_stress = dataset.get_data(
            field="Stress", component=component, data_type="GaussPoint"
        )
        data_strain = dataset.get_data(
            field="Strain", component=component, data_type="GaussPoint"
        )
        stress_array[i, k + 1] = data_stress[element]
        strain_array[i, k + 1] = data_strain[element]

plot_data_6D(stress_array, strain_array, time)
# print(stress_array)
# print(strain_array)

# data_stress = dataset.get_data(
#     field="Fext(MeanStrain)", component=component, data_type="GaussPoint"
# )
# data_MeanStrain = dataset.get_data(
#     field="MeanStrain", component=component, data_type="GaussPoint"
# )
# mean_stress_array[i, k + 1] = data_stress[0]
# mean_strain_array[i, k + 1] = data_MeanStrain[0]
