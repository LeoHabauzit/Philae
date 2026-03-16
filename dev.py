from pathlib import Path
import numpy as np
from tools_homogeneisation import *
from Umat.loi_SMAUT_props import umat_smaut
from simuEF.tools_fea import process_data_fea, mises_strain_fea
import matplotlib.pyplot as plt

cell = "Cuboctahedron40"
props = vect_props_smaut_test(
    load_variable_props(f"results_params/params_smaut_{cell}.txt")
)
typesim = "tension"

# umat_smaut(props, typesim)
# outputfile_global = f"Umat/results_smaut/results_{typesim}_global-0.txt"

# (
#     e11,
#     e22,
#     e33,
#     e12,
#     e13,
#     e23,
#     s11,
#     s22,
#     s33,
#     s12,
#     s13,
#     s23,
#     xi,
#     et11,
#     et22,
#     et33,
#     et12,
#     et13,
#     et23,
# ) = np.loadtxt(
#     outputfile_global,
#     usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 25, 26, 27, 28, 29, 30, 31),
#     unpack=True,
# )
# E = np.vstack([et11, et22, et33, et12, et13, et23]).T
# et = np.array([mises_strain(v) for v in E])
# typesim_to_loads = {
#     # "tension",
#     "biaxial_tension",
#     "compression",
#     "biaxial_compression",
#     "tencomp",
#     "shear",
# }
# for typesim in typesim_to_loads:
#     process_data_fea(typesim, "Cuboctahedron40")

# stress_exp = np.loadtxt(
#     f"datas_simu/Cuboctahedron40/SXX/data_tension/Stress_tension.txt"
# )
# strain_mises_exp = np.loadtxt(
#     f"datas_simu/Cuboctahedron40/SXX/data_tension/Transformation_strain_tension.txt"
# )
# strain_labels = [
#     "et11",
#     "et22",
#     "et33",
#     "et12",
#     "et13",
#     "et23",
# ]
# strain_mises = np.zeros(6)
# for k, label in enumerate(strain_labels):
#     plt.figure()
#     et = np.loadtxt(f"{label}.txt")
#     strain_mises[k] = et[-1]
#     plt.plot(et, stress_exp, label=f"{label}")

# print(mises_strain_fea(strain_mises))
# plt.figure()
# plt.plot(strain_mises_exp, stress_exp)
# # plt.plot(et, s11)
# plt.legend()
# plt.show()

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_strain_loads(props, cell=cell, axs=axes_strain)
fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_mises_strain_loads(props, cell, axes_strain)
plt.show()
