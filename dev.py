from pathlib import Path
import numpy as np
from tools_homogeneisation import *
from Umat.loi_SMAUT_props import umat_smaut
from simuEF.tools_fea import process_data_fea, mises_strain_fea
import matplotlib.pyplot as plt

# cell = "Cuboctahedron40"
# process_data_fea("compression", cell)
# props = vect_props_smaut_test(
#     load_variable_props(f"results_params/params_smaut_{cell}.txt")
# )
# typesim = "tension"
# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_stress_strain_loads(props, cell=cell, axs=axes_strain)
# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_stress_mises_strain_loads(props, cell, axes_strain)
# plt.show()
cell = "Cuboctahedron40"
props_smaut = vect_props_smaut_test(
    load_variable_props(f"results_params/params_smaut_{cell}.txt")
)
xi_modif = right_artificial_xi(props_smaut, cell)

props_smaut = load_variable_props(f"results_params/params_smaut_{cell}.txt")
props_smani = load_variable_props(f"results_params/params_smani_{cell}.txt")
final_props = vect_props_smani(props_smani, props_smaut)
fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_xi_stress(final_props, cell=cell, axs=axes_strain)

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_strain_loads(final_props, cell=cell, axs=axes_strain)
plt.show()
