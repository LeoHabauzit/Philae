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
from tools_homogeneisation import *

cell = "RhombicCuboctahedron40"
typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}


props_cubic = run_linear_homogenization(f"{cell}")
props_var_smaac = load_variable_props(f"results_params/params_smaac_{cell}.txt")
finalprops_smaac = vect_props_smaac(props_var_smaac, props_cubic)

props_var_smadi = load_variable_props(f"results_params/params_smadi_{cell}.txt")
finalprops_smadi = vect_props_smadi(props_var_smadi)
# xi_values = np.arange(0.05, 0.15, 0.02)
xi_values = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# fig, axes_iso = plt.subplots(1, 2, figsize=(12, 5))
# # plot_xi_gradient(finalprops_smaac, fig, nx=100, ny=100, axes=axes_iso)
# plot_isosurface_strut_material(finalprops_smaac, xi_values, cell=cell, axes=axes_iso)
# handles, labels = axes_iso[0].get_legend_handles_labels()
# axes_iso[0].set_xlim(-260, 260)
# axes_iso[0].set_ylim(-260, 260)
# axes_iso[1].set_xlim(-260, 260)
# axes_iso[1].set_ylim(-260, 260)

# fig.legend(
#     dict(zip(labels, handles)).values(),
#     dict(zip(labels, handles)).keys(),
#     fontsize=8,
#     loc="lower center",
# )

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_strain_loads(
    finalprops_smaac, UMAT="SMAAC", cell=cell, axs=axes_strain, color="orange"
)
print(
    "erreur_finale=",
    calc_cost_smaac(
        props_var_smaac,
        list_typesim=typesim_to_loads,
        cell=cell,
        props_cubic=props_cubic,
    ),
)
# plot_stress_strain_loads(
#     finalprops_smadi, UMAT="SMADI", cell=cell, axs=axes_strain, color="red"
# )
# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))

print(
    "erreur_finale=",
    calc_cost_smadi(props_var_smadi, list_typesim=typesim_to_loads, cell=cell),
)

# fig, axes_xi = plt.subplots(2, 3, figsize=(12, 8))
# plot_xi_stress(finalprops_smaac, cell=cell, axs=axes_xi)
# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_xi_stress(finalprops, cell=cell, axs=axes_strain)
plt.show()
