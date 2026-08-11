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
from plot_criteria import *
from tools_homogeneisation import *

cell = "Cuboctahedron40"
typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}


props_var_smadi = load_variable_props("results_params/params_smadi_bulk.txt")
finalprops_smadi = vect_props_smadi(props_var_smadi)
# xi_values = np.arange(0.05, 0.15, 0.02)
xi_values = []
# fig, axes_iso = plt.subplots(1, 2, figsize=(12, 5))
# plot_xi_gradient(finalprops, fig, nx=100, ny=100, axes=axes_iso)
# plot_isosurface_strut_material(finalprops, xi_values, cell=cell, axes=axes_iso)
# handles, labels = axes_iso[0].get_legend_handles_labels()
# fig.legend(
#     dict(zip(labels, handles)).values(),
#     dict(zip(labels, handles)).keys(),
#     fontsize=8,
#     loc="lower center",
# )

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_stress_strain_loads(
#     finalprops_smaac, UMAT="SMAAC", cell=cell, axs=axes_strain, color="orange"
# )
# print(
#     "erreur_finale=",
#     calc_cost_smaac(
#         props_var_smaac,
#         list_typesim=typesim_to_loads,
#         cell=cell,
#         props_cubic=props_cubic,
#     ),
# )
plot_stress_strain_loads(
    finalprops_smadi, UMAT="SMADI", cell=cell, axs=axes_strain, color="red"
)
# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))

print(
    "erreur_finale=",
    calc_cost_smadi(props_var_smadi, list_typesim=typesim_to_loads, cell=cell),
)
fig, axes_iso = plt.subplots(1, 1, figsize=(10, 10))
# plot_drucker_contour(
#     axes_iso,
#     props=finalprops_smadi,
#     xi=0.1,
#     T=300,
#     plane="s11-s22",
#     npts=400,
# )


axes_iso.grid(True)
axes_iso.spines["left"].set_position("zero")
axes_iso.spines["right"].set_color("none")
axes_iso.spines["bottom"].set_position("zero")
axes_iso.spines["top"].set_color("none")
axes_iso.set_aspect("equal", adjustable="box")
axes_iso.set_xlabel("S11 [MPa]", loc="right")
axes_iso.set_ylabel("S22 [MPa]", loc="top")

fig, axes_iso = plt.subplots(1, 1, figsize=(10, 10))
# plot_drucker_contour(
#     axes_iso,
#     props=finalprops_smadi,
#     xi=0.1,
#     T=300,
#     plane="s11-s12",
#     npts=400,
# )


axes_iso.grid(True)
axes_iso.spines["left"].set_position("zero")
axes_iso.spines["right"].set_color("none")
axes_iso.spines["bottom"].set_position("zero")
axes_iso.spines["top"].set_color("none")
axes_iso.set_aspect("equal", adjustable="box")
axes_iso.set_xlabel("S11 [MPa]", loc="right")
axes_iso.set_ylabel("S12 [MPa]", loc="top")
plt.show()
