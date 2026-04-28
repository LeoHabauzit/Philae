from tools_homogeneisation import *
# import matplotlib as mpl

cell = "RhombicCuboctahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)
# xi_values = np.arange(0.05, 0.15, 0.02)
xi_values = []
fig, axes_iso = plt.subplots(1, 2, figsize=(12, 5))
plot_xi_gradient(finalprops, fig, nx=100, ny=100, axes=axes_iso)
plot_isosurface_strut_material(finalprops, xi_values, cell=cell, axes=axes_iso)
handles, labels = axes_iso[0].get_legend_handles_labels()
fig.legend(
    dict(zip(labels, handles)).values(),
    dict(zip(labels, handles)).keys(),
    fontsize=8,
    loc="lower center",
)

# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_stress_strain_loads(finalprops, cell=cell, axs=axes_strain)

# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_xi_stress(finalprops, cell=cell, axs=axes_strain)
plt.show()
