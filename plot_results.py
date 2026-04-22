from tools_homogeneisation import *

cell = "RhombicCuboctahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)
# xi_values = np.arange(0.05, 0.15, 0.02)
xi_values = [0.35, 0.4, 0.45]
fig, axes_iso = plt.subplots(1, 2, figsize=(12, 5))
plot_isosurface_strut_material(finalprops, xi_values, cell=cell, axes=axes_iso)
# # plot_isosurface_strut_material(
# #     final_props, xi_modif=0.01, cell=cell, axes=axes_iso, i=0
# # )

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_strain_loads(finalprops, cell=cell, axs=axes_strain)

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_xi_stress(finalprops, cell=cell, axs=axes_strain)
plt.show()
