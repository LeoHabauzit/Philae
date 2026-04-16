from tools_homogeneisation import *

cell = "Cuboctahedron40"
# props_smadi = vect_props_smadi_test(
#     load_variable_props(f"results_params/params_smadi_{cell}.txt")
# )
# xi_modif = right_artificial_xi(props_smadi, cell)

# props_smadi = load_variable_props(f"results_params/params_smadi_{cell}.txt")
# props_smani = load_variable_props(f"results_params/params_smani_{cell}.txt")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")

props_cubic = run_linear_homogenization(f"{cell}")
final_props = vect_props_smaac(props_var, props_cubic)

# fig, axes_iso = plt.subplots(1, 2, figsize=(12, 5))
# # plot_isosurface_strut_material(final_props, xi_modif=xi_modif, cell=cell, axes=axes_iso)
# plot_isosurface_strut_material(
#     final_props, xi_modif=0.01, cell=cell, axes=axes_iso, i=0
# )

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_strain_loads(final_props, cell=cell, axs=axes_strain)

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_xi_stress(final_props, cell=cell, axs=axes_strain)
plt.show()
