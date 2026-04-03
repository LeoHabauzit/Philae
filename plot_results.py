from tools_homogeneisation import *

cell = "Cuboctahedron40"
props_smadi = vect_props_smadi_test(
    load_variable_props(f"results_params/params_smadi_{cell}.txt")
)
xi_modif = right_artificial_xi(props_smadi, cell)

props_smadi = load_variable_props(f"results_params/params_smadi_{cell}.txt")
props_smani = load_variable_props(f"results_params/params_smani_{cell}.txt")
final_props = vect_props_smani(props_smani, props_smadi)
final_props2 = vect_props_smadi_test(props_smadi)

fig, axes_iso = plt.subplots(1, 2, figsize=(12, 5))
# plot_isosurface_strut_material(final_props, xi_modif=xi_modif, cell=cell, axes=axes_iso)
plot_isosurface_strut_material(
    final_props, xi_modif=xi_modif, cell=cell, axes=axes_iso, i=0
)

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
plot_stress_strain_loads(final_props, cell=cell, axs=axes_strain)
plt.show()
