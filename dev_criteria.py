from plot_criteria import *
from tools_homogeneisation import *
import matplotlib.pyplot as plt

cell = "RhombicCuboctahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)
fig, ax = plt.subplots()
plot_drucker_ani_radius(ax, finalprops, 0.4, 300)
print(radius_for_drucker_ani(finalprops, 0.4, 300, theta=np.pi))
plt.show()
