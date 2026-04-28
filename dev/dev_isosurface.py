import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import matplotlib.lines as mlines
from continum_mech import *
from criteria import *
from matplotlib.colors import LinearSegmentedColormap

# grille x, y de -10 à 10
from tools_homogeneisation import *

cell = "RhombicCuboctahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)


def xi_for_drucker_ani(props, r, T, theta, plane="s11-s22"):
    """Rayon r tel que σ_drucker = sigma_y pour un angle theta"""
    # print(theta)

    def f(xi):
        v = stress_vector_from_polar(r, theta, plane)
        return get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=1)

    a = 1e-9
    b = 1.0
    fa = f(a)
    fb = f(b)

    # max_expand = 20
    # i = 0
    # while fa * fb > 0 and i < max_expand and a > 1e-10:
    #     a *= 0.5
    #     fa = f(a)

    # while fa * fb > 0 and i < max_expand and b < 1e8:
    #     b *= 2
    #     fb = f(b)

    # if fa * fb > 0:
    #     return fa * fb * 10
    # else:
    sol = root_scalar(f, bracket=[a, b], method="brentq")
    return sol.root


fig, axes_iso = plt.subplots(1, 2, figsize=(20, 15))
fig.suptitle("Carte de xi critique", fontsize=14)
# xi_values = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
xi_values = []

plot_xi_gradient(finalprops, axes_iso)
plot_isosurface_strut_material(finalprops, xi_values, cell=cell, axes=axes_iso)
# axes_iso[0].set_title("Plan \sigma_{11}-s22")
# axes_iso[1].set_title("Plan s11-s12")

plt.show()
# print(get_Phi_forward_SMA(props=finalprops, v=v, xi=xi, T=T, ani=1))
# print(radius_for_drucker_ani(finalprops, r, T, theta, plane="s11-s22"))
# r_found, theta_found = get_polar_coords(100, -100)
# print(radius_for_drucker_ani(finalprops, r_found, T, theta_found, plane="s11-s22"))
