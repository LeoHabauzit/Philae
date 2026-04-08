from pathlib import Path
import numpy as np
from tools_homogeneisation import *
from Umat.loi_sma import umat_sma

from simuEF.tools_fea import run_linear_homogenization
import matplotlib.pyplot as plt

cell = "RhombicCuboctahedron40"
typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
# E = props_var[0]
# C = props_var[1]
Hmax = props_var[0]
sigmacrit = props_var[1]
# dT = props_var[4]
# sigmacaliber = props_var[5]
b_prager = props_var[2]
n_prager = props_var[3]
Mf0_test = props_var[4]
Dsf = props_var[5]

props_cubic = run_linear_homogenization(f"{cell}")
E = props_cubic[0]
nu = props_cubic[1]
G = props_cubic[2]
# nu = 0.42
# G = E / (2 * (1 + nu))
print(E, nu, G)
alpha = 1.0e-4
E_A = E
E_M = E

G_A = G
G_M = G
nu_A = nu
nu_M = nu
alphaA = alpha
alphaM = alpha
flagT = 0.0
Hmin = 0.00
k1 = 0.008
C_A = 5.36
C_M = 5.36
Ms0 = 270
Mf0 = 250
As0 = 280
Af0 = 300
n1 = 0.05
n2 = 0.05
n3 = 0.05
n4 = 0.05
sigmacaliber = 300.0
b_prager = b_prager
n_prager = n_prager
c_lambda = 1.0e-6
p0_lambda = 1.0e-3
n_lambda = 1.0
alpha_lambda = 1.0e8
F = 0.5
L = 1.5
K = 0
F_dfa = F
G_dfa = F_dfa
H_dfa = F_dfa
L_dfa = L
M_dfa = L_dfa
N_dfa = L_dfa
K_dfa = K

props_di = np.array(
    [
        flagT,
        E_A,
        E_M,
        nu_A,
        nu_M,
        alphaA,
        alphaM,
        Hmin,
        Hmax,
        k1,
        sigmacrit,
        C_A,
        C_M,
        Ms0,
        Mf0,
        As0,
        Af0,
        n1,
        n2,
        n3,
        n4,
        sigmacaliber,
        b_prager,
        n_prager,
        c_lambda,
        p0_lambda,
        n_lambda,
        alpha_lambda,
    ]
)

# Mf0_test = Mf0
# Mf = 253
# Dsf = 20
# Ms0_test, As0_test, Af0_test = get_martensite_temp(Mf=Mf, Dsf=Dsf, T0=Mf + Dsf + 5)
# print("Ms=", Ms0_test, "As=", As0_test, "Af=", Af0_test)

props_test = vect_props_smaac(props_var, props_cubic)
print(calc_cost_strain(props_var, typesim_to_loads, cell, props_cubic))

fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
for i, typesim in enumerate(sorted(typesim_to_loads)):
    losses = []
    row = i // 3
    col = i % 3
    ax = axes_strain[row, col]
    umat_sma(props_test, typesim, "SMAAC")

    outputfile_global = f"Umat/results_SMAAC/results_{typesim}_global-0.txt"

    e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
        outputfile_global,
        usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
        unpack=True,
    )
    if typesim == "shear":
        ax.plot(e12, s12, label=f"{typesim} AC", color="orange")
        ax.legend()

    else:
        ax.plot(e11, s11, label=f"{typesim} AC", color="orange")
        ax.legend()

    umat_sma(props_di, typesim, "SMADI")

    outputfile_global = f"Umat/results_SMADI/results_{typesim}_global-0.txt"

    e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
        outputfile_global,
        usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
        unpack=True,
    )
    if typesim == "shear":
        ax.plot(e12, s12, label=f"{typesim} DI", color="red")
        ax.legend()

    else:
        ax.plot(e11, s11, label=f"{typesim} DI", color="red")
        ax.legend()
    data_simu_dir = f"simuEF/datas_simu/{cell}"
    results_dir = typesim
    if typesim == "shear":
        strain_num = e12
        stress_num = s12
        stress_exp = np.loadtxt(
            f"{data_simu_dir}/SXY/data_{results_dir}/Stress_{results_dir}.txt"
        )
        strain_exp = np.loadtxt(
            f"{data_simu_dir}/SXY/data_{results_dir}/MeanStrain_{results_dir}.txt"
        )
    else:
        stress_exp = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
        )
        strain_exp = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/MeanStrain_{results_dir}.txt"
        )
    ax.plot(strain_exp, stress_exp, color="blue", label="cellule")
    ax.legend()

# plt.plot(e11, s11)

props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)
# fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
# plot_xi_stress(finalprops, cell=cell, axs=axes_strain)
plt.show()
plt.close("all")
