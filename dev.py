from pathlib import Path
import numpy as np
from tools_homogeneisation import *
from Umat.loi_sma import umat_sma

from simuEF.tools_fea import process_data_fea, mises_strain_fea
import matplotlib.pyplot as plt

cell = "Cuboctahedron40"
typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}
props_var = load_variable_props(f"results_params/params_smadi_{cell}.txt")
E = props_var[0]
C = props_var[1]
Hmax = props_var[2]
sigmacrit = props_var[3]
dT = props_var[4]
sigmacaliber = props_var[5]
b_prager = props_var[6]
n_prager = props_var[7]

nu = 0.42
alpha = 1.0e-4
E_A = E
E_M = E
G_A = 25000
G_M = 25000
nu_A = nu
nu_M = nu
alphaA = alpha
alphaM = alpha
flagT = 0.0
Hmin = 0.00
k1 = 0.008
C_A = C
C_M = C
Ms0 = 273
Mf0 = 253
As0 = 283
Af0 = 303
n1 = 0.05
n2 = 0.05
n3 = 0.05
n4 = 0.05
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

Mf0_test = Mf0
Ms0_test, As0_test, Af0_test = get_martensite_temp(Mf=253, Dsf=20, T0=278)
print("Ms=", Ms0_test, "As=", As0_test, "Af=", Af0_test)
props_ac = np.array(
    [
        flagT,
        E_A,
        E_M,
        nu_A,
        nu_M,
        G_A,
        G_M,
        alphaA,
        alphaM,
        Hmin,
        Hmax,
        k1,
        sigmacrit,
        C_A,
        C_M,
        Ms0_test,
        Mf0_test,
        As0_test,
        Af0_test,
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
        F_dfa,
        G_dfa,
        H_dfa,
        L_dfa,
        M_dfa,
        N_dfa,
        K_dfa,
    ]
)


fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
for i, typesim in enumerate(sorted(typesim_to_loads)):
    losses = []
    row = i // 3
    col = i % 3
    ax = axes_strain[row, col]
    umat_sma(props_ac, typesim, "SMAAC")

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
# plt.plot(e11, s11)
plt.show()
