import numpy as np
import os, sys
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from plot_criteria import *
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, Bounds
from functools import partial
import pandas as pd
# parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# sys.path.append(parent_dir)
from Umat.loi_SMAUT_props import umat_smaut


def setup_ax(ax, xlabel, ylabel):
    ax.spines["left"].set_position("zero")
    ax.spines["bottom"].set_position("zero")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.text(xlim[1], 0, xlabel, fontsize=12, ha="right", va="bottom")
    ax.text(0, ylim[1], ylabel, fontsize=12, ha="left", va="top")


def load_variable_props(filepath):
    """
    Lit un fichier texte et retourne un vecteur numpy
    contenant toutes les valeurs numériques trouvées.

    Ignore automatiquement les lignes non convertibles en float.
    """

    values = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                values.append(float(line))
            except ValueError:
                # Ignore les lignes texte comme "traction-compression"
                continue

    return np.array(values)


def calc_cost_smaut(props_var, list_typesim, cell):
    """Calcule la fonction de coût pour un jeu de paramètres donné et un type de simulation,
    en comparant les résultats numériques et expérimentaux (erreur quadratique moyenne)

    Args:
        props_var (np.array): vecteur complet des propriétés matériaux
        typesim (_type_): _description_

    Returns:
        float : erreur liée a la simulation calculée avec la méthode des moindrs carrées
    """
    losses = []
    data_simu_dir = f"datas_simu/{cell}"
    for typesim in list_typesim:
        results_dir = typesim
        props = vect_props_smaut(props_var)
        umat_smaut(props, typesim)
        outputfile_global = f"Umat/results_smaut/results_{typesim}_global-0.txt"

        e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
            outputfile_global,
            usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
            unpack=True,
        )

        strain_num = e11
        stress_num = s11
        stress_exp = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
        )
        strain_exp = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/MeanStrain_{results_dir}.txt"
        )

        interp = interp1d(
            strain_exp, stress_exp, kind="linear", fill_value="extrapolate"
        )
        stress_exp_interp = interp(strain_num)

        losses.append(np.mean((stress_exp_interp - stress_num) ** 2))
    loss = np.mean(losses)
    return loss


def vect_props_smaut(props_var):
    """Transforme un vecteur de paramètres variables (E, Hmax, sigmacrit, C, dT, sigmacaliber)
    en un vecteur complet des propriétés attendu par la simulation UMAT.

    Args:
        props_var (_type_): vecteur props des propriétés a optimiser par l'évolution différentielle

    Returns:
        full_props: vecteur props complet utilisée dans la simulation SMAUT
    """

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
    nu_A = nu
    nu_M = nu
    alphaA = alpha
    alphaM = alpha
    flagT = 0.0
    Hmin = 0.00
    k1 = 0.008
    C_A = C
    C_M = C
    Ms0 = 250 + dT
    Mf0 = 230 + dT
    As0 = 240 + dT
    Af0 = 260 + dT
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

    full_props = np.array(
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
    return full_props


def vect_props_smaut_test(props_var):
    """Transforme un vecteur de paramètres variables (E, Hmax, sigmacrit, C, dT, sigmacaliber)
    en un vecteur complet des propriétés attendu par la simulation UMAT.

    Args:
        props_var (_type_): vecteur props des propriétés a optimiser par l'évolution différentielle

    Returns:
        full_props: vecteur props complet utilisée dans la simulation SMAUT
    """

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
    nu_A = nu
    nu_M = nu
    alphaA = alpha
    alphaM = alpha
    flagT = 0.0
    Hmin = 0.00
    k1 = 0.008
    C_A = C
    C_M = C
    Ms0 = 250 + dT
    Mf0 = 230 + dT
    As0 = 240 + dT
    Af0 = 260 + dT
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

    full_props = np.array(
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
            F_dfa,
            G_dfa,
            H_dfa,
            L_dfa,
            M_dfa,
            N_dfa,
            K_dfa,
        ]
    )
    return full_props


def find_first_stress_at_xi_limit(typesim, xi_lim, cell):
    data_simu_dir = f"datas_simu/{cell}"
    results_dir = typesim
    # xi_init = get_xi_init_sma(props, typesim)
    xi_init = 0
    if typesim == "shear":
        xi = np.loadtxt(f"{data_simu_dir}/SXY/data_{results_dir}/Xi_{results_dir}.txt")
    else:
        xi = np.loadtxt(f"{data_simu_dir}/SXX/data_{results_dir}/Xi_{results_dir}.txt")
    xi = xi + xi_init
    xi_lim = xi_lim + xi_init
    # print(xi_lim)

    idx = np.argmax(xi > xi_lim)  # premier True

    if xi[idx] <= xi_lim:
        raise ValueError(f"{typesim} xi ne dépasse jamais la limite")

    if typesim == "shear":
        s12 = np.loadtxt(
            f"{data_simu_dir}/SXY/data_{results_dir}/Stress_{results_dir}.txt"
        )
        s11 = np.zeros(np.shape(s12))
        s22 = np.zeros(np.shape(s12))
    elif typesim == "tension" or typesim == "compression":
        s11 = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
        )
        s22 = np.zeros(np.shape(s11))
        s12 = np.zeros(np.shape(s11))
    else:
        s11 = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
        )
        s22 = np.loadtxt(
            f"{data_simu_dir}/SYY/data_{results_dir}/Stress_{results_dir}.txt"
        )
        s12 = np.zeros(np.shape(s11))

    return s11[idx], s22[idx], s12[idx]


def vect_props_smani(props_var, props_smaut):
    """Transforme un vecteur de paramètres variables (E, Hmax, sigmacrit, C, dT, sigmacaliber)
    en un vecteur complet des propriétés attendu par la simulation UMAT.

    Args:
        props_var (_type_): vecteur props des propriétés a optimiser par l'évolution différentielle

    Returns:
        full_props: vecteur props complet utilisée dans la simulation SMAUT
    """

    b_prager = props_var[0]
    n_prager = props_var[1]
    F = props_var[2]
    L = props_var[3]
    K = props_var[4]

    E = props_smaut[0]
    C = props_smaut[1]
    Hmax = props_smaut[2]
    sigmacrit = props_smaut[3]
    dT = props_smaut[4]
    sigmacaliber = props_smaut[5]

    nu = 0.42
    alpha = 1.0e-4
    E_A = E
    E_M = E
    nu_A = nu
    nu_M = nu
    alphaA = alpha
    alphaM = alpha
    flagT = 0.0
    Hmin = 0.00
    k1 = 0.008
    C_A = C
    C_M = C
    Ms0 = 250 + dT
    Mf0 = 230 + dT
    As0 = 240 + dT
    Af0 = 260 + dT
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
    F_dfa = F
    G_dfa = F_dfa
    H_dfa = F_dfa
    L_dfa = L
    M_dfa = L_dfa
    N_dfa = L_dfa
    K_dfa = K

    full_props = np.array(
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
            F_dfa,
            G_dfa,
            H_dfa,
            L_dfa,
            M_dfa,
            N_dfa,
            K_dfa,
        ]
    )
    return full_props


def calc_cost_smani(props_var, list_typesim, xi_modif, cell):
    losses = []
    props_smaut = load_variable_props(f"results_params/params_smaut_{cell}.txt")
    props = vect_props_smani(props_var, props_smaut)
    # print(props_var)
    for typesim, theta in list_typesim.items():
        # props = vect_props_smani(props_var)
        s11_exp, s22_exp, s12_exp = find_first_stress_at_xi_limit(
            typesim, xi_lim=0.01, cell=cell
        )
        if typesim == "shear":
            r_umat = radius_for_drucker_ani(
                props=props, xi=xi_modif, T=300.0, theta=theta, plane="s11-s12"
            )
            r_exp = np.sqrt(s11_exp**2 + s12_exp**2)
        else:
            r_umat = radius_for_drucker_ani(
                props=props, xi=xi_modif, T=300.0, theta=theta, plane="s11-s22"
            )
            r_exp = np.sqrt(s11_exp**2 + s22_exp**2)

        if typesim == "compression" or typesim == "traction":
            losses.append(10 * ((r_exp - r_umat) ** 2))
        else:
            losses.append((r_exp - r_umat) ** 2)
    loss = np.mean(losses)
    return loss


def right_artificial_xi(props, cell):
    s11_exp, s22_exp, s12_exp = find_first_stress_at_xi_limit(
        "compression", xi_lim=0.01, cell=cell
    )
    r_exp = np.sqrt(s11_exp**2 + s22_exp**2)

    def f(xi_test):
        r_umat = radius_to_find_xi_lim(
            props=props,
            xi=xi_test,
            T=300.0,
            theta=0.0,
            plane="s11-s22",
        )

        return r_umat - r_exp

    a = 0.01
    b = 1
    # print(f(a), f(b))
    sol = root_scalar(f, bracket=[a, b], method="brentq")
    return sol.root


def plot_isosurface_strut_material(full_props, xi_modif, cell, axes, i=0):
    ax = axes[0]
    # -------------------------------------------Partie dans le plan (S11,S22)--------------------------------------
    typesim_to_loads = {
        "tension",
        "biaxial_tension",
        "compression",
        "biaxial_compression",
        "tencomp",
        # "shear",
    }
    if i == 1:
        plot_drucker_radius(ax, full_props, xi=xi_modif, T=300.0, plane="s11-s22")
    else:
        plot_drucker_ani_radius(ax, full_props, xi=xi_modif, T=300.0, plane="s11-s22")
    # plot_drucker_ani_radius(ax, props, xi=0.01, T=300.0, plane="s11-s22")

    X_points = []
    Y_points = []

    for typesim in typesim_to_loads:
        coords_typesim = find_first_stress_at_xi_limit(typesim, xi_lim=0.01, cell=cell)

        if typesim in ["tension", "compression", "tencomp"]:
            X_points.append(coords_typesim[0])
            Y_points.append(coords_typesim[1])

            X_points.append(coords_typesim[1])
            Y_points.append(coords_typesim[0])

        else:
            X_points.append(coords_typesim[0])
            Y_points.append(coords_typesim[1])

    ax.scatter(X_points, Y_points, label=f"{cell}")

    setup_ax(ax, r"$\sigma_{11} [MPa]$", r"$\sigma_{22} [MPa]$")
    ax.legend(fontsize=8)

    # -------------------------------------------Partie dans le plan (S11,S12)--------------------------------------

    typesim_to_loads = {
        "tension",
        # "biaxial_tension",
        "compression",
        # "biaxial_compression",
        # "tencomp",
        "shear",
    }
    ax = axes[1]
    # plot_drucker_radius(ax, full_props, xi=xi_modif, T=300.0, plane="s11-s12")
    if i == 1:
        plot_drucker_radius(ax, full_props, xi=xi_modif, T=300.0, plane="s11-s12")
    else:
        plot_drucker_ani_radius(ax, full_props, xi=xi_modif, T=300.0, plane="s11-s12")

    # for typemat in materials_to_load:
    X_points = []
    Y_points = []

    for typesim in typesim_to_loads:
        coords_typesim = find_first_stress_at_xi_limit(typesim, xi_lim=0.01, cell=cell)

        if typesim in ["shear"]:
            X_points.append(coords_typesim[0])
            Y_points.append(coords_typesim[2])

            X_points.append(coords_typesim[0])
            Y_points.append(-coords_typesim[2])

        else:
            X_points.append(coords_typesim[0])
            Y_points.append(coords_typesim[2])

    ax.scatter(X_points, Y_points, label=f"{cell}")

    setup_ax(ax, r"$\sigma_{11} [MPa]$", r"$\sigma_{12} [MPa]$")
    ax.legend(fontsize=8)


def plot_stress_strain_loads(full_props, cell, axs):
    data_simu_dir = f"datas_simu/{cell}"
    typesim_to_loads = {
        "tension",
        "biaxial_tension",
        "compression",
        "biaxial_compression",
        "tencomp",
        "shear",
    }

    for i, typesim in enumerate(sorted(typesim_to_loads)):
        losses = []
        row = i // 3
        col = i % 3

        ax = axs[row, col]
        results_dir = typesim
        umat_smaut(full_props, typesim)

        outputfile_global = f"Umat/results_smaut/results_{typesim}_global-0.txt"

        e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23, xi = np.loadtxt(
            outputfile_global,
            usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 25),
            unpack=True,
        )

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
            strain_num = e11
            stress_num = s11
            stress_exp = np.loadtxt(
                f"{data_simu_dir}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
            )
            strain_exp = np.loadtxt(
                f"{data_simu_dir}/SXX/data_{results_dir}/MeanStrain_{results_dir}.txt"
            )

        ax.plot(
            strain_exp,
            stress_exp,
            label=typesim,
        )
        ax.plot(strain_num, stress_num, c="red", label="UMAT SMA")
        ax.set_xlabel("E11 [%]")
        ax.set_ylabel("S11 [MPa]")
        ax.grid()
        ax.legend()
    plt.suptitle(f"{cell}")

    plt.tight_layout()

    # plt.title(f"Plot {typesim}")
    plt.legend(loc="upper left", fontsize=8)
    plt.grid(True)
    # plt.xlabel("E11[%]")
    # plt.ylabel("S11 [MPa]")


def evol_diff_smaut(bounds, cell, n_iter):
    all_params = {}
    typesim_to_loads = {
        "tension",
        # "biaxial_tension",
        "compression",
        "biaxial_compression",
        # "tencomp",
        # "shear",
    }

    loss = partial(calc_cost_smaut, list_typesim=typesim_to_loads, cell=cell)
    result = differential_evolution(
        func=loss,
        bounds=bounds,
        maxiter=n_iter,
        tol=1e-4,
        disp=True,
    )

    all_params["smaut"] = result.x
    df = pd.DataFrame(all_params)
    df.to_csv(
        f"results_params/params_smaut_{cell}.txt",
        index=False,
        sep=" ",
        float_format="%.8e",
    )


def evol_diff_smani(bounds, cell, xi_modif, n_iter):
    all_params_ani = {}
    props_smaut = vect_props_smaut_test(
        load_variable_props(f"results_params/params_smaut_{cell}.txt")
    )

    xi_modif = right_artificial_xi(props_smaut, cell)
    typesim_to_loads = {
        "tension": 0,
        "compression": np.pi,
        "biaxial_compression": -3 * np.pi / 4,
        "tencomp": 3 * np.pi / 4,
        "biaxial_tension": np.pi / 4,
        "shear": np.pi / 2,
    }

    loss = partial(
        calc_cost_smani, list_typesim=typesim_to_loads, xi_modif=xi_modif, cell=cell
    )
    result = differential_evolution(
        loss,
        bounds,
        maxiter=n_iter,
        popsize=25,
        tol=1e-6,
        mutation=(0.7, 1.5),
        recombination=0.9,
        disp=True,
    )

    all_params_ani["smani"] = result.x
    df = pd.DataFrame(all_params_ani)
    df.to_csv(
        f"results_params/params_smani_{cell}.txt",
        index=False,
        sep=" ",
        float_format="%.8e",
    )


def run_homogeneisation(cell):
    bounds = [
        (5000, 12000),  # E
        (1.0, 11.0),  # C_A C_M
        (0.02, 0.12),  # Hmax
        (0, 100),  # sigma crit
        (0, 80),  # dT
        (0, 400),  # sigma caliber
        (-2.0, 2.0),
        (0.1, 5),  # b  # n
    ]

    evol_diff_smaut(bounds, cell=cell, n_iter=40)
    ######################################################Iso_surface################################################""
    props_smaut = vect_props_smaut_test(
        load_variable_props(f"results_params/params_smaut_{cell}.txt")
    )
    xi_modif = right_artificial_xi(props_smaut, cell)

    bounds = [
        (-2.0, 3.0),  # b
        (0.1, 3.0),  # n
        (0, 10),  # F
        (0, 10),  # L
        (0, 15),  # K
    ]

    evol_diff_smani(bounds, cell, xi_modif=xi_modif, n_iter=100)
