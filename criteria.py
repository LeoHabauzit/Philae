import numpy as np
from continum_mech import *


def mises_stress(v):
    # v = np.asarray(v, dtype=float)
    vdev = dev(v)
    vdev2 = vdev.copy()
    for i in range(3, 6):
        vdev2[i] = 2 * vdev2[i]  # notation ingénieur
    return np.sqrt(3 / 2 * np.sum(vdev * vdev2))


def drucker_stress(v, b, n):
    v = np.asarray(v, dtype=float)
    vm = mises_stress(v)

    J2_3_2 = pow(J2_stress(v), 3 / 2)
    m = 1 / n
    return vm * pow(1 + b * J3_stress(v) / J2_3_2, m)


def find_P_dfa(params):
    F = params[0]
    G = params[1]
    H = params[2]
    L = params[3]
    M = params[4]
    N = params[5]
    K = params[6]
    P = np.zeros((6, 6))

    P[0, 0] = G + H + K / 9
    P[0, 1] = -H + K / 9
    P[0, 2] = -G + K / 9

    P[1, 0] = -H + K / 9
    P[1, 1] = F + H + K / 9
    P[1, 2] = -F + K / 9

    P[2, 0] = -G + K / 9
    P[2, 1] = -F + K / 9
    P[2, 2] = F + G + K / 9

    # composantes de cisaillement
    P[3, 3] = 2 * N
    P[4, 4] = 2 * M
    P[5, 5] = 2 * L

    return P


def dfa_stress(v, params):
    P_dfa = find_P_dfa(params)
    return np.sqrt(np.sum(v * (P_dfa @ v)))


def drucker_ani_stress(v, b, n, params):
    v = np.asarray(v, dtype=float)
    v_dfa = dfa_stress(v, params)
    J2_3_2 = pow(J2_stress(v), 3 / 2)
    m = 1 / n
    return v_dfa * pow(1 + b * J3_stress(v) / J2_3_2, m)


def get_Phi_forward_SMA(props, v, xi, T, ani=0):
    """

    Si ani = 0 -> critère de drucker
    Si ani = 1 -> drucker anisotrope (dfa)
    """
    flagT = 0
    E_A = props[1]
    E_M = props[2]
    nu_A = props[3]
    nu_M = props[4]
    alphaA_iso = props[5]
    alphaM_iso = props[6]
    Hmin = props[7]
    Hmax = props[8]
    k1 = props[9]
    sigmacrit = props[10]

    C_A = props[11]
    C_M = props[12]
    Ms0 = props[13]
    Mf0 = props[14]
    As0 = props[15]
    Af0 = props[16]

    n1 = props[17]
    n2 = props[18]
    n3 = props[19]
    n4 = props[20]
    sigmacaliber = props[21]
    prager_b = props[22]
    prager_n = props[23]

    c_lambda = props[24]
    p0_lambda = props[25]
    n_lambda = props[26]
    alpha_lambda = props[27]

    F_dfa = props[28]
    G_dfa = props[29]
    H_dfa = props[30]
    L_dfa = props[31]
    M_dfa = props[32]
    N_dfa = props[33]
    K_dfa = props[34]
    DFA_params = np.array([F_dfa, G_dfa, H_dfa, L_dfa, M_dfa, N_dfa, K_dfa])

    Dalpha = (alphaM_iso - alphaA_iso) * np.array((1, 1, 1, 0, 0, 0))

    K_A = E_A / (3.0 * (1.0 - 2 * nu_A))
    mu_A = E_A / (2.0 * (1.0 + nu_A))

    K_M = E_M / (3.0 * (1.0 - 2 * nu_M))
    mu_M = E_M / (2.0 * (1.0 + nu_M))
    M_A = M_iso_Kmu(K_A, mu_A)
    M_M = M_iso_Kmu(K_M, mu_M)
    DM = M_M - M_A

    DM_sig = DM @ v

    if sigmacaliber > sigmacrit:
        sigmastar = sigmacaliber - sigmacrit
    else:
        sigmastar = 0.0

    ###Axif
    Hcurstar = Hmin + (Hmax - Hmin) * (1.0 - np.exp(-k1 * sigmastar))
    dHcurstar = (Hmax - Hmin) * k1 * np.exp(-k1 * sigmastar)

    if flagT == 0:
        MsSmooth = 0.5 * Ms0 * (
            1.0 + (n1 + 1) * pow(2.0, -n1) + (n2 - 1) * pow(2.0, -n2)
        ) / (n1 * pow(2.0, -n1) + n2 * pow(2.0, -n2)) + 0.5 * Mf0 * (
            -1.0 + (n1 - 1) * pow(2.0, -n1) + (n2 + 1) * pow(2.0, -n2)
        ) / (n1 * pow(2.0, -n1) + n2 * pow(2.0, -n2))
        MfSmooth = 0.5 * Ms0 * (
            -1.0 + (n1 + 1) * pow(2.0, -n1) + (n2 - 1) * pow(2.0, -n2)
        ) / (n1 * pow(2.0, -n1) + n2 * pow(2.0, -n2)) + 0.5 * Mf0 * (
            1.0 + (n1 - 1) * pow(2.0, -n1) + (n2 + 1) * pow(2.0, -n2)
        ) / (n1 * pow(2.0, -n1) + n2 * pow(2.0, -n2))
        AsSmooth = 0.5 * As0 * (
            1.0 + (n3 - 1) * pow(2.0, -n3) + (n4 + 1) * pow(2.0, -n4)
        ) / (n3 * pow(2.0, -n3) + n4 * pow(2.0, -n4)) + 0.5 * Af0 * (
            -1.0 + (n3 + 1) * pow(2.0, -n3) + (n4 - 1) * pow(2.0, -n4)
        ) / (n3 * pow(2.0, -n3) + n4 * pow(2.0, -n4))
        AfSmooth = 0.5 * As0 * (
            -1.0 + (n3 - 1) * pow(2.0, -n3) + (n4 + 1) * pow(2.0, -n4)
        ) / (n3 * pow(2.0, -n3) + n4 * pow(2.0, -n4)) + 0.5 * Af0 * (
            1.0 + (n3 + 1) * pow(2.0, -n3) + (n4 - 1) * pow(2.0, -n4)
        ) / (n3 * pow(2.0, -n3) + n4 * pow(2.0, -n4))

    else:
        MsSmooth = Ms0
        MfSmooth = Mf0
        AsSmooth = As0
        AfSmooth = Af0

    rhoDs0 = (
        -2.0
        * C_M
        * C_A
        * (Hcurstar + sigmacaliber * (dHcurstar + (1 / E_M - 1 / E_A)))
        / (C_M + C_A)
    )

    rhoDE0 = 0.5 * rhoDs0 * (MsSmooth + AfSmooth)
    D = (
        (C_M - C_A)
        * (Hcurstar + sigmacaliber * (dHcurstar + (1 / E_M - 1 / E_A)))
        / ((C_A + C_M) * (Hcurstar + sigmacaliber * dHcurstar))
    )
    a1 = rhoDs0 * (MfSmooth - MsSmooth)
    a2 = rhoDs0 * (AsSmooth - AfSmooth)
    a3 = -0.25 * a1 * (1 + 1.0 / (n1 + 1.0) - 1.0 / (n2 + 1.0)) + 0.25 * a2 * (
        1.0 + 1.0 / (n3 + 1.0) - 1.0 / (n4 + 1.0)
    )
    Y0t = 0.5 * rhoDs0 * (MsSmooth - AfSmooth) - a3

    if n1 == 1.0 and n2 == 1.0:
        HfF = a1 * xi + a3
    else:
        if xi > 0.0 and 1.0 - xi > 0.0:
            HfF = 0.5 * a1 * (1.0 + pow(xi, n1) - pow(1.0 - xi, n2)) + a3
        elif xi <= 0.0:
            HfF = 0.5 * a1 - 1.0 + a3

        elif xi >= 1.0:
            HfF = a1 + a3

    A_xiF = rhoDs0 * T - rhoDE0 + 0.5 * sum(v * DM_sig) - HfF

    ###Lambda 1
    h = 1.0 - c_lambda
    if xi <= h:
        lambda1 = p0_lambda * xi * pow(1.0 - xi, -n_lambda)
    else:
        lambda1 = p0_lambda * (
            (pow(1.0 - h, -n_lambda) + n_lambda * h * pow(1.0 - h, -n_lambda - 1.0))
            * (xi - h)
            + h * pow(1.0 - h, -n_lambda)
        ) + alpha_lambda * pow(xi - h, 2.0)

    if mises_stress(v) > sigmacrit:
        sigmastar = mises_stress(v) - sigmacrit
    else:
        sigmastar = 0.0

    Hcur = Hmin + (Hmax - Hmin) * (1.0 - np.exp(-1.0 * k1 * sigmastar))
    YtF = Y0t + D * Hcur * mises_stress(v)

    ##Calcul du critère
    if ani == 0:
        PhihatF = Hcur * drucker_stress(
            v,
            prager_b,
            prager_n,
        )
    else:
        PhihatF = Hcur * drucker_ani_stress(
            v,
            prager_b,
            prager_n,
            DFA_params,
        )

    return PhihatF + A_xiF - lambda1 - YtF
