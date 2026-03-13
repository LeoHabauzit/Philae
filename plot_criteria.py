import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import matplotlib.lines as mlines
from continum_mech import *
from criteria import *


def radius_for_von_mises(theta, sigma_y, plane):
    """Rayon r tel que σ_VM = sigma_y pour un angle theta"""

    def f(r):
        v = stress_vector_from_polar(r, theta, plane)
        return mises_stress(v) - sigma_y

    sol = root_scalar(f, bracket=[0, 10 * sigma_y], method="brentq")
    return sol.root


def radius_for_drucker(props, xi, T, theta, plane):
    """Rayon r tel que σ_drucker = sigma_y pour un angle theta"""
    # print(theta)

    def f(r):
        v = stress_vector_from_polar(r, theta, plane)
        return get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=0)

    a = 1e-3
    b = 10.0
    fa = f(a)
    fb = f(b)

    max_expand = 20
    i = 0
    while fa * fb > 0 and i < max_expand and a > 1e-10:
        a *= 0.5
        fa = f(a)

    while fa * fb > 0 and i < max_expand and b < 1e8:
        b *= 2
        fb = f(b)

    if fa * fb > 0:
        return fa * fb

    else:
        sol = root_scalar(f, bracket=[a, b], method="brentq")
        return sol.root


def radius_for_drucker_ani(props, xi, T, theta, plane):
    """Rayon r tel que σ_drucker = sigma_y pour un angle theta"""
    # print(theta)

    def f(r):
        v = stress_vector_from_polar(r, theta, plane)
        return get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=1)

    a = 1e-3
    b = 10.0
    fa = f(a)
    fb = f(b)

    max_expand = 20
    i = 0
    while fa * fb > 0 and i < max_expand and a > 1e-10:
        a *= 0.5
        fa = f(a)

    while fa * fb > 0 and i < max_expand and b < 1e8:
        b *= 2
        fb = f(b)

    if fa * fb > 0:
        return fa * fb * 10
    else:
        sol = root_scalar(f, bracket=[a, b], method="brentq")
        return sol.root


def radius_to_find_xi_lim(props, xi, T, theta, plane):
    """Rayon r tel que σ_drucker = sigma_y pour un angle theta"""
    # print(theta)

    def f(r):
        v = stress_vector_from_polar(r, theta, plane)
        return get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=1)

    a = 1e-3
    b = 10.0
    fa = f(a)
    fb = f(b)

    max_expand = 20
    i = 0
    while fa * fb > 0 and i < max_expand and a > 1e-10:
        a *= 0.5
        fa = f(a)

    while fa * fb > 0 and i < max_expand and b < 1e8:
        b *= 2
        fb = f(b)

    if fa * fb > 0 and f(a) > 0:
        return -fa * fb * 10
    else:
        sol = root_scalar(f, bracket=[a, b], method="brentq")
        return sol.root


def plot_von_mises_radius(ax, sigma_y, plane="s11-s22", npts=400):
    """Trace le critère de Von Mises dans le plan choisi"""

    theta = np.linspace(0, 2 * np.pi, npts)
    r = np.array([radius_for_von_mises(t, sigma_y, plane) for t in theta])

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    ax.plot(x, y, label="Von Mises", color="red")


def plot_drucker_radius(ax, props, xi, T, plane="s11-s22", npts=400):
    theta = np.linspace(0, 2 * np.pi, npts)
    r = np.array(
        [
            radius_for_drucker(props=props, xi=xi, T=T, theta=t, plane=plane)
            for t in theta
        ]
    )

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    ax.plot(x, y, label="Drucker", color="green")


def plot_drucker_ani_radius(ax, props, xi, T, plane="s11-s22", npts=400):
    theta = np.linspace(0, 2 * np.pi, npts)
    r = np.array(
        [
            radius_for_drucker_ani(props=props, xi=xi, T=T, theta=t, plane=plane)
            for t in theta
        ]
    )

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    ax.plot(x, y, label="Drucker Ani", color="orange")


def plot_dfa_contour(ax, sigma_y, params, plane="s11-s22", npts=400):
    """Trace le critère de Von Mises dans le plan choisi"""

    s = np.linspace(-800, 800, 600)
    X, Y = np.meshgrid(s, s)

    F = np.zeros_like(X)
    if plane == "s11-s22":
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                v = np.array([X[i, j], Y[i, j], 0, 0, 0, 0])
                F[i, j] = dfa_stress(v, params) - sigma_y
    elif plane == "s11-s12":
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                v = np.array([X[i, j], 0, 0, Y[i, j], 0, 0])
                F[i, j] = dfa_stress(v, params) - sigma_y

    cs = ax.contour(X, Y, F, levels=[0], colors="green")
    ax.clabel(cs, fmt="Dfa")


def plot_drucker_contour(ax, props, xi, T, plane="s11-s22", npts=400):
    """Trace le critère de Von Mises dans le plan choisi"""

    s = np.linspace(-800, 800, 600)
    X, Y = np.meshgrid(s, s)

    F = np.zeros_like(X)
    if plane == "s11-s22":
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                v = np.array([X[i, j], Y[i, j], 0, 0, 0, 0])
                F[i, j] = get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=0)
    elif plane == "s11-s12":
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                v = np.array([X[i, j], 0, 0, Y[i, j], 0, 0])
                F[i, j] = get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=0)

    cs = ax.contour(X, Y, F, levels=[0], colors="k")
    ax.clabel(cs, fmt="Drucker")


def plot_drucker_ani_contour(ax, props, xi, T, plane="s11-s22", npts=400):
    """Trace le critère de Von Mises dans le plan choisi"""

    s = np.linspace(-800, 800, 600)
    X, Y = np.meshgrid(s, s)

    F = np.zeros_like(X)
    if plane == "s11-s22":
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                v = np.array([X[i, j], Y[i, j], 0, 0, 0, 0])
                F[i, j] = get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=1)
    elif plane == "s11-s12":
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                v = np.array([X[i, j], 0, 0, Y[i, j], 0, 0])
                F[i, j] = get_Phi_forward_SMA(props=props, v=v, xi=xi, T=T, ani=1)

    # tracé
    cs = ax.contour(X, Y, F, levels=[0], colors="orange")
    # ax.clabel(cs, fmt="Dfa_Drucker")
