import numpy as np
from scipy.interpolate import PchipInterpolator


def dev(v):
    v = np.asarray(v, dtype=float)
    mean_stress = (v[0] + v[1] + v[2]) / 3
    vdev = v.copy()
    vdev[0:3] -= mean_stress
    return vdev


def mises_strain(v):
    vdev = dev(v)
    vdev2 = vdev.copy()

    vdev2[3:] *= 0.5
    return np.sqrt(2.0 / 3.0 * np.sum(vdev * vdev2))


def v2t(v):
    v = np.asarray(v, dtype=float)
    m = np.zeros((3, 3))
    for i in range(3):
        m[i, i] = v[i]
        for j in range(i + 1, 3):
            m[i, j] = v[i + j + 2]
            m[j, i] = v[i + j + 2]
    return m


def J2_stress(v):
    v = np.asarray(v, dtype=float)
    vdev = dev(v)
    vdev2 = vdev.copy()
    for i in range(3, 6):
        vdev2[i] = 2 * vdev2[i]
    return 0.5 * np.sum(vdev * vdev2)


def J3_stress(v):
    v = np.asarray(v, dtype=float)
    vdev = dev(v)
    mat_vdev = v2t(vdev)
    mat_cub = np.linalg.matrix_power(mat_vdev, 3)

    return (1 / 3) * (mat_cub[0, 0] + mat_cub[1, 1] + mat_cub[2, 2])


def Ivol():
    """Tenseur volumique identité en notation Voigt 6x6"""
    Iv = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            Iv[i, j] = 1 / 3
    return Iv


def Idev2():
    """Tenseur deviatorique identité en notation Voigt 6x6"""
    Id = np.zeros((6, 6))
    # partie normale (indices 0,1,2)
    for i in range(3):
        for j in range(3):
            if i == j:
                Id[i, j] = 2 / 3
            else:
                Id[i, j] = -1 / 3
    # partie cisaillement (indices 3,4,5)
    for i in range(3, 6):
        Id[i, i] = 2
    return Id


def M_iso_Kmu(C1, C2):
    K = C1
    mu = C2
    return 1 / (3.0 * K) * Ivol() + 1 / (2.0 * mu) * Idev2()


def stress_vector_from_polar(r, theta, plane):
    """
    Génère le vecteur de contraintes à partir de (r, theta)
    plane :
        's11-s22' -> (σ11, σ22)
        's11-s12' -> (σ11, σ12)
    """
    if plane == "s11-s22":
        return np.array([r * np.cos(theta), r * np.sin(theta), 0, 0, 0, 0])
    elif plane == "s11-s12":
        return np.array([r * np.cos(theta), 0, 0, r * np.sin(theta), 0, 0])
    else:
        raise ValueError("plane doit être 's11-s22' ou 's11-s12'")


def prepare_interp(stress, xi):
    # tri
    idx = np.argsort(stress)
    stress = stress[idx]
    xi = xi[idx]

    # suppression des doublons
    stress_unique, indices = np.unique(stress, return_index=True)
    xi_unique = xi[indices]

    return PchipInterpolator(stress_unique, xi_unique)
