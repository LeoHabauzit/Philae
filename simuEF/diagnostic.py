"""diagnostic_smaut_multigp.py"""

# ── Fix OpenMP deadlock simcoon / macOS ARM64 ────────────────────────────────
import os

os.environ["OMP_NUM_THREADS"] = "1"  # avant tout import de simcoon/fedoo
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
# ─────────────────────────────────────────────────────────────────────────────
try:
    from scipy.sparse.linalg._dsolve.linsolve import useUmfpack as _uu

    if not hasattr(_uu, "u"):
        _uu.u = True
except Exception:
    pass

import numpy as np
from simcoon import simmit as sim

props = np.asfortranarray(
    np.c_[
        np.array(
            [
                0,
                67538.0,
                67538.0,
                0.349,
                0.349,
                1e-6,
                1e-6,
                0.0,
                0.0418,
                0.021,
                0.0,
                10.0,
                10.0,
                250.0,
                230.0,
                260.0,
                280.0,
                0.2,
                0.2,
                0.2,
                0.2,
                300.0,
                1.4,
                2.0,
                1e-6,
                1e-3,
                1.0,
                1e8,
            ]
        )
    ]
)  # shape (28, 1)

# ── Test avec N points de Gauss croissant ────────────────────────────────────
for n_gp in [1, 2, 4, 8]:
    print(f"Test n_gp = {n_gp} ...", end=" ", flush=True)
    try:
        statev = np.zeros((50, n_gp), order="F")
        zeros_6 = np.zeros((6, n_gp), order="F")
        DR = np.empty((3, 3, n_gp), order="F")
        DR[...] = np.eye(3).reshape(3, 3, 1)
        F = np.array([])  # nlgeom=False → F vide
        Wm = np.zeros((4, n_gp), order="F")
        temp = 300.0 * np.ones(n_gp)  # shape (n_gp,)

        result = sim.umat(
            "SMAUT",
            zeros_6,
            zeros_6,  # strain_start, dstrain
            F,
            F,  # F0, F1 (vides car nlgeom=False)
            zeros_6,  # stress_start
            DR,
            props,
            statev,
            0,
            0,  # time, dtime
            Wm,
            temp,
            ndi=3,
        )
        print(f"✅ OK  (TangentMatrix shape: {result[3].shape})")

    except Exception as e:
        print(f"❌ CRASH : {type(e).__name__}: {e}")
        break

print("Diagnostic terminé.")
