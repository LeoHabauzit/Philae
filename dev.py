from tools_homogeneisation import *

cell = "RhombicDodecahedron40"
typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}

bounds = [
    (0.02, 0.12),  # Hmax
    (0, 100),  # sigma crit
    (-2.0, 2.0),  # b
    (0.1, 5),  # n
    (250, 280),  # Mf0
    (20, 50),  # dsf
    (0, 10),  # F
    (0, 10),  # L
    (0, 10),  # K
]


evol_diff_strain(bounds, cell, n_iter=40)
