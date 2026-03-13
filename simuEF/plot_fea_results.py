from pathlib import Path
import numpy as np

# import fedoo as fd
# import pyvista as pv
import matplotlib.pyplot as plt
from typing import NamedTuple
from simcoon import simmit as sim
import numpy.typing as npt
import os

typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}
list_cell={
    # 'RhombicDodecahedron40',
    'Cuboctahedron40',
    # 'Gyroid40'
    # 'Cuboctahedron40_finemesh'
}

fig, axs = plt.subplots(2, 3, figsize=(10, 8))
for i, typesim in enumerate(sorted(typesim_to_loads)):
    row = i // 3
    col = i % 3

    ax = axs[row, col]
    results_dir = typesim   
    for cellule in list_cell:
        if typesim=='shear':
            stress_array = np.loadtxt(
                f"datas_simu/{cellule}/SXY/data_{results_dir}/Stress_{results_dir}.txt"
            )
            strain_array = np.loadtxt(
                f"datas_simu/{cellule}/SXY/data_{results_dir}/MeanStrain_{results_dir}.txt"
            )

        else:
            stress_array = np.loadtxt(
                f"datas_simu/{cellule}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
            )
            strain_array = np.loadtxt(
                f"datas_simu/{cellule}/SXX/data_{results_dir}/MeanStrain_{results_dir}.txt"
            )
        xi_array = np.loadtxt(f"datas_simu/{cellule}/SXX/data_{results_dir}/Xi_{results_dir}.txt")

        ax.plot(
            xi_array,
            stress_array,
            label=f"{typesim}/{cellule}",
        )
    ax.set_title(f"Plot {typesim}")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True)
    ax.set_xlabel("E11[%]")
    ax.set_ylabel("S11 [MPa]")

plt.tight_layout()

plt.title(f"Plot {typesim}")
plt.legend(loc="upper left", fontsize=8)
plt.grid(True)
plt.xlabel("E11[%]")
plt.ylabel("S11 [MPa]")
plt.show()

