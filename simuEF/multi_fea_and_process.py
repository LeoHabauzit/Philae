import os
import fedoo as fd

# from fedoo.core.boundary_conditions import ListBC, BoundaryCondition

# import Pat

from tools_fea import *

np.float_ = np.float64
cell = "RhombicDodecahedron40"
# meshfile = f"cellules/{cell}.vtk"

material_law = "SMAUT"

props = read_props("simuEF/params_sma_init.txt")

typesim_to_loads = define_typesim_to_loads(0.05)

for typesim in typesim_to_loads.keys():
    load = typesim_to_loads.get(typesim)
    cell_fea(props, material_law, typesim, load, cell)
    process_data_fea(typesim, cell)
    erase_fea_file(typesim)

plot_results_fea(cell, typesim_to_loads)
plt.show()
