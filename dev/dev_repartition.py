from pathlib import Path
import numpy as np
from tools_homogeneisation import *

import os
import fedoo as fd

# from fedoo.core.boundary_conditions import ListBC, BoundaryCondition

# import Pat

from simuEF.tools_fea import *

cell = "Cuboctahedron40"

typesim = "biaxial_compression"
process_element_repartition(typesim, cell)
