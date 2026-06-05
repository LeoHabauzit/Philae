import os
from tools_fea import *
import fedoo as fd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
import sys
from pathlib import Path
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
from tools_homogeneisation import (
    run_linear_homogenization,
    load_variable_props,
    vect_props_smaac,
)
from simuEF.tools_fea import read_props


fd.ModelingSpace("3D")

cell = "RhombicDodecahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)
props_init = read_props("simuEF/params_sma_init.txt")
print("props=", props_init)

os.remove("simuEF/train_fea.csv")
for j in range(1):
    mesh = fd.mesh.box_mesh(
        nx=2,
        ny=2,
        nz=2,  # 2 nodes/dir → 1 hex8 element, 8 corner nodes
        x_min=0,
        x_max=1,
        y_min=0,
        y_max=1,
        z_min=0,
        z_max=1,
        elm_type="hex8",
        name="Domain",
    )

    stress_lim = 40
    stress_target = np.random.uniform(-stress_lim * 10, 10 * stress_lim, (1, 6))
    stress_zeros = np.zeros((1, 6))
    print(stress_target)

    material = fd.constitutivelaw.Simcoon("SMAUT", props_init)

    # material = fd.constitutivelaw.ElasticIsotrop(E, nu, name="Material")
    wf = fd.weakform.StressEquilibrium(material, nlgeom=False)
    assembly = fd.Assembly.create(wf, mesh, name="Assembly")

    temp = 300.0
    if isinstance(temp, float):
        assembly.sv["Temp"] = temp * np.ones(assembly.n_gauss_points, order="F")

    pb = fd.problem.NonLinear(assembly)
    pb.set_solver("direct")
    pb.bc.add(fd.constraint.PeriodicBC(periodicity_type="small_strain"))

    ref_node = mesh.nearest_node(mesh.bounding_box.center)
    pb.bc.add("Dirichlet", ref_node, "Disp", 0)

    volume = mesh.bounding_box.volume  # = 1.0 mm³ for the unit cube

    pb.bc.add("Neumann", "E_xx", stress_target[0, 0] * volume)
    pb.bc.add("Neumann", "E_yy", stress_target[0, 1] * volume)
    pb.bc.add("Neumann", "E_zz", stress_target[0, 2] * volume)
    pb.bc.add(
        "Neumann", "E_xy", stress_target[0, 3] * volume
    )  # conjugate of γ_xy is σ_xy
    pb.bc.add("Neumann", "E_xz", stress_target[0, 4] * volume)
    pb.bc.add("Neumann", "E_yz", stress_target[0, 5] * volume)

    results = pb.add_output(
        "fea_1element",
        assembly,
        ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
    )
    pb.nlsolve(dt=0.2, tmax=1, t0=0, update_dt=True, print_info=1, interval_output=0.02)

    for k in range(6):
        pb.bc.remove(-1)
    pb.bc.add("Neumann", "E_xx", 0)
    pb.bc.add("Neumann", "E_yy", 0)
    pb.bc.add("Neumann", "E_zz", 0)
    pb.bc.add("Neumann", "E_xy", 0)  # conjugate of γ_xy is σ_xy
    pb.bc.add("Neumann", "E_xz", 0)
    pb.bc.add("Neumann", "E_yz", 0)

    pb.nlsolve(dt=0.2, tmax=2, t0=1, update_dt=True, print_info=2, interval_output=0.02)
    res = pb.get_results(
        # "test",
        assembly,
        ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
    )

    dataset = fd.read_data("fea_1element.fdz")
    n_iter = dataset.n_iter
    time = np.linspace(0, 2, n_iter + 1)
    stress_array = np.zeros((6, n_iter + 1))
    strain_array = np.zeros((6, n_iter + 1))
    mean_strain_array = np.zeros((6, n_iter + 1))
    mean_stress_array = np.zeros((6, n_iter + 1))

    all_components = ["XX", "YY", "ZZ", "XY", "XZ", "YZ"]
    for k in range(dataset.n_iter):
        dataset.load(k)

        for i, component in enumerate(all_components, start=0):
            data_stress = dataset.get_data(
                field="Fext(MeanStrain)", component=component, data_type="GaussPoint"
            )
            data_MeanStrain = dataset.get_data(
                field="MeanStrain", component=component, data_type="GaussPoint"
            )
            mean_stress_array[i, k + 1] = data_stress[0]
            mean_strain_array[i, k + 1] = data_MeanStrain[0]
    print(mean_strain_array.shape)
    save_data_csv(mean_stress_array, mean_strain_array, time, 1)

# plot_data_6D(mean_stress_array, mean_strain_array, time)

# n_simulations = 100
# save_data_csv(mean_stress_array, mean_strain_array, time, 1)
