import os

os.environ["OMP_NUM_THREADS"] = "1"  # avant tout import de simcoon/fedoo
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

try:
    from scipy.sparse.linalg._dsolve.linsolve import useUmfpack as _scipy_uu

    if not hasattr(_scipy_uu, "u"):
        _scipy_uu.u = True  # active umfpack, évite le crash dans base.py
except Exception:
    pass
import fedoo as fd
import numpy as np

import sys
from pathlib import Path
import time
from tools_fea import *

sys.path.append(str(Path(__file__).resolve().parent.parent))
from tools_homogeneisation import (
    run_linear_homogenization,
    load_variable_props,
    vect_props_smaac,
)
from simuEF.tools_fea import read_props, plot_data_6D


def save_stress_target(stress_target, sim_id, csv_file="simuEF/stress_target.csv"):
    s11 = stress_target[:, 0]
    s22 = stress_target[:, 1]
    s33 = stress_target[:, 2]
    s12 = stress_target[:, 3]
    s13 = stress_target[:, 4]
    s23 = stress_target[:, 5]
    # print(s11, s22, s33, s12, s13, s23)
    df = pd.DataFrame(
        {
            "stress_xx": s11,
            "stress_yy": s22,
            "stress_zz": s33,
            "stress_xy": s12,
            "stress_xz": s13,
            "stress_yz": s23,
            # "timestep": time,
            "simulation_load_id": sim_id,
        }
    )
    if not os.path.isfile(csv_file):
        df.to_csv(csv_file, index=False)
    else:
        # append sans réécrire les colonnes
        df.to_csv(csv_file, mode="a", header=False, index=False)


# t1 = time.time()
fd.ModelingSpace("3D")

cell = "RhombicDodecahedron40"
meshfile = f"simuEF/cellules/{cell}.vtk"
n_simulations = 1
os.remove("simuEF/train_fea.csv") if os.path.exists("simuEF/train_fea.csv") else None
os.remove("simuEF/stress_target.csv") if os.path.exists(
    "simuEF/stress_target.csv"
) else None
sim_ids = np.random.choice(np.arange(1, 10001), size=n_simulations, replace=False)

for j in range(n_simulations):
    print(f"Running random {j + 1} FE computation")
    # props_cubic = run_linear_homogenization(f"{cell}")
    # props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
    # finalprops = vect_props_smaac(props_var, props_cubic)
    props_init = read_props("simuEF/params_sma_init.txt")

    mesh = fd.Mesh.read(meshfile)

    stress_lim = 150
    sign = np.random.choice([-1, 1], size=(1, 6))
    magnitude = np.random.uniform(20, stress_lim, size=(1, 6))

    stress_target = sign * magnitude
    stress_zeros = np.zeros((1, 6))
    save_stress_target(stress_target, sim_ids[j])

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
        "simuEF/fea_RandCell",
        assembly,
        ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
    )
    pb.set_nr_criterion("Displacement", tol=1e-4, max_subiter=20)
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
    # t2 = time.time()
    # print("Time =", t2 - t1, "s")
    res = pb.get_results(
        # "test",
        assembly,
        ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
    )

    dataset = fd.read_data("simuEF/fea_RandCell.fdz")
    n_iter = dataset.n_iter
    time = np.linspace(0, 1, n_iter + 1)
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
    save_data_csv(mean_stress_array, mean_strain_array, time, sim_ids[j])

# plot_data_6D(mean_stress_array, mean_strain_array, time)
