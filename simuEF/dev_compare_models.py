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
import pandas as pd
from tools_fea import *


sys.path.append(str(Path(__file__).resolve().parent.parent))
from tools_homogeneisation import *
from lstm.tools_database import generer_path_txt, add_zeros_to_arrays
from Umat.loi_sma import umat_sma_random


cell = "RhombicDodecahedron40"
props_cubic = run_linear_homogenization(f"{cell}")
props_var = load_variable_props(f"results_params/params_strain_{cell}.txt")
finalprops = vect_props_smaac(props_var, props_cubic)


def simu_fea(stress_target):
    # A LANCER AVEC FEDOO 0.7.0 (ou plus récent) et SIMCOON 1.10.0 (ou plus récent)
    fd.ModelingSpace("3D")
    cell = "RhombicDodecahedron40"
    meshfile = f"simuEF/cellules/{cell}.vtk"
    n_simulations = 1
    # os.remove("simuEF/train_fea.csv") if os.path.exists("simuEF/train_fea.csv") else None
    # os.remove("simuEF/stress_target.csv") if os.path.exists(
    #     "simuEF/stress_target.csv"
    # ) else None
    sim_ids = np.random.choice(np.arange(1, 10001), size=n_simulations, replace=False)
    failed_sims = []
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
        # save_stress_target(stress_target, sim_ids[j])

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
            "simuEF/fdz_files/fea_compare",
            assembly,
            ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
        )
        pb.set_nr_criterion("Displacement", tol=1e-4, max_subiter=20)
        try:
            pb.nlsolve(
                dt=0.2, tmax=1, t0=0, update_dt=True, print_info=1, interval_output=0.02
            )
        except (RuntimeError, NameError) as e:
            print(f"[WARNING] Simulation {j + 1} (id={sim_ids[j]}) échouée : {e}")
            failed_sims.append({"index": j, "sim_id": sim_ids[j], "error": str(e)})
            continue

        for k in range(6):
            pb.bc.remove(-1)
        pb.bc.add("Neumann", "E_xx", 0)
        pb.bc.add("Neumann", "E_yy", 0)
        pb.bc.add("Neumann", "E_zz", 0)
        pb.bc.add("Neumann", "E_xy", 0)  # conjugate of γ_xy is σ_xy
        pb.bc.add("Neumann", "E_xz", 0)
        pb.bc.add("Neumann", "E_yz", 0)

        try:
            pb.nlsolve(
                dt=0.2, tmax=2, t0=1, update_dt=True, print_info=2, interval_output=0.02
            )
        except (RuntimeError, NameError) as e:
            print(
                f"[WARNING] Simulation {j + 1} (id={sim_ids[j]}) échouée lors de la deuxième résolution : {e}"
            )
            failed_sims.append({"index": j, "sim_id": sim_ids[j], "error": str(e)})
            continue
        # t2 = time.time()
        # print("Time =", t2 - t1, "s")
        res = pb.get_results(
            # "test",
            assembly,
            ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
        )

        dataset = fd.read_data("simuEF/fdz_files/fea_compare.fdz")
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
                    field="Fext(MeanStrain)",
                    component=component,
                    data_type="GaussPoint",
                )
                data_MeanStrain = dataset.get_data(
                    field="MeanStrain", component=component, data_type="GaussPoint"
                )
                mean_stress_array[i, k + 1] = data_stress[0]
                mean_strain_array[i, k + 1] = data_MeanStrain[0]
        save_data_csv(
            mean_stress_array,
            mean_strain_array,
            time,
            sim_ids[j],
            "simuEF/csv_files/compare_fea.csv",
        )


def simu_umat(
    stress_target,
    props=finalprops,
    umat_name="SMAAC",
    csv_file="simuEF/csv_files/compare_umat.csv",
):
    generer_path_txt(
        filename="Umat/data/path_random.txt",
        temperature_initiale=300,
        Dn_inc=0.02,
        stress_lim=1,
        stress_target=stress_target,
    )
    umat_sma_random(props, umat_name)
    os.remove("Umat/data/path_random.txt")
    outputfile_global = f"Umat/results_{umat_name}/results_random_global-0.txt"

    time, e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
        outputfile_global,
        usecols=(4, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
        unpack=True,
    )
    sim_id = np.full(len(time), 1)
    (time, e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23, sim_id) = (
        add_zeros_to_arrays(
            [
                time,
                e11,
                e22,
                e33,
                e12,
                e13,
                e23,
                s11,
                s22,
                s33,
                s12,
                s13,
                s23,
                sim_id,
            ]
        )
    )

    df = pd.DataFrame(
        {
            "total_strain_xx": e11,
            "total_strain_yy": e22,
            "total_strain_zz": e33,
            "total_strain_xy": e12,
            "total_strain_xz": e13,
            "total_strain_yz": e23,
            "stress_xx": s11,
            "stress_yy": s22,
            "stress_zz": s33,
            "stress_xy": s12,
            "stress_xz": s13,
            "stress_yz": s23,
            "timestep": time,
            "simulation_load_id": sim_id,
        }
    )

    if not os.path.isfile(csv_file):
        df.to_csv(csv_file, index=False)
    else:
        # append sans réécrire les colonnes
        df.to_csv(csv_file, mode="a", header=False, index=False)


stress_target = np.zeros((1, 6))

stress_target[0, 0] = 65.243535
stress_target[0, 1] = -94.699172
stress_target[0, 2] = -123.286809
stress_target[0, 3] = -64.226
stress_target[0, 4] = -149.743054
stress_target[0, 5] = 89.158459


# simu_fea(stress_target)
# simu_umat(stress_target)


def plot_compare(list_csv, i=0):
    j = 101 * i
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    axs[0, 0].set_title("ε11-σ11")
    axs[0, 0].set_ylabel("σ [MPa]")

    axs[0, 0].grid(True)
    axs[0, 1].set_title("ε22-σ22")
    axs[0, 1].set_ylabel("σ [MPa]")
    axs[0, 1].grid(True)

    axs[0, 2].set_title("ε33-σ33")
    axs[0, 2].set_ylabel("σ [MPa]")
    axs[0, 2].grid(True)

    axs[1, 0].set_xlabel("ε[-]")
    axs[1, 0].set_ylabel("σ [MPa]")
    axs[1, 0].grid(True)

    axs[1, 1].set_title("ε13-σ13")
    axs[1, 1].set_xlabel("ε[-]")
    axs[1, 1].grid(True)

    axs[1, 2].set_title("ε23-σ23")
    axs[1, 2].set_xlabel("ε[-]")
    axs[1, 2].grid(True)

    for filename in list_csv:
        print(filename)
        if filename.endswith("test_fea_dataset_23.csv"):
            label = "FEA"
            color = "tab:blue"
        elif filename.endswith("compare_umat.csv"):
            label = "UMAT"
            color = "tab:orange"

        elif filename.endswith("compare_lstm_tuned.csv"):
            label = "Test LSTM"
            color = "tab:green"
        df = pd.read_csv(filename)

        e11 = df["total_strain_xx"].values[j : j + 101]
        e22 = df["total_strain_yy"].values[j : j + 101]
        e33 = df["total_strain_zz"].values[j : j + 101]
        e12 = df["total_strain_xy"].values[j : j + 101]
        e13 = df["total_strain_xz"].values[j : j + 101]
        e23 = df["total_strain_yz"].values[j : j + 101]
        s11 = df["stress_xx"].values[j : j + 101]
        s22 = df["stress_yy"].values[j : j + 101]
        s33 = df["stress_zz"].values[j : j + 101]
        s12 = df["stress_xy"].values[j : j + 101]
        s13 = df["stress_xz"].values[j : j + 101]
        s23 = df["stress_yz"].values[j : j + 101]
        # time = df["timestep"].values[j : j + 101]

        axs[0, 0].plot(e11, s11, color=color, label=label)

        axs[0, 1].plot(e22, s22, color=color, label=label)

        axs[0, 2].plot(e33, s33, color=color, label=label)

        axs[1, 0].plot(e12, s12, color=color, label=label)

        axs[1, 1].plot(e13, s13, color=color, label=label)

        axs[1, 2].plot(e23, s23, color=color, label=label)

        fig.suptitle("Courbes contrainte-déformation", fontsize=16)

    axs[0, 0].legend()
    axs[0, 1].legend()
    axs[0, 2].legend()
    axs[1, 0].legend()
    axs[1, 1].legend()
    axs[1, 2].legend()
    plt.tight_layout()
    plt.show()


plot_compare(
    list_csv=[
        "lstm/dataset/test_fea_dataset_23.csv",
        "simuEF/csv_files/compare_umat.csv",
        "simuEF/csv_files/compare_lstm_tuned.csv",
    ],
    i=0,
)
