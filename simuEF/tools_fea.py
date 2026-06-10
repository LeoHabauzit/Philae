import numpy as np
import fedoo as fd
import os
from pathlib import Path
import matplotlib.pyplot as plt
from simcoon import simmit as sim
import pandas as pd


def dev_fea(v):
    v = np.asarray(v, dtype=float)
    mean_stress = (v[0] + v[1] + v[2]) / 3
    vdev = v.copy()
    vdev[0:3] -= mean_stress
    return vdev


def mises_strain_fea(v):
    vdev = dev_fea(v)
    vdev2 = vdev.copy()

    vdev2[3:] *= 0.5
    return np.sqrt(2.0 / 3.0 * np.sum(vdev * vdev2))


def read_props(filename):
    values = {}

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # ignorer lignes vides ou commentaires
            if not line or line.startswith("#"):
                continue

            name, expr = [x.strip() for x in line.split("=", 1)]

            # évaluer l'expression en utilisant les variables déjà définies
            values[name] = eval(expr, {}, values)

    props = np.array(
        [
            values["flagT"],
            values["E_A"],
            values["E_M"],
            values["nu_A"],
            values["nu_M"],
            values["alphaA"],
            values["alphaM"],
            values["Hmin"],
            values["Hmax"],
            values["k1"],
            values["sigmacrit"],
            values["C_A"],
            values["C_M"],
            values["Ms0"],
            values["Mf0"],
            values["As0"],
            values["Af0"],
            values["n1"],
            values["n2"],
            values["n3"],
            values["n4"],
            values["sigmacaliber"],
            values["b_prager"],
            values["n_prager"],
            values["c_lambda"],
            values["p0_lambda"],
            values["n_lambda"],
            values["alpha_lambda"],
        ]
    )

    return props


def define_typesim_to_loads(strain_value):
    tensile_load = [strain_value, 0, 0]  # 4% tensile strain along x
    biaxial_tension_load = [
        strain_value,
        strain_value,
        0,
    ]  # 4% tensile strain along x and y
    compression_load = [-strain_value, 0, 0]  # 5% compression strain along x
    biaxial_compression_load = [
        -strain_value,
        -strain_value,
        0,
    ]  # 5% compression strain along x and y
    tension_compression_load = [
        strain_value,
        -strain_value,
        0,
    ]  # 5% tensile strain along x compression along y
    shear_load = [0, 0, strain_value]  # 5% shear strain in plane (x,y) (/!\ 2*gamma!!!)

    typesim_to_loads = {
        "tension": tensile_load,
        "biaxial_tension": biaxial_tension_load,
        "compression": compression_load,
        "biaxial_compression": biaxial_compression_load,
        "tencomp": tension_compression_load,
        "shear": shear_load,
    }
    return typesim_to_loads


def cell_fea(props, material_law, typesim, load_typesim, cell):
    meshfile = f"simuEF/cellules/{cell}.vtk"

    print("Running " + typesim + " FE computation")
    results_dir = str(Path(__file__).parent / typesim)
    output_file = typesim
    if not (os.path.isdir(results_dir)):
        os.mkdir(results_dir)

    output_file_ext: str = "fdz"
    temp = 300.0

    fd.ModelingSpace("3D")

    mesh = fd.Mesh.read(meshfile)
    bounds = mesh.bounding_box

    center = mesh.nearest_node(mesh.bounding_box.center)
    print(material_law)
    material = fd.constitutivelaw.Simcoon(material_law, props)
    weakform = fd.weakform.StressEquilibrium(material, nlgeom=False)
    assembly = fd.Assembly.create(weakform, mesh)
    if isinstance(temp, float):
        assembly.sv["Temp"] = temp * np.ones(assembly.n_gauss_points, order="F")

    pb = fd.problem.NonLinear(assembly)
    pb.set_solver("direct")
    pb.add_output(
        results_dir + "/" + output_file,
        assembly,
        [
            "Disp",
            "Stress",
            "Strain",
            "Fext",
            "Statev",
            "MeanStrain",
            "Fext(MeanStrain)",
        ],
        file_format=output_file_ext,
        compressed=True,
    )

    pb.bc.add(fd.constraint.PeriodicBC(periodicity_type="small_strain", dim=3))

    strain_components = ["E_xx", "E_yy", "E_xy"]
    load = load_typesim
    for comp, value in zip(strain_components, load):
        if value != 0:
            print(f"pb.bc.add(Dirichlet, {comp}, {value})")
            pb.bc.add("Dirichlet", comp, value)

    center = mesh.nearest_node(mesh.bounding_box.center)
    pb.bc.add("Dirichlet", center, "Disp", 0)

    pb.nlsolve(dt=0.2, tmax=1, update_dt=True, print_info=1, interval_output=0.01)


def process_element_repartition(typesim, cell):
    base_dir = Path("datas_simu") / cell

    results_dir = typesim
    dataset = fd.read_data(f"simuEF/{typesim}/{typesim}.fdz")
    if typesim == "shear":
        list_component = {"XY"}
    elif typesim == "tencomp":
        list_component = {"XX", "YY"}
    else:
        list_component = {"XX"}
    for component in list_component:
        mesh_volume = dataset.mesh.to_pyvista().volume
        rve_volume = dataset.mesh.bounding_box.volume
        density = mesh_volume / rve_volume

        n_iter = dataset.n_iter
        # stress_array = np.zeros(n_iter + 1)
        xi_array = np.zeros(n_iter + 1)
        # meanStrain_array = np.zeros(n_iter + 1)

        # transformation_strain_array = np.zeros(n_iter + 1)
        # et_arrays = np.zeros((6, n_iter + 1))

        print(n_iter)

        data_dir = base_dir / f"S{component}" / f"data_{results_dir}"
        data_dir.mkdir(parents=True, exist_ok=True)

        # list_iter = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        fig, axes = plt.subplots(5, 5, figsize=(10, 10))
        axes = axes.flatten()
        k = 0
        for i in range(n_iter):
            dataset.load(i)
            print(i)

            statev = dataset.get_data(field="Statev", data_type="GaussPoint")[:2]

            cell_volume = density / mesh_volume

            xi_array[i + 1] = cell_volume * dataset.mesh.integrate_field(
                statev[1], "GaussPoint"
            )

            if (i - 1) % 4 == 0:  # sélectionne 1 itération sur 4
                ax = axes[k]
                k += 1
                fractions = statev[1]
                bins = np.linspace(0, 1, 101)
                ax.hist(fractions, bins=bins, edgecolor="black")

        plt.xlabel("Fraction volumique")
        plt.ylabel("Nombre d'éléments")
        plt.title("Distribution des transformations")
        plt.show()
        # np.savetxt(
        #     data_dir / f"Transformation_strain_{results_dir}.txt",
        #     transformation_strain_array,
        # )
        # strain_labels = ["et11", "et22", "et33", "et12", "et13", "et23"]

        # for k, label in enumerate(strain_labels):
        #     np.savetxt(f"{label}.txt", et_arrays[k])


def process_data_fea(typesim, cell):
    base_dir = Path("datas_simu") / cell

    results_dir = typesim
    dataset = fd.read_data(f"simuEF/{typesim}/{typesim}.fdz")
    if typesim == "shear":
        list_component = {"XY"}
    elif typesim == "tencomp":
        list_component = {"XX", "YY"}
    else:
        list_component = {"XX"}
    for component in list_component:
        mesh_volume = dataset.mesh.to_pyvista().volume
        rve_volume = dataset.mesh.bounding_box.volume
        density = mesh_volume / rve_volume

        n_iter = dataset.n_iter
        stress_array = np.zeros(n_iter + 1)
        xi_array = np.zeros(n_iter + 1)
        meanStrain_array = np.zeros(n_iter + 1)

        # transformation_strain_array = np.zeros(n_iter + 1)
        # et_arrays = np.zeros((6, n_iter + 1))

        print(n_iter)

        data_dir = base_dir / f"S{component}" / f"data_{results_dir}"
        data_dir.mkdir(parents=True, exist_ok=True)

        # list_iter = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        k = 0
        fig, axes = plt.subplots(5, 5, figsize=(10, 10))
        axes = axes.flatten()
        for i in range(n_iter):
            dataset.load(i)
            print(i)

            data_stress = dataset.get_data(
                field="Stress", component=component, data_type="GaussPoint"
            )
            vol_avg_stress = (density / mesh_volume) * dataset.mesh.integrate_field(
                field=data_stress, type_field="GaussPoint"
            )
            stress_array[i + 1] = vol_avg_stress

            data_Mstrain = dataset.get_data(
                field="MeanStrain",
                component=component,
                data_type="GaussPoint",
            )
            meanStrain_array[i + 1] = data_Mstrain[0]

            statev = dataset.get_data(field="Statev", data_type="GaussPoint")[:2]

            cell_volume = density / mesh_volume

            xi_array[i + 1] = cell_volume * dataset.mesh.integrate_field(
                statev[1], "GaussPoint"
            )

            if (i - 1) % 4 == 0:  # sélectionne 1 itération sur 4
                ax = axes[k]
                k += 1
                fractions = statev[1]
                bins = np.linspace(0, 1, 101)
                ax.hist(fractions, bins=bins, edgecolor="black")

            # vol_avg_strain = np.array(
            #     [
            #         cell_volume * dataset.mesh.integrate_field(statev[j], "GaussPoint")
            #         for j in range(2, 8)
            #     ]
            # )

            # et_arrays[:, i + 1] = vol_avg_strain

            # transformation_strain_array[i + 1] = mises_strain_fea(vol_avg_strain)

        np.savetxt(data_dir / f"Stress_{results_dir}.txt", stress_array)
        np.savetxt(data_dir / f"Xi_{results_dir}.txt", xi_array)
        np.savetxt(data_dir / f"MeanStrain_{results_dir}.txt", meanStrain_array)

        plt.xlabel("Fraction volumique")
        plt.ylabel("Nombre d'éléments")
        plt.title("Distribution des transformations")

        if not os.path.exists(f"fig_repartition/{cell}"):
            os.makedirs(f"fig_repartition/{cell}")

        plt.savefig(f"fig_repartition/{cell}/{cell}_{typesim}.png")
        # np.savetxt(
        #     data_dir / f"Transformation_strain_{results_dir}.txt",
        #     transformation_strain_array,
        # )
        # strain_labels = ["et11", "et22", "et33", "et12", "et13", "et23"]

        # for k, label in enumerate(strain_labels):
        #     np.savetxt(f"{label}.txt", et_arrays[k])


def erase_fea_file(typesim):
    file = Path(f"simuEF/{typesim}/{typesim}.fdz")
    if file.exists():
        file.unlink()


def plot_results_fea(cellule, typesim_to_loads):
    fig, axs = plt.subplots(2, 3, figsize=(10, 8))
    for i, typesim in enumerate(sorted(typesim_to_loads)):
        row = i // 3
        col = i % 3

        ax = axs[row, col]
        results_dir = typesim

        if typesim == "shear":
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
        # xi_array = np.loadtxt(
        #     f"datas_simu/{cellule}/SXX/data_{results_dir}/Xi_{results_dir}.txt"
        # )

        ax.plot(
            strain_array,
            stress_array,
            label=f"{typesim}/{cellule}",
        )
        ax.set_title(f"Plot {typesim}")
        ax.legend(loc="upper left", fontsize=8)
        ax.grid(True)
        ax.set_xlabel("E11[%]")
        ax.set_ylabel("S11 [MPa]")

    plt.tight_layout()

    plt.title(f"Plot {cellule}")
    plt.legend(loc="upper left", fontsize=8)
    plt.grid(True)
    plt.xlabel("E11[%]")
    plt.ylabel("S11 [MPa]")


def run_linear_homogenization(
    cell: str, young_modulus: float = 67538, poisson_ratio: float = 0.42
):
    linear_path = f"results_params/parametres_lineaires/E_nu_G_{cell}.txt"
    if os.path.exists(linear_path):
        values = []

        with open(linear_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    values.append(float(line))
                except ValueError:
                    # Ignore les lignes texte comme "traction-compression"
                    continue

        return np.array(values)
    else:
        print("Run linear homogeneisation")
        fd.ModelingSpace("3D")
        mesh = fd.Mesh.read(f"simuEF/cellules/{cell}.vtk")
        material = fd.constitutivelaw.ElasticIsotrop(young_modulus, poisson_ratio)
        weakform = fd.weakform.StressEquilibrium(material, nlgeom=False)
        assembly = fd.Assembly.create(weakform, mesh, mesh.elm_type, name="Assembly")

        effective_stiffness_tensor = fd.homogen.get_homogenized_stiffness(assembly)

        props_cubic = sim.L_cubic_props(effective_stiffness_tensor)
        props_cubic = props_cubic.flatten()

        with open(linear_path, "w") as f:
            f.write("cubic params\n")
            for val in props_cubic:
                f.write(f"{val:.8e}\n")

    return props_cubic


def plot_data_6D(stress_array, strain_array, time):
    e11 = strain_array[0]
    e22 = strain_array[1]
    e33 = strain_array[2]
    e12 = strain_array[3]
    e13 = strain_array[4]
    e23 = strain_array[5]
    s11 = stress_array[0]
    s22 = stress_array[1]
    s33 = stress_array[2]
    s12 = stress_array[3]
    s13 = stress_array[4]
    s23 = stress_array[5]

    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)

    # ---------- Contraintes s ----------
    axs[0, 0].plot(time, s11, color="tab:blue")
    axs[0, 0].set_title("σ11")
    axs[0, 0].set_ylabel("[MPa]")
    axs[0, 0].grid(True)

    axs[0, 1].plot(time, s22, color="tab:orange")
    axs[0, 1].set_title("σ22")
    axs[0, 1].grid(True)

    axs[0, 2].plot(time, s33, color="tab:green")
    axs[0, 2].set_title("σ33")
    axs[0, 2].grid(True)

    axs[1, 0].plot(time, s12, color="tab:red")
    axs[1, 0].set_title("σ12 ")
    axs[1, 0].set_xlabel("Temps (s)")
    axs[1, 0].set_ylabel("[MPa]")
    axs[1, 0].grid(True)

    axs[1, 1].plot(time, s13, color="tab:purple")
    axs[1, 1].set_title("σ13")
    axs[1, 1].set_xlabel("Temps (s)")
    axs[1, 1].grid(True)

    axs[1, 2].plot(time, s23, color="tab:brown")
    axs[1, 2].set_title("σ23")
    axs[1, 2].set_xlabel("Temps (s)")
    axs[1, 2].grid(True)

    fig.suptitle("Contraintes σ", fontsize=16)
    plt.tight_layout()
    plt.show()

    # ---------- Déformations e ----------’
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)

    axs[0, 0].plot(time, e11, color="tab:blue")
    axs[0, 0].set_title("ε11")
    axs[0, 0].set_ylabel("[–]")
    axs[0, 0].grid(True)

    axs[0, 1].plot(time, e22, color="tab:orange")
    axs[0, 1].set_title("ε22")
    axs[0, 1].grid(True)

    axs[0, 2].plot(time, e33, color="tab:green")
    axs[0, 2].set_title("ε33")
    axs[0, 2].grid(True)

    axs[1, 0].plot(time, e12, color="tab:red")
    axs[1, 0].set_title("ε12")
    axs[1, 0].set_xlabel("Temps (s)")
    axs[1, 0].set_ylabel("[–]")
    axs[1, 0].grid(True)

    axs[1, 1].plot(time, e13, color="tab:purple")
    axs[1, 1].set_title("ε13")
    axs[1, 1].set_xlabel("Temps (s)")
    axs[1, 1].grid(True)

    axs[1, 2].plot(time, e23, color="tab:brown")
    axs[1, 2].set_title("ε23")
    axs[1, 2].set_xlabel("Temps (s)")
    axs[1, 2].grid(True)

    fig.suptitle("Déformations ε", fontsize=16)

    # ---------- Déformations e-s ----------’
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)

    axs[0, 0].plot(e11, s11, color="tab:blue")
    axs[0, 0].set_title("ε11-σ11")
    axs[0, 0].set_ylabel("contrainte[MPa]")
    axs[0, 0].grid(True)

    axs[0, 1].plot(e22, s22, color="tab:orange")
    axs[0, 1].set_title("ε22-σ22")
    axs[0, 1].set_ylabel("contrainte[MPa]")
    axs[0, 1].grid(True)

    axs[0, 2].plot(e33, s33, color="tab:green")
    axs[0, 2].set_title("ε33-σ33")
    axs[0, 2].set_ylabel("contrainte[MPa]")
    axs[0, 2].grid(True)

    axs[1, 0].plot(e12, s12, color="tab:red")
    axs[1, 0].set_title("ε12-σ12")
    axs[1, 0].set_xlabel("deformation[-]")
    axs[1, 0].set_ylabel("contrainte[MPa]")
    axs[1, 0].grid(True)

    axs[1, 1].plot(e13, s13, color="tab:purple")
    axs[1, 1].set_title("ε13-σ13")
    axs[1, 1].set_xlabel("deformation[-]")
    axs[1, 1].grid(True)

    axs[1, 2].plot(e23, s23, color="tab:brown")
    axs[1, 2].set_title("ε23-σ23")
    axs[1, 2].set_xlabel("deformation[-]")
    axs[1, 2].grid(True)

    fig.suptitle(
        "Courbes contrainte-déformation sur un élément (aléatoire)", fontsize=16
    )
    plt.tight_layout()
    plt.show()

    return time, stress_array, strain_array


def save_data_csv(
    mean_stress_array,
    mean_strain_array,
    time,
    sim_id,
    csv_file="simuEF/train_fea.csv",
):
    df = pd.DataFrame(
        {
            "total_strain_xx": mean_strain_array[0],
            "total_strain_yy": mean_strain_array[1],
            "total_strain_zz": mean_strain_array[2],
            "total_strain_xy": mean_strain_array[3],
            "total_strain_xz": mean_strain_array[4],
            "total_strain_yz": mean_strain_array[5],
            "stress_xx": mean_stress_array[0],
            "stress_yy": mean_stress_array[1],
            "stress_zz": mean_stress_array[2],
            "stress_xy": mean_stress_array[3],
            "stress_xz": mean_stress_array[4],
            "stress_yz": mean_stress_array[5],
            "timestep": time,
            "simulation_load_id": sim_id,
        }
    )

    if not os.path.isfile(csv_file):
        df.to_csv(csv_file, index=False)
    else:
        # append sans réécrire les colonnes
        df.to_csv(csv_file, mode="a", header=False, index=False)
