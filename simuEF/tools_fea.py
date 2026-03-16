import numpy as np
import fedoo as fd
import os
from pathlib import Path
import matplotlib.pyplot as plt


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

        transformation_strain_array = np.zeros(n_iter + 1)
        et_arrays = np.zeros((6, n_iter + 1))

        print(n_iter)

        data_dir = base_dir / f"S{component}" / f"data_{results_dir}"
        data_dir.mkdir(parents=True, exist_ok=True)

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

            statev = dataset.get_data(field="Statev", data_type="GaussPoint")[:8]

            cell_volume = density / mesh_volume

            xi_array[i + 1] = cell_volume * dataset.mesh.integrate_field(
                statev[1], "GaussPoint"
            )

            vol_avg_strain = np.array(
                [
                    cell_volume * dataset.mesh.integrate_field(statev[j], "GaussPoint")
                    for j in range(2, 8)
                ]
            )

            et_arrays[:, i + 1] = vol_avg_strain

            transformation_strain_array[i + 1] = mises_strain_fea(vol_avg_strain)

        np.savetxt(data_dir / f"Stress_{results_dir}.txt", stress_array)
        np.savetxt(data_dir / f"Xi_{results_dir}.txt", xi_array)
        np.savetxt(data_dir / f"MeanStrain_{results_dir}.txt", meanStrain_array)
        np.savetxt(
            data_dir / f"Transformation_strain_{results_dir}.txt",
            transformation_strain_array,
        )
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
