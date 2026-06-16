import os
import numpy as np
import pandas as pd
import sys
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent.parent))

from Umat.loi_sma import umat_sma_random


def generer_path_txt(
    filename,
    n_steps=2,
    Dn_inc=0.001,
    temperature_initiale=300,
    stress_lim=1000,
    stress_target=None,
):
    if stress_target is None:
        stress = np.random.uniform(-stress_lim, stress_lim, (1, 6))
    else:
        stress = stress_target

    with open(filename, "w", encoding="utf-8") as f:
        f.write("#Initial_temperature\n")
        f.write(f"{temperature_initiale}\n")

        f.write("#Number_of_blocks\n")
        f.write("1\n\n")

        f.write("#Block\n")
        f.write("1\n")

        f.write("#Loading_type\n")
        f.write("1\n")

        f.write("#Control_type(NLGEOM)\n")
        f.write("1\n")

        f.write("#Repeat\n")
        f.write("1\n")

        f.write("#Steps\n")
        f.write(f"{n_steps}\n\n")

        f.write("#Mode\n")
        f.write("1\n")

        f.write("#Dn_init 1.\n")
        f.write("#Dn_mini 1.\n")
        f.write(f"#Dn_inc {Dn_inc}\n")

        f.write("#time\n")
        f.write("50\n")

        f.write("#Consigne\n")

        # Valeurs des contraintes
        f.write(f"S {stress[0, 0]}\n")
        f.write(f"S {stress[0, 3]} S {stress[0, 1]}\n")
        f.write(f"S {stress[0, 4]} S {stress[0, 5]} S {stress[0, 2]}\n")

        f.write("#Consigne_T\n")
        f.write(f"T {temperature_initiale}\n\n")

        f.write("#Mode\n")
        f.write("1\n")

        f.write("#Dn_init 1.\n")
        f.write("#Dn_mini 1.\n")
        f.write(f"#Dn_inc {Dn_inc}\n")

        f.write("#time\n")
        f.write("50\n")

        f.write("#Consigne\n")

        # Valeurs des contraintes
        f.write("S 0\n")
        f.write("S 0 S 0\n")
        f.write("S 0 S 0 S 0\n")

        f.write("#Consigne_T\n")
        f.write(f"T {temperature_initiale}\n\n")

    # print(f"Fichier '{filename}' généré avec succès.")


def add_zeros_to_arrays(variables):
    return [np.insert(v, 0, 0) for v in variables]


def generate_data_csv(
    props,
    n_simulations,
    umat_name="SMAAC",
    csv_file="lstm/dataset/train_dataset.csv",
    stress_lim=150,
):
    os.remove(csv_file) if os.path.exists(csv_file) else None

    sim_ids = np.random.choice(np.arange(1, 10001), size=n_simulations, replace=False)
    for k in range(n_simulations):
        print(k)
        generer_path_txt(
            filename="Umat/data/path_random.txt",
            temperature_initiale=300,
            Dn_inc=0.02,
            stress_lim=stress_lim,
        )
        umat_sma_random(props, umat_name)
        os.remove("Umat/data/path_random.txt")
        outputfile_global = f"Umat/results_{umat_name}/results_random_global-0.txt"

        time, e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23 = np.loadtxt(
            outputfile_global,
            usecols=(4, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19),
            unpack=True,
        )
        sim_id = np.full(len(time), sim_ids[k])
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

    print("Toutes les simulations ont été ajoutées.")


def read_data(filename, i=0):
    j = 101 * i
    # e11, e22, e12, s11, s22, s12, time, sim_id = np.loadtxt(
    #     filename,
    #     delimiter=",",
    #     skiprows=1,
    #     max_rows=102,
    #     usecols=(0, 1, 2),
    #     unpack=True,
    # )
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
    time = df["timestep"].values[j : j + 101]

    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)

    # ---------- Contraintes s ----------
    axs[0, 0].plot(time, s11, color="tab:blue")
    axs[0, 0].set_title("σ11")
    axs[0, 0].set_ylabel("σ [MPa]")
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
    axs[1, 0].set_ylabel("σ [MPa]")
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
    axs[0, 0].set_ylabel("ε [–]")
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
    axs[1, 0].set_ylabel("ε [–]")
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
    axs[0, 0].set_ylabel("σ [MPa]")
    axs[0, 0].grid(True)

    axs[0, 1].plot(e22, s22, color="tab:orange")
    axs[0, 1].set_title("ε22-σ22")
    axs[0, 1].set_ylabel("σ [MPa]")
    axs[0, 1].grid(True)

    axs[0, 2].plot(e33, s33, color="tab:green")
    axs[0, 2].set_title("ε33-σ33")
    axs[0, 2].set_ylabel("σ [MPa]")
    axs[0, 2].grid(True)

    axs[1, 0].plot(e12, s12, color="tab:red")
    axs[1, 0].set_xlabel("ε[-]")
    axs[1, 0].set_ylabel("σ [MPa]")
    axs[1, 0].grid(True)

    axs[1, 1].plot(e13, s13, color="tab:purple")
    axs[1, 1].set_title("ε13-σ13")
    axs[1, 1].set_xlabel("ε[-]")
    axs[1, 1].grid(True)

    axs[1, 2].plot(e23, s23, color="tab:brown")
    axs[1, 2].set_title("ε23-σ23")
    axs[1, 2].set_xlabel("ε[-]")
    axs[1, 2].grid(True)

    fig.suptitle("Courbes contrainte-déformation", fontsize=16)
    plt.tight_layout()
    plt.show()


def read_multiple_data(dict_file=None, i=0):
    """
    dict_file={
    "simulation_ref.csv": {
        "color": "black",
        "linestyle": "-",
        "label": "Référence"
    },
    "simulation_lstm.csv": {
        "color": "red",
        "linestyle": "--",
        "label": "LSTM"
    }
    }
    """
    j = 101 * i
    components = [
        ("xx", 0, 0),
        ("yy", 0, 1),
        ("zz", 0, 2),
        ("xy", 1, 0),
        ("xz", 1, 1),
        ("yz", 1, 2),
    ]

    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for filename, params in dict_file.items():
        df = pd.read_csv(filename)

        for comp, row, col in components:
            strain = df[f"total_strain_{comp}"].values[j : j + 101]
            stress = df[f"stress_{comp}"].values[j : j + 101]

            ax = axs[row, col]

            ax.plot(
                strain,
                stress,
                color=params["color"],
                linestyle=params["linestyle"],
                label=params["label"],
            )

            ax.set_title(f"{comp}")
            ax.set_xlabel("ε [-]")
            ax.set_ylabel("σ [MPa]")
            ax.grid(True)
            ax.legend()

        fig.suptitle("Courbes contrainte-déformation", fontsize=16)
        plt.tight_layout()

    plt.show()


def split_csv(csv_file, n_train, block_size=101):

    df = pd.read_csv(csv_file)

    split_idx = n_train * block_size

    train_df = df.iloc[:split_idx]
    test_df = df.iloc[split_idx:]
    os.remove(
        f"lstm/dataset/train_fea_dataset_{int(len(train_df) / block_size)}.csv"
    ) if os.path.exists(
        f"lstm/dataset/train_fea_dataset_{int(len(train_df) / block_size)}.csv"
    ) else None
    os.remove(
        f"lstm/dataset/test_fea_dataset_{int(len(test_df) / block_size)}.csv"
    ) if os.path.exists(
        f"lstm/dataset/test_fea_dataset_{int(len(test_df) / block_size)}.csv"
    ) else None

    train_df.to_csv(
        f"lstm/dataset/train_fea_dataset_{int(len(train_df) / block_size)}.csv",
        index=False,
    )
    test_df.to_csv(
        f"lstm/dataset/test_fea_dataset_{int(len(test_df) / block_size)}.csv",
        index=False,
    )

    print(f"Train : {len(train_df) // block_size} blocs")
    print(f"Test  : {len(test_df) // block_size} blocs")


def save_data_csv(
    mean_stress_array,
    mean_strain_array,
    time=None,
    sim_id=None,
    csv_file="simuEF/csv_files/default_csv_file.csv",
):
    os.remove(csv_file) if os.path.exists(csv_file) else None
    if time is None:
        time = np.linspace(0, 1, len(mean_strain_array))

    if sim_id is None:
        sim_id = np.full(len(time), 1)
    print(mean_strain_array[0, 0])
    if mean_strain_array[0, 0] != 0:
        add_zeros_to_arrays(
            [
                mean_strain_array[:, 0],
                mean_strain_array[:, 1],
                mean_strain_array[:, 2],
                mean_strain_array[:, 3],
                mean_strain_array[:, 4],
                mean_strain_array[:, 5],
                mean_stress_array[:, 0],
                mean_stress_array[:, 1],
                mean_stress_array[:, 2],
                mean_stress_array[:, 3],
                mean_stress_array[:, 4],
                mean_stress_array[:, 5],
                time,
                sim_id,
            ]
        )
    df = pd.DataFrame(
        {
            "total_strain_xx": mean_strain_array[:, 0],
            "total_strain_yy": mean_strain_array[:, 1],
            "total_strain_zz": mean_strain_array[:, 2],
            "total_strain_xy": mean_strain_array[:, 3],
            "total_strain_xz": mean_strain_array[:, 4],
            "total_strain_yz": mean_strain_array[:, 5],
            "stress_xx": mean_stress_array[:, 0],
            "stress_yy": mean_stress_array[:, 1],
            "stress_zz": mean_stress_array[:, 2],
            "stress_xy": mean_stress_array[:, 3],
            "stress_xz": mean_stress_array[:, 4],
            "stress_yz": mean_stress_array[:, 5],
            "timestep": time,
            "simulation_load_id": sim_id,
        }
    )

    if not os.path.isfile(csv_file):
        df.to_csv(csv_file, index=False)
    else:
        # append sans réécrire les colonnes

        df.to_csv(csv_file, mode="a", header=False, index=False)
