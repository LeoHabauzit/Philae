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
):
    stress = np.random.uniform(-stress_lim, stress_lim, (1, 6))

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


def generate_data_csv(props, n_simulations, umat_name="SMAAC"):
    os.remove("lstm/train_dataset.csv") if os.path.exists(
        "lstm/train_dataset.csv"
    ) else None
    csv_file = "lstm/train_dataset.csv"
    sim_ids = np.random.choice(np.arange(1, 10001), size=n_simulations, replace=False)
    for k in range(n_simulations):
        generer_path_txt(
            filename="Umat/data/path_random.txt",
            temperature_initiale=300,
            Dn_inc=0.02,
            stress_lim=200,
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

    e11 = df["e11"].values[j : j + 102]
    e22 = df["e22"].values[j : j + 102]
    e33 = df["e33"].values[j : j + 102]
    e12 = df["e12"].values[j : j + 102]
    e13 = df["e13"].values[j : j + 102]
    e23 = df["e23"].values[j : j + 102]
    s11 = df["s11"].values[j : j + 102]
    s22 = df["s22"].values[j : j + 102]
    s33 = df["s33"].values[j : j + 102]
    s12 = df["s12"].values[j : j + 102]
    s13 = df["s13"].values[j : j + 102]
    s23 = df["s23"].values[j : j + 102]
    time = df["timestep"].values[j : j + 102]

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
    plt.tight_layout()
    plt.show()
