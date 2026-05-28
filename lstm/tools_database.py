import os
import numpy as np
import pandas as pd
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from Umat.loi_sma import umat_sma_random


def generer_path_txt(
    filename,
    n_steps=2,
    Dn_inc=0.001,
    temperature_initiale=300,
):
    stress = np.random.uniform(-1000, 1000, (1, 6))

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

    print(f"Fichier '{filename}' généré avec succès.")


def add_zeros_to_arrays(variables):
    return [np.insert(v, 0, 0) for v in variables]


def generate_data_csv(props, n_simulations):
    os.remove("lstm/train_dataset.csv") if os.path.exists(
        "lstm/train_dataset.csv"
    ) else None
    csv_file = "lstm/train_dataset.csv"
    sim_ids = np.random.choice(np.arange(1, 10001), size=n_simulations, replace=False)
    for k in range(n_simulations):
        generer_path_txt(
            filename="Umat/data/path_random.txt", temperature_initiale=300, Dn_inc=0.02
        )
        umat_sma_random(props, "SMADI")
        os.remove("Umat/data/path_random.txt")
        outputfile_global = "Umat/results_SMADI/results_random_global-0.txt"

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
                "e11": e11,
                "e22": e22,
                "e12": e12,
                "s11": s11,
                "s22": s22,
                "s12": s12,
                "timestep": time,
                "simulation_id": sim_id,
            }
        )

        if not os.path.isfile(csv_file):
            df.to_csv(csv_file, index=False)
        else:
            # append sans réécrire les colonnes
            df.to_csv(csv_file, mode="a", header=False, index=False)

    print("Toutes les simulations ont été ajoutées.")
