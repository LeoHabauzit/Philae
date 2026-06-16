import pandas as pd
from pathlib import Path

filename = "lstm/dataset/train_fea_dataset_80.csv"


def swap_cols(df, case="x-y"):
    df2 = df.copy()
    if case == "x-y":
        swaps = [
            ("total_strain_xx", "total_strain_yy"),
            ("total_strain_xz", "total_strain_yz"),
            ("stress_xx", "stress_yy"),
            ("stress_xz", "stress_yz"),
        ]
    elif case == "x-z":
        swaps = [
            ("total_strain_xx", "total_strain_zz"),
            ("total_strain_xy", "total_strain_yz"),
            ("stress_xx", "stress_zz"),
            ("stress_xy", "stress_yz"),
        ]

    elif case == "y-z":
        swaps = [
            ("total_strain_yy", "total_strain_zz"),
            ("total_strain_xy", "total_strain_xz"),
            ("stress_yy", "stress_zz"),
            ("stress_xy", "stress_xz"),
        ]
    for a, b in swaps:
        df2[a], df2[b] = df2[b].copy(), df2[a].copy()

    return df2


def data_augmentation_symmetry(
    filename,
):

    name = Path(filename).with_suffix("")
    components = {"x-y", "x-z", "y-z"}
    df = pd.read_csv(filename)
    df.to_csv(f"{name}_augmented.csv", mode="a", header=True, index=False)
    for i in range(0, len(df), 101):
        bloc = df.iloc[i : i + 101]
        for component in components:
            df2 = swap_cols(bloc, component)
            df2.to_csv(f"{name}_augmented.csv", mode="a", header=False, index=False)


data_augmentation_symmetry(filename)
