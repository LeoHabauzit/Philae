from tools_homogeneisation import *
from Umat.loi_sma import umat_sma

cell = "Cuboctahedron40"
props_var = load_variable_props(f"results_params/params_smani_{cell}.txt")
finalprops = vect_props_smadi(props_var)
print(finalprops)

data_simu_dir = f"simuEF/datas_simu/{cell}"
typesim_to_loads = {
    "tension",
    "biaxial_tension",
    "compression",
    "biaxial_compression",
    "tencomp",
    "shear",
}
fig, axes_strain = plt.subplots(2, 3, figsize=(10, 8))
for i, typesim in enumerate(sorted(typesim_to_loads)):
    losses = []
    row = i // 3
    col = i % 3

    ax = axes_strain[row, col]
    results_dir = typesim
    umat_sma(finalprops, typesim, "SMADI")

    outputfile_global = f"Umat/results_SMADI/results_{typesim}_global-0.txt"

    e11, e22, e33, e12, e13, e23, s11, s22, s33, s12, s13, s23, xi = np.loadtxt(
        outputfile_global,
        usecols=(8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 25),
        unpack=True,
    )

    if typesim == "shear":
        strain_num = e12
        stress_num = s12
        stress_exp = np.loadtxt(
            f"{data_simu_dir}/SXY/data_{results_dir}/Stress_{results_dir}.txt"
        )
        strain_exp = np.loadtxt(
            f"{data_simu_dir}/SXY/data_{results_dir}/MeanStrain_{results_dir}.txt"
        )
    else:
        strain_num = e11
        stress_num = s11
        stress_exp = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/Stress_{results_dir}.txt"
        )
        strain_exp = np.loadtxt(
            f"{data_simu_dir}/SXX/data_{results_dir}/MeanStrain_{results_dir}.txt"
        )

    ax.plot(
        strain_exp,
        stress_exp,
        label=typesim,
    )
    ax.plot(strain_num, stress_num, c="red", label="UMAT SMA")
    ax.set_xlabel("E11 [%]")
    ax.set_ylabel("S11 [MPa]")
    ax.grid()
    ax.legend()
plt.suptitle(f"{cell}")

plt.tight_layout()

# plt.title(f"Plot {typesim}")
plt.legend(loc="upper left", fontsize=8)
plt.grid(True)
plt.show()
