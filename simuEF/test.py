import fedoo as fd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab

fd.ModelingSpace("3D")


def plot_data(stress_array, strain_array, time):
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
    plt.tight_layout()
    plt.show()

    return time, stress_array, strain_array


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

stress_lim = 150
stress_target = np.random.uniform(-stress_lim, stress_lim, (1, 6))
stress_zeros = np.zeros((1, 6))
print(stress_target)
# stress_target = np.zeros((1, 6))
# stress_target[0, 1] = 100
# print(stress_target)

E = 200
nu = 0.3

material = fd.constitutivelaw.ElasticIsotrop(E, nu, name="Material")
wf = fd.weakform.StressEquilibrium(material)
assembly = fd.Assembly.create(wf, mesh, name="Assembly")

pb = fd.problem.NonLinear(assembly)

pb.bc.add(fd.constraint.PeriodicBC(periodicity_type="small_strain"))

print(mesh.bounding_box)
ref_node = mesh.nearest_node(mesh.bounding_box.center)
print(ref_node)

pb.bc.add("Dirichlet", 0, "Disp", 0)

volume = mesh.bounding_box.volume  # = 1.0 mm³ for the unit cube

pb.bc.add("Neumann", "E_xx", stress_target[0, 0] * volume)
pb.bc.add("Neumann", "E_yy", stress_target[0, 1] * volume)
pb.bc.add("Neumann", "E_zz", stress_target[0, 2] * volume)
pb.bc.add("Neumann", "E_xy", stress_target[0, 3] * volume)  # conjugate of γ_xy is σ_xy
pb.bc.add("Neumann", "E_xz", stress_target[0, 4] * volume)
pb.bc.add("Neumann", "E_yz", stress_target[0, 5] * volume)
print(pb.bc)
results = pb.add_output(
    "test",
    assembly,
    ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
)
pb.nlsolve(dt=0.2, tmax=1, update_dt=True, print_info=1, interval_output=0.02)

for k in range(6):
    pb.bc.remove(-1)
pb.bc.add("Neumann", "E_xx", 0)
pb.bc.add("Neumann", "E_yy", 0)
pb.bc.add("Neumann", "E_zz", 0)
pb.bc.add("Neumann", "E_xy", 0)  # conjugate of γ_xy is σ_xy
pb.bc.add("Neumann", "E_xz", 0)
pb.bc.add("Neumann", "E_yz", 0)

pb.nlsolve(dt=0.2, tmax=1, update_dt=True, print_info=1, interval_output=0.02)
res = pb.get_results(
    # "test",
    assembly,
    ["Stress", "Strain", "Disp", "MeanStrain", "Fext(MeanStrain)"],
)

dataset = fd.read_data("test.fdz")
n_iter = dataset.n_iter
time = np.linspace(0, 1, 101)
# dataset.load(-1)
# print(dataset.n_iter)
stress_array = np.zeros((6, n_iter + 1))
print(stress_array.shape)
strain_array = np.zeros((6, n_iter + 1))
all_components = ["XX", "YY", "ZZ", "XY", "XZ", "YZ"]
for k in range(dataset.n_iter):
    dataset.load(k)
    for i, component in enumerate(all_components, start=0):
        data_stress = dataset.get_data(
            field="Stress", component=component, data_type="GaussPoint"
        )
        data_strain = dataset.get_data(
            field="Strain", component=component, data_type="GaussPoint"
        )
        stress_array[i, k + 1] = data_stress[i]
        strain_array[i, k + 1] = data_strain[i]

# plt.plot(time, stress_array)
plot_data(stress_array, strain_array, time)
# fd.viewer(res)
