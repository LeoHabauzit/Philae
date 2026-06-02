import fedoo as fd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab

fd.ModelingSpace("3D")

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
all_components = {"XX", "YY", "ZZ", "XY", "XZ", "YZ"}
for k in range(dataset.n_iter):
    dataset.load(k)
    for i, component in enumerate(sorted(all_components)):
        print(i, component)
        data_stress = dataset.get_data(
            field="Stress", component="XX", data_type="GaussPoint"
        )
        data_strain = dataset.get_data(
            field="Strain", component="XX", data_type="GaussPoint"
        )
        stress_array[i, k + 1] = data_stress[0]
        strain_array[i, k + 1] = data_strain[0]

    # stress_array[k + 1] = data_stress[0]
    # strain_array[k + 1] = data_strain[0]
    print(data_stress[0])
# plt.plot(time, stress_array)
plt.plot(time, strain_array)
plt.show()
# fd.viewer(res)
