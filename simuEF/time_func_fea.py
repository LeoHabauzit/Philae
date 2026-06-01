import fedoo as fd
import numpy as np

fd.ModelingSpace("3D")

mesh = fd.mesh.box_mesh(
    nx=10,
    ny=5,
    nz=5,
    x_min=0,
    x_max=100,
    y_min=0,
    y_max=10,
    z_min=0,
    z_max=10,
    elm_type="hex8",
    name="Domain",
)

fd.constitutivelaw.ElasticIsotrop(200e3, 0.3, name="Steel")
wf = fd.weakform.StressEquilibrium("Steel", nlgeom=True)
assembly = fd.Assembly.create(wf, mesh, name="Assembly")
pb = fd.problem.NonLinear("Assembly")

nodes_left = mesh.find_nodes("X", mesh.bounding_box.xmin)
nodes_right = mesh.find_nodes("X", mesh.bounding_box.xmax)

# ── Condition bloquée (sans time_func) ────────────────────────────────────
pb.bc.add("Dirichlet", nodes_left, "Disp", 0)

# ── Évolution LINÉAIRE (comportement par défaut si time_func=None) ─────────
pb.bc.add("Dirichlet", nodes_right, "DispX", 5.0)
# équivalent explicite :
# pb.bc.add("Dirichlet", nodes_right, "DispX", 5.0, time_func=lambda t: t)


# ── Évolution en MARCHE (step function) ───────────────────────────────────
# def step_func(t_fact):
#     """Applique immédiatement la valeur complète dès le début du pas."""
#     if t_fact == 0:
#         return 0.0
#     return 1.0


# pb.bc.add("Dirichlet", nodes_right, "DispY", 2.0, time_func=step_func)


# ── Évolution SINUSOÏDALE (ex: chargement cyclique) ───────────────────────
# def sinus_func(t_fact):
#     """Suit sin(π/2 * t) : part de 0, atteint 1 à t=1 (quart de cycle)."""
#     return np.sin(np.pi / 2 * t_fact)


# pb.bc.add("Dirichlet", nodes_right, "DispZ", 1.0, time_func=sinus_func)


# ── Évolution QUADRATIQUE (rampe douce, accélération progressive) ─────────
def quad_func(t_fact):
    """Évolution en t² : démarrage doux, arrivée rapide."""
    return t_fact**2


# ── start_value : valeur initiale explicite (sinon = état courant) ─────────
# Utile pour enchaîner plusieurs nlsolve avec des valeurs différentes
pb.bc.add(
    "Dirichlet", nodes_right, "DispX", 10.0, time_func=quad_func, start_value=0.0
)  # forcer le départ à 0 même si état courant != 0

# ── Résolution non-linéaire avec gestion du temps ─────────────────────────
pb.set_nr_criterion("Displacement", tol=1e-3, max_subiter=10)

# t_fact évolue de 0 à 1 à l'intérieur de chaque appel nlsolve
# (décomposé en sous-incréments via dt)
pb.nlsolve(dt=0.1, tmax=1.0, update_dt=True, print_info=2)

# ── Enchaîner un deuxième pas (déchargement) ──────────────────────────────
# pb.bc.remove(-1)  # retire la dernière BC
# pb.bc.add(
#     "Dirichlet",
#     nodes_right,
#     "DispX",
#     0.0,
#     time_func=lambda t: t,  # retour linéaire vers 0
#     start_value=None,
# )  # repart de l'état courant (fin du pas précédent)

# pb.nlsolve(dt=0.1, tmax=1.0, update_dt=True, print_info=1)
