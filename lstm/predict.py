"""
predict.py
Charge le modèle entraîné (model.pth) et affiche la prédiction
pour une simulation choisie du jeu de test.
"""

from pathlib import Path

import torch

import example_utils
from lstm_example import RNNModel

# ─────────────────────────────────────────────
# PARAMÈTRES À MODIFIER SI BESOIN
# ─────────────────────────────────────────────
TRAIN_CSV = Path("lstm") / "dataset" / "train_dataset.csv"
# TEST_CSV = Path("lstm") / "dataset" / "test_dataset.csv"
TEST_CSV = Path("lstm") / "dataset" / "test_fea.csv"
MODEL_PTH = "model_2000.pth"

# Index de la simulation à afficher (0 = première simulation du jeu de test)
SAMPLE_INDEX = 0

# Architecture : doit être IDENTIQUE à celle utilisée à l'entraînement
INPUT_SIZE = 6
HIDDEN_SIZE = 64
OUTPUT_SIZE = 6
NUM_LAYERS = 2

COMPONENT_LABELS = ["xx", "yy", "zz", "xy", "xz", "zy"]
# ─────────────────────────────────────────────


def main() -> None:
    # 1. Device
    if torch.backends.mps.is_available():
        device = "mps"
    elif torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"
    print(f"Using device = {device}")

    # 2. Chargement des données
    #    On a besoin du train pour recalculer les stats de normalisation
    print("Chargement des données...")
    x_train, y_train, _ = example_utils.load_data(TRAIN_CSV.as_posix())
    x_test, y_test, sim_ids_test = example_utils.load_data(TEST_CSV.as_posix())

    # 3. Statistiques de normalisation (calculées sur le train, comme à l'entraînement)
    x_mean, x_std = x_train.mean(), x_train.std()
    y_mean, y_std = y_train.mean(), y_train.std()

    # 4. Normalisation du jeu de test
    x_test_norm = (x_test - x_mean) / x_std
    y_test_norm = (y_test - y_mean) / y_std

    # 5. Chargement du modèle
    model = RNNModel(
        input_features_size=INPUT_SIZE,
        hidden_state_size=HIDDEN_SIZE,
        output_size=OUTPUT_SIZE,
        num_layers=NUM_LAYERS,
    )
    # model.load_state_dict(torch.load(MODEL_PTH, map_location=device, weights_only=True))
    checkpoint = torch.load(MODEL_PTH, map_location=device, weights_only=True)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    model.eval()
    print(f"Modèle chargé depuis '{MODEL_PTH}'")

    # 6. Prédiction pour la simulation choisie
    strain_norm = x_test_norm[SAMPLE_INDEX].unsqueeze(0).to(device)  # (1, T, 6)
    true_stress_norm = y_test_norm[SAMPLE_INDEX]  # (T, 6)

    with torch.no_grad():
        pred_stress_norm, _ = model(strain_norm)

    pred_stress_norm = pred_stress_norm.squeeze(0).cpu()  # (T, 6)

    # 7. Dénormalisation
    strain_sample = x_test[SAMPLE_INDEX]  # déjà en valeurs réelles
    true_stress = true_stress_norm * y_std + y_mean
    pred_stress = pred_stress_norm * y_std + y_mean

    sim_id = sim_ids_test[SAMPLE_INDEX].item()
    print(f"Affichage de la simulation id={sim_id} (index {SAMPLE_INDEX})")

    # 8. Graphique
    example_utils.plot_stress_strain_sample_with_prediction(
        strain_sequence=strain_sample,
        true_stress_sequence=true_stress,
        predicted_stress_sequence=pred_stress,
        sim_id=int(sim_id),
        component_labels=COMPONENT_LABELS,
    )


if __name__ == "__main__":
    main()
