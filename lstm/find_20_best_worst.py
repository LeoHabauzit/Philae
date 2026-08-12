"""
find_best_worst.py
Charge le modèle entraîné (model.pth), prédit sur TOUTES les simulations
du jeu de test, calcule l'erreur (MSE) par simulation, puis affiche :
  - la MEILLEURE prédiction (MSE la plus faible)
  - la PIRE prédiction     (MSE la plus élevée)
"""

from pathlib import Path

import example_utils
import matplotlib.pyplot as plt
import numpy as np
import torch
from first_training import RNNModel
from matplotlib.backends.backend_pdf import PdfPages
from torch import randn
from torchmetrics.regression import WeightedMeanAbsolutePercentageError


def wmape_verif(y_true, y_pred):
    return torch.sum(torch.abs(y_true - y_pred)) / torch.sum(torch.abs(y_true))


# ─────────────────────────────────────────────
# PARAMÈTRES À MODIFIER SI BESOIN
# ─────────────────────────────────────────────

shuffle_id = 6329

TRAIN_CSV = Path("lstm") / "dataset" / "train_dataset.csv"
# TEST_CSV = Path("lstm") / "dataset" / "test_dataset.csv"
TEST_CSV = (
    Path("lstm")
    / "dataset"
    / f"datasets_{shuffle_id}"
    / f"test_shuffled_{shuffle_id}.csv"
)

# TEST_CSV = Path("lstm") / "dataset" / "test_fea_dataset_23_augmented.csv"

MODEL_PTH = f"lstm/models_cuboctahedron40/model_finetuned_test_{shuffle_id}.pth"
# MODEL_PTH = "lstm/models_cuboctahedron40/model_5000_HS32.pth"

# Architecture : doit être IDENTIQUE à celle utilisée à l'entraînement
INPUT_SIZE = 6
HIDDEN_SIZE = 32
OUTPUT_SIZE = 6
NUM_LAYERS = 2

COMPONENT_LABELS = ["xx", "yy", "zz", "xy", "xz", "zy"]

# ─────────────────────────────────────────────


def predict_all(model, x_norm, device):
    """Prédit toutes les simulations, renvoie un tenseur (N, T, 6)."""
    model.eval()
    all_preds = []
    with torch.no_grad():
        for i in range(x_norm.shape[0]):
            strain = x_norm[i].unsqueeze(0).to(device)  # (1, T, 6)
            pred, _ = model(strain)
            all_preds.append(pred.squeeze(0).cpu())  # (T, 6)
    return torch.stack(all_preds)  # (N, T, 6)


def mse_per_simulation(preds, targets):
    """Calcule le MSE pour chaque simulation. Retourne un tenseur (N,)."""
    # preds et targets : (N, T, 6)
    diff = preds - targets  # (N, T, 6)
    mse = (diff**2).mean(dim=[1, 2])  # (N,)  une valeur par simulation
    return mse


def save_simulation_to_pdf(
    pdf,
    idx,
    sim_ids_test,
    x_test,
    y_test_norm,
    preds_norm,
    y_std,
    y_mean,
    mse,
    category_label,
):
    """Génère le graphique d'une simulation et l'enregistre dans le PDF."""
    sim_id = sim_ids_test[idx].item()

    strain_sample = x_test[idx]  # valeurs réelles
    true_stress = y_test_norm[idx] * y_std + y_mean  # dénormalisé
    pred_stress = preds_norm[idx] * y_std + y_mean  # dénormalisé

    title = f"[{category_label}] idx={idx} – WMAPE={mse[idx]:.6f}"
    print(f"Enregistrement : {title} (sim_id={sim_id})")

    example_utils.plot_stress_strain_sample_with_prediction(
        strain_sequence=strain_sample,
        true_stress_sequence=true_stress,
        predicted_stress_sequence=pred_stress,
        sim_id=int(sim_id),
        component_labels=COMPONENT_LABELS,
        title=title,
        pdf_file=pdf,  # on ne sauvegarde pas séparément, on utilise PdfPages
    )

    # pdf.savefig(plt.gcf())
    plt.close()


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
    print("Chargement des données...")
    x_train, y_train, _ = example_utils.load_data(TRAIN_CSV.as_posix())
    x_test, y_test, sim_ids_test = example_utils.load_data(TEST_CSV.as_posix())

    # 3. Statistiques de normalisation (sur le train)
    x_mean, x_std = x_train.mean(), x_train.std()
    y_mean, y_std = y_train.mean(), y_train.std()

    # 4. Normalisation
    x_test_norm = (x_test - x_mean) / x_std
    y_test_norm = (y_test - y_mean) / y_std

    # 5. Chargement du modèle
    model = RNNModel(
        input_features_size=INPUT_SIZE,
        hidden_state_size=HIDDEN_SIZE,
        output_size=OUTPUT_SIZE,
        num_layers=NUM_LAYERS,
    )
    checkpoint = torch.load(MODEL_PTH, map_location=device, weights_only=True)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    print(f"Modèle chargé depuis '{MODEL_PTH}'")

    # 6. Prédictions sur tout le jeu de test
    print(f"Prédiction sur {x_test.shape[0]} simulations...")
    preds_norm = predict_all(model, x_test_norm, device)  # (N, T, 6) normalisé
    print(preds_norm.shape[0], x_test_norm.shape)
    true_stress = y_test_norm * y_std + y_mean  # dénormalisé
    pred_stress = preds_norm * y_std + y_mean  # dénormalisé

    wmape = WeightedMeanAbsolutePercentageError()
    values = []
    plt.figure(figsize=(10, 6))
    for k in range(pred_stress.shape[0]):
        # print(pred_stress[k].shape, true_stress[k].shape)
        # print(values)
        values.append(wmape(pred_stress[k], true_stress[k]))
        # wmape.update(preds_norm[k], x_test_norm[k])
    values = np.array(values)
    wmape_val = torch.tensor(values)
    bins = np.linspace(np.min(values), np.max(values), 50)
    plt.hist(values, bins=bins, edgecolor="black")
    plt.xlabel("Erreur (WMAPE)")
    plt.ylabel("Prédictions")
    plt.title(
        f"Distribution de l'erreur WMAPE sur le jeu de test shuffle {shuffle_id} - n=60"
    )
    # wmape(preds_norm, x_test_norm)
    # fig_, ax_ = wmape.plot(values)
    # 7. Calcul du MSE par simulation (en espace normalisé)
    mse = mse_per_simulation(preds_norm, y_test_norm)  # (N,)

    sorted_idx = torch.argsort(wmape_val)
    best_idx = sorted_idx[:10].tolist()
    worst_idx = sorted_idx[-10:].tolist()
    worst_idx.reverse()
    mse = np.array(mse)
    # print(mse)
    plt.figure(figsize=(10, 6))
    bins = np.linspace(np.min(mse), np.max(mse), 50)
    plt.hist(mse, bins=bins, edgecolor="black")
    plt.xlabel("WMAPE")
    plt.ylabel("Prédictions")
    plt.title("Distribution des MSE sur le jeu de test - n=92")
    plt.show()
    print(np.max(mse), np.min(mse))
    # print("─────────────────────────────────────────────────────\n")

    # 8. Dénormalisation et enregistrement dans le PDF
    with PdfPages("lstm/best_worst/20Best_predi.pdf") as pdf:
        fig_title = plt.figure(figsize=(10, 4))
        plt.axis("off")
        plt.text(
            0.5,
            0.6,
            "Rapport de prédictions LSTM",
            ha="center",
            va="center",
            fontsize=20,
            fontweight="bold",
        )
        plt.text(
            0.5,
            0.4,
            f"Modèle : {MODEL_PTH}  |  20 meilleures simulations",
            ha="center",
            va="center",
            fontsize=12,
        )
        pdf.savefig(fig_title)
        plt.close()

        # 20 meilleures prédictions
        print("=== 20 MEILLEURES prédictions ===")
        for idx in best_idx:
            save_simulation_to_pdf(
                pdf,
                idx,
                sim_ids_test,
                x_test,
                y_test_norm,
                preds_norm,
                y_std,
                y_mean,
                wmape_val,
                category_label="MEILLEURE",
            )

        # 20 pires prédictions

    with PdfPages("lstm/best_worst/20_Worst_predi.pdf") as pdf:
        fig_title = plt.figure(figsize=(10, 4))
        plt.axis("off")
        plt.text(
            0.5,
            0.6,
            "Rapport de prédictions LSTM",
            ha="center",
            va="center",
            fontsize=20,
            fontweight="bold",
        )
        plt.text(
            0.5,
            0.4,
            f"Modèle : {MODEL_PTH}  |  20 meilleures simulations",
            ha="center",
            va="center",
            fontsize=12,
        )
        pdf.savefig(fig_title)
        plt.close()
        print("\n=== 20 PIRES prédictions ===")
        for idx in worst_idx:
            save_simulation_to_pdf(
                pdf,
                idx,
                sim_ids_test,
                x_test,
                y_test_norm,
                preds_norm,
                y_std,
                y_mean,
                wmape_val,
                category_label="PIRE",
            )


if __name__ == "__main__":
    main()
