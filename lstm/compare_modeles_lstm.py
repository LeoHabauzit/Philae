from pathlib import Path
from torch import randn
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import torch
from torchmetrics.regression import WeightedMeanAbsolutePercentageError
import numpy as np
import example_utils
from lstm_example import RNNModel
import seaborn as sns

sns.set_style("white")


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


TRAIN_CSV = Path("lstm") / "dataset" / "train_dataset.csv"
TEST_CSV = Path("lstm") / "dataset" / "test_fea_dataset_23_augmented.csv"

# Architecture : doit être IDENTIQUE à celle utilisée à l'entraînement
INPUT_SIZE = 6
# HIDDEN_SIZE = 64
OUTPUT_SIZE = 6
NUM_LAYERS = 2
# fig, ((ax0, ax1, ax2, ax3)) = plt.subplots(figsize=(15, 10), nrows=4, ncols=1)
fig, ax = plt.subplots(figsize=(8, 6))

# Boxplots


files = {
    # "model_5000.pth": {
    #     "hidden_size": 64,
    #     "color": "green",
    #     "linestyle": "--",
    #     "label": "Sans fine tuning",
    #     "face color": (1, 0, 1, 0.5),
    #     "edgecolor": "pink",
    #     # "ax": ax0,
    # },
    # "model_finetuned_5000.pth": {
    #     "hidden_size": 64,
    #     "color": "green",
    #     "linestyle": "--",
    #     "label": "Avec fine tuning",
    #     "face color": (0, 0, 1, 0.5),
    #     "edgecolor": "blue",
    #     # "ax": ax0,
    # },
    # "model_finetuned_augmented_HS16.pth": {
    #     "hidden_size": 16,
    #     "color": "black",
    #     "linestyle": "-",
    #     "label": "Data augmentation - Hidden Size 16",
    #     "face color": (1, 0, 0, 0.5),
    #     "edgecolor": "red",
    #     # "ax": ax1,
    # },
    "model_finetuned_augmented_HS32.pth": {
        "hidden_size": 32,
        "color": "blue",
        "linestyle": "-",
        "label": " Hidden Size 32",
        "face color": (0, 0, 1, 0.5),
        "edgecolor": "blue",
        # "ax": ax3,
    },
    "model_finetuned_augmented.pth": {
        "hidden_size": 64,
        "color": "red",
        "linestyle": "--",
        "label": " Hidden Size 64",
        # "face color": (0, 0, 0, 0.5),
        "face color": "none",
        "edgecolor": "red",
        # "ax": ax2,
    },
}

# fig = plt.figure(figsize=(8, 6))
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
data = []
for i, (model_name, params) in enumerate(files.items(), start=1):
    HIDDEN_SIZE = params["hidden_size"]
    # color = params["color"]
    model = RNNModel(
        input_features_size=INPUT_SIZE,
        hidden_state_size=HIDDEN_SIZE,
        output_size=OUTPUT_SIZE,
        num_layers=NUM_LAYERS,
    )
    checkpoint = torch.load(model_name, map_location=device, weights_only=True)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    print(f"Modèle chargé depuis '{model_name}'")

    # 6. Prédictions sur tout le jeu de test
    print(f"Prédiction sur {x_test.shape[0]} simulations...")
    preds_norm = predict_all(model, x_test_norm, device)  # (N, T, 6) normalisé
    true_stress = y_test_norm * y_std + y_mean  # dénormalisé
    pred_stress = preds_norm * y_std + y_mean  # dénormalisé

    wmape = WeightedMeanAbsolutePercentageError()
    values = []
    for k in range(pred_stress.shape[0]):
        values.append(wmape(pred_stress[k], true_stress[k]))

    moyenne = np.mean(values)
    ecart_type = np.std(values, ddof=1)
    values = np.array(values)
    data.append(values)
    bins = np.linspace(0, 0.1, 60)
    # ax.boxplot(
    #     values,
    #     positions=[i],
    #     patch_artist=True,
    #     boxprops=dict(facecolor=params["color_supp"]),
    #     medianprops=dict(color="black", linewidth=2),
    # )
    plt.hist(
        values,
        bins=bins,
        histtype="stepfilled",
        facecolor=params["face color"],
        edgecolor=params["edgecolor"],
        ls=params["linestyle"],
        # facecolor=params["face color"],
        # density=True,
        label=params["label"]
        + f" (accu={(1 - moyenne) * 100:.2f}%)"
        + f" (écart-type={(1 - ecart_type) * 100:.2f}%)",
        lw=1.5,
    )
    # params["ax"].hist(
    #     values,
    #     bins=bins,
    #     edgecolor="black",
    #     color=params["color_supp"],
    #     alpha=0.7,
    #     label=params["label"],
    # )
# plt.set_xlabel("Erreur WMAPE (%)")
# plt.set_ylabel("Prédictions")
# params["ax"].set_title("Distribution des WMAPE sur le jeu de test")


# params["ax"].grid()
# params["ax"].set_ylim(0, 41)
plt.title("Comparaison des modèles LSTM sur le jeu de test")
plt.legend()
# for model_name, params in files.items():
ax.set_xlabel("Erreur WMAPE")
ax.set_ylabel("Nombre de prédictions")
# # Boxplots
# ax.boxplot(
#     data,
#     patch_artist=True,
#     boxprops=dict(facecolor="lightblue"),
#     medianprops=dict(color="black", linewidth=2),
# )

# # Boucle pour ajouter moyenne et écart-type
# for i, valeurs in enumerate(data, start=1):
#     moyenne = np.mean(valeurs)
#     ecart_type = np.std(valeurs, ddof=1)

#     # ax.errorbar(i, moyenne, yerr=ecart_type, fmt="o", color="red", capsize=5)

# ax.set_ylim(90, 100)
# ax.set_ylabel("Valeur")
# ax.set_title("Boîtes à moustaches avec moyenne ± écart-type")
# ax.grid(axis="y", linestyle="--", alpha=0.5)

# plt.tight_layout()
# plt.tight_layout()
plt.show()
