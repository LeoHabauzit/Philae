import torch
import pandas as pd
from lstm_example import RNNModel  # ou là où ta classe est définie
import example_utils
from pathlib import Path

# Choisir le bon device
if torch.backends.mps.is_available():
    device = "mps"
else:
    device = "cpu"

inputs_headers = [
    "total_strain_xx",
    "total_strain_yy",
    "total_strain_zz",
    "total_strain_xy",
    "total_strain_xz",
    "total_strain_yz",
]
TRAIN_CSV = Path("lstm") / "dataset" / "train_dataset.csv"
# TEST_CSV = Path("lstm") / "dataset" / "test_dataset.csv"
TEST_CSV = Path("lstm") / "dataset" / "test_fea_dataset_23.csv"
x_train, y_train, _ = example_utils.load_data(TRAIN_CSV.as_posix())
x_test, y_test, sim_ids_test = example_utils.load_data(TEST_CSV.as_posix())

x_mean, x_std = x_train.mean(), x_train.std()
y_mean, y_std = y_train.mean(), y_train.std()

# 4. Normalisation du jeu de test
x_test_norm = (x_test - x_mean) / x_std
y_test_norm = (y_test - y_mean) / y_std
SAMPLE_INDEX = 0  # index de la simulation à prédire
# Recréer le modèle avec la MÊME architecture qu'à l'entraînement
model = RNNModel(
    input_features_size=6,
    hidden_state_size=64,
    output_size=6,
    num_layers=2,
)

# Charger les poids sauvegardés
checkpoint = torch.load("model_finetuned.pth", map_location=device, weights_only=True)
model.load_state_dict(checkpoint["model_state_dict"])
model.to(device)


# Mettre en mode "prédiction" (désactive le dropout etc.)
model.eval()


# Récupérer les résultats en numpy
strain_norm = x_test_norm[SAMPLE_INDEX].unsqueeze(0).to(device)  # (1, T, 6)
true_stress_norm = y_test_norm[SAMPLE_INDEX]  # (T, 6)

with torch.no_grad():
    pred_stress_norm, _ = model(strain_norm)

pred_stress_norm = pred_stress_norm.squeeze(0).cpu()  # (T, 6)

# 7. Dénormalisation
strain_sample = x_test[SAMPLE_INDEX]  # déjà en valeurs réelles
true_stress = true_stress_norm * y_std + y_mean
pred_stress = pred_stress_norm * y_std + y_mean

print(strain_sample.shape, pred_stress.shape)

# stress_array = pred_stress[:, i].detach().numpy()
targets_headers = [
    "total_strain_xx",
    "total_strain_yy",
    "total_strain_zz",
    "total_strain_xy",
    "total_strain_xz",
    "total_strain_yz",
    "stress_xx",
    "stress_yy",
    "stress_zz",
    "stress_xy",
    "stress_xz",
    "stress_yz",
]

data = torch.cat([strain_sample, pred_stress], dim=1)  # [100, 12]

# Conversion vers pandas
df = pd.DataFrame(data.numpy(), columns=targets_headers)

# Sauvegarde CSV
df.to_csv("simuEF/csv_files/compare_lstm_tuned.csv", index=False)

# df_predictions = pd.DataFrame(predictions, columns=targets_headers)
# df_predictions.to_csv("simuEF/csv_files/compare_lstm_tuned.csv", index=False)
# print("Sauvegardé !")
