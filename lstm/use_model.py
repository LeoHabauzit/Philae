import torch
import pandas as pd
from lstm_example import RNNModel  # ou là où ta classe est définie

# Choisir le bon device
if torch.backends.mps.is_available():
    device = "mps"
else:
    device = "cpu"

# Recréer le modèle avec la MÊME architecture qu'à l'entraînement
model = RNNModel(
    input_features_size=6,
    hidden_state_size=64,
    output_size=6,
    num_layers=2,
)

# Charger les poids sauvegardés
model.load_state_dict(torch.load("model_save.pth", map_location=device))
model.to(device)

# Mettre en mode "prédiction" (désactive le dropout etc.)
model.eval()


df = pd.read_csv("lstm/dataset/test_dataset.csv")

inputs_headers = [
    "total_strain_xx",
    "total_strain_yy",
    "total_strain_zz",
    "total_strain_xy",
    "total_strain_xz",
    "total_strain_yz",
]

# Préparer les données : même format qu'à l'entraînement
data = torch.tensor(df[inputs_headers].values, dtype=torch.float32)

# Le LSTM attend la forme : (séquence, batch, features)
# Ici on met batch=1 car on prédit une seule séquence
data = data.unsqueeze(1).to(device)  # shape : (T, 1, 6)

# Prédire !
with torch.no_grad():  # pas besoin de calculer les gradients
    predictions = model(data)

# Récupérer les résultats en numpy
predictions, _ = predictions
predictions = predictions.squeeze(1).cpu().numpy()
print(predictions.shape)
targets_headers = [
    "stress_xx",
    "stress_yy",
    "stress_zz",
    "stress_xy",
    "stress_xz",
    "stress_zy",
]

df_predictions = pd.DataFrame(predictions, columns=targets_headers)
df_predictions.to_csv("predictions.csv", index=False)
print("Sauvegardé !")
