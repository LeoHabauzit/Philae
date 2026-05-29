import torch
from example_utils import RNNModel

device = "mps" if torch.backends.mps.is_available() else "cpu"

# 1. Charger le modèle pré-entraîné (pareil qu'avant)
model = RNNModel(
    input_features_size=6,
    hidden_state_size=64,
    output_size=6,
    num_layers=2,
)
model.load_state_dict(torch.load("model.pth", map_location=device))
model.to(device)

# 2. PAS de model.eval() cette fois ! On reste en mode entraînement
model.train()

# 3. Utiliser un learning rate PLUS PETIT qu'à l'entraînement initial
#    (pour ne pas "effacer" ce que le modèle a déjà appris)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)  # au lieu de 1e-3 par exemple

# 4. Entraîner normalement avec ton nouveau dataset
# ... (même boucle d'entraînement qu'avant)

# 5. Sauvegarder le modèle affiné
torch.save(model.state_dict(), "model_finetuned.pth")
