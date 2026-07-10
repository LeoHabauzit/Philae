import torch
from pathlib import Path
import example_utils
from lstm_example import *
from tools_database import write_shuffled_repartition_dataset

# ─────────────────────────────────────────────
# PARAMÈTRES À MODIFIER SI BESOIN
# ─────────────────────────────────────────────
TRAIN_1_CSV = Path("lstm") / "dataset" / "train_dataset.csv"

write_shuffled_repartition_dataset("lstm/dataset/all_fea_cuboctahedron_augmented.csv")
TRAIN_FT_CSV = Path("lstm") / "dataset" / "train_shuffled.csv"
# VALIDATION_CSV = Path("lstm") / "dataset" / "test_dataset.csv"
VALIDATION_CSV = Path("lstm") / "dataset" / "validation_shuffled.csv"


MODEL_PTH = "model_2000_HS16.pth"

# Architecture : doit être IDENTIQUE à celle utilisée à l'entraînement
# INPUT_SIZE = 6
# HIDDEN_SIZE = 32
# OUTPUT_SIZE = 6
# NUM_LAYERS = 2

device = "mps" if torch.backends.mps.is_available() else "cpu"
print(f"Using device = {device}")
train_dataset_path = TRAIN_FT_CSV
test_dataset_path = VALIDATION_CSV
# Load data
x_train1, y_train1, sim_ids_train1 = example_utils.load_data(TRAIN_1_CSV.as_posix())
x_train, y_train, sim_ids_train = example_utils.load_data(train_dataset_path.as_posix())
x_test, y_test, sim_ids_test = example_utils.load_data(test_dataset_path.as_posix())
# for i in range(3):
#     example_utils.plot_stress_strain_sample(
#         x_train, y_train, sim_ids_train, sim_index=i
#     )
# Normalize data using standardization
x_mean, x_std = x_train1.mean(), x_train1.std()
y_mean, y_std = y_train1.mean(), y_train1.std()
x_train = (x_train - x_mean) / x_std
y_train = (y_train - y_mean) / y_std
x_test = (x_test - x_mean) / x_std
y_test = (y_test - y_mean) / y_std
# Datasets and loaders
train_dataset = TensorDataset(x_train, y_train)
test_dataset = TensorDataset(x_test, y_test)
batch_size = 32

train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
test_dataloader = DataLoader(test_dataset, batch_size=batch_size)


# # 1. Charger le modèle pré-entraîné (pareil qu'avant)
model = RNNModel_droppout(
    input_features_size=6,
    hidden_state_size=16,
    output_size=6,
    num_layers=2,
    dropout_p=0.0,
)
# freeze les params de du LSTM
# for param in model.rnn.parameters():
#     param.requires_grad = False
criterion = torch.nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001, weight_decay=1e-5)
checkpoint = torch.load(MODEL_PTH, map_location=device, weights_only=True)
model.load_state_dict(checkpoint["model_state_dict"])
model.to(device)

num_epochs = 5000
train_loss_curve, test_loss_curve = train(
    model,
    criterion,
    optimizer,
    train_dataloader,
    test_dataloader,
    num_epochs=num_epochs,
    device=device,
)

weights_path = "model_finetuned_augmented_HS16.pth"
torch.save(
    {
        "model_state_dict": model.state_dict(),
        "train_losses": train_loss_curve,
        "test_losses": test_loss_curve,
    },
    weights_path,
)
