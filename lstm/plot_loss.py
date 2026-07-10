"""
plot_loss.py
Recharge et affiche les courbes de loss depuis le fichier model.pth
(nécessite que model.pth ait été sauvegardé en format checkpoint complet)
"""

import torch
import example_utils

MODEL_PTH = "model_NO_SMAAC.pth"


def main() -> None:
    # Charger le checkpoint
    checkpoint = torch.load(MODEL_PTH, map_location="cpu", weights_only=True)

    # Vérifier que les courbes sont bien présentes
    if "train_losses" not in checkpoint:
        print("Erreur : ce fichier model ne contient pas les courbes de loss.")

        return

    train_losses = checkpoint["train_losses"]
    test_losses = checkpoint["test_losses"]

    print(f"Nombre d'époques : {len(train_losses)}")
    print(f"Loss finale train : {train_losses[-1]:.6f}")
    print(f"Loss finale test  : {test_losses[-1]:.6f}")

    # Afficher le graphique (utilise la fonction déjà existante dans example_utils)
    example_utils.plot_loss_curves(train_losses, test_losses)


if __name__ == "__main__":
    main()
