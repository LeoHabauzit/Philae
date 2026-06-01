"""
StressLSTM
Copyright (C) 2025 Manuel Ricardo GUEVARA GARBAN

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

from pathlib import Path
from typing import Callable, Optional, Tuple

import torch
from torch import Tensor
from torch.utils.data import DataLoader, TensorDataset

import example_utils


class RNNModel(torch.nn.Module):
    def __init__(
        self,
        input_features_size: int,
        hidden_state_size: int,
        output_size: int,
        num_layers: int,
    ):
        super().__init__()
        self.input_features_size = input_features_size
        self.hidden_state_size = hidden_state_size
        self.output_size = output_size
        self.num_layers = num_layers
        self.rnn = torch.nn.LSTM(
            input_features_size,
            hidden_state_size,
            num_layers,
            batch_first=True,
        )
        self.mlp_output_layer = torch.nn.Linear(hidden_state_size, output_size)

    def forward(
        self,
        input_sequence: Tensor,
        hidden_states: Optional[Tensor] = None,
        cell_state: Optional[Tensor] = None,
    ) -> Tuple[Tensor, Tuple[Tensor, Tensor]]:
        if hidden_states is None or cell_state is None:
            output, (hidden_states, cell_state) = self.rnn(input_sequence)
        else:
            output, (hidden_states, cell_state) = self.rnn(
                input_sequence, (hidden_states, cell_state)
            )
        preds = self.mlp_output_layer(output)
        return preds, (hidden_states, cell_state)


def train(
    model: torch.nn.Module,
    criterion: Callable[[Tensor, Tensor], Tensor],
    optimizer: torch.optim.Optimizer,
    train_dataloader: DataLoader,
    test_dataloader: DataLoader,
    num_epochs: int = 100,
    device: str = "cuda",
) -> tuple[list[float], list[float]]:
    model.to(device)
    train_losses: list[float] = []
    test_losses: list[float] = []
    for epoch in range(num_epochs):
        model.train()
        train_loss_total = 0.0

        for input_sequence, target_sequence in train_dataloader:
            input_sequence, target_sequence = (
                input_sequence.to(device),
                target_sequence.to(device),
            )
            optimizer.zero_grad()
            predicted_sequence, _ = model.forward(input_sequence)

            loss = criterion(predicted_sequence, target_sequence)
            loss.backward()
            optimizer.step()
            train_loss_total += loss.item()

        # Evaluation
        model.eval()
        test_loss_total = 0.0
        with torch.no_grad():
            for input_sequence, target_sequence in test_dataloader:
                input_sequence, target_sequence = (
                    input_sequence.to(device),
                    target_sequence.to(device),
                )
                predicted_sequence, _ = model(input_sequence)
                loss = criterion(predicted_sequence, target_sequence)
                test_loss_total += loss.item()

        avg_train_loss = train_loss_total / len(train_dataloader)
        avg_test_loss = test_loss_total / len(test_dataloader)

        print(
            f"Epoch [{epoch + 1}/{num_epochs}] "
            f"Train Loss: {avg_train_loss:.6f} "
            f"Test Loss: {avg_test_loss:.6f}"
        )
        train_losses.append(avg_train_loss)
        test_losses.append(avg_test_loss)
    return train_losses, test_losses


def main() -> None:
    if torch.backends.mps.is_available():
        device = "mps"
    elif torch.cuda.is_available():
        device = "cuda"
    else:
        device = "cpu"

    print(f"Using device = {device}")
    train_dataset_path = Path("lstm") / Path("dataset") / Path("train_dataset.csv")
    test_dataset_path = Path("lstm") / Path("dataset") / Path("test_dataset.csv")
    # Load data
    x_train, y_train, sim_ids_train = example_utils.load_data(
        train_dataset_path.as_posix()
    )
    x_test, y_test, sim_ids_test = example_utils.load_data(test_dataset_path.as_posix())
    for i in range(3):
        example_utils.plot_stress_strain_sample(
            x_train, y_train, sim_ids_train, sim_index=i
        )
    # Normalize data using standardization
    x_mean, x_std = x_train.mean(), x_train.std()
    y_mean, y_std = y_train.mean(), y_train.std()
    x_train = (x_train - x_mean) / x_std
    y_train = (y_train - y_mean) / y_std
    x_test = (x_test - x_mean) / x_std
    y_test = (y_test - y_mean) / y_std
    # Datasets and loaders
    train_dataset = TensorDataset(x_train, y_train)
    test_dataset = TensorDataset(x_test, y_test)
    batch_size = 512
    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size)

    print(f"{x_train.shape=}")
    print(f"{y_train.shape=}")
    print(f"{sim_ids_train.shape=}")
    print(f"{x_test.shape=}")
    print(f"{y_test.shape=}")
    print(f"{sim_ids_test.shape=}")

    # Model, loss, optimizer
    model = RNNModel(
        input_features_size=6,
        hidden_state_size=64,
        output_size=6,
        num_layers=2,
    )
    criterion = torch.nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

    # Train model
    num_epochs = 500
    train_loss_curve, test_loss_curve = train(
        model,
        criterion,
        optimizer,
        train_dataloader,
        test_dataloader,
        num_epochs=num_epochs,
        device=device,
    )
    example_utils.plot_loss_curves(train_loss_curve, test_loss_curve)
    weights_path = "model.pth"
    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "train_losses": train_loss_curve,
            "test_losses": test_loss_curve,
        },
        weights_path,
    )
    # Test the model with a sample from test dataset
    with torch.no_grad():
        model.load_state_dict(torch.load(weights_path, weights_only=True))
        model.eval()

        # Choose a sample from test set
        sample_index = 0
        strain_sample = x_test[sample_index].unsqueeze(0)  # Add batch dimension
        true_stress_sample = y_test[sample_index]

        # Predict
        model.eval()
        with torch.no_grad():
            # No need to normalize the input data because it is already normalized
            predicted_stress_sample, _ = model(strain_sample.to(device))

        # Remove batch dimension for plotting and send to cpu
        predicted_stress_sample = predicted_stress_sample.squeeze(0).cpu()
        strain_sample = strain_sample.squeeze(0)
        # Denormalize data
        strain_sample = (strain_sample * x_std) + x_mean
        true_stress_sample = (true_stress_sample * y_std) + y_mean
        predicted_stress_sample = (predicted_stress_sample * y_std) + y_mean

        # Plot comparison
        example_utils.plot_stress_strain_sample_with_prediction(
            strain_sequence=strain_sample,
            true_stress_sequence=true_stress_sample,
            predicted_stress_sequence=predicted_stress_sample,
            component_labels=["xx", "yy", "zz", "xy", "xz", "yz"],
        )


if __name__ == "__main__":
    main()
