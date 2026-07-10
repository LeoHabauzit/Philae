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

from typing import Callable, List, Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import torch
from torch import Tensor
import numpy as np
import os


def plot_stress_strain_curve(
    strain_sequence: Tensor,
    stress_sequence: Tensor,
    component_labels: List[str],
    title: str = "Stress–Strain Curves",
    sim_id: Optional[int] = None,
) -> None:
    """
    Plot each stress–strain component in a 2x3 grid of subplots.
    """
    n_components = strain_sequence.shape[1]

    n_rows = 2
    n_cols = 3

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(5 * n_cols, 4 * n_rows),
        sharex=True,
    )

    axes = axes.flatten()

    for i in range(n_components):
        ax = axes[i]
        ax.plot(
            strain_sequence[:, i].numpy(),
            stress_sequence[:, i].numpy(),
            label=component_labels[i],
            linewidth=2,
            color="blue",
        )
        ax.set_xlabel("Strain", fontsize=12)
        ax.set_ylabel("Stress", fontsize=12)
        ax.grid(True)
        ax.legend(fontsize=10)
        ax.set_title(f"Component {component_labels[i]}", fontsize=12)

    # Cache les axes inutilisés si n_components < 6
    for j in range(n_components, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        title if sim_id is None else f"{title} (Simulation {sim_id})",
        fontsize=14,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()


def plot_stress_strain_sample(
    X: Tensor,
    Y: Tensor,
    sim_ids: Tensor,
    sim_index: int = 0,
) -> None:
    """
    Plot true stress–strain curves in separate subplots.
    """
    strain = X[sim_index]
    stress = Y[sim_index]
    sim_id = sim_ids[sim_index].item()
    plot_stress_strain_curve(
        strain,
        stress,
        component_labels=["xx", "yy", "zz", "xy", "xz", "yz"],
        sim_id=sim_id,
    )


def plot_stress_strain_sample_with_prediction(
    strain_sequence: Tensor,
    true_stress_sequence: Tensor,
    predicted_stress_sequence: Tensor,
    sim_id: Optional[int] = None,
    component_labels: List[str] = ["xx", "yy", "zz", "xy", "xz", "yz"],
    title: str = "Stress–Strain Prediction vs Ground Truth",
    pdf_file=None,
) -> None:
    """
    Plot predicted vs true stress–strain curves in separate subplots.
    """
    n_components = strain_sequence.shape[1]
    n_rows = 3
    n_cols = 2
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows), sharex=True
    )
    axes = axes.flatten()
    for i in range(n_components):
        ax = axes[i]
        ax.plot(
            strain_sequence[:, i].numpy(),
            true_stress_sequence[:, i].numpy(),
            label="True",
            linewidth=2,
            color="blue",
        )
        ax.plot(
            strain_sequence[:, i].numpy(),
            predicted_stress_sequence[:, i].detach().numpy(),
            linestyle="--",
            linewidth=2,
            color="red",
            label="Predicted",
        )
        ax.set_ylabel("Stress [MPa]", fontsize=12)
        ax.set_xlabel("Strain [-]")
        ax.grid(True)
        ax.legend(fontsize=10)
        ax.set_title(f"Component {component_labels[i]}", fontsize=12)

    axes[-1].set_xlabel("Strain", fontsize=12)
    fig.suptitle(
        f"{title}" if sim_id is None else f"{title} (Simulation {sim_id})",
        fontsize=14,
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if pdf_file is None:
        plt.show()
    else:
        pdf_file.savefig(plt.gcf())


def plot_loss_curves(
    train_losses: List[float],
    val_losses: Optional[List[float]] = None,
    title: str = "Training and Validation Loss",
    xlabel: str = "Epoch",
    ylabel: str = "Loss",
    legend_labels: Optional[List[str]] = None,
    save_path: Optional[str] = None,
) -> None:
    """
    Plot training and (optional) validation loss curves.

    Args:
        train_losses (List[float]): Training loss per epoch.
        val_losses (Optional[List[float]]): Validation loss per epoch (same length as train_losses or shorter).
        title (str): Plot title.
        xlabel (str): X-axis label.
        ylabel (str): Y-axis label.
        legend_labels (Optional[List[str]]): Custom labels for the legend.
        save_path (Optional[str]): If provided, saves the figure to this path.
    """
    plt.figure(figsize=(8, 5))
    epochs = range(1, len(train_losses) + 1)

    plt.plot(
        epochs,
        train_losses,
        label=legend_labels[0] if legend_labels else "Train",
        color="black",
        linewidth=2,
    )

    plt.yscale("log")

    if val_losses is not None:
        val_epochs = range(1, len(val_losses) + 1)
        plt.plot(
            val_epochs,
            val_losses,
            label=legend_labels[1] if legend_labels else "Validation",
            linestyle="--",
            color="red",
            linewidth=2,
        )

    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(fontsize=10)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()


def _get_input_and_target_seq_from_df_sample(
    df_sample: pd.DataFrame,
    inputs_headers: list[str],
    targets_headers: list[str],
) -> Tuple[Tensor, Tensor]:
    total_strain = torch.tensor(df_sample[inputs_headers].values, dtype=torch.float32)
    stress = torch.tensor(df_sample[targets_headers].values, dtype=torch.float32)
    return total_strain, stress


def _get_inputs_targets_and_sims_ids_data_from_dataframe(
    df: pd.DataFrame,
    inputs_headers: list[str],
    targets_headers: list[str],
) -> Tuple[Tensor, Tensor, Tensor]:
    simulation_ids = list(set(df["simulation_load_id"]))
    max_seq_length: int = df["simulation_load_id"].value_counts().max()
    n_sequences = len(simulation_ids)

    # Determine feature dimensions
    tmp_seq_input, tmp_seq_target = _get_input_and_target_seq_from_df_sample(
        df[df["simulation_load_id"] == simulation_ids[0]],
        inputs_headers,
        targets_headers,
    )
    num_input_features = tmp_seq_input.shape[1]
    num_target_features = tmp_seq_target.shape[1]

    padded_input_sequences = torch.zeros(
        n_sequences, max_seq_length, num_input_features
    )
    padded_target_sequences = torch.zeros(
        n_sequences, max_seq_length, num_target_features
    )

    for i, sim_id in enumerate(simulation_ids):
        sample_df = df[df["simulation_load_id"] == sim_id]
        seq_input, seq_target = _get_input_and_target_seq_from_df_sample(
            sample_df, inputs_headers, targets_headers
        )
        seq_len = seq_input.shape[0]
        padded_input_sequences[i, :seq_len, :] = seq_input
        padded_target_sequences[i, :seq_len, :] = seq_target

    return (
        padded_input_sequences,
        padded_target_sequences,
        torch.tensor(simulation_ids),
    )


def load_data(csv_path: str) -> Tuple[Tensor, Tensor, Tensor]:
    dataframe = pd.read_csv(csv_path)

    # Filter and select headers
    inputs_headers = [
        "total_strain_xx",
        "total_strain_yy",
        "total_strain_zz",
        "total_strain_xy",
        "total_strain_xz",
        "total_strain_yz",
    ]
    targets_headers = [
        "stress_xx",
        "stress_yy",
        "stress_zz",
        "stress_xy",
        "stress_xz",
        "stress_yz",
    ]
    dataframe = dataframe[dataframe["timestep"] != 0]

    return _get_inputs_targets_and_sims_ids_data_from_dataframe(
        dataframe, inputs_headers, targets_headers
    )
