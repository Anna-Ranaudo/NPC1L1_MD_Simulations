#!/usr/bin/env python3
"""Create a readable five-replica summary from protein-cholesterol H-bond CSVs."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


REPLICA_NAMES = ["md1", "md2", "md3", "md4", "md5"]


def parse_args():
    script_dir = Path(__file__).resolve().parent
    default_input = script_dir.parent / "final_data/open/bound/prot_chol_hbond"
    default_output = script_dir.parent / "final_plots/open/bound"

    parser = argparse.ArgumentParser(
        description="Plot protein-cholesterol hydrogen-bond presence for five replicas."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input,
        help="Directory containing md1..md5_hbonds_presence_final.csv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output,
        help="Directory in which the combined PNG is written.",
    )
    parser.add_argument(
        "--output-name",
        default="protein-cholesterol-hbonds_5replicas.png",
        help="Name of the combined PNG file.",
    )
    parser.add_argument(
        "--title",
        default="Protein-cholesterol hydrogen bonds | open bound",
        help="Figure title.",
    )
    return parser.parse_args()


def load_replica_csvs(input_dir):
    data = {}
    missing = []

    for replica in REPLICA_NAMES:
        csv_path = input_dir / f"{replica}_hbonds_presence_final.csv"
        if not csv_path.is_file():
            missing.append(str(csv_path))
            continue

        frame_data = pd.read_csv(csv_path)
        if "Time (ns)" not in frame_data.columns:
            raise ValueError(f"Missing 'Time (ns)' column in {csv_path}")
        if len(frame_data.columns) < 2:
            raise ValueError(f"No hydrogen-bond columns found in {csv_path}")
        data[replica] = frame_data

    if missing:
        missing_paths = "\n  ".join(missing)
        raise FileNotFoundError(f"Missing replica CSV file(s):\n  {missing_paths}")

    return data


def plot_replicas(replica_data, output_path, title):
    labels = list(replica_data[REPLICA_NAMES[0]].columns[1:])
    for replica, frame_data in replica_data.items():
        replica_labels = list(frame_data.columns[1:])
        if replica_labels != labels:
            raise ValueError(
                f"Hydrogen-bond columns differ between replicas: {replica}"
            )

    n_bonds = len(labels)
    figure_height = max(9.5, 2.2 * len(REPLICA_NAMES))
    figure, axes = plt.subplots(
        len(REPLICA_NAMES),
        1,
        figsize=(15, figure_height),
        sharex=True,
        sharey=True,
        gridspec_kw={"hspace": 0.28},
    )
    axes = np.atleast_1d(axes)

    cmap = ListedColormap(["#f1f3f2", "#147d86"])
    max_time = max(frame_data["Time (ns)"].max() for frame_data in replica_data.values())

    for axis, replica in zip(axes, REPLICA_NAMES):
        frame_data = replica_data[replica]
        time = frame_data["Time (ns)"].to_numpy(dtype=float)
        matrix = frame_data.iloc[:, 1:].to_numpy(dtype=float).T
        matrix = np.where(matrix > 0, 1, 0)

        # Cell edges preserve the exact binary presence information without
        # producing the dense black vertical strokes of the old line plot.
        if len(time) > 1:
            step = np.median(np.diff(time))
        else:
            step = 1.0
        time_edges = np.r_[time - step / 2, time[-1] + step / 2]
        bond_edges = np.arange(n_bonds + 1)
        axis.pcolormesh(
            time_edges,
            bond_edges,
            matrix,
            cmap=cmap,
            vmin=0,
            vmax=1,
            shading="flat",
            rasterized=True,
        )

        axis.set_title(replica.upper(), loc="left", fontsize=12, weight="bold", pad=7)
        axis.set_yticks(np.arange(n_bonds) + 0.5)
        axis.set_yticklabels(labels, fontsize=9)
        axis.grid(axis="x", color="white", linewidth=0.7, alpha=0.8)
        axis.set_axisbelow(False)
        axis.spines[["top", "right"]].set_visible(False)

    axes[-1].set_xlabel("Time (ns)", fontsize=11, labelpad=8)
    for axis in axes:
        axis.set_xlim(0, max_time + 5)
        axis.set_xticks(np.arange(0, max_time + 1, 50))
        axis.tick_params(axis="x", labelsize=9)
        axis.tick_params(axis="y", length=0)

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#f1f3f2", edgecolor="#b7c2c0", label="Absent"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#147d86", label="Present"),
    ]
    figure.suptitle(
        title,
        x=0.28,
        y=0.995,
        ha="left",
        fontsize=17,
        weight="bold",
    )
    figure.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.995),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    figure.subplots_adjust(left=0.28, right=0.98, top=0.94, bottom=0.08)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, facecolor="white")
    plt.close(figure)


def main():
    args = parse_args()
    replica_data = load_replica_csvs(args.input_dir)
    output_path = args.output_dir / args.output_name
    plot_replicas(replica_data, output_path, args.title)
    print(f"Saved combined plot: {output_path}")


if __name__ == "__main__":
    main()