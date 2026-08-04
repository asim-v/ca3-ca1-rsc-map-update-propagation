#!/usr/bin/env python3
"""Plot mouse-level traversal-scale cross-area prediction effects."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ORDER = [("CA3", "CA1"), ("CA1", "CA3"), ("CA1", "RSC"), ("RSC", "CA1")]
COLORS = {"forward": "#277DA1", "reverse": "#F3722C"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mouse", type=Path, default=Path("outputs/cross_area/mouse_effects.csv")
    )
    parser.add_argument(
        "--summary", type=Path, default=Path("outputs/cross_area/summary.csv")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/figures/cross_area_prediction.png"),
    )
    args = parser.parse_args()
    mouse = pd.read_csv(args.mouse)
    summary = pd.read_csv(args.summary)

    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.6), sharey=True, constrained_layout=True)
    for lag, ax in enumerate(axes):
        for position, (source, target) in enumerate(ORDER):
            subset = mouse[
                (mouse["source"] == source)
                & (mouse["target"] == target)
                & (mouse["lag_same_direction_traversals"] == lag)
            ]
            direction = "forward" if (source, target) in (("CA3", "CA1"), ("CA1", "RSC")) else "reverse"
            values = subset["mouse_median_delta_mse"].to_numpy(float)
            jitter = np.linspace(-0.08, 0.08, len(values))
            ax.scatter(
                position + jitter,
                values,
                s=46,
                color=COLORS[direction],
                alpha=0.82,
                edgecolor="white",
                linewidth=0.5,
                zorder=3,
            )
            ax.plot(
                [position - 0.16, position + 0.16],
                [np.mean(values), np.mean(values)],
                color="#222222",
                linewidth=2.2,
            )
            result = summary[
                (summary["source"] == source)
                & (summary["target"] == target)
                & (summary["lag_same_direction_traversals"] == lag)
            ].iloc[0]
            if bool(result["supported_positive_gain"]):
                ax.text(position, 0.225, "q<0.05", ha="center", va="top", fontsize=8)
        ax.axhline(0, color="#333333", linewidth=1)
        ax.set(
            title="Same traversal" if lag == 0 else "Previous same-direction traversal",
            xticks=np.arange(4),
            xticklabels=[f"{source}->{target}" for source, target in ORDER],
            xlabel="Source added to target-history model",
            ylim=(-0.09, 0.235),
        )
        ax.tick_params(axis="x", rotation=25)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("LOOCV MSE reduction (positive = predictive gain)")
    figure.suptitle(
        "Reference-independent cross-area prediction is heterogeneous and not lagged",
        fontsize=13,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
