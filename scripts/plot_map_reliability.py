#!/usr/bin/env python3
"""Plot the region-blind late-map reliability audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, nargs="?", default=Path("outputs/quality/map_reliability.csv"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/figures/map_reliability_audit.png"),
    )
    args = parser.parse_args()

    table = pd.read_csv(args.input)
    regions = ["CA3", "CA1", "RSC"]
    colors = {"CA3": "#577590", "CA1": "#43AA8B", "RSC": "#F3722C"}
    rng = np.random.default_rng(20260803)

    figure, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    values = [
        table.loc[table["region"] == region, "spearman_brown_reliability"].dropna().to_numpy()
        for region in regions
    ]
    boxes = ax.boxplot(values, positions=np.arange(1, 4), widths=0.48, patch_artist=True, showfliers=False)
    for patch, region in zip(boxes["boxes"], regions, strict=True):
        patch.set(facecolor=colors[region], alpha=0.20, edgecolor=colors[region], linewidth=1.5)
    for element in ("whiskers", "caps", "medians"):
        for artist in boxes[element]:
            artist.set(color="#333333", linewidth=1.25)

    for position, (region, region_values) in enumerate(zip(regions, values, strict=True), start=1):
        jitter = rng.uniform(-0.13, 0.13, len(region_values))
        ax.scatter(
            position + jitter,
            region_values,
            s=23,
            color=colors[region],
            alpha=0.72,
            edgecolor="white",
            linewidth=0.35,
        )
    ax.axhline(0.15, color="#6D597A", linestyle="--", linewidth=1.2, label="0.15 sensitivity cutoff")
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set(
        title="Late-map reliability before onset analysis",
        ylabel="Spearman-Brown split-half reliability",
        xticks=[1, 2, 3],
        xticklabels=regions,
        ylim=(-0.03, 0.92),
    )
    ax.legend(frameon=False, loc="lower center")
    ax.spines[["top", "right"]].set_visible(False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=180)
    plt.close(figure)


if __name__ == "__main__":
    main()
