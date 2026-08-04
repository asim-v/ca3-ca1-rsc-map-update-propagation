#!/usr/bin/env python3
"""Plot the audited unit and traversal support without neural effect estimates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", type=Path, default=Path("data/metadata/navigation_session_ledger.csv"))
    parser.add_argument("--output", type=Path, default=Path("outputs/figures/dataset_audit.png"))
    args = parser.parse_args()
    table = pd.read_csv(args.ledger)
    table["label"] = table.subject + "\n" + table.date.astype(str).str[4:]
    x = np.arange(len(table))
    colors = {"CA1": "#355C7D", "CA3": "#C06C84", "RSC": "#6C5B7B"}

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    axes[0].plot(x, table.n_ca1_pyramidal, "o-", color=colors["CA1"], label="CA1 pyramidal")
    axes[0].plot(x, table.n_ca3_pyramidal, "o-", color=colors["CA3"], label="CA3 pyramidal")
    axes[0].plot(x, table.n_rsc_non_narrow, "o-", color=colors["RSC"], label="RSC non-narrow")
    axes[0].set_ylabel("Audited units")
    axes[0].legend(frameon=False, ncol=3)

    bar_colors = np.where(table.is_two_maze, "#F67280", "#99B898")
    axes[1].bar(x, table.n_traversals, color=bar_colors)
    axes[1].set_ylabel("Detected traversals")
    axes[1].set_xticks(x, table.label, rotation=45, ha="right")
    axes[1].set_xlabel("Subject and session date (MMDD)")
    axes[1].text(0.99, 0.95, "coral = familiar two-maze control", transform=axes[1].transAxes, ha="right", va="top")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#DDDDDD", linewidth=0.7, zorder=0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=200)
    plt.close(fig)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
