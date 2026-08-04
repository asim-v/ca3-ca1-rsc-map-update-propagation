#!/usr/bin/env python3
"""Summarize all regional map-maturity curves without cell/trial pseudoreplication."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


KEYS = ["subject", "file", "block", "direction"]
REGIONS = ["CA3", "CA1", "RSC"]
STATUS_ORDER = {
    "unresolved_nonmonotonic_or_weak": 0,
    "present_from_first_observed_traversals": 1,
    "detectable_transition": 2,
}


def date_from_filename(filename: str) -> str:
    match = re.search(r"ses-(\d{8})", filename)
    return match.group(1)[4:] if match else filename


def paired_transition_summary(table: pd.DataFrame, status_column: str, midpoint_column: str) -> dict:
    status = table.pivot(index=KEYS, columns="region", values=status_column)
    midpoint = table.pivot(index=KEYS, columns="region", values=midpoint_column)
    all_three = status.eq("detectable_transition").all(axis=1)
    ca3_ca1 = status["CA3"].eq("detectable_transition") & status["CA1"].eq("detectable_transition")
    ca1_rsc = status["CA1"].eq("detectable_transition") & status["RSC"].eq("detectable_transition")
    return {
        "all_three_transition_curves": int(all_three.sum()),
        "ca3_ca1_comparable_curves": int(ca3_ca1.sum()),
        "ca3_minus_ca1_midpoints": [
            float(value) for value in (midpoint.loc[ca3_ca1, "CA3"] - midpoint.loc[ca3_ca1, "CA1"])
        ],
        "ca1_rsc_comparable_curves": int(ca1_rsc.sum()),
        "ca1_minus_rsc_midpoints": [
            float(value) for value in (midpoint.loc[ca1_rsc, "CA1"] - midpoint.loc[ca1_rsc, "RSC"])
        ],
    }


def make_figure(fits: pd.DataFrame, controlled: pd.DataFrame, output: Path) -> None:
    fits = fits.copy()
    fits["early_registered"] = (
        fits["early_identity_null_p"].le(0.05) & fits["early_spatial_null_p"].le(0.05)
    )
    row_table = fits[KEYS].drop_duplicates().sort_values(KEYS).reset_index(drop=True)
    matrix = np.zeros((len(row_table), len(REGIONS)), dtype=int)
    early = np.zeros_like(matrix, dtype=bool)
    row_labels = []
    for row_index, key_row in row_table.iterrows():
        key_mask = np.logical_and.reduce([fits[key] == key_row[key] for key in KEYS])
        subset = fits.loc[key_mask].set_index("region")
        for region_index, region in enumerate(REGIONS):
            matrix[row_index, region_index] = STATUS_ORDER[subset.loc[region, "fit_status"]]
            early[row_index, region_index] = bool(subset.loc[region, "early_registered"])
        arrow = "L->R" if key_row["direction"] == "LtoR" else "R->L"
        row_labels.append(f"{key_row['subject']} {date_from_filename(key_row['file'])} {arrow}")

    figure = plt.figure(figsize=(11.8, 8.3), constrained_layout=True)
    grid = figure.add_gridspec(1, 2, width_ratios=[1.05, 1.25])
    ax_heat = figure.add_subplot(grid[0, 0])
    cmap = ListedColormap(["#D9D9D9", "#4C956C", "#F4A261"])
    ax_heat.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=2.5)
    for row, column in np.argwhere(early):
        ax_heat.scatter(column, row, marker="o", s=26, facecolor="white", edgecolor="#222222", linewidth=0.7)
    ax_heat.set(
        title="Every session and direction",
        xticks=np.arange(3),
        xticklabels=REGIONS,
        yticks=np.arange(len(row_labels)),
        yticklabels=row_labels,
    )
    ax_heat.tick_params(axis="y", labelsize=7.6)
    ax_heat.set_xticks(np.arange(-0.5, 3, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=1.0)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    legend = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#D9D9D9", markersize=10, label="Unresolved/irregular"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#4C956C", markersize=10, label="Present from first traversals"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor="#F4A261", markersize=10, label="Detectable maturation"),
        Line2D([0], [0], marker="o", color="#222222", markerfacecolor="white", markersize=6, linewidth=0, label="Early registered identity > both nulls"),
    ]
    ax_heat.legend(handles=legend, frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.045))

    right = grid[0, 1].subgridspec(2, 1, height_ratios=[1.1, 1.0])
    ax_counts = figure.add_subplot(right[0, 0])
    early_counts = fits.groupby("region")["early_registered"].sum().reindex(REGIONS)
    raw_counts = (
        fits["fit_status"].eq("detectable_transition").groupby(fits["region"]).sum().reindex(REGIONS)
    )
    adjusted_counts = (
        controlled["adjusted_fit_status"]
        .eq("detectable_transition")
        .groupby(controlled["region"])
        .sum()
        .reindex(REGIONS)
    )
    x = np.arange(3)
    width = 0.23
    ax_counts.bar(x - width, early_counts, width, color="#277DA1", label="Early registered identity")
    ax_counts.bar(x, raw_counts, width, color="#F4A261", label="Raw maturation")
    ax_counts.bar(x + width, adjusted_counts, width, color="#9C6644", label="Behavior-controlled maturation")
    ax_counts.axhline(24, color="#555555", linestyle=":", linewidth=1)
    ax_counts.set(
        title="Descriptive counts (24 curves per region)",
        ylabel="Session x direction curves",
        xticks=x,
        xticklabels=REGIONS,
        ylim=(0, 25.5),
    )
    ax_counts.legend(frameon=False, fontsize=8)
    ax_counts.spines[["top", "right"]].set_visible(False)

    ax_mouse = figure.add_subplot(right[1, 0])
    fits["early_margin"] = fits["early_five_mean"] - fits[
        ["early_identity_null_q95", "early_spatial_null_q95"]
    ].max(axis=1)
    mouse = fits.groupby(["subject", "region"])["early_margin"].median().unstack()
    region_colors = {"CA3": "#577590", "CA1": "#43AA8B", "RSC": "#F3722C"}
    offsets = {"CA3": -0.16, "CA1": 0.0, "RSC": 0.16}
    for region in REGIONS:
        ax_mouse.scatter(
            np.arange(len(mouse)) + offsets[region],
            mouse[region],
            s=55,
            color=region_colors[region],
            label=region,
            zorder=3,
        )
    ax_mouse.axhline(0, color="#333333", linewidth=1)
    ax_mouse.set(
        title="Mouse-level early identity margin",
        ylabel="Median early score - stricter null q95",
        xlabel="Primary biological replicate",
        xticks=np.arange(len(mouse)),
        xticklabels=mouse.index,
    )
    ax_mouse.legend(frameon=False, ncol=3, fontsize=8)
    ax_mouse.spines[["top", "right"]].set_visible(False)

    figure.suptitle("Late-map identity is often registered early; maturation is heterogeneous", fontsize=14)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", type=Path, default=Path("outputs/maturity/map_maturity_fits.csv"))
    parser.add_argument(
        "--controlled",
        type=Path,
        default=Path("outputs/maturity/behavior_controlled_fits.csv"),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("outputs/maturity/descriptive_summary.json"),
    )
    parser.add_argument(
        "--output-figure",
        type=Path,
        default=Path("outputs/figures/map_maturity_overview.png"),
    )
    args = parser.parse_args()

    fits = pd.read_csv(args.fits)
    controlled = pd.read_csv(args.controlled)
    fits["early_registered"] = (
        fits["early_identity_null_p"].le(0.05) & fits["early_spatial_null_p"].le(0.05)
    )
    fits["early_margin"] = fits["early_five_mean"] - fits[
        ["early_identity_null_q95", "early_spatial_null_q95"]
    ].max(axis=1)
    mouse_margins = (
        fits.groupby(["subject", "region"])["early_margin"].median().reset_index()
    )
    summary = {
        "n_animals": int(fits["subject"].nunique()),
        "n_sessions": int(fits["file"].nunique()),
        "n_curves_per_region": int(len(fits) / len(REGIONS)),
        "early_registered_counts": fits.groupby("region")["early_registered"].sum().astype(int).to_dict(),
        "raw_status_counts": (
            fits.groupby(["region", "fit_status"]).size().rename("n").reset_index().to_dict(orient="records")
        ),
        "behavior_controlled_status_counts": (
            controlled.groupby(["region", "adjusted_fit_status"])
            .size()
            .rename("n")
            .reset_index()
            .to_dict(orient="records")
        ),
        "mouse_median_early_margin": mouse_margins.to_dict(orient="records"),
        "raw_pairing": paired_transition_summary(fits, "fit_status", "midpoint"),
        "behavior_controlled_pairing": paired_transition_summary(
            controlled, "adjusted_fit_status", "adjusted_midpoint"
        ),
        "median_cross_validated_behavior_r2": float(controlled["behavior_cv_r2"].median()),
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    make_figure(fits, controlled, args.output_figure)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
