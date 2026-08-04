#!/usr/bin/env python3
"""Calibrate candidate traversal-window map-maturity scores on synthetic data.

This script does not read the DANDI neural recordings. It simulates a known
transition from an initially unstable spatial map to a stable late map, then
asks whether simple candidate scores recover the imposed midpoint without
acquiring a systematic advantage from larger recorded populations.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


@dataclass(frozen=True)
class SimulationConfig:
    n_traversals: int = 80
    n_position_bins: int = 24
    transition_midpoint: float = 20.0
    transition_scale: float = 4.0
    late_start: int = 55
    rate_scale: float = 1.8
    baseline_rate: float = 0.08


def logistic(x: np.ndarray, midpoint: float, scale: float) -> np.ndarray:
    """Increasing logistic curve constrained to the unit interval."""
    return 1.0 / (1.0 + np.exp(-(x - midpoint) / scale))


def rowwise_correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Correlation after centering every cell's spatial profile separately."""
    a_centered = a - a.mean(axis=1, keepdims=True)
    b_centered = b - b.mean(axis=1, keepdims=True)
    a_flat = a_centered.ravel()
    b_flat = b_centered.ravel()
    denominator = np.linalg.norm(a_flat) * np.linalg.norm(b_flat)
    if denominator == 0:
        return np.nan
    return float(np.dot(a_flat, b_flat) / denominator)


def simulate_session(
    rng: np.random.Generator,
    n_cells: int,
    config: SimulationConfig,
) -> np.ndarray:
    """Return Poisson spike counts shaped traversal x cell x position."""
    position = np.arange(config.n_position_bins)
    centers_late = rng.uniform(0, config.n_position_bins, n_cells)
    centers_early = rng.uniform(0, config.n_position_bins, n_cells)
    widths = rng.uniform(2.0, 4.5, n_cells)
    amplitudes = rng.lognormal(mean=0.0, sigma=0.35, size=n_cells)

    def circular_fields(centers: np.ndarray) -> np.ndarray:
        distance = np.abs(position[None, :] - centers[:, None])
        distance = np.minimum(distance, config.n_position_bins - distance)
        return amplitudes[:, None] * np.exp(-0.5 * (distance / widths[:, None]) ** 2)

    early_map = circular_fields(centers_early)
    late_map = circular_fields(centers_late)
    maturity = logistic(
        np.arange(config.n_traversals),
        config.transition_midpoint,
        config.transition_scale,
    )
    expected = (
        (1.0 - maturity[:, None, None]) * early_map[None, :, :]
        + maturity[:, None, None] * late_map[None, :, :]
        + config.baseline_rate
    )
    return rng.poisson(config.rate_scale * expected).astype(float)


def window_maps(counts: np.ndarray, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Average consecutive traversals and return their center indices."""
    maps = []
    centers = []
    for start in range(0, counts.shape[0] - width + 1):
        stop = start + width
        maps.append(counts[start:stop].mean(axis=0))
        centers.append((start + stop - 1) / 2)
    return np.stack(maps), np.asarray(centers)


def score_maps(
    maps: np.ndarray,
    centers: np.ndarray,
    config: SimulationConfig,
) -> np.ndarray:
    """Compute cell-centered correlation with an independent late template."""
    late = maps[centers >= config.late_start]
    late_a = late[::2].mean(axis=0)
    late_b = late[1::2].mean(axis=0)
    late_template = 0.5 * (late_a + late_b)
    return np.asarray([rowwise_correlation(item, late_template) for item in maps])


def estimate_midpoint(x: np.ndarray, y: np.ndarray) -> float:
    """Fit a four-parameter sigmoid and return its midpoint."""
    finite = np.isfinite(y)
    x_fit = x[finite]
    y_fit = y[finite]
    if len(x_fit) < 8:
        return np.nan

    def sigmoid(t: np.ndarray, low: float, high: float, midpoint: float, scale: float) -> np.ndarray:
        return low + (high - low) * logistic(t, midpoint, scale)

    try:
        parameters, _ = curve_fit(
            sigmoid,
            x_fit,
            y_fit,
            p0=(np.quantile(y_fit, 0.1), np.quantile(y_fit, 0.9), 20.0, 5.0),
            bounds=([-2.0, -2.0, 0.0, 0.5], [2.0, 2.0, 60.0, 30.0]),
            maxfev=20_000,
        )
    except (RuntimeError, ValueError):
        return np.nan
    return float(parameters[2])


def calibrate(
    seed: int,
    replicates: int,
    cell_counts: list[int],
    window_widths: list[int],
    config: SimulationConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    estimates = []
    example_curves = []
    for n_cells in cell_counts:
        for replicate in range(replicates):
            counts = simulate_session(rng, n_cells, config)
            for window_width in window_widths:
                maps, centers = window_maps(counts, window_width)
                values = score_maps(maps, centers, config)
                midpoint = estimate_midpoint(centers, values)
                estimates.append(
                    {
                        "n_cells": n_cells,
                        "replicate": replicate,
                        "window_width": window_width,
                        "estimated_midpoint": midpoint,
                        "error": midpoint - config.transition_midpoint,
                    }
                )
                if replicate == 0:
                    example_curves.extend(
                        {
                            "n_cells": n_cells,
                            "window_width": window_width,
                            "traversal": center,
                            "score": value,
                        }
                        for center, value in zip(centers, values, strict=True)
                    )
    return pd.DataFrame(estimates), pd.DataFrame(example_curves)


def plot_calibration(
    estimates: pd.DataFrame,
    examples: pd.DataFrame,
    config: SimulationConfig,
    output_path: Path,
) -> None:
    window_widths = sorted(int(value) for value in estimates["window_width"].unique())
    palette = ["#277DA1", "#43AA8B", "#F3722C"]
    colors = {width: palette[index % len(palette)] for index, width in enumerate(window_widths)}
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))

    largest = int(examples["n_cells"].max())
    for window_width in window_widths:
        subset = examples[
            (examples["n_cells"] == largest)
            & (examples["window_width"] == window_width)
        ]
        axes[0].plot(
            subset["traversal"],
            subset["score"],
            color=colors[window_width],
            linewidth=2,
            label=f"{window_width}-traversal window",
        )
    axes[0].axvline(config.transition_midpoint, color="#333333", linestyle="--", linewidth=1.5)
    axes[0].set(
        title=f"Synthetic session ({largest} cells)",
        xlabel="Traversal",
        ylabel="Late-map similarity",
    )
    axes[0].legend(frameon=False, fontsize=9)

    offsets = dict(zip(window_widths, np.linspace(-0.3, 0.3, len(window_widths)), strict=True))
    for window_width in window_widths:
        subset = estimates[estimates["window_width"] == window_width]
        summary = subset.groupby("n_cells")["estimated_midpoint"].agg(["mean", "std"]).reset_index()
        axes[1].errorbar(
            summary["n_cells"] + offsets[window_width],
            summary["mean"],
            yerr=summary["std"],
            color=colors[window_width],
            marker="o",
            capsize=3,
            linewidth=1.8,
            label=f"{window_width}-traversal window",
        )
    axes[1].axhline(config.transition_midpoint, color="#333333", linestyle="--", linewidth=1.5)
    axes[1].set(
        title="Recovery of the imposed transition",
        xlabel="Recorded cells",
        ylabel="Estimated midpoint (mean +/- SD)",
        xticks=sorted(estimates["n_cells"].unique()),
    )
    axes[1].legend(frameon=False, fontsize=9)
    figure.suptitle("Estimator calibration without inspecting regional data", fontsize=13, y=0.98)
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--cell-counts", type=int, nargs="+", default=[12, 24, 48, 96])
    parser.add_argument("--window-widths", type=int, nargs="+", default=[1, 5, 9])
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/calibration"))
    args = parser.parse_args()

    config = SimulationConfig()
    estimates, examples = calibrate(
        args.seed, args.replicates, args.cell_counts, args.window_widths, config
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    estimates.to_csv(args.output_dir / "estimator_recovery.csv", index=False)
    plot_calibration(estimates, examples, config, args.output_dir / "estimator_calibration.png")

    summary = (
        estimates.groupby(["window_width", "n_cells"])["error"]
        .agg(n="count", mean_error="mean", median_absolute_error=lambda x: np.median(np.abs(x)), sd_error="std")
        .reset_index()
    )
    payload = {
        "seed": args.seed,
        "replicates": args.replicates,
        "configuration": config.__dict__,
        "summary": summary.to_dict(orient="records"),
    }
    (args.output_dir / "calibration_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
