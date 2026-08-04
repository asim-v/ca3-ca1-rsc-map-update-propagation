#!/usr/bin/env python3
"""Score traversal-resolved similarity to a held-out late spatial map.

The implementation follows the locked endpoint in docs/estimator_calibration.md.
It analyzes each region and running direction separately, equalizes cell count
within session, and does not force a midpoint when an increasing transition is
not supported.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO
from scipy.optimize import curve_fit

from audit_map_reliability import (
    REGION_RULES,
    aggregate_rate,
    centered_map_correlation,
    detect_traversals,
    trial_counts_and_occupancy,
)
from extract_traversals import linearize_position


def increasing_sigmoid(
    x: np.ndarray,
    low: float,
    amplitude: float,
    midpoint: float,
    scale: float,
) -> np.ndarray:
    return low + amplitude / (1.0 + np.exp(-(x - midpoint) / scale))


def fit_maturity_curve(
    x: np.ndarray,
    y: np.ndarray,
    early_null_q95: float,
    early_null_p: float,
) -> dict[str, float | str | bool]:
    finite = np.isfinite(x) & np.isfinite(y)
    x = x[finite]
    y = y[finite]
    if len(x) < 15:
        return {"fit_status": "insufficient_scores", "midpoint_usable": False}

    constant_rss = float(np.sum((y - y.mean()) ** 2))
    try:
        parameters, _ = curve_fit(
            increasing_sigmoid,
            x,
            y,
            p0=(
                float(np.quantile(y, 0.10)),
                float(max(np.quantile(y, 0.90) - np.quantile(y, 0.10), 0.10)),
                float(np.median(x)),
                5.0,
            ),
            bounds=(
                [-1.0, 0.0, float(x.min()), 0.5],
                [1.0, 2.0, float(x.max()), 30.0],
            ),
            maxfev=30_000,
        )
        fitted = increasing_sigmoid(x, *parameters)
    except (RuntimeError, ValueError):
        return {"fit_status": "fit_failed", "midpoint_usable": False}

    sigmoid_rss = float(np.sum((y - fitted) ** 2))
    epsilon = np.finfo(float).tiny
    constant_aic = len(x) * np.log(max(constant_rss / len(x), epsilon)) + 2
    sigmoid_aic = len(x) * np.log(max(sigmoid_rss / len(x), epsilon)) + 8
    low, amplitude, midpoint, scale = (float(value) for value in parameters)
    r_squared = 1.0 - sigmoid_rss / constant_rss if constant_rss > 0 else np.nan
    midpoint_interior = (midpoint > x.min() + 1) and (midpoint < x.max() - 1)
    midpoint_usable = bool(amplitude >= 0.10 and sigmoid_aic <= constant_aic - 2 and midpoint_interior)

    early_mean = float(np.mean(y[x <= np.sort(x)[min(4, len(x) - 1)]]))
    late_mean = float(np.mean(y[x >= np.sort(x)[max(0, len(x) - 5)]]))
    if midpoint_usable:
        status = "detectable_transition"
    elif early_mean >= late_mean - 0.10 and early_mean > early_null_q95:
        status = "present_from_first_observed_traversals"
    else:
        status = "unresolved_nonmonotonic_or_weak"

    return {
        "fit_status": status,
        "midpoint_usable": midpoint_usable,
        "low": low,
        "amplitude": amplitude,
        "midpoint": midpoint,
        "scale": scale,
        "r_squared": float(r_squared),
        "constant_aic": float(constant_aic),
        "sigmoid_aic": float(sigmoid_aic),
        "delta_aic_sigmoid_minus_constant": float(sigmoid_aic - constant_aic),
        "early_five_mean": early_mean,
        "late_five_mean": late_mean,
        "late_minus_early": late_mean - early_mean,
        "early_identity_null_q95": early_null_q95,
        "early_identity_null_p": early_null_p,
        "n_scores": int(len(x)),
    }


def score_subset(map_a: np.ndarray, map_b: np.ndarray, cell_indices: np.ndarray) -> float:
    return centered_map_correlation(map_a[cell_indices], map_b[cell_indices])


def analyze_session(
    path: Path,
    seed: int,
    n_bins: int,
    late_traversals: int,
    speed_threshold: float,
    n_subsamples: int,
    n_nulls: int,
) -> tuple[list[dict], list[dict]]:
    with NWBHDF5IO(str(path), "r", load_namespaces=True) as io:
        nwb = io.read()
        subject = nwb.subject.subject_id if nwb.subject is not None else path.name.split("_")[0]
        behavior = nwb.processing["behavior"]
        position_series = behavior["AnimalPosition"]["Position"]
        position = linearize_position(np.asarray(position_series.data[:]).squeeze())
        timestamps = np.asarray(position_series.timestamps[:]).squeeze()
        speed = np.asarray(behavior["Speed"].data[:]).squeeze()
        units = nwb.units.to_dataframe()
        traversals = detect_traversals(position, timestamps, 0.15, 0.85, 0.70, 120.0)

        region_units = {region: units.loc[selector(units)] for region, selector in REGION_RULES.items()}
        n_equal = min(len(table) for table in region_units.values())
        rng = np.random.default_rng(seed)
        subsets = {
            region: np.stack(
                [rng.choice(len(table), size=n_equal, replace=False) for _ in range(n_subsamples)]
            )
            for region, table in region_units.items()
        }

        score_rows: list[dict] = []
        fit_rows: list[dict] = []
        for region, table in region_units.items():
            spike_times = [np.asarray(value) for value in table["spike_times"]]
            trial_counts = []
            trial_occupancy = []
            for traversal in traversals:
                counts, occupancy = trial_counts_and_occupancy(
                    traversal,
                    spike_times,
                    position,
                    timestamps,
                    speed,
                    n_bins,
                    speed_threshold,
                )
                trial_counts.append(counts)
                trial_occupancy.append(occupancy)
            counts_array = np.stack(trial_counts)
            occupancy_array = np.stack(trial_occupancy)
            single_maps = np.stack(
                [
                    aggregate_rate(counts_array, occupancy_array, np.asarray([index]))
                    for index in range(len(traversals))
                ]
            )

            for block in sorted({int(item["block"]) for item in traversals}):
                for direction in ("LtoR", "RtoL"):
                    matching = np.asarray(
                        [
                            index
                            for index, item in enumerate(traversals)
                            if int(item["block"]) == block and item["direction"] == direction
                        ],
                        dtype=int,
                    )
                    if len(matching) < late_traversals:
                        continue
                    late = matching[-late_traversals:]
                    full_reference = aggregate_rate(counts_array, occupancy_array, late)
                    direction_scores = []
                    for ordinal, trial_index in enumerate(matching, start=1):
                        reference_indices = late[late != trial_index] if trial_index in late else late
                        reference = aggregate_rate(counts_array, occupancy_array, reference_indices)
                        values = np.asarray(
                            [
                                score_subset(single_maps[trial_index], reference, subset)
                                for subset in subsets[region]
                            ]
                        )
                        finite_values = values[np.isfinite(values)]
                        score = float(np.median(finite_values)) if len(finite_values) else np.nan
                        direction_scores.append(score)
                        score_rows.append(
                            {
                                "subject": subject,
                                "file": path.name,
                                "region": region,
                                "block": block,
                                "direction": direction,
                                "direction_traversal": ordinal,
                                "global_traversal": int(trial_index + 1),
                                "map_similarity": score,
                                "subsample_q025": float(np.quantile(finite_values, 0.025)) if len(finite_values) else np.nan,
                                "subsample_q975": float(np.quantile(finite_values, 0.975)) if len(finite_values) else np.nan,
                                "is_late_reference_traversal": bool(trial_index in late),
                                "n_cells_available": len(table),
                                "n_cells_equalized": n_equal,
                            }
                        )
                    observed_early = float(np.nanmean(direction_scores[:5]))
                    null_early = []
                    for _ in range(n_nulls):
                        subset = rng.choice(len(table), size=n_equal, replace=False)
                        permuted_reference = rng.permutation(subset)
                        null_values = [
                            centered_map_correlation(
                                single_maps[trial_index][subset],
                                full_reference[permuted_reference],
                            )
                            for trial_index in matching[:5]
                        ]
                        null_early.append(float(np.nanmean(null_values)))
                    finite_null = np.asarray(null_early)[np.isfinite(null_early)]
                    null_q95 = float(np.quantile(finite_null, 0.95))
                    null_p = float((1 + np.sum(finite_null >= observed_early)) / (len(finite_null) + 1))
                    fit = fit_maturity_curve(
                        np.arange(1, len(direction_scores) + 1, dtype=float),
                        np.asarray(direction_scores),
                        null_q95,
                        null_p,
                    )
                    fit_rows.append(
                        {
                            "subject": subject,
                            "file": path.name,
                            "region": region,
                            "block": block,
                            "direction": direction,
                            "n_cells_available": len(table),
                            "n_cells_equalized": n_equal,
                            **fit,
                        }
                    )
    return score_rows, fit_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nwb", type=Path, nargs="+")
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--n-bins", type=int, default=24)
    parser.add_argument("--late-traversals", type=int, default=10)
    parser.add_argument("--speed-threshold", type=float, default=2.5)
    parser.add_argument("--n-subsamples", type=int, default=200)
    parser.add_argument("--n-nulls", type=int, default=500)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/maturity"))
    args = parser.parse_args()

    all_scores = []
    all_fits = []
    for file_number, path in enumerate(args.nwb):
        print(f"Scoring {path.name}", flush=True)
        score_rows, fit_rows = analyze_session(
            path,
            args.seed + file_number,
            args.n_bins,
            args.late_traversals,
            args.speed_threshold,
            args.n_subsamples,
            args.n_nulls,
        )
        all_scores.extend(score_rows)
        all_fits.extend(fit_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores = pd.DataFrame(all_scores)
    fits = pd.DataFrame(all_fits)
    scores.to_csv(args.output_dir / "map_maturity_scores.csv", index=False)
    fits.to_csv(args.output_dir / "map_maturity_fits.csv", index=False)
    summary = {
        "parameters": {
            "seed": args.seed,
            "n_bins": args.n_bins,
            "late_traversals": args.late_traversals,
            "speed_threshold": args.speed_threshold,
            "n_subsamples": args.n_subsamples,
            "n_nulls": args.n_nulls,
        },
        "n_files": int(scores["file"].nunique()),
        "n_subjects": int(scores["subject"].nunique()),
        "n_score_rows": int(len(scores)),
        "fit_status_counts": fits.groupby(["region", "fit_status"]).size().rename("n").reset_index().to_dict(orient="records"),
    }
    (args.output_dir / "map_maturity_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
