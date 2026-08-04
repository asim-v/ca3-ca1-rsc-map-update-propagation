#!/usr/bin/env python3
"""Remove nonlinear trial-level measured-behavior predictions from map scores."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from score_map_maturity import fit_maturity_curve


KEYS = ["subject", "file", "region", "block", "direction"]


def behavior_features(table: pd.DataFrame) -> np.ndarray:
    speed = table["moving_speed_mean"].to_numpy(float)
    duration = table["duration_seconds"].to_numpy(float)
    moving = table["moving_sample_fraction"].to_numpy(float)
    return np.column_stack(
        [
            speed,
            speed**2,
            duration,
            duration**2,
            moving,
            moving**2,
            speed * duration,
        ]
    )


def loocv_ridge_predictions(features: np.ndarray, outcome: np.ndarray, alpha: float) -> np.ndarray:
    predictions = np.full(len(outcome), np.nan)
    finite = np.isfinite(outcome) & np.all(np.isfinite(features), axis=1)
    indices = np.flatnonzero(finite)
    for held_out in indices:
        training = indices[indices != held_out]
        x_train = features[training]
        y_train = outcome[training]
        mean_x = x_train.mean(axis=0)
        sd_x = x_train.std(axis=0)
        sd_x[sd_x == 0] = 1.0
        standardized = (x_train - mean_x) / sd_x
        mean_y = float(y_train.mean())
        gram = standardized.T @ standardized + alpha * np.eye(standardized.shape[1])
        beta = np.linalg.solve(gram, standardized.T @ (y_train - mean_y))
        predictions[held_out] = mean_y + ((features[held_out] - mean_x) / sd_x) @ beta
    return predictions


def permutation_trend(curve: np.ndarray, rng: np.random.Generator, n_nulls: int) -> tuple[float, float]:
    finite = np.isfinite(curve)
    ranks = rankdata(curve[finite])
    traversal = np.arange(len(ranks), dtype=float)
    rho = float(np.corrcoef(traversal, ranks)[0, 1])
    null = np.asarray(
        [np.corrcoef(traversal, rng.permutation(ranks))[0, 1] for _ in range(n_nulls)]
    )
    p_value = float((1 + np.sum(null >= rho)) / (len(null) + 1))
    return rho, p_value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scores", type=Path, default=Path("outputs/maturity/map_maturity_scores.csv")
    )
    parser.add_argument(
        "--fits", type=Path, default=Path("outputs/maturity/map_maturity_fits.csv")
    )
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--n-nulls", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument(
        "--output", type=Path, default=Path("outputs/maturity/behavior_controlled_fits.csv")
    )
    args = parser.parse_args()

    scores = pd.read_csv(args.scores)
    raw_fits = pd.read_csv(args.fits).set_index(KEYS)
    rng = np.random.default_rng(args.seed)
    rows = []
    for key, group in scores.groupby(KEYS, sort=True):
        group = group.sort_values("direction_traversal")
        outcome = group["map_similarity"].to_numpy(float)
        predictions = loocv_ridge_predictions(behavior_features(group), outcome, args.alpha)
        mean_outcome = float(np.nanmean(outcome))
        adjusted = outcome - predictions + mean_outcome
        finite = np.isfinite(outcome) & np.isfinite(predictions)
        denominator = float(np.sum((outcome[finite] - np.mean(outcome[finite])) ** 2))
        behavior_cv_r2 = (
            1.0 - float(np.sum((outcome[finite] - predictions[finite]) ** 2)) / denominator
            if denominator > 0
            else np.nan
        )
        trend_rho, trend_p = permutation_trend(adjusted, rng, args.n_nulls)
        raw = raw_fits.loc[key]
        early_correction = float(np.nanmean(predictions[:5]) - mean_outcome)
        adjusted_fit = fit_maturity_curve(
            group["direction_traversal"].to_numpy(float),
            adjusted,
            float(raw["early_identity_null_q95"] - early_correction),
            float(raw["early_identity_null_p"]),
            float(raw["early_spatial_null_q95"] - early_correction),
            float(raw["early_spatial_null_p"]),
            trend_rho,
            trend_p,
        )
        rows.append(
            {
                **dict(zip(KEYS, key, strict=True)),
                "raw_fit_status": raw["fit_status"],
                "raw_midpoint": raw.get("midpoint", np.nan),
                "behavior_cv_r2": behavior_cv_r2,
                "behavior_ridge_alpha": args.alpha,
                **{f"adjusted_{name}": value for name, value in adjusted_fit.items()},
            }
        )

    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    summary = {
        "parameters": {"alpha": args.alpha, "n_nulls": args.n_nulls, "seed": args.seed},
        "n_curves": int(len(output)),
        "cross_validated_behavior_r2_median": float(output["behavior_cv_r2"].median()),
        "status_counts": (
            output.groupby(["region", "adjusted_fit_status"])
            .size()
            .rename("n")
            .reset_index()
            .to_dict(orient="records")
        ),
        "raw_transitions": int((output["raw_fit_status"] == "detectable_transition").sum()),
        "behavior_adjusted_transitions": int(
            (output["adjusted_fit_status"] == "detectable_transition").sum()
        ),
    }
    args.output.with_suffix(".summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
