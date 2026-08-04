#!/usr/bin/env python3
"""Test traversal-scale cross-area prediction with circular-shift nulls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


CURVE_KEYS = ["subject", "file", "block", "direction"]
ROW_KEYS = CURVE_KEYS + ["direction_traversal", "global_traversal"]
CONNECTIONS = [
    ("CA3", "CA1", "forward"),
    ("CA1", "CA3", "reverse"),
    ("CA1", "RSC", "forward"),
    ("RSC", "CA1", "reverse"),
]


def loocv_ridge_mse(features: np.ndarray, outcome: np.ndarray, alpha: float) -> float:
    predictions = np.full(len(outcome), np.nan)
    for held_out in range(len(outcome)):
        training = np.arange(len(outcome)) != held_out
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
    return float(np.mean((outcome - predictions) ** 2))


def base_features(table: pd.DataFrame, target_previous: np.ndarray) -> np.ndarray:
    current = table.iloc[1:]
    speed = current["moving_speed_mean"].to_numpy(float)
    duration = current["duration_seconds"].to_numpy(float)
    moving = current["moving_sample_fraction"].to_numpy(float)
    traversal = current["direction_traversal"].to_numpy(float)
    return np.column_stack(
        [
            target_previous,
            speed,
            speed**2,
            duration,
            duration**2,
            moving,
            moving**2,
            speed * duration,
            traversal,
            traversal**2,
        ]
    )


def analyze_curve(
    table: pd.DataFrame,
    source: str,
    target: str,
    lag: int,
    alpha: float,
) -> tuple[dict, np.ndarray] | None:
    table = table.sort_values("direction_traversal")
    source_values = table[source].to_numpy(float)
    target_values = table[target].to_numpy(float)
    if len(table) < 18:
        return None
    target_now = target_values[1:]
    target_previous = target_values[:-1]
    source_predictor = source_values[1:] if lag == 0 else source_values[:-1]
    reduced = base_features(table, target_previous)
    finite = (
        np.isfinite(target_now)
        & np.isfinite(source_predictor)
        & np.all(np.isfinite(reduced), axis=1)
    )
    target_now = target_now[finite]
    source_predictor = source_predictor[finite]
    reduced = reduced[finite]
    if len(target_now) < 15 or np.std(target_now) == 0:
        return None
    target_now = (target_now - np.mean(target_now)) / np.std(target_now)
    reduced_mse = loocv_ridge_mse(reduced, target_now, alpha)
    full = np.column_stack([reduced, source_predictor])
    full_mse = loocv_ridge_mse(full, target_now, alpha)
    observed = reduced_mse - full_mse

    null_effects = []
    for shift in range(2, len(table) - 1):
        shifted = np.roll(source_values, shift)
        shifted_predictor = shifted[1:] if lag == 0 else shifted[:-1]
        shifted_predictor = shifted_predictor[finite]
        null_full = np.column_stack([reduced, shifted_predictor])
        null_effects.append(reduced_mse - loocv_ridge_mse(null_full, target_now, alpha))
    null_array = np.asarray(null_effects)
    circular_p = float((1 + np.sum(null_array >= observed)) / (len(null_array) + 1))
    return (
        {
            "n_traversals_modeled": int(len(target_now)),
            "reduced_loocv_mse": reduced_mse,
            "full_loocv_mse": full_mse,
            "delta_mse_source_added": observed,
            "curve_circular_shift_p": circular_p,
            "n_unique_circular_shifts": int(len(null_array)),
        },
        null_array,
    )


def mouse_level_null(
    observed: pd.DataFrame,
    nulls: dict[int, np.ndarray],
    rng: np.random.Generator,
    iterations: int,
) -> tuple[pd.DataFrame, float, float]:
    mouse_observed = (
        observed.groupby("subject")["delta_mse_source_added"].median().rename("mouse_median_delta_mse").reset_index()
    )
    grand_observed = float(mouse_observed["mouse_median_delta_mse"].mean())
    null_grand = []
    subject_rows = {
        subject: group.index.to_numpy() for subject, group in observed.groupby("subject")
    }
    for _ in range(iterations):
        subject_values = []
        for indices in subject_rows.values():
            sampled = [float(rng.choice(nulls[int(index)])) for index in indices]
            subject_values.append(float(np.median(sampled)))
        null_grand.append(float(np.mean(subject_values)))
    null_grand = np.asarray(null_grand)
    p_value = float((1 + np.sum(null_grand >= grand_observed)) / (len(null_grand) + 1))
    return mouse_observed, grand_observed, p_value


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scores", type=Path, default=Path("outputs/maturity/map_maturity_scores.csv")
    )
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--mouse-null-iterations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/cross_area"))
    args = parser.parse_args()

    scores = pd.read_csv(args.scores)
    wide = scores.pivot(index=ROW_KEYS, columns="region", values="map_similarity").reset_index()
    behavior_columns = [
        "moving_speed_mean",
        "duration_seconds",
        "moving_sample_fraction",
    ]
    behavior = scores[ROW_KEYS + behavior_columns].drop_duplicates(ROW_KEYS)
    wide = wide.merge(behavior, on=ROW_KEYS, validate="one_to_one")

    rows = []
    nulls: dict[int, np.ndarray] = {}
    for key, curve in wide.groupby(CURVE_KEYS, sort=True):
        for source, target, anatomical_direction in CONNECTIONS:
            for lag in (0, 1):
                result = analyze_curve(curve, source, target, lag, args.alpha)
                if result is None:
                    continue
                values, curve_null = result
                row_id = len(rows)
                rows.append(
                    {
                        "row_id": row_id,
                        **dict(zip(CURVE_KEYS, key, strict=True)),
                        "source": source,
                        "target": target,
                        "anatomical_direction": anatomical_direction,
                        "lag_same_direction_traversals": lag,
                        **values,
                    }
                )
                nulls[row_id] = curve_null

    curve_results = pd.DataFrame(rows).set_index("row_id", drop=False)
    rng = np.random.default_rng(args.seed)
    summaries = []
    mouse_rows = []
    for (source, target, lag), subset in curve_results.groupby(
        ["source", "target", "lag_same_direction_traversals"], sort=True
    ):
        mouse, grand_effect, p_value = mouse_level_null(
            subset, nulls, rng, args.mouse_null_iterations
        )
        mouse["source"] = source
        mouse["target"] = target
        mouse["lag_same_direction_traversals"] = lag
        mouse_rows.append(mouse)
        summaries.append(
            {
                "source": source,
                "target": target,
                "lag_same_direction_traversals": int(lag),
                "n_curves": int(len(subset)),
                "n_mice": int(subset["subject"].nunique()),
                "mouse_level_mean_delta_mse": grand_effect,
                "mouse_level_circular_shift_p": p_value,
                "mice_positive": int((mouse["mouse_median_delta_mse"] > 0).sum()),
            }
        )
    summary = pd.DataFrame(summaries)
    summary["bh_q_across_8_tests"] = benjamini_hochberg(
        summary["mouse_level_circular_shift_p"].to_numpy(float)
    )
    mouse_results = pd.concat(mouse_rows, ignore_index=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    curve_results.to_csv(args.output_dir / "curve_effects.csv", index=False)
    mouse_results.to_csv(args.output_dir / "mouse_effects.csv", index=False)
    summary.to_csv(args.output_dir / "summary.csv", index=False)
    payload = {
        "parameters": {
            "alpha": args.alpha,
            "mouse_null_iterations": args.mouse_null_iterations,
            "seed": args.seed,
            "lag_definition": "one same-direction traversal, approximately two physical traversals",
        },
        "results": summary.to_dict(orient="records"),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
