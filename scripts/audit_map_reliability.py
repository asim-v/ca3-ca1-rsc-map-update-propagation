#!/usr/bin/env python3
"""Audit late-map occupancy and split-half reliability without onset analysis.

The script intentionally never scores early traversals and never compares
regional onset. It supplies the final quality gate needed before inspecting the
CA3-versus-CA1-versus-RSC ordering hypothesis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO

from extract_traversals import behavior_blocks, endpoint_events, linearize_position


REGION_RULES = {
    "CA1": lambda table: (table["cell_area"] == "CA1") & (table["cell_type"] == "Pyramidal Cell"),
    "CA3": lambda table: (table["cell_area"] == "CA3") & (table["cell_type"] == "Pyramidal Cell"),
    "RSC": lambda table: (table["cell_area"] == "RSC") & (table["cell_type"] != "Narrow Interneuron"),
}


def detect_traversals(
    position: np.ndarray,
    timestamps: np.ndarray,
    low: float,
    high: float,
    min_coverage: float,
    max_duration: float,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for block_number, (block_start, block_stop) in enumerate(behavior_blocks(timestamps), start=1):
        block_position = position[block_start : block_stop + 1]
        lo, hi = np.nanpercentile(block_position, [1, 99])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            continue
        normalized = np.clip((block_position - lo) / (hi - lo), 0, 1)
        events = endpoint_events(normalized, low, high)
        trial_number = 0
        for (start_local, side_a), (stop_local, side_b) in zip(events[:-1], events[1:]):
            if side_a == side_b:
                continue
            start = block_start + start_local
            stop = block_start + stop_local
            duration = float(timestamps[stop] - timestamps[start])
            coverage = float(
                np.nanmax(normalized[start_local : stop_local + 1])
                - np.nanmin(normalized[start_local : stop_local + 1])
            )
            if duration <= 0 or duration > max_duration or coverage < min_coverage:
                continue
            trial_number += 1
            rows.append(
                {
                    "block": block_number,
                    "trial": trial_number,
                    "direction": f"{side_a}to{side_b}",
                    "start_sample": start,
                    "stop_sample": stop,
                    "start_time": float(timestamps[start]),
                    "stop_time": float(timestamps[stop]),
                    "position_lo": float(lo),
                    "position_hi": float(hi),
                }
            )
    return rows


def trial_counts_and_occupancy(
    trial: dict[str, float | int | str],
    spike_times: list[np.ndarray],
    position: np.ndarray,
    timestamps: np.ndarray,
    speed: np.ndarray,
    n_bins: int,
    speed_threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    start = int(trial["start_sample"])
    stop = int(trial["stop_sample"]) + 1
    sample_time = timestamps[start:stop]
    sample_speed = speed[start:stop]
    normalized = np.clip(
        (position[start:stop] - float(trial["position_lo"]))
        / (float(trial["position_hi"]) - float(trial["position_lo"])),
        0,
        1,
    )
    moving = np.isfinite(normalized) & np.isfinite(sample_speed) & (sample_speed > speed_threshold)
    bin_edges = np.linspace(0, 1, n_bins + 1)
    dt = float(np.median(np.diff(sample_time))) if len(sample_time) > 1 else 0.0
    occupancy = np.histogram(normalized[moving], bins=bin_edges)[0].astype(float) * dt
    counts = np.zeros((len(spike_times), n_bins), dtype=float)
    for cell_index, cell_spikes in enumerate(spike_times):
        lo = np.searchsorted(cell_spikes, float(trial["start_time"]), side="left")
        hi = np.searchsorted(cell_spikes, float(trial["stop_time"]), side="right")
        selected = cell_spikes[lo:hi]
        if len(selected) == 0:
            continue
        spike_position = np.interp(selected, sample_time, normalized)
        spike_speed = np.interp(selected, sample_time, sample_speed)
        valid = np.isfinite(spike_position) & np.isfinite(spike_speed) & (spike_speed > speed_threshold)
        counts[cell_index] = np.histogram(spike_position[valid], bins=bin_edges)[0]
    return counts, occupancy


def aggregate_rate(counts: np.ndarray, occupancy: np.ndarray, indices: np.ndarray) -> np.ndarray:
    total_counts = counts[indices].sum(axis=0)
    total_occupancy = occupancy[indices].sum(axis=0)
    rate = np.full(total_counts.shape, np.nan, dtype=float)
    occupied = total_occupancy > 0
    rate[:, occupied] = total_counts[:, occupied] / total_occupancy[occupied]
    return np.sqrt(rate)


def centered_map_correlation(a: np.ndarray, b: np.ndarray) -> float:
    valid_cells = np.isfinite(a).any(axis=1) & np.isfinite(b).any(axis=1)
    if not np.any(valid_cells):
        return np.nan
    a = a[valid_cells]
    b = b[valid_cells]
    a = a - np.nanmean(a, axis=1, keepdims=True)
    b = b - np.nanmean(b, axis=1, keepdims=True)
    valid = np.isfinite(a) & np.isfinite(b)
    if valid.sum() < 3:
        return np.nan
    x = a[valid]
    y = b[valid]
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / denominator) if denominator > 0 else np.nan


def audit_file(path: Path, n_bins: int, late_traversals: int, speed_threshold: float) -> list[dict]:
    with NWBHDF5IO(str(path), "r", load_namespaces=True) as io:
        nwb = io.read()
        behavior = nwb.processing["behavior"]
        position_series = behavior["AnimalPosition"]["Position"]
        position = linearize_position(np.asarray(position_series.data[:]).squeeze())
        timestamps = np.asarray(position_series.timestamps[:]).squeeze()
        speed = np.asarray(behavior["Speed"].data[:]).squeeze()
        units = nwb.units.to_dataframe()
        traversals = detect_traversals(position, timestamps, 0.15, 0.85, 0.70, 120.0)

        rows = []
        for region, selector in REGION_RULES.items():
            region_units = units.loc[selector(units)]
            spike_times = [np.asarray(value) for value in region_units["spike_times"]]
            all_counts = []
            all_occupancy = []
            for trial in traversals:
                counts, occupancy = trial_counts_and_occupancy(
                    trial,
                    spike_times,
                    position,
                    timestamps,
                    speed,
                    n_bins,
                    speed_threshold,
                )
                all_counts.append(counts)
                all_occupancy.append(occupancy)
            counts_array = np.stack(all_counts)
            occupancy_array = np.stack(all_occupancy)

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
                        rows.append(
                            {
                                "file": path.name,
                                "region": region,
                                "block": block,
                                "direction": direction,
                                "n_cells": len(region_units),
                                "n_traversals": len(matching),
                                "late_traversals": 0,
                                "occupied_bin_fraction": np.nan,
                                "split_half_r": np.nan,
                                "spearman_brown_reliability": np.nan,
                                "status": "insufficient_late_traversals",
                            }
                        )
                        continue
                    late = matching[-late_traversals:]
                    odd = late[::2]
                    even = late[1::2]
                    map_odd = aggregate_rate(counts_array, occupancy_array, odd)
                    map_even = aggregate_rate(counts_array, occupancy_array, even)
                    correlation = centered_map_correlation(map_odd, map_even)
                    corrected = 2 * correlation / (1 + correlation) if np.isfinite(correlation) and correlation > -1 else np.nan
                    occupied = occupancy_array[late].sum(axis=0) > 0
                    rows.append(
                        {
                            "file": path.name,
                            "region": region,
                            "block": block,
                            "direction": direction,
                            "n_cells": len(region_units),
                            "n_traversals": len(matching),
                            "late_traversals": len(late),
                            "occupied_bin_fraction": float(occupied.mean()),
                            "split_half_r": correlation,
                            "spearman_brown_reliability": corrected,
                            "status": "ok",
                        }
                    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nwb", type=Path, nargs="+")
    parser.add_argument("--n-bins", type=int, default=24)
    parser.add_argument("--late-traversals", type=int, default=10)
    parser.add_argument("--speed-threshold", type=float, default=2.5)
    parser.add_argument("--output", type=Path, default=Path("outputs/quality/map_reliability.csv"))
    args = parser.parse_args()

    rows = []
    for path in args.nwb:
        print(f"Auditing {path.name}", flush=True)
        rows.extend(audit_file(path, args.n_bins, args.late_traversals, args.speed_threshold))
    table = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    summary = {
        "n_files": int(table["file"].nunique()),
        "n_rows": int(len(table)),
        "parameters": {
            "n_bins": args.n_bins,
            "late_traversals": args.late_traversals,
            "speed_threshold": args.speed_threshold,
        },
        "status_counts": table["status"].value_counts().to_dict(),
        "occupied_bin_fraction_minimum": float(table["occupied_bin_fraction"].min()),
        "reliability_by_region": (
            table.groupby("region")["spearman_brown_reliability"]
            .agg(["count", "min", "median", "max"])
            .reset_index()
            .to_dict(orient="records")
        ),
    }
    summary_path = args.output.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
