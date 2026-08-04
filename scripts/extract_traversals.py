#!/usr/bin/env python3
"""Extract alternating linear-track traversals from behavior only.

This script intentionally does not load spikes. Its output is a development
ledger for choosing defensible traversal rules before inspecting neural map
maturity.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO


def behavior_blocks(timestamps: np.ndarray) -> list[tuple[int, int]]:
    delta = np.diff(timestamps)
    threshold = max(1.0, float(np.median(delta) * 10))
    split_after = np.flatnonzero(delta > threshold)
    starts = np.r_[0, split_after + 1]
    stops = np.r_[split_after, len(timestamps) - 1]
    return [(int(start), int(stop)) for start, stop in zip(starts, stops)]


def linearize_position(position: np.ndarray) -> np.ndarray:
    position = np.asarray(position).squeeze()
    if position.ndim == 1:
        return position
    if position.ndim != 2:
        raise ValueError(f"Unsupported position shape: {position.shape}")
    valid = np.all(np.isfinite(position), axis=1)
    center = np.nanmean(position[valid], axis=0)
    _, _, vh = np.linalg.svd(position[valid] - center, full_matrices=False)
    return (position - center) @ vh[0]


def endpoint_events(normalized: np.ndarray, low: float, high: float) -> list[tuple[int, str]]:
    events: list[tuple[int, str]] = []
    state: str | None = None
    for index, value in enumerate(normalized):
        new_state = "L" if value <= low else "R" if value >= high else None
        if new_state is not None and new_state != state:
            events.append((index, new_state))
            state = new_state
    return events


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("nwb", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/traversal_audit"))
    parser.add_argument("--low", type=float, default=0.15)
    parser.add_argument("--high", type=float, default=0.85)
    parser.add_argument("--min-coverage", type=float, default=0.70)
    parser.add_argument("--max-duration", type=float, default=120.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with NWBHDF5IO(str(args.nwb), "r", load_namespaces=True) as io:
        nwb = io.read()
        behavior = nwb.processing["behavior"]
        series = behavior["AnimalPosition"]["Position"]
        position_xy = np.asarray(series.data[:]).squeeze()
        position = linearize_position(position_xy)
        timestamps = np.asarray(series.timestamps[:]).squeeze()
        speed = np.asarray(behavior["Speed"].data[:]).squeeze()

    rows: list[dict[str, float | int | str]] = []
    for block_number, (block_start, block_stop) in enumerate(behavior_blocks(timestamps), start=1):
        pos = position[block_start : block_stop + 1]
        lo, hi = np.nanpercentile(pos, [1, 99])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            continue
        normalized = np.clip((pos - lo) / (hi - lo), 0, 1)
        events = endpoint_events(normalized, args.low, args.high)
        trial_number = 0
        for (start_local, side_a), (stop_local, side_b) in zip(events[:-1], events[1:]):
            if side_a == side_b:
                continue
            start = block_start + start_local
            stop = block_start + stop_local
            duration = float(timestamps[stop] - timestamps[start])
            coverage = float(np.nanmax(normalized[start_local : stop_local + 1]) - np.nanmin(normalized[start_local : stop_local + 1]))
            if duration <= 0 or duration > args.max_duration or coverage < args.min_coverage:
                continue
            trial_number += 1
            rows.append(
                {
                    "file": args.nwb.name,
                    "block": block_number,
                    "trial": trial_number,
                    "direction": f"{side_a}to{side_b}",
                    "start_sample": start,
                    "stop_sample": stop,
                    "start_time": float(timestamps[start]),
                    "stop_time": float(timestamps[stop]),
                    "duration_seconds": duration,
                    "normalized_coverage": coverage,
                    "samples": int(stop - start + 1),
                    "moving_fraction": float(np.mean(speed[start : stop + 1] > 2.5)),
                }
            )

    table = pd.DataFrame(rows)
    destination = args.output_dir / f"{args.nwb.stem}_traversals.csv"
    table.to_csv(destination, index=False)

    relative_minutes = (timestamps - timestamps[0]) / 60
    fig, ax = plt.subplots(figsize=(11, 3.8), constrained_layout=True)
    ax.plot(relative_minutes, position, color="#355C7D", linewidth=0.5)
    colors = {"LtoR": "#F67280", "RtoL": "#6C5B7B"}
    for row in rows:
        start_min = (float(row["start_time"]) - timestamps[0]) / 60
        stop_min = (float(row["stop_time"]) - timestamps[0]) / 60
        ax.axvspan(start_min, stop_min, color=colors[str(row["direction"])], alpha=0.12)
    ax.set(title=args.nwb.name, xlabel="Minutes", ylabel="Linearized position")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(args.output_dir / f"{args.nwb.stem}_traversals.png", dpi=180)
    plt.close(fig)

    if table.empty:
        print("No traversals detected")
    else:
        print(table.groupby(["block", "direction"]).size().rename("n").to_string())
        print(f"\nTotal traversals: {len(table)}")


if __name__ == "__main__":
    main()
