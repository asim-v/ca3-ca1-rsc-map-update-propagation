#!/usr/bin/env python3
"""Audit behavioral time coverage without evaluating the neural hypothesis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pynwb import NWBHDF5IO


def as_array(value) -> np.ndarray:
    return np.asarray(value[:])


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("nwb", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/navigation_audit"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with NWBHDF5IO(str(args.nwb), "r", load_namespaces=True) as io:
        nwb = io.read()
        behavior = nwb.processing["behavior"]
        position_series = behavior["AnimalPosition"]["Position"]
        speed_series = behavior["Speed"]
        position_xy = as_array(position_series.data).squeeze()
        position = linearize_position(position_xy)
        timestamps = as_array(position_series.timestamps).squeeze()
        speed = as_array(speed_series.data).squeeze()

    if timestamps.ndim != 1 or len(timestamps) < 2:
        raise SystemExit("Position timestamps are missing or malformed")
    if len(position) != len(timestamps) or len(speed) != len(timestamps):
        raise SystemExit(
            f"Behavior length mismatch: position={len(position)}, speed={len(speed)}, "
            f"timestamps={len(timestamps)}"
        )

    delta = np.diff(timestamps)
    order = np.argsort(delta)[::-1]
    top_gaps = [
        {
            "after_sample": int(i),
            "gap_seconds": float(delta[i]),
            "time_before": float(timestamps[i]),
            "time_after": float(timestamps[i + 1]),
        }
        for i in order[: min(10, len(order))]
    ]
    duration_minutes = float((timestamps[-1] - timestamps[0]) / 60)
    gap_threshold = max(1.0, float(np.median(delta) * 10))
    split_after = np.flatnonzero(delta > gap_threshold)
    starts = np.r_[0, split_after + 1]
    stops = np.r_[split_after, len(timestamps) - 1]
    blocks = [
        {
            "block": int(index + 1),
            "start_sample": int(start),
            "stop_sample": int(stop),
            "start_time": float(timestamps[start]),
            "stop_time": float(timestamps[stop]),
            "duration_minutes": float((timestamps[stop] - timestamps[start]) / 60),
            "samples": int(stop - start + 1),
        }
        for index, (start, stop) in enumerate(zip(starts, stops))
    ]
    two_maze_candidate = len(blocks) == 2 and all(x["duration_minutes"] >= 10 for x in blocks)
    summary = {
        "file": args.nwb.name,
        "samples": int(len(timestamps)),
        "timestamp_start": float(timestamps[0]),
        "timestamp_stop": float(timestamps[-1]),
        "duration_minutes": duration_minutes,
        "median_sample_interval_seconds": float(np.median(delta)),
        "maximum_gap_seconds": float(np.max(delta)),
        "position_min": float(np.nanmin(position)),
        "position_max": float(np.nanmax(position)),
        "speed_median": float(np.nanmedian(speed)),
        "speed_fraction_above_2_5": float(np.mean(speed > 2.5)),
        "behavior_blocks": blocks,
        "two_maze_gap_candidate": two_maze_candidate,
        "top_timestamp_gaps": top_gaps,
    }
    stem = args.nwb.stem
    (args.output_dir / f"{stem}_behavior_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    relative_minutes = (timestamps - timestamps[0]) / 60
    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex=True, constrained_layout=True)
    axes[0].plot(relative_minutes, position, color="#355C7D", linewidth=0.55)
    axes[0].set_ylabel("Linearized position")
    axes[0].set_title(args.nwb.name)
    axes[1].plot(relative_minutes, speed, color="#C06C84", linewidth=0.45)
    axes[1].axhline(2.5, color="black", linestyle="--", linewidth=0.8, label="2.5 cm/s")
    axes[1].set_xlabel("Minutes from first behavioral sample")
    axes[1].set_ylabel("Speed")
    axes[1].legend(frameon=False, loc="upper right")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        for split_index in split_after:
            split_minute = float((timestamps[split_index] - timestamps[0]) / 60)
            ax.axvline(split_minute, color="#6C5B7B", linestyle=":", linewidth=1.1)
    fig.savefig(args.output_dir / f"{stem}_behavior_overview.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
