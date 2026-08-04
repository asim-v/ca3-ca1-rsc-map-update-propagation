#!/usr/bin/env python3
"""Build the non-effect navigation-session eligibility ledger."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pynwb import NWBHDF5IO

from extract_traversals import behavior_blocks, endpoint_events, linearize_position


def count_traversals(position: np.ndarray, timestamps: np.ndarray, speed: np.ndarray) -> tuple[int, list[int]]:
    counts: list[int] = []
    for block_start, block_stop in behavior_blocks(timestamps):
        pos = position[block_start : block_stop + 1]
        lo, hi = np.nanpercentile(pos, [1, 99])
        normalized = np.clip((pos - lo) / (hi - lo), 0, 1)
        events = endpoint_events(normalized, 0.15, 0.85)
        accepted = 0
        for (start_local, side_a), (stop_local, side_b) in zip(events[:-1], events[1:]):
            if side_a == side_b:
                continue
            start, stop = block_start + start_local, block_start + stop_local
            duration = timestamps[stop] - timestamps[start]
            coverage = np.nanmax(normalized[start_local : stop_local + 1]) - np.nanmin(normalized[start_local : stop_local + 1])
            if 0 < duration <= 120 and coverage >= 0.70 and np.mean(speed[start : stop + 1] > 2.5) > 0:
                accepted += 1
        counts.append(accepted)
    return sum(counts), counts


def inspect(path: Path) -> dict[str, object]:
    with NWBHDF5IO(str(path), "r", load_namespaces=True) as io:
        nwb = io.read()
        units = nwb.units.to_dataframe()
        behavior = nwb.processing["behavior"]
        position_series = behavior["AnimalPosition"]["Position"]
        position = linearize_position(np.asarray(position_series.data[:]))
        timestamps = np.asarray(position_series.timestamps[:]).squeeze()
        speed = np.asarray(behavior["Speed"].data[:]).squeeze()

    areas = units["cell_area"].astype(str)
    types = units["cell_type"].astype(str)
    blocks = behavior_blocks(timestamps)
    total_traversals, block_counts = count_traversals(position, timestamps, speed)
    gaps = np.diff(timestamps)
    subject_match = re.search(r"sub-(M\d+)", path.name)
    date_match = re.search(r"ses-(\d{8})", path.name)
    return {
        "subject": subject_match.group(1) if subject_match else "",
        "date": date_match.group(1) if date_match else "",
        "file": path.name,
        "session_description": nwb.session_description,
        "n_units_total": len(units),
        "n_ca1_pyramidal": int(((areas == "CA1") & (types == "Pyramidal Cell")).sum()),
        "n_ca3_pyramidal": int(((areas == "CA3") & (types == "Pyramidal Cell")).sum()),
        "n_rsc_non_narrow": int(((areas == "RSC") & (types != "Narrow Interneuron")).sum()),
        "behavior_samples": len(timestamps),
        "behavior_duration_minutes": round(float((timestamps[-1] - timestamps[0]) / 60), 4),
        "maximum_timestamp_gap_seconds": round(float(gaps.max()), 4),
        "n_behavior_blocks": len(blocks),
        "is_two_maze": len(blocks) == 2 and all((timestamps[stop] - timestamps[start]) / 60 >= 10 for start, stop in blocks),
        "n_traversals": total_traversals,
        "traversals_by_block": ";".join(map(str, block_counts)),
        "speed_fraction_above_2_5": round(float(np.mean(speed > 2.5)), 4),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--output", type=Path, default=Path("data/metadata/navigation_session_ledger.csv"))
    args = parser.parse_args()

    files = sorted(args.raw_dir.glob("*behavior+ecephys.nwb"))
    if not files:
        raise SystemExit(f"No behavior NWB files found in {args.raw_dir}")
    rows = [inspect(path) for path in files]
    table = pd.DataFrame(rows).sort_values(["subject", "date"])
    table["eligible_region_counts"] = (
        (table.n_ca1_pyramidal > 0)
        & (table.n_ca3_pyramidal > 0)
        & (table.n_rsc_non_narrow > 0)
    )
    table["eligible_behavior"] = table.n_traversals >= 20
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)

    summary = {
        "sessions": len(table),
        "subjects": sorted(table.subject.unique().tolist()),
        "single_maze_sessions": int((~table.is_two_maze).sum()),
        "two_maze_sessions": int(table.is_two_maze.sum()),
        "sessions_with_three_regions": int(table.eligible_region_counts.sum()),
        "sessions_with_at_least_20_traversals": int(table.eligible_behavior.sum()),
        "units": {
            region: {"min": int(table[column].min()), "median": float(table[column].median()), "max": int(table[column].max())}
            for region, column in {
                "CA1 pyramidal": "n_ca1_pyramidal",
                "CA3 pyramidal": "n_ca3_pyramidal",
                "RSC non-narrow": "n_rsc_non_narrow",
            }.items()
        },
        "traversals": {"min": int(table.n_traversals.min()), "median": float(table.n_traversals.median()), "max": int(table.n_traversals.max())},
    }
    args.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
