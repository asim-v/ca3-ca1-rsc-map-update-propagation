#!/usr/bin/env python3
"""Consolidate strict, reference-excluded cross-area sensitivity results."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def main() -> None:
    root = Path("outputs/sensitivity")
    sources = {
        "primary": Path("outputs/cross_area/exclude_reference/summary.csv"),
        "16 position bins": root / "bins16/cross_area_exclude_reference/summary.csv",
        "32 position bins": root / "bins32/cross_area_exclude_reference/summary.csv",
        "speed > 5 cm/s": root / "speed5/cross_area_exclude_reference/summary.csv",
        "8 late-reference traversals": root / "late8/cross_area_exclude_reference/summary.csv",
        "12 late-reference traversals": root / "late12/cross_area_exclude_reference/summary.csv",
    }
    rows = []
    for configuration, path in sources.items():
        table = pd.read_csv(path)
        for _, row in table.iterrows():
            rows.append(
                {
                    "configuration": configuration,
                    "connection": f"{row['source']}->{row['target']}",
                    "lag_same_direction_traversals": int(row["lag_same_direction_traversals"]),
                    "mouse_level_mean_delta_mse": float(row["mouse_level_mean_delta_mse"]),
                    "mice_positive": int(row["mice_positive"]),
                    "bh_q_positive_gain": float(row["bh_q_positive_gain_across_8_tests"]),
                    "supported_positive_gain": bool(row["supported_positive_gain"]),
                }
            )
    output = pd.DataFrame(rows)
    destination = Path("outputs/sensitivity/cross_area_summary.csv")
    destination.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(destination, index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
