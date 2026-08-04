#!/usr/bin/env python3
"""Refresh the majority-mouse support flag in existing cross-area summaries."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def main() -> None:
    summary_paths = sorted(
        path
        for path in Path("outputs").glob("**/summary.csv")
        if any(part.startswith("cross_area") for part in path.parts)
    )
    for csv_path in summary_paths:
        table = pd.read_csv(csv_path)
        table["supported_positive_gain"] = (
            table["mouse_level_mean_delta_mse"].gt(0)
            & table["bh_q_positive_gain_across_8_tests"].le(0.05)
            & table["mice_positive"].ge(3)
        )
        table.to_csv(csv_path, index=False)
        json_path = csv_path.with_suffix(".json")
        if json_path.exists():
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            payload["results"] = table.to_dict(orient="records")
            json_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print(csv_path)


if __name__ == "__main__":
    main()
