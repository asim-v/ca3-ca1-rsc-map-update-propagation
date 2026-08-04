#!/usr/bin/env python3
"""Consolidate prespecified map-maturity sensitivity runs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REGIONS = ["CA3", "CA1", "RSC"]
KEYS = ["subject", "file", "block", "direction"]


def summarize(name: str, path: Path) -> dict:
    table = pd.read_csv(path)
    metadata = json.loads((path.parent / "map_maturity_summary.json").read_text(encoding="utf-8"))
    table["early_registered"] = table["early_identity_null_p"].le(0.05) & table[
        "early_spatial_null_p"
    ].le(0.05)
    table["early_margin"] = table["early_five_mean"] - table[
        ["early_identity_null_q95", "early_spatial_null_q95"]
    ].max(axis=1)
    status = table.pivot(index=KEYS, columns="region", values="fit_status")
    mouse_margin = table.groupby(["subject", "region"])["early_margin"].median()
    row = {
        "configuration": name,
        "n_subsamples": int(metadata["parameters"]["n_subsamples"]),
        "n_nulls": int(metadata["parameters"]["n_nulls"]),
        "all_three_transition_curves": int(status.eq("detectable_transition").all(axis=1).sum()),
        "paired_ca3_ca1_transition_curves": int(
            (status["CA3"].eq("detectable_transition") & status["CA1"].eq("detectable_transition")).sum()
        ),
    }
    for region in REGIONS:
        subset = table[table["region"] == region]
        row[f"early_registered_{region}"] = int(subset["early_registered"].sum())
        row[f"detectable_transition_{region}"] = int(
            subset["fit_status"].eq("detectable_transition").sum()
        )
        row[f"mice_positive_early_margin_{region}"] = int(
            mouse_margin.xs(region, level="region").gt(0).sum()
        )
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/sensitivity/summary.csv"))
    args = parser.parse_args()
    root = Path("outputs/sensitivity")
    sources = {
        "primary": Path("outputs/maturity/map_maturity_fits.csv"),
        "16 position bins": root / "bins16/map_maturity_fits.csv",
        "32 position bins": root / "bins32/map_maturity_fits.csv",
        "speed > 5 cm/s": root / "speed5/map_maturity_fits.csv",
        "8 late-reference traversals": root / "late8/map_maturity_fits.csv",
        "12 late-reference traversals": root / "late12/map_maturity_fits.csv",
    }
    rows = [summarize(name, path) for name, path in sources.items()]
    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    print(output.to_string(index=False))


if __name__ == "__main__":
    main()
