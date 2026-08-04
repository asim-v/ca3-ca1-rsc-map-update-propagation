#!/usr/bin/env python3
"""Create a compact, non-effect NWB inventory for one pilot session."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
from pynwb import NWBHDF5IO


def shape_of(value: Any) -> list[int] | None:
    shape = getattr(value, "shape", None)
    return list(shape) if shape is not None else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("nwb", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/pilot_inventory"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with NWBHDF5IO(str(args.nwb), mode="r", load_namespaces=True) as io:
        nwb = io.read()
        units = nwb.units.to_dataframe() if nwb.units is not None else pd.DataFrame()
        group_columns = [x for x in ("cell_area", "cell_type") if x in units.columns]
        if group_columns:
            counts = units.groupby(group_columns, dropna=False).size().rename("n_units").reset_index()
        else:
            counts = pd.DataFrame({"n_units": [len(units)]})
        counts.to_csv(args.output_dir / "units_by_region_cell_type.csv", index=False)

        processing: dict[str, Any] = {}
        for module_name, module in nwb.processing.items():
            processing[module_name] = {}
            for name, interface in module.data_interfaces.items():
                item = {"type": type(interface).__name__}
                if hasattr(interface, "data"):
                    item["data_shape"] = shape_of(interface.data)
                if hasattr(interface, "time_series"):
                    item["time_series"] = {
                        ts_name: {
                            "type": type(ts).__name__,
                            "data_shape": shape_of(getattr(ts, "data", None)),
                            "timestamps_shape": shape_of(getattr(ts, "timestamps", None)),
                            "rate": getattr(ts, "rate", None),
                        }
                        for ts_name, ts in interface.time_series.items()
                    }
                if hasattr(interface, "spatial_series"):
                    item["spatial_series"] = {
                        ts_name: {
                            "type": type(ts).__name__,
                            "data_shape": shape_of(getattr(ts, "data", None)),
                            "timestamps_shape": shape_of(getattr(ts, "timestamps", None)),
                            "rate": getattr(ts, "rate", None),
                        }
                        for ts_name, ts in interface.spatial_series.items()
                    }
                processing[module_name][name] = item

        intervals = {}
        for name, table in nwb.intervals.items():
            intervals[name] = {"rows": len(table), "columns": list(table.colnames)}

        inventory = {
            "path": str(args.nwb.resolve()),
            "identifier": nwb.identifier,
            "session_id": nwb.session_id,
            "session_description": nwb.session_description,
            "session_start_time": nwb.session_start_time.isoformat(),
            "n_units": len(units),
            "unit_columns": list(units.columns),
            "acquisition": {name: type(obj).__name__ for name, obj in nwb.acquisition.items()},
            "processing": processing,
            "intervals": intervals,
        }
        (args.output_dir / "nwb_inventory.json").write_text(
            json.dumps(inventory, indent=2, default=str), encoding="utf-8"
        )
        print(json.dumps({"session": inventory["session_description"], "n_units": len(units), "counts": counts.to_dict("records")}, indent=2))


if __name__ == "__main__":
    main()
