#!/usr/bin/env python3
"""Download one audited DANDI asset and verify size and SHA-256."""

from __future__ import annotations

import argparse
import csv
import hashlib
import urllib.request
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--asset-id", required=True)
    parser.add_argument("--ledger", type=Path, default=Path("data/metadata/assets.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("data/raw"))
    args = parser.parse_args()

    with args.ledger.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    matches = [row for row in rows if row["asset_id"] == args.asset_id]
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one ledger row for {args.asset_id}; found {len(matches)}")
    row = matches[0]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    destination = args.output_dir / Path(row["path"]).name
    partial = destination.with_suffix(destination.suffix + ".partial")
    request = urllib.request.Request(row["download_url"], headers={"User-Agent": "ca3-ca1-rsc-download/0.1"})
    digest = hashlib.sha256()
    downloaded = 0
    with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as f:
        while chunk := response.read(1024 * 1024):
            f.write(chunk)
            digest.update(chunk)
            downloaded += len(chunk)

    expected_size = int(row["size_bytes"])
    expected_sha = row["sha256"].lower()
    observed_sha = digest.hexdigest()
    if downloaded != expected_size or observed_sha != expected_sha:
        partial.unlink(missing_ok=True)
        raise SystemExit(
            f"Verification failed: size {downloaded}/{expected_size}, "
            f"sha256 {observed_sha}/{expected_sha}"
        )
    partial.replace(destination)
    print(destination.resolve())


if __name__ == "__main__":
    main()
