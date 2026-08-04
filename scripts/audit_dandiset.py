#!/usr/bin/env python3
"""Metadata-only audit of a published DANDI version.

Uses only Python's standard library so the audit can run before installing the
NWB analysis environment. It does not download neural data.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

API = "https://api.dandiarchive.org/api"


def get_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"User-Agent": "ca3-ca1-rsc-audit/0.1"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def session_record(asset: dict[str, Any]) -> dict[str, Any]:
    session = next(
        (x for x in asset.get("wasGeneratedBy", []) if x.get("schemaKey") == "Session"),
        {},
    )
    participant = next(iter(asset.get("wasAttributedTo", [])), {})
    description = session.get("description", "")
    lower = description.lower()
    path = asset.get("path", "")
    content_urls = asset.get("contentUrl", [])
    digest = asset.get("digest", {})
    return {
        "asset_id": asset.get("identifier", ""),
        "subject": participant.get("identifier", ""),
        "path": path,
        "size_bytes": asset.get("contentSize", 0),
        "size_mb": round(asset.get("contentSize", 0) / 1_000_000, 3),
        "session_start": session.get("startDate", ""),
        "session_end": session.get("endDate", ""),
        "session_description": description,
        "has_behavior": "behavior+ecephys" in path,
        "is_maze": "maze" in lower,
        "is_novel_maze": "novel maze" in lower,
        "is_two_maze": bool(re.search(r"(?:2|two)[ -]?maze", lower)),
        "has_pre_sleep": "pre" in lower and "sleep" in lower,
        "has_post_sleep": "post" in lower and "sleep" in lower,
        "sha256": digest.get("dandi:sha2-256", ""),
        "download_url": content_urls[0] if content_urls else "",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dandiset", default="001695")
    parser.add_argument("--version", default="0.260319.2023")
    parser.add_argument("--output-dir", type=Path, default=Path("data/metadata"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    version_url = f"{API}/dandisets/{args.dandiset}/versions/{args.version}/"
    metadata = get_json(version_url)

    page_url = f"{version_url}assets/?page_size=100"
    listed: list[dict[str, Any]] = []
    while page_url:
        page = get_json(page_url)
        listed.extend(page.get("results", []))
        page_url = page.get("next")

    asset_ids = [x.get("asset_id") or x.get("identifier") for x in listed]
    detail_urls = [f"{API}/assets/{asset_id}/" for asset_id in asset_ids]
    with ThreadPoolExecutor(max_workers=8) as pool:
        details = list(pool.map(get_json, detail_urls))

    rows = sorted((session_record(x) for x in details), key=lambda x: x["path"])
    fields = list(rows[0]) if rows else []
    with (args.output_dir / "assets.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    (args.output_dir / "dandiset_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    descriptions: dict[str, int] = {}
    for row in rows:
        descriptions[row["session_description"]] = descriptions.get(row["session_description"], 0) + 1
    summary = {
        "dandiset": args.dandiset,
        "version": args.version,
        "doi": metadata.get("doi"),
        "asset_count": len(rows),
        "total_size_bytes": sum(int(x["size_bytes"]) for x in rows),
        "subjects": sorted({x["subject"] for x in rows if x["subject"]}),
        "behavior_assets": sum(bool(x["has_behavior"]) for x in rows),
        "maze_assets": sum(bool(x["is_maze"]) for x in rows),
        "novel_maze_assets": sum(bool(x["is_novel_maze"]) for x in rows),
        "two_maze_assets": sum(bool(x["is_two_maze"]) for x in rows),
        "session_descriptions": descriptions,
    }
    (args.output_dir / "audit_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    lines = [
        "# DANDI metadata audit",
        "",
        f"- Dandiset: `{args.dandiset}`",
        f"- Version: `{args.version}`",
        f"- Assets: {summary['asset_count']}",
        f"- Total size: {summary['total_size_bytes'] / 1_000_000_000:.3f} GB",
        f"- Subjects: {', '.join(summary['subjects'])}",
        f"- Assets with behavior: {summary['behavior_assets']}",
        f"- Maze assets: {summary['maze_assets']}",
        f"- Novel-maze assets: {summary['novel_maze_assets']}",
        f"- Two-maze assets: {summary['two_maze_assets']}",
        "",
        "## Session descriptions",
        "",
    ]
    lines.extend(f"- {count} × {description or '(missing)'}" for description, count in sorted(descriptions.items()))
    (args.output_dir / "AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
