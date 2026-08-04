# Dataset audit findings

Audit date: 2026-08-03
Source: DANDI:001695, version `0.260319.2023`

## What is actually available

The published Dandiset contains 22 NWB assets (3.088 GB) from six identifiers.
Fifteen assets contain freely moving behavior and simultaneous CA3, CA1, and
RSC recordings from four mice (M01, M02, M03, M05). The remaining assets are
non-equivalent head-fixed recordings and are outside the spatial endpoint.

The fifteen navigation sessions divide cleanly into:

- 12 single-maze novel-learning sessions from four mice;
- 3 familiar two-maze context-switch sessions from M01, M02, and M03.

The two-maze classification was verified from behavioral timestamps rather
than filenames. Each of those sessions has two long behavioral blocks separated
by a conspicuous gap (22.2–128.2 s). The original protocol states that both
mazes had been experienced previously, although room enclosure and context
changed. They therefore cannot support a claim about formation of a novel map.

## Audited support

All 15 navigation files contain all three target regions and at least 20
behaviorally detected traversals.

| Quantity | Minimum | Median | Maximum |
|---|---:|---:|---:|
| CA1 pyramidal units | 44 | 85 | 173 |
| CA3 pyramidal units | 11 | 37 | 141 |
| RSC non-narrow units | 7 | 55 | 116 |
| Detected traversals | 71 | 89 | 116 |

The RSC definition follows the released demonstration notebook: cells not
classified as narrow interneurons. Region/cell-type definitions will be
checked against the final paper and sensitivity-tested before inference.

## Consequences for the project

1. **Primary question:** emergence of a stable late-session spatial map during
   twelve novel-maze sessions.
2. **Secondary control:** rapid switching between two familiar maps in three
   mice.
3. **Primary biological replicate:** mouse. The twelve sessions are repeated
   observations from four mice, not twelve independent animals.
4. **Cell-count imbalance:** regional comparisons require within-session
   equal-unit resampling or a reliability-normalized statistic.
5. **No effect has been inspected:** the completed work covers metadata,
   behavior, unit counts, traversal detection, and late-map reliability only.

## Effect-blind late-map gate

Every navigation session supplied ten late traversals per direction and at
least 87.5% occupied position bins. Spearman-Brown corrected odd/even
late-map reliability was positive in every region/direction/block. Across 36
direction/block observations per region, median reliability was 0.59 in CA1,
0.48 in CA3, and 0.41 in RSC. These values are quality measurements, not tests
of regional onset.

All 12 novel sessions remain in the primary cohort. A prespecified sensitivity
requires a session-level minimum reliability of 0.15; this gate was fixed
without computing any early-traversal score.

## Frozen decision gate

The map-maturity estimator now:

- treats running directions separately;
- uses only movement samples;
- learns a late-session template without evaluating it on the same traversals;
- equalizes cell-count support across regions;
- returns one value per traversal and ultimately one contrast per mouse;
- is calibrated on label/time-shift nulls and same-session split reliability.

The consolidated ledger is in
`data/metadata/navigation_session_ledger.csv`; the versioned audit figure is
`outputs/figures/dataset_audit.png`.
