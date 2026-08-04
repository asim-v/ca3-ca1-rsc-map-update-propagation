# CA3–CA1–RSC map-update propagation

This repository tests where a new spatial-map state first becomes detectable
and how much of the CA1 change is predictable from upstream CA3 activity and
downstream RSC activity.

## Scientific question

> When a mouse enters a different maze, does CA3 express the new map before
> CA1, and does CA1 preserve an upstream component while adding a component
> that better predicts RSC?

The project is deliberately narrower than a “whole hippocampus model.” Its
first target is the simultaneously recorded CA3–CA1–RSC navigation subset of
[DANDI:001695](https://dandiarchive.org/dandiset/001695/0.260319.2023).
DG and CA2 recordings in the release are not assumed to be comparable to the
freely navigating three-region sessions.

## Claim boundary

The observational recordings can establish temporal ordering and
cross-validated prediction, not synaptic causality. Sessions are repeated
within animals, so animal—not neuron, time bin, lap, or session—is the primary
unit of population inference.

## Planned endpoint

For each area, learn a position- and direction-conditioned late-session spatial
template from held-out stable traversals. Apply the frozen template to early
traversals and estimate:

1. the first traversal at which the late-session map geometry is reliably
   expressed during novel-maze learning;
2. whether CA3 activity predicts later CA1 map-change score beyond CA1 history
   and measured behavior;
3. whether the CA1 residual not predicted by CA3 predicts RSC;
4. whether a separate three-animal, familiar two-maze subset shows a comparable
   ordering during rapid context switches;
5. whether these relations survive position, direction, speed, occupancy,
   circular-shift, and within-position label nulls.

The primary development set contains 12 novel-maze sessions from four mice.
The context-switch subset contains three two-maze sessions from three mice;
both mazes had been experienced previously, so this subset cannot establish
learning of a new map.

## Repository layout

- `scripts/audit_dandiset.py`: metadata-only DANDI audit.
- `scripts/download_asset.py`: checksum-verified download of a selected NWB.
- `scripts/inspect_nwb.py`: NWB inventory and unit counts by area/cell type.
- `scripts/calibrate_map_maturity.py`: synthetic, effect-blind estimator check.
- `scripts/audit_map_reliability.py`: late-map occupancy and reliability gate.
- `scripts/score_map_maturity.py`: frozen traversal-resolved regional score.
- `scripts/control_behavior_maturity.py`: nonlinear trial-level speed and
  movement-duration control.
- `scripts/cross_area_prediction.py`: traversal-scale forward/reverse prediction
  with mouse-level circular-shift inference.
- `data/metadata/`: small, versioned audit tables and source metadata.
- `data/raw/`: ignored public NWB files.
- `outputs/`: derived audit and analysis artifacts.
- `docs/`: decisions, milestones, and analysis specifications.

## Reproduce the metadata audit

```powershell
python scripts\audit_dandiset.py
```

To inspect one pilot file after installing the analysis dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
python scripts\download_asset.py --asset-id <DANDI_ASSET_ID>
.\.venv\Scripts\python.exe scripts\inspect_nwb.py data\raw\<FILE>.nwb
```

## Data provenance

Gonzalez, J. & Voroslakos, M. (2026). *High-density extracellular
recordings (Neuropixels/SiNAPS) of Hippocampal-Cortical dynamics during
spatial behavior*. DANDI Archive, version `0.260319.2023`.
[https://doi.org/10.48324/dandi.001695/0.260319.2023](https://doi.org/10.48324/dandi.001695/0.260319.2023)

No raw data are redistributed by this repository.

## Current status

The complete navigation subset has been downloaded and audited locally. No
CA3-versus-CA1-versus-RSC neural ordering has been inspected. The
single-traversal map-maturity estimator, its synthetic calibration, and the
effect-blind occupancy/reliability gate are documented in
`docs/estimator_calibration.md`. All sessions passed the primary gate. The
frozen-cohort result and sensitivity analysis are in
`docs/preliminary_findings.md`: early CA1 map identity is common, while a
universal CA3-to-CA1-to-RSC onset sequence is not supported. Cross-area lag,
reverse-direction, reference-exclusion, and familiar-map controls likewise do
not support a robust traversal-scale propagation claim.
