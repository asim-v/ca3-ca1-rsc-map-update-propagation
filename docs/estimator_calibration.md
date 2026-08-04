# Map-maturity estimator: calibration and provisional lock

## Purpose

Regional ordering will not be inspected until the traversal-scale endpoint is
specified. `scripts/calibrate_map_maturity.py` tests the estimator on synthetic
place-cell populations whose transition midpoint is known. The simulation does
not read spikes from DANDI:001695.

## Locked primary endpoint

For each session, direction, and region:

1. Retain movement samples above 2.5 cm/s and map every traversal to 24
   normalized position bins.
2. Form square-root-transformed occupancy-normalized rate maps.
3. Use each traversal as the primary temporal unit.
4. Construct the late reference from the final ten traversals of that
   direction, splitting odd and even traversals so reliability can be measured.
5. Center each cell's spatial profile and correlate the flattened window map
   with the late reference.
6. Equalize the number of cells across CA3, CA1, and RSC within a session by
   repeated subsampling. The regional statistic is the median across 200 fixed
   subsamples.
7. Fit one bounded sigmoid per session and region. The descriptive timing
   endpoint is its midpoint; the confirmatory contrast is the within-session
   difference in midpoint between regions, summarized within mouse.

Five- and nine-traversal windows remain required sensitivities. Neither window
width is interpreted as a biological timescale.

## Synthetic calibration result

Across 100 simulations at 12, 24, 48, and 96 cells, cell count did not create a
systematic shift in the recovered midpoint. Single-traversal scores recovered
the imposed midpoint about 0.6 traversals early, compared with about 1.5 for
five-traversal windows and 1.9 for nine-traversal windows. Wider windows blurred
the transition without improving midpoint recovery in this calibration.

The single-traversal score is therefore primary. The remaining offset is an
expected property of the nonlinear map-to-correlation relationship; absolute
midpoint values are descriptive, while paired regional differences are the
target estimand.

## Effect-blind quality gate

All 15 navigation sessions supplied at least ten late traversals per direction,
at least 87.5% occupied spatial bins, and positive late-map split-half
reliability in every region/direction/block. All 12 novel-maze sessions
therefore enter the primary analysis. No session will be excluded based on its
CA3-versus-CA1-versus-RSC ordering.

The primary analysis carries measurement uncertainty without a hard
reliability exclusion. A prespecified sensitivity requires every regional,
directional late-map reliability in a session to be at least 0.15. This cutoff
was set before early-traversal or regional-onset scores were computed.

Generated artifacts:

- `outputs/calibration/estimator_calibration.png`
- `outputs/calibration/estimator_recovery.csv`
- `outputs/calibration/calibration_summary.json`
- `outputs/quality/map_reliability.csv`
- `outputs/quality/map_reliability.summary.json`
