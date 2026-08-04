# Analysis specification — development version

## Biological estimand

The primary object is traversal-resolved similarity to a held-out late-session
spatial template in each simultaneously recorded area. It is not a claim that
the same neurons or axes exist across regions.

## Stage 1: dataset gate

A session can enter development only if it has:

- behavioral position and timestamps spanning identifiable maze epochs;
- CA3 and CA1 pyramidal units and RSC principal/non-narrow units;
- enough movement samples and repeated traversals in each compared maze;
- no unresolved mismatch between neural and behavioral time bases.

The confirmatory cohort and minimum unit/traversal thresholds will be frozen
only after a metadata-and-quality audit that does not inspect the proposed
cross-area effect.

## Stage 2: regional map-change score

Within each session and region:

1. bin spiking on a common time grid;
2. restrict to movement and define traversal, direction, position, and maze;
3. fit a cross-validated classifier/regression score for maze identity using
   late, stable traversals while matching position, direction, speed, and
   occupancy;
4. freeze the model and project all early traversals;
5. aggregate time bins to one score per traversal before regional comparisons.

The principal onset statistic is the earliest traversal after maze entry whose
animal-level score reliably approaches the late-session spatial template. The
exact similarity measure and sequential criterion must be frozen before effect
inspection. The three familiar two-maze sessions are a context-switch control,
not evidence about novel-map learning.

## Stage 3: propagation test

Test whether CA3 score at traversal `k` predicts CA1 score at `k` or `k+1`
beyond:

- past CA1 score;
- position/direction coverage;
- speed and occupancy summaries;
- elapsed time and traversal number.

Then test whether the component of CA1 unexplained by CA3 predicts RSC. These
are predictive temporal relations, not causal mediation estimates.

## Required nulls

- maze labels shuffled within position × direction strata;
- circular time shift of one region relative to the others;
- traversal-order shuffle within maze;
- cell-label and equal-unit-count resampling;
- fixed-condition or same-maze recurrence where available;
- alternative lags and bin widths;
- leave-one-animal-out sensitivity.

## Inference

- Primary unit: animal.
- Sessions and traversals: repeated measurements, not independent replicates.
- Report every animal and session alongside hierarchical summaries.
- No neuron-count-inflated p-values.
