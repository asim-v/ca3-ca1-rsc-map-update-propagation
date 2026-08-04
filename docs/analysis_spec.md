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

The confirmatory cohort was frozen after a metadata-and-quality audit that did
not inspect early-traversal scores or the proposed cross-area ordering. All 12
novel-maze sessions passed: each direction supplied ten late traversals, at
least 87.5% of bins were occupied, and late-map reliability was positive.

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

For traversal-scale coupling, fit fixed-alpha ridge models separately by
session and direction. Compare leave-one-traversal-out error for a reduced
model containing target history, traversal number and its square, and measured
behavior with a full model that adds either current or previous source-region
score. Standardize the target so error reduction is comparable across curves.
Use every unique circular shift of the source except shifts 0 and +/-1 as the
curve-level null. Aggregate curve effects by median within mouse and then mean
across mice; obtain the group null by sampling one circular null per curve.
Test CA3-to-CA1 and CA1-to-RSC as planned directions and CA1-to-CA3 and
RSC-to-CA1 as reverse-direction controls. Correct the eight exploratory tests
with Benjamini-Hochberg. Label a positive gain as supported only when the
mouse-level mean is positive, corrected `q <= 0.05`, and at least three of four
mouse medians are positive. Circular-shift p-values are conditional tests in
the recorded animals, not substitutes for population replication.

As an explicitly exploratory localization, repeat the coupling analysis on the
first 20 and last 20 traversals of each direction. These equal-length phase
windows are not additional confirmatory tests and may overlap in the shortest
sessions.

Because late-reference traversals participate in template construction even
under leave-one-out scoring, the cross-area coupling claim must survive a
strict analysis that removes all final ten reference traversals. Results that
exist only inside the reference block are treated as estimator dependence, not
evidence of biological propagation.

Before interpreting a regional transition, remove the component of its
traversal-level similarity score predicted out of sample by mean running speed,
duration, moving-sample fraction, their squared terms, and a speed-by-duration
interaction. Use fixed-alpha (`alpha = 1`) leave-one-traversal-out ridge
prediction. Report both raw and behavior-controlled classifications; measured
behavior is a control, not a latent behavioral-state estimate.

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
