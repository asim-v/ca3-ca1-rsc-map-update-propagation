# Preliminary findings after frozen-cohort analysis

Analysis date: 2026-08-03

## Result in one sentence

The data do not reveal a reproducible traversal-scale CA3-to-CA1-to-RSC map
construction sequence; instead, late-session spatial identity is often already
registered during the first five traversals, especially in CA1, while the
strength of that registered map can continue to mature heterogeneously.

## Primary novel-maze cohort

The analysis contains 12 sessions from four mice and two directions per
session, giving 24 repeated session-direction curves per region. These curves
are descriptive repeated measurements; the biological replicate remains the
mouse.

Early activity matched the cell-position identity of the late map beyond both
a cell-identity permutation null and an independent within-cell spatial-shift
null, after Benjamini-Hochberg correction across all 72 regional curves, in:

- CA1: 21/24 curves;
- CA3: 13/24 curves;
- RSC: 15/24 curves.

The median early-minus-null margin was positive for CA1 in all four mice, for
CA3 in two of four mice, and for RSC in three of four. With only four animals,
this supports a descriptive pattern rather than a population-level significance
claim.

Detectable monotonic maturation occurred in 8/24 CA1, 4/24 CA3, and 2/24 RSC
curves. All eight CA1 maturation curves already showed registered late-map
identity in the first five traversals. Thus early identity and later maturation
are not mutually exclusive: CA1 can select a stable spatial arrangement early
and subsequently strengthen or stabilize it.

## No supported serial propagation claim

No raw session-direction curve showed detectable maturation in all three
regions. Only three curves allowed a CA3-versus-CA1 midpoint comparison, and
their CA3-minus-CA1 differences were +12.39, -0.73, and -7.39 traversals. The
signs are mixed and the curves come from only two mice.

After nonlinear leave-one-traversal-out adjustment for measured speed,
duration, moving-sample fraction, squared terms, and a speed-by-duration
interaction, only one curve showed transitions in all three regions. Five
CA3-versus-CA1 comparisons were available, all from M03, and again had mixed
signs. Measured behavior predicted map similarity poorly out of sample
(median cross-validated R-squared = -0.11).

These results do not support a universal CA3-first propagation sequence. They
also do not exclude fast propagation occurring within a traversal or at a
timescale finer than this analysis.

## Traversal-scale cross-area prediction

An autoregressive ridge model tested whether a source-region map-similarity
score improved leave-one-traversal-out prediction of another region beyond
target history, traversal number, and nonlinear measured-behavior terms.
Circular shifts preserved each source trajectory while destroying its alignment
to the target.

Using every traversal initially produced positive same-traversal prediction in
both directions for CA3-CA1 and CA1-RSC, but no positive gain from the previous
same-direction traversal. The forward directions were not stronger than their
reverse controls. Temporal localization then showed that CA1-RSC coupling was
absent in the first 20 traversals and strongest near the end of the session.

Because the last ten traversals define the late template, the strict analysis
removed them completely. CA3-CA1 prediction disappeared. Same-traversal
CA1-to-RSC prediction remained in the primary parameterization, but only one to
three of four mouse medians were positive across five sensitivity settings.
No connection had a positive, corrected, majority-mouse gain in every
configuration. The familiar-map control likewise had no corrected predictive
gain.

Therefore cross-area coupling is retained as a heterogeneous exploratory
observation, not evidence for propagation. Conditional circular-shift p-values
can be small within the recorded animals, but with four biological replicates
they do not establish population generalization.

## Familiar-map control

The three two-maze sessions contain previously experienced maps and are a
context-switch control, not novel learning. Early identity exceeded both nulls
in 9/12 CA1, 10/12 CA3, and 8/12 RSC block-direction curves. CA1 nevertheless
showed detectable maturation in 7/12 raw curves. Therefore a later increase in
map similarity is not specific evidence for learning a novel map; it can also
occur while reinstating a familiar representation after a context switch.

## Sensitivity

Across 16 versus 24 versus 32 position bins, a 5 cm/s movement threshold, and
8 versus 10 versus 12 late-reference traversals:

- CA1 early identity occurred in 18-21 of 24 curves and had a positive
  mouse-level median in 4/4 animals in every configuration;
- CA3 early identity occurred in 11-15 of 24 curves;
- RSC early identity occurred in 14-17 of 24 curves;
- simultaneous three-region maturation occurred in zero curves in five of six
  configurations and one curve with 16 bins;
- only three or four curves per configuration allowed CA3-versus-CA1 midpoint
  comparison.

The CA1 early-identity pattern is the most stable observation. Exact maturation
counts are more parameter-sensitive and should remain secondary.

## Claim boundary and next test

The strongest justified statement is that late-map cell-position identity is
detectable very early in CA1 and that subsequent map-strength changes are
heterogeneous and not uniquely associated with novelty. This is not evidence
of innateness, synaptic causality, or an instantaneous whole-circuit map.

The next decisive step should not add more models to these same four animals.
It should seek an independent simultaneous CA3-CA1-RSC dataset or a new
experiment with finer within-traversal timing. The present dataset supports an
early-CA1-identity result, but does not have enough independent replication to
establish a circuit-wide propagation mechanism.
