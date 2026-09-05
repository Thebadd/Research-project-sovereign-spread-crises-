# Project writing standard

All prose in this project — `EMPIRICAL_ANALYSIS.md`, `RESULTS_SECTION_DRAFT.md`,
`METHODOLOGY.md`, `PAPER_FRAMING.md`, do-file header comments, table notes,
figure notes, commit messages — is written as a researcher writing an academic
economics paper, consistent with the empirical literature this project is
modelled on (Asonuma, Chamon, Erce & Sasahara 2024; Jordà 2005; Jordà & Taylor
2016).

## What that requires

**Present the methodology, then justify it.** Every specification choice —
estimator, control set, fixed effects, inference, sample restriction — is
stated and then defended on economic or statistical grounds. A choice that
cannot be defended is a choice to reconsider, not to leave unexplained. Where
the project departs from the reference paper, say so and give the reason (see
`METHODOLOGY.md`).

**Present the results, then interpret them.** Report the estimate first and
the reading second, so the reader can disagree with the interpretation while
accepting the number.

**Explain why each robustness check exists.** A robustness table is only
informative if the reader knows what would have counted as failure. State the
concern the check addresses, then whether the result survives it — including
when it does not.

## The binding constraint: the narrative follows the evidence

**The narrative must emerge from the actual empirical results, not from a
story decided in advance.** Do not force results to fit a preconceived
hypothesis or mechanism. Let the evidence determine the argument, while
keeping the exposition coherent and economically meaningful.

In practice this means:

- **No sign priors written into a design.** `13b_exposure_heterogeneity.do`
  originally asserted that a negative interaction "means the channel
  transmits." That was never estimated; it was an assumption about what the
  exposure variable measures, and most exposures came back positive. Define
  in advance what each coefficient *would* mean under competing readings, and
  let the data pick.
- **A null is a result.** Report it as such. The regional-contagion predictor
  never worked in any first-stage column; the pooled `claims_govt` channel
  never moves. Both are reported, not quietly dropped.
- **A result that dies under a better specification is reported as having
  died.** The `claims_govt`/GDP interaction was significant on a legacy
  control set and is not on the common core. The write-up says so and explains
  why the assets-scaled measure is the estimable one.
- **State what the estimate cannot bear.** Sample sizes, multiple-testing
  arithmetic, and cells too thin to support a claim belong next to the claim,
  not only in a limitations section.
- **Distinguish co-movement, statistical absorption, and identified
  amplification.** These are three different strengths of evidence and the
  language must not blur them. "Is associated with" is not "causes"; a
  Gelbach share is not mediation.
- **Do not present the same object estimated twice as two findings.** Where a
  question is answered by more than one design, organise them as a ladder —
  effect, robustness, mechanism — and say which is doing the identifying work
  (see `EMPIRICAL_ANALYSIS.md` §8 for the pattern).

## Tone

Plain declarative academic prose. Continuous argument rather than bullet
lists in the paper-facing documents. No overclaiming, no hedging past the
point of informativeness, and no rhetorical flourishes standing in for
evidence. Numbers carry the argument; adjectives do not.
