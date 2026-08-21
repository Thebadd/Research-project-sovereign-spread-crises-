# Methodology

> Why this paper uses several estimators, what question each one answers,
> and exactly where the design departs from Asonuma, Chamon, Erce &
> Sasahara (2024). `EMPIRICAL_ANALYSIS.md` reports the results;
> `RESULTS_SECTION_DRAFT.md` is the paper-ready prose. This file is the
> specification reference behind both.

---

## 1. The estimator ladder

The paper runs five estimators. They are not competing refinements of one
another, and the later ones are not "better versions" of the earlier ones.
Each exists because the previous stage leaves a specific question
unanswered. Read in order, they are three questions:

**Question 1 — What is the average dynamic response to a spread crisis?**

Jordà (2005) local projections, estimated horizon by horizon
(`02_lp_all.do` for all onsets, `03_lp_resolution.do` for the non-default /
default-linked split). This produces Figure 1, Figure 2 and Table 2 — the
descriptive core of the paper. Everything downstream is a response to
something this specification cannot settle on its own.

**Question 2 — Is the comparison group legitimate, or is this selection?**

The local projection compares crisis country-years to tranquil
country-years conditional on controls. But countries that go on to default
were already different: more indebted, weaker external positions, further
into a banking crisis. If those same differences also drive future output,
the LP coefficient blends the effect of the crisis with the effect of
having been the kind of country that has one.

Two estimators address this. Plain inverse-probability weighting
(`08_ipw_lp.do`) models the probability of onset and reweights the control
group to resemble the treated on observables. The doubly-robust AIPW
estimator (`08b_aipw.do`, following Jordà & Taylor 2016 and Asonuma et
al.'s Eq. 3) adds the outcome model back on top of the weighting, so the
estimate is consistent if *either* the propensity model or the outcome
model is correctly specified rather than requiring both. AIPW is the
preferred selection-corrected specification; plain IPW is retained as the
intermediate step and a robustness row.

**Question 3 — Through what, and for whom?**

Having established a cost and shown it is not an artifact of selection, the
paper opens it up. The same projection is re-estimated on six intermediate
outcomes (`11_channels.do`, split by resolution type in
`12_channels_resolution.do`); a Gelbach (2016) decomposition
(`12b_gelbach_decomposition.do`) asks how much of the default-linked GDP
coefficient each channel statistically absorbs; and interaction and
median-split designs (`13b`, `13c`, `13d`) ask whether the cost depends on
pre-crisis bank exposure to the sovereign.

So the sequence is: **establish it → check it is not selection → open it
up.** Three questions, not one question estimated five ways.

**Question 4 — What does it cost to be *in* one?**

The three questions above all treat a crisis as an event. `20_lp_flow.do` asks
the state question instead: with treatment set to 1 in every year of an
episode, what is output doing while a country is in a spread crisis? This is a
different estimand, not a refinement of the first three, and it sits on a
different rung. Questions 1–3 identify, subject to the selection correction;
Question 4 describes, because no propensity model for a time-varying treatment
is well-posed here (§5.2). Where the two disagree, the onset tier carries the
claim and the disagreement is itself reported.

### AIPW here is for selection, not nonlinearity

This distinction matters for how the estimator is justified in writing, and
it is easy to get wrong because the applied local-projection literature
often motivates AIPW the other way.

Where the treatment is a *continuous* shock, researchers commonly define
treatment by quantile — a large improvement is the top quartile, a large
deterioration the bottom, with the middle 50% as the control group — and
use AIPW to recover nonlinear, state-dependent effects. The benchmark in
that design is the middle of the shock distribution, not a zero-shock
counterfactual.

**None of that applies here.** The treatment in this paper is binary: a
spread crisis either begins in a given country-year or it does not. There
is no distribution to split, and the control group is tranquil years — a
genuine zero-treatment counterfactual. The reason for AIPW is entirely the
selection problem set out under Question 2. A justification framed around
nonlinearity or state dependence would not describe what this estimator is
doing in this paper.

### Why Eq. (3) is hand-coded rather than `teffects aipw`

Stata's canned `teffects aipw` was considered and rejected for three
reasons: it does not admit country fixed effects in the outcome model; it
cannot produce a paired bootstrap of the *difference* between two treatment
types, which is the quantity this paper actually tests; and it would not
reproduce the reference paper's estimator, foreclosing a like-for-like
comparison. The `_aipw` program in `08b_aipw.do` implements Eq. (3)
directly, with the IPW-weighted outcome regression supplying the
conditional means.

---

## 2. Fixed effects, by stage

The split is by **estimator stage, not by Act**. This is the single most
misread aspect of the design, so the table is authoritative:

| Stage | Files | Fixed effects |
|---|---|---|
| Single-stage local projection | `02`, `03`, `06`, `07`, `12b`, `13b`, and the OLS halves of `11`, `11b`, `12` | country **and year** |
| Two-stage IPW / AIPW | `08`, `08b`, `13c`, `13d`, `21`, and the IPW halves of `11`, `11b`, `12` | country **only** |
| Flow tier, single-stage | `20` | country **and year** — the rule applies |

**The flow tier splits across the rule, and that has a consequence worth
stating.** `20_lp_flow.do` is single-stage and keeps year fixed effects;
`21_aipw_flow.do` is two-stage and drops them, matching the reference paper,
whose `$convar` carries `c1-c74` and no year dummies. So Figure 9b and
Figure 10 — the analogues of their Fig. 3 (OLS) and Fig. 4 (AIPW) — differ in
the **estimator and in the year fixed effects at the same time**, and the
difference between them cannot be attributed to the estimator alone.

`20` therefore reports a **`r_noyearfe`** robustness row, dropping year FE for
both the pooled coefficient and the def−nd difference. That row is the
specification `21` runs, so it is what isolates the estimator, and the write-up
should quote it whenever the two figures are set side by side.

Expect the row to differ materially from the baseline. The year dummies are
economically real — 2020 at −8.95, 2019 at −9.92, 2009 at −4.01 — and
default-linked episodes cluster in 2019–2022, so without them COVID and the GFC
can load on the treatment. The onset analogue is on record below: the
default-linked coefficient moves from −5.53/−6.75 at Years 1–2 with year FE to
−7.58/−9.48 without.

One consequence to carry forward: §1 of `EMPIRICAL_ANALYSIS.md` justifies
leaving VIX and UST10Y out of the outcome equation on the grounds that year
fixed effects absorb every shock common to all countries. That justification
holds in `20` but **not** in `21`. The global-push terms remain excluded there,
by the same omission the reference paper uses rather than by absorption.

Year fixed effects in the single-stage projections are a deliberate
improvement on the reference paper, whose `$convar` carries country dummies
(`c1-c74`) and no year dummies at all. They absorb every shock common to
all countries in a given year — the global financial cycle, in particular —
directly, rather than relying on a proxy.

The two-stage estimators drop them to match the reference design exactly,
so the doubly-robust results are comparable like-for-like with theirs.

**What the year fixed effects are worth, in coefficients.** This is not a
cosmetic difference and the two designs do not produce the same number, so
a reader comparing tables across the paper needs the reason. On identical
data and controls, the default-linked coefficient is **−5.53 at Year 1 and
−6.75 at Year 2 with year FE** (`03`, `xtscc`) and **−7.58 and −9.48
without them** (`08`, `areg`, country FE only). The pooled Act-1 cost moves
similarly, from −2.98 to −4.42 at Year 2. Year fixed effects are therefore
worth roughly two percentage points of estimated output cost.

The interpretation is straightforward: year dummies absorb the global
financial cycle, so a country-FE-only specification leaves the common
component of a crisis year inside the treatment coefficient. Since spread
crises cluster in global risk-off episodes, that inflates the estimate. The
paper's headline is the smaller, year-FE number — the conservative one —
while the two-stage estimates are reported on the reference paper's design
so the doubly-robust comparison is like-for-like. Neither is wrong; they
answer slightly different questions, and the difference between them is
itself informative about how much of the measured cost is global rather
than country-specific.

Note in particular that `03_lp_resolution.do` — the non-default /
default-linked comparison behind Table 2 and Figure 2 — is **single-stage
and does carry year fixed effects**. "The resolution split has no year FE"
is not an accurate summary; the resolution split estimated *by the
two-stage estimators* has none.

---

## 3. Inference, by stage

Three layers, each protecting against something different.

**Driscoll-Kraay** (`xtscc`, lag length `max(1, h+1)`) in every
single-stage local projection. Clustering by country permits arbitrary
correlation within a country over time but requires independence *across*
countries within a year — close to untenable in a panel of emerging-market
sovereign spreads, where a common global risk factor means a bad year for
one borrower is systematically a bad year for the others. Driscoll-Kraay
corrects for that cross-sectional dependence as well as for the serial
correlation that grows mechanically with the horizon, which is why the lag
length is tied to `h`. The full argument, with its qualifications, is in
`EMPIRICAL_ANALYSIS.md` §1.

This is a departure from the reference paper, which uses `vce(robust)` for
its descriptive projections and `cluster(wdicode)` for its doubly-robust
estimates. It is also a departure from the common recommendation to treat
clustered errors as the baseline and Driscoll-Kraay as a strict robustness
check. Leading with the conservative estimator is defensible here precisely
because the headline result **survives** it: when the estimator that grants
the fewest assumptions still rejects, reporting it as the baseline is a
strength rather than a risk. Where results do weaken under Driscoll-Kraay —
the Year 4–5 horizons — that is reported as weakening, not suppressed.

**Cluster-robust** (`vce(cluster cid)`) in the two-stage estimators. This
is mechanical rather than elective: `xtscc` does not accept probability
weights, so an inverse-probability-weighted outcome regression cannot be
estimated with it (see the note at `08_ipw_lp.do:186`). Those stages fall
back to `areg ... absorb(cid) vce(cluster cid)`, which has the incidental
benefit of matching the reference paper's inference for the same stage.

**Stratified cluster bootstrap** (`bsample`, in `08b`, `13c`, `13d`) for
the AIPW confidence intervals, and — the important case — for the
default-minus-non-default *difference*, which is the quantity the paper's
central claim rests on. Bootstrapping the difference on paired draws gives
it an interval of its own, rather than inferring it from the visual gap
between two separately estimated bands. With roughly twenty default-linked
episodes, resampling is more defensible than asymptotic standard errors.

One deviation worth recording: the reference paper's `bsample` carries no
`cluster()` option and resamples individual rows within treatment-type
strata. Ours resamples whole countries (`bsample, cluster(cid) strata(...)
idcluster(_bid)`), which is the stricter choice given that observations
within a country are plainly not independent.

**Flow specification (`20_lp_flow.do`).** Driscoll-Kraay again, but at
`lag(max(2, h+3))` rather than `max(1, h+1)`. The onset rule covers the $h+1$
overlap of successive outcome windows. Flow coding introduces a second source
of dependence that onset coding does not have — the *regressor* is serially
correlated within an episode — so the lag carries an additional 2 for
persistence at the median episode duration of two years. The onset lag rule is
reported as a second row, so that a flow-versus-onset comparison is not
confounded by a simultaneous change of inference, and country-clustered
standard errors (`areg … vce(cluster cid)`, 52 clusters) as a third.

The more important discipline here is not the lag but the counting. Flow coding
produces 234 treated rows, but they still carry only **61 episodes in 52
countries**, and a handful of chronic cases supply a large share of them. Every
flow table reports episodes and countries alongside N, because a table showing
only N would imply far more independent information than the design contains.

---

## 4. Controls: the deviation map

The reference paper's control object is:

```stata
global convar gdpg2 gov_exp2 open2 banking_duration2 credit_bank2 ///
              hyperinf_dummy ex_dum1-ex_dum5 c1-c74
```

plus `g'v'_0`, the own-outcome pre-trend, added by hand in every regression
and *not* inside `$convar`. Ours is:

```stata
global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp ///
                  l_open l_credit_bank l_hyperinfl"
```

Six of the eight terms map one-to-one:

| Ours | Theirs | Content |
|---|---|---|
| `l1_gdpg` | `gdpg2` | lagged real GDP growth |
| `l_govexp` | `gov_exp2` | government expenditure / GDP, t−1 |
| `l_open` | `open2` | trade openness, t−1 |
| `l_banking_duration` | `banking_duration2` | years a Laeven–Valencia banking crisis has run, t−1 |
| `l_credit_bank` | `credit_bank2` | bank credit to the private sector / GDP, t−1 |
| `l_hyperinfl` | `hyperinf_dummy` | hyperinflation flag, t−1 |

Three differences, each a choice with a reason:

**We add `l_debt` and `l_ca`; they have neither.** Deliberate. Their
treatment is a restructuring; ours is a market-priced spread crossing a
threshold, and public debt and the external position are precisely what
that spread is priced off. Omitting them would leave the most direct
determinants of the treatment out of the conditioning set.

**We lack their `ex_dum1`–`ex_dum5`; they have five.** Their bins are
quintiles of `ln(1+L.exchange) − ln(1+L2.exchange)`, cut at the p5/p25/p50/
p75/p95 of the restructuring sample. This is a nominal-depreciation
control, and depreciation is bundled with default in ways that matter. This
is the one **genuine gap** in the control set, and it is a data limitation
rather than a design choice: the panel carries no nominal exchange-rate
series. It should be stated as such rather than passed over.

**The own-outcome pre-trend plays the same role in both.** Their `g'v'_0`
is rebuilt per outcome as `L.var'v' - L2.var'v'`; our per-channel
`pre_<var>` is the same object. For the GDP equation theirs coincides with
`gdpg2` and is therefore entered twice; ours enters once via `l1_gdpg`.

One further detail: their `hyperinf_dummy` is built as
`replace hyperinf_dummy = 1 if L.inflation > 50` with no missing-value
guard, so it reads 0 wherever inflation is missing. Ours is conditioned on
`!missing(L.infl)`. The two are close but not identical, and ours is the
more careful construction.

### Excluded predictors

Both papers exclude a set of predictors from the outcome equation so they
can identify the propensity model. Theirs is
`federal_funds2 cont_all2 past_preemp2`; ours is the lagged fed funds rate,
a regional contagion measure, and a count of the country's own past onsets
(`past_onsets`, or `past_def_onsets` for the resolution-type stage).

The exclusion mechanism differs by stage in our design and should be
described accordingly. In the **two-stage** files the predictors are
excluded by simple omission from the outcome equation — the same mechanism
the reference paper uses. In the **single-stage** projections there is no
first stage at all, so nothing is being excluded there; the global cycle is
absorbed by year fixed effects, which is why a pure time-series regressor
would be redundant in those specifications rather than excluded from them.

---

## 5. Treatment definition

The project estimates two treatments, because the question "what follows the
start of a spread crisis" and the question "what is the output cost of being
in one" are not the same question and cannot be answered by the same variable.

### 5.1 Onset coding — the identifying design

In `02`, `03`, `06`–`08b` and `11`–`13f`, each episode contributes **exactly
one treated observation**, its onset year; continuation years are dropped from
the estimation sample entirely (`18_transforms.do`, `sample`). The multi-year
dynamics of an episode are traced through the horizon of the outcome, not
through additional treated rows.

Horizons are labelled so the crisis year itself is **Year 1**, matching the
reference paper's convention. Year 0 is a hard-coded pre-crisis baseline,
always zero and never estimated. A single genuinely estimable pre-trend
placebo is reported at Year −1. Internal variable names (`dy_0`…`dy_4`) are
unchanged — only the displayed axis is shifted.

This convention is what makes the propensity-score tier coherent: `Pr(onset)`
is a hazard model on predetermined covariates, so the two-stage estimators in
`08`/`08b`/`13c`/`13d` have a well-posed first stage. Onset coding therefore
carries the paper's identification claim.

### 5.2 Flow coding — the state, estimated in `20_lp_flow.do`

`in_crisis` is 1 for **every year of an episode**, onset and continuation
alike: 234 treated country-years rather than 61, with tranquil years as
controls. Treatment is *episode membership*, not the annual criterion flag
`crisis_any`, which is 0 on 13 rows the episode-dating rule places inside an
episode and would push mid-episode years into the control pool.

The horizon convention is **the same as §5.1** — Year 0 the baseline, Year 1 the
crisis year. `dy_h` is differenced against the row's own $t-1$, so $h=-1$ gives
$y_{t-1}-y_{t-1}=0$ by construction for continuation rows exactly as for onset
rows; the Year-0 zero is a genuine normalisation under either treatment. What
differs is what the baseline year *is*: "the year before this crisis year"
rather than "the year before the crisis started".

At Year 1 the outcome is growth **during** that crisis year, which is the object
the research question names. Beyond it, the estimate is the cumulative change
from $t-1$ associated with being in a crisis year, averaged over how long the
country has *already* been in crisis — it is **not** "years after onset",
because the row generating it is itself treated, and elapsed duration is an
outcome of the crisis. Sharing the axis lets Figure 9 be laid directly over
Figure 1, but the two Year-1 coefficients are **not the same object**: the onset
one is the first year of every episode, the flow one pools the first year of one
episode with the fourth year of another. No pre-trend placebo is written at
Year −1 (below).

**Controls are dated at the episode's entry year.** Under flow coding, every
element of `$ctrl_core` measured at a continuation row's own $t-1$ is an
outcome of the crisis that row is already in — lagged growth most obviously,
but debt, the current account, bank credit, government spending and
banking-crisis duration equally. Conditioning on them is conditioning on
mediators, on the right-hand side. `$ctrl_flow` (built in `18_transforms.do`)
therefore evaluates the common core at the year the episode began and holds it
fixed across the episode; tranquil rows keep their own $t-1$, which is
predetermined for them trivially.

This is not a departure from the reference paper — it *reproduces* their
timing. Their treated rows are always onsets, so their `L.`/`L2.` controls are
predetermined by construction and the question never arises for them; row-dated
controls are the thing that would depart. Onset rows are unaffected here, since
for them $t-1$ *is* the entry year, so the onset and flow specifications share
an identical control set at Year 1 and their Year-1 coefficients are directly
comparable. The construction follows Callaway & Sant'Anna (2021), who measure
covariates in each cohort's pre-treatment period for the same reason. The
row-dated set is reported as a robustness column so the choice is visible.

**The flow tier stops at the single-stage projection, deliberately.** For a
time-varying treatment the strongest predictor of being in crisis at *t* is
being in crisis at *t−1* — a post-treatment outcome of the same episode — so
unconfoundedness given X is false whether or not that term enters the first
stage. Include it and the propensity model becomes a persistence model whose
scores approach 1 on continuation rows, which the `[0.01, 0.99]` trim in
`08_ipw_lp.do` then deletes: precisely the rows flow coding exists to add. The
correct estimator would be a marginal structural model with
inverse-probability-of-treatment weights over the treatment *history* (Robins,
Hernán & Brumback 2000), which is out of scope. Flow results are accordingly
reported as **magnitude and state conditional on selection into crisis**, not
as a second identified effect.

No pre-trend placebo is estimable under flow coding: `dy_m2` differences
$t-2$ against $t-1$, and for a row three or more years into an episode both
endpoints are crisis years, so regressing it on `in_crisis` would reject by
construction. The valid restriction is the onset sample, whose placebo `02`
and `03` already report.

### 5.3 What each design can and cannot say about duration

Onset coding identifies the coefficient equally from a five-year crisis and a
one-year one, so on its own it cannot speak to whether protracted crises are
more costly. This matters more than it first appears, because the two
resolution groups differ systematically in duration — non-default episodes
average 2.8 years and 113 crisis-years across 40 episodes, default-linked
episodes 5.8 years and 121 crisis-years across 21. At Years 3–5 most
non-default episodes have ended while most default-linked ones are still
running, so part of the measured resolution gap reflects that difference.
Flow coding weights episodes by their length instead, which is the appropriate
weighting for a question about the *state* but concentrates the treatment in a
handful of chronic cases; `20_lp_flow.do` reports the treated share of each
country's panel years and a leave-Venezuela-out variant for that reason.
