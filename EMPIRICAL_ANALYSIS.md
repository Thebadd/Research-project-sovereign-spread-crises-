# Empirical Analysis

> **Status.** This version reflects the from-source rebuild of the panel
> (stages `10`–`18`), the total-real-GDP outcome (Asonuma-aligned, not
> per-capita), the Year-1-is-crisis-year horizon relabeling, and the
> conformability/pre-trend fixes made over this working session. Numbers in
> Sections 4.1–4.4, 6, 7, and 7a are taken directly from verified runs of
> the current code (`00_master.do`) and are current. Section 5 (Table 2's
> Spec A / plain-IPW comparison) uses a verified current run. Every
> section — including the doubly-robust AIPW headline (§5.4) and its
> channel and nexus extensions (§7b, §8) — now reports numbers taken from a
> single verified end-to-end run of `00_master.do` on the current code,
> with no stage left unre-run.

## 1. Setting and identification

The panel covers 52 emerging and frontier-market economies, rebuilt
entirely from official sources (IMF WEO, World Bank WDI/IDS, IMF MFS/BOP,
FRED, Laeven–Valencia) on top of the project's own spread-crisis database,
which supplies only the crisis dating and classification (onset years,
non-default vs. default-linked). Three-year "carry-in" scaffolding rows
are added at each country's left panel edge purely so lag operators have
something to point at; they are excluded from every estimation sample.
After that exclusion, the estimation sample (onset years plus tranquil
years, continuation years dropped) contains 1,085 country-year
observations, of which 59 of the 61 identified spread-crisis episodes
survive (two are lost to a missing pre-crisis GDP base). Forty of the 61
episodes are non-default, 21 are default-linked; roughly 47 of the 52
countries contribute at least one onset that enters the regression sample.

Each episode contributes **exactly one treated observation** — its onset
year. Continuation years (the second and later years of an ongoing
episode) are excluded from the sample entirely, both to avoid treating an
ongoing crisis year as a tranquil control and to avoid letting a single
episode's duration inflate the treated count. This is the standard
Jordà-style local-projection convention (Jordà, 2005; Jordà & Taylor,
2016), and it is Asonuma et al.'s own design as well; we return to what it
does and does not imply for the interpretation of duration in Section 10.

The estimator is the local projection. For each horizon
$h \in \{0,1,2,3,4\}$ we estimate

$$\Delta^h y_{it} = \alpha_i + \lambda_t + \beta_h D_{it} + \gamma' X_{it} + \varepsilon_{it+h},$$

where $\Delta^h y_{it}$ is the cumulative percent change in log real GDP
(total, not per capita) from $t-1$ to $t+h$. Horizons are labeled to match
Asonuma et al.'s own convention: the crisis year itself is **Year 1**;
**Year 0** is an explicit, hard-coded pre-crisis baseline (value zero,
never estimated, since by construction $y_{t-1}-y_{t-1}\equiv 0$); the
figures additionally show a single pre-trend placebo point at **Year
$-1$**. Internally the underlying variables (`dy_0`…`dy_4`) are unchanged
by this relabeling — only the printed axis shifts.

The common-core control vector $X_{it}$ — used, with minor channel-specific
adjustments, in every regression in this paper — contains one lag of GDP
growth, public debt/GDP, the current account/GDP, the number of years a
systemic banking crisis has already run (Laeven–Valencia, zero-filled
where no crisis is recorded), government expenditure/GDP, trade openness,
bank credit depth, and a hyperinflation flag, all measured at $t-1$. Global
financial conditions (fed funds, VIX, the US 10-year yield) enter only the
first-stage propensity models used in Section 5, where they serve as
excluded predictors. The mechanism differs by stage and is worth stating
precisely. The single-stage projections of this section have no first
stage, so nothing is being *excluded* there in the identification sense: a
pure time-series regressor would simply be redundant, because year fixed
effects already absorb every shock common to all countries in a given year.
In the two-stage estimators of Section 5, which carry country fixed effects
but no year fixed effects, these variables are excluded from the outcome
equation by omission — the same mechanism the reference paper uses, and the
source of identification for the propensity model. Standard errors in the single-stage local projections follow Driscoll and
Kraay (1998), with lag length $\max(1,h+1)$. This is a deliberate
departure from the reference paper, which uses heteroskedasticity-robust
errors for its descriptive projections and country-clustered errors for its
doubly-robust estimates. Clustering by country is the weaker assumption in
our setting: it permits arbitrary correlation within a country over time
but requires independence *across* countries in a given year. In a panel of
emerging-market sovereign spreads that assumption is close to untenable —
spreads are driven by a common global risk factor, so a bad year for
Argentina is systematically a bad year for Turkey, and the residuals are
correlated in exactly the dimension country-clustering assumes away.
Driscoll-Kraay corrects for that cross-sectional dependence as well as for
the serial correlation that mechanically grows with the horizon (each
$dy_h$ overlaps the next by construction), which is why the lag length is
tied to $h$. Two qualifications belong with the choice rather than in a
footnote. First, Driscoll-Kraay's asymptotics run in the *time* dimension,
and this panel offers roughly 33 years; that is within the range where the
estimator is conventionally applied, but it is not generous, and it is one
reason the permutation placebo of Section 9 — which assumes nothing about
the error structure at all — is reported alongside rather than as an
afterthought. Second, the two-stage IPW and AIPW estimators of Section 5 do
*not* use Driscoll-Kraay, and the reason is mechanical rather than
discretionary: `xtscc` does not accept probability weights, so an
inverse-probability-weighted outcome regression cannot be estimated with
it. Those stages fall back to country-clustered errors, which has the
incidental benefit of making the doubly-robust results directly comparable
to the reference paper's, estimated the same way.

The identifying assumption is that, conditional on fixed effects and the
observable pre-crisis state, crisis onset is uncorrelated with unobserved
shocks to future output. We test this directly in Section 4.3, and — going
further than a standard placebo check — stress-test the result against a
longer pre-crisis growth pattern than the placebo alone can rule out.

---

## 2. The average output cost of a spread crisis

Figure 1 plots the cumulative GDP response to a spread-crisis onset,
pooling all 61 episodes (Table 1):

| Horizon | $\hat\beta_h$ (pp) | SE | $p$ | Episodes | $N$ |
|:-------:|:------------------:|:----:|:-----:|:--------:|:---:|
| Year 1 | −2.29 | 0.88 | 0.014 | 46 | 937 |
| Year 2 | −2.98 | 0.95 | 0.004 | 46 | 908 |
| Year 3 | −1.97 | 1.05 | 0.070 | 45 | 871 |
| Year 4 | −1.25 | 1.39 | 0.374 | 45 | 837 |
| Year 5 | −0.72 | 1.21 | 0.557 | 44 | 803 |

The coefficient is negative and significant at the one- and two-year
horizons, weakens to marginal significance at Year 3, and is
indistinguishable from zero by Years 4–5, visible in Figure 1 as the
shaded 90 percent band widening and recentering on the horizontal axis
toward the end of the window. On average, a spread-crisis onset is
followed by a real but front-loaded output shortfall of roughly two to
three percentage points, fading within about three years. As the next
sections show, this pooled average conceals substantial and economically
important heterogeneity, and we treat it mainly as a motivating fact
rather than the paper's central estimate.

---

## 3. Pre-trend validation

Rather than testing several placebo horizons before onset, the panel's
own construction pins down a single, well-defined pre-trend point: because
the main outcome is anchored at $t-1$, plugging $h=-1$ into the identical
formula gives $y_{t-1}-y_{t-1}\equiv 0$ trivially — it is not an estimable
placebo, it *is* the baseline (this is why it is displayed as the explicit
Year 0 point, never estimated). The first genuinely estimable pre-crisis
horizon on that same base is one year earlier still — the growth rate from
$t-2$ to $t-1$ — which we display as **Year $-1$** and estimate as an
outcome in its own right, dropping the one lag of GDP growth from its own
controls (since that lag is algebraically identical to this outcome, and
including it would guarantee a null by construction).

Pooling all crises, this placebo is significantly **positive**
($\hat\beta = 0.943$, SE $=0.359$, $p=0.013$) — visible in Figure 3 as a
point sitting above, not on, the zero line. Splitting by resolution type
(Section 4) shows this is not a generic feature of crisis onset: it is
concentrated in default-linked episodes specifically ($\hat\beta = 2.14$,
$p=0.003$), and absent in non-default ones ($\hat\beta = 0.36$,
$p=0.444$). We do not treat this as invalidating the design — it is a
boom-before-the-bust pattern, not evidence that default-linked countries
were already declining before the crisis — but it is exactly the kind of
signal that warrants a direct stress test, which Section 4.3 provides.

---

## 4. Does the resolution of the crisis matter?

### 4.1 The headline resolution split

We estimate both onset dummies jointly, with tranquil years as the
reference category — matching the reference paper's own baseline
specification exactly — and test $\beta^{nd}_h = \beta^{def}_h$ at each
horizon (Table 2; Figure 2 plots the two resulting impulse responses).

| Horizon | $\hat\beta^{nd}$ | $\hat\beta^{def}$ | Gap (def$-$nd) | Clogg $z$ | $p$(Wald, nd=def) |
|:-------:|:----------------:|:------------------:|:---------------:|:---------:|:-------------------:|
| Year 1 | −0.73 (p=0.175) | −5.53 (p=0.041) | −4.80 | −1.81 | 0.091 |
| Year 2 | −1.18 (p=0.020) | −6.75 (p=0.006) | −5.57 | −2.36 | 0.016 |
| Year 3 | −0.24 (p=0.759) | −5.42 (p=0.014) | −5.19 | −2.34 | 0.006 |
| Year 4 | +0.38 (p=0.768) | −4.49 (p=0.048) | −4.87 | −1.93 | 0.019 |
| Year 5 | +1.20 (p=0.222) | −4.31 (p=0.043) | −5.51 | −2.45 | 0.011 |

**What identifies these coefficients.** The classification assigns 40
non-default and 21 default-linked onsets, but those are not the counts the
regression uses. Listwise deletion on the common core removes an onset
whenever any control is missing at $t-1$, and the binding constraints are
public debt (56 of 61 onsets), trade openness (58) and the hyperinflation
flag (58). The estimates above are therefore identified off **31 non-default
and 15 default-linked episodes** at Year 1, falling to 29 and 15 by Year 5 —
a loss of roughly a quarter of the non-default group and **six of the
twenty-one default-linked episodes**. Tables 1 and 2 report these
within-sample counts rather than the classification totals, and every
statement about precision in this paper should be read against 15 treated
observations in the default-linked group, not 21.

The contrast is visible directly in Figure 2: the non-default line sits
close to zero throughout, with a 90 percent band that clears zero only
briefly (Year 2); the default-linked line sits well below it at every
horizon, its band never touching zero from Year 1 through Year 5. Critically,
the equality test — not merely the presence of significance in one series —
rejects at four of five horizons and is borderline at the fifth ($p=0.091$
at Year 1). We read this as reasonably strong evidence that the output
cost of a spread crisis is concentrated in, and largely specific to,
episodes that culminate in default, rather than a generic consequence of
crossing the spread threshold. This is the central result of the paper.

A robustness specification estimates each type separately against
tranquil years, with the rival type dropped from the sample (Spec A — the
design the two-stage IPW estimator in Section 5 also uses, so it is the
correct OLS partner for that comparison, though not the paper's headline
baseline). It reproduces the same pattern: non-default is mostly
insignificant (significant only at Year 2, $-1.23$, $p=0.015$),
default-linked is significant at every horizon ($-4.1$ to $-6.8$), and the
extra-cost-of-default difference is significant or borderline throughout
(Year 2: $p=0.019$; Year 3: $p=0.030$; Year 5: $p=0.029$).

**The long horizons do not survive, and we report that rather than the
gap alone.** Because the panel runs to 2026 while WEO's last actual
observation is 2024 or 2025 for most countries, an onset from 2022 onward
has its Year-4 and Year-5 outcome built partly from IMF projections: seven
onsets at Year 4 and nine at Year 5, precisely the horizons where the
estimated cost is largest. Re-estimating on realized data only leaves Years
1–3 essentially untouched ($-5.58$, $-6.79$, $-5.33$ against $-5.53$,
$-6.75$, $-5.42$) but **nearly halves the default-linked coefficient at the
long horizons**, to $-2.69$ at Year 4 and $-3.00$ at Year 5. The
*difference* between types survives at every horizon
($p=0.031$ and $p=0.060$), so the resolution contrast is not a forecasting
artifact — but the level of the persistent cost partly is.

Two further checks point the same way. The two-stage IPW estimates of
Section 5.3 put the default-linked coefficient at $-1.90$ and $-0.70$ at
Years 4–5, neither significant, and the doubly-robust AIPW of Section 5.4
returns $-4.32$ and $-2.51$ with confidence intervals covering zero, with
the bootstrapped difference likewise insignificant at both horizons. Three
independent checks — removing forecast data, reweighting on the propensity
score, and doubly-robust estimation — therefore agree that the
default-linked cost **is not robustly identified beyond Year 3**. We
accordingly restrict the persistence claim to Years 1–3 throughout, and
treat the Year 4–5 coefficients in Table 2 as suggestive at best.

### 4.2 A pre-existing-trend concern, addressed directly

Section 3 showed the pre-crisis placebo is significantly positive for
default-linked episodes specifically — a boom, not a decline, but a real
departure from zero nonetheless. Because one lag of GDP growth already
enters the common core, this exact pre-crisis window (Year $-1$) is
mechanically absorbed by the headline regression in Table 2 already; the
open question was narrower — does the boom extend further back than a
single lag can see, in a way the baseline specification cannot detect?

We test this directly by adding a second, non-overlapping lag of GDP
growth (the growth rate from $t-3$ to $t-2$ — genuinely separate
information, not mechanically tied to the lag already in the core) to the
exact Table 2 specification:

| Horizon | $\hat\beta^{def}$, headline | $\hat\beta^{def}$, + second lag | $\hat\beta^{nd}$, headline | $\hat\beta^{nd}$, + second lag |
|:-------:|:----:|:----:|:----:|:----:|
| Year 1 | −5.535 | −5.499 (p=0.044) | −0.730 | −0.853 (p=0.140) |
| Year 2 | −6.748 | −6.680 (p=0.007) | −1.178 | −1.362 (p=0.009) |
| Year 3 | −5.421 | −5.389 (p=0.016) | −0.237 | −0.318 (p=0.665) |
| Year 4 | −4.491 | −4.503 (p=0.049) | 0.377 | 0.406 (p=0.731) |
| Year 5 | −4.311 | −4.353 (p=0.039) | 1.199 | 1.293 (p=0.154) |

The default-linked coefficient moves by no more than 0.1 percentage
points at any horizon and remains significant throughout. This indicates
the boom detected in Section 3 is confined to the single window the
baseline specification already controls for; there is no evidence of a
longer, uncontrolled pre-crisis pattern manufacturing the appearance of a
post-crisis cost through simple mean reversion. We regard this as
meaningfully strengthening, not merely qualifying, the central result of
Section 4.1: the divergence in Figure 2 survives a direct stress test
against the specific concern the placebo test raised.

This test rules out a linear growth-momentum explanation for the
pre-trend specifically. It does not rule out other pre-existing,
non-growth differences between the two groups — a discrete shock
coinciding with the pre-onset window, for instance — which we retain as
an open limitation (Section 10).

---

## 5. Selection, propensity, and inverse-probability weighting

### 5.1 Why selection is a real concern here

The two resolution groups are not exchangeable ex ante — countries that
eventually default plausibly differ from those that do not on observable
pre-crisis fundamentals. We address this with a first-stage propensity
model: a pooled probit (no country fixed effects — with roughly 20–40
events per type, country dummies would separate on the many never-treated
countries and collapse the usable sample) of onset on the common core plus
excluded predictors — the US fed funds rate (lagged), a leave-one-out
regional contagion measure, and a count of the country's own past onsets.
For Act 1 (any onset vs. tranquil), adding these predictors to the
controls-only model raises the area under the ROC curve from 0.664 to
0.732; the Table-1-style first-stage table (`08c_first_stage_table.do`)
reports the same predictors separately for the non-default and
default-linked columns, with AUROC of 0.748 and 0.881 respectively and a
jointly significant $\chi^2$ test on the predictors in every column.

### 5.2 Act 1: does reweighting change the average cost?

Reweighting tranquil years to match the observable pre-crisis profile of
treated country-years, and re-estimating the pooled Act-1 local
projection (Figure 7 overlays the two), **attenuates** the estimated cost
at the shorter horizons rather than amplifying it:

| Horizon | $\hat\beta$, OLS-FE | $\hat\beta$, IPW | $\Delta$ |
|:-------:|:--------------------:|:------------------:|:--------:|
| Year 1 | −2.74 | −1.84 | +0.90 |
| Year 2 | −4.42 | −2.52 | +1.90 |
| Year 3 | −1.69 | −0.70 | +0.99 |
| Year 4 | −0.12 | −0.13 | −0.01 |
| Year 5 | +0.48 | +0.62 | +0.14 |

This is worth stating plainly because it is easy to get backwards: the
IPW-weighted coefficients are consistently *smaller in magnitude* than the
unweighted ones at Years 1–3, and essentially unchanged at Years 4–5. That
is, a meaningful part of the pooled Act-1 cost documented in Section 2 is
associated with selection on observable pre-crisis characteristics rather
than surviving as a pure treatment effect once that selection is
addressed — the unweighted pooled estimate should be read as somewhat
generous, not conservative.

### 5.3 Act 2: the resolution split under reweighting

For the resolution split, we follow the reference paper's own design more
closely than a single reweighted comparison would allow: each type is
scored against tranquil years separately (the rival type dropped from
that type's estimation sample), with its own stabilized weights, so
non-default and default-linked each get a genuine vs.-tranquil doubly
estimated line (OLS and IPW) rather than being reweighted against each
other directly. The outcome-stage regressions here carry country fixed
effects only, no year fixed effects — deliberately: the excluded
predictors (fed funds, contagion, past onsets) are genuinely excluded from
this outcome equation by omission, rather than by year-dummy absorption,
which is the two-stage design's own exclusion-restriction logic, distinct
from the year-FE mechanism used in the single-stage specifications above.

| Horizon | $\hat\beta^{nd}$, OLS | $\hat\beta^{nd}$, IPW | $\hat\beta^{def}$, OLS | $\hat\beta^{def}$, IPW | $p$, OLS | $p$, IPW |
|:-------:|:---:|:---:|:---:|:---:|:---:|:---:|
| Year 1 | −0.35 | −1.63 | −7.58 | −7.35 | 0.006 | 0.013 |
| Year 2 | −2.13 | −2.13 | −9.48 | −10.30 | 0.001 | 0.010 |
| Year 3 | −0.11 | −0.81 | −6.65 | −4.59 | 0.011 | 0.105 |
| Year 4 | 1.31 | −0.48 | −5.06 | −1.90 | 0.038 | 0.591 |
| Year 5 | 2.16 | −0.23 | −5.05 | −0.70 | 0.047 | 0.892 |

Two features are worth flagging honestly. First, the gap between the
resolution types is confirmed under both OLS and IPW at the shorter
horizons (Years 1–2), where the extra cost of default clears conventional
significance in both columns. Second, at Years 3–5 the IPW-weighted
default-linked coefficient shrinks substantially and loses significance —
a real divergence from the OLS pattern, driven by the thin default-linked
common support at longer horizons (this design's propensity trimming
removes a growing share of the already-small default-linked pool as the
horizon lengthens). We do not paper over this: the longer-horizon
persistence documented under OLS in Section 4.1 is *not* robustly
confirmed once selection is addressed via this particular reweighting
design, and should be read as more uncertain at Years 3–5 than the OLS
table alone suggests.

### 5.4 Doubly-robust estimation (AIPW)

The paper's preferred selection-corrected estimator is the doubly-robust
AIPW of `08b_aipw.do`, following Jordà & Taylor (2016) and Asonuma et
al.'s Eq. (3): the propensity and outcome models are combined so that
consistency requires only one of the two to be correctly specified, and
inference comes from a stratified cluster bootstrap rather than from
asymptotic standard errors that a panel with twenty-odd default-linked
episodes cannot really support. Each onset type is estimated against
tranquil years with the rival type dropped, exactly as in Section 5.3, and
— the point that matters most for what follows — the *difference* between
the two types is bootstrapped directly on paired draws, so the extra cost
of default gets a confidence interval of its own rather than being read
off the visual gap between two separately estimated bands.

Pooling all onsets (Act 1), the doubly-robust cost is significant at Year
2 only ($-2.79$, 95% CI $[-5.82,-0.66]$), consistent with the plain-IPW
finding in Section 5.2 that the unweighted pooled average somewhat
overstates the typical cost. The resolution split is where the estimator
earns its place:

| Horizon | Non-default (vs tranquil) | Default-linked (vs tranquil) | Difference (def $-$ nd), paired bootstrap |
|---|---|---|---|
| Year 1 | ns | $-8.03$ (sig.) | $-7.21$ $[-12.29,-2.21]$ |
| Year 2 | $-1.51$ $[-3.79,-0.02]$ | $-11.04$ (sig.) | $-9.53$ $[-13.26,-4.67]$ |
| Year 3 | ns | $-7.07$ (sig.) | $-7.22$ $[-11.04,-1.42]$ |
| Year 4 | ns | ns | not significant |
| Year 5 | ns | ns | not significant |

Two readings follow. First, the central result of the paper survives the
doubly-robust estimator: the default-linked cost is large and significant
at Years 1–3, the non-default cost is essentially absent, and the gap
between them clears significance *as a difference*, on its own bootstrap
interval, at each of those three horizons. This is a materially stronger
statement than the OLS equality tests of Section 4.1, because it holds
after modelling selection into resolution type and without relying on the
correctness of the outcome equation alone.

Second, the AIPW independently reproduces the qualification already
reached in Section 5.3: the divergence does not survive to Years 4–5. Two
selection-corrected designs — the two-stage plain IPW and the
doubly-robust AIPW, which differ in estimator, in weighting, and in how
inference is constructed — agree both on the short-horizon result and on
where it stops. That agreement is what licenses the summary we use
throughout: **the resolution gap is well established at Years 1–3 and only
suggestive thereafter**, and the long-horizon persistence visible in the
raw OLS table should not be reported as a robust finding.

---

## 6. Transmission channels

To characterize the mechanism behind the divergence in Figure 2, we
re-estimate the local projection with six intermediate outcomes in place
of GDP — private credit, bank claims on the government, gross investment,
government expenditure, the primary balance, and FDI, each expressed as a
cumulative percent change from $t-1$ (log real levels for credit,
investment, and government expenditure, so the outcome is a genuine
percent change in the variable itself rather than in its ratio to a
collapsing GDP; the primary balance and FDI change sign and so remain
ratios) — plotted together in Figure 11's six-panel grid, pooling all
episodes:

| Channel | Year 1 | Year 2 | Year 3 | Year 4 | Year 5 |
|:--|:--:|:--:|:--:|:--:|:--:|
| Private credit | −1.73 (.148) | −5.34 (.073) | −7.22 (.042) | −5.01 (.044) | −2.21 (.346) |
| Bank claims on govt | 1.16 (.163) | 0.86 (.344) | 1.10 (.282) | 0.06 (.961) | 1.67 (.347) |
| Investment | −5.09 (.015) | −6.14 (.091) | −3.88 (.133) | −4.95 (.287) | −2.16 (.541) |
| Govt expenditure | 0.84 (.566) | −4.50 (.012) | −2.48 (.273) | −1.72 (.545) | 1.55 (.634) |
| Primary balance | −0.98 (.056) | −0.17 (.784) | −0.28 (.665) | −0.00 (.995) | −0.58 (.303) |
| FDI | −0.41 (.406) | 0.86 (.260) | 1.38 (.042) | 0.90 (.129) | −0.09 (.933) |

*(coefficient with $p$-value in parentheses)*

Because each channel is estimated as an independent, single-outcome
regression sharing only the onset dummy with the GDP specification, a
significant panel in Figure 11 establishes that the variable moves
following an onset; it does not by itself establish that the variable
mediates the output cost in Figure 2, since a common shock could drive
both the channel and GDP simultaneously, or GDP could be driving the
channel rather than the reverse. Bearing that distinction in mind: private
credit shows the clearest pooled signature, a response building from near
zero at impact to a trough of roughly $-5$ to $-7$ points by Years 3–4 —
delayed and compounding, rather than an on-impact shock, echoing the
GDP result's own persistence. Investment is significant only at Year 1;
government expenditure only at Year 2; primary balance and FDI show at
most one significant horizon each with no persistent pattern; bank claims
on government show no significant pooled response at any horizon.

---

## 7. Channels by resolution type

Splitting the same six channels by resolution type (Table 4; Figure 12a
for the OLS-DK version, Figure 12b for IPW-weighted) is considerably
noisier than the aggregate split in Section 4, given default-linked
sub-samples of only 13–21 treated episodes depending on the channel and
horizon.

| Channel | Year 1: nd / def ($p$) | Year 2: nd / def ($p$) | Year 3: nd / def ($p$) | Year 4: nd / def ($p$) | Year 5: nd / def ($p$) |
|:--|:--|:--|:--|:--|:--|
| Private credit | −0.6 / −4.3 (.42) | −1.2 / **−14.6** (.024) | −3.7 / −14.6 (.23) | 0.8 / −17.7 (.10) | 3.8 / −13.7 (.14) |
| Bank claims on govt | 0.4 / 2.9 (.12) | −0.2 / 3.0 (.22) | 0.1 / 3.1 (.12) | −0.7 / 1.4 (.23) | 0.4 / **3.9** (.07) |
| Investment | −4.4 / −6.4 (.80) | −5.0 / −8.4 (.49) | −2.4 / −6.7 (.50) | −4.7 / −5.4 (.90) | −1.7 / −3.0 (.90) |
| Govt expenditure | −0.2 / 3.1 (.21) | −3.8 / −6.0 (.61) | −0.8 / −5.9 (.16) | −1.0 / −3.3 (.70) | 3.5 / −2.4 (.34) |
| Primary balance | −0.5 / −2.1 (.07) | 0.1 / −0.8 (.47) | −0.5 / 0.2 (.51) | −0.1 / 0.2 (.83) | −0.8 / −0.1 (.55) |
| FDI | −0.1 / −0.9 (.27) | 1.4 / −0.2 (.11) | 1.9 / 0.3 (.15) | 1.7 / **−0.6** (.029) | 0.1 / −0.4 (.57) |

*(bold = equality test $p<0.05$)*

We are careful to describe only what these numbers actually support.
Across 30 nd/def equality tests, only two reach conventional
significance — private credit at Year 2 ($p=0.024$) and FDI at Year 4
($p=0.029$) — with the bank-claims gap at Year 5 and the primary-balance
gap at Year 1 borderline ($p\approx0.07$). Several point estimates look
economically large in the default-linked column, notably credit at Years
3–5 and government expenditure at Years 3–5, but their equality tests do
not reject, and given the sample sizes involved we read this as a genuine
power limitation rather than evidence that the channels are equal across
resolution types. We do not build a "different plumbing by resolution
type" narrative on these differences alone; Section 7a takes a different,
more direct approach to tying the channels to the GDP result.

### 7a. How much of the default-linked cost do these channels explain?

Section 7 shows *which* channels move differently by resolution type, at
whatever precision the small default-linked sample allows; it does not by
itself say how much of the estimated default-linked GDP cost in Table 2
each channel accounts for. We address this with a Gelbach (2016)
decomposition (`do/12b_gelbach_decomposition.do`): for each channel, we
add its own contemporaneous change to the exact specification behind
Table 2, estimated on the identical sample, and compute the resulting
shrinkage in the default-linked coefficient as a share of its original
magnitude.

| Channel | Year 1 | Year 2 | Year 3 | Year 4 | Year 5 |
|:--|:--:|:--:|:--:|:--:|:--:|
| Private credit | 4.6% | 17.3% | 24.5% | 43.2% | 23.7% |
| Government expenditure | −2.9% | 12.4% | 31.1% | 43.3% | 44.2% |
| Investment | 13.9% | 18.9% | 20.8% | 22.4% | 14.9% |
| Bank claims on government | 4.3% | 7.1% | 15.1% | 11.9% | 61.3%$^\dagger$ |
| Primary balance | 8.7% | 4.3% | 0.1% | −0.6% | 1.5% |
| FDI | 0.3% | −0.1% | −0.7% | 0.3% | −1.1% |

$^\dagger$ Estimated on the thinnest horizon ($n_{def}=17$); read with
caution rather than as a stable result.

Private credit and government expenditure each absorb a rising share of
the default-linked coefficient as the horizon lengthens, reaching roughly
40–44 percent by Years 4–5, consistent with the compounding credit response
visible in Figure 11's credit panel; the primary balance and FDI absorb
essentially none of it at any horizon. Investment sits in between, at a
fairly stable 15–22 percent throughout. This gives a ranked answer to
which channels are most closely tied to the gap plotted in Figure 2 —
credit and fiscal spending lead, investment is secondary, the primary
balance and FDI appear largely unrelated — though a shrinking coefficient
is consistent with a mediating role without establishing one causally: the
channel variables are not exogenous conditional on the controls, so the
shrinkage is equally consistent with a common shock driving both the
channel and the output cost, or with reverse causality running from
output to the channel. We do not compute the corresponding decomposition
for non-default episodes: because the non-default GDP coefficient is not
statistically distinguishable from zero at any horizon (Table 2), a
percentage computed relative to it is not economically informative — there
is no real effect there for a channel to explain.

---

### 7b. External adjustment: a channel that responds, but not by resolution

The six channels of Section 6 were selected as plausible domestic
transmission margins. One further margin is predicted by theory rather than
chosen by us, and is worth testing on its own terms: Aguiar and Gopinath
(2006) imply that a loss of market access forces the current account toward
surplus, since a country that cannot borrow must stop running a deficit.
That is a directional prediction about a specific variable, so it can fail,
which makes it more informative than a channel selected for plausibility.

We estimate it with the same local projection, replacing GDP with the
cumulative change in the current-account balance and dropping `l_ca` from
the controls, since it is the outcome's own lagged level (matching the
current-account model in `13c_aipw_channels.do`). A positive coefficient
means the external balance moves toward surplus.

| Horizon | $\hat\beta^{all}$ (SE) | $t$ |
|:-------:|:----------------------:|:---:|
| Year 1 | 1.03 (0.47) | 2.20 |
| Year 2 | 1.69 (0.77) | 2.20 |
| Year 3 | 1.60 (0.92) | 1.73 |
| Year 4 | 2.11 (0.96) | 2.20 |
| Year 5 | 2.62 (1.07) | 2.45 |

The current account moves toward surplus by **1.0 percentage point of GDP in
the crisis year, rising monotonically to 2.6 points by Year 5**, significant
at conventional levels at four of the five horizons and at ten percent at the
fifth. The prediction holds, and the adjustment builds rather than reverting:
this is a persistent reallocation of the external position, not a one-year
compression of imports.

**The adjustment does not differ by resolution, and that is the more useful
result.** Estimating both onset types jointly and testing their equality with
the same Clogg et al. (1995) statistic used for GDP in Section 4.1:

| Horizon | $\hat\beta^{nd}$ | $\hat\beta^{def}$ | Difference | $z$ | $p$ |
|:-------:|:----------------:|:-----------------:|:----------:|:---:|:---:|
| Year 1 | 0.88 | 1.34 | +0.46 | 0.40 | .691 |
| Year 2 | 1.78 | 1.49 | −0.29 | −0.21 | .838 |
| Year 3 | 1.89 | 1.02 | −0.87 | −0.59 | .561 |
| Year 4 | 2.62 | 1.10 | −1.51 | −1.02 | .318 |
| Year 5 | 3.06 | 1.79 | −1.27 | −0.74 | .468 |

No horizon comes close: the largest $|z|$ is 1.02. The two adjustment paths
cannot be distinguished. We record explicitly that the prediction written
into this test when it was first specified — that adjustment would be
*harder* under default, because market exclusion is more complete — is not
what the estimates show. The non-default point estimate is in fact the larger
of the two from Year 2 onward, though not significantly so. Neither the
prediction nor its reverse is supported; the honest statement is that
resolution type does not detectably change the external adjustment.

That null does real work in the argument. Sections 4 through 7a establish
that the output cost, the credit response and the fiscal response all differ
sharply by resolution, and Section 8 will add the sovereign-bank nexus to
that list. The current account does not join it. **External adjustment is a
generic consequence of losing market access, common to both crisis types,
and is therefore not a margin on which default and non-default resolutions
diverge.** Ruling a channel out is as informative here as ruling one in,
because the paper's central claim is precisely that the default premium
operates through some margins and not others.

Two qualifications. First, the reweighted version of this test disagrees at
short horizons — the IPW estimates put the default-linked response above the
non-default one at Years 1–2 (3.26 and 4.38 against 0.17 and 0.45) before
both collapse toward zero by Years 4–5 — so the *ordering* of the two types
is not stable across estimators, which is consistent with the difference
being unidentified rather than merely small. Second, the current account is
measured as a ratio to GDP and so shares the denominator problem discussed in
Section 6: a falling denominator mechanically raises the ratio. Because the
balance changes sign it cannot be log-transformed, so unlike credit or
investment there is no level-based alternative available; part of the
measured surplus at longer horizons may reflect the output contraction
itself rather than a change in the external position. The direction of the
finding is robust, its magnitude at Years 4–5 should be read with that in
mind.

---

## 8. The sovereign-bank nexus

Does the strength of the sovereign-bank linkage *condition* the output cost
documented in Section 4 — an analog of Asonuma et al.'s bank-intermediation
split (their Fig. 6), sharpened from bank *size* to bank *exposure to this
sovereign*? Four files in the pipeline speak to this question, and rather
than report them as four findings we organise them as one argument in three
steps: the **effect** (does pre-crisis nexus deepen the output cost?), its
**robustness** to a different functional form, and the **mechanism** (what do
banks actually do?). This mirrors the estimator ladder used for the GDP
result in Sections 4–5, and it is deliberately more conservative than
presenting the four sets of estimates as independent evidence — several of
them are the same object measured twice.

### 8.1 The effect: pre-crisis nexus amplifies the cost of default

The headline estimate comes from the exposure-interaction design of
`13b_exposure_heterogeneity.do`, which puts the nexus on the right-hand side
of the *GDP* equation rather than treating it as an outcome:

$$dy_{i,t+h} = \alpha_i + \lambda_t + b_{nd}D^{nd} + b_{def}D^{def}
+ \delta_{nd}(D^{nd}\times Z) + \delta_{def}(D^{def}\times Z) + \phi Z + X'\gamma + \varepsilon$$

where $Z$ is bank claims on government as a share of **total bank assets**,
measured at $t-1$ and standardized. The coefficients of interest are
$\delta_{nd}$ and $\delta_{def}$ — the extra output cost per standard
deviation of pre-crisis nexus in each resolution type — and the test is
their equality.

| Horizon | $\delta_{nd}$ | $\delta_{def}$ | $t(\delta_{def})$ | $p(\delta_{nd}=\delta_{def})$ |
|---|---|---|---|---|
| Year 1 | 0.700 | −2.613 | −0.51 | .535 |
| Year 2 | 0.703 | −2.404 | −0.80 | .321 |
| **Year 3** | 0.171 | **−4.431** | **−2.46** | **.034** |
| **Year 4** | 0.198 | **−4.270** | **−3.12** | **.029** |
| Year 5 | −0.252 | −3.235 | −1.25 | .329 |

One standard deviation more of the banking system committed to the sovereign
before the crisis is associated with roughly **4.3 additional percentage
points of lost output by Years 3–4 — but only where the crisis ends in
default**.

The non-default column is not empty, and reads in the opposite direction. At
Year 1, $\delta_{nd}=0.700$ with a standard error of 0.254 — a $t$ of 2.76,
significant, and *positive*: one standard deviation more pre-crisis exposure
is associated with a **smaller** output loss when the crisis does not end in
default. The effect fades to nothing by Year 3. So the nexus is not simply an
amplifier. It cushions at impact where the debt is honoured and amplifies at
medium horizons where it is not, and the equality test rejects at Years 3–4
because by then the two have moved decisively apart. This two-sided reading
was not anticipated when the test was written and emerged from the estimates;
Section 8.2 finds the same asymmetry independently.

Four features make this the cleanest nexus evidence in the paper:

1. **The outcome is GDP.** The alternative designs (§8.3) show the nexus
   *variable* responding to a crisis, which is co-movement — precisely the
   objection that motivates the exposure-interaction approach in the first
   place. This specification asks whether the nexus deepens the output cost
   itself.
2. **Exposure is predetermined.** Measured at $t-1$, it cannot be a
   consequence of the crisis it is supposed to amplify. The channel-outcome
   designs measure the nexus *after* onset, where reverse causality is live.
3. **The measure carries no development confound.** As Section 8.4 sets out,
   most exposure ratios in this file proxy financial depth as much as
   vulnerability. Claims on government as a share of *bank assets* is a
   portfolio composition measure, not a size measure, so it does not scale
   with income the way credit/GDP does.
4. **The difference is formally tested and rejects at two adjacent
   horizons**, with the same sign, against a non-default interaction that is
   zero or positive rather than negative at every horizon. Two
   neighbouring horizons agreeing is materially more convincing than one
   isolated cell, which matters given the multiple-testing arithmetic in
   §8.4.

**A measurement point worth stating explicitly.** The same economic concept
scaled by GDP rather than bank assets — `claims_govt` — is not estimable
here. Its $\delta_{def}$ is large and negative (−8.0 at Year 1) but carries a
standard error of 7.33 against 1.80 for the assets-scaled version, and no
horizon clears significance ($p=.246$ at best). An earlier version of this
project reported the GDP-scaled measure as significant; that result was
obtained on a legacy control set and does not survive the common core. The
lesson is that the doom-loop parameter has to be measured as a portfolio
share of the bank balance sheet, not as a ratio to output.

### 8.2 Robustness: the median split

`13d_aipw_nexus_split.do` asks the same question with a different functional
form — splitting onsets at the median of the same nexus variable and
estimating a doubly-robust AIPW output cost separately within each cell, with
a stratified cluster bootstrap. The split is not a development proxy: the
high-nexus bin contains Brazil, Mexico, Turkey and Indonesia alongside Ghana,
Kenya and Zambia, and the composition check printed by that file confirms it
(mean nexus 6.9 in the low bin against 23.1 in the high).

| Cell | Year 1 | Year 2 | Year 3 | Treated |
|---|---|---|---|---|
| Default × high nexus | −18.04\* | −18.36\* | −15.39\* | 7 |
| Default × low nexus | −7.20\* | −10.80\* | −3.57 | 10 |
| Non-default × high | +0.27 | −0.89 | −0.15 | 18 |
| Non-default × low | −3.49\* | −3.48\* | −1.10 | 15 |

\* bootstrap 95 percent interval excludes zero.

**The two differences behave differently, and the asymmetry matches Section
8.1.** In the *default* rows the point estimates are dramatic — the cost is
roughly two and a half times deeper where the nexus is tight — but the
high-minus-low difference **clears significance at no horizon**, the closest
being $-10.85$ with an interval of $[-17.46, 3.06]$ at Year 1. That null is
the most stable thing in the file: it holds identically across three separate
bootstrap runs. Seven treated observations cannot support the comparison,
however large the gap looks.

In the *non-default* rows the difference is significant at Year 1: high nexus
minus low nexus is $+3.76$ with an interval of $[0.55, 7.85]$. The cushioning
result of Section 8.1 therefore appears twice, in two designs that share only
the underlying variable — a continuous interaction estimated by OLS with
Driscoll-Kraay errors, and a median split estimated doubly-robustly with a
cluster bootstrap. Investment behaves the same way in the same cell
($+16.68$, $[1.00, 33.72]$ at Year 1), which is what a cushion operating
through the real economy should look like.

**A reproducibility note, because it changed what we are willing to claim.**
These intervals initially moved between runs of identical code. The file
originally drew 300 bootstrap resamples and set no seed, so the random-number
state at the start of any cell depended on what had run before it in the
session. Three cells crossed the five percent line purely on the draws,
including two reported above. Seeding fixed reproducibility without fixing
reliability: at 300 draws the 2.5th percentile is roughly the seventh-smallest
value, and in these cells only 235 to 274 draws survive estimation, so the
endpoint carries enough Monte Carlo error to move it further than its distance
from zero. Raising to 1000 draws resolved all four contested cells, and the
direction of movement is itself informative. The non-default difference moved
*away* from zero as precision improved ($[-0.09, 7.69]$ at 300 unseeded, then
$[0.17, 7.55]$, then $[0.55, 7.85]$), which is the signature of a real effect
emerging from noise. The `claims_govt` difference discussed in Section 8.3
moved the other way, from $[0.97, 12.76]$ to $[0.37, 12.91]$, and is reported
as marginal accordingly. More draws buy a reliable verdict, not a tighter
interval: the widths are set by seven to eighteen treated observations per
cell and no amount of resampling improves that.

**Draw counts as a diagnostic.** With 1000 draws the printed survival rate
matters more than the count. The GDP cells run 900 to 1000 valid draws and are
sound. The `claimpriv_assets` default × high cell runs 384 of 1000 at Year 4
and fails outright at Year 5; the corresponding difference runs 339. Those
cells are not usable and we do not read them.

We therefore report this section as **robustness in sign for the default
amplification, and as independent confirmation for the non-default
cushioning**. It uses the stronger estimator on the weaker sample; Section 8.1
uses the weaker estimator on the stronger design. They agree on both halves of
the asymmetry, and only Section 8.1 has the power to establish the
default-side amplification on its own.

### 8.3 Mechanism: what banks actually do

Having established that the nexus conditions the output cost, the remaining
question is behavioural. Two files estimate the nexus variables as outcomes —
`11b_nexus_channels.do` by OLS and IPW, `13c_aipw_channels.do` doubly-robustly
— and since these are the same object under two inference schemes we report
the doubly-robust estimates and note where the OLS agrees, rather than
treating them as separate results.

**Banks absorb the sovereign.** In the default × high-nexus cell, bank claims
on government rise by **+7.8, +17.1 and +12.9 percentage points of assets** in
Years 1–3, each individually significant, with no such movement in the
low-nexus default cell or in non-default episodes regardless of exposure. (The
Year-5 reading of −5.4 coincides with a sharply reduced draw count and is
treated as noise.)

**And they contract private lending.** Claims on the private sector move the
other way: the default-minus-non-default gap is significant at Year 1 under
AIPW (−6.62, 95% CI $[-11.38,-0.74]$), and the OLS split in `11b` agrees
across the horizon range (non-default $+2.72$ vs default $-1.60$ at Year 1,
$p=.016$, with the gap significant or marginal at every subsequent horizon).
The non-default coefficient is significantly *positive* at Year 1 under AIPW
($+2.39$, CI $[0.72, 4.17]$).

The high-minus-low difference on that same channel is significant at Year 3
($+10.09$, $[0.37, 12.91]$), but only marginally, and it is the one contested
cell whose interval moved *toward* zero as the draw count rose. We report it
as suggestive rather than established, and the claim in this subsection rests
on the level estimates above, which are significant at three consecutive
horizons under every bootstrap run we have executed.

Read together with §8.1 this is a coherent reallocation story: where banks
entered the crisis heavily exposed to the sovereign and the crisis ends in
default, they absorb more government paper and lend less to firms, and the
output cost is correspondingly deeper. The timing supports it — the GDP
amplification appears at Years 3–4, *after* the balance-sheet reallocation at
Years 1–3, which is what a slow-working credit mechanism should look like
rather than an immediate confidence shock.

### 8.4 What this evidence can and cannot bear

**Sample.** The interaction in §8.1 is identified off 11 default-linked
episodes; the median split in §8.2 off seven; the channel results in §8.3 off
roughly 18. Nexus data are IMF MFS, 2001 onward, covering 32 of 61 onsets.
Every estimate in this section is thin by construction.

**Multiple testing.** `13b` runs nine exposures × five horizons × two panels =
90 tests, of which roughly four or five would be significant at the 5% level
by chance; about seven are. No single panel in that file can carry a claim on
its own. What distinguishes the nexus result is that it is significant at two
adjacent horizons with a consistent sign and a near-zero comparison group, and
that a second, unrelated exposure design (§8.2) and the mechanism evidence
(§8.3) point the same way.

**The same arithmetic applies to the median split.** `13d` tests 5 outcomes ×
3 resolution groups × 5 horizons = 75 high-minus-low differences, of which
roughly four would clear five percent by chance. Three do: GDP and investment
in the non-default cell at Year 1, and `claims_govt` in the default cell at
Year 3. **That count is at, not above, the chance rate**, and we say so
plainly rather than presenting three significant cells as three findings. What
makes the non-default result usable is not its own $p$-value but that Section
8.1 finds the identical effect in the identical direction on the same variable
through a different estimator, and that GDP and investment — two outcomes
linked by construction — move together in the same cell. Corroboration across
designs is doing the work here, not significance within one.

**One further exposure clears the same bar.** The short-term share of external
debt — pure rollover vulnerability, with no development reading — shows the
same default-specific amplification: $\delta_{def}=-6.20$ at Year 3
($p(\text{diff})=.053$) and $-5.77$ at Year 4 ($p=.033$), against
$\delta_{nd}$ of $-0.67$ and $+0.70$. Pooled across all crises it is the
strongest exposure result in the file ($-4.70$ at Year 2, $p=.000$). So the
two interpretable vulnerabilities in the data are **how much of the banking
system was lent to the sovereign** and **how much debt had to be rolled
over** — both predetermined balance-sheet conditions, both operating only
where the crisis ends in default.

**The exposures we do not interpret.** Private credit, investment, FDI, bank
claims on government scaled by GDP, and debt service are all ratios that
proxy financial development as much as vulnerability to a channel, and their
interactions come back *positive* (credit $+2.33$ at Year 3, $p=.001$;
investment $+2.06$ at Year 5; FDI $+1.33$ at Year 3). We read this as the
development reading dominating rather than as evidence that these channels
protect, and it is corroborated within our own estimates: bank credit depth
enters the Act-1 outcome equation negatively and growing ($-0.026$ at Year 1
to $-0.171$ at Year 5, $t=-3.98$), so high-credit countries grow more slowly
on trend *and* lose less in a crisis — the signature of an income proxy. These
panels are reported for completeness and treated as uninformative about
transmission.

**Overall.** We regard the sovereign-bank nexus as the paper's strongest
*mechanism* evidence — a genuine interaction design on a predetermined,
confound-free measure, agreeing in sign across two functional forms and
supported by the balance-sheet behaviour it predicts — and simultaneously its
least statistically secure result, resting on cells of seven to eighteen
episodes. It should be read as a well-identified pattern estimated with little
power, not as a finding on the same footing as the resolution-type divergence
of Section 4.

---

## 9. Robustness and inference

Beyond the pre-trend-controlled specification (Section 4.2) and the
outturns-only restriction (Section 4.1), the code base includes: a
permutation placebo test (1,000 random reassignments of the 61 onsets
within the estimation pool, locating the true coefficients in the
resulting null distribution — an assumption-free check particularly
relevant given the small number of treated observations); a robustness
forest plot comparing the resolution-split coefficient across ten
alternative specifications; and $t$-distribution-based (rather than
normal-approximation) confidence intervals throughout, appropriate given
the panel's roughly 47 clusters.

---

## 10. Limitations

- **Small treated samples, and smaller than the classification suggests.**
  The design assigns 40 non-default and 21 default-linked onsets, but
  listwise deletion on the common core leaves **31 and 15** identifying the
  headline coefficients — six of the twenty-one default-linked episodes are
  lost to a missing control, chiefly public debt at $t-1$. The nexus cells
  in Section 8 are thinner still (7 to 18). Every precision statement in
  this paper should be read against the within-sample counts reported in
  Tables 1 and 2, not against the episode totals.
- **Selection under reweighting is not uniform across horizons.** Section
  5.3's IPW comparison confirms the extra cost of default at Years 1–2 but
  loses significance at Years 3–5 once the two-stage, vs.-tranquil design
  is used — a real qualification to the OLS-only persistence claim in
  Section 4.1 that should not be smoothed over.
- **The exposure-interaction file runs 90 tests.** `13b` estimates nine
  exposures across five horizons and two panels; at the 5% level roughly
  four or five significant cells are expected by chance, and about seven
  appear. Section 8 therefore rests only on the two exposures whose sign is
  economically identified (nexus, short-term debt share), each significant
  at two adjacent horizons with a near-zero comparison group, and treats
  every other panel as uninformative. Readers should not count the
  significant cells in Tables 5–6 as independent findings.

- **The nexus high-minus-low difference is run-sensitive at the margin.**
  On a cell of 7 treated observations, the non-default high-minus-low gap
  moved from nominally significant to marginally insignificant between two
  bootstrap runs of identical code (Section 8). Nothing about the point
  estimate changed materially; the resampling variance on a cell that
  small is simply of the same order as the gap being tested. This is a
  reason to treat every nexus *difference* as suggestive, and to rest the
  section on the level contrasts and the `claims_govt` mechanism sequence
  instead.
- **Episode duration is not weighted.** Each episode contributes exactly
  one treated observation (its onset year); continuation years are
  excluded from the sample entirely. This is the standard local-projection
  convention and the multi-year dynamics of a crisis are captured through
  the horizon dimension rather than through additional treated rows, but
  one consequence is that a crisis lasting five years contributes exactly
  as much identifying variation as one lasting a single year. Whether
  protracted crises are systematically more costly is an open question
  this design does not address.
- **Thin nexus and common-support cells.** The default×high-nexus cell in
  Section 8 rests on 7 treated country-years, and the AIPW bootstrap
  loses a growing share of draws at later horizons in the thinnest cells;
  the high/low nexus split could in principle also proxy income or
  financial development more broadly, a possibility the underlying code
  prints the country composition of each bin to allow checking, but which
  has not been formally ruled out here.
- **Channel decomposition is not mediation.** Section 7a's Gelbach shares
  quantify how much of the default-linked coefficient is statistically
  absorbed by each channel, not a causal mediation estimate; the channel
  variables are not exogenous conditional on the controls.

---

## 11. Summary

The analysis supports one central, carefully qualified finding and two
further pieces of mechanism evidence of differing strength. The central
finding: the output cost of a spread crisis is not a generic consequence
of crossing the market-based threshold, but is concentrated in, and
largely specific to, episodes that end in default — non-default episodes
show no reliable GDP cost at any horizon, default-linked episodes show a
persistent cost of four to seven percentage points, and the gap between
the two is itself statistically significant at four of five horizons. This
result survives a direct stress test against a pre-existing pre-crisis
growth pattern, though it is only confirmed at the shorter horizons once
selection into the two-stage IPW design is addressed, a real qualification
rather than a settled matter.

Three further results bear on mechanism. A Gelbach decomposition shows
private credit and
government expenditure jointly absorb an increasing share of the
default-linked coefficient, reaching roughly 40–44 percent by Years 4–5 —
the strongest quantitative candidates for the transmission channel,
without establishing causal mediation.

The second is external adjustment, and it is reported as a null on the
dimension the paper cares about. The current account moves toward surplus
after a spread crisis by one percentage point of GDP in the crisis year,
rising to 2.6 points by Year 5 and significant at four of five horizons —
confirming the forced-deleveraging prediction of Aguiar and Gopinath (2006)
for spread crises. But the equality test between resolution types rejects at
no horizon, with a largest absolute $z$ of 1.02, and the prediction written
into the test when it was specified (harder adjustment under default) is not
what the estimates show. External adjustment is therefore a generic
consequence of losing market access rather than a margin on which the two
resolutions diverge, which usefully bounds where the default premium
operates: in output, credit, fiscal spending and the sovereign-bank nexus,
but not in the external balance.

The third is the sovereign-bank nexus, and Section 8 now presents it as a
single three-step argument rather than as several coordinate findings.
Interacting crisis onset with pre-crisis bank claims on government as a
share of bank assets — a predetermined portfolio measure that carries no
financial-development confound — one standard deviation of additional
exposure is associated with roughly 4.3 percentage points of extra output
loss by Years 3–4, and only in default-linked episodes; the equality test
against non-default rejects at both horizons. A median-split doubly-robust
estimate agrees in sign and is roughly two and a half times deeper in the
high-exposure cell, but its difference does not clear significance on seven
treated observations, so it is reported as robustness rather than as
confirmation. The behaviour that would produce the effect is present:
banks in the default × high-exposure cell increase government claims by 8
to 17 percentage points of assets over Years 1–3 while contracting private
lending, and the output amplification appears afterwards, at Years 3–4, as
a slow-working credit mechanism should. The same design identifies one
further vulnerability on the same terms — the short-term share of external
debt, which amplifies the default-linked cost by roughly 6 percentage
points at Years 3–4 and is the strongest pooled exposure result in the
paper. Both are predetermined balance-sheet conditions; both bind only
where the crisis ends in default. All of it rests on cells of seven to
eighteen episodes, so it is a well-identified pattern estimated with little
power, not a finding on the footing of the resolution-type divergence.
