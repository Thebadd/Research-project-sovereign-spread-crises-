# Empirical Analysis — Structure and Regression Equations

This document describes every regression run in the project, in the order
the do-files execute them, with the exact specification and the reason each
one is there.

---

## Data construction (`01_build_panel.do`, `01b_merge_new_controls.do`)

The panel covers 52 emerging and frontier market economies over 1994–2025
at annual frequency. Continuation years (years 2, 3, … of an episode, after
the onset year) are dropped from every estimation sample. The identifying
variation is the **onset year** of each episode only.

**Outcome variable — same construction used in every LP:**

$$dy_{i,t+h} = \bigl(\ln\text{GDPpc}_{i,t+h} - \ln\text{GDPpc}_{i,t-1}\bigr) \times 100$$

Cumulative percentage change in real GDP per capita from one year before
the crisis to $h$ years after, for $h \in \{0,1,2,3,4\}$. For channel
regressions the same formula is applied to the channel variable instead
of GDP.

**Controls used in every baseline LP:**

| Variable | Content |
|---|---|
| `l1_gdpg`, `l2_gdpg` | GDP growth at $t-1$ and $t-2$ (captures pre-crisis momentum) |
| `ca` | Current account balance / GDP (external vulnerability) |
| `debt` | General government debt / GDP (fiscal space) |
| `infl` | Inflation (macro instability) |
| `imf` | IMF program dummy (crisis management) |
| `vix` | VIX index (global risk appetite) |
| `ust10y` | US 10-year rate (global financial conditions) |

---

## Step 0 — Pre-crisis balance (`05_balance_table.do`)

**What it does.** Before any regression, compare the observable
characteristics of the two groups at crisis onset:

- **Non-default episodes** (40 onsets)
- **Default-linked episodes** (21 onsets)

A two-sample $t$-test with unequal variances is run for each of:
`gdpg`, `l1_gdpg`, `debt`, `ca`, `infl`, `l_spr_mean`, `imf`, `fdi`,
`govexp`.

**Why.** If the two groups look statistically different before the crisis
on observables (higher debt, worse current account, etc. in the
default-linked group), then a naive comparison of their output losses
conflates the crisis *effect* with selection into default. The balance
table tells the reader how much of a problem this is and motivates
the IPW correction in later steps.

---

## Act 1 — Baseline LP: all crisis episodes (`02_lp_all.do`)

**Research question.** What is the average output cost of a sovereign
spread crisis, pooling all 61 onset episodes?

**Regression.** For each horizon $h \in \{-2,-1,0,1,2,3,4\}$:

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta_h \cdot D_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

- $\alpha_i$ = country fixed effect
- $\lambda_t$ = year fixed effect
- $D_{it}$ = `onset_all` (= 1 in first year of any spread crisis)
- SE: Driscoll-Kraay with bandwidth $= \max(1, h+1)$

**Horizons $h = -2$ and $h = -1$ (pre-trend placebo).** The outcome is
*past* GDP growth, so $\beta_{-2}$ and $\beta_{-1}$ should be
statistically indistinguishable from zero. A significant pre-trend
would indicate that the crisis dummy is picking up a pre-existing
trend rather than a causal effect.

**Output.** The sequence $\{\hat\beta_0, \hat\beta_1, \hat\beta_2,
\hat\beta_3, \hat\beta_4\}$ is the impulse response function (IRF) of
output to a spread crisis onset.

---

## Act 2 — LP by resolution type (`03_lp_resolution.do`)

**Research question.** Does the output cost differ between spread crises
resolved without default and those linked to debt restructuring?

### Spec A — Separate regressions (one group at a time)

**Non-default episodes** ($D^{nd}_{it}$ = `onset_nd`, 40 episodes):

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{nd}_h \cdot D^{nd}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

**Default-linked episodes** ($D^{def}_{it}$ = `onset_def`, 21 episodes):

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{def}_h \cdot D^{def}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

These give separate IRFs for each group, estimated on the full panel
(both groups serve as controls for each other within year and country).

### Spec B — Joint regression with equality test

Both dummies enter the same regression simultaneously:

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{nd}_h \cdot D^{nd}_{it} + \beta^{def}_h \cdot D^{def}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

Then test $H_0 : \beta^{nd}_h = \beta^{def}_h$ at each horizon with an
$F$-test.

**Why the joint regression matters.** The separate regressions cannot
tell you whether the *gap* between the two IRFs is statistically
significant. The joint regression does. If the gap widens and becomes
significant at $h = 2, 3, 4$, that is direct evidence that the
resolution outcome affects the depth and persistence of the output loss —
the key empirical finding of the paper.

---

## Act 3 — Role of IMF programs (`09_lp_imf.do`)

**Research question.** Does an IMF program at crisis onset attenuate the
output cost?

**Spec A — Separate LP:**

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{imf}_h \cdot D^{imf}_{it} + \gamma' X^{-imf}_{it} + \varepsilon_{i,t+h}$$

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{no\_imf}_h \cdot D^{no\_imf}_{it} + \gamma' X^{-imf}_{it} + \varepsilon_{i,t+h}$$

Note: `imf` is excluded from the controls since it is now the treatment
modifier.

**Spec B — Joint regression** with equality test $H_0:
\beta^{imf}_h = \beta^{no\_imf}_h$.

**Spec C — IPW correction.** Countries that get IMF programs are
systematically different (higher debt, worse current account, lower growth)
from those that do not. A probit of $\Pr(\text{IMF} = 1 \mid
\text{crisis onset}, X)$ generates propensity scores, which are used to
construct stabilized IPW weights that reweight non-IMF episodes to look
like IMF episodes on observables. The IPW-weighted LP is then estimated
with `areg` (country FE) and clustered standard errors.

**Why.** Without IPW, a finding that IMF episodes have smaller output
losses would be confounded by the fact that the IMF typically intervenes
in countries that are already in severe distress — the comparison is not
apples-to-apples.

---

## Act 4 — Transmission channels, all crises (`11_channels.do`)

**Research question.** Through which channels does a spread crisis
transmit to output? The LP outcome is replaced by six channel variables.

For each channel variable $Z$ and each horizon $h \in \{0,1,2,3,4\}$:

$$\Delta^h Z_{it} = \alpha_i + \lambda_t + \beta^Z_h \cdot D_{it} + \gamma^Z{}' X^Z_{it} + \varepsilon_{i,t+h}$$

where $\Delta^h Z_{it} = Z_{i,t+h} - Z_{i,t-1}$ (same anchoring as the
GDP outcome).

| Channel | $Z$ | What it tests | Controls |
|---|---|---|---|
| 1. Credit | Private credit / GDP | Does the banking sector contract? | `l1_gdpg l2_gdpg debt infl ca banking_crisis` |
| 2. Sovereign-bank nexus | Bank claims on govt / GDP | Do banks accumulate sovereign bonds (doom loop)? | `L.claims_govt L.credit pb banking_crisis` |
| 3. Investment | Gross investment / GDP | Does investment fall (capital channel)? | `l1_gdpg l2_gdpg debt ca L.credit banking_crisis` |
| 4. Govt expenditure | Govexp / GDP | Is there fiscal austerity? | `L.govexp debt revenue_gdp` |
| 5. Primary balance | PB / GDP | Does fiscal adjustment occur? | `l1_gdpg l2_gdpg debt ca L.pb banking_crisis` |
| 6. FDI | FDI / GDP | Does external financing stop? | `l1_gdpg L.fdi infl reer_chg` |

Each channel is also re-estimated with **IPW weights** (same probit
first stage as in `08_ipw_lp.do`) to check that the channel responses
are not driven by selection on pre-crisis observables.

**Why channel-specific controls?** VIX and UST10y are pure time-series
variables with zero cross-sectional variation — they are fully absorbed
by year fixed effects and would cause collinearity. They are dropped from
all channel regressions except FDI, where the US rate captures the carry
trade incentive of foreign investors and adds genuine identification
content not in the year FE.

---

## Act 5 — Channels by resolution type (`12_channels_resolution.do`)

**Research question.** Do the six transmission channels operate
differently in episodes that end in default versus those that do not?

**Spec A — OLS joint LP** (both dummies simultaneously, for each channel):

$$\Delta^h Z_{it} = \alpha_i + \lambda_t + \beta^{nd,Z}_h \cdot D^{nd}_{it} + \beta^{def,Z}_h \cdot D^{def}_{it} + \gamma^Z{}' X^Z_{it} + \varepsilon_{i,t+h}$$

Equality test $H_0: \beta^{nd,Z}_h = \beta^{def,Z}_h$ at each horizon.

**Spec B — IPW-corrected LP.** A probit of
$\Pr(\text{default-linked} = 1 \mid \text{crisis onset}, \text{debt}, \text{CA})$
generates Act 2 stabilized weights. Non-default episodes are reweighted
to match default-linked episodes on pre-crisis observables. Estimated
with `areg` and clustered SE.

**Why.** This is where the model predictions become testable:
- The credit-to-output ratio should be higher in non-default episodes
  (Proposition 2): credit contracts without the extra capital depletion
  that depresses output further in default episodes.
- Investment should fall more in default-linked episodes through the
  external financing stop (equation 6 of the model, $\chi > 0$).
- Bank sovereign bond accumulation (doom loop) may differ if banks
  respond to default risk rather than spread level alone.

---

## Act 6 — Mechanistic tests (`13_mechanisms.do`)

**Research question.** Are the credit and investment channel responses
causal mechanisms, or correlations driven by common factors?

### Test 1 — Supply vs. demand in the credit channel

Run the credit LP twice:

**Baseline:** $\Delta^h\text{Credit}_{it} = \alpha_i + \lambda_t + \beta_h \cdot D_{it} + \gamma' X^{cr} + \varepsilon$

**With sovereign bond stock:** same + $L.\text{ClaimsGovt}_{it}$

If adding $L.\text{ClaimsGovt}$ significantly shrinks $\hat\beta_h$,
sovereign bond accumulation absorbs bank balance sheet capacity
(portfolio substitution — supply-side contraction). If $\hat\beta_h$
is unchanged, the credit contraction is demand-driven or operates
through a different supply channel (funding costs, risk aversion).

### Test 2 — Credit as mediator of investment contraction

Run the investment LP twice:

**Without credit:** $\Delta^h\text{Inv}_{it} = \alpha_i + \lambda_t + \beta^{total}_h \cdot D_{it} + \gamma' X^{inv,-cr} + \varepsilon$

**With credit:** same + $L.\text{Credit}_{it}$

The **mediation share** at horizon $h$ is:

$$\text{Mediation share}_h = \frac{\hat\beta^{total}_h - \hat\beta^{direct}_h}{\hat\beta^{total}_h}$$

A mediation share close to 1 means almost all of the investment decline
passes through the credit channel. A small share means investment falls
directly (through the lending rate / working-capital channel) even
independently of the credit quantity contraction.

---

## Robustness (`06_robustness.do`)

All robustness checks replicate the Act 1 LP (`onset_all`) across ten
alternative specifications. Results are collected in a single matrix and
exported to `robustness_summary.csv`.

| Label | What changes |
|---|---|
| Baseline | Main spec (DK SE, all 61 onsets) |
| Crit1 only | Restrict treatment to 1000 bps threshold episodes |
| Crit2 only | Restrict treatment to HRT 99th-percentile threshold episodes |
| Drop Argentina | Removes the single country with the most extreme episodes |
| Drop Venezuela | Removes hyperinflation outlier |
| Drop ARG + VEN | Both together |
| Drop 2020–2021 | Excludes COVID pandemic years |
| Drop 2008–2009 | Excludes Global Financial Crisis |
| Drop 1997–1999 | Excludes Asian / EM crisis cluster |
| Country cluster | Replaces DK SE with standard within-country cluster SE |

**Why.** Each check addresses a specific identification concern: (i)
threshold sensitivity — does the result depend on which crisis definition
you use? (ii) outlier sensitivity — do a handful of extreme episodes
mechanically drive the IRF? (iii) common-shock sensitivity — are the
results explained by one global episode that hit many countries at once?
(iv) SE sensitivity — does the inference depend on the Driscoll-Kraay
correction?

---

## Placebo test (`07_placebo.do`)

**Procedure.** The 61 real crisis onsets are randomly reassigned to 61
country-year observations drawn from the estimation pool. The LP is run
on the placebo dummy. This is repeated 1,000 times, generating an
empirical null distribution of $\hat\beta_h$ under no treatment effect.

The true $\hat\beta_h$ from Act 1 is then located in this null
distribution and an empirical $p$-value is computed.

**Why.** The Driscoll-Kraay $t$-statistic relies on asymptotic
approximations that may be unreliable with a small number of treated
units ($N_{\text{treat}} = 61$). The permutation test is
assumption-free and provides an exact finite-sample $p$-value.

---

## IPW baseline LP (`08_ipw_lp.do`)

**Why a separate IPW file in addition to `11_channels.do`?** This file
applies IPW to the baseline *output* LP — not to the channels. The
procedure is:

1. **First-stage probit:** $\Pr(\text{onset}_{it} = 1 \mid X_{it-1})$ on
   the full set of macro controls, estimated on the estimation sample.
2. **Propensity scores** trimmed at $[0.01, 0.99]$ to avoid explosive
   weights.
3. **Stabilized IPW weights:**
   $w_{it} = \bar{p} / \hat{p}_{it}$ for treated, $w_{it} = (1-\bar{p})/(1-\hat{p}_{it})$ for control.
4. **Weighted LP** with `areg` (country FE) and clustered SE.

The IPW-LP estimate is plotted against the unweighted baseline. If they
are close, selection on observables is not driving the main result. If
they diverge, the IPW estimate is preferred.

---

## Heterogeneity (`10_heterogeneity.do`)

Two dimensions of heterogeneity, each with a joint regression and
equality test.

**Block A — Market classification.**

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{fr}_h \cdot D^{frontier}_{it} + \beta^{em}_h \cdot D^{EM}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

Test $H_0: \beta^{fr}_h = \beta^{em}_h$. Frontier markets have thinner
domestic financial systems and less investor base — the transmission
channels may be stronger or weaker.

**Block B — Episode duration.**

Episodes are split at the median into short ($\leq 2$ years) and long
($\geq 3$ years):

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{short}_h \cdot D^{short}_{it} + \beta^{long}_h \cdot D^{long}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

Test $H_0: \beta^{short}_h = \beta^{long}_h$. Longer episodes may
generate more capital depletion through sustained investment stop,
consistent with the model's prediction that output losses in the default
regime widen monotonically with horizon.

---

## Model calibration and IRF comparison (`14_calibration.do`, `15_solve_default.do`, `16_model_irf.do`)

These three files bridge the empirical and theoretical sections.

**`14_calibration.do`** sets all structural parameters:
- *Externally calibrated* from the EM-DSGE literature (Arellano 2008,
  Neumeyer-Perri 2005, Cruces-Trebesch 2013): $\beta, \beta', \sigma,
  \alpha, \delta, R^*, \theta, \mu, \lambda$.
- *Data-determined* from tranquil-period means in the panel (years with
  `onset_all == 0` and `continuation == 0`): steady-state ratios
  $b^B_{ss}/N_{ss}$, $\ell_{ss}/N_{ss}$, $s_{ss}$, $R^L_{ss}$.
- *Estimated from the panel*: spread persistence $\rho_s$, income
  process $(\rho_y, \sigma_y)$, working-capital share $\xi$, sovereign
  ceiling pass-through $\gamma$, external investment share $\chi$.

**`15_solve_default.do`** solves the sovereign default block (Block 1)
by value function iteration on a discretized $(B, y)$ grid. It generates
the equilibrium spread paths for both crisis types.

**`16_model_irf.do`** simulates the transmission block (Block 2)
conditional on the equilibrium spread paths from Block 1. It produces
model-implied IRFs for output, credit, and investment across both crisis
regimes and overlays them on the empirical IRFs from Acts 1 and 2.

**Why.** The comparison between empirical and model-implied IRFs is
the core test of the theoretical framework. The model is not fit to
match the IRFs — its parameters are calibrated independently. The
out-of-sample comparison then tells you whether the structural
mechanisms (balance-sheet channel, working-capital channel, capital
depletion, external financing stop) are quantitatively consistent with
the reduced-form local projection estimates.

---

## Execution order

```
00_master.do
│
├── 01_build_panel.do          ← data construction + LP outcome variables
├── 01b_merge_new_controls.do  ← merge additional macro controls
│
├── 02_lp_all.do               ← Act 1: baseline IRF, all crises
├── 03_lp_resolution.do        ← Act 2: non-default vs. default-linked
├── 04_graphs.do               ← all main figures
├── 05_balance_table.do        ← pre-crisis balance check
│
├── 06_robustness.do           ← 10 robustness specifications
├── 07_placebo.do              ← permutation test (1,000 draws)
├── 08_ipw_lp.do               ← IPW-corrected baseline LP
├── 09_lp_imf.do               ← Act 3: IMF program heterogeneity
├── 10_heterogeneity.do        ← market type + episode duration
│
├── 11_channels.do             ← Act 4: 6 transmission channels
├── 12_channels_resolution.do  ← Act 5: channels by resolution type
├── 13_mechanisms.do           ← Act 6: supply/demand + mediation tests
│
├── 14_calibration.do          ← structural parameter calibration
├── 15_solve_default.do        ← Block 1 VFI solution
└── 16_model_irf.do            ← Block 2 simulation + model vs. data
```
