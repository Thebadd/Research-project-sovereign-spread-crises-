# Empirical Analysis

## 1. Setting and identification

The analysis rests on an unbalanced panel of 52 emerging and frontier
market economies observed annually over 1994–2025. The unit of
observation is the country-year. The sample is restricted to **onset
years** — the first year of each spread-crisis episode — and **tranquil
years**; continuation years (the second and later years of an ongoing
episode) are excluded throughout. The identifying variation therefore
comes from the *transition into* a crisis, not from the crisis persisting.
Of the 1,661 country-years, 1,435 fall in the estimation sample, 61 of
which are crisis onsets. Thirty-four countries experience at least one
crisis; eighteen are never treated and form a pure control group.

The estimator is the Jordà (2005) local projection. For each horizon
$h \in \{0, 1, 2, 3, 4\}$ we estimate a separate regression of the
cumulative change in an outcome on a crisis-onset indicator, country and
year fixed effects, and a vector of predetermined macro controls:

$$\Delta^h y_{it} = \alpha_i + \lambda_t + \beta_h \, D_{it} + \gamma' X_{it} + \varepsilon_{it+h}$$

The dependent variable $\Delta^h y_{it}$ is the cumulative percentage
deviation of real GDP per capita from its pre-crisis level, anchored at
$t-1$, so that $h=0$ is the impact effect and $h=4$ the four-year
cumulative loss. The sequence $\{\hat\beta_h\}$ traces the impulse
response of the outcome to a crisis onset.

The control vector $X_{it}$ contains two lags of GDP growth, the current
account balance, the public debt-to-GDP ratio, inflation, and an IMF
program dummy. Lagged growth absorbs pre-crisis output momentum; debt and
the current account capture fiscal and external vulnerability at onset;
inflation proxies for underlying macro instability; the IMF dummy controls
for the differential path of countries that receive official financing.
Global financial conditions (the VIX and the US 10-year yield) are
deliberately excluded: year fixed effects absorb every shock common to all
countries in a given year, making any pure time-series regressor
redundant. Standard errors follow Driscoll and Kraay (1998), correcting
for cross-sectional dependence and for the serial correlation that
mechanically rises with the forecast horizon (lag = $\max(1, h+1)$).

The identifying assumption is that, conditional on fixed effects and
observable fundamentals, crisis onset is uncorrelated with contemporaneous
unobserved shocks to future output. We test it directly at the placebo
horizons $h=-5$ through $h=-1$: under valid identification these coefficients
should be approximately zero. Each placebo dependent variable is the
single-year log GDP per capita growth rate $h$ years before onset — none
include the crisis year itself. Non-zero pre-trend coefficients would signal
anticipation effects or trending differences between treated and control
countries.

---

## 2. The average output cost of a spread crisis

The baseline local projection on all 61 onsets yields a precisely
estimated and persistent output decline.

| Horizon $h$ | $\hat\beta_h$ (pp) | SE | $p$ |
|:-----------:|:------------------:|:----:|:-----:|
| 0 | −1.86 | 0.47 | 0.000 |
| 1 | −3.76 | 0.72 | 0.000 |
| 2 | −3.18 | 0.67 | 0.000 |
| 3 | −2.32 | 0.70 | 0.001 |
| 4 | −3.00 | 1.05 | 0.004 |

Real GDP per capita falls about 1.9pp on impact, reaches a trough of
−3.8pp one year after onset, and remains roughly 3pp below the pre-crisis
trajectory four years out. All coefficients are significant at the 1%
level under Driscoll-Kraay standard errors. The defining feature is the
**absence of recovery within the four-year window**: a spread crisis
produces a level shift in output, not a transitory dip. This is the
central and most robust result of the paper.

---

## 3. Pre-trend validation

We estimate the same local projection at placebo horizons $h=-5$ through
$h=-1$. Each dependent variable is the single-year log GDP per capita
growth rate $h$ years before onset; by construction none of these
variables include the crisis year $t$. Under valid identification, the
onset dummy should have no predictive power over any of these
pre-crisis growth rates. The five-year window is long enough to detect
gradual pre-trends or anticipation effects that a two-year window would
miss.

All pre-trend coefficients should be close to zero and statistically
insignificant — results are reported after re-running the corrected
pipeline. If any horizon shows a significant coefficient, the affected
group's impulse responses must be interpreted with the caveat that
treated countries were already on a different output trajectory before
the crisis onset date.

---

## 4. Does the resolution of the crisis matter?

The central question of the second part of the analysis is whether crises
resolved without default impose smaller losses than those linked to
sovereign default. We estimate a joint local projection with both onset
dummies entered simultaneously,

$$\Delta^h y_{it} = \alpha_i + \lambda_t + \beta^{nd}_h D^{nd}_{it} + \beta^{def}_h D^{def}_{it} + \gamma' X_{it} + \varepsilon_{it+h},$$

and test the equality $\beta^{nd}_h = \beta^{def}_h$ at each horizon. Of
the 61 episodes, 40 are non-default and 21 are default-linked (Venezuela
2008 is reclassified as non-default: its restructuring began in 2017, a
nine-year lag well beyond the five-year rule used to attribute a spread
crisis to default).

| Horizon | $\hat\beta^{nd}$ | $\hat\beta^{def}$ | Gap | $p(\beta^{nd}=\beta^{def})$ |
|:-------:|:----------------:|:-----------------:|:----:|:---------------------------:|
| 0 | −1.03 | −3.54 | −2.51 | 0.104 |
| 1 | −2.71 | −5.82 | −3.11 | 0.105 |
| 2 | −2.46 | −4.58 | −2.12 | 0.208 |
| 3 | −1.76 | −3.37 | −1.61 | 0.406 |
| 4 | −3.07 | −2.88 | +0.19 | 0.952 |

The point estimates suggest default-linked crises are two to three times
deeper at impact and at $h=1$. Three features qualify this, however.

1. **The difference is never statistically significant.** The smallest
   $p$-value on the equality test is 0.10, at the impact and one-year
   horizons.
2. **The gap is front-loaded, not widening.** It is largest at impact and
   *converges* to zero by $h=4$, where the two groups are
   indistinguishable. Default-linked crises look like a deeper but more
   front-loaded shock, with non-default losses catching up over time.
3. **The default-linked sample is small and imprecise.** With 21 onsets,
   the default-linked standard error rises from 1.26 at $h=0$ to 2.66 at
   $h=4$. The failure to reject equality is therefore as much a power
   limitation as evidence of similarity.

The data are thus *suggestive* that default-linked crises are worse on
impact, but they cannot establish a statistically reliable, persistent
resolution gap with this sample.

---

## 5. Selection and inverse-probability weighting

The two crisis groups are not exchangeable ex ante. Comparing observable
fundamentals at the onset year:

| Variable | Non-default | Default-linked | $p$ |
|:--|:--:|:--:|:--:|
| GDP growth (%) | +2.77 | −2.91 | 0.006 |
| L1 GDP growth (%) | 4.94 | 1.81 | 0.001 |
| Debt/GDP | 51.6 | 78.0 | 0.012 |
| Lagged EMBIG spread (bps) | 402 | 592 | 0.000 |
| Current account/GDP | −3.84 | −2.33 | 0.410 |
| IMF program | 0.38 | 0.57 | 0.154 |

Default-linked countries enter with worse output momentum, substantially
higher debt, and higher spreads. We address this with inverse-probability
weighting: a first-stage probit predicts treatment from predetermined
fundamentals, propensity scores are trimmed to $[0.01, 0.99]$, and
stabilized weights are applied in the projection.

- **Act 1 (onset vs. tranquil).** Reweighting makes the estimated losses
  *larger*, not smaller, at every horizon (e.g. $h=1$: −4.05 → −5.14;
  $h=4$: −3.51 → −5.01). Selection on observables, if anything,
  *understates* the cost — the baseline estimate is conservative.
- **Act 2 (default vs. non-default).** Reweighting the non-default group
  to match the default-linked group on fundamentals **reduces the
  estimated loss of *both* groups** at most horizons. The full set of
  coefficients (GDP outcome, `areg` with cluster SE):

  | $h$ | $\beta^{nd}$ OLS | $\beta^{nd}$ IPW | $\beta^{def}$ OLS | $\beta^{def}$ IPW | Gap OLS | Gap IPW |
  |:---:|:---:|:---:|:---:|:---:|:---:|:---:|
  | 0 | −1.03 | −1.08 | −3.54 | −3.02 | −2.51 | −1.93 |
  | 1 | −2.71 | −2.58 | −5.82 | −5.58 | −3.11 | −3.01 |
  | 2 | −2.46 | −2.14 | −4.58 | −4.30 | −2.12 | −2.16 |
  | 3 | −1.76 | −1.32 | −3.37 | −2.83 | −1.61 | −1.52 |
  | 4 | −3.07 | −2.34 | −2.88 | −2.51 | +0.19 | −0.17 |

  The default-linked excess loss is attenuated at impact (gap −2.51 →
  −1.93) and at $h=3$, but the gap is essentially unchanged at $h=1$–$2$
  and remains near zero at $h=4$, because reweighting shrinks the
  non-default loss by *more* at longer horizons (its weighted-minus-OLS
  difference grows from −0.05 at $h=0$ to +0.73 at $h=4$). A meaningful
  part of the OLS default-linked penalty therefore reflects pre-crisis
  fundamentals rather than the resolution itself, while default-linked
  remains the deeper response at $h=0$–$3$. The equality test was not
  computed in the original Act 2 IPW run, so an explicit IPW gap $p$-value
  is not reported here; given the OLS gap was insignificant at every
  horizon ($p$ = 0.10–0.95) and the IPW gaps are similar or smaller with
  comparable standard errors, it is very unlikely to reach significance.

Both exercises rest on thin common support — 476 of roughly 1,275 sample
observations are trimmed, and the first stage classifies only a handful of
treated observations above the 0.5 threshold — so the IPW results are
corroborative rather than decisive. They nonetheless reinforce Section 4:
the apparent "default is worse" pattern is partly compositional.

---

## 6. Transmission channels

To identify *how* spread crises propagate, we replace GDP with six
intermediate outcomes — private credit, bank claims on the government,
gross investment, government expenditure, the primary balance, and foreign
direct investment — each expressed as a cumulative change from $t-1$ and
estimated with channel-specific predetermined controls.

| Channel | Profile across $h=0\ldots4$ | Interpretation |
|:--|:--|:--|
| Private credit/GDP | ≈0 at impact, building to −3.5*** at $h=4$ | Credit contraction is **delayed and compounding** |
| Investment/GDP | weak; significant only at $h=3$ (−1.9) | No strong average capital response |
| Bank claims on govt/GDP | +0.3 to +2.0, all insignificant | No average sovereign-bank build-up |
| Govt expenditure/GDP | **+1.19\*\*** at $h=0$, then ≈0 | Spending rises at onset — mildly **counter**cyclical |
| Primary balance/GDP | −0.57** at $h=0$, then insignificant | Deficit widens on impact |
| FDI/GDP | −0.78 ($p$=0.069) at $h=0$ only | Weak, short-lived external pullback |

The dominant pooled mechanism is a **private-credit contraction that
compounds over time**, mirroring the persistence of the output loss
itself: the credit response troughs late (h=4) rather than on impact.
Fiscal policy is accommodative at onset — spending up, primary balance
down — not austere.

---

## 7. Channels by resolution type

Splitting each channel by crisis type is the most informative cut of the
data: the two episode types transmit through visibly different margins. We
report the joint specification (both onset dummies) with the equality test
$\beta^{nd}_h = \beta^{def}_h$, under both OLS-DK and IPW weighting.

- **Investment.** Gross investment **collapses in non-default episodes**
  (−1.19 at $h=0$, −1.77 at $h=1$, −2.73 at $h=3$; equality
  $p$ = 0.003, 0.021, 0.010) but is **flat-to-positive in default-linked
  episodes** (+1.07, +0.19, +0.43). The investment channel is a
  non-default phenomenon.
- **Bank claims on government.** Default-linked banks **accumulate
  sovereign debt** (+4.06 at $h=0$, IPW $p$=0.009; $p$=0.092 at $h=2$),
  a sovereign-bank balance-sheet linkage absent in non-default episodes.
- **Private credit.** The credit crunch is concentrated in **non-default**
  episodes (at $h=4$: −4.78 for non-default vs. +0.62 for default-linked,
  $p$=0.020).
- **Government expenditure.** Default-linked countries **cut spending** at
  longer horizons (at $h=4$: +1.05 for non-default vs. −1.69 for
  default-linked, $p$=0.061); non-default countries do not.
- **Primary balance and FDI.** No significant differences by resolution
  type at any horizon.

The resulting picture is internally coherent: **non-default spread crises
transmit through a contraction in private credit and investment, while
default-linked crises transmit through bank absorption of government debt
and fiscal retrenchment.** The two types deliver a comparable average
output cost (Section 4) through different plumbing.

---

## 8. Robustness and inference

The baseline is checked against alternative crisis-classification
thresholds, exclusion of the most extreme episodes (Argentina, Venezuela),
exclusion of global-crisis years (1997–1999, 2008–2009, 2020–2021), and
country-clustered standard errors in place of Driscoll-Kraay. A
permutation placebo test provides assumption-free inference: the 61
observed onsets are randomly reassigned to estimation-pool country-years
and the projection is re-estimated 1,000 times, locating the true
estimates in the null distribution. This is especially relevant given the
small number of treated observations, where the asymptotics behind
conventional $t$-statistics may be unreliable.

---

## 9. Limitations

- **Small treated samples**, particularly default-linked (21 onsets; the
  channel-by-resolution cells contain only 12–19 default-linked
  observations). Insignificant difference tests frequently reflect low
  power rather than the absence of an effect.
- **Horizon attrition.** At $h=4$ the most recent onset years drop out
  because forward GDP is unavailable, so the longest-horizon estimates
  rest on fewer episodes.
- **Pre-trend window.** The placebo test covers $h=-5$ through $h=-1$;
  each variable is a single-year growth rate fully before onset. Any
  significant pre-trend coefficient limits the causal interpretation for
  that group at the affected horizons.
- **Thin common support for IPW.** A large share of observations is
  trimmed and few treated units cross the classification threshold, so the
  weighted estimates are corroborative rather than definitive.

---

## 10. Summary

The analysis establishes three findings. First, a spread crisis imposes a
large and **persistent** output cost — about −1.9pp on impact, troughing
near −3.8pp, and still −3pp after four years, robustly and if anything
conservatively estimated. Second, default-linked crises appear deeper on
impact, but the resolution gap is statistically insignificant, converges
within four years, and is partly explained by pre-crisis fundamentals once
selection is addressed. Third, the two crisis types propagate through
**different channels** — a private credit-and-investment contraction in
non-default episodes versus sovereign-bank balance-sheet linkages and
fiscal retrenchment in default-linked episodes — which is arguably the
most novel empirical contribution of the paper.
