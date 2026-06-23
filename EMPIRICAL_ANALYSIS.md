# Empirical Strategy

## 1. Setting and identification

The empirical analysis rests on a panel of 52 emerging and frontier market
economies over 1994–2025. The unit of observation is the country-year, and
the sample is restricted to onset years — the first year of each crisis
episode — and tranquil years. Continuation years are excluded throughout,
so the identifying variation comes from the transition into crisis, not
from the crisis persisting.

The estimator is the Jordà (2005) local projection. For each horizon
$h \in \{0, 1, 2, 3, 4\}$, we estimate a separate regression of the
cumulative change in an outcome variable on a crisis onset indicator,
country and year fixed effects, and a vector of pre-determined macro
controls:

$$\Delta^h y_{it} = \alpha_i + \lambda_t + \beta_h \cdot D_{it} + \gamma' X_{it} + \varepsilon_{it+h}$$

The dependent variable $\Delta^h y_{it}$ is the cumulative percentage
deviation of real GDP per capita from its pre-crisis level — anchored at
$t-1$ so that $h = 0$ captures the impact effect and $h = 4$ the
four-year cumulative loss. The sequence $\{\hat\beta_h\}$ traces the
impulse response of the outcome to a crisis onset.

The control vector $X_{it}$ includes two lags of GDP growth, the
current account balance, the public debt-to-GDP ratio, inflation, and
an IMF program dummy. Lagged growth controls for pre-crisis output
momentum that would independently predict future trajectories; debt and
the current account capture the fiscal and external vulnerability of
the economy at onset, which jointly predict both crisis incidence and
subsequent output performance; inflation proxies for underlying macro
instability; and the IMF program dummy absorbs the differential
trajectory of countries that obtain official financing during the
episode. Global financial conditions — the VIX and the US long-term
rate — are not included in the specification because year fixed effects
already absorb all shocks common to all countries in a given year,
rendering any pure time-series variable redundant.
Standard errors follow Driscoll and Kraay (1998), correcting for
cross-sectional dependence and for the serial correlation that
mechanically increases with the forecast horizon.

The identifying assumption is that, conditional on fixed effects and
observable macro fundamentals, the onset of a spread crisis is
uncorrelated with contemporaneous unobserved shocks to future output. We
test this directly by estimating the regression at placebo horizons
$h = -2$ and $h = -1$: if pre-crisis output trends predict crisis
incidence, the assumption fails. Non-zero pre-trend coefficients would
indicate anticipation effects or trending differences between treated and
control countries that bias the main estimates.

---

## 2. The central comparison: crisis resolution and output losses

The paper's core empirical question is whether the output cost of a spread
crisis depends on its resolution — specifically, whether episodes resolved
without default inflict smaller and less persistent losses than those
culminating in debt restructuring. To answer this, we estimate a joint
local projection in which the onset dummies for the two crisis types enter
simultaneously:

$$\Delta^h y_{it} = \alpha_i + \lambda_t + \beta^{nd}_h \cdot D^{nd}_{it} + \beta^{def}_h \cdot D^{def}_{it} + \gamma' X_{it} + \varepsilon_{it+h}$$

The coefficient $\hat\beta^{nd}_h$ captures the average output loss $h$
years after the onset of a spread crisis resolved without default;
$\hat\beta^{def}_h$ captures the same for default-linked episodes. The
joint specification is essential: it tests whether the gap
$\hat\beta^{def}_h - \hat\beta^{nd}_h$ is statistically significant at
each horizon, and whether it widens as $h$ increases.

The theoretical model predicts that this gap is strictly negative for
$h \geq 1$ and non-decreasing in $h$ (Proposition 1). The intuition is
that in default-linked episodes, access to external financing collapses
when the sovereign enters financial autarky, depressing investment and
eroding the capital stock at a rate that compounds with each year of
exclusion. This capital-depletion mechanism has no counterpart in
non-default crises, where the transmission operates exclusively through
lending rates and bank balance sheets. A widening, statistically
significant gap in the empirical impulse response functions is the
direct test of this mechanism.

A legitimate concern is that the two groups are not exchangeable ex ante:
countries that eventually default tend to enter the crisis with higher
debt ratios, worse current account positions, and weaker output momentum.
The balance table in Section 4.1 quantifies these differences. We address
them through inverse probability weighting, constructing propensity scores
from a first-stage probit of crisis type on pre-crisis fundamentals and
reweighting the non-default group to match the default-linked group on
observables. Stability of the gap estimates across OLS and IPW
specifications is evidence that selection on observable characteristics
does not drive the result.

---

## 3. Transmission channels

Establishing that output losses differ across crisis types is a necessary
but not sufficient contribution: the paper also aims to identify the
structural mechanisms behind the difference. We do this by replacing GDP
with a set of intermediate outcome variables — private credit, bank
sovereign bond holdings, gross investment, government expenditure, the
primary fiscal balance, and foreign direct investment — and running the
same local projection for each.

The credit response tests the balance-sheet channel. When sovereign
spreads widen, banks holding government bonds suffer mark-to-market
losses that erode their net worth. Through the leverage constraint,
this forces a contraction in lending that is proportional to the equity
loss. A negative $\hat\beta^{cr}_h$ confirms that this channel is
empirically active; the ratio of the credit response to the output
response across the two crisis types then tests Proposition 2 of the
model, which predicts that this ratio is larger in non-default episodes,
where output does not carry the additional capital depletion term.

The investment response is the empirical counterpart of the capital
accumulation channel. A persistent decline in investment — one that
remains negative at horizons $h = 3$ and $h = 4$, well after the spread
itself has likely normalised — is evidence that the crisis leaves a
lasting imprint on the capital stock, not merely a temporary
working-capital disruption. In default-linked episodes specifically,
a larger investment decline at medium horizons relative to non-default
episodes would point to the external financing stop mechanism embedded in
the model: when the sovereign loses market access, the share of investment
normally financed by external flows disappears, generating capital
depletion beyond what lending rate movements alone would predict.

Bank sovereign bond holdings test whether banks actively accumulate
government debt during crises — the doom loop. The theoretical model
abstracts from this margin, assuming fixed sovereign bond portfolios. A
significant increase in bank claims on the sovereign following a crisis
onset would indicate that this abstraction is costly and that the
balance-sheet channel is amplified beyond the model's prediction.

The fiscal responses — government expenditure and the primary balance —
capture the role of procyclical fiscal policy. If governments cut spending
or improve their fiscal balance during a crisis, they amplify the output
contraction through a demand channel that operates independently of the
financial transmission. The foreign direct investment response directly
measures the external financing stop, providing empirical content for the
parameter $\chi$ in the model, which governs what fraction of investment
requires access to international capital markets.

For each channel, we estimate the regression separately by crisis type and
test whether the impulse response functions differ, applying the same IPW
correction to account for pre-crisis differences between the two groups.

---

## 4. Mechanistic tests

Two additional regressions sharpen the identification of the credit
channel. The first asks whether the credit contraction is supply-driven or
demand-driven. We add pre-crisis bank sovereign bond holdings to the credit
regression. If the onset coefficient shrinks materially, it means that
banks are substituting sovereign bonds for firm loans — a supply
contraction driven by portfolio reallocation rather than by a fall in
credit demand. The second test asks whether the credit contraction
mediates the investment decline. We estimate the investment regression with
and without the credit stock as a control. The fraction of the total
investment response absorbed by the credit control — the mediation share
— identifies what portion of the investment decline passes through the
quantity of bank lending, as opposed to the direct effect of higher
borrowing costs on the return to new capital.

---

## 5. Robustness and placebo tests

The main results are verified across a battery of robustness checks. We
replicate the central regression under alternative crisis classification
thresholds, dropping outlier countries with the most extreme episodes
(Argentina, Venezuela), excluding global crisis years (2008–2009,
1997–1999, 2020–2021), and replacing Driscoll-Kraay standard errors with
country-clustered standard errors. Stability across these specifications
addresses the concern that the baseline estimate is sensitive to the
definition of the treatment, to a handful of influential observations, or
to the particular SE correction.

The permutation placebo test provides assumption-free inference. We
randomly reassign the 61 observed crisis onsets to country-year
observations drawn from the estimation pool and re-estimate the local
projection 1,000 times. The empirical distribution of the placebo
estimates under the null of no treatment effect is then compared to the
true estimates. This test is particularly relevant given the relatively
small number of treated observations, where asymptotic approximations
underlying standard $t$-statistics may be unreliable.

---

## 6. Connecting theory and empirics

The model is calibrated independently of the local projection estimates.
Structural parameters are drawn from the quantitative sovereign default
literature, from balance-of-payments data over tranquil periods, and from
panel regressions that exploit only the pre-crisis variation in the data.
The model is then simulated to produce impulse response functions for
output, credit, and investment under each crisis regime. These
model-implied responses are overlaid on the empirical estimates.

The comparison is out-of-sample by construction: the model parameters are
not chosen to match the IRFs, so agreement between the two is not
mechanical. If the model-implied and empirical impulse responses coincide
in sign, magnitude, and dynamic profile, it provides structural
interpretation for the reduced-form estimates — the output losses are not
mere statistical patterns but the predictable consequence of the
mechanisms the model isolates. Divergence, by contrast, would point to
margins the model does not capture and would guide future theoretical work.
