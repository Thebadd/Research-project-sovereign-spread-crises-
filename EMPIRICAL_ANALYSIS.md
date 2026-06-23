# Empirical Strategy

---

## 1. Research questions

The paper addresses three related questions.

**First question.** What is the output cost of a sovereign spread crisis
in emerging and frontier market economies? Spread crises — episodes in
which the sovereign borrowing spread rises sharply relative to the
risk-free rate — are common in emerging markets but their macroeconomic
consequences are poorly measured. Existing work either focuses on outright
defaults (Reinhart and Rogoff 2009, Cruces and Trebesch 2013) or on
sudden stops (Calvo et al. 2006), leaving the broader category of spread
crises without default largely undocumented.

**Second question.** Does the output cost differ depending on whether the
episode is eventually resolved without default or culminates in debt
restructuring? This is the central comparative question of the paper.
The conjecture is that default-linked crises generate larger and more
persistent output losses, for two reasons: the sovereign ceiling channels
a more severe financing disruption to domestic firms, and external
capital flows stop abruptly when the country enters financial autarky,
depleting the capital stock at a rate that has no counterpart in
non-default crises.

**Third question.** Through which channels does the transmission operate?
Is it primarily a credit contraction (the banking sector channel), a
collapse in investment (the capital accumulation channel), a withdrawal
of external financing (the sudden stop channel), a fiscal consolidation,
or some combination? And do these channels operate symmetrically across
the two crisis types?

The theoretical model in Section 3 makes sharp predictions about the sign,
magnitude, and relative pattern of the output responses across regimes and
channels. The empirical strategy tests those predictions without the
analysis being designed to fit them.

---

## 2. Data and sample

The sample covers **52 emerging and frontier market economies** over
**1994–2025** at annual frequency. Crisis episodes are identified from
a purpose-built database that classifies sovereign spread episodes using
two complementary criteria: a threshold criterion (spreads exceeding
1,000 basis points) and a relative criterion (spreads above the 99th
percentile of the country-specific historical distribution following
Hadzi-Vaskov et al.). For each episode, the first year is designated
the **onset year**; subsequent years within the episode are designated
**continuation years**. The analysis is entirely based on onset years —
continuation years are dropped from every estimation sample. This
prevents the outcome from reflecting the persistence of the crisis itself
rather than its causal impact at inception.

Episodes are classified by resolution outcome: **40 episodes** end
without a sovereign default; **21 episodes** are linked to debt
restructuring. This classification is the key dimension of heterogeneity
in the paper.

---

## 3. Outcome variable and identification approach

### 3.1 The local projection estimator

The empirical backbone of the paper is the **Jordà (2005) local
projection** method. For each forecast horizon $h$, a separate
regression is run with the cumulative change in the outcome variable as
the dependent variable:

$$dy_{i,t+h} = \left(\ln Y_{i,t+h} - \ln Y_{i,t-1}\right) \times 100$$

where $Y_{i,t}$ is real GDP per capita. The outcome at horizon $h$
measures the **cumulative percentage deviation of output from its
pre-crisis level**, with the pre-crisis year $t-1$ as the anchor. At
$h = 0$ this is the impact effect; at $h = 4$ it is the total output
loss four years after onset. The local projection approach has a key
advantage over a VAR: it imposes no cross-horizon restrictions on the
shape of the impulse response function, which is important here because
the two crisis types are expected to generate qualitatively different
dynamic profiles.

### 3.2 Treatment variable and identification

The treatment variable $D_{it}$ equals one in the first year of a crisis
episode and zero in all tranquil years. The identifying assumption is
that, conditional on country and year fixed effects and the set of
observable controls, the onset of a crisis is unrelated to
contemporaneous unobserved shocks that also affect future output growth.

Two features of the design reinforce this assumption. First, the
sample includes **pre-trend placebo horizons** ($h = -2$ and $h = -1$):
the onset dummy is regressed on past output growth. A significant
pre-trend coefficient would indicate that the crisis dummy picks up a
pre-existing trend rather than a causal effect. Second, the sample
restriction to onset years only avoids the mechanical contamination
that would arise if continuation years — in which the crisis is ongoing
by construction — were included in the control group.

### 3.3 Fixed effects and controls

Every regression includes **country fixed effects** $\alpha_i$ and
**year fixed effects** $\lambda_t$. Country fixed effects absorb all
time-invariant country characteristics (institutional quality, geography,
financial development) that could jointly predict both crisis incidence
and output performance. Year fixed effects absorb all global shocks that
hit multiple countries simultaneously — the Global Financial Crisis, the
COVID pandemic, global commodity price cycles — preventing common-shock
contamination of the estimated crisis effect.

The control vector $X_{it}$ includes lagged GDP growth (two lags),
the current account balance, the public debt-to-GDP ratio, inflation,
an IMF program dummy, the VIX index, and the US 10-year Treasury rate.
These controls serve two purposes. First, they account for pre-crisis
macroeconomic conditions that independently predict future output growth
and are also correlated with crisis incidence — omitting them would
produce biased estimates. Second, they partial out the role of global
financial conditions (VIX, US rate) beyond what year fixed effects
capture, which matters because some crisis episodes are partly triggered
by external tightening.

Standard errors follow Driscoll and Kraay (1998), with bandwidth
$\max(1, h+1)$, correcting simultaneously for cross-sectional
dependence (multiple countries hit by a common shock) and serial
correlation in the residuals (which mechanical increases with $h$
because the overlapping forecast windows introduce MA structure).

---

## 4. Act 1 — Average output cost of spread crises

### Specification

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta_h \cdot D_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}, \quad h \in \{-2,-1,0,1,2,3,4\}$$

### What is estimated

The coefficient $\hat\beta_h$ is the average cumulative output loss
$h$ years after a spread crisis onset, relative to tranquil periods
in the same country and year. The sequence $\{\hat\beta_h\}_{h=0}^{4}$
traces the impulse response function of output to a crisis onset.

### What is tested

This is the baseline characterisation. There is no theoretical prior
on the exact magnitude, but the model predicts that $\hat\beta_h < 0$
for all $h \geq 0$ and that the loss persists beyond the spread episode
itself — due to the capital accumulation channel (the capital stock
depresses output even after spreads normalise).

### How to read the results

A coefficient of $\hat\beta_2 = -4$ means that, on average, a country
that experiences a spread crisis onset has real GDP per capita four
percentage points lower two years later than a comparable country in a
tranquil year, after accounting for country and year fixed effects and
observable macro fundamentals. The impulse response function plots this
estimate with 90% and 95% confidence bands for each horizon.

---

## 5. Act 2 — The central comparison: non-default vs. default-linked crises

This is the core empirical contribution of the paper.

### Specification A — Separate regressions

For each group and each horizon $h$:

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{nd}_h \cdot D^{nd}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{def}_h \cdot D^{def}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

### Specification B — Joint regression with equality test

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{nd}_h \cdot D^{nd}_{it} + \beta^{def}_h \cdot D^{def}_{it} + \gamma' X_{it} + \varepsilon_{i,t+h}$$

followed by a test of $H_0: \beta^{nd}_h = \beta^{def}_h$ at each $h$.

### What is estimated

The separate regressions produce two impulse response functions — one
for each crisis type. The joint regression produces both simultaneously
and tests whether the gap between them is statistically significant.

### What is tested

The model makes three predictions about the comparative dynamics.

**Proposition 1** (ranking and divergence): the output loss in the
default regime should be weakly larger than in the non-default regime
at every horizon ($\beta^{def}_h \leq \beta^{nd}_h$), with the gap
widening over time. The divergence is the signature of the
capital-depletion mechanism — external investment stops when the
country enters financial autarky, and this compounds geometrically
with each year of exclusion.

**Proposition 3** (persistence): the persistence should be governed
by different parameters in each regime. In the non-default regime,
persistence comes from the capital stock's slow recovery (governed by
the depreciation rate and the investment sensitivity to the lending
rate). In the default regime, persistence comes from market exclusion
(governed by the re-entry probability $\mu$). These different sources
of inertia imply different shapes for the two IRFs, not merely
different levels.

### Why the joint regression matters

The separate regressions show the level of output losses for each
group but cannot address whether the *gap* is statistically
significant. The joint regression provides the formal test. A finding
that the gap is significant and widens from $h = 1$ onward is the
key empirical result of the paper.

### A note on selection

Countries that eventually default are not randomly drawn from the set
of countries that experience spread crises. Default-linked episodes
tend to start with higher debt, larger current account deficits, and
lower output growth. This pre-crisis imbalance means that part of the
larger output loss in default-linked episodes may reflect worse
initial conditions rather than a causal effect of the default itself.
The balance table (Step 0) quantifies this selection. The IPW
correction (Section 7) addresses it.

---

## 6. Act 3 — Role of IMF programs

### Specification

$$dy_{i,t+h} = \alpha_i + \lambda_t + \beta^{imf}_h \cdot D^{imf}_{it} + \beta^{no\text{-}imf}_h \cdot D^{no\text{-}imf}_{it} + \gamma' X^{-imf}_{it} + \varepsilon_{i,t+h}$$

Test: $H_0: \beta^{imf}_h = \beta^{no\text{-}imf}_h$.

### What is tested

Among spread crisis episodes, does access to IMF financing at onset
attenuate the output loss? The IMF programs directly affect the fiscal
channel (by providing financing that avoids abrupt austerity) and
potentially the sovereign ceiling channel (by providing an implicit
backstop that limits spread contagion to the banking sector).

The IMF dummy is excluded from the control vector in this regression
because it is itself the treatment modifier — including it as a control
in other regressions (where it does appear) is appropriate because its
effect there is absorbed into the baseline, not the treatment.

### Selection concern and IPW correction

Countries that obtain IMF programs at crisis onset are systematically
different from those that do not: they tend to have higher debt ratios,
worse current account positions, and lower growth. Comparing the two
groups directly would confound the IMF program effect with these
pre-existing differences. A first-stage probit of
$\Pr(\text{IMF} = 1 \mid \text{crisis onset}, X)$ generates propensity
scores, which are used to reweight non-IMF episodes to match IMF
episodes on observables. The IPW-corrected LP gives the counterfactual
effect of the IMF program purged of selection.

---

## 7. Act 4 — Transmission channels

### Why channel regressions?

The output LP tells you that crises are costly. It does not tell you
*why*. The channel regressions decompose the transmission mechanism by
replacing the GDP outcome with six intermediate variables, each
corresponding to a specific theoretical channel. This is a direct test
of the model's structural mechanisms.

### Common structure

For each channel variable $Z$ and each horizon $h$:

$$\Delta^h Z_{it} = \alpha_i + \lambda_t + \beta^Z_h \cdot D_{it} + \gamma^Z{}' X^Z_{it} + \varepsilon_{i,t+h}$$

where $\Delta^h Z_{it} = Z_{i,t+h} - Z_{i,t-1}$ is the cumulative
change in the channel variable from one year before the crisis to $h$
years after. The controls are channel-specific rather than uniform,
because each channel is driven by a different set of confounders.

### Channel 1 — Private credit / GDP

**What is estimated.** The cumulative change in private credit relative
to GDP following a crisis onset.

**Theoretical link.** This tests equation (2) of the theoretical model:
a spread widening reduces bank net worth (equation 1), which forces a
proportional contraction in lending through the leverage constraint.
The credit response is the empirical counterpart of
$\hat\ell_h = (\lambda N_{ss}/\ell_{ss})\hat n_h$.

**Economic interpretation.** A negative $\hat\beta^{cr}_h$ means the
banking sector actively contracts lending following a crisis — a supply
contraction, not just a demand response. Combined with the output LP,
it reveals the credit-to-output ratio, which is the key quantity in
Proposition 2.

### Channel 2 — Bank claims on sovereign / GDP

**What is estimated.** The cumulative change in the share of bank assets
held in government bonds following a crisis onset.

**Theoretical link.** This tests the doom loop mechanism. In the model,
banks hold a fixed portfolio ($b^B_{ss}$ constant). If the empirical
estimate is significantly positive, it indicates that banks are
*increasing* their sovereign bond exposure during the crisis —
the gambling-for-resurrection or carry-trade motive that the baseline
model abstracts from. A null result supports the symmetric bank
behaviour assumption.

**Economic interpretation.** An increase in bank claims on the sovereign
during a crisis is consistent with portfolio substitution: banks reduce
lending to firms and replace it with sovereign bonds, amplifying the
credit contraction and deepening the doom loop.

### Channel 3 — Gross investment / GDP

**What is estimated.** The cumulative change in gross fixed capital
formation relative to GDP.

**Theoretical link.** This directly tests the capital accumulation
channel, the most persistent transmission mechanism in the model.
Equation (4) predicts that a spread increase depresses investment
at $h=0$, which lowers the capital stock at $h=1$, which depresses
output at $h=1$ even after the spread has normalised.

**Economic interpretation.** The investment response is the bridge
between the contemporaneous lending rate shock and the persistence of
the output loss. A coefficient that remains negative at $h = 3$ or
$h = 4$ — after the spread itself has likely returned to normal —
is evidence that the capital channel is quantitatively important.

In the default regime specifically, the external investment share $\chi$
(the fraction of investment financed by foreign capital flows) matters:
when the country enters autarky, external financing disappears and
investment falls by more than the credit contraction alone would imply.
This is the mechanism in equation (6) of the model.

### Channel 4 — Government expenditure / GDP

**What is estimated.** The cumulative change in government spending
relative to GDP.

**Economic interpretation.** A negative coefficient means fiscal
austerity — governments cut spending when the crisis hits, either
because market access deteriorates (forcing fiscal tightening) or
because IMF conditionality requires it. Procyclical fiscal consolidation
amplifies the output contraction beyond what the financial channel alone
produces.

### Channel 5 — Primary balance / GDP

**What is estimated.** The cumulative change in the primary fiscal
balance following a crisis onset.

**Economic interpretation.** An improvement in the primary balance
(positive coefficient) during a crisis means the government is
tightening fiscal policy at precisely the moment when automatic
stabilisers would normally imply a deterioration. This is classic
procyclical fiscal policy, which deepens the recession. Combined with
Channel 4, it identifies whether fiscal austerity is a quantitatively
significant transmission mechanism beyond the financial channels.

### Channel 6 — FDI / GDP

**What is estimated.** The cumulative change in foreign direct
investment relative to GDP.

**Theoretical link.** This directly tests the $\chi$ parameter of the
model — the external investment share. Equation (6) assumes that a
fraction $\chi$ of steady-state replacement investment requires access
to international capital markets. A significant decline in FDI following
a crisis onset provides empirical support for $\chi > 0$ and identifies
the external financing stop that underlies the capital depletion
mechanism in the default regime.

**Why the US rate appears in this regression only.** The US 10-year
Treasury rate captures the carry trade incentive of foreign investors —
when US rates are high, the opportunity cost of investing in emerging
markets rises and FDI declines. This variable adds genuine identification
content for foreign capital flows but is irrelevant for domestic
real variables (credit, investment, fiscal aggregates) that are
driven by domestic rather than foreign financing conditions.

---

## 8. Act 5 — Channels by resolution type

### What is estimated

For each of the six channel variables, a joint regression is run with
both the non-default and the default-linked onset dummies simultaneously:

$$\Delta^h Z_{it} = \alpha_i + \lambda_t + \beta^{nd,Z}_h \cdot D^{nd}_{it} + \beta^{def,Z}_h \cdot D^{def}_{it} + \gamma^Z{}' X^Z_{it} + \varepsilon_{i,t+h}$$

with an equality test $H_0: \beta^{nd,Z}_h = \beta^{def,Z}_h$ at each
horizon.

### What is tested

This is the cleanest test of the model's structural predictions.

**Credit-to-output ratio (Proposition 2).** The model predicts that the
ratio of the credit contraction to the output contraction is *larger*
in non-default episodes. In default episodes, output carries an
additional autonomous depletion term from the capital stock — absent
from the credit equation — which makes output fall more relative to
credit. Empirically, $|\hat\beta^{nd,cr}_h / \hat\beta^{nd,y}_h| >
|\hat\beta^{def,cr}_h / \hat\beta^{def,y}_h|$ would confirm this.

**Investment channel asymmetry.** Investment should fall more in
default-linked episodes through the external financing stop, even after
controlling for the credit contraction. This is $\hat\beta^{def,inv}_h
< \hat\beta^{nd,inv}_h$ at horizons $h \geq 1$.

**Sovereign-bank nexus.** If banks respond to default risk by
concentrating sovereign bond holdings (doom loop), the claims-on-
government response should be more positive in default-linked episodes.
A null result across both groups supports the symmetric bank behaviour
assumption of the model.

### IPW correction for selection into default

Even among crisis episodes, the subset that ends in default is not a
random draw. Default-linked episodes have systematically higher pre-crisis
debt and worse external positions. A first-stage probit of
$\Pr(\text{default-linked} = 1 \mid \text{crisis onset}, \text{debt},
\text{CA})$ generates propensity scores that reweight non-default
episodes to match default-linked episodes on observables. The IPW
correction is applied to every channel regression, producing two
estimates — unweighted OLS and IPW — for each channel at each horizon.
Stability of the estimates across the two methods signals that the
channel asymmetries are genuine effects of the resolution outcome, not
artefacts of selection.

---

## 9. Act 6 — Mechanistic tests

### Test 1 — Is the credit contraction supply-driven or demand-driven?

**The question.** A crisis reduces credit. But is this because banks are
*willing* to lend less (supply contraction) or because firms want to
borrow less (demand contraction)? The distinction matters: a supply
contraction implies active balance-sheet management by banks; a demand
contraction may require a different policy response.

**The approach.** Add the lagged stock of bank sovereign bond holdings
to the credit LP. If the crisis-onset coefficient on credit shrinks
significantly when this variable is included, it means that part of the
credit decline is explained by banks substituting sovereign bonds for
firm loans — a portfolio reallocation that contracts credit supply
independently of demand. If the coefficient is unchanged, the credit
decline is driven by demand or by funding-cost mechanisms that are
orthogonal to sovereign bond accumulation.

**Interpretation.** The mediation through sovereign bond accumulation
is the empirical analogue of the doom loop: the balance-sheet hit from
rising spreads depletes net worth, forcing deleveraging, and banks
offset the loss of yield by loading up on high-carry sovereign bonds
rather than rolling over loans to firms. If this mechanism is
quantitatively important, the model's abstraction from endogenous
portfolio choice — the fixed $b^B_{ss}$ assumption — understates the
amplification.

### Test 2 — Does credit mediate the investment contraction?

**The question.** Reduced investment could occur because the lending
rate rose (the direct channel, through the firm's zero-profit condition
for capital), or because credit availability fell (the quantity channel,
through the leverage constraint). These are distinct mechanisms.

**The approach.** Run the investment LP twice: once without controlling
for the credit stock, and once including the lagged credit-to-GDP ratio.
The coefficient on the crisis dummy in the first regression captures the
**total effect** of the crisis on investment. In the second regression,
it captures only the **direct effect** net of the credit quantity
contraction. The difference — expressed as a share of the total effect —
is the mediation share attributable to the credit channel.

**Interpretation.** A high mediation share (close to one) means that
if you prevented the credit contraction from happening, you would also
prevent most of the investment decline. A low mediation share means
that investment falls even when credit is held constant, which is
consistent with the direct lending-rate channel: firms cut investment
because the cost of new capital rises, not because they are rationed
out of the credit market.

This test has direct implications for policy: a credit-mediated
investment contraction is better addressed by bank recapitalisation or
credit guarantees; a directly rate-mediated contraction requires spread
compression, which is the role of IMF programs or central bank
interventions.

---

## 10. Identification challenges and corrections

### Selection on observables — IPW

The central concern is that default-linked episodes are not randomly
assigned. Countries that default are systematically weaker along
observable dimensions (higher debt, worse current account, lower growth)
before the crisis begins. This means a naive comparison of output losses
between the two groups confounds the causal effect of default with
pre-existing fragility. The inverse probability weighting correction
reweights the control group to match the treated group on observables,
recovering the average treatment effect on the treated under the
assumption of selection on observables.

### Selection on unobservables — placebo test

IPW corrects for selection on measured characteristics. Selection on
unmeasured characteristics — a country is more likely to default because
of a political crisis, institutional fragility, or other unobservables
that also affect output growth — cannot be addressed by reweighting.
The permutation placebo test provides a finite-sample, assumption-free
$p$-value: the 61 real crisis onsets are randomly reassigned 1,000 times,
and the LP is re-estimated each time. The true estimate's position in
the resulting null distribution measures how unusual the result is
under the null of no causal effect, regardless of the distributional
assumptions underlying the Driscoll-Kraay $t$-statistic.

### Robustness across specifications

Ten alternative specifications are tested: two alternative crisis
thresholds (testing whether the classification drives the result),
four country and time exclusions (Argentina, Venezuela, GFC years,
Asian crisis years), and an alternative standard error approach
(country clustering). Stability of the main coefficients across
these specifications — particularly the gap between $\hat\beta^{nd}_h$
and $\hat\beta^{def}_h$ — is the key test of robustness.

---

## 11. Connection between theory and empirics

The theoretical model in Section 3 is calibrated independently of the
LP estimates — its parameters come from the prior EM-DSGE literature,
from balance-of-payments data (for the external investment share $\chi$),
and from tranquil-period panel means (for the steady-state ratios). The
model is then simulated to generate impulse response functions for output,
credit, and investment under each crisis regime. These model-implied IRFs
are overlaid on the empirical IRFs from Acts 1 and 2.

This comparison is out-of-sample in the following sense: the model
parameters are not chosen to match the LP estimates. The comparison
therefore tests whether the transmission mechanisms embedded in the model
— the balance-sheet channel, the working-capital channel, the capital
accumulation channel, the external financing stop — are, taken together,
quantitatively consistent with what the local projections measure
reduced-form.

Three specific predictions connect the model to the empirical results.

**Prediction 1 (output divergence)** corresponds to the significance
and widening of the gap between $\hat\beta^{def}_h$ and
$\hat\beta^{nd}_h$ in the joint regression of Act 2. The model
explains this divergence through the autonomous capital depletion term
in equation (7), which has no counterpart in the non-default path.

**Prediction 2 (credit-to-output ratio)** corresponds to the relative
magnitude of the credit response (Channel 1) versus the output response
in Act 5. The model explains the asymmetry through the difference in
loadings: credit loads only on net worth, while output loads on both
net worth and the capital stock.

**Prediction 3 (persistence)** corresponds to the shape of the IRFs at
long horizons ($h = 3, 4$). The model predicts that output remains below
trend even after spreads normalise, through the capital stock channel
in the non-default regime and through market exclusion in the default
regime. The persistence of the empirical IRFs is the direct test of
this prediction.

If the model-implied IRFs match the empirical IRFs quantitatively — not
just qualitatively — it provides structural content to the causal
interpretation of the local projection estimates: the output losses are
not mere correlations but the expected consequences of the structural
mechanisms the model isolates.
