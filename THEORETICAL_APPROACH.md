# 3. Theoretical Framework

To interpret the empirical findings and discipline their quantitative
magnitudes, we develop a two-block model of sovereign risk and real
transmission. The first block is a nonlinear sovereign default model in the
tradition of Arellano (2008), which endogenously generates the sovereign
spread and classifies crisis episodes as either default or non-default based
on the government's optimal default decision. The second block is a
log-linearized banking and production system, adapted from Bocola (2016),
which maps the sovereign spread path into output and credit responses through
three distinct channels: a balance-sheet channel operating through bank net
worth, a working-capital channel operating through the cost of intermediate
inputs, and a capital-accumulation channel operating through firms' investment
decisions. The two blocks are deliberately kept separate: the default block
generates the spread paths that feed into the transmission block, but the
transmission block does not feed back into the default decision. This
recursive structure keeps the model analytically tractable and ensures that
each block can be separately identified.

Crucially, the model is constructed without reference to the local projection
estimates it will be compared against. All parameters are calibrated either
from the prior literature or from moments of the data that are orthogonal to
the impulse response functions — specifically, long-run averages of sovereign
spreads, default frequencies, debt-to-GDP ratios, and bank balance sheet
ratios computed over tranquil periods. The comparison of model-implied impulse
responses to the LP estimates in Section 4 then constitutes a genuine
out-of-sample test: the model makes predictions that could be falsified.

---

## 3.1 Economic Environment

Time is discrete and indexed by $t$. The economy is populated by a
representative household, a government, a unit measure of competitive banks,
and a representative firm. The economy is small and open, facing a
world risk-free rate $R^*$. The household supplies labor inelastically,
consumes, and saves through bank deposits. The government issues
one-period discount bonds to foreign investors to smooth consumption, and
retains the option to default on its obligations. Banks intermediate between
the sovereign bond market and the domestic firm sector, holding both
government bonds and issuing loans to firms. Firms combine capital and labor
using a Cobb-Douglas technology and borrow working capital from banks to
finance their wage bill.

The economy is subject to two types of aggregate shock. The first is a
persistent endowment shock $y_t$, which follows a log-AR(1) process
estimated from the panel data (Section 3.9). This shock drives the
government's incentive to default: a sufficiently bad income draw makes
autarky preferable to debt repayment. The second is the endogenous sovereign
spread $s_t$, which rises in response to deteriorating fundamentals and
transmits to the real economy through the banking and production blocks.

---

## 3.2 The Sovereign's Problem

The government is benevolent but discounts the future at rate $\beta' < \beta$,
where $\beta$ is the household's discount factor. This impatience — standard
in the sovereign default literature (Arellano 2008; Aguiar and Gopinath 2006)
— generates a realistic default frequency without requiring implausibly large
output costs. The government issues one-period discount bonds $B_t$ with face
value normalized to one. A bond sold at price $q_t$ raises $q_t \cdot B_t$
units of resources today and promises to repay $B_t$ next period.

**Preferences.** The government maximizes the expected discounted utility of
the representative household:

$$V_t = \mathbb{E}_t \sum_{j=0}^{\infty} (\beta')^j u(c_{t+j})$$

where $u(c) = c^{1-\sigma}/(1-\sigma)$ with $\sigma = 2$ is a constant
relative risk aversion utility function. The CRRA specification is standard
in the quantitative sovereign default literature and implies that the
government has a stronger incentive to default when income is low, since the
marginal utility of consumption is high at precisely those times.

**Repayment.** If the government repays, it chooses next-period debt $B_{t+1}$
to maximize continuation value, subject to the budget constraint:

$$c_t = y_t - B_t + q(B_{t+1}, y_t) \cdot B_{t+1}$$

where $q(B_{t+1}, y_t)$ is the equilibrium bond price, which is endogenous
and depends on the borrowing level and current income because both affect the
probability of future default. The value of repayment is:

$$V^R(B_t, y_t) = \max_{B_{t+1}} \Bigl\lbrace  u\left(y_t - B_t + q(B_{t+1}, y_t) \cdot B_{t+1}\right) + \beta' \,\mathbb{E}_{y'|y}\left[V(B_{t+1}, y')\right] \Bigr\rbrace$$

**Default.** If the government defaults, it receives income $y_t^{def}$
in the current period and is excluded from international capital markets.
Following Arellano (2008), the default output cost is asymmetric: the
autarky income is $y_t^{def} = \min(y_t,\, \bar{h} \cdot \bar{y})$, where
$\bar{y} = \mathbb{E}[y]$ is mean income and $\bar{h} = 0.965$. This
formulation captures the empirical regularity that defaults occur precisely
when output is already low, so the output cost of exclusion is small for
countries in severe distress but non-trivial for those with moderate income.
The excluded sovereign re-enters international markets with probability $\mu$
each period, at which point it restarts with zero debt. The value of
default is:

$$V^D(y_t) = u(y_t^{def}) + \beta' \Bigl\lbrace  \mu \cdot \mathbb{E}_{y'|y}\left[V(0, y')\right] + (1-\mu) \cdot \mathbb{E}_{y'|y}\left[V^D(y')\right] \Bigr\rbrace$$

Note that $V^D$ does not depend on $B_t$: once the government defaults, the
stock of outstanding debt is irrelevant because it is being repudiated.

**Default decision.** The government defaults whenever the value of default
exceeds the value of repayment. The optimal policy is:

$$d(B_t, y_t) = \mathbf{1}\Bigl\lbrace  V^D(y_t) > V^R(B_t, y_t) \Bigr\rbrace$$

The default set $\mathcal{D}(B_t)$ is the set of income realizations for
which default is optimal given debt level $B_t$:

$$\mathcal{D}(B_t) = \bigl\lbrace y_t \;:\; V^D(y_t) > V^R(B_t, y_t) \bigr\rbrace$$

Because $V^R$ is increasing in $y_t$ (higher income makes repayment easier)
while $V^D$ is non-decreasing in $y_t$ (default always relieves the debt
burden), the default set takes the form of a lower interval
$\mathcal{D}(B_t) = (-\infty, \hat{y}(B_t)]$ for some threshold $\hat{y}$.
The government defaults on bad income draws and repays on good ones.

The value function satisfies $V(B_t, y_t) = \max( V^R(B_t, y_t),\, V^D(y_t) )$.

---

## 3.3 Bond Pricing and the Sovereign Spread

Foreign investors are risk-neutral and have access to the world risk-free
rate $R^*$. Competition drives bond prices to the zero-profit level. A bond
that promises to pay one unit next period delivers 1 unit if the sovereign
repays and $\theta \in (0,1)$ units if it defaults, where $\theta$ is the
recovery rate (one minus the haircut). The bond price $q(B_{t+1}, y_t)$
satisfies the zero-profit condition:

$$q(B_{t+1}, y_t) = \frac{1 - \delta(B_{t+1}, y_t)\cdot(1 - \theta)}{1 + R^*}$$

where $\delta(B_{t+1}, y_t)$ is the next-period default probability:

$$\delta(B_{t+1}, y_t) = \mathbb{E}_{y'|y}\left[ d(B_{t+1}, y') \right]$$

This is evaluated at the borrowing level $B_{t+1}$ and current income $y_t$,
which determines the distribution of $y'$ via the transition matrix. This pricing equation
embeds three features of the data. First, higher debt $B_{t+1}$ raises
default risk, lowering the bond price and raising the cost of borrowing —
the government faces an upward-sloping supply curve for credit. Second, lower
current income $y_t$ shifts the transition distribution toward low future
income, also raising default risk and the spread. Third, partial recovery
$\theta > 0$ limits the spread to a finite level even when the default
probability approaches one, consistent with the observation that sovereign
bonds trade at positive prices even during debt crises.

The sovereign spread is the excess yield on government bonds over the
risk-free rate:

$$s_t = \frac{1}{q(B_{t+1}, y_t)} - 1 - R^*$$

When $\delta = 0$ (no default risk), $q = 1/(1+R^*)$ and $s = 0$. When
default risk is positive, the price falls below $1/(1+R^*)$ and the spread
becomes positive. At the calibrated parameters, the model generates a mean
spread of 200-400 basis points and a standard deviation consistent with the
EMBIG data, without targeting these moments directly in the transmission block.

**Equilibrium.** An equilibrium consists of value functions $V^R$, $V^D$,
$V$, a bond price schedule $q(B', y)$, default policy $d(B, y)$, and debt
policy $B'(B, y)$ such that: (i) the default policy is optimal given the
bond price; (ii) the debt policy is optimal given the bond price and the
default option; (iii) the bond price reflects competitive foreign investors
pricing in the endogenous default probability. We solve for the equilibrium
by value function iteration over a discretized state space (Section 3.9).

---

## 3.4 Banking Sector and Financial Intermediation

Banks are the key link between the sovereign bond market and the domestic
real economy. Each bank holds two assets: sovereign bonds $b^B_t$ and loans
to firms $\ell_t$. These are funded by equity (net worth) $N_t$ and
short-term deposits $d_t$. The bank balance sheet identity is:

$$b^B_t + \ell_t = N_t + d_t$$

**Net worth.** Bank net worth evolves according to realized returns on
assets minus the cost of liabilities. Sovereign bonds pay off $q_{t+1}$
tomorrow per unit of face value held today. When the sovereign spread widens,
the market value of sovereign bonds held by banks falls, inflicting a
capital loss on bank equity. Specifically, a bank that entered the period
with sovereign bond holdings $b^B_{t-1}$ at purchase price $q_{t-1}$ suffers
a mark-to-market loss of $b^B_{t-1} \cdot \Delta q_t$ when the bond price
falls by $\Delta q_t = q_t - q_{t-1}$. In log-linear terms, the deviation of
net worth from its steady-state value $N_{ss}$ satisfies:

$$\hat{n}_t = \Phi_N \hat{n}_{t-1} - B_{sens} \hat{s}_t - \Phi^O \gamma \hat{s}_{t-1} \quad\text{(LL.5a)}$$

where $\hat{n}_t = (N_t - N_{ss})/N_{ss}$ denotes the proportional deviation
of net worth from steady state, and $\hat{s}_t = s_t - s_{ss}$ is the
spread deviation measured in decimal units. The coefficient
$B_{sens} = (b^B_{ss}/N_{ss})/(1 + s_{ss})^2$ captures the balance-sheet
sensitivity of net worth to spread changes: it is proportional to the
ratio of sovereign bond holdings to equity, $b^B_{ss}/N_{ss}$, which
measures how leveraged the bank is with respect to sovereign risk. A bank
with large sovereign bond holdings relative to equity suffers a large net
worth loss for any given spread increase. The parameter $\Phi_N \in (0,1)$
governs the persistence of net worth deviations — it reflects retained
earnings, recapitalization frictions, and the gradual adjustment of bank
balance sheets. The term $\Phi^O \gamma \hat{s}_{t-1}$ captures delayed
pass-through from wholesale funding costs, discussed further below.

**Leverage constraint.** A binding leverage constraint limits the volume of
loans a bank can extend to a multiple $\lambda$ of its net worth:

$$\ell_t = \lambda \cdot N_t \quad\text{(IR.4)}$$

This constraint — standard in the bank-capital literature (Gertler and
Karadi 2011; Bocola 2016) — is the key mechanism linking sovereign risk to
private credit. When a spread widening reduces bank net worth, the
leverage constraint forces a proportional contraction in lending:
$\hat{\ell}_t = \hat{n}_t$. The parameter $\lambda$ is the leverage
multiplier; a higher $\lambda$ amplifies the credit response to any given
net worth shock. We calibrate $\lambda$ directly from the data as the ratio
of total bank assets (sovereign bonds plus loans) to bank equity in tranquil
periods (Section 3.9).

**Lending rate and the sovereign ceiling.** Banks fund firm investment at the
lending rate $R^L_t$. Following the sovereign ceiling doctrine — the empirical
regularity that banks in emerging markets cannot borrow more cheaply than
their own sovereign (Borensztein et al. 2006) — the bank lending rate is tied
to the sovereign spread through a pass-through coefficient $\gamma$:

$$R^L_t = R^* + \gamma \cdot s_t \quad\text{(SC)}$$

where $\gamma \in (0,1]$ measures how much of the sovereign spread is
transmitted to the private lending rate. We set $\gamma = 0.80$ following
Bocola (2016), consistent with the empirical pass-through estimates for
emerging markets. In the steady state, the lending rate $R^L_{ss}$ satisfies
$R^L_{ss} = R^* + \gamma s_{ss}$, exceeding the world risk-free
rate by the steady-state spread premium. An additional balance-sheet
component links the lending rate to bank net worth: when banks are
undercapitalized, they raise lending rates to preserve margins, adding the
term $\Omega \hat{n}_{t}{}$ to the effective financing cost faced by firms, where
$\Omega = \phi \cdot b^B_{ss} N_{ss}^{-1}{}$ is the balance-sheet pass-through
scaled by the sovereign exposure ratio.

---

## 3.5 Firm Production and the Working-Capital Channel

The representative firm produces output $Y_t$ using capital $K_t$ and
labor $L_t$ with a Cobb-Douglas production function:

$$Y_t = A \cdot K_t^\alpha \cdot L_t^{1-\alpha}$$

where $\alpha$ is the capital share and $A$ is total factor productivity
(normalized to one in the steady state). Labor is supplied inelastically by
the household at a competitive wage $w_t = (1-\alpha) Y_t / L_t$.

**Working-capital requirement.** Following Neumeyer and Perri (2005), a
fraction $\xi$ of the wage bill must be financed in advance through bank
loans, before production takes place. This requires the firm to borrow
$\xi \cdot w_t \cdot L_t$ at the bank lending rate $R^L_t$. The total
cost of labor, including financing costs, is:

$$\tilde{w}_t = w_t \cdot \left[1 + \xi \cdot (R^L_t - 1)\right]$$

An increase in the sovereign spread raises $R^L_t$, which raises the
effective cost of labor and reduces the firm's optimal employment and output.
Log-linearizing around the steady state, the elasticity of output with
respect to the lending rate is:

$$\varepsilon_p = \frac{\xi \cdot \frac{1-\alpha}{\alpha} \cdot R^L_{ss}}{1 + \xi(R^L_{ss} - 1)}$$

This elasticity is increasing in $\xi$ (a higher working-capital share
makes output more sensitive to interest rates) and in $R^L_{ss}$ (a higher
steady-state lending rate amplifies the effect of any given rate change).

**Investment and capital accumulation.** Firms also invest in physical
capital. Investment $I_t$ is financed through bank loans and must satisfy
the zero-profit condition for capital:

$$\alpha \frac{Y_t}{K_t} = R^L_t - (1-\delta)$$

where the right-hand side is the net return to capital (gross lending rate
minus one plus depreciation net of the replacement cost). Log-linearizing
this condition and combining with the capital accumulation identity
$K_{t+1} = (1-\delta)K_t + I_t$, we obtain the capital recursion:

$$\hat{k}_{t+1} = (1-\delta)\hat{k}_t - \eta \left(\gamma \hat{s}_t - \Omega \hat{n}_t\right) \quad\text{(LL.3a)}$$

where $\hat{k}_t = (K_t - K_{ss})/K_{ss}$ is the proportional capital
deviation and the parameter:

$$\eta = \frac{\delta}{(1-\alpha)(R^L_{ss} - (1-\delta))}$$

is the capital-investment sensitivity. The derivation of $\eta$ follows
directly from the capital accumulation identity. Log-linearizing
$K_{t+1} = (1-\delta)K_t + I_t$ around the steady state yields:

$$\hat{k}_{t+1} = (1-\delta)\hat{k}_t + \delta \hat{i}_t$$

Log-linearizing the investment demand condition then pins the sensitivity of
investment $\hat{i}_{t}{}$ to the lending rate change, yielding the factor
$\delta / [(1-\alpha)(R^L_{ss} - (1-\delta))]$. The presence of
$\delta$ in the numerator is critical: it reflects the fact that only the
fraction $\delta$ of the capital stock is replaced each period, so a
change in the lending rate affects the flow of new investment, not the
existing stock of capital, proportionally to the depreciation rate. For our
calibrated parameters ($\alpha = 0.33$, $\delta = 0.10$, $R^L_{ss} \approx 1.08$),
we obtain $\eta \approx 0.67$, implying that a 100 basis point increase in
the lending rate reduces the capital stock by approximately
0.67 percent after one year.

---

## 3.6 Log-Linearized Transmission: The Full System

Combining the banking and production blocks, the complete log-linear
system that maps a sovereign spread path $\hat{s}_h$ (for $h = 0, 1, \ldots, H$) into
output and credit responses consists of four equations, indexed by horizon
$h$ measured in years after crisis onset.

**Horizon $h = 0$ (impact).**
At the moment of onset, the capital stock has not yet adjusted
($\hat{k}_0 = 0$). The impact responses of net worth and output are:

$$\hat{n}_0 = -B_{sens} \cdot \hat{s}_0 \quad\text{(1)}$$

$$\hat{y}_0 = -\varepsilon_p \left(\gamma \hat{s}_0 - \Omega \hat{n}_0\right) \quad\text{(2)}$$

Equation (1) states that a spread increase at onset reduces bank net worth
proportionally to sovereign bond exposure. Equation (2) states that output
falls through the working-capital channel: higher lending rates raise
effective labor costs, contracting output. The term $\Omega \hat{n}_0$
partially offsets the rate increase — when banks are also losing net worth,
they raise rates further, amplifying the output contraction.

**Horizons $h \geq 1$ (dynamic propagation).**
For subsequent periods, the full recursion governs the dynamic evolution:

$$\hat{n}_h = \Phi_N \hat{n}_{h-1} - B_{sens} \hat{s}_h - \Phi^O \gamma \hat{s}_{h-1} \quad\text{(LL.5a)}$$

$$\hat{k}_h = (1-\delta)\hat{k}_{h-1} - \eta\left(\gamma \hat{s}_{h-1} - \Omega \hat{n}_{h-1}\right) \quad\text{(LL.3a)}$$

$$\hat{y}_h = \alpha \hat{k}_h - \varepsilon_p\left(\gamma \hat{s}_h - \Omega \hat{n}_h\right) \quad\text{(LL.2a)}$$

$$\hat{\ell}_h = \frac{\lambda N_{ss}}{\ell_{ss}} \hat{n}_h \quad\text{(IR.4)}$$

The recursion in LL.5a propagates bank net worth losses forward: the
parameter $\Phi_N$ governs how quickly banks repair their balance sheets
through retained earnings and new equity. A high $\Phi_N$ means that a net
worth loss at $h=0$ persists for many periods, generating prolonged credit
tightening. The additional term $\Phi^O \gamma \hat{s}_{h-1}$ captures a
lagged effect from wholesale funding costs — banks that rely on short-term
market funding also see their liabilities repriced when sovereign spreads
rise, compressing net interest margins with a one-period delay.

The capital recursion LL.3a shows how the lending-rate wedge depresses
investment and erodes the capital stock over time. Unlike the net worth
channel (which operates through the balance sheet), the capital channel
is persistent by construction: a higher lending rate at $h=0$ reduces
investment at $h=0$, which lowers the capital stock at $h=1$, which
lowers output at $h=1$ even if the spread has already returned to its
steady-state level. This inertia is the reason why spread crises can
generate persistent output losses even when the spread itself is transitory.

**Output is reported in percentage points:** $y_h = \hat{y}_h \times 100$.

**Credit response.** The private credit response (IR.4) is proportional to
the net worth response, scaled by the leverage multiplier
$\lambda N_{ss}/\ell_{ss}$. This amplification factor is greater than one
when leverage is high (bank assets are large relative to equity), meaning
that a given percentage decline in net worth translates into a larger
percentage decline in credit. Crucially, credit and output respond to the
same underlying net worth shock but with different loadings: credit loads
exclusively on $\hat{n}_h$ (through the leverage constraint), while output
loads on both $\hat{k}_h$ (through the capital stock) and on
$\varepsilon_p(\gamma\hat{s}_{h} - \Omega\hat{n}_{h}){}$ (through working
capital costs). This difference in loadings generates the divergence in
credit-to-GDP ratios across the two types of episode, as we discuss in
Section 3.8.

---

## 3.7 The Default Regime: Capital Depletion and Autarky

When the government defaults, the country enters a period of financial
autarky. It loses access to international capital markets, which eliminates
the ability to issue new sovereign bonds and — through the sovereign ceiling
(SC) — also disrupts domestic financial intermediation. We model the autarky
period through two distinct mechanisms.

**Investment stop.** In autarky, domestic investment depends on internal
funds only. In the aggregate, the sudden stop of capital inflows forces
$I_t = 0$: domestic savings are insufficient to replace depreciating capital.
The capital law of motion then simplifies to:

$$K_{t+1} = (1-\delta) K_t \quad\text{(3)}$$

Log-linearizing with $I_{ss} = \delta K_{ss}$ (the steady-state investment
required to maintain the capital stock), the autarky capital path satisfies:

$$\hat{k}^{excl}_{h+1} = (1-\delta)\hat{k}^{excl}_h - \delta \quad\text{(4)}$$

The subtraction of $\delta$ in equation (4) reflects that zero investment
falls short of the steady-state investment level by exactly $\delta K_{ss}$,
generating a persistent decline in the capital stock below its steady-state
level. Starting from $\hat{k}^{excl}_0 = 0$ at the moment of default (the
capital stock has not yet been affected), the autarky capital path is:

$$\hat{k}^{excl}_h = -\delta \sum_{j=0}^{h-1} (1-\delta)^j = -(1 - (1-\delta)^h)$$

Capital falls by $\delta$ percent in the first year of exclusion, and by a
cumulative $(1-(1-\delta)^h)$ percent after $h$ years. For $\delta = 0.10$,
capital is 10 percent below its steady state after one year and 34 percent
below after four years of uninterrupted exclusion.

**Autarky lending rate wedge.** In addition to the investment stop, the
sovereign ceiling implies that domestic firms face an elevated financing
cost during the default episode: without the government's ability to
borrow cheaply, the private lending rate is stuck at an autarky level
$R^{L,aut} = R^L_{ss} + \Delta R^L_{aut}$, where $\Delta R^L_{aut} > 0$
is the autarky wedge. This wedge captures the disruption of financial
intermediation that accompanies sovereign default — including bank
undercapitalization, credit rationing, and the loss of trade finance —
and generates an immediate output contraction through the working-capital
channel even before capital depletion sets in.

**Survival-weighted output path.** The empirical local projection estimates
the average output response across all country-year observations following
a crisis onset, including both those that remain in default and those that
have already re-entered markets. To match this average, we weight the
output path of still-excluded countries by the survival probability
$(1-\mu)^h$ — the probability of being excluded for at least $h$ consecutive
periods — and assume that re-entered countries (fraction $1-(1-\mu)^h$)
have output approximately back at the steady-state level:

$$\bar{y}^{def}_h = (1-\mu)^h \cdot \left[\alpha \hat{k}^{excl}_h - \varepsilon_p \cdot \Delta R^L_{aut}\right] \times 100 \quad\text{(5)}$$

The autarky wedge $\Delta R^L_{aut}$ is the only parameter not set
independently: it is pinned to match the $h=0$ output observation for
default episodes, which captures the instantaneous disruption of financial
intermediation at the moment of default. All subsequent horizons ($h \geq 1$)
are genuine out-of-sample predictions.

---

## 3.8 Identifying Non-Default and Pre-Default Episodes from the Model

A central requirement for the empirical test to be valid is that the two
crisis types — non-default and default-linked — be identified from the model
itself, not from the data. We achieve this through the simulation of the
Block 1 economy.

We simulate the Arellano model for $T = 100{,}000$ periods (with a 1,000
period burn-in) to generate time series of spreads $\{s_t\}$, default
decisions $\{d_t\}$, income $\{y_t\}$, and market access indicators
$\{m_t\}$. We then collect all periods in which the sovereign spread crosses
the 1,000 basis point threshold from below while the country is in good
standing ($m_t = 1$, $s_t \geq 0.10$, $s_{t-1} < 0.10$). Each such
crossing constitutes a model crisis onset. We classify each onset as:

- **Non-default episode:** the country does not default in any of the four
  years following the onset ($d_{t+k} = 0$ for $k = 0, 1, 2, 3, 4$).
- **Pre-default episode:** the country defaults in at least one of the four
  years following the onset ($\exists k \in \{0,1,2,3,4\} : d_{t+k} = 1$).

For each type, we compute the average spread path over the event window
$h \in \{-2, -1, 0, 1, 2, 3, 4\}$ by averaging the spread $s_{t+h}$ across
all onsets of that type. This yields two model-implied spread paths,
$\hat{s}^{nd}_h$ and $\hat{s}^{def}_h$, which serve as the inputs
to the Block 2 transmission system. The model generates these two paths
endogenously from the same calibrated economy, so the difference between
them — pre-default episodes have a higher and more persistent spread, since
the country's fundamentals are worse — is a model prediction, not an
assumption.

---

## 3.9 Calibration

The model is calibrated at an annual frequency. We divide parameters into
three groups: externally calibrated constants drawn from the prior literature;
data-determined steady-state ratios computed from tranquil periods in our
panel (country-years without a crisis onset or continuation); and stochastic
process parameters estimated directly from the panel data.

**Group 1: Externally calibrated parameters.**
The household discount factor $\beta = 0.96$ implies an annual real interest
rate of approximately 4 percent, consistent with the world risk-free rate
$R^* = 0.04$. The government discount factor $\beta' = 0.95$ introduces
a small degree of government impatience relative to the household, in line
with the calibration of Arellano (2008) and Aguiar and Gopinath (2006).
The CRRA coefficient $\sigma = 2$ is standard in the international
macroeconomics literature (Aguiar and Gopinath 2007). The capital share
$\alpha = 0.33$ and annual depreciation rate $\delta = 0.10$ follow
the long-run averages for emerging market economies in the Penn World
Tables. The recovery rate $\theta = 0.62$ corresponds to one minus the
average haircut of 38 percent reported by Cruces and Trebesch (2013) for
sovereign restructurings in emerging markets. The re-entry probability
$\mu = 0.22$ implies an average exclusion period of approximately 4.5 years,
consistent with the post-default exclusion estimates of Gelos et al. (2011).
The bank leverage multiplier $\lambda = 10$ (assets-to-equity ratio) is
standard for emerging market banking sectors (Bocola 2016). The
sovereign-ceiling pass-through coefficient $\gamma = 0.80$ is taken from
the empirical estimates of Bocola (2016) for the transmission of sovereign
spreads to private lending rates in emerging markets.

**Group 2: Data-determined steady-state ratios.**
These parameters are computed from the subset of country-years without a
crisis onset or continuation in our panel — which we call "tranquil" periods
— ensuring that the steady-state calibration reflects normal economic
conditions rather than crisis dynamics. The steady-state sovereign spread
$s_{ss}$ is set equal to the cross-country median of the EMBIG spread during
tranquil periods, converted from basis points to decimal. The bank sovereign
bond exposure ratio $b^B_{ss}/Y_{ss}$ is set to the tranquil-period mean
of the BIS locational banking statistics measure of bank claims on the
government, normalized by GDP; where this variable is unavailable,
we use $b^B_{ss}/Y_{ss} = 0.15$, consistent with the IMF Financial
Soundness Indicators for emerging markets. Private credit-to-GDP
$\ell_{ss}/Y_{ss}$ is set to the tranquil-period mean of private sector
credit from the World Bank Financial Development Database. Bank net worth
is derived as $N_{ss}/Y_{ss} = (b^B_{ss} + \ell_{ss})/(\lambda \cdot Y_{ss})$,
using the balance sheet identity and the leverage constraint. Expenditure
shares ($s_C$, $s_I$, $s_G$, $s_{CA}$) and the tax rate ($\tau$) are set
to tranquil-period sample means. These data-determined parameters ensure
that the model's steady state is internally consistent with the average
characteristics of the 52 economies in our sample.

**Group 3: Estimated stochastic process parameters.**
The persistence $\rho_y$ and innovation standard deviation $\sigma_\varepsilon$
of the income process are estimated from a panel AR(1) regression of
country-specific detrended log real GDP per capita, where the trend is
removed by country-specific linear regression prior to estimation.
The spread persistence $\rho_s$ is estimated from a fixed-effects AR(1)
regression of the EMBIG spread on its one-period lag.

**SMM transmission parameters.**
The three transmission parameters $(\xi, \phi, \Phi_N)$ — governing the
working-capital share, the balance-sheet pass-through, and net worth
persistence — are jointly calibrated by Simulated Method of Moments (SMM)
to minimize the sum of squared deviations between the model-implied
non-default output path and the LP estimates at horizons $h = 0, 1, 2, 3, 4$:

$$(\xi^{\ast},\, \phi^{\ast},\, \Phi_N^{\ast}) \;=\; \arg\min_{\xi,\,\phi,\,\Phi_N} \sum_{h=0}^{4} \Bigl(\hat{y}^{nd}_{h} - \hat{\beta}^{nd}_{h}\Bigr)^2 \quad\text{(SMM)}$$

This step is the only place in the calibration where the empirical LP
estimates enter. The default path, the credit responses for both episode
types, and all qualitative predictions about the credit-to-GDP ratio are
then fully out-of-sample.

| Parameter | Value | Source |
|-----------|-------|--------|
| $\beta$ (household discount factor) | 0.96 | Standard EM calibration |
| $\beta'$ (government discount factor) | 0.95 | Arellano (2008) |
| $\sigma$ (CRRA coefficient) | 2.0 | Aguiar-Gopinath (2007) |
| $\alpha$ (capital share) | 0.33 | Penn World Tables |
| $\delta$ (depreciation rate) | 0.10 | Penn World Tables |
| $R^*$ (world risk-free rate) | 0.04 | US Treasury average |
| $\theta$ (recovery rate) | 0.62 | Cruces-Trebesch (2013) |
| $\mu$ (re-entry probability) | 0.22 | Gelos et al. (2011) |
| $\lambda$ (bank leverage) | 10 | Bocola (2016) |
| $\gamma$ (sovereign ceiling pass-through) | 0.80 | Bocola (2016) |
| $\bar{h}$ (default output cost) | 0.965 | Arellano (2008) |
| $s_{ss}$ (steady-state spread) | Median tranquil EMBIG | Our sample |
| $b^B_{ss}/Y_{ss}$ (bank sovereign exposure) | Tranquil mean | BIS / IMF FSI |
| $\ell_{ss}/Y_{ss}$ (credit-to-GDP) | Tranquil mean | World Bank |
| $\rho_y$, $\sigma_\varepsilon$ (income process) | Panel AR(1) | Our sample |
| $\xi$, $\phi$, $\Phi_N$ (transmission) | SMM | Non-default LP (Section 4) |

---

## 3.10 Testable Predictions

The model generates three predictions that can be compared to the empirical
estimates without any further tuning.

**Prediction 1 (output ranking).** Default-linked crises produce a deeper
and more persistent output contraction than non-default crises. This follows
from the superposition of channels: pre-default episodes operate through
all three transmission mechanisms (balance-sheet, working-capital, and
capital-depletion), while non-default episodes operate only through the
first two. The capital-depletion mechanism is absent in non-default crises
because investment continues — the country retains market access and can
borrow to finance capital accumulation. The model does not assume this
ordering: if the balance-sheet and working-capital channels are large enough
relative to the capital-depletion channel, non-default crises could produce
comparably large output contractions. The predicted ordering emerges from
the calibrated parameter values.

**Prediction 2 (credit-to-GDP paradox).** Non-default crises produce a
larger decline in the private credit-to-GDP ratio than default-linked
crises. In non-default episodes, the credit contraction is the primary
transmission channel: bank net worth falls through the balance-sheet
mechanism, the leverage constraint binds, and private credit contracts
faster than GDP. The credit-to-GDP ratio therefore falls. In default
episodes, GDP falls faster than credit: the dominant channel is capital
depletion (investment stops, the capital stock erodes, output falls
persistently), while credit contracts but the GDP denominator contracts
even faster. The credit-to-GDP ratio therefore falls by less, or may
even rise. This prediction is an implication of the model's structure —
the relative magnitude of the balance-sheet channel versus the
capital-depletion channel — and is not imposed by the calibration.

**Prediction 3 (default path is out-of-sample).** With transmission
parameters fixed by the SMM step on non-default episodes, the model
generates a predicted default path $\{\bar{y}^{def}_h\}$ for $h \geq 1$
(only $h=0$ is pinned by the autarky wedge). The shape — an initial
sharp contraction deepening at $h=1$ as capital depletion accumulates,
followed by partial recovery as countries exit autarky at rate $\mu$ —
is a structural prediction. The correspondence between this path and the
empirical LP estimate for default-linked episodes in Section 4 is the
primary test of the model.

---

## References

Aguiar, M. and Gopinath, G. (2006). Defaultable Debt, Interest Rates and
the Current Account. *Journal of International Economics*, 69(1), 64–83.

Aguiar, M. and Gopinath, G. (2007). Emerging Market Business Cycles: The
Cycle Is the Trend. *Journal of Political Economy*, 115(1), 69–102.

Arellano, C. (2008). Default Risk and Income Fluctuations in Emerging
Economies. *American Economic Review*, 98(3), 690–712.

Bocola, L. (2016). The Pass-Through of Sovereign Risk. *Journal of
Political Economy*, 124(4), 879–926.

Borensztein, E., Cowan, K. and Valenzuela, P. (2006). Sovereign Ceilings
"Lite"? The Impact of Sovereign Ratings on Corporate Ratings in Emerging
Market Economies. IMF Working Paper 06/75.

Cruces, J. J. and Trebesch, C. (2013). Sovereign Defaults: The Price of
Haircuts. *American Economic Journal: Macroeconomics*, 5(3), 85–117.

Gelos, G., Sahay, R. and Sandleris, G. (2011). Sovereign Borrowing by
Developing Countries: What Determines Market Access? *Journal of
International Economics*, 83(2), 243–254.

Gertler, M. and Karadi, P. (2011). A Model of Unconventional Monetary
Policy. *Journal of Monetary Economics*, 58(1), 17–34.

Jordà, O. (2005). Estimation and Inference of Impulse Responses by Local
Projections. *American Economic Review*, 95(1), 161–182.

Neumeyer, P. A. and Perri, F. (2005). Business Cycles in Emerging
Economies: The Role of Interest Rates. *Journal of Monetary Economics*,
52(2), 345–380.
