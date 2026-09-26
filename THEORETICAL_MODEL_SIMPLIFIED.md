# Theoretical Appendix

> **Scope note.** The structural model is **not part of the first version of
> the paper**, which is purely empirical. This document, its companions, and
> `do/14_calibration.do`–`do/16_model_irf.do` are retained for a later
> version; the corresponding lines in `do/00_master.do` are commented out.
> Nothing in `EMPIRICAL_ANALYSIS.md` or `RESULTS_SECTION_DRAFT.md` depends
> on any of it.

## A Two-Block Model of Sovereign Spread Crises

---

## Overview

The model is organised in two blocks that communicate sequentially. **Block 1** is a fully
nonlinear sovereign default model in the tradition of Arellano (2008): it takes the income
process and government preferences as primitives and solves — by value-function iteration —
for the endogenous default policy, the bond-price schedule, and the model-implied spread path
around crisis onsets. **Block 2** is a log-linearised banking and production model: it takes
the spread path from Block 1 as a forcing variable and maps it into output, investment, and
credit responses through the sovereign-ceiling and balance-sheet channels. This two-block
decomposition follows Bocola (2016) and has two advantages: the nonlinear object (the default
decision) is solved exactly, while the transmission object (the banking sector) is solved
analytically in closed form, making the mechanisms transparent.

The model has two regimes — a **non-default spread-crisis path** and a **default-linked path**
— that differ structurally in three ways: whether investment is possible ($I_t > 0$ vs.
$I_t = 0$), whether private credit is positive ($\ell_t > 0$ vs. $\ell_t = 0$), and whether
the sovereign-ceiling channel is active (on vs. severed by market exclusion). These structural
differences generate the two distinct transmission mechanisms that the empirical local
projections identify.

---

## Part I — Households and Firms

### 1. Households

A representative household maximises expected lifetime utility over consumption $C_t$ and
labour supply $L_t^s$:

$$\mathbb{E}_0 \sum_{t=0}^{\infty} \beta^t \left[\frac{C_t^{1-\sigma}}{1-\sigma} - \chi \frac{(L_t^s)^{1+\psi}}{1+\psi}\right] \qquad\qquad (1)$$

where $\beta \in (0,1)$ is the household's **subjective discount factor** (calibrated to
$\beta = 0.96$ at annual frequency), $\sigma > 0$ is the **inverse of the intertemporal
elasticity of substitution** (CRRA risk-aversion parameter, calibrated to $\sigma = 2$),
$\chi > 0$ is the **weight on the labour disutility**, and $\psi > 0$ is the **inverse of
the Frisch elasticity of labour supply**. The household earns wage income $w_t L_t^s$,
receives firm profits $\Pi_t$, and saves in domestic bank deposits that pay gross return
$R^D_t$. The household budget constraint is:

$$C_t + D_t = w_t L_t^s + \Pi_t + R^D_{t} D_{t-1} \qquad\qquad (2)$$

where $D_t$ denotes the stock of bank deposits held at the end of period $t$ and $R^D_t$ is
the **gross deposit rate** — the return the household earns on its savings. Optimisation over
consumption and labour supply yields the standard conditions: the **consumption Euler equation**
$C_t^{-\sigma} = \beta R^D_{t+1} \mathbb{E}_t[C_{t+1}^{-\sigma}]$ and the **labour-supply
condition** $\chi (L_t^s)^\psi = w_t C_t^{-\sigma}$.

### 2. Firms

A representative firm produces output $Y_t$ using capital $K_t$ and labour $L_t$ according to
a **Cobb-Douglas production function**:

$$Y_t = A K_t^\alpha L_t^{1-\alpha} \qquad\qquad (3)$$

where $A > 0$ is **total factor productivity** (normalised to one at the steady state),
$\alpha \in (0,1)$ is the **capital share** (calibrated to $\alpha = 0.33$ following the
standard macro literature on emerging markets), and $1-\alpha$ is the labour share. The firm
takes the wage $w_t$ and the domestic lending rate $R^L_t$ as given.

**The working-capital constraint.** Before production takes place, the firm must
**pre-finance a fraction $\xi \in (0,1)$ of its total wage bill** $w_t L_t$ by borrowing from
domestic banks at the lending rate $R^L_t$. The working-capital loan $\xi w_t L_t$ must be
repaid at rate $R^L_t$ at the end of the period. The effective unit cost of labour for the
firm is therefore $(w_t + \xi w_t (R^L_t - 1)) = w_t(1 + \xi(R^L_t - 1))$, which is
increasing in $R^L_t$: any rise in the lending rate makes labour more expensive. The firm's
profit is:

$$\Pi_t = Y_t - w_t(1 + \xi(R^L_t - 1))L_t - R^L_t K_t \qquad\qquad (4)$$

Maximising (4) over $L_t$ and $K_t$ yields the **labour demand condition**:

$$(1-\alpha) A K_t^\alpha L_t^{-\alpha} = w_t(1 + \xi(R^L_t-1)) \qquad\qquad (5)$$

which equates the marginal product of labour to its effective cost, and the **capital demand
condition**:

$$\alpha A K_t^{\alpha-1} L_t^{1-\alpha} = R^L_t \qquad\qquad (6)$$

which equates the marginal product of capital to the lending rate. Substituting (5) into
the production function and solving for $Y_t$ yields the **output supply function**:

$$Y_t = A K_t^\alpha \cdot \Psi(R^L_t), \qquad\qquad \Psi(R^L_t) \equiv \left[\frac{1-\alpha}{1 + \xi(R^L_t - 1)}\right]^{n}, \quad n = \frac{1-\alpha}{\alpha} \qquad\qquad (7)$$

where $\Psi(R^L_t)$ is the **working-capital distortion factor** — a decreasing function
of the lending rate that captures how the cost of pre-financing the wage bill compresses
output below its frictionless level. When $\xi = 0$ (no working-capital requirement), $\Psi = 1$
and output is purely determined by capital and TFP. When $\xi > 0$, any rise in $R^L_t$
reduces $\Psi$ and therefore $Y_t$ even for a given capital stock. Log-linearising (7) around
the non-crisis steady state (where $R^L_t = R^{L,ss}$, $K_t = K^{ss}$):

$$\hat{y}_t = \alpha\hat{k}_t - \varepsilon_p(R^L_t - R^{L,ss}) \qquad\qquad (8)$$

where $\hat{y}_t \equiv \log(Y_t/Y^{ss})$ and $\hat{k}_t \equiv \log(K_t/K^{ss})$ are
log-deviations from steady state, and:

$$\varepsilon_p \equiv \frac{\xi\, n\, R^{L,ss}}{1 + \xi(R^{L,ss}-1)} > 0 \qquad\qquad (9)$$

is the **working-capital output elasticity** — the percentage-point fall in output for each
unit rise in the lending rate, through the channel of more expensive working-capital finance.
Equation (8) reveals two channels: the **capital channel** $\alpha\hat{k}_t$ (zero at impact
since capital is predetermined, accumulating over time as investment responds to higher
lending rates) and the **working-capital channel** $-\varepsilon_p(R^L_t - R^{L,ss})$
(operating immediately at $h=0$ without any delay).

### 3. Capital Accumulation and Investment

The capital stock evolves according to:

$$K_{t+1} = (1-\delta)K_t + I_t \qquad\qquad (10)$$

where $\delta \in (0,1)$ is the **annual depreciation rate** (calibrated to $\delta = 0.10$)
and $I_t \geq 0$ is gross investment. The firm's investment decision follows from the
**user-cost condition**: capital is accumulated until its net return equals the lending rate.
Log-linearising around the steady state:

$$\hat{k}_{t+1} = (1-\delta)\hat{k}_t - \eta(R^L_t - R^{L,ss}) \qquad\qquad (11)$$

where $\eta \equiv 1/[(1-\alpha)(R^{L,ss} - (1-\delta))] > 0$ is the **interest
semi-elasticity of the capital stock** — the fall in next period's capital for each unit rise
in the lending rate. Equation (11) is the capital accumulation equation used throughout the
transmission block: a rise in $R^L_t$ reduces investment today, which reduces $K_{t+1}$ and
thereby $Y_{t+1}$ through the capital channel in (8). This intertemporal link generates the
**progressive deepening** of the output response over horizons $h=1,2,3,4$ even as the
spread itself decays.

---

## Part II — Banks and the Financial Sector

### 4. Bank Balance Sheet

A representative bank intermediates between the household (which holds deposits $D_t$),
foreign wholesale funding markets (which supply funds $F_t$ at the international rate $R^O_t$),
and the real sector (private loans $\ell_t$) and the government (sovereign bonds $b^B_t$).
The **balance-sheet identity** is:

$$b^B_t + \ell_t = D_t + F_t + N_t \qquad\qquad (12)$$

where $N_t \geq 0$ is **bank net worth** (equity capital, the residual claim after all
liabilities are repaid). Total assets on the left equal total liabilities plus equity on the
right. The bank holds sovereign bonds $b^B_t$ at market price $q_t$ (so the market value of
holdings is $q_t \cdot \tilde{b}^B_t$, where $\tilde{b}^B_t$ is face value; we use
$b^B_t$ to denote market value throughout) and extends private loans $\ell_t$ at rate $R^L_t$.

### 5. The Leverage Constraint

A regulatory or market-imposed **leverage constraint** limits total bank assets to a multiple
$\lambda > 1$ of equity:

$$b^B_t + \ell_t \leq \lambda N_t \qquad\qquad (13)$$

where $\lambda = 10$ is the **bank leverage multiplier** (calibrated from EM bank balance-sheet
data; typical EM commercial banks operate at assets-to-equity ratios of 8–12). The constraint
(13) is **always binding in equilibrium**: profit-maximising banks fully deploy their
balance-sheet capacity because the return on assets exceeds the cost of liabilities. With the
constraint binding with equality, total bank assets are exactly $\lambda N_t$, and the
composition between $b^B_t$ and $\ell_t$ is the bank's choice variable.

### 6. Net Worth Dynamics

Bank net worth evolves through **retained earnings** minus **capital losses on the sovereign
bond portfolio**. When the sovereign spread $s_t$ rises, the market value of sovereign bonds
$q_t$ falls, generating a mark-to-market **capital loss** proportional to the bank's
pre-existing holdings. If the bank holds $b^{B,ss}$ of sovereign bonds at the steady state,
a unit rise in the spread causes a capital loss of approximately
$b^{B,ss}/(1+s^{ss})^2$ on the asset side of the balance sheet, which flows directly through
to net worth. Log-linearising net worth around the steady state:

$$\hat{n}_{t+1} = \Phi^N \hat{n}_t - \mathrm{B}\,\hat{s}_{t+1} \qquad\qquad (14)$$

where $\hat{n}_t \equiv \log(N_t/N^{ss})$ is the **net-worth log-deviation** from steady state,
$\Phi^N \in (0,1)$ is the **net-worth persistence parameter** — the fraction of last period's
net-worth deviation that survives into this period through retained earnings, after paying
interest on liabilities and dividends — and:

$$\mathrm{B} \equiv \frac{b^{B,ss}}{N^{ss}(1+s^{ss})^2} > 0 \qquad\qquad (15)$$

is the **balance-sheet sensitivity**: how much net worth falls (in log-deviation terms) per
unit rise in the spread, proportional to the steady-state ratio of sovereign bond holdings
to equity $b^{B,ss}/N^{ss}$. A bank with large pre-existing sovereign exposure loses more
net worth for the same spread shock — this is the **diabolic loop** between sovereign stress
and banking-sector fragility documented in Brunnermeier et al. (2016). The term
$-\mathrm{B}\hat{s}_{t+1}$ captures this: the fresh capital loss at each period $t+1$ is
proportional to the spread at that period and to the pre-existing exposure ratio.

### 7. Private Credit Supply

With the leverage constraint (13) binding with equality and sovereign bond holdings
$b^B_t$ partially sticky in the short run (banks cannot instantaneously unwind illiquid
sovereign bond positions), private lending adjusts as the **residual** on the bank's
balance sheet. Log-linearising the binding constraint and holding $\hat{b}^B_t$ fixed at zero
in the non-default path (no systematic change in sovereign bond holdings during non-default
episodes):

$$\hat{\ell}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\,\hat{n}_t \qquad\qquad (16)$$

where $\hat{\ell}_t \equiv \log(\ell_t/\ell^{ss})$ is the log-deviation of private credit
from steady state. The coefficient $\lambda N^{ss}/\ell^{ss} > 1$ is a **leverage amplifier**:
since total assets are $\lambda N^{ss}$ and private loans are only $\ell^{ss} < \lambda N^{ss}$
(the rest being sovereign bonds), a 1% fall in net worth produces a more-than-1% fall in
private lending. This amplification is larger when banks hold more sovereign bonds relative
to private loans at the steady state — the same factor that determines $\mathrm{B}$ in (15).

### 8. The Sovereign Ceiling and the Domestic Lending Rate

Banks fund their sovereign bond holdings and private loans partly through **international
wholesale markets** — interbank borrowing, bond issuance, repo. The sovereign ceiling is the
empirical regularity that no private borrower in a country can obtain external funding at a
rate below its own government's borrowing rate. This arises because foreign investors view the
sovereign as the ultimate guarantor of the private sector's obligations, so private default
risk cannot be lower than sovereign default risk. Formally, the **international wholesale
funding rate** $R^O_t$ faced by domestic banks is bounded below by the sovereign rate:

$$R^O_t = R^* + \gamma s_t \qquad\qquad (17)$$

where $R^* > 0$ is the **international risk-free rate** (calibrated to $R^* = 0.04$,
approximately the US 10-year Treasury yield over the sample), $s_t \geq 0$ is the **sovereign
EMBIG spread** in decimal form (the excess yield on sovereign Eurobonds over US Treasuries),
and $\gamma \in (0,1]$ is the **sovereign-ceiling pass-through parameter** — the fraction of
the sovereign spread that is transmitted to the bank's wholesale funding cost. When $\gamma = 1$,
the full sovereign spread passes through to banks; when $\gamma < 1$, domestic-currency
mismatch or regulatory buffers dampen the transmission.

Domestic banks price private loans at a **mark-up** over their funding cost that also reflects
the balance-sheet risk premium arising from their sovereign bond exposure. A bank with a high
ratio $b^B_t/N_t$ of sovereign bonds to net worth charges a higher lending rate because it
is more exposed to sovereign default risk and must compensate lenders for this risk. The
**domestic lending rate** is:

$$R^L_t = R^* + \gamma s_t + \varphi\, \frac{b^B_t}{N_t} \qquad\qquad (18)$$

where $\varphi > 0$ is the **balance-sheet pass-through parameter** — the sensitivity of the
lending rate to the bank's sovereign exposure ratio $b^B_t/N_t$. Log-linearising (18) around
the steady state (where $b^B_t/N_t = b^{B,ss}/N^{ss}$):

$$R^L_t - R^{L,ss} = \gamma\hat{s}_t - \Omega\hat{n}_t \qquad\qquad (19)$$

where $\Omega \equiv \varphi\, b^{B,ss}/N^{ss} > 0$ is the **balance-sheet amplification
coefficient** — the sensitivity of the lending rate to net-worth erosion, proportional to the
steady-state sovereign exposure ratio. Equation (19) has two terms: (i) the **direct
sovereign-ceiling channel** $\gamma\hat{s}_t > 0$ — the rise in wholesale funding costs
driven by the higher sovereign spread, operating on all banks regardless of their individual
sovereign exposure; and (ii) the **indirect balance-sheet channel** $-\Omega\hat{n}_t > 0$
(since $\hat{n}_t < 0$ during a crisis) — the additional lending-rate premium that arises as
net-worth erosion raises the bank's effective leverage in sovereign bonds, increasing the risk
premium charged to private borrowers.

### 9. The Complete Transmission System (Non-Default Path)

Substituting (19) into (8) and (11) and combining with (14) and (16), the complete
log-linearised transmission system for the **non-default path** is:

$$\hat{y}_t = \alpha\hat{k}_t - \varepsilon_p\left(\gamma\hat{s}_t - \Omega\hat{n}_t\right) \qquad\qquad \text{(LL.2a)}$$

$$\hat{k}_{t+1} = (1-\delta)\hat{k}_t - \eta\left(\gamma\hat{s}_t - \Omega\hat{n}_t\right) \qquad\qquad \text{(LL.3a)}$$

$$\hat{n}_{t+1} = \Phi^N\hat{n}_t - \mathrm{B}\,\hat{s}_{t+1} \qquad\qquad \text{(LL.5a)}$$

$$\hat{\ell}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\,\hat{n}_t \qquad\qquad \text{(LL.6)}$$

These four equations — driven solely by the exogenous spread path $\{\hat{s}_t\}$ from
Block 1 — constitute the complete non-default transmission system. The system has three
free parameters $(\gamma, \varphi, \Phi^N)$ that are calibrated by SMM (Section 10).
All other parameters ($\alpha$, $\delta$, $\xi$, $\lambda$, $\mathrm{B}$, $\Omega$, $\varepsilon_p$,
$\eta$) are determined externally from the literature or from steady-state data moments.

---

## Part III — Block 1: The Nonlinear Default and Spread Block

### 10. Income Process

The country's income $y_t$ follows a stationary **AR(1) process in logarithms**:

$$\log y_{t+1} = \rho_y \log y_t + \varepsilon_{t+1}, \qquad \varepsilon_{t+1} \overset{iid}{\sim} \mathcal{N}(0,\sigma_\varepsilon^2) \qquad\qquad (20)$$

where $\rho_y \in (0,1)$ is the **income persistence** and $\sigma_\varepsilon > 0$ is the
**innovation standard deviation**. Both are estimated from the panel data by an AR(1)
regression on country-detrended log GDP per capita (country fixed effects plus country-specific
linear trends are removed before estimation). The process is discretised using the
**Tauchen (1986)** method into $N_y = 21$ income states $\{y_1, \ldots, y_{N_y}\}$ with
$N_y \times N_y$ **transition matrix** $\mathbf{P}$, where $P_{ij} = \Pr(y_{t+1} = y_j \mid y_t = y_i)$.
The grid spans $\pm 3$ unconditional standard deviations of $\log y$.

### 11. The Government's Default Problem

The government issues **one-period discount bonds** with face value $B_{t+1}$ (the amount
owed next period) and receives price $q_t$ per unit today, so it raises $q_t B_{t+1}$ in
revenue from bond issuance. It finances primary expenditure $G_t$ through tax revenue $\tau Y_t$
(where $\tau$ is the tax rate) and bond issuance, subject to the **government budget
constraint**:

$$G_t = \tau Y_t - B_t + q_t(B_{t+1}, y_t)\, B_{t+1} \qquad\qquad (21)$$

where $B_t$ is the face value of debt maturing this period (which must be repaid in full if
the government does not default) and $q_t(B_{t+1}, y_t)$ is the bond price — endogenous
because it reflects the probability that next period's government will choose to default on
$B_{t+1}$.

The government maximises the discounted utility of primary expenditure:

$$\mathbb{E}_0 \sum_{t=0}^{\infty} \beta_p^t U(G_t) \qquad\qquad (22)$$

where $U(G) = G^{1-\sigma}/(1-\sigma)$ is a **CRRA utility function** over government
spending and $\beta_p \in (0,1)$ is the **government's discount factor**, which may differ
from the household's $\beta$. A lower $\beta_p$ captures political short-termism — governments
place excessive weight on current expenditure relative to future debt repayment, creating a
**temptation to default** even when repayment is feasible. The government's problem is
recursive in the **state variables** $(B_t, y_t)$. The **value of repayment** satisfies:

$$V^R(B_t, y_t) = \max_{G_t,\; B_{t+1}\, \leq\, \bar{B}(y_t)} \left\{ U(G_t) + \beta_p\, \mathbb{E}_t\left[V(B_{t+1}, y_{t+1})\right] \right\} \qquad\qquad (23)$$

subject to (21), where $\bar{B}(y_t)$ is the **natural borrowing limit** (the maximum debt
the government can commit to repay) and $V(B_{t+1}, y_{t+1}) = \max\{V^R(B_{t+1}, y_{t+1}),\,V^D(y_{t+1})\}$
is the next-period **continuation value** — the government will optimally choose repayment or
default next period depending on which delivers higher utility.

### 12. The Value of Default and Autarky

Upon default, the government is **excluded from international capital markets**: it cannot
issue new bonds, so it cannot smooth consumption across periods through debt. During exclusion,
primary expenditure is financed entirely from within-period tax revenue, subject to a
**direct output cost of exclusion** $\phi > 0$ that captures the real disruption to trade
finance, imports of intermediate inputs, and foreign direct investment that typically
accompanies a sovereign default:

$$G^{\text{aut}}_t = (1-\phi)\tau Y_t \qquad\qquad (24)$$

The country **re-enters international capital markets with probability $\mu$ each period**
(calibrated to $\mu = 0.22$, corresponding to an average exclusion period of approximately
4.5 years, as documented by Gelos et al. 2011). Upon re-entry, the country's debt is reset
to zero ($B = 0$) following the haircut settlement. The **value of default** satisfies:

$$V^D(y_t) = U\left((1-\phi)\tau Y_t\right) + \beta_p\, \mathbb{E}_t\left[\mu\, V^R(0, y_{t+1}) + (1-\mu)\, V^D(y_{t+1})\right] \qquad\qquad (25)$$

where the continuation value combines: with probability $\mu$, the country re-enters markets
with $B = 0$ and receives $V^R(0, y_{t+1})$; with probability $(1-\mu)$, it remains excluded
for another period and receives $V^D(y_{t+1})$. Note that $V^D$ depends only on $y_t$ (not
$B_t$) because debt obligations are suspended during default.

### 13. The Optimal Default Decision

The government defaults whenever the value of default strictly exceeds the value of repayment.
The **default set** $\mathcal{D}(B_t)$ is the set of income realisations at which default is
optimal given debt level $B_t$:

$$\mathcal{D}(B_t) = \left\{y_t : V^D(y_t) \geq V^R(B_t, y_t)\right\} \qquad\qquad (26)$$

The **default indicator** is $d_t = \mathbf{1}[y_t \in \mathcal{D}(B_t)]$. Two comparative
statics follow immediately. First, $\mathcal{D}(B_t)$ is **increasing in $B_t$**: higher debt
raises the cost of repayment $V^R$ (through tighter budget constraints) without affecting $V^D$,
making default relatively more attractive. Second, $\mathcal{D}(B_t)$ is **decreasing in $y_t$**:
higher income raises the value of repayment (through higher primary expenditure while remaining
in markets) and raises the cost of exclusion (through the $h_{\text{def}}\bar{y}$ threshold
in the autarky income function, see below), making default less attractive. These two properties
together generate defaults that occur at high debt during income downturns — the empirical
pattern documented across EM sovereign crises.

The default cost function follows Arellano (2008): autarky income is $y^{\text{def}}_t = \min(y_t, h_{\text{def}}\bar{y})$,
where $h_{\text{def}} < 1$ and $\bar{y} = \mathbb{E}[y]$ is mean income. When income is below
$h_{\text{def}}\bar{y}$, the output cost of default is zero (the country is already depressed);
when income is above the threshold, the government suffers a real output loss $h_{\text{def}}\bar{y} - y_t < 0$
if it defaults. This **asymmetric cost** creates a stronger incentive to default in bad times
(when the cost is zero) than in good times (when the cost is positive), amplifying the
countercyclicality of default.

### 14. The Bond Price Fixed Point

Foreign investors are **risk-neutral** and can invest at the world risk-free rate $R^*$.
Zero-profit on bond holdings requires the bond price $q_t(B_{t+1}, y_t)$ to equal the
expected discounted repayment. If the government defaults next period (with probability
$\pi_{t+1} = \Pr(d_{t+1}=1 \mid B_{t+1}, y_t)$), the investor recovers a fraction
$\theta \in (0,1)$ of face value — the **recovery rate** (calibrated to $\theta = 0.62$,
corresponding to a mean haircut of 38% as documented by Cruces and Trebesch 2013 for a
broad sample of EM restructurings). The **bond price schedule** is:

$$q_t(B_{t+1}, y_t) = \frac{1}{1+R^*}\, \mathbb{E}_t\left[(1-d_{t+1})\cdot 1 + d_{t+1}\cdot\theta\right] \qquad\qquad (27)$$

$$= \frac{1}{1+R^*}\left[(1 - \pi_{t+1}) + \pi_{t+1}\theta\right] = \frac{1 - \pi_{t+1}(1-\theta)}{1+R^*}$$

The **sovereign spread** implied by the bond price is the excess return over the risk-free rate:

$$s_t(B_{t+1}, y_t) = \frac{1}{q_t(B_{t+1}, y_t)} - 1 - R^* \qquad\qquad (28)$$

This spread is **endogenous and nonlinear**: it rises sharply as $B_{t+1}$ approaches the
natural debt limit $\bar{B}(y_t)$ (because default probability approaches one) and declines
as income $y_t$ rises (because repayment becomes more attractive). Equations (23)–(28)
constitute a **fixed-point problem**: the bond price depends on the default policy (equation 26),
the default policy depends on the value functions (equations 23 and 25), and the value
functions depend on the bond price through the budget constraint (21). This fixed point
is solved by value-function iteration.

### 15. Numerical Solution: Value-Function Iteration

The model is solved on a discrete state space: a uniform debt grid
$\{B_1, \ldots, B_{N_B}\}$ with $N_B = 200$ points on $[0,\, 1.5\bar{y}]$, and the
Tauchen-discretised income grid $\{y_1, \ldots, y_{N_y}\}$ with $N_y = 21$ points. Starting
from an initial guess for the value functions and bond price schedule
(typically $V^R = V^D = 0$ and $q = 1/(1+R^*)$), the algorithm iterates the following steps
until convergence:

**Step 1 — Update value of default** (equation 25): for each income state $y_j$, compute
$V^D_{\text{new}}(y_j)$ using the current guess for $V$ and the transition matrix $\mathbf{P}$:
$$V^D_{\text{new}}(y_j) = U(y^{\text{def}}_j) + \beta_p \sum_{k=1}^{N_y} P_{jk}\left[\mu\, V(B_1, y_k) + (1-\mu)\, V^D(y_k)\right]$$

**Step 2 — Update value of repayment** (equation 23): for each $(B_i, y_j)$, compute the
consumption matrix $c(B_i, B'_m, y_j) = \tau y_j - B_i + q(B'_m, y_j)\cdot B'_m$ for all
candidate borrowing choices $B'_m$, apply CRRA utility, add discounted continuation
$\beta_p\sum_k P_{jk} V(B'_m, y_k)$, and take the maximum over $m$:
$$V^R_{\text{new}}(B_i, y_j) = \max_{m \in \{1,\ldots,N_B\}} \left\{U(c(B_i,B'_m,y_j)) + \beta_p\sum_{k=1}^{N_y}P_{jk}V(B'_m,y_k)\right\}$$

**Step 3 — Update default policy** (equation 26): set $d_{\text{new}}(B_i, y_j) = 1$ if
$V^D_{\text{new}}(y_j) \geq V^R_{\text{new}}(B_i, y_j)$, else 0.

**Step 4 — Update bond price** (equation 27): compute
$q_{\text{new}}(B'_m, y_j) = \frac{1}{1+R^*}\sum_k P_{jk}[(1-d_{\text{new}}(B'_m,y_k)) + d_{\text{new}}(B'_m,y_k)\cdot\theta]$.
Apply dampening $q \leftarrow 0.5\,q_{\text{new}} + 0.5\,q_{\text{old}}$ for stability.

**Step 5 — Check convergence**: if $\max(||V_{\text{new}}-V_{\text{old}}||_\infty,\,||q_{\text{new}}-q_{\text{old}}||_\infty) < 10^{-6}$,
stop. Otherwise return to Step 1.

Convergence is typically achieved in fewer than 500 iterations.

### 16. Simulation and Crisis Event Study

After solving the fixed point, the model is **simulated for $T = 100{,}000$ periods**
(discarding the first 1,000 burn-in periods) to compute model moments and the crisis spread
path. The simulation draws income $y_t$ from the Tauchen transition matrix, applies the
optimal default policy, and records the spread $s_t(B'_t, y_t)$ on newly issued bonds each
period. A **spread crisis onset** is defined as the first period the spread crosses the
1,000 bps threshold from below while the country is in good standing ($d_t = 0$). Around
each onset at date $t^*$, the spread path at horizons $h = -2, \ldots, 4$ is recorded.
Averaging across all crisis events yields the **model-implied spread event path**
$\{\bar{s}^{\text{model}}_h\}_{h=-2}^{4}$, which is then passed to Block 2 as the forcing
variable in (LL.2a)–(LL.6).

**Model moments targeted:**

| Moment | Model target | Data source |
|---|---|---|
| Annual default frequency | $\approx 3\%$ | Sample mean, 52 EM countries 1994–2025 |
| Mean spread (bps) | $\approx 350$ bps | Tranquil-period median EMBIG |
| Spread standard deviation | $\approx 200$ bps | Panel standard deviation |
| Mean debt/GDP in market | $\approx 0.45$ | Tranquil-period mean |

Parameters $\beta_p$ and $h_{\text{def}}$ are jointly tuned to match default frequency and
mean spread. The other Block 1 parameters ($\theta = 0.62$, $\mu = 0.22$, $\rho_y$,
$\sigma_\varepsilon$) are fixed prior to simulation.

---

## Part IV — Impulse Responses and the Two Paths

### 17. Non-Default Path: Closed-Form Impulse Responses

The spread path from Block 1 follows an **AR(1) process** around the onset:
$\hat{s}_t = \rho_s^t \hat{s}_0$, where $\rho_s$ is the spread persistence estimated from the
panel (see Section 3) and $\hat{s}_0 = s_0 - s^{ss} > 0$ is the onset-period spread
deviation from steady state. Starting from the non-crisis steady state ($\hat{k}_0 = 0$,
$\hat{n}_0 = -\mathrm{B}\hat{s}_0$ from equation LL.5a at $h=0$), iterating forward yields:

**Net worth at horizon $h$:**

$$\hat{n}^{nd}_h = (\Phi^N)^h \hat{n}_0 - \sum_{h'=1}^{h}(\Phi^N)^{h-h'}\,\mathrm{B}\,\rho_s^{h'}\hat{s}_0 \qquad\qquad \text{(IR.1)}$$

The first term $(\Phi^N)^h \hat{n}_0$ is the **decaying legacy** of the initial balance-sheet
shock; the summation is the **fresh capital loss** accumulated at each period $h' = 1,\ldots,h$
as the still-elevated spread continues to erode net worth. Net worth reaches its trough when
the retained-earnings recovery term $\Phi^N|\hat{n}_{h-1}|$ first outweighs the fresh loss
$\mathrm{B}\rho_s^h\hat{s}_0$.

**Capital stock at horizon $h$** (iterating LL.3a with $\hat{k}_0 = 0$):

$$\hat{k}^{nd}_h = -\eta\gamma\sum_{h'=1}^{h}(1-\delta)^{h-h'}\rho_s^{h'-1}\hat{s}_0 + \eta\Omega\sum_{h'=1}^{h}(1-\delta)^{h-h'}\hat{n}_{h'-1} \qquad\qquad \text{(IR.2)}$$

The capital stock deviation accumulates as the discounted sum of past lending-rate increases,
weighted by the depreciation factor $(1-\delta)$. The first sum is driven by the
**sovereign-ceiling channel** ($\eta\gamma$ term) and the second by the **balance-sheet
amplification** ($\eta\Omega$ term, negative since $\hat{n}_t < 0$).

**Output at horizon $h$:**

$$\hat{y}^{nd}_h = \underbrace{\alpha\hat{k}^{nd}_h}_{\text{capital channel}} - \underbrace{\varepsilon_p\gamma\rho_s^h\hat{s}_0}_{\text{working-capital channel}} + \underbrace{\varepsilon_p\Omega\hat{n}^{nd}_h}_{\text{balance-sheet amplification}} \qquad\qquad \text{(IR.3)}$$

All three terms are negative ($\hat{k}^{nd}_h < 0$, $\hat{s}_0 > 0$, $\hat{n}^{nd}_h < 0$).
The working-capital channel ($h=0$ impact: $-\varepsilon_p\gamma\hat{s}_0$) is immediate;
the capital and balance-sheet channels deepen over time, generating the trough at $h=1$
(consistent with the empirical estimate of $-3.76$ pp at $h=1$).

**Private credit at horizon $h$:**

$$\hat{\ell}^{nd}_h = \frac{\lambda N^{ss}}{\ell^{ss}}\,\hat{n}^{nd}_h \qquad\qquad \text{(IR.4)}$$

Credit is negligible at $h=0$ (since $\hat{n}_0 = -\mathrm{B}\hat{s}_0$ is small when
$b^{B,ss}/N^{ss}$ is moderate) and deepens progressively, reaching the empirical estimate
of $-4.78$ pp at $h=4$. **In the non-default path, private credit $\hat{\ell}_t$ falls
faster than output $\hat{y}_t$, so the ratio private credit/GDP falls substantially — the
credit channel is the empirically dominant transmission mechanism.**

### 18. Default Path: Three Structural Breaks

Under default-linked episodes, three simultaneous structural breaks modify the transmission
system (LL.2a)–(LL.6). These breaks are not choices made by the model; they are the
real-world consequences of market exclusion that are imposed at the onset date.

**Break 1 — Investment suspended ($I_t = 0$).** Market exclusion eliminates the country's
access to external capital flows. With domestic savings insufficient to finance the pre-crisis
level of investment, and with the IMF or bilateral creditors typically conditioning any support
on fiscal adjustment (which further crowds out private investment), gross investment falls to
zero during the exclusion period. Capital therefore depletes at the depreciation rate with no
offsetting investment. From equation (10):

$$K^{\text{def}}_{t+h} = (1-\delta)^h K_0 \quad \Rightarrow \quad \hat{k}^{\text{def}}_h = h\ln(1-\delta) \qquad\qquad \text{(IR.5a)}$$

This is a **cumulative and strictly worsening** mechanism: each period without investment
reduces the capital stock by a further 10%, dragging output down proportionally through the
capital channel in (LL.2a). At $h=1$ the capital loss is $\approx 10\%$, at $h=4$ it is
$\approx 34\%$, generating an output loss through the capital channel of
$\alpha \times 34\% \approx 11$ pp by $h=4$ — even before accounting for the working-capital
channel.

**Break 2 — Private lending collapses: gambling for resurrection ($\ell^{\text{def}}_t = 0$).**
When the government defaults, the value of sovereign bonds $b^B_t$ on bank balance sheets
falls to the recovery value $\theta \times$ face value. Banks' net worth $N_t$ is eroded.
However, rather than reducing total assets (which would require repaying wholesale funding
or accepting equity dilution), banks optimally **redirect their entire balance-sheet capacity
toward sovereign bonds**: $b^{B,\text{def}}_t = \lambda N_t$, $\ell^{\text{def}}_t = 0$.
This is profitable when the expected recovery on sovereign bonds (probability $\mu\theta$ per
period) exceeds the net return on private lending, which happens precisely during a crisis
when private default risk is also elevated. The strategy is termed **gambling for resurrection**
because the bank's value is maximised by concentrating risk — if the country recovers and
sovereign bonds are repaid at par, the bank survives; if the country restructures at $\theta$,
the bank would have failed anyway.

The consequence is that **private credit collapses to zero**, and firms cannot pre-finance
their wage bill. The working-capital distortion is maximal: $\Psi(R^{L,\text{aut}}) \ll 1$,
where $R^{L,\text{aut}}$ is the autarky lending rate (set by the domestic closed-economy
equilibrium). The output equation under default becomes:

$$\hat{y}^{\text{def}}_h = \alpha \cdot h\ln(1-\delta) - \varepsilon_p\,\Delta r^{L,\text{aut}} \qquad\qquad \text{(IR.5)}$$

where $\Delta r^{L,\text{aut}} \equiv R^{L,\text{aut}} - R^{L,ss} > 0$ is the **autarky
lending-rate wedge** — the excess of the autarky rate above the steady-state rate. This
wedge is determined in the closed-economy equilibrium by the capital market clearing
condition: the marginal product of capital equals the lending rate. Since capital is already
depleted ($K^{\text{def}} < K^{ss}$), its marginal product is above its steady-state value,
so $R^{L,\text{aut}} > R^{L,ss}$. The wedge $\Delta r^{L,\text{aut}}$ is pinned to the
**single $h=0$ empirical default moment** from the local projections, leaving
$h = 1, \ldots, 4$ as **genuine out-of-sample predictions**.

**Break 3 — Sovereign ceiling severed.** With wholesale funding markets closed ($F_t = 0$),
the international rate $R^O_t$ is no longer relevant: banks cannot borrow externally at any
rate. The sovereign-ceiling channel $\gamma\hat{s}_t$ in equation (LL.2a) is therefore
**inoperative during exclusion**. The domestic lending rate is set entirely by the autarky
equilibrium, independently of the sovereign spread. This is why the deepening of the default
recession is driven by capital depletion (equation IR.5a), not by the spread rising further —
the spread channel is severed at the moment of default.

### 19. The Credit/GDP Paradox

A central empirical finding is that **private credit/GDP contracts more in non-default
episodes** despite private lending collapsing to zero under default ($\ell^{\text{def}}_t = 0$,
equation LL.6D). This is resolved by the **GDP-denominator effect**.

The empirical outcome variable is the **ratio** $\ell_t/Y_t$. Under default, the output
collapse is driven by capital depletion (IR.5a) — a cumulative mechanism that adds
approximately $\alpha\delta = 3.3$ pp of output loss per year. By $h=3$, capital depletion
alone has pushed GDP down by $\alpha \times |3\ln(1-\delta)| \approx 9$ pp. The **GDP
denominator therefore shrinks at least as fast as the private-credit numerator**, leaving
the ratio $\ell_t/Y_t$ approximately unchanged. Under non-default, the reverse holds: the
credit channel dominates, $\hat{\ell}_t$ falls much faster than $\hat{y}_t$, and the ratio
falls substantially.

$$\underbrace{|\Delta(\ell/Y)|^{nd} > |\Delta(\ell/Y)|^{def}}_{\text{empirical finding}} \qquad \text{because} \qquad \begin{cases} \text{ND:} & |\Delta\ell| \gg |\Delta Y| \quad (\text{credit channel dominates}) \\ \text{DEF:} & |\Delta Y| \gg |\Delta\ell| \quad (\text{capital-depletion channel dominates}) \end{cases}$$

**Proposition (Credit/GDP Paradox).** A smaller measured credit/GDP response in
default-linked episodes is not evidence that domestic banks are lending more during default
than during non-default spread crises. It is evidence that GDP is collapsing faster than
credit, through the capital-depletion channel (IR.5a) that is absent in non-default episodes.
The two paths differ not in the intensity of the same mechanism but in the identity of the
dominant mechanism.

---

## Part V — Calibration and Fit

### 20. Parameter Groups

**Group 1 — Externally calibrated.** Standard values from the EM-DSGE literature:

| Parameter | Value | Source |
|---|---|---|
| $\beta$ | 0.96 | Household discount factor (annual) |
| $\beta_p$ | 0.95 | Government discount factor (tuned to default frequency) |
| $\sigma$ | 2.0 | CRRA risk aversion |
| $\alpha$ | 0.33 | Capital share |
| $\delta$ | 0.10 | Annual depreciation |
| $R^*$ | 0.04 | World risk-free rate |
| $\theta$ | 0.62 | Recovery rate (Cruces-Trebesch 2013) |
| $\mu$ | 0.22 | Re-entry probability (Gelos et al. 2011) |
| $\lambda$ | 10 | Bank leverage |

**Group 2 — Data-determined steady-state ratios.** Computed from tranquil-period means
(country-years with no spread crisis onset and no continuation) in the panel:

| Moment | Value | Interpretation |
|---|---|---|
| $s^{ss}$ | Median tranquil spread | Steady-state EMBIG spread |
| $b^{B,ss}/Y^{ss}$ | Mean claims on govt / GDP | Sovereign bond exposure |
| $\ell^{ss}/Y^{ss}$ | Mean private credit / GDP | Steady-state credit ratio |
| $N^{ss}/Y^{ss}$ | $(b^{B,ss} + \ell^{ss})/(\lambda Y^{ss})$ | Implied steady-state net worth |

**Group 3 — SMM-calibrated transmission parameters.** The three free parameters
$(\gamma, \varphi, \Phi^N)$ are calibrated by minimising the sum of squared deviations
between the model-implied and empirical IRFs for output and private credit, targeting the
**non-default path only**:

$$\min_{\gamma,\,\varphi,\,\Phi^N} \sum_{h=0}^{4}\left[\left(\hat{y}^{nd,\text{model}}_h - \hat{y}^{nd,\text{LP}}_h\right)^2 + \left(\hat{\ell}^{nd,\text{model}}_h - \hat{\ell}^{nd,\text{LP}}_h\right)^2\right] \qquad\qquad (29)$$

The **default path is not used in calibration**: it is validated out-of-sample. Only
$\Delta r^{L,\text{aut}}$ is pinned to the single $h=0$ default LP moment; horizons
$h=1,\ldots,4$ under default are zero-free-parameter predictions of equation (IR.5), testing
whether the capital-depletion mechanism quantitatively matches the empirical deepening of
default recessions.

### 21. Model vs. Data Summary

| Variable / Finding | Model prediction | Empirical estimate |
|---|---|---|
| ND output at $h=0$ | $-\varepsilon_p\gamma\hat{s}_0 < 0$ | $-1.86$ pp |
| ND output trough | $h=1$–$2$ (balance-sheet deepening) | $-3.76$ pp at $h=1$ |
| ND credit/GDP at $h=4$ | $-(\lambda N^{ss}/\ell^{ss})|\hat{n}_4^{nd}|$ | $-4.78$ pp |
| DEF output at $h=0$ | $-\varepsilon_p\Delta r^{L,\text{aut}}$ (pinned) | $-3.42$ pp |
| DEF output at $h=1$ | $\alpha\ln(1-\delta) - \varepsilon_p\Delta r^{L,\text{aut}}$ (oos) | $-5.50$ pp |
| DEF credit/GDP | $\approx 0$ (GDP-denominator effect) | $\approx 0$ (not significant) |
| Doom loop (bank sov. bonds) | $+\lambda N_t$ under default | $+4$–$6$ pp (default only) |

---

## 22. Notation Summary

| Symbol | Definition |
|---|---|
| $C_t$, $D_t$, $L_t$ | Consumption, deposits, labour |
| $\beta$, $\sigma$ | Household discount factor; CRRA risk aversion |
| $Y_t$, $K_t$, $I_t$ | Output, capital stock, investment |
| $A$, $\alpha$, $\delta$ | TFP; capital share ($0.33$); depreciation ($0.10$) |
| $\xi$, $\varepsilon_p$ | Working-capital share of wage bill; working-capital output elasticity |
| $n$ | Labour exponent: $(1-\alpha)/\alpha$ |
| $\Psi(R^L_t)$ | Working-capital distortion factor: decreasing in $R^L_t$ |
| $\eta$ | Interest semi-elasticity of capital stock |
| $b^B_t$, $\ell_t$, $N_t$ | Bank sovereign bond holdings; private loans; net worth |
| $\lambda$ | Bank leverage multiplier ($= 10$) |
| $\mathrm{B}$ | Balance-sheet sensitivity: $b^{B,ss}/[N^{ss}(1+s^{ss})^2]$ |
| $\Phi^N$ | Net-worth persistence parameter (SMM-calibrated) |
| $R^*$ | International risk-free rate ($= 0.04$) |
| $s_t$ | Sovereign EMBIG spread (decimal) |
| $\gamma$ | Sovereign-ceiling pass-through parameter (SMM-calibrated) |
| $R^O_t$ | International wholesale funding rate: $R^* + \gamma s_t$ |
| $R^L_t$, $R^{L,ss}$ | Domestic lending rate; its steady-state value |
| $\varphi$, $\Omega$ | Balance-sheet pass-through; amplification coefficient $\varphi b^{B,ss}/N^{ss}$ |
| $\hat{x}_t$ | Log-deviation of variable $x_t$ from steady state: $\log(x_t/x^{ss})$ |
| $h$ | Horizon (years after crisis onset) |
| $y_t$ | Government income (Tauchen-discretised AR(1)) |
| $\rho_y$, $\sigma_\varepsilon$ | Income AR(1) persistence; innovation standard deviation |
| $B_t$ | Government face-value debt |
| $q_t$ | Bond price: endogenous, from fixed point (27) |
| $d_t$ | Default indicator: $\mathbf{1}[y_t \in \mathcal{D}(B_t)]$ |
| $\theta$ | Recovery rate on restructured sovereign debt ($= 0.62$) |
| $\mu$ | Probability of market re-entry after default ($= 0.22$) |
| $\beta_p$ | Government discount factor ($< \beta$) |
| $h_{\text{def}}$ | Arellano asymmetric default cost parameter |
| $\mathbf{P}$ | Tauchen transition matrix ($N_y \times N_y$) |
| $\Delta r^{L,\text{aut}}$ | Autarky lending-rate wedge: $R^{L,\text{aut}} - R^{L,ss}$ |
| $\pi_{t+1}$ | Default probability: $\Pr(d_{t+1}=1 \mid B_{t+1}, y_t)$ |
| $\rho_s$ | Spread AR(1) persistence (estimated from panel) |
