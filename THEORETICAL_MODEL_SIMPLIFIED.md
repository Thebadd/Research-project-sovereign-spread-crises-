# Theoretical Model — Simplified Representation
## Keeping the Nonlinear Default Block

This appendix presents a parsimonious model that retains the **fully nonlinear sovereign
default block** (endogenous default decision, bond-price fixed point, value-function iteration)
while stripping away secondary elements — the household Euler equation, government recursive
optimisation, and the wholesale-funding income share $\Phi^O$ — that are not needed to
reproduce the nine empirical findings. The model therefore has two blocks: a **nonlinear
default and spread block** (Block 1) that generates the endogenous spread path, and a
**log-linear transmission block** (Block 2) that maps that spread path into output and credit
responses. This two-block design follows Bocola (2016).

---

## Block 1 — Nonlinear Default and Spread Determination

### 1.1 Income Process

The country's endowment income $y_t$ follows a stationary AR(1) process in logs:

$$\log y_{t+1} = \rho_y \log y_t + \varepsilon_{t+1}, \qquad \varepsilon_{t+1} \sim \mathcal{N}(0, \sigma_\varepsilon^2) \qquad\qquad (1)$$

where $\rho_y \in (0,1)$ is the income persistence and $\sigma_\varepsilon > 0$ is the innovation
standard deviation. Both are estimated from the panel data by an AR(1) regression on
country-detrended log GDP per capita (see Section 3 of the paper). The process is discretised
using the Tauchen (1986) method into $N_y = 21$ states with transition matrix $\mathbf{P}$,
so the state space for income is the finite grid $\{y_1, \ldots, y_{N_y}\}$.

### 1.2 Government's Default Problem

The government issues one-period discount bonds with face value $B_{t+1}$ and price $q_t$.
It maximises the discounted utility of primary expenditure $G_t$ subject to its budget
constraint. The **state variables** are current debt $B_t$ and current income $y_t$. The
**value of repayment** $V^R(B_t, y_t)$ satisfies the Bellman equation:

$$V^R(B_t, y_t) = \max_{G_t,\; B_{t+1}} \left\{ U(G_t) + \beta_p\, \mathbb{E}_t\left[V(B_{t+1}, y_{t+1})\right] \right\} \qquad\qquad (2)$$

subject to the budget constraint:

$$G_t = y_t - B_t + q_t(B_{t+1}, y_t)\, B_{t+1} \qquad\qquad (3)$$

where $U(G_t) = G_t^{1-\sigma}/(1-\sigma)$ is the government's CRRA utility with risk-aversion
parameter $\sigma$, $\beta_p \in (0,1)$ is the **government's discount factor** (set below the
household's $\beta$ to capture political short-termism and the temptation to default), and
$V(B_{t+1}, y_{t+1}) = \max\{V^R(B_{t+1}, y_{t+1}),\, V^D(y_{t+1})\}$ is the next-period
continuation value — the government will optimally choose repayment or default next period.
The bond price $q_t(B_{t+1}, y_t)$ depends on the borrowing choice $B_{t+1}$ and current
income $y_t$ because foreign investors price in the probability of future default (see §1.4).

Upon default the government is **excluded from international capital markets** and receives
autarky income $y^{\text{def}}_t = \min(y_t,\, h_{\text{def}}\, \bar{y})$, where $h_{\text{def}} < 1$
imposes an asymmetric output cost — when income is above the threshold $h_{\text{def}}\bar{y}$,
the government bears a real cost of default; when income is already below the threshold, the
cost is zero (Arellano 2008). Re-entry to markets occurs with probability $\mu$ each period.
The **value of default** $V^D(y_t)$ satisfies:

$$V^D(y_t) = U(y^{\text{def}}_t) + \beta_p\, \mathbb{E}_t\left[\mu\, V^R(0, y_{t+1}) + (1-\mu)\, V^D(y_{t+1})\right] \qquad\qquad (4)$$

where $\mu = 0.22$ (calibrated to the average EM exclusion period of 4.5 years, Gelos et al.
2011) and re-entry occurs with debt reset to $B = 0$ after the post-haircut settlement.

### 1.3 Optimal Default Decision

The government defaults whenever the value of default exceeds the value of repayment. The
**default set** $\mathcal{D}(B_t)$ — the set of income realisations at which default is
optimal given debt level $B_t$ — is:

$$\mathcal{D}(B_t) = \left\{y_t : V^D(y_t) \geq V^R(B_t, y_t)\right\} \qquad\qquad (5)$$

The default indicator is $d_t = \mathbf{1}[y_t \in \mathcal{D}(B_t)]$: the government defaults
if and only if income falls into the default set. The default set is **increasing in $B_t$**
(higher debt makes default more attractive for any given income) and **decreasing in $y_t$**
(higher income makes repayment easier), which generates the empirical pattern of defaults
occurring during income downturns at elevated debt levels.

### 1.4 Bond Price Fixed Point

Foreign investors are risk-neutral and can invest at the world risk-free rate $R^*$. Zero
profit requires the bond price $q_t(B_{t+1}, y_t)$ to equal the expected repayment discounted
at $R^*$. If the country defaults (probability $\pi_{t+1}$), the investor recovers a fraction
$\theta \in (0,1)$ of face value (the recovery rate, calibrated to $\theta = 0.62$ following
Cruces and Trebesch 2013, who document a mean EM haircut of approximately 38%). The bond price
is therefore:

$$q_t(B_{t+1}, y_t) = \frac{1}{1+R^*}\, \mathbb{E}_t\left[(1-d_{t+1})\cdot 1 + d_{t+1}\cdot\theta\right] \qquad\qquad (6)$$

where the expectation is taken over next-period income $y_{t+1}$ given today's $y_t$. The
**sovereign spread** implied by this price is:

$$s_t(B_{t+1}, y_t) = \frac{1}{q_t(B_{t+1}, y_t)} - 1 - R^* \qquad\qquad (7)$$

This spread is **endogenous**: it rises when debt is high relative to income, because foreign
investors assign a higher probability of future default. Equations (2)–(7) form a fixed-point
problem: the bond price depends on the default policy, which depends on the value functions,
which depend on the bond price. This fixed point is solved by **value-function iteration**
(see §1.5).

### 1.5 Numerical Solution: Value-Function Iteration

The fixed point is solved on a discrete grid over $(B_t, y_t)$: a uniform debt grid with
$N_B = 200$ points on $[0, 1.5\,\bar{y}]$ and the Tauchen-discretised income grid with
$N_y = 21$ points. The algorithm iterates until convergence (sup-norm tolerance $10^{-6}$,
maximum 2000 iterations):

1. **Value of default** (equation 4): update $V^D(y_i)$ for each income state $y_i$ using
   the current value function $V$ and transition matrix $\mathbf{P}$.

2. **Value of repayment** (equation 2): for each $(B_i, y_j)$, compute the consumption matrix
   $c(B_i, B'_k, y_j) = y_j - B_i + q(B'_k, y_j)\cdot B'_k$ for all candidate debt choices
   $B'_k$, apply CRRA utility, add discounted continuation $\beta_p\, \mathbb{E}[V(B'_k, y')]$,
   and take the maximum over $B'_k$.

3. **Default policy** (equation 5): set $d(B_i, y_j) = 1$ if $V^D(y_j) \geq V^R(B_i, y_j)$,
   else 0.

4. **Bond price update** (equation 6): compute $q^{\text{new}}(B'_k, y_j)$ using the updated
   default policy and transition probabilities. Apply dampening
   $q \leftarrow 0.5\,q^{\text{new}} + 0.5\,q^{\text{old}}$ for numerical stability.

5. **Check convergence**: if $\sup|V^{\text{new}} - V^{\text{old}}| < 10^{-6}$ and
   $\sup|q^{\text{new}} - q^{\text{old}}| < 10^{-6}$, stop; else return to step 1.

### 1.6 Crisis Event Study: Model-Implied Spread Path

After solving the fixed point, the model is **simulated for 100,000 periods** (discarding
1,000 burn-in periods). A spread crisis **onset** is defined as the first period the spread
crosses the 1,000 bps threshold from below while the country is in good standing. Around each
onset, the spread path at horizons $h = -2, \ldots, 4$ is recorded and averaged across all
crisis events in the simulation. This yields the **model-implied spread event path**
$\{s^{\text{model}}_h\}_{h=-2}^{4}$, which feeds directly into Block 2 as the forcing variable.

**Model moments** compared to data targets:

| Moment | Model target | Data |
|---|---|---|
| Default frequency (annual) | $\approx 3\%$ | $2$–$4\%$ (EM sample) |
| Mean spread (bps) | $\approx 350$ bps | $s^{ss}$ from tranquil-period median |
| SD spread (bps) | $\approx 200$ bps | Panel AR(1) residuals |
| Mean debt/GDP | $\approx 0.45$ | Tranquil-period mean |

Parameters $\beta_p$ and $h_{\text{def}}$ are tuned to jointly match default frequency and mean
spread.

---

## Block 2 — Log-Linear Transmission

Block 2 takes the model-implied spread path $\{s^{\text{model}}_h\}$ from Block 1 as given
and maps it into output and credit responses through the banking and production sectors.
This block is **log-linearised** around the non-crisis steady state, so it is computationally
trivial and analytically transparent.

### 2.1 Firms and the Working-Capital Channel

A representative firm produces output $Y_t = A K_t^\alpha L_t^{1-\alpha}$ and must borrow a
fraction $\xi$ of its wage bill from domestic banks at rate $R^L_t$ before production takes
place. Optimisation yields the log-linearised output equation:

$$\hat{y}_t = \alpha\hat{k}_t - \varepsilon_p(R^L_t - R^{L,ss}) \qquad\qquad (8)$$

where $\hat{k}_t \equiv \log(K_t/K^{ss})$ is the capital stock log-deviation,
$\varepsilon_p \equiv \xi n R^{L,ss}/[1+\xi(R^{L,ss}-1)] > 0$ is the **working-capital output
elasticity**, and $n = (1-\alpha)/\alpha$. Capital accumulates as
$K_{t+1} = (1-\delta)K_t + I_t$, where the investment equation (from the firm's user-cost
condition) gives:

$$\hat{k}_{t+1} = (1-\delta)\hat{k}_t - \eta(R^L_t - R^{L,ss}) \qquad\qquad (9)$$

with $\eta \equiv 1/[(1-\alpha)(R^{L,ss}-(1-\delta))] > 0$ the **interest semi-elasticity of
capital**: a 1 pp rise in the lending rate reduces next period's capital stock by $\eta$ pp.

### 2.2 Banks, Net Worth, and the Leverage Constraint

A bank holds sovereign bonds $b^B_t$ and private loans $\ell_t$, financed by net worth $N_t$
and debt. The **leverage constraint** $b^B_t + \ell_t \leq \lambda N_t$ binds with equality
($\lambda = 10$). Net worth evolves as bank equity absorbs mark-to-market losses on sovereign
bonds whenever the spread rises:

$$\hat{n}_{t+1} = \Phi^N \hat{n}_t - \mathrm{B}\,\hat{s}_{t+1} \qquad\qquad (10)$$

where $\hat{n}_t \equiv \log(N_t/N^{ss})$, $\Phi^N \in (0,1)$ is the **net-worth persistence**
(fraction of net worth preserved through retained earnings each period), and
$\mathrm{B} \equiv b^{B,ss}/[N^{ss}(1+s^{ss})^2] > 0$ is the **balance-sheet sensitivity**
(net-worth loss per unit rise in the spread, proportional to steady-state sovereign exposure).

Private credit contracts proportionally to net worth, amplified by leverage:

$$\hat{\ell}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\,\hat{n}_t \qquad\qquad (11)$$

### 2.3 Sovereign Ceiling and the Domestic Lending Rate

No domestic borrower can obtain external funding below the sovereign rate. The **sovereign
ceiling** pins banks' international wholesale funding cost to:

$$R^O_t = R^* + \gamma s_t \qquad\qquad (12)$$

where $\gamma \in (0,1]$ is the **sovereign-ceiling pass-through**. Banks price private loans
at a mark-up that also reflects their balance-sheet leverage in sovereign bonds:

$$R^L_t - R^{L,ss} = \gamma\hat{s}_t - \Omega\hat{n}_t \qquad\qquad (13)$$

where $\Omega \equiv \varphi\, b^{B,ss}/N^{ss} > 0$ is the **balance-sheet amplification
coefficient** ($\varphi$ is the sensitivity of the lending rate to the sovereign-bond
leverage ratio). Substituting (13) into (8) and (9) gives the complete system:

$$\hat{y}_t = \alpha\hat{k}_t - \varepsilon_p\left(\gamma\hat{s}_t - \Omega\hat{n}_t\right) \qquad\qquad \text{(LL.2a)}$$

$$\hat{k}_{t+1} = (1-\delta)\hat{k}_t - \eta\left(\gamma\hat{s}_t - \Omega\hat{n}_t\right) \qquad\qquad \text{(LL.3a)}$$

$$\hat{n}_{t+1} = \Phi^N\hat{n}_t - \mathrm{B}\,\hat{s}_{t+1} \qquad\qquad \text{(LL.5a)}$$

These three equations — driven by the spread path $\{\hat{s}_h\}$ from Block 1 — constitute
the complete non-default transmission system. There are three free parameters:
$\gamma$, $\varphi$ (or equivalently $\Omega$), and $\Phi^N$.

### 2.4 Default Path Modifications

Under default, three structural breaks are imposed:

**Break 1 — Investment suspended.** Market exclusion cuts external financing, so $I_t = 0$.
Capital declines at the depreciation rate:

$$\hat{k}^{\text{def}}_h = h\ln(1-\delta) \qquad\qquad \text{(IR.5)}$$

This term is **zero at $h=0$** and **strictly more negative at each subsequent horizon**,
generating the deepening output loss observed empirically.

**Break 2 — Private lending collapses (gambling for resurrection).** When the option value
of sovereign-bond recovery $\theta$ exceeds the return from private lending, banks optimally
set $\ell^{\text{def}}_t = 0$ and $b^{B,\text{def}}_t = \lambda N_t$. The autarky lending
rate $R^{L,\text{aut}} > R^{L,ss}$ produces the impact output loss at $h=0$:

$$\hat{y}^{\text{def}}_h = \alpha\cdot h\ln(1-\delta) - \varepsilon_p\,\Delta r^{L,\text{aut}} \qquad\qquad \text{(IR.5 full)}$$

where $\Delta r^{L,\text{aut}} \equiv R^{L,\text{aut}} - R^{L,ss} > 0$ is pinned to the
single $h=0$ default LP moment (leaving $h=1,\ldots,4$ as **genuine out-of-sample
predictions**).

**Break 3 — Sovereign ceiling severed.** With $F_t = 0$, equation (12) no longer applies.
The working-capital channel in the non-default system (LL.2a) is replaced by the autarky
wedge $\Delta r^{L,\text{aut}}$.

---

## 3. The Credit/GDP Paradox

A central empirical finding is that **private credit/GDP falls more in non-default episodes**
despite private lending collapsing to zero under default ($\ell^{\text{def}}_t = 0$). The
resolution is the **GDP-denominator effect**.

The empirical outcome is the **ratio** $\ell_t / Y_t$. Under default, output collapses via
capital depletion (IR.5) — a cumulative mechanism adding $\approx 3.3$ pp of output loss per
year. By $h = 3$, GDP has fallen by $\approx 9$ pp from capital depletion alone. The GDP
denominator shrinks at least as fast as the private-credit numerator, leaving the ratio
approximately unchanged.

Under non-default, the credit channel dominates: $\hat{\ell}_t$ (equation 11) falls much
faster than $\hat{y}_t$ (equation LL.2a), so the ratio $\ell_t / Y_t$ falls substantially.

$$\underbrace{|\Delta\ell/Y|^{nd} > |\Delta\ell/Y|^{def}}_{\text{empirical finding}}
\quad \text{because} \quad
\begin{cases}
\text{ND:} & |\Delta\ell_t| \gg |\Delta Y_t| \quad \text{(credit channel dominates)} \\
\text{DEF:} & |\Delta Y_t| \gg |\Delta\ell_t| \quad \text{(capital depletion dominates)}
\end{cases}$$

**A smaller credit/GDP response in default is not evidence that banks lend more — it is
evidence that GDP collapses faster than credit, driven by a completely different mechanism.**

---

## 4. SMM Calibration

**Block 1 parameters** ($\beta_p$, $h_{\text{def}}$, $\theta$, $\mu$, $\rho_y$, $\sigma_\varepsilon$)
are calibrated from the EM literature and panel AR(1) estimates (see Section 3 of the paper).

**Block 2 parameters** ($\gamma$, $\varphi$, $\Phi^N$) are calibrated by **Simulated Method
of Moments** targeting the **non-default empirical IRFs** for output and private credit
at horizons $h = 0, \ldots, 4$:

$$\min_{\gamma,\,\varphi,\,\Phi^N} \sum_{h=0}^{4}\left[(\hat{y}^{nd,\text{model}}_h - \hat{y}^{nd,\text{LP}}_h)^2 + (\hat{\ell}^{nd,\text{model}}_h - \hat{\ell}^{nd,\text{LP}}_h)^2\right] \qquad\qquad (14)$$

The **default path is validated out-of-sample**: only the single $h=0$ default moment pins
$\Delta r^{L,\text{aut}}$; the paths at $h = 1, \ldots, 4$ are zero-free-parameter predictions
of equation (IR.5 full), providing a genuine test of the capital-depletion mechanism.

---

## 5. Summary of Mechanisms

| Variable | Non-default | Default | Dominant mechanism |
|---|---|---|---|
| **Output** | $-3.76$ pp at $h=1$ | $-5.50$ pp at $h=1$, deepening | ND: working capital + balance sheet; DEF: capital depletion (IR.5) |
| **Private credit/GDP** | $-4.78$ pp at $h=4$ | $\approx 0$ measured | ND: $\ell_t$ falls faster than $Y_t$; DEF: $Y_t$ falls faster than $\ell_t$ |
| **Investment** | Reduced ($I_t > 0$) | Suspended ($I_t = 0$) | DEF: market exclusion |
| **Sovereign bonds at banks** | $\approx 0$ | $+4$–$6$ pp | DEF: gambling for resurrection |
| **Current account** | Gradual surplus | Front-loaded surge | DEF: sudden stop ($F_{-1} + (1-\delta^B)B_0$) |
| **Govt expenditure** | Positive throughout | Boom then bust | DEF: windfall at $h=0$, conditionality at $h \geq 2$ |

---

## 6. Notation Summary

| Symbol | Definition |
|---|---|
| $y_t$ | Endowment income (Tauchen-discretised AR(1)) |
| $\rho_y$, $\sigma_\varepsilon$ | Income persistence and innovation SD |
| $B_t$ | Government debt (face value) |
| $q_t$ | Bond price (endogenous) |
| $s_t$ | Sovereign EMBIG spread: $1/q_t - 1 - R^*$ |
| $d_t$ | Default indicator: $\mathbf{1}[y_t \in \mathcal{D}(B_t)]$ |
| $\theta$ | Recovery rate on restructured debt ($= 0.62$) |
| $\mu$ | Probability of re-entry after default ($= 0.22$) |
| $\beta_p$ | Government discount factor |
| $h_{\text{def}}$ | Asymmetric default cost parameter |
| $Y_t$, $K_t$, $L_t$ | Output, capital, labour |
| $A$, $\alpha$, $\delta$ | TFP, capital share ($0.33$), depreciation ($0.10$) |
| $\xi$, $\varepsilon_p$ | Working-capital share; working-capital output elasticity |
| $\eta$ | Interest semi-elasticity of capital |
| $b^B_t$, $\ell_t$, $N_t$ | Bank sovereign bonds, private loans, net worth |
| $\lambda$ | Bank leverage ($= 10$) |
| $\mathrm{B}$ | Balance-sheet sensitivity of net worth to spread |
| $\Phi^N$ | Net-worth persistence parameter |
| $R^*$ | International risk-free rate ($= 0.04$) |
| $\gamma$ | Sovereign-ceiling pass-through |
| $\varphi$, $\Omega$ | Balance-sheet pass-through; amplification coefficient |
| $\Delta r^{L,\text{aut}}$ | Autarky lending-rate wedge |
| $\hat{x}_t$ | Log-deviation of $x_t$ from steady state |
