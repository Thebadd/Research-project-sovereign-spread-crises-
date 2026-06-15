# A Stylized Model of Sovereign Spread Crises With and Without Default

**Theoretical Appendix** — *Companion to "The Aftermath of Sovereign Spread Crises With and Without Default"*

---

## Overview

This appendix develops a small open economy model that rationalizes the central empirical findings of the paper: (i) sovereign spread crises generate large and persistent output losses *even when no default occurs*; (ii) default-linked crises are roughly twice as costly and front-loaded; (iii) the private-credit and investment contractions are concentrated in *non-default* episodes; (iv) bank accumulation of sovereign debt (the doom loop) is exclusive to *default-linked* episodes; and (v) the current-account adjustment is larger and earlier under default.

The model's defining feature is a **separation between the international sovereign bond market** — where the government issues dollar-denominated external debt priced by foreign investors — **and the domestic credit market** — where banks lend to firms in local currency. Transmission from the former to the latter runs through two channels: a **wholesale funding-cost channel** (the sovereign ceiling, by which international lenders cap domestic bank borrowing costs at the sovereign's own rate) and a **balance-sheet channel** (sovereign bond price declines erode bank net worth and tighten lending capacity). This dual structure delivers both a non-default transmission mechanism and a distinct default mechanism operating through fundamentally different channels.

The solution uses a **two-block decomposition** following Bocola (2016): a nonlinear sovereign-risk block that determines the endogenous default decision and the equilibrium spread (solved by value-function iteration in `do/15_solve_default.do`), and a log-linear transmission block that maps the spread into output, credit, investment, and the current account (solved analytically in closed form in `do/16_model_irf.do`). Nonlinearity is retained where it matters — the default choice — and linearized where it does not.

---

## 1. Environment and Agents

Consider a small open economy operating in discrete time $t = 0, 1, 2, \dots$ Five types of agents interact: a representative household, a representative firm, a domestic banking sector, a government, and a continuum of competitive foreign investors.

---

### 1.1 Households

The representative household maximizes lifetime expected utility over consumption:

$$\max \; \mathbb{E}_0 \sum_{t=0}^{\infty} \beta^t \, \frac{C_t^{1-\sigma}}{1-\sigma} \qquad\qquad (1)$$

where:
- $C_t > 0$ is real consumption in period $t$
- $\beta \in (0,1)$ is the subjective discount factor (how much the household values future consumption relative to today; $\beta = 0.96$ means one unit of future utility is worth 0.96 today)
- $\sigma > 0$ is the inverse of the intertemporal elasticity of substitution (IES); higher $\sigma$ means stronger preference for smooth consumption across time; $\sigma = 2$ is standard for emerging markets
- $\mathbb{E}_0[\cdot]$ denotes the mathematical expectation conditional on information available at $t=0$

The household's period budget constraint is:

$$C_t + D_{t+1} = R^D_t \, D_t + w_t + \Pi_t \qquad\qquad (2)$$

where:
- $D_t \ge 0$ is the stock of domestic bank deposits held at the beginning of period $t$
- $D_{t+1}$ is new deposits saved for next period
- $R^D_t \ge 1$ is the gross return on deposits (the deposit interest rate; $R^D_t = 1.05$ means 5% annual interest)
- $w_t > 0$ is the equilibrium wage rate per unit of labor
- $\Pi_t$ are lump-sum profits from bank ownership, received each period

The household supplies labor inelastically ($L_t = 1$ in all periods) and holds only domestic deposits — it has **no direct access to international capital markets**, consistent with the financial constraints typical of emerging market households.

Maximizing (1) subject to (2) yields the **consumption Euler equation**:

$$C_t^{-\sigma} = \beta \, R^D_{t+1} \, \mathbb{E}_t\!\left[C_{t+1}^{-\sigma}\right] \qquad\qquad (3)$$

This equation says the household equates the marginal utility of consuming today ($C_t^{-\sigma}$) to the discounted expected marginal utility of consuming tomorrow, scaled by the deposit return $R^D_{t+1}$. A higher deposit rate $R^D_{t+1}$ induces the household to save more and consume less today — the **intertemporal substitution effect** that is central to the model's transmission mechanism.

---

### 1.2 Firms

A representative competitive firm produces output using capital and labor with a Cobb-Douglas technology:

$$Y_t = A \, K_t^{\alpha} \, L_t^{1-\alpha}, \qquad \alpha \in (0,1) \qquad\qquad (4)$$

where:
- $Y_t > 0$ is real output (GDP) in period $t$
- $A > 0$ is total factor productivity (TFP), assumed constant
- $K_t > 0$ is the stock of physical capital (machinery, buildings, equipment) available at the beginning of period $t$
- $L_t = 1$ is inelastic labor supply (normalized to one)
- $\alpha \in (0,1)$ is the capital share of income (calibrated at $\alpha = 0.33$, so capital earns one-third of output)

Capital evolves over time according to:

$$K_{t+1} = (1-\delta) K_t + I_t \qquad\qquad (5)$$

where:
- $K_{t+1}$ is next period's capital stock
- $\delta \in (0,1)$ is the annual physical depreciation rate (calibrated at $\delta = 0.10$, meaning 10% of capital wears out each year)
- $I_t \ge 0$ is gross investment (new capital installed in period $t$)

**Working-capital constraint.** The firm must pre-finance a fraction $\xi \in (0,1)$ of its wage bill with short-term bank loans *before* production takes place. This reflects the reality that firms in emerging markets need credit to pay workers before revenues arrive. The firm also borrows from banks to finance new investment. Total credit demand from the firm is therefore:

$$\ell_t = \xi \, w_t + I_t \qquad\qquad (6)$$

where:
- $\ell_t \ge 0$ is total short-term bank lending to the firm
- $\xi \in (0,1)$ is the working-capital share (the fraction of the wage bill that must be pre-financed; calibrated by SMM)
- $\xi w_t$ is the working-capital loan (to pre-finance wages)
- $I_t$ is the investment loan

The working-capital cost raises the effective cost of labor from $w_t$ to $[1 + \xi(R^L_t - 1)]w_t$, where $R^L_t \ge 1$ is the domestic lending rate charged by banks. Firm profit maximization yields the labor demand condition:

$$(1-\alpha) A K_t^{\alpha} L_t^{-\alpha} = [1 + \xi(R^L_t - 1)] \, w_t \qquad\qquad (7)$$

The left side is the marginal product of labor; the right side is the effective wage cost inclusive of pre-financing. Combining (4) and (7), equilibrium output can be written purely as a function of the lending rate and the capital stock:

$$Y_t = A K_t^{\alpha} \cdot \Psi(R^L_t) \qquad\qquad (8)$$

where the function $\Psi$ captures the effect of the lending rate on effective labor input:

$$\Psi(R^L_t) \equiv \left[\frac{1-\alpha}{1+\xi(R^L_t-1)}\right]^{n}, \qquad n = \frac{1-\alpha}{\alpha}$$

- $\Psi(R^L_t) > 0$ is a decreasing function of $R^L_t$: higher lending rates raise the cost of pre-financing wages, reduce effective labor demand, and lower output
- $n = (1-\alpha)/\alpha$ is the ratio of the labor share to the capital share (approximately $2$ when $\alpha = 0.33$)
- $\Psi'(R^L_t) < 0$ for all $\xi > 0$, confirming that **output is strictly decreasing in the lending rate**

**This is the key result of the firm block**: any upstream shock that raises $R^L_t$ — whether through the sovereign ceiling, through bank balance-sheet losses, or through market exclusion — generates an *immediate within-period output contraction*, independently of and in addition to the longer-run contraction from reduced capital accumulation.

The firm equates the marginal product of capital to the user cost of capital to determine investment:

$$\alpha A K_{t+1}^{\alpha-1} = R^L_t - (1-\delta) \qquad\qquad (9)$$

where the left side is the marginal product of capital next period, and the right side is the user cost (the lending rate minus the undepreciated value of capital). A rise in $R^L_t$ raises the user cost, reduces investment demand, and slows capital accumulation — the **investment channel** that generates medium-run output losses.

---

### 1.3 Banking Sector with Dual Funding Structure

The banking sector is the central transmission mechanism in the model. Banks differ from standard frameworks in two respects: they hold sovereign bonds priced in international markets, exposing their balance sheets to foreign investor sentiment; and they fund themselves through a dual structure of domestic deposits and international wholesale funding, creating a direct channel from the sovereign spread to domestic credit conditions.

#### 1.3.1 Balance Sheet

The bank's balance sheet identity is:

$$b^B_t + \ell_t = D_t + F_t + N_t \qquad\qquad (10)$$

where:
- $b^B_t \ge 0$ is the bank's holdings of dollar-denominated sovereign bonds, valued at the international market price $q_t$ (so the market value of the portfolio is $q_t \cdot b^B_t$)
- $\ell_t \ge 0$ is total domestic lending to firms (working capital plus investment, from equation 6)
- $D_t \ge 0$ is deposits received from households — a domestic liability
- $F_t \ge 0$ is international wholesale funding — funds borrowed from foreign interbank markets in dollars
- $N_t \ge 0$ is bank net worth (equity capital, the residual claim of bank shareholders)

The left side is total bank assets (sovereign bonds plus domestic loans); the right side is total liabilities plus equity.

#### 1.3.2 International Wholesale Funding and the Sovereign Ceiling

Banks can borrow in international interbank markets to supplement domestic deposit funding. Foreign lenders charge domestic banks a gross funding rate that reflects both the international risk-free rate $R^*$ and the sovereign's credit risk. This is the **sovereign ceiling mechanism**: because a sovereign default would disrupt the legal and payments infrastructure that cross-border contracts depend on, foreign lenders will not lend to domestic banks at a rate *below* the sovereign's own cost of external borrowing:

$$R^O_t = R^* + \gamma s_t, \qquad \gamma \in (0,1] \qquad\qquad (11)$$

where:
- $R^O_t \ge 1$ is the gross international wholesale funding rate charged to domestic banks
- $R^* > 0$ is the international risk-free interest rate (the rate on US Treasury bills, approximately 4% annually)
- $s_t \ge 0$ is the sovereign EMBIG spread — the excess yield demanded by international investors over the risk-free rate to hold the government's external bonds, expressed in decimal form (e.g., $s_t = 0.10$ means 1000 basis points)
- $\gamma \in (0,1]$ is the **sovereign ceiling pass-through parameter**: when $\gamma = 1$ the ceiling binds fully and domestic banks cannot borrow internationally at a rate below the sovereign; when $\gamma < 1$ banks may partially escape sovereign risk through collateral or institutional standing

At the interior equilibrium where banks are indifferent between domestic and international funding at the margin, arbitrage between the two sources forces:

$$R^D_t = R^O_t = R^* + \gamma s_t \qquad\qquad (12)$$

where $R^D_t$ is the domestic deposit rate that households receive.

**Equation (12) is the sovereign ceiling's key implication**: the domestic deposit rate is pinned to the international sovereign spread. Any rise in $s_t$ immediately raises $R^D_{t+1}$, which via the Euler equation (3) induces households to defer consumption and reduce demand. Substituting (12) into (3):

$$C_t^{-\sigma} = \beta \, (R^* + \gamma s_{t+1}) \, \mathbb{E}_t\!\left[C_{t+1}^{-\sigma}\right] \qquad\qquad (3')$$

This modified Euler equation makes explicit that **a sovereign spread spike reduces consumption today** through the intertemporal substitution channel, even before any change in credit quantities or investment. This channel is absent from models that treat the spread as only a government borrowing cost.

#### 1.3.3 Bank Net Worth and the Balance-Sheet Channel

Bank net worth evolves as the return on assets minus the cost of liabilities:

$$N_{t+1} = R^L_t \ell_t + (R^B_t + \Delta q_{t+1}) b^B_t - R^D_t D_t - R^O_t F_t \qquad\qquad (13)$$

where:
- $R^L_t$ is the gross domestic lending rate charged to firms
- $R^B_t$ is the coupon yield on sovereign bonds
- $\Delta q_{t+1} = q_{t+1} - q_t$ is the **capital gain or loss** on the sovereign bond portfolio between period $t$ and $t+1$: if the bond price falls, the bank suffers a mark-to-market loss on its sovereign holdings

Using the bond pricing relationship $q_t = 1/(1+s_t)$, a spread increase $\Delta s_t > 0$ generates an immediate capital loss:

$$\Delta q_t \approx -\frac{\Delta s_t}{(1+s_t)^2} < 0 \qquad\qquad (14)$$

The balance-sheet shock to bank net worth upon a spread increase is therefore:

$$\Delta N_t = b^B_{t-1} \cdot \Delta q_t = -\frac{b^B_{t-1} \, \Delta s_t}{(1+s_t)^2} < 0 \qquad\qquad (15)$$

where $b^B_{t-1}$ is the **pre-existing stock of sovereign bonds** held by the bank entering the crisis. Two important properties follow from (15):

1. The balance-sheet loss is **proportional to pre-existing sovereign exposure** $b^B_{t-1}$: banks that entered the crisis heavily invested in sovereign bonds suffer larger losses for the same spread shock.
2. The balance-sheet channel and the sovereign-ceiling channel (11)–(12) are **additive and simultaneous**: the ceiling raises the cost of new funding while the balance-sheet loss reduces the capacity to fund existing assets.

#### 1.3.4 Leverage Constraint and the Domestic Lending Rate

Banks face a regulatory or market-imposed leverage constraint limiting total assets to a multiple of net worth:

$$b^B_t + \ell_t \le \lambda N_t, \qquad \lambda > 1 \qquad\qquad (16)$$

where $\lambda > 1$ is the **maximum leverage ratio** (assets per unit of equity; calibrated at $\lambda = 10$, consistent with EM banking sector data). When $N_t$ falls through the balance-sheet channel (15), the binding constraint (16) forces a contraction in total assets — either sovereign bond holdings, or domestic lending, or both.

The domestic lending rate charged to firms reflects both the bank's weighted-average funding cost and a risk premium that increases with the degree of balance-sheet impairment:

$$R^L_t = \omega R^D_t + (1-\omega) R^O_t + \varphi \, \frac{b^B_{t-1}}{N_t} \qquad\qquad (17)$$

where:
- $\omega \in (0,1)$ is the **share of domestic deposits** in total bank funding (the remainder $1-\omega$ is international wholesale funding)
- $\varphi > 0$ is the **balance-sheet pass-through parameter**: how strongly the ratio of sovereign bond holdings to net worth amplifies the lending rate above the bank's average funding cost; calibrated by SMM
- $b^B_{t-1}/N_t$ is the **sovereign-bond-to-net-worth ratio**: this rises when the spread increases (because $N_t$ falls via equation 15), further raising lending costs

Substituting (11) and (12) into (17):

$$\boxed{R^L_t = R^* + \gamma s_t + \varphi \, \frac{b^B_{t-1}}{N_t}} \qquad\qquad (17')$$

**Equation (17') is the central transmission equation of the model.** It decomposes the domestic lending rate into three components:

1. **Baseline** $R^*$: the international risk-free rate, the floor cost of funds in the absence of sovereign risk
2. **Sovereign-ceiling term** $\gamma s_t$: the EMBIG spread transmitted to all domestic borrowers equally through the bank's international funding cost — this operates contemporaneously for all banks regardless of their sovereign bond holdings
3. **Balance-sheet amplification** $\varphi(b^B_{t-1}/N_t)$: an additional spread that grows as net worth erodes, amplifying the basic ceiling effect in proportion to pre-existing sovereign exposure — this is heterogeneous across banks and deepens over time as net worth falls

> **Remark on identification.** Equation (17') reveals why the empirical finding that the doom-loop channel is *absorbed by the lagged stock of bank claims on government* is structurally meaningful. The sovereign-ceiling term $\gamma s_t$ affects all banks equally and is fully captured by year fixed effects in the local projections. The balance-sheet amplification term $\varphi(b^B_{t-1}/N_t)$ is proportional to the *lagged* sovereign bond position $b^B_{t-1}$: banks that entered the crisis with higher sovereign exposure suffer larger balance-sheet losses for a given spread shock. Including $L.b^B$ (lagged bank claims on government) as a control in the empirical specification absorbs this term, isolating crisis-period *new accumulation* from pre-existing balance-sheet exposure. The empirical finding that controlling for $L.b^B$ renders the claims-on-government coefficient insignificant is therefore structurally predicted by the model — it is not a spurious result.

---

## 2. International Sovereign Bond Market and the Endogenous Spread

### 2.1 Foreign Investor Zero-Profit Condition

A continuum of risk-neutral foreign investors can purchase the government's external bonds at price $q_t$. They do not purchase domestic bank liabilities or firm equity — they operate exclusively in the sovereign bond market. Their participation requires that expected returns equal the risk-free rate $R^*$:

$$q_t (1+R^*) = \mathbb{E}_t\!\left[(1-\pi_{t+1}) \cdot 1 + \pi_{t+1} \cdot \theta\right] \qquad\qquad (18)$$

where:
- $q_t \in (0,1)$ is the current price of a sovereign bond with face value of 1 (a discount bond)
- $R^* > 0$ is the gross risk-free rate (see above)
- $\pi_{t+1} \in [0,1]$ is the **rational expectation** of the probability that the government defaults in period $t+1$, formed using all information available at $t$
- $\theta \in [0,1)$ is the **recovery rate**: the fraction of face value that bondholders recover in the event of default/restructuring (calibrated at $\theta = 0.62$, based on Cruces and Trebesch 2013's database of EM haircuts — a 38% average haircut implies $\theta = 0.62$)
- The right side of (18) is the expected payoff: the bond pays 1 in full with probability $1-\pi_{t+1}$, and pays $\theta < 1$ (the recovery value) with probability $\pi_{t+1}$

### 2.2 Equilibrium Spread

Using the bond pricing relationship $q_t = 1/(1+s_t)$ where $s_t$ is the spread, rearranging (18) yields the **equilibrium EMBIG spread**:

$$s_t = \frac{(1+R^*)(1-\theta)\,\pi_{t+1}}{1 - \pi_{t+1}(1-\theta)} \qquad\qquad (19)$$

The spread $s_t$ is:
- **Strictly increasing** in the default probability $\pi_{t+1}$: the higher the perceived probability of default, the larger the risk premium demanded by investors
- **Strictly decreasing** in the recovery rate $\theta$: the more investors expect to recover in case of default, the smaller the required spread
- For small default probabilities, the approximation $s_t \approx (1+R^*)(1-\theta)\pi_{t+1}$ holds

**Important:** $s_t$ here is exclusively the **sovereign's external bond spread** — the excess yield on dollar-denominated government bonds relative to US Treasuries. It is *not* the domestic lending-deposit spread or a domestic interbank rate. The transmission of $s_t$ to domestic credit conditions operates entirely through the banking sector via equations (11)–(12) and (15)–(17').

### 2.3 Spread Crisis Definition

A spread crisis episode is defined as any period in which the EMBIG spread crosses the threshold $\bar{s} = 1000$ basis points (10 percentage points) from below:

$$\Omega_t = \mathbf{1}\!\left[s_t \ge \bar{s}\right] \qquad\qquad (20)$$

where $\Omega_t = 1$ denotes a spread crisis onset and $\mathbf{1}[\cdot]$ is the indicator function. This threshold is consistent with the Pescatori-Sy (2007) criterion used in the empirical analysis. By inverting (19), $s_t \ge \bar{s}$ if and only if the default probability exceeds a threshold $\bar{\pi}$ derived from $\bar{s}$ and the structural parameters.

---

## 3. Government and the Default Decision

### 3.1 Budget Constraint

Each period the government must finance primary expenditure and repay maturing debt, either from tax revenue or by issuing new bonds:

$$G_t + (1-\delta^B) B_t = \tau Y_t + q_t \!\left[B_{t+1} - (1-\delta^B)B_t\right] \qquad\qquad (21)$$

where:
- $G_t \ge \bar{G}$ is primary government expenditure, subject to a **minimum spending floor** $\bar{G} > 0$ (essential services the government cannot cut below)
- $B_t > 0$ is the **face value of outstanding external debt** (the total amount the government owes to foreign bondholders)
- $\delta^B \in (0,1)$ is the fraction of total debt that matures each period (calibrated at $\delta^B = 0.22$, corresponding to an average debt maturity of about 4.5 years, consistent with EM sovereign debt)
- $(1-\delta^B)B_t$ is therefore the **amount of debt coming due** this period (principal repayment)
- $\tau \in (0,1)$ is the **proportional tax rate** applied to output; $\tau Y_t$ is total government revenue
- $q_t$ is the international bond price from (19)
- $B_{t+1}$ is new debt issued this period at price $q_t$, so $q_t[B_{t+1}-(1-\delta^B)B_t]$ is the **net proceeds from bond issuance** (positive when the government rolls over more than it repays)
- The borrowing limit is $B_{t+1} \le \bar{B}(s_t)$, which tightens as $s_t$ rises, reflecting the difficulty of placing new bonds when spreads are elevated

### 3.2 Recursive Default Problem

The government maximizes the discounted value of primary expenditure utility subject to the budget constraint. Let the **state variables** be $(B_t, K_t)$ — the current debt level and the current capital stock (which determines tax capacity $\tau Y_t$ through production). The **value of repayment** satisfies:

$$V^R_t(B_t,K_t) = \max_{G_t,\; B_{t+1}}\left\{ U(G_t) + \beta_p \, \mathbb{E}_t\!\left[V_{t+1}(B_{t+1},K_{t+1})\right] \right\} \qquad\qquad (22)$$

where:
- $U(G_t)$ is the government's per-period utility from primary expenditure (strictly increasing and concave)
- $\beta_p \in (0,1)$ is the **government's discount factor**, which may differ from the household's $\beta$ — a lower $\beta_p$ captures political short-termism and a stronger temptation to default
- $V_{t+1}(B_{t+1},K_{t+1}) = \max\{V^R_{t+1}, V^D_{t+1}\}$ is the continuation value — the government will optimally choose repayment or default next period
- Maximization is subject to: the budget constraint (21), the capital accumulation equation (5), the spending floor $G_t \ge \bar{G}$, and the borrowing limit $B_{t+1} \le \bar{B}(s_t)$

Upon default, the government is **excluded from international capital markets** with probability $1-\mu$ each period (it re-enters with probability $\mu$, calibrated at $\mu = 0.22$ based on Gelos et al. 2011 who find an average exclusion period of approximately 4.5 years). The **value of default** is:

$$V^D_t(K_t) = \max_{G_t}\left\{ U(G_t) + \beta_p \, \mathbb{E}_t\!\left[\mu \, V^R_{t+1}(0, K_{t+1}) + (1-\mu) V^D_{t+1}(K_{t+1})\right] \right\} \qquad\qquad (23)$$

where:
- $V^D$ depends only on $K_t$ (not $B_t$) because debt obligations are suspended during default
- The continuation value combines: with probability $\mu$, the country **re-enters markets** with debt reset to $B = 0$ after the post-haircut settlement; with probability $1-\mu$, it **remains excluded** for another period
- The autarky budget constraint is $G_t = (1-\phi)\tau Y_t$, where $\phi \in (0,1)$ is the **direct output cost of market exclusion** (capturing trade finance disruption, import-input collapse, and loss of foreign direct investment)
- Investment is zero during exclusion ($I_t = 0$) because firms cannot access international funding for capital expenditure; capital therefore decumulates at rate $\delta$

### 3.3 Default Decision and Equilibrium

The government defaults when:

$$V^D_t(K_t) \ge V^R_t(B_t, K_t) \qquad\qquad (24)$$

This inequality defines a **default set** characterized by a threshold $B^*_t(K_t)$: the government defaults if and only if current debt exceeds the threshold. The threshold is **decreasing in $s_t$** (higher spreads make repayment more costly, expanding the default region) and **decreasing in $K_t$** (lower capital reduces tax revenues, making repayment harder). Foreign investors form rational expectations of the default probability:

$$\pi_{t+1} = \Pr\!\left(B_{t+1} \ge B^*_{t+1}(K_{t+1}) \mid \mathcal{F}_t\right) \qquad\qquad (25)$$

where $\mathcal{F}_t$ is the information set at time $t$ (all variables observed up to and including period $t$).

Equations (19) and (25) form the **fixed-point system** that closes the model: the EMBIG spread reflects rational beliefs about default, and the default decision reflects the borrowing costs and output consequences generated by that spread through (11)–(17') and (8). This is solved numerically in `do/15_solve_default.do` using value-function iteration.

---

## 4. The Non-Default Path: Transmission Mechanisms

Suppose the government repays its debt in all periods ($B_t < B^*_t(K_t)$ for all $t$). A spread crisis occurs at $t=0$ because foreign investors revise $\pi_{t+1}$ upward — perhaps due to a global risk-off shock, a deterioration in fiscal fundamentals, or contagion from another EM country — driving $s_0$ past $\bar{s} = 1000$ bps. Four mechanisms activate sequentially.

---

**Mechanism 1 — Impact: Sovereign Ceiling on Consumption ($h=0$)**

The spread spike immediately raises the bank's international wholesale funding cost through (11): $R^O_0 = R^* + \gamma s_0$. By the arbitrage condition (12), the domestic deposit rate also rises: $R^D_0 = R^* + \gamma s_0$. Substituting into the household Euler equation (3'):

The rise in $R^D_0$ induces the household to defer consumption — future consumption is now more attractive relative to today because deposits pay more. This generates a **contemporaneous demand-side output reduction** through goods-market clearing (39), even before any bank credit quantity effect operates.

---

**Mechanism 2 — Impact: Sovereign Ceiling on Working Capital ($h=0$)**

Simultaneously, the domestic lending rate $R^L_0$ rises at onset through (17'): $R^L_0 = R^* + \gamma s_0 + \varphi(b^B_{-1}/N_0)$. From the output equation (8), the immediate output contraction is:

$$\frac{dY_0}{dR^L_0} = A K_0^{\alpha} \, \Psi'(R^L_0) < 0, \qquad \Psi'(R^L_0) = -\frac{\xi n \, \Psi(R^L_0)}{1+\xi(R^L_0-1)} \qquad\qquad (26)$$

where $\Psi'(R^L_0)$ is the derivative of the $\Psi$ function with respect to $R^L_0$, which is always negative. Higher lending rates raise the effective cost of pre-financing wages, reduce equilibrium labor demand and employment, and compress output within the same period. This **supply-side channel** combines with the demand-side consumption contraction above to generate the impact output loss of $-1.86$ pp at $h=0$.

---

**Mechanism 3 — Progressive: Credit Contraction and Investment ($h=1$ to $h=4$)**

As the spread remains elevated beyond $t=0$, bank net worth erodes progressively. From equation (13), each period brings a fresh mark-to-market loss on the sovereign bond portfolio (via equation 15) and compressed net interest margins (higher funding costs from equation 11 not fully passed through to lending rates). The leverage constraint (16) therefore forces a progressive contraction in total domestic lending:

$$\ell_t = \lambda N_t - b^B_t \qquad\qquad (27)$$

where $b^B_t \approx b^B_{-1}$ in the non-default path (no doom-loop incentive absent gambling-for-resurrection motives). As $N_t$ falls, $\ell_t$ falls proportionally by the factor $\lambda$. Private credit contracts gradually: empirically $-2.25$ pp at $h=2$, deepening to $-3.48$ pp at $h=4$. The rising lending rate (17') simultaneously raises the user cost (equation 9), deterring investment and slowing capital accumulation — the channel generating persistent output losses at $h=1$ to $h=4$.

---

**Mechanism 4 — Fiscal Response and Current Account ($h=0$ to $h=4$)**

The government initially finances expenditure by rolling over external bonds, exploiting residual borrowing capacity $B_{t+1} < B^*_t(K_t)$. As the spread persists and the debt stock grows, the borrowing limit tightens and expenditure gradually declines toward the floor $\bar{G}$. The **current account adjusts partially**: international wholesale funding to banks contracts as $R^O_t$ rises, reducing capital inflows, but the government continues rolling its bonds so the adjustment is gradual rather than abrupt.

Log-linearizing the goods-market clearing condition (39) around the non-crisis steady state, with $s_C$, $s_I$, $s_G$, $s_{CA}$ denoting the steady-state shares of consumption, investment, government expenditure, and the current account in output:

$$\hat{y}_t = s_C \hat{c}_t + s_I \iota_t + s_G \breve{g}_t + s_{CA} \, ca_t \qquad\qquad (28)$$

where $\hat{y}_t = \log(Y_t/Y^{ss})$, $\hat{c}_t = \log(C_t/C^{ss})$, $\iota_t = I_t/Y^{ss}$, $\breve{g}_t = G_t/Y^{ss}$, and $ca_t = CA_t/Y^{ss}$ are the respective log-deviations or normalized levels relative to the non-crisis steady state. Under the non-default path: $\hat{c}_0 < 0$ (sovereign ceiling consumption effect), $\iota_0 < 0$ (working capital and user cost effects), $\breve{g}_0 > 0$ (debt-financed fiscal expansion), and $ca_0 > 0$ partially (declining wholesale inflows). The aggregate output loss therefore begins on impact and deepens progressively as the capital stock declines.

---

## 5. The Default-Linked Path: Transmission Mechanisms

Suppose the government defaults at $t=0$: $B_0 \ge B^*_0(K_0)$. Default severs the sovereign-ceiling link entirely — with no bond issuance in international markets, the spread $s_t$ is no longer defined for the defaulting sovereign, and international wholesale funding to domestic banks collapses. Three mechanisms operate simultaneously and jointly generate an output collapse that is larger, more front-loaded, and operates through fundamentally different channels than the non-default path.

---

### 5.1 Sudden Stop: Collapse of International Wholesale Funding

Upon default, foreign interbank lenders immediately withdraw all wholesale funding from domestic banks ($F_t = 0$), since contractual enforceability across borders collapses and credit risk is unquantifiable. The bank's balance sheet identity (10) under autarky becomes:

$$b^B_t + \ell_t = D_t + N_t, \qquad F_t = 0 \qquad\qquad (29)$$

The collapse of $F_t$ from its pre-crisis level $F_{-1} > 0$ to zero constitutes a **sudden stop in international capital inflows to the banking sector**. From the current-account identity — the current account equals output minus domestic absorption:

$$CA_t = Y_t - C_t - I_t - G_t = F_{-1} + (1-\delta^B) B_0 \qquad\qquad (30)$$

where:
- $F_{-1}$ is the stock of international wholesale funding that banks lose when excluded from interbank markets
- $(1-\delta^B)B_0$ is the debt service (principal and interest) that the government *stops paying* upon default

Both components improve the current account simultaneously and immediately: the banking sector no longer attracts foreign wholesale inflows (a capital-account improvement), and the government suspends external debt service (effectively a forced current-account improvement). This **two-component sudden stop** is larger and more abrupt than the gradual adjustment under the non-default path — matching the empirical estimates of $+1.44$ pp at $h=0$ and $+2.61$ pp at $h=1$.

With $F_t = 0$, domestic banks can no longer arbitrage between domestic and international funding, severing the sovereign-ceiling transmission channel (12). The domestic deposit rate is now determined by domestic deposit-market equilibrium under autarky (see §5.4 below):

$$C_0^{-\sigma} = \beta \, R^{D,\text{aut}}_0 \, \mathbb{E}_0\!\left[C_1^{-\sigma}\right], \qquad R^{D,\text{aut}}_0 > R^{D,ss} \qquad\qquad (31)$$

where $R^{D,\text{aut}}_0$ is the autarky deposit rate, which *exceeds* the non-crisis steady-state rate because domestic savings are scarce under financial autarky. This induces a **large front-loaded consumption contraction** that contributes to the immediate output collapse.

---

### 5.2 Doom Loop: Gambling for Resurrection

With international wholesale funding eliminated ($F_t = 0$), the bank's asset portfolio is constrained to $\lambda N_t$ (from the leverage constraint 16) and must be funded entirely from domestic deposits and net worth. The government simultaneously turns to **domestic bond issuance** to finance expenditure, offering high-yield domestic sovereign bonds backed by implicit government guarantees on the banking sector.

An undercapitalized bank — whose net worth $N_0$ has been severely eroded by the bond price collapse (equation 15) — faces a portfolio choice under the leverage constraint between private loans and domestic sovereign bonds:

$$\ell_t + b^B_t \le \lambda N_t, \qquad \ell_t \ge 0, \qquad b^B_t \ge 0 \qquad\qquad (32)$$

Under autarky, private firm loans carry a high probability of non-performance $p^{RDL} \in (0,1)$ (firms cannot import inputs or access trade finance), so the expected return from private lending is $R^L_t(1-p^{RDL}) \approx 0$. The expected return from domestic sovereign bonds (priced at a discount reflecting restructuring uncertainty, with coupon rate $R^{BDK}_t$) is:

$$\mathbb{E}_t[r^{BDK}] = \mu \, \theta \, R^{BDK}_t + (1-\mu) \cdot 0 = \mu \theta R^{BDK}_t$$

where the second term reflects that if the country remains excluded next period (probability $1-\mu$), bonds pay nothing; if the country re-enters (probability $\mu$), the haircut settlement pays fraction $\theta$ of face value. The bank optimally sets $\ell_t = 0$ and saturates the leverage constraint with sovereign bonds when:

$$\mu \theta R^{BDK}_t + \gamma^B N_t \ge R^L_t (1 - p^{RDL}) \qquad\qquad (33)$$

where $\gamma^B > 0$ is the implicit government guarantee subsidy per unit of net worth. When (33) holds, the net change in sovereign bond holdings is:

$$\Delta b^B_0 = \lambda N_0 - b^B_{-1} > 0 \qquad \text{when} \qquad \lambda N_0 > b^B_{-1} \qquad\qquad (34)$$

generating the empirically observed $+4$ to $+6$ pp accumulation in bank sovereign bond holdings at $h=0$ to $h=2$ in default-linked episodes.

**The critical implication:** since private lending collapses to zero ($\ell_t = 0$), **the private-credit channel to firms is absent under the default path** — not because banks are deleveraging (gross assets remain at $\lambda N_t$), but because they are redirecting their constrained balance-sheet capacity entirely toward sovereign bonds. Firms cannot access working-capital finance, the wage bill cannot be pre-funded, and the working-capital output contraction (equation 8) is maximized.

---

### 5.3 Fiscal Boom-Bust

Upon default, the government suspends external debt service. The budget constraint under autarky (no new external bond issuance, $q_t B_{t+1} = 0$) becomes:

$$G_0 = (1-\phi)\tau Y_0 + (1-\delta^B) B_0 \qquad\qquad (35)$$

where $(1-\phi)\tau Y_0$ is tax revenue net of the direct output cost of market exclusion $\phi$, and $(1-\delta^B)B_0$ is the **debt-service suspension windfall** — the principal and interest that the government no longer has to pay. This permits a **large initial fiscal expansion**, matching the empirical estimate of $+2.19$ pp at $h=0$.

Post-restructuring, the government re-enters international markets with post-haircut debt $B^\circ = \theta B_0$, subject to IMF program conditionality or creditor-imposed conditions requiring a primary surplus $\tau Y_t - G_t \ge \rho B^\circ$ (where $\rho > 0$ is the required debt service ratio on restructured bonds). This forces **sharp expenditure retrenchment** — the "austerity phase" — generating the empirically observed boom-bust pattern: $-1.39$ pp at $h=2$ and $-1.69$ pp at $h=4$.

---

### 5.4 The Autarky Deposit Rate

Under default ($F_t = 0$, post-restructuring $b^B = 0$), the bank funding constraint becomes:

$$\ell^{\text{aut}}_t = D^{\text{aut}}_t + N^{\text{aut}}_t$$

With no sovereign-bond crowding out of private lending ($\varphi \cdot 0/N = 0$), bank optimization gives $R^{L,\text{aut}}_t = R^{D,\text{aut}}_t$ — the lending rate equals the deposit rate. The **closed-economy investment first-order condition** (9) then pins the autarky rate:

$$R^{D,\text{aut}}_t = \alpha A (K^{\text{aut}}_t)^{\alpha-1} + (1-\delta) \qquad\qquad (36)$$

Since investment is suspended during exclusion ($I_t = 0$), capital decumulates: $K^{\text{aut}}_t < K^{ss}$. The marginal product of capital $\alpha A (K^{\text{aut}}_t)^{\alpha-1}$ is therefore *higher* than in the non-crisis steady state. Hence:

$$R^{D,\text{aut}} = \alpha A (K^{\text{aut}})^{\alpha-1} + (1-\delta) > R^{D,ss} = R^* + \gamma s^{ss} \qquad\qquad (37)$$

**The autarky lending rate strictly exceeds the sovereign-ceiling rate under non-default.** The wedge $\Delta r^{L,\text{aut}} \equiv R^{D,\text{aut}} - R^{D,ss} > \gamma \hat{s}_0$ is the structural object that makes default-path output losses larger than non-default on impact. In the calibration (`do/16_model_irf.do`), $\Delta r^{L,\text{aut}}$ is pinned by the $h=0$ default data point; all subsequent horizons $h=1$–$4$ are then zero-free-parameter out-of-sample predictions.

---

## 6. Propositions

**Proposition 1 (Non-Default Channel).** Following a spread crisis with no default ($B_t < B^*_t$ for all $t$), the EMBIG spread $s_t$ is transmitted to domestic output through three sequential mechanisms: (i) the **sovereign-ceiling channel** — equations (11)–(12) — raises domestic deposit and lending rates contemporaneously, contracting household consumption through (3') and firm output through the working-capital channel (8) at $h=0$; (ii) the **balance-sheet channel** — equations (14)–(17') — progressively erodes bank net worth in proportion to pre-existing sovereign bond exposure $b^B_{-1}$, tightening the leverage constraint on domestic lending (27) and raising the user cost of investment (9) at $h=1$ through $h=4$; (iii) the **current-account adjustment** is gradual as international wholesale bank funding contracts but government bond rollover continues. The fiscal response is mildly expansionary until borrowing capacity is exhausted.

**Proposition 2 (Default-Linked Channel).** Following a spread crisis with default ($B_0 \ge B^*_0$), international wholesale funding collapses simultaneously with the suspension of sovereign external debt service, generating a **two-component sudden stop** (30) that forces an immediate and large current-account adjustment. The sovereign-ceiling transmission channel is severed, but the **autarky deposit rate rises sharply** (36)–(37), inducing front-loaded household consumption contraction. Banks **gamble for resurrection** (33)–(34), redirecting their full leverage capacity toward sovereign bonds and generating $+4$ to $+6$ pp of GDP in new sovereign bond accumulation without a conventional private credit squeeze. The fiscal response exhibits a **boom-bust pattern** (35) driven by the debt-service suspension windfall and subsequent restructuring conditionality.

**Corollary (The Central Paradox and the Role of the International Market).** The model provides a structural explanation for two findings that are paradoxical under single-channel frameworks. *First*, non-default spread crises generate large and persistent output losses without default: the sovereign-ceiling channel transmits the international spread to domestic credit conditions regardless of whether the government ultimately repays — the crisis itself is the shock, not the default. *Second*, default-linked episodes do not exhibit larger *measured* private-credit contractions despite the doom loop: banks do not deleverage (gross assets remain at $\lambda N_t$); they merely redirect capacity from private loans to sovereign bonds. The econometrician observes bank credit to government rising while private credit is unchanged in aggregate, because the reallocation leaves total assets constant. The key structural difference is that under non-default the sovereign-ceiling channel $\gamma s_t$ persists and deepens over the crisis horizon, while under default this channel is severed by market exclusion and replaced by the more acute sudden-stop and autarky mechanisms.

---

## 7. Competitive Equilibrium

A competitive equilibrium consists of sequences of quantities and prices

$$\{C_t,\, I_t,\, K_t,\, L_t,\, D_t,\, F_t,\, N_t,\, b^B_t,\, \ell_t,\, G_t,\, B_t,\, s_t,\, q_t,\, R^L_t,\, R^D_t,\, R^O_t,\, w_t\}_{t=0}^{\infty}$$

such that: (i) households maximize (1) subject to (2), satisfying Euler equation (3); (ii) firms maximize profits subject to the working-capital constraint, satisfying (7)–(9); (iii) banks optimize subject to the balance sheet (10) and leverage constraint (16), with funding costs (11)–(12) and lending rate (17'); (iv) the government satisfies the budget constraint (21) and the default condition (24); (v) foreign investors satisfy the participation constraint (18), yielding the endogenous spread (19); (vi) labor, domestic deposit, and international wholesale funding markets all clear; and (vii) the **goods market clears**:

$$Y_t = C_t + I_t + G_t + CA_t \qquad\qquad (39)$$

where $CA_t = Y_t - C_t - I_t - G_t$ is the current-account balance (positive = current-account surplus = net exporter of goods).

**Non-crisis steady state.** Setting all variables constant and assuming $CA^{ss} = 0$ (balanced current account in normal times), the steady state satisfies: $\beta(R^* + \gamma s^{ss}) = 1$ (from the Euler equation, pinning $s^{ss}$); output $Y^{ss} = A(K^{ss})^{\alpha}\Psi(R^{L,ss})$ (from equation 8); capital $\alpha A(K^{ss})^{\alpha-1} = R^{L,ss}-(1-\delta)$ (from equation 9); and lending rate $R^{L,ss} = R^* + \gamma s^{ss} + \varphi\, b^{B,ss}/N^{ss}$ (from equation 17'). The six conditions SS.1–SS.6 in `do/14_calibration.do` jointly determine the steady state $\{R^{D,ss}, R^{L,ss}, K^{ss}, Y^{ss}, \ell^{ss}, b^{B,ss}, N^{ss}, s^{ss}\}$ as functions of the structural parameters.

---

## 8. Log-Linearization and Analytical Solution

### 8.1 Conventions

**Log-deviation notation.** For any variable $X_t$ with steady-state value $X^{ss}$, define $\hat{x}_t \equiv \log(X_t/X^{ss}) \approx (X_t - X^{ss})/X^{ss}$. For the EMBIG spread (which enters additively, not multiplicatively), define the **level deviation** $\hat{s}_t \equiv s_t - s^{ss}$, so $\hat{s}_t = 0$ in normal times and $\hat{s}_t > 0$ during a crisis.

**Spread process.** The spread shock is modeled as a stationary AR(1) process:

$$\hat{s}_t = \rho_s \hat{s}_{t-1} + \varepsilon_t, \qquad \rho_s \in (0,1) \qquad\qquad \text{(LL.0)}$$

where $\rho_s$ is the **spread persistence parameter** — the fraction of the spread elevation that carries over to the next period — and $\varepsilon_t$ is a mean-zero innovation (the surprise component of the spread shock). $\rho_s$ is estimated directly from the data using a panel AR(1) regression on EMBIG spreads in `do/14_calibration.do`. The crisis-onset event corresponds to $\varepsilon_0 = \hat{s}_0 > 0$; subsequent periods evolve as $\hat{s}_h = \rho_s^h \hat{s}_0$. As $\rho_s \to 1$, the spread never returns to steady state and output losses are permanently larger.

### 8.2 Log-Linearized Equilibrium Equations

**Lending rate** (from 17'), with balance-sheet amplification $\Omega \equiv \varphi\, b^{B,ss}/N^{ss}$ and no new sovereign bond accumulation in the non-default path ($\hat{n}^{-1}_{t-1} = 0$):

$$\Delta r^{L,nd}_t = \gamma \hat{s}_t - \Omega \acute{n}_t \qquad\qquad \text{(LL.1a)}$$

where $\acute{n}_t \equiv \log(N_t/N^{ss})$ is the log-deviation of bank net worth from steady state (negative during a crisis), and $\Delta r^{L,nd}_t \equiv R^L_t - R^{L,ss}$ is the lending-rate deviation. The first term $\gamma \hat{s}_t > 0$ captures the sovereign-ceiling channel; the second term $-\Omega \acute{n}_t > 0$ (since $\acute{n}_t < 0$) captures the balance-sheet amplification that raises lending costs further as net worth erodes.

**Output** (from 8), where $\varepsilon_p \equiv \xi n R^{L,ss}/[1+\xi(R^{L,ss}-1)] > 0$ is the **elasticity of output with respect to the lending rate** through the working-capital channel (also called the working-capital output elasticity):

$$\hat{y}^{nd}_t = \alpha \hat{k}_t - \varepsilon_p\!\left(\gamma \hat{s}_t - \Omega \acute{n}_t\right) \qquad\qquad \text{(LL.2a)}$$

where $\hat{k}_t \equiv \log(K_t/K^{ss})$ is the capital stock deviation. Output deviates from steady state through two channels: the **capital accumulation channel** $\alpha \hat{k}_t$ (zero at $h=0$ since capital is predetermined, growing in magnitude thereafter) and the **working-capital channel** $-\varepsilon_p \Delta r^{L,nd}_t$ (operating immediately at $h=0$).

**Capital** (from 5 and 9), with $\eta \equiv 1/[(1-\alpha)(R^{L,ss}-(1-\delta))] > 0$ being the **interest semi-elasticity of the capital stock** with respect to the lending rate:

$$\hat{k}^{nd}_{t+1} = (1-\delta)\hat{k}_t - \eta\!\left(\gamma \hat{s}_t - \Omega \acute{n}_t\right) \qquad\qquad \text{(LL.3a)}$$

The capital stock evolves as a first-order linear difference equation driven by both the current spread (via the sovereign-ceiling channel on investment) and current net worth (via the balance-sheet amplification). Since $(1-\delta) < 1$, the system is stable.

**Bank net worth** (from 13 and 15), with:
- $\mathrm{B} \equiv b^{B,ss}/[N^{ss}(1+s^{ss})^2] > 0$ the **balance-sheet sensitivity**: how much net worth falls per unit of spread increase, proportional to the steady-state sovereign bond exposure relative to net worth
- $\Phi^N \equiv \Phi^L \lambda N^{ss}/\ell^{ss} + \Phi^B - \Phi^D - \Phi^O$ the **net-worth multiplier**: the fraction of last period's net-worth deviation that persists into this period through retained earnings; requires $|\Phi^N| < 1$ for stability
- $\Phi^O \equiv R^{O,ss}F^{ss}/N^{ss}$ the **wholesale-funding income share**: the share of bank income consumed by international wholesale funding costs

$$\acute{n}^{nd}_{t+1} = \Phi^N \acute{n}_t - \mathrm{B}\,\hat{s}_{t+1} - \Phi^O \gamma \hat{s}_t \qquad\qquad \text{(LL.5a)}$$

Net worth is driven down by: (i) the **bond price capital loss** $\mathrm{B}\hat{s}_{t+1}$ — the balance-sheet channel operating through the mark-to-market value of the sovereign bond portfolio; and (ii) the **wholesale funding cost increase** $\Phi^O \gamma \hat{s}_t$ — the sovereign-ceiling channel operating on the liability side of the bank's income statement, compressing net interest margins.

**Private credit** (from 16 and the balance sheet, with $\hat{n}^B_t = 0$ in the non-default path):

$$\hat{\ell}^{nd}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\,\acute{n}_t \qquad\qquad \text{(LL.6)}$$

Credit contracts in proportion to the net-worth deviation, amplified by the leverage ratio $\lambda N^{ss}/\ell^{ss} > 1$. Because $\acute{n}_0 = -\mathrm{B}\hat{s}_0$ is small when $b^{B,ss}/N^{ss}$ is modest, credit contraction is negligible at $h=0$ and deepens progressively as net worth erodes — matching the empirical pattern.

**Consumption** (from 3', the modified Euler equation with the sovereign ceiling substituted in):

$$\hat{c}_t = \mathbb{E}_t[\hat{c}_{t+1}] - \frac{\gamma}{\sigma \, R^{D,ss}}\,\hat{s}_{t+1} \qquad\qquad \text{(LL.7)}$$

This is a forward-looking IS-type equation: current consumption falls relative to future consumption whenever $\hat{s}_{t+1} > 0$. Since the spread is persistent ($\rho_s < 1$), the household anticipates elevated future spreads and accelerates intertemporal substitution toward the future, generating a contemporaneous demand contraction. Solving (LL.7) forward:

$$\hat{c}_t = -\frac{\gamma}{\sigma \, R^{D,ss}} \sum_{h=t}^{\infty} \rho_s^{h-t} \mathbb{E}_t[\hat{s}_h]$$

The impact contraction is $\hat{c}_0 = -(\gamma/\sigma R^{D,ss}) \hat{s}_0/(1-\rho_s)$, larger when the spread is more persistent ($\rho_s \to 1$), the ceiling pass-through is stronger ($\gamma \to 1$), or the household is less willing to substitute intertemporally ($\sigma \to 0$).

### 8.3 Closed-Form Impulse Responses (Non-Default)

Given the AR(1) process $\hat{s}_t = \rho_s^t \hat{s}_0$ and the initial impact balance-sheet loss $\acute{n}_0 = -\mathrm{B}\hat{s}_0$ (from LL.5a with $\acute{n}_{-1} = 0$ and $\hat{s}_{-1} = 0$), iterating forward yields:

**Net worth at horizon $h$:**

$$\acute{n}^{nd}_h = (\Phi^N)^h \acute{n}_0 \;-\; \sum_{h'=1}^{h} (\Phi^N)^{h-h'}\!\left(\mathrm{B}\,\rho_s^{h'} + \Phi^O \gamma\,\rho_s^{h'-1}\right)\hat{s}_0 \qquad\qquad \text{(IR.1)}$$

where the first term $(\Phi^N)^h \acute{n}_0$ is the **decaying legacy** of the initial net-worth shock, and the summation is the **cumulative fresh loss** from the still-elevated spread at each subsequent period.

> **Correction note.** An earlier draft wrote this as $-(\mathrm{B}+\Phi^O\gamma)\sum(\Phi^N)^{h-h'}\rho_s^{h'}\hat{s}_0$, which factors incorrectly. The $\mathrm{B}$ term multiplies $\rho_s^{h'}$ (the spread *at* period $h'$, entering through $-\mathrm{B}\hat{s}_{h'+1}$ in LL.5a) while the $\Phi^O\gamma$ term multiplies $\rho_s^{h'-1}$ (the spread *one period earlier*, entering through $-\Phi^O\gamma\hat{s}_{h'}$). These only coincide when $\rho_s = 1$. In `do/16_model_irf.do`, LL.5a is iterated directly period by period, which is numerically exact and avoids the summation.

**Capital stock at horizon $h$** (from iterating LL.3a with $\hat{k}_0 = 0$):

$$\hat{k}^{nd}_h = -\eta\gamma \sum_{h'=1}^{h}(1-\delta)^{h-h'}\rho_s^{h'-1}\hat{s}_0 \;+\; \eta\Omega \sum_{h'=1}^{h}(1-\delta)^{h-h'}\acute{n}_{h'-1} \qquad\qquad \text{(IR.2)}$$

The capital stock deviation is the **discounted cumulative sum** of past lending-rate increases, weighted by the depreciation factor $(1-\delta)$. Two forces drive the capital decline: the direct effect of elevated spreads on investment demand (the $\eta\gamma$ term, proportional to the sovereign ceiling) and the indirect balance-sheet amplification (the $\eta\Omega$ term, proportional to $\Omega = \varphi\, b^{B,ss}/N^{ss}$, the steady-state sovereign exposure ratio).

**Output at horizon $h$:**

$$\hat{y}^{nd}_h = \underbrace{\alpha \hat{k}^{nd}_h}_{\substack{\text{capital} \\ \text{channel}}} \;-\; \underbrace{\varepsilon_p\gamma\rho_s^h \hat{s}_0}_{\substack{\text{working-capital} \\ \text{channel}}} \;+\; \underbrace{\varepsilon_p\Omega\,\acute{n}^{nd}_h}_{\substack{\text{balance-sheet} \\ \text{amplification}}} \qquad\qquad \text{(IR.3)}$$

All three terms are negative (since $\hat{k}^{nd}_h < 0$, $\hat{s}_0 > 0$, and $\acute{n}^{nd}_h < 0$), and their sum is negative at all horizons:

$$\hat{y}^{nd}_h < 0 \quad \forall\, h \ge 0 \quad \text{whenever} \quad \hat{s}_0 > 0 \qquad\qquad \text{(IR.3a)}$$

The **trough** occurs where the deepening balance-sheet and capital effects outweigh the decaying spread shock — between $h=1$ and $h=3$ for plausible parameters, consistent with the empirical trough of $-3.76$ pp at $h=1$.

**Private credit at horizon $h$:**

$$\hat{\ell}^{nd}_h = \frac{\lambda N^{ss}}{\ell^{ss}}\,\acute{n}^{nd}_h \qquad\qquad \text{(IR.4)}$$

insignificant at $h=0$ (when $\acute{n}_0 = -\mathrm{B}\hat{s}_0$ is small if $b^{B,ss}/N^{ss}$ is modest), deepening to $-3.48$ pp at $h=4$ as net worth erodes.

**Lemma (Persistence Without Mean Reversion).** Let $\rho_s, \Phi^N \in (0,1)$. Then: (i) $\hat{y}^{nd}_h < 0$ for all $h \ge 0$; (ii) the output trough occurs at finite $h^* > 0$; (iii) the cumulative output loss diverges as $\rho_s \to 1$, because when the spread does not return to steady state, net worth never recovers and capital remains permanently depressed.

---

### 8.4 The Default Path: Modified Log-Linearized System

Under default, three structural breaks modify the transmission equations.

**Break 1 — Collapse of international wholesale funding.** Setting $F_t = 0$ severs the sovereign-ceiling link. The domestic deposit rate is now set by closed-economy equilibrium (§5.4) rather than by international arbitrage. Equation LL.7 is replaced by the autarky Euler equation with $R^{D,\text{aut}}_t > R^{D,ss}$, inducing a larger front-loaded consumption contraction than under non-default.

**Break 2 — Gambling for resurrection.** When condition (33) holds, the bank's optimal policy is $\ell^{\text{def}}_t = 0$ and $b^{B,\text{def}}_t = \lambda N_t$. Private lending is zero, so the log-deviation of credit diverges:

$$\hat{\ell}^{\text{def}}_t \approx -\infty \qquad\text{(private lending collapses to zero)} \qquad\qquad \text{(LL.6D)}$$

> **Proposition 3 (Credit Paradox).** For all $h \ge 0$, private lending under the default path is strictly lower than under the non-default path ($\hat{\ell}^{\text{def}}_h \ll \hat{\ell}^{nd}_h < 0$). Yet the private-credit contraction in aggregate data may appear smaller under default: banks do not reduce total assets (gross assets remain $\lambda N_t$), they merely redirect them from private loans to sovereign bonds. The econometrician measuring *total bank credit* (including government) observes a reallocation, not a deleveraging, so the *private-credit* squeeze is invisible in credit aggregates. This is why default-linked episodes show no statistically larger measured private-credit contraction despite the doom loop.

**Break 3 — Doom-loop sovereign bond accumulation.** Log-linearizing $\Delta b^B_0 = \lambda N_0 - b^B_{-1}$:

$$\hat{n}^{B,\text{def}}_t = \left(\frac{\lambda N^{ss}}{b^{B,ss}} - 1\right) - \frac{\lambda N^{ss}}{b^{B,ss}}\,\mathrm{B}\,\hat{s}_0 \qquad\qquad \text{(LL.10)}$$

The first term (positive when $\lambda N^{ss} > b^{B,ss}$, i.e., banks had spare leverage capacity pre-crisis) represents the **reallocation gain**; the second term represents the **net-worth erosion cost**. The net accumulation is positive precisely when the gambling-for-resurrection condition (33) holds.

**Output under default** (from equation 8 with $\ell^{\text{def}}_t = 0$ and $I_t = 0$). With investment suspended, the capital stock declines at the depreciation rate: $K^{\text{def}}_h = (1-\delta)^h K_0$, so $\hat{k}^{\text{def}}_h = h \ln(1-\delta) < 0$ (a negative number growing in magnitude). Output:

$$\hat{y}^{\text{def}}_h = \alpha \cdot h \ln(1-\delta) \;-\; \varepsilon_p \, \Delta r^{L,\text{aut}} \qquad\qquad \text{(IR.5)}$$

where:
- $\alpha \cdot h \ln(1-\delta) < 0$ is the **capital-depletion channel**: each period without investment the capital stock shrinks by factor $(1-\delta)$, reducing output by $\alpha$ times the capital log-deviation. This term is zero at $h=0$ and becomes more negative at $h=1,2,3,4$ — generating the deepening output losses under default even after the initial impact
- $-\varepsilon_p \Delta r^{L,\text{aut}} < 0$ is the **autarky working-capital channel**: with private lending at zero, firms cannot pre-finance wages, and output contracts by $\varepsilon_p$ times the autarky lending-rate wedge. This is the **impact loss at $h=0$**

> **Correction note.** An earlier draft wrote the capital term as $\alpha(1-\delta)^h \hat{k}_0$ which equals zero since $\hat{k}_0 = 0$ (capital is predetermined at the onset). This incorrectly implies no capital-depletion effect, contradicting the stated deepening of output losses under default. The correct derivation uses the log of the depreciated capital stock $\hat{k}^{\text{def}}_h = h\ln(1-\delta)$, which accumulates over time and drives the progressive worsening of default-path outcomes.

---

## 9. Mapping to the Empirical Findings

### Table 1 — Model Sign Predictions by Channel, Path, and Horizon

| Variable | Path | $h=0$ | $h=1$–$2$ | $h=3$–$4$ | Mechanism |
|---|---|:---:|:---:|:---:|---|
| **Output** | Non-default | − | − (deepens) | −/≈0 | working capital (IR.3) + balance sheet |
| | Default | −− | −− | −− | autarky working capital + capital depletion (IR.5) |
| | *Comparison* | nd > d | nd > d | nd > d | default always deeper |
| **Private credit** | Non-default | − (small) | − (growing) | − (deep, −3.5 pp) | net-worth erosion via IR.4 |
| | Default | ≈0 | ≈0 | ≈0 | banks redirect to sovereign bonds (LL.6D) |
| | *Comparison* | nd < d | nd < d | nd < d | ND has worse *measured* credit |
| **Investment** | Non-default | − | − | −/≈0 | user cost rises via (9) |
| | Default | ≈0 | − | −− | $I=0$ during exclusion; capital depletes |
| | *Comparison* | nd > d | nd > d | nd > d | |
| **Claims on govt** | Non-default | ≈0 | ≈0 | ≈0 | no doom-loop incentive ($\hat{n}^B=0$) |
| | Default | + | ++ (+4–6 pp) | + | gambling for resurrection (33)–(34) |
| | *Comparison* | nd < d | nd < d | nd < d | default exclusively |
| **Current account** | Non-default | + (gradual) | + (fading) | ≈0 | wholesale cost rises, flows continue |
| | Default | ++ | ++ (+2.6 pp) | + | two-component sudden stop (30) |
| | *Comparison* | nd < d | nd < d | nd ≈ d | default front-loaded |
| **Govt expenditure** | Non-default | + | + | + | retained market access, gradual tightening |
| | Default | ++ (+2.2 pp) | − | −− (−1.7 pp) | windfall then conditionality (35) |
| | *Comparison* | nd < d | nd > d | nd > d | paths cross near $h \approx 1$ |

---

### Table 2 — Empirical Findings, Model Equations, and Parameter Conditions

| Finding | Model equation(s) | Parameter condition required |
|---|---|---|
| **(a)** ND output loss (−1.86 pp at $h=0$, −3.76 pp at $h=1$) | LL.2a, IR.3 | $\varepsilon_p > 0$ (requires $\xi > 0$); $\gamma > 0$; $\hat{s}_0 > 0$ |
| **(b)** ND trough at $h=1$, not $h=0$ | IR.3 combined with IR.1–IR.2 | $\Phi^N \in (0,1)$ and $\Omega > 0$: balance-sheet deepening overtakes spread decay at $h=1$ |
| **(c)** Default deeper and front-loaded (−3.42 pp at $h=0$, −5.50 pp at $h=1$) | IR.5 and eq. (37) | $\Delta r^{L,\text{aut}} > \gamma\hat{s}_0$; guaranteed by $K^{\text{aut}} < K^{ss}$ from capital depletion |
| **(d)** Credit contraction in ND only (−4.78 pp at $h=4$; ≈0 default) | IR.4 (ND) and LL.6D (default) | ND: $\acute{n} < 0$; default: gambling-for-resurrection (33) holds, $\ell = 0$ |
| **(e)** Doom loop +4–6 pp in default only | Eqs. (33)–(34) and LL.10 | $\lambda N_0 > b^B_{-1}$ (spare leverage capacity); $\mu\theta R^{BDK} + \gamma^B N \ge R^L(1-p^{RDL})$ |
| **(f)** Diabolic loop absorbed by lagged $b^B$ | Eqs. (14)–(15) and (17') | Amplification $\propto b^B_{t-1}$: controlling $L.b^B$ in LP absorbs $\mathrm{B}$ in LL.5a |
| **(g)** CA adjustment larger and earlier in default (+2.61 pp at $h=1$) | Eq. (30) vs. eq. (28) | $F_{-1} > 0$ (pre-crisis wholesale funding) and $\delta^B > 0$ (maturing debt = service flow) |
| **(h)** Govexp boom-bust in default (+2.19 pp at $h=0$, −1.69 pp at $h=4$) | Eq. (35) | $\delta^B > 0$ (windfall); exclusion $> 2$ periods so conditionality bites at $h=2$–$4$ |
| **(i)** Govexp sustained in ND (positive throughout) | Eq. (21) + ceiling (11) | $\gamma < \bar{\gamma}$: spread rise insufficient to trigger solvency crisis; market access retained |

**Theorem (Joint Consistency).** Under the conditions in Table 2, the log-linearized model under both paths — non-default IR.1–IR.4 and default IR.5, LL.6D, LL.10 — jointly reproduces the qualitative and directional predictions of all nine empirical findings. No individual condition contradicts another; there exists an open set $\Theta \subset \mathbb{R}^{11}_+$ of parameter values under which all findings hold simultaneously. The model is therefore internally consistent: no result requires a parameter restriction that invalidates another.

---

## 10. Calibration and Numerical Solution

### 10.1 Parameter Groups

**Group 1 — Externally calibrated from the literature.** Standard annual EM values: discount factor $\beta = 0.96$; risk aversion $\sigma = 2$; capital share $\alpha = 0.33$; depreciation $\delta = 0.10$; world risk-free rate $R^* = 0.04$; recovery rate $\theta = 0.62$ (Cruces and Trebesch 2013, mean EM haircut $\approx 38\%$); re-entry probability $\mu = 0.22$ (Gelos et al. 2011, mean exclusion $\approx 4.5$ years); maturing fraction $\delta^B = 0.22$ (average EM sovereign debt maturity $\approx 4.5$ years); bank leverage $\lambda = 10$.

**Group 2 — Data-determined steady-state ratios.** Computed from *tranquil-period means* (non-crisis, non-continuation country-years: `onset_all==0 & continuation==0`) in `panel_lp.dta`: background spread $s^{ss}$; expenditure shares $s_C, s_I, s_G, s_{CA}$; tax rate $\tau$ (from `revenue_gdp`); bank sovereign exposure $b^{B,ss}/Y$ (from `claims_govt`); bank net worth $N^{ss}/Y = (b^{B,ss}/Y + \ell^{ss}/Y)/\lambda$; and hence the structural elasticities $\varepsilon_p$, $\eta$, $\mathrm{B}$, and amplification $\Omega$. Implemented in `do/14_calibration.do`.

**Group 3 — Estimated stochastic processes.** Spread persistence $\rho_s$ from a panel fixed-effects AR(1) on EMBIG spreads (`xtreg spr_mean L.spr_mean, fe`); income process parameters $(\rho_y, \sigma_y)$ from a country-trend-detrended AR(1) on log GDP per capita. The deep transmission parameters $\{\xi, \varphi, \Phi^N\}$ are calibrated by **Simulated Method of Moments (SMM) targeting the empirical non-default IRF only**.

### 10.2 Two-Block Solution

**Block 1 — Nonlinear default block** (`do/15_solve_default.do`). Canonical Arellano (2008) value-function iteration (VFI) over states $(B, y)$ in Mata (Stata's matrix programming language): Tauchen (1986)-discretized income process with 21 income states; the default decision (25); and the bond-price fixed point (18)–(19) with recovery $\theta$ and re-entry $\mu$. Simulation of 100,000 periods (with 1,000 burn-in) yields model moments (annual default frequency, mean/SD of the spread, mean debt/GDP ratio) and a **crisis-event spread path** at horizons $h=-2$ to $h=4$ (averages across all threshold-crossing events in the simulation) that feeds Block 2.

**Block 2 — Log-linear transmission block** (`do/16_model_irf.do`). The model spread path enters equations LL.5a, LL.3a, LL.2a, and IR.4, iterated directly period by period. The SMM grid search over $\{\xi, \varphi, \Phi^N\}$ minimizes:

$$\min_{\xi,\,\varphi,\,\Phi^N} \sum_{h=0}^{4} \left[\left(\hat{y}^{nd,\text{model}}_h - \hat{y}^{nd,\text{LP}}_h\right)^2 + \left(\hat{\ell}^{nd,\text{model}}_h - \hat{\ell}^{nd,\text{LP}}_h\right)^2\right]$$

where $\hat{y}^{nd,\text{LP}}_h$ and $\hat{\ell}^{nd,\text{LP}}_h$ are the empirical LP point estimates from `irf_nd.dta` and `11_channels.do`.

### 10.3 Identification and Out-of-Sample Validation

The calibration uses **only the non-default IRF** to pin $\{\xi, \varphi, \Phi^N\}$ — six moment conditions (output and credit at five horizons) identifying three parameters. The **default path is then validated out of sample**: the single autarky wedge $\Delta r^{L,\text{aut}}$ is pinned to the $h=0$ default output data point; $h=1$ through $h=4$ of the default IRF are then **zero-free-parameter predictions** generated by IR.5. The degree to which these match the empirical default LP constitutes the model's central test, converting the framework from a fitting exercise into a falsifiable structural prediction.

---

## 11. Scope and Limitations

The sovereign-ceiling pass-through $\gamma$ and the wholesale-funding share $\omega$ are treated as constants, whereas in practice both vary with the degree of financial integration, the currency composition of bank liabilities, and the institutional credibility of the sovereign guarantee. Modeling $\gamma$ as endogenous would introduce richer cross-country predictions but at the cost of additional identification requirements.

The default decision abstracts from renegotiation dynamics, long-maturity debt with duration risk, and partial default. TFP $A$ is constant, abstracting from the endogenous productivity feedback modeled in Mendoza and Yue (2012), whereby firms in default episodes lose access to imported intermediate inputs and productivity falls directly.

These are deliberate simplifications consistent with the model's purpose — analytical clarity on the transmission channels identified in the empirical analysis — rather than quantitative discipline. A full quantitative extension would integrate the dual funding structure with a Bocola (2016)-type balance-sheet model for the non-default channel and a Gennaioli, Martin, and Rossi (2014) framework for the default-linked channel, calibrated to the full 52-country emerging-market sample.

---

## Notation Summary

| Symbol | Definition |
|---|---|
| $C_t$ | Household real consumption in period $t$ |
| $D_t$ | Household deposits at domestic banks |
| $\beta$ | Household subjective discount factor |
| $\sigma$ | Inverse of the intertemporal elasticity of substitution (CRRA coefficient) |
| $R^D_t$ | Gross domestic deposit interest rate |
| $w_t$ | Equilibrium real wage |
| $\Pi_t$ | Lump-sum bank profits distributed to households |
| $Y_t$ | Real output (GDP) |
| $A$ | Total factor productivity (constant) |
| $K_t$ | Physical capital stock at beginning of period $t$ |
| $\alpha$ | Capital share in production ($\approx 0.33$) |
| $\delta$ | Annual physical capital depreciation rate ($= 0.10$) |
| $I_t$ | Gross investment |
| $\xi$ | Working-capital share: fraction of wage bill pre-financed by bank loans |
| $\ell_t$ | Total domestic bank lending to firms |
| $n = (1-\alpha)/\alpha$ | Ratio of labor share to capital share |
| $\Psi(R^L_t)$ | Working-capital efficiency function; decreasing in $R^L_t$ |
| $\varepsilon_p$ | Elasticity of output with respect to the lending rate (working-capital channel) |
| $\eta$ | Interest semi-elasticity of the capital stock |
| $R^L_t$ | Gross domestic bank lending rate charged to firms |
| $b^B_t$ | Bank holdings of dollar-denominated sovereign bonds |
| $F_t$ | Bank international wholesale funding (from foreign interbank markets) |
| $N_t$ | Bank net worth (equity) |
| $\lambda$ | Maximum bank leverage ratio (assets/equity, $= 10$) |
| $\omega$ | Share of domestic deposits in total bank funding |
| $\varphi$ | Balance-sheet pass-through: sensitivity of lending rate to sovereign-bond/net-worth ratio |
| $R^O_t$ | Gross international wholesale funding rate paid by domestic banks |
| $R^*$ | International risk-free interest rate (US Treasury, $\approx 0.04$) |
| $s_t$ | EMBIG sovereign spread (decimal; $s_t = 0.10$ means 1000 bps) |
| $\gamma$ | Sovereign-ceiling pass-through parameter ($\in (0,1]$) |
| $q_t = 1/(1+s_t)$ | Sovereign bond price |
| $\pi_{t+1}$ | Rational expectation of sovereign default probability in $t+1$ |
| $\theta$ | Recovery rate upon restructuring ($= 0.62$, i.e., 38% haircut) |
| $G_t$ | Primary government expenditure |
| $B_t$ | Face value of outstanding sovereign external debt |
| $\delta^B$ | Fraction of external debt maturing per period ($= 0.22$) |
| $\tau$ | Proportional tax rate on output |
| $\beta_p$ | Government discount factor (may differ from household $\beta$) |
| $\phi$ | Direct output cost of market exclusion during default |
| $\mu$ | Probability of re-entering international markets after default ($= 0.22$) |
| $\hat{x}_t$ | Log-deviation of variable $X_t$ from non-crisis steady state |
| $\hat{s}_t$ | Level deviation of spread from steady state: $\hat{s}_t = s_t - s^{ss}$ |
| $\rho_s$ | Spread persistence (AR(1) coefficient, estimated from EMBIG panel data) |
| $\rho_y$ | Income persistence (AR(1) coefficient, estimated from GDP data) |
| $\Omega = \varphi b^{B,ss}/N^{ss}$ | Balance-sheet amplification in the log-linear system |
| $\mathrm{B} = b^{B,ss}/[N^{ss}(1+s^{ss})^2]$ | Balance-sheet sensitivity of net worth to the spread |
| $\Phi^N$ | Net-worth multiplier in LL.5a (fraction of net-worth shock that persists) |
| $\Phi^O = R^{O,ss}F^{ss}/N^{ss}$ | Wholesale-funding income share in bank net worth |
| $\acute{n}_t = \log(N_t/N^{ss})$ | Log-deviation of bank net worth |
| $\hat{k}_t = \log(K_t/K^{ss})$ | Log-deviation of capital stock |
| $\Delta r^{L,\text{aut}}$ | Autarky lending-rate wedge under default: $R^{L,\text{aut}} - R^{L,ss}$ |
| $CA_t$ | Current-account balance ($= Y_t - C_t - I_t - G_t$; positive = surplus) |
| $s_C, s_I, s_G, s_{CA}$ | Steady-state expenditure shares of consumption, investment, government, current account |

---

*Implementation: `do/14_calibration.do` (Groups 1–3 calibration) · `do/15_solve_default.do` (Arellano VFI, nonlinear default block) · `do/16_model_irf.do` (log-linear transmission, SMM, overlay figure `fig_model_vs_data.pdf`)*
