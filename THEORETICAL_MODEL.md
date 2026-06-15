# A Stylized Model of Sovereign Spread Crises With and Without Default

**Theoretical Appendix** — *Companion to "The Aftermath of Sovereign Spread Crises With and Without Default"*

---

## Overview

This appendix develops a small open economy model that rationalizes the central empirical findings of the paper: (i) sovereign spread crises generate large and persistent output losses *even when no default occurs*; (ii) default-linked crises are roughly twice as costly and front-loaded; (iii) the private-credit and investment contractions are concentrated in *non-default* episodes; (iv) bank accumulation of sovereign debt (the doom loop) is exclusive to *default-linked* episodes; and (v) the current-account adjustment is larger and earlier under default.

The model's defining feature is a **separation between the international sovereign bond market** — where the government issues dollar-denominated external debt priced by foreign investors — **and the domestic credit market** — where banks lend to firms in local currency. Transmission from the former to the latter runs through two channels: a **wholesale funding-cost channel** (the sovereign ceiling) and a **balance-sheet channel** (bond-price losses erode bank net worth). This dual structure delivers both a non-default transmission mechanism and a distinct default mechanism.

The solution uses a **two-block decomposition** following Bocola (2016): a nonlinear sovereign-risk block that determines the endogenous default decision and the equilibrium spread (solved by value-function iteration), and a log-linear transmission block that maps the spread into output, credit, investment, and the current account (solved in closed form).

---

## 1. Environment and Agents

Consider a small open economy in discrete time, $t = 0, 1, 2, \dots$ Five agents: a representative household, a representative firm, a domestic banking sector with access to both domestic and international wholesale funding, a government, and a continuum of competitive foreign investors.

---

### 1.1 Households

The representative household maximizes:

$$\max \; \mathbb{E}_0 \sum_{t=0}^{\infty} \beta^t \, \frac{C_t^{1-\sigma}}{1-\sigma} \qquad\qquad (1)$$

where $\beta \in (0,1)$ and $\sigma > 0$ is the inverse EIS. Budget constraint:

$$C_t + D_{t+1} = R^D_t D_t + w_t + \Pi_t \qquad\qquad (2)$$

First-order condition (Euler equation):

$$C_t^{-\sigma} = \beta \, R^D_{t+1} \, \mathbb{E}_t\!\left[C_{t+1}^{-\sigma}\right] \qquad\qquad (3)$$

The household holds only domestic deposits $D_t$ and has no direct access to international capital markets.

---

### 1.2 Firms

Cobb-Douglas production:

$$Y_t = A K_t^{\alpha} L_t^{1-\alpha}, \qquad \alpha \in (0,1) \qquad\qquad (4)$$

Capital accumulation ($\delta$ = depreciation rate):

$$K_{t+1} = (1-\delta)K_t + I_t \qquad\qquad (5)$$

**Working-capital constraint:** a fraction $\xi \in (0,1)$ of the wage bill must be pre-financed by short-term bank loans. Total credit demand:

$$\ell_t = \xi w_t + I_t \qquad\qquad (6)$$

Labor demand with the effective cost of labor raised to $[1+\xi(R^L_t-1)]w_t$:

$$(1-\alpha)A K_t^{\alpha} L_t^{-\alpha} = [1+\xi(R^L_t-1)]\,w_t \qquad\qquad (7)$$

Equilibrium output as a function of the domestic lending rate $R^L_t$:

$$Y_t = A K_t^{\alpha} \cdot \Psi(R^L_t) \qquad\qquad (8)$$

where

$$\Psi(R^L_t) \equiv \left[\frac{1-\alpha}{1+\xi(R^L_t-1)}\right]^{\!n}, \qquad n = \frac{1-\alpha}{\alpha}$$

Since $\Psi'(R^L_t) < 0$, output is **strictly decreasing in the lending rate** for any $\xi > 0$. This is the key result of the firm block: any upstream rise in $R^L_t$ generates an *immediate* within-period output contraction.

Investment demand (user cost condition):

$$\alpha A K_{t+1}^{\alpha-1} = R^L_t - (1-\delta) \qquad\qquad (9)$$

---

### 1.3 Banking Sector with Dual Funding Structure

#### 1.3.1 Balance Sheet

$$b^B_t + \ell_t = D_t + F_t + N_t \qquad\qquad (10)$$

where $b^B_t$ = sovereign bond holdings (marked to market), $\ell_t$ = domestic lending, $D_t$ = deposits, $F_t$ = international wholesale funding, $N_t$ = net worth.

#### 1.3.2 International Wholesale Funding — The Sovereign Ceiling

Foreign lenders charge domestic banks a funding rate reflecting the sovereign's credit risk:

$$R^O_t = R^* + \gamma s_t, \qquad \gamma \in (0,1] \qquad\qquad (11)$$

where $s_t$ is the sovereign EMBIG spread. At the interior equilibrium where both funding sources are active, arbitrage gives:

$$R^D_t = R^O_t = R^* + \gamma s_t \qquad\qquad (12)$$

**Equation (12) pins the domestic deposit rate to the international spread.** Substituting into the Euler equation:

$$C_t^{-\sigma} = \beta\,(R^* + \gamma s_{t+1})\,\mathbb{E}_t\!\left[C_{t+1}^{-\sigma}\right] \qquad\qquad (3')$$

A spread spike raises the household's effective discount rate, generating a contemporaneous consumption contraction *before any balance-sheet effect*.

#### 1.3.3 Net Worth and the Balance-Sheet Channel

$$N_{t+1} = R^L_t \ell_t + (R^B_t + \Delta q_{t+1}) b^B_t - R^D_t D_t - R^O_t F_t \qquad\qquad (13)$$

Using $q_t = 1/(1+s_t)$, a spread increase generates a capital loss on outstanding bonds:

$$\Delta q_t \approx -\frac{\Delta s_t}{(1+s_t)^2} < 0 \qquad\qquad (14)$$

so the balance-sheet shock to net worth is:

$$\Delta N_t = b^B_{t-1} \cdot \Delta q_t = -\frac{b^B_{t-1}\,\Delta s_t}{(1+s_t)^2} < 0 \qquad\qquad (15)$$

The loss is **proportional to pre-existing exposure** $b^B_{t-1}$. The two channels are additive: the sovereign-ceiling channel raises the *price* of credit; the balance-sheet channel reduces the *quantity* banks can supply.

#### 1.3.4 Leverage Constraint and the Lending Rate

$$b^B_t + \ell_t \le \lambda N_t, \qquad \lambda > 1 \qquad\qquad (16)$$

The domestic lending rate combines weighted funding costs and a balance-sheet risk premium:

$$R^L_t = \omega R^D_t + (1-\omega) R^O_t + \varphi\,\frac{b^B_{t-1}}{N_t} \qquad\qquad (17)$$

Substituting (11)–(12):

$$\boxed{R^L_t = R^* + \gamma s_t + \varphi\,\frac{b^B_{t-1}}{N_t}} \qquad\qquad (17')$$

**Equation (17') is the central transmission equation.** It decomposes the lending rate into: a baseline $R^*$; a **sovereign-ceiling term** $\gamma s_t$ affecting all borrowers contemporaneously; and a **balance-sheet term** $\varphi(b^B_{t-1}/N_t)$ that amplifies the effect in proportion to pre-existing sovereign exposure.

> **Remark on identification.** The amplification term $\varphi(b^B_{t-1}/N_t)$ is proportional to the *lagged* position $b^B_{t-1}$. This is why the empirical finding that the doom-loop channel is absorbed by the lagged stock of bank claims on government is structurally meaningful: controlling for $L.b^B$ in the local projections isolates crisis-period accumulation from the pre-existing exposure that already drives the balance-sheet transmission.

---

## 2. International Sovereign Bond Market and the Endogenous Spread

### 2.1 Foreign Investor Zero-Profit Condition

$$q_t(1+R^*) = \mathbb{E}_t\!\left[(1-\pi_{t+1}) + \pi_{t+1}\,\theta\right] \qquad\qquad (18)$$

where $\pi_{t+1} \in [0,1]$ is the rational default probability and $\theta \in [0,1)$ is the recovery rate.

### 2.2 Equilibrium Spread

Using $q_t = 1/(1+s_t)$:

$$s_t = \frac{(1+R^*)(1-\theta)\,\pi_{t+1}}{1 - \pi_{t+1}(1-\theta)} \qquad\qquad (19)$$

The spread is strictly increasing in $\pi_{t+1}$ and decreasing in $\theta$. For small $\pi_{t+1}$, $s_t \approx (1+R^*)(1-\theta)\pi_{t+1}$.

### 2.3 Spread Crisis Definition

$$\Omega_t = \mathbf{1}\!\left[s_t \ge \bar{s}\right], \qquad \bar{s} = 1000\text{ bps} \qquad\qquad (20)$$

consistent with Pescatori-Sy (2007). The onset indicator $\Omega_t = 1$ in the first period the threshold is crossed.

---

## 3. Government and the Default Decision

### 3.1 Budget Constraint

$$G_t + (1-\delta^B)B_t = \tau Y_t + q_t\!\left[B_{t+1} - (1-\delta^B)B_t\right] \qquad\qquad (21)$$

where $G_t$ is primary expenditure, $B_t$ is external debt face value, $\delta^B$ is the fraction maturing per period, $\tau Y_t$ is tax revenue.

### 3.2 Recursive Default Problem

**Value of repayment:**

$$V^R_t(B_t,K_t) = \max_{G_t,\, B_{t+1}}\left\{ U(G_t) + \beta_p\,\mathbb{E}_t\!\left[V_{t+1}(B_{t+1},K_{t+1})\right] \right\} \qquad\qquad (22)$$

subject to (21), (5), $G_t \ge \bar{G}$, $B_{t+1} \le \bar{B}(s_t)$.

**Value of default** (excluded from markets, re-enters with probability $\mu$):

$$V^D_t(K_t) = \max_{G_t}\left\{ U(G_t) + \beta_p\,\mathbb{E}_t\!\left[\mu V^R_{t+1}(0,K_{t+1}) + (1-\mu)V^D_{t+1}(K_{t+1})\right] \right\} \qquad\qquad (23)$$

subject to autarky budget $G_t = (1-\phi)\tau Y_t$, with $I_t = 0$ during exclusion.

### 3.3 Default Threshold

Default occurs when $V^D_t(K_t) \ge V^R_t(B_t, K_t)$. This defines a threshold $B^*_t(K_t)$ decreasing in $s_t$ and $K_t$. Foreign investors form rational expectations:

$$\pi_{t+1} = \Pr\!\left(B_{t+1} \ge B^*_{t+1}(K_{t+1}) \mid \mathcal{F}_t\right) \qquad\qquad (24)$$

Equations (19) and (24) form the **fixed-point system** closing the model.

---

## 4. The Non-Default Path

Suppose the government repays in all periods. A crisis occurs at $t=0$ because investors revise $\pi_{t+1}$ upward, driving $s_0$ past $\bar{s}$.

**$h=0$ — Sovereign ceiling on consumption.** The spread spike raises $R^D_0 = R^* + \gamma s_0$ and, via (3'), the household defers consumption. This demand contraction operates before any credit-quantity effect.

**$h=0$ — Sovereign ceiling on working capital.** The lending rate rises via (17'). From (8):

$$\frac{dY_0}{dR^L_0} = A K_0^{\alpha}\,\Psi'(R^L_0) < 0, \qquad \Psi'(R^L_0) = -\frac{\xi n\,\Psi(R^L_0)}{1+\xi(R^L_0-1)} \qquad\qquad (25)$$

**$h=1$ to $4$ — Progressive credit and investment contraction.** Bank net worth erodes continuously. The leverage constraint forces a progressive contraction in lending:

$$\ell_t = \lambda N_t - b^B_t \qquad\qquad (26)$$

with $b^B_t \approx b^B_{-1}$ (no doom-loop accumulation). Credit contracts from $-2.25$ pp at $h=2$ to $-3.48$ pp at $h=4$.

**Fiscal and current-account dynamics.** Log-linearizing goods-market clearing:

$$\hat{y}_t = s_C \hat{c}_t + s_I \iota_t + s_G \breve{g}_t + s_{CA}\, ca_t \qquad\qquad (27)$$

Under non-default: $\hat{c}_0 < 0$, $\iota_0 < 0$, $\breve{g}_0 > 0$ (debt-financed expansion), $ca_0 > 0$ partially. The output loss deepens as capital declines.

---

## 5. The Default-Linked Path

Default occurs at $t=0$: $B_0 \ge B^*_0(K_0)$. International wholesale funding collapses. Three mechanisms operate simultaneously.

### 5.1 Sudden Stop

The balance sheet under autarky: $b^B_t + \ell_t = D_t + N_t$ with $F_t = 0$. The current-account identity gives:

$$CA_t = Y_t - C_t - I_t - G_t = F_{-1} + (1-\delta^B)B_0 \qquad\qquad (28)$$

The improvement combines **withdrawal of wholesale funding** $F_{-1}$ and **suspended debt service** $(1-\delta^B)B_0$. Both are forced immediately — a larger and more abrupt adjustment than the gradual non-default path, matching $+1.44$ pp at $h=0$ and $+2.61$ pp at $h=1$.

With $F_t = 0$, the deposit rate is set by domestic equilibrium (see §5.4):

$$C_0^{-\sigma} = \beta\,R^{D,\text{aut}}_0\,\mathbb{E}_0\!\left[C_1^{-\sigma}\right], \qquad R^{D,\text{aut}}_0 > R^{D,ss} \qquad\qquad (29)$$

### 5.2 Doom Loop: Gambling for Resurrection

Under autarky, expected returns from private lending $R^L_t(1-p^{RDL}) \to 0$. Banks optimally set $\ell_t = 0$ and saturate the leverage constraint with sovereign bonds when:

$$\mu\theta R^{BDK}_t + \gamma^B N_t \ge R^L_t(1-p^{RDL}) \qquad\qquad (30)$$

The net accumulation of sovereign bonds is:

$$\Delta b^B_0 = \lambda N_0 - b^B_{-1} > 0 \qquad \text{when} \qquad \lambda N_0 > b^B_{-1} \qquad\qquad (31)$$

generating the empirical $+4$ to $+6$ pp accumulation. Since private lending collapses to zero, **the private-credit channel is absent under default** — not because banks deleverage, but because they redirect balance-sheet capacity entirely toward sovereign bonds.

### 5.3 Fiscal Boom-Bust

The debt-service suspension provides a fiscal windfall:

$$G_0 = (1-\phi)\tau Y_0 + (1-\delta^B)B_0 \qquad\qquad (32)$$

permitting a large initial expansion ($+2.19$ pp at $h=0$). Post-restructuring conditionality requires $\tau Y_t - G_t \ge \rho B^\circ$, forcing retrenchment ($-1.39$ pp at $h=2$, $-1.69$ pp at $h=4$).

### 5.4 The Autarky Deposit Rate

Under default, $F_t = 0$ and post-restructuring $b^B = 0$, so the bank funding constraint becomes $\ell^{\text{aut}}_t = D^{\text{aut}}_t + N^{\text{aut}}_t$. With no sovereign-bond crowding, bank optimization gives $R^{L,\text{aut}}_t = R^{D,\text{aut}}_t$. The investment first-order condition (9) then pins the autarky rate:

$$R^{D,\text{aut}}_t = \alpha A (K^{\text{aut}}_t)^{\alpha-1} + (1-\delta) \qquad\qquad (33)$$

Since $I_t = 0$ during exclusion, the capital stock falls below its steady state ($K^{\text{aut}} < K^{ss}$), raising the marginal product. Hence:

$$R^{D,\text{aut}} = \alpha A (K^{\text{aut}})^{\alpha-1} + (1-\delta) > R^{D,ss} = R^* + \gamma s^{ss} \qquad\qquad (34)$$

**The autarky lending rate exceeds the sovereign-ceiling rate.** This is the wedge $\Delta r^{L,\text{aut}} > \gamma \hat{s}_0$ that makes default-path output losses larger than non-default on impact.

---

## 6. Propositions

**Proposition 1 (Non-Default Channel).** Following a spread crisis with no default, the spread transmits to output through three sequential mechanisms: the **sovereign-ceiling channel** (11)–(12) raises rates contemporaneously, contracting consumption via (3') and output via (8) at $h=0$; the **balance-sheet channel** (14)–(17') progressively erodes net worth in proportion to $b^B_{-1}$, tightening lending (26) and the user cost (9) at $h=1$–$4$; the current account adjusts gradually and the fiscal response is mildly expansionary until borrowing capacity is exhausted.

**Proposition 2 (Default-Linked Channel).** Following default, wholesale funding collapses simultaneously with debt-service suspension, generating a **two-component sudden stop** (28). The sovereign-ceiling channel is severed, but the autarky deposit rate rises sharply (29). Banks **gamble for resurrection** (30)–(31), accumulating $4$–$6$ pp of sovereign bonds without a conventional private-credit squeeze. The fiscal response is **boom-bust** (32).

**Corollary (The Central Paradox).** Non-default crises generate large output losses without default because the sovereign-ceiling channel transmits the spread regardless of repayment. Default episodes do not exhibit larger private-credit contractions despite the doom loop: banks redirect capacity toward sovereign bonds rather than contracting. Under non-default, the channel $\gamma s_t$ persists and deepens; under default, it is severed and replaced by the more acute sudden-stop mechanism.

---

## 7. Competitive Equilibrium

A competitive equilibrium consists of sequences

$$\{C_t,\, I_t,\, K_t,\, D_t,\, F_t,\, N_t,\, b^B_t,\, \ell_t,\, G_t,\, B_t,\, s_t,\, q_t,\, R^L_t,\, R^D_t,\, R^O_t,\, w_t\}$$

such that: (i) households satisfy (1)–(3); (ii) firms satisfy (7)–(9); (iii) banks optimize subject to (10) and (16) with costs (11)–(12) and lending rate (17'); (iv) the government satisfies (21) and (24); (v) foreign investors satisfy (18), yielding (19); (vi) labor, deposit, and wholesale funding markets clear; (vii) goods market clears:

$$Y_t = C_t + I_t + G_t + CA_t \qquad\qquad (39)$$

In the **non-crisis steady state**: $\beta(R^* + \gamma s^{ss}) = 1$, output $Y^{ss} = A(K^{ss})^{\alpha}\Psi(R^{L,ss})$ with $R^{L,ss} = R^* + \gamma s^{ss} + \varphi\, b^{B,ss}/N^{ss}$, and $\alpha A(K^{ss})^{\alpha-1} = R^{L,ss}-(1-\delta)$.

---

## 8. Log-Linearization and Analytical Solution

### 8.1 Conventions

For variable $X_t$: $\hat{x}_t \equiv \log(X_t/X^{ss})$. For the spread: $\hat{s}_t \equiv s_t - s^{ss}$ (level deviation). Spread AR(1):

$$\hat{s}_t = \rho_s \hat{s}_{t-1} + \varepsilon_t, \qquad \rho_s \in (0,1) \qquad\qquad \text{(LL.0)}$$

### 8.2 Log-Linearized Equations

**Lending rate** — with $\Omega \equiv \varphi\, b^{B,ss}/N^{ss}$, non-default path ($\hat{n}^{-1}_{t-1}=0$):

$$\Delta r^{L,nd}_t = \gamma \hat{s}_t - \Omega \acute{n}_t \qquad\qquad \text{(LL.1a)}$$

**Output** — with $\varepsilon_p \equiv \xi n R^{L,ss}/[1+\xi(R^{L,ss}-1)]$:

$$\hat{y}^{nd}_t = \alpha \hat{k}_t - \varepsilon_p\!\left(\gamma \hat{s}_t - \Omega \acute{n}_t\right) \qquad\qquad \text{(LL.2a)}$$

**Capital** — with $\eta \equiv 1/[(1-\alpha)(R^{L,ss}-(1-\delta))]$:

$$\hat{k}^{nd}_{t+1} = (1-\delta)\hat{k}_t - \eta\!\left(\gamma \hat{s}_t - \Omega \acute{n}_t\right) \qquad\qquad \text{(LL.3a)}$$

**Bank net worth** — with $\mathrm{B} \equiv b^{B,ss}/[N^{ss}(1+s^{ss})^2]$:

$$\acute{n}^{nd}_{t+1} = \Phi^N \acute{n}_t - \mathrm{B}\,\hat{s}_{t+1} - \Phi^O \gamma \hat{s}_t \qquad\qquad \text{(LL.5a)}$$

Driven down by the bond-price loss ($\mathrm{B}\hat{s}_{t+1}$, balance-sheet channel) and higher wholesale funding costs ($\Phi^O \gamma \hat{s}_t$, sovereign-ceiling channel). Stability requires $|\Phi^N| < 1$.

**Private credit:**

$$\hat{\ell}^{nd}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\,\acute{n}_t \qquad\qquad \text{(LL.6)}$$

**Consumption:**

$$\hat{c}_t = \mathbb{E}_t[\hat{c}_{t+1}] - \frac{\gamma}{\sigma\, R^{D,ss}}\,\hat{s}_{t+1} \qquad\qquad \text{(LL.7)}$$

### 8.3 Closed-Form Impulse Responses (Non-Default)

Iterating LL.5a from $\acute{n}_0 = -\mathrm{B}\hat{s}_0$ with $\hat{s}_t = \rho_s^t \hat{s}_0$:

$$\acute{n}^{nd}_h = (\Phi^N)^h \acute{n}_0 \;-\; \sum_{h'=1}^{h} (\Phi^N)^{h-h'}\!\left(\mathrm{B}\,\rho_s^{h'} + \Phi^O \gamma\,\rho_s^{h'-1}\right)\hat{s}_0 \qquad\qquad \text{(IR.1)}$$

> **Note.** An earlier draft factored this as $-(\mathrm{B}+\Phi^O\gamma)\sum(\Phi^N)^{h-h'}\rho_s^{h'}\hat{s}_0$, which holds only when $\rho_s=1$. The correct form keeps the two sums distinct: $\mathrm{B}$ multiplies $\rho_s^{h'}$ (the spread at period $h'$) while $\Phi^O\gamma$ multiplies $\rho_s^{h'-1}$ (one period earlier). Qualitative results are unaffected. In `16_model_irf.do` we iterate LL.5a directly and avoid the summation.

Iterating LL.3a from $\hat{k}_0 = 0$:

$$\hat{k}^{nd}_h = -\eta\gamma \sum_{h'=1}^{h}(1-\delta)^{h-h'}\rho_s^{h'-1}\hat{s}_0 \;+\; \eta\Omega \sum_{h'=1}^{h}(1-\delta)^{h-h'}\acute{n}_{h'-1} \qquad\qquad \text{(IR.2)}$$

Output at horizon $h$:

$$\hat{y}^{nd}_h = \underbrace{\alpha \hat{k}^{nd}_h}_{\text{capital}} \;-\; \underbrace{\varepsilon_p\gamma\rho_s^h \hat{s}_0}_{\text{working capital}} \;+\; \underbrace{\varepsilon_p\Omega\,\acute{n}^{nd}_h}_{\text{balance sheet}} \qquad\qquad \text{(IR.3)}$$

The capital channel is zero at $h=0$ and grows thereafter; the working-capital channel is largest at $h=0$ and decays at rate $\rho_s$; the balance-sheet term deepens as net worth erodes. The sum is negative at all horizons:

$$\hat{y}^{nd}_h < 0 \quad \forall\, h \ge 0 \quad \text{whenever} \quad \hat{s}_0 > 0 \qquad\qquad \text{(IR.3a)}$$

Private credit:

$$\hat{\ell}^{nd}_h = \frac{\lambda N^{ss}}{\ell^{ss}}\,\acute{n}^{nd}_h \qquad\qquad \text{(IR.4)}$$

Small at $h=0$, deepening to $-3.48$ pp at $h=4$.

**Lemma (Persistence Without Mean Reversion).** Let $\rho_s, \Phi^N \in (0,1)$. Then (i) $\hat{y}^{nd}_h < 0$ for all $h \ge 0$; (ii) the trough occurs at finite $h^* > 0$; (iii) the cumulative output loss diverges as $\rho_s \to 1$.

### 8.4 The Default Path

**Break 1 — Wholesale funding collapse.** $F_t = 0$ replaces LL.7 with the autarky Euler equation using $R^{D,\text{aut}} > R^{D,ss}$ (§5.4), generating a larger front-loaded consumption contraction.

**Break 2 — Gambling for resurrection.** When (30) holds, $\hat{\ell}^{\text{def}}_t \approx -\infty$ while $\hat{\ell}^{nd}_t \in (-\infty, 0)$.

> **Proposition 3 (Credit Paradox).** For all $h \ge 0$, private credit is lower under default than non-default, yet the *gross* bank balance sheet does not shrink — banks redirect the full capacity $\lambda N_t$ toward sovereign bonds. The econometrician observes bank credit to government rising while private credit falls, with aggregate bank assets roughly unchanged. This is why default-linked episodes show **no statistically larger private-credit contraction** despite the doom loop.

**Break 3 — Doom-loop accumulation:**

$$\hat{n}^{B,\text{def}}_t = \left(\frac{\lambda N^{ss}}{b^{B,ss}} - 1\right) - \frac{\lambda N^{ss}}{b^{B,ss}}\,\mathrm{B}\,\hat{s}_0 \qquad\qquad \text{(LL.10)}$$

positive when the reallocation gain exceeds net-worth erosion.

**Output under default.** With $I_t=0$, capital depletes at rate $\delta$, so $K^{\text{def}}_h = (1-\delta)^h K_0$ and $\hat{k}^{\text{def}}_h = h\ln(1-\delta) < 0$. Output:

$$\hat{y}^{\text{def}}_h = \alpha \cdot h\ln(1-\delta) - \varepsilon_p\,\Delta r^{L,\text{aut}} \qquad\qquad \text{(IR.5)}$$

> **Note.** An earlier draft wrote the capital term as $\alpha(1-\delta)^h \hat{k}_0 = 0$ (since $\hat{k}_0=0$), which is identically zero and contradicts the stated deepening of losses. The correct form uses the cumulative depreciation $h\ln(1-\delta)$, which is negative and growing in magnitude.

The **autarky wedge** $\Delta r^{L,\text{aut}}$ is pinned by the $h=0$ default data point; $h=1$–$4$ are then zero-free-parameter out-of-sample predictions (see §10.3).

---

## 9. Mapping to the Empirical Findings

### Table 1 — Model Sign Predictions by Channel, Path, and Horizon

| Variable | Path | $h=0$ | $h=1$–$2$ | $h=3$–$4$ | Mechanism |
|---|---|:---:|:---:|:---:|---|
| **Output** | Non-default | − | − (deepens) | −/≈0 | working capital + balance sheet |
| | Default | −− | −− | −− | autarky working capital + capital depletion |
| | *Comparison* | nd > d | nd > d | nd > d | default always deeper |
| **Private credit** | Non-default | − (small) | − (growing) | − (deep) | net-worth channel, IR.4 |
| | Default | ≈0 | ≈0 | ≈0 | banks redirect to sovereign bonds |
| | *Comparison* | nd < d | nd < d | nd < d | nd has worse measured credit |
| **Investment** | Non-default | − | − | −/≈0 | user cost rises, eq. (9) |
| | Default | ≈0 | − | −− | $I=0$, capital depletes |
| | *Comparison* | nd > d | nd > d | nd > d | |
| **Claims on govt** | Non-default | ≈0 | ≈0 | ≈0 | no doom-loop incentive |
| | Default | + | ++ | + | gambling for resurrection, eq. (31) |
| | *Comparison* | nd < d | nd < d | nd < d | +4–6 pp accumulation |
| **Current account** | Non-default | + (gradual) | + (fading) | ≈0 | wholesale cost rises, flows continue |
| | Default | ++ | ++ | + | two-component sudden stop, eq. (28) |
| | *Comparison* | nd < d | nd < d | nd ≈ d | default front-loaded |
| **Govt expenditure** | Non-default | + | + | + | retained market access |
| | Default | ++ | − | −− | boom-bust: windfall then conditionality |
| | *Comparison* | nd < d | nd > d | nd > d | paths cross near $h\approx1$ |

---

### Table 2 — Empirical Findings, Model Equations, and Parameter Conditions

| Finding | Model equation(s) | Parameter condition |
|---|---|---|
| **(a)** ND output loss (−1.86 pp at $h=0$, −3.76 pp at $h=1$) | LL.2a, IR.3 | $\varepsilon_p>0$ (req. $\xi>0$), $\gamma>0$, $\hat{s}_0>0$ |
| **(b)** ND trough at $h=1$, not $h=0$ | IR.3 + IR.2 | $\Phi^N\in(0,1)$, $\Omega>0$: net-worth deepening overtakes spread decay at $h=1$ |
| **(c)** Default deeper and front-loaded (−3.42 pp at $h=0$, −5.50 pp at $h=1$) | IR.5, eq. (34) | $\Delta r^{L,\text{aut}}>\gamma\hat{s}_0$; guaranteed by $K^{\text{aut}}<K^{ss}$ |
| **(d)** Credit contraction in ND only (−4.78 pp at $h=4$; ≈0 default) | IR.4, Prop. 3 | ND: $\acute{n}<0$; default: GFR condition (30) $\Rightarrow \ell=0$ |
| **(e)** Doom loop +4–6 pp in default only | eqs. (30)–(31), LL.10 | $\lambda N_0 > b^B_{-1}$; $\mu\theta R^{BDK}+\gamma^B N \ge R^L(1-p^{RDL})$ |
| **(f)** Diabolic loop absorbed by lagged $b^B$ | eqs. (14)–(15), (17') | Amplification $\propto b^B_{t-1}$: controlling $L.b^B$ absorbs $\mathrm{B}$ in LL.5a |
| **(g)** CA larger and earlier in default (+2.61 pp at $h=1$) | eq. (28) vs. eq. (27) | $F_{-1}>0$, $\delta^B>0$: sudden stop >> gradual adjustment |
| **(h)** Govexp boom-bust in default (+2.19 pp at $h=0$, −1.69 pp at $h=4$) | eq. (32) | $\delta^B>0$ (windfall); exclusion $>2$ periods (conditionality bites) |
| **(i)** Govexp sustained in ND (positive throughout) | eq. (21) + ceiling (11) | $\gamma<\bar{\gamma}$: spread rise insufficient to trigger solvency crisis |

**Theorem (Joint Consistency).** Under the conditions in Table 2, the log-linearized model under both paths jointly reproduces the qualitative and directional predictions of all nine empirical findings. No individual condition contradicts another; there exists an open set $\Theta \subset \mathbb{R}^{11}_+$ under which all findings hold simultaneously.

---

## 10. Calibration and Numerical Solution

### 10.1 Parameter Groups

**Group 1 — Literature values.** $\beta=0.96$, $\sigma=2$, $\alpha=0.33$, $\delta=0.10$, $R^*=0.04$; recovery $\theta=0.62$ (Cruces-Trebesch 2013); re-entry $\mu=0.22$ (Gelos et al. 2011, mean exclusion 4.5 years); maturing fraction $\delta^B=0.22$; leverage $\lambda=10$.

**Group 2 — Data steady-state ratios.** Computed from *tranquil-period* means (non-crisis, non-continuation) in the panel: background spread $s^{ss}$; expenditure shares $s_C, s_I, s_G, s_{CA}$; tax rate $\tau$; bank sovereign exposure $b^{B,ss}/Y$; net worth $N^{ss}/Y$; and hence $\varepsilon_p$, $\eta$, $\mathrm{B}$, $\Omega$. Implemented in `do/14_calibration.do`.

**Group 3 — Estimated processes.** $\rho_s$ from a panel AR(1) on EMBIG spreads; $(\rho_y, \sigma_y)$ from a detrended AR(1) on log GDP per capita. Both estimated within `do/14_calibration.do`.

The deep transmission parameters $\{\xi, \varphi, \Phi^N\}$ are calibrated by **Simulated Method of Moments to the non-default IRF only**.

### 10.2 Two-Block Solution

**Block 1 — Nonlinear default block** (`do/15_solve_default.do`). Canonical Arellano (2008) VFI over states $(B, y)$ in Mata: Tauchen-discretized income, the default decision (24), and the bond-price fixed point (18)–(19) with recovery $\theta$ and re-entry $\mu$. Simulation of 100,000 periods yields model moments (default frequency, mean/SD spread, debt/GDP) and a **crisis-event spread path** ($h=-2$ to $4$) that feeds Block 2.

**Block 2 — Log-linear transmission block** (`do/16_model_irf.do`). The model spread path enters LL.5a, LL.3a, LL.2a, and IR.4 (iterated directly), producing model IRFs for output, net worth, capital, and credit. SMM grid search over $\{\xi, \varphi, \Phi^N\}$ minimizes squared deviations from the empirical non-default LP at $h=0$–$4$:

$$\min_{\xi,\,\varphi,\,\Phi^N} \sum_{h=0}^{4} \left[\left(\hat{y}^{nd,\text{model}}_h - \hat{y}^{nd,\text{LP}}_h\right)^2 + \left(\hat{\ell}^{nd,\text{model}}_h - \hat{\ell}^{nd,\text{LP}}_h\right)^2\right]$$

### 10.3 Identification and Out-of-Sample Validation

$\{\xi, \varphi, \Phi^N\}$ are pinned entirely by the **non-default IRF**. The **default path is then validated out of sample**: the single autarky wedge $\Delta r^{L,\text{aut}}$ is pinned to the $h=0$ default data point, after which $h=1$–$4$ are zero-free-parameter predictions from IR.5. The degree to which these match the empirical default LP constitutes the model's central test.

---

## 11. Scope and Limitations

The sovereign-ceiling pass-through $\gamma$ and wholesale-funding share $\omega$ are constants; in practice both vary with financial integration and the currency composition of bank liabilities. The default decision abstracts from renegotiation dynamics and long-maturity debt. TFP $A$ is constant, abstracting from the productivity feedback in Mendoza-Yue (2012). These simplifications serve analytical clarity rather than quantitative discipline. A full quantitative extension would combine a Bocola (2016)-type balance-sheet model for the non-default channel with a Gennaioli-Martin-Rossi (2014) framework for the default-linked channel.

---

## Notation Summary

| Symbol | Meaning | Symbol | Meaning |
|---|---|---|---|
| $C_t$, $D_t$ | consumption, deposits | $s_t$ | EMBIG spread |
| $K_t$, $I_t$ | capital, investment | $q_t$ | bond price $= 1/(1+s_t)$ |
| $\ell_t$, $b^B_t$ | bank lending, sovereign bond holdings | $\pi_{t+1}$ | default probability |
| $N_t$, $F_t$ | bank net worth, wholesale funding | $\theta$, $\mu$ | recovery rate, re-entry prob. |
| $R^D_t$, $R^O_t$, $R^L_t$ | deposit, wholesale, lending rates | $\gamma$, $\varphi$ | ceiling pass-through, balance-sheet pass-through |
| $G_t$, $B_t$ | govt expenditure, external debt | $\xi$, $\lambda$ | working-capital share, bank leverage |
| $\Omega$, $\mathrm{B}$ | amplification, balance-sheet sensitivity | $\Phi^N$, $\Phi^O$ | net-worth multiplier, wholesale income share |
| $\varepsilon_p$, $\eta$ | output and capital elasticities to $R^L$ | $\rho_s$, $\rho_y$ | spread and income persistence |

---

*Empirical implementation: `do/14_calibration.do` (calibration), `do/15_solve_default.do` (nonlinear default block, Arellano VFI in Mata), `do/16_model_irf.do` (log-linear transmission, SMM, model-vs-data overlay figure).*
