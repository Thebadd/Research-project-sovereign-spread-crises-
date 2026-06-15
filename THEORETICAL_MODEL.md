# A Stylized Model of Sovereign Spread Crises With and Without Default

**Theoretical Appendix**

*Companion to "The Aftermath of Sovereign Spread Crises With and Without Default"*

---

## Overview

This appendix develops a small open economy model that rationalizes the central empirical findings of the paper: (i) sovereign spread crises generate large and persistent output losses *even when no default occurs*; (ii) default-linked crises are roughly twice as costly and front-loaded; (iii) the private-credit and investment contractions are concentrated in *non-default* episodes; (iv) bank accumulation of sovereign debt (the "doom loop") is exclusive to *default-linked* episodes; and (v) the current-account adjustment is larger and earlier under default.

The model's defining feature is a **separation between the international sovereign bond market**, where the government issues dollar-denominated external debt priced by foreign investors, **and the domestic credit market**, where banks lend to firms in local currency. Transmission from the former to the latter runs through two channels: a **wholesale funding-cost channel** (the sovereign ceiling), and a **balance-sheet channel** (bond-price losses erode bank net worth). This dual structure is what allows the model to deliver both a non-default transmission mechanism and a distinct default mechanism.

The solution strategy is a **two-block decomposition** following Bocola (2016): a nonlinear sovereign-risk block that determines the endogenous default decision and the equilibrium spread (solved by value-function iteration), and a log-linear transmission block that maps the spread into output, credit, investment, and the current account (solved in closed form). Nonlinearity is retained where it is essential — the default choice — and linearized where it is not.

---

## 1. Environment and Agents

Consider a small open economy in discrete time, $t = 0, 1, 2, \dots$ The economy consists of five agents: a representative household, a representative firm, a domestic banking sector with access to both domestic and international wholesale funding, a government, and a continuum of competitive foreign investors active in the international sovereign bond market.

### 1.1 Households

The representative household maximizes lifetime utility over consumption $C_t$:

$$\max \; \mathbb{E}_0 \sum_{t=0}^{\infty} \beta^t \, u(C_t), \qquad u(C_t) = \frac{C_t^{1-\sigma}}{1-\sigma} \tag{1}$$

where $\beta \in (0,1)$ is the discount factor and $\sigma > 0$ is the inverse of the intertemporal elasticity of substitution. The household holds one-period domestic bank deposits $D_t$ at gross return $R^D_t$, earns labor income $w_t$, and receives lump-sum bank profits $\Pi_t$. The period budget constraint is:

$$C_t + D_{t+1} = R^D_t D_t + w_t + \Pi_t \tag{2}$$

The first-order condition yields the consumption Euler equation:

$$C_t^{-\sigma} = \beta R^D_{t+1} \, \mathbb{E}_t\!\left[C_{t+1}^{-\sigma}\right] \tag{3}$$

The household holds only domestic deposits and has no direct access to international capital markets, consistent with the financial-development constraints characteristic of emerging market economies. The deposit rate $R^D_{t+1}$ is determined in the domestic banking equilibrium (Section 1.3).

### 1.2 Firms

A representative competitive firm produces with Cobb-Douglas technology:

$$Y_t = A K_t^{\alpha} L_t^{1-\alpha}, \qquad \alpha \in (0,1) \tag{4}$$

where $A>0$ is total factor productivity and $L_t = 1$ is inelastic labor supply. Capital depreciates at rate $\delta \in (0,1)$:

$$K_{t+1} = (1-\delta)K_t + I_t \tag{5}$$

Firms face a **working-capital constraint**: a fraction $\xi \in (0,1)$ of the wage bill must be pre-financed by short-term bank loans before production. Total firm credit demand is:

$$\ell_t = \ell_t^w + \ell_t^I = \xi w_t + I_t \tag{6}$$

The working-capital cost raises the effective cost of labor from $w_t$ to $[1+\xi(R^L_t-1)]w_t$. Profit maximization yields:

$$(1-\alpha)A K_t^{\alpha} L_t^{-\alpha} = [1+\xi(R^L_t-1)]\,w_t \tag{7}$$

and equilibrium output as a function of the domestic lending rate:

$$Y_t = A K_t^{\alpha} \cdot \Psi(R^L_t), \qquad \Psi(R^L_t) \equiv \left[\frac{1-\alpha}{1+\xi(R^L_t-1)}\right]^{n}, \quad n=\frac{1-\alpha}{\alpha} \tag{8}$$

with $\Psi'(R^L_t)<0$: output is strictly decreasing in the lending rate for any $\xi>0$. **Equation (8) is the key result of the firm block**: a rise in $R^L_t$ from any upstream shock generates an *immediate within-period output contraction*, over and above the deferred loss from reduced capital accumulation. Investment demand equates the marginal product of capital to the user cost:

$$\alpha A K_{t+1}^{\alpha-1} = R^L_t - (1-\delta) \tag{9}$$

### 1.3 Banking Sector with Dual Funding Structure

The banking sector is the central transmission mechanism. Banks hold sovereign bonds priced in international markets, so their balance-sheet values depend on foreign investor sentiment; and they fund through a dual structure of domestic deposits and international wholesale funding.

#### 1.3.1 Balance Sheet

$$b^B_t + \ell_t = D_t + F_t + N_t \tag{10}$$

where $b^B_t$ is holdings of dollar-denominated sovereign bonds (marked to market at international price $q_t$), $\ell_t$ is total domestic lending, $D_t$ is deposits, $F_t$ is international wholesale funding, and $N_t$ is bank net worth.

#### 1.3.2 International Wholesale Funding and the Sovereign Ceiling

Foreign lenders charge domestic banks a gross funding rate reflecting both the risk-free rate $R^*$ and the sovereign's credit risk — the **sovereign ceiling**:

$$R^O_t = R^* + \gamma s_t, \qquad \gamma \in (0,1] \tag{11}$$

where $s_t$ is the sovereign's EMBIG spread (eq. 19) and $\gamma$ is the pass-through parameter. When $\gamma=1$ the ceiling binds fully; $\gamma<1$ allows banks to partially escape sovereign risk. In the interior equilibrium where both funding sources are active, banks equate marginal funding costs:

$$R^D_t = R^O_t = R^* + \gamma s_t \tag{12}$$

**Equation (12) pins the domestic deposit rate to the international spread.** Substituting into the Euler equation:

$$C_t^{-\sigma} = \beta\,(R^* + \gamma s_{t+1})\,\mathbb{E}_t\!\left[C_{t+1}^{-\sigma}\right] \tag{3$'$}$$

A spread spike raises the rate at which the household discounts future consumption, generating a contemporaneous consumption contraction *before any balance-sheet or credit-quantity effect operates*. This channel is absent from models that treat the spread purely as a government borrowing cost.

#### 1.3.3 Net Worth and the Balance-Sheet Channel

$$N_{t+1} = R^L_t \ell_t + (R^B_t + \Delta q_{t+1}) b^B_t - R^D_t D_t - R^O_t F_t \tag{13}$$

A spread increase generates an immediate capital loss on outstanding bonds through $q_t = 1/(1+s_t)$:

$$\Delta q_t \approx -\frac{\Delta s_t}{(1+s_t)^2} < 0 \tag{14}$$

so the balance-sheet shock to net worth is:

$$\Delta N_t = b^B_{t-1}\cdot \Delta q_t = -\,\frac{b^B_{t-1}\,\Delta s_t}{(1+s_t)^2} < 0 \tag{15}$$

**proportional to pre-existing exposure $b^B_{t-1}$.** The two channels are additive: the funding-cost channel (11)–(12) raises the *price* of credit; the balance-sheet channel (15) reduces the *quantity* banks can supply.

#### 1.3.4 Leverage Constraint and the Lending Rate

Banks face a leverage constraint:

$$b^B_t + \ell_t \le \lambda N_t, \qquad \lambda > 1 \tag{16}$$

The domestic lending rate combines weighted funding costs and a balance-sheet risk premium:

$$R^L_t = \omega R^D_t + (1-\omega) R^O_t + \varphi\,\frac{b^B_{t-1}}{N_t} \tag{17}$$

Substituting (11)–(12):

$$\boxed{\,R^L_t = R^* + \gamma s_t + \varphi\,\frac{b^B_{t-1}}{N_t}\,} \tag{17$'$}$$

**Equation (17$'$) is the central transmission equation.** It decomposes the lending rate into: a baseline $R^*$; a **sovereign-ceiling term** $\gamma s_t$ affecting all borrowers contemporaneously; and a **balance-sheet term** $\varphi(b^B_{t-1}/N_t)$ that amplifies the spread effect in proportion to pre-existing sovereign exposure.

> **Remark on identification.** Equation (17$'$) explains why the empirical finding that the doom-loop channel is *absorbed by the lagged stock of bank claims on government* is structurally meaningful. The sovereign-ceiling term $\gamma s_t$ affects all banks equally and is captured by time fixed effects in the local projections. The amplification term $\varphi(b^B_{t-1}/N_t)$ is proportional to the *lagged* sovereign position $b^B_{t-1}$ — exactly the control included in the empirical specification to isolate crisis-period accumulation.

---

## 2. International Sovereign Bond Market and the Endogenous Spread

### 2.1 Foreign Investor Problem

A continuum of risk-neutral foreign investors purchase external bonds at price $q_t$. Their zero-profit participation constraint requires:

$$q_t(1+R^*) = \mathbb{E}_t\!\left[(1-\pi_{t+1})\cdot 1 + \pi_{t+1}\cdot\theta\right] \tag{18}$$

where $\pi_{t+1} \in [0,1]$ is the rational expectation of default and $\theta \in [0,1)$ is the recovery rate.

### 2.2 Equilibrium Spread

Using $q_t = 1/(1+s_t)$:

$$s_t = \frac{(1+R^*)(1-\theta)\pi_{t+1}}{1 - \pi_{t+1}(1-\theta)} \tag{19}$$

The spread is strictly increasing in $\pi_{t+1}$ and decreasing in $\theta$. For small $\pi_{t+1}$, $s_t \approx (1+R^*)(1-\theta)\pi_{t+1}$. **$s_t$ is exclusively the sovereign's external bond spread** — not the domestic lending-deposit spread. Transmission to domestic credit runs entirely through the banking block, equations (11)–(12) and (15)–(17$'$).

### 2.3 Spread Crisis Definition

A spread crisis is an episode where the EMBIG spread crosses the threshold $\bar s = 1000$ bps from below:

$$\Omega_t = \mathbf{1}\!\left[s_t \ge \bar s\right] \tag{20}$$

consistent with the Pescatori-Sy (2007) criterion used empirically. By inverting (19), $s_t \ge \bar s$ iff the expected default probability exceeds a threshold $\bar\pi$.

---

## 3. Government and the Default Decision

### 3.1 Budget Constraint

$$G_t + (1-\delta^B)B_t = \tau Y_t + q_t\!\left[B_{t+1} - (1-\delta^B)B_t\right] \tag{21}$$

where $G_t$ is primary expenditure, $B_t$ is the face value of external debt, $\delta^B \in (0,1)$ is the fraction maturing each period, $\tau Y_t$ is tax revenue, and $q_t$ is the bond price from (19). The government faces an expenditure floor $\bar G$ and a borrowing limit $B_{t+1} \le \bar B(s_t)$ that tightens as the spread rises.

### 3.2 Recursive Default Problem

With state $(B_t, K_t)$, the value of repayment is:

$$V^R_t(B_t,K_t) = \max_{G_t, B_{t+1}}\Big\{\, U(G_t) + \beta_p\,\mathbb{E}_t[V_{t+1}(B_{t+1},K_{t+1})]\,\Big\} \tag{22}$$

subject to (21), (5), $G_t \ge \bar G$, $B_{t+1} \le \bar B(s_t)$, where $\beta_p$ is the government discount factor and $V_{t+1}=\max\{V^R_{t+1},V^D_{t+1}\}$. Upon default the government is excluded from markets, re-entering with probability $\mu$ each period:

$$V^D_t(K_t) = \max_{G_t}\Big\{\, U(G_t) + \beta_p\,\mathbb{E}_t\!\left[\mu V^R_{t+1}(0,K_{t+1}) + (1-\mu)V^D_{t+1}(K_{t+1})\right]\Big\} \tag{23}$$

subject to the autarky budget $G_t = (1-\phi)\tau Y_t$, where $\phi$ is the direct output cost of market exclusion. Capital accumulation stalls during autarky ($I_t = 0$). Re-entry resets debt to the post-haircut level $\theta B_{t-1}$.

### 3.3 Default Threshold and Equilibrium

The government defaults when $V^D_t(K_t) \ge V^R_t(B_t,K_t)$, defining a default set with threshold $B^*_t(K_t)$ decreasing in $s_t$ and $K_t$. Foreign investors form rational expectations:

$$\pi_{t+1} = \Pr\!\left(B_{t+1} \ge B^*_{t+1}(K_{t+1}) \mid \mathcal F_t\right) \tag{24}$$

Equations (19) and (24) form the **fixed-point system** that closes the model.

---

## 4. The Non-Default Path: Transmission Mechanisms

Suppose the government repays in all periods. A spread crisis occurs at $t=0$ because investors revise $\pi_{t+1}$ upward, driving $s_0$ past $\bar s$.

### 4.1 Impact: Sovereign Ceiling on Consumption ($h=0$)

The spread spike raises $R^O_0 = R^* + \gamma s_0$ and, by arbitrage, $R^D_0 = R^* + \gamma s_0$. From (3$'$), the household defers consumption — a contemporaneous demand contraction *before any credit-quantity effect*.

### 4.2 Impact: Sovereign Ceiling on Working Capital ($h=0$)

The lending rate rises at onset through (17$'$). From (8), the immediate output contraction is:

$$\frac{dY_0}{dR^L_0} = A K_0^{\alpha}\,\Psi'(R^L_0) < 0, \qquad \Psi'(R^L_0) = -\frac{\xi n\,\Psi(R^L_0)}{1+\xi(R^L_0-1)} \tag{25}$$

### 4.3 Progressive: Credit and Investment ($h=1$ to $h=4$)

As the spread persists, net worth erodes (eq. 13). The leverage constraint forces a progressive contraction in lending:

$$\ell_t = \lambda N_t - b^B_t \tag{26}$$

with $b^B_t \approx b^B_{-1}$ (no doom-loop accumulation absent gambling-for-resurrection). Credit contracts gradually over $h=1$ to $h=4$, matching the empirical $-2.25$ pp at $h=2$ deepening to $-3.48$ pp at $h=4$.

### 4.4 Fiscal Response and Current Account

The government issues additional bonds, exploiting residual borrowing capacity; as the spread persists, expenditure declines toward $\bar G$. The current account adjusts *partially* — wholesale funding contracts but the government keeps rolling its bonds, so the adjustment is gradual. Log-linearizing goods-market clearing (eq. 39):

$$\hat y_t = s_C \hat c_t + s_I \iota_t + s_G \breve g_t + s_{CA}\, ca_t \tag{27}$$

Under non-default: $\hat c_0<0$ (sovereign ceiling), $\iota_0<0$ (working capital + user cost), $\breve g_0>0$ (debt-financed expansion), $ca_0>0$ partially. The output loss begins on impact and deepens as capital declines.

---

## 5. The Default-Linked Path: Transmission Mechanisms

Suppose the government defaults at $t=0$: $B_0 \ge B^*_0(K_0)$. Default severs the sovereign-ceiling link, and international wholesale funding collapses.

### 5.1 Sudden Stop: Collapse of Wholesale Funding

The balance sheet under autarky becomes $b^B_t + \ell_t = D_t + N_t$ with $F_t = 0$. From the current-account identity:

$$CA_t = Y_t - C_t - I_t - G_t = F_{-1} + (1-\delta^B)B_0 \tag{28}$$

The current-account improvement combines **two forced components**: withdrawal of wholesale bank funding $F_{-1}$ and suspended debt service $(1-\delta^B)B_0$. This two-component sudden stop is larger and more abrupt than the gradual non-default adjustment, matching the empirical $+1.44$ pp at $h=0$ and $+2.61$ pp at $h=1$. With $F_t=0$, the deposit rate is set by domestic equilibrium (see §5.4):

$$C_0^{-\sigma} = \beta\,R^{D,\text{aut}}_0\,\mathbb{E}_0\!\left[C_1^{-\sigma}\right], \qquad R^{D,\text{aut}}_0 > R^{D,ss} \tag{29}$$

### 5.2 Doom Loop: Gambling for Resurrection

Under autarky, private loans carry a high non-performance probability $p^{RDL}$, so expected returns from private lending $R^L_t(1-p^{RDL})$ approach zero. The expected return from domestic sovereign bonds is $\mathbb{E}_t[r^{BDK}] = \mu\theta R^{BDK}_t$. The bank optimally sets $\ell_t=0$ and saturates the leverage constraint with sovereign bonds when:

$$\mu\theta R^{BDK}_t + \gamma^B N_t \ge R^L_t(1-p^{RDL}) \tag{30}$$

where $\gamma^B$ is the implicit guarantee subsidy. The net accumulation is:

$$\Delta b^B_0 = \lambda N_0 - b^B_{-1} > 0 \quad\text{when}\quad \lambda N_0 > b^B_{-1} \tag{31}$$

generating the empirical $+4$ to $+6$ pp accumulation in default episodes. **Since private lending collapses to zero, the private-credit channel to firms is absent under default** — not because banks deleverage, but because they redirect balance-sheet capacity entirely toward sovereign bonds.

### 5.3 Fiscal Boom-Bust

The debt-service suspension provides a fiscal windfall:

$$G_0 = (1-\phi)\tau Y_0 + (1-\delta^B)B_0 \tag{32}$$

permitting a large initial expansion ($+2.19$ pp at $h=0$). Post-restructuring, the government re-enters with $B^\circ = \theta B_0$ subject to conditionality requiring a primary surplus $\tau Y_t - G_t \ge \rho B^\circ$, forcing retrenchment ($-1.39$ pp at $h=2$, $-1.69$ pp at $h=4$).

### 5.4 The Autarky Deposit Rate

Under default, $F_t = 0$ and post-restructuring $b^B = 0$, so the bank funding constraint is $\ell^{\text{aut}}_t = D^{\text{aut}}_t + N^{\text{aut}}_t$. With no sovereign-bond crowding ($\varphi\cdot 0/N = 0$), bank optimization gives $R^{L,\text{aut}}_t = R^{D,\text{aut}}_t$. The autarky rate then solves the closed-economy version of the investment first-order condition (9):

$$R^{D,\text{aut}}_t = \alpha A (K^{\text{aut}}_t)^{\alpha-1} + (1-\delta) \tag{33}$$

Because investment is suspended during exclusion, the capital stock falls ($K^{\text{aut}} < K^{ss}$), raising the marginal product of capital. Hence:

$$R^{D,\text{aut}}_{ss} = \alpha A (K^{\text{aut}})^{\alpha-1} + (1-\delta) > R^{D,ss} = R^* + \gamma s^{ss} \tag{34}$$

**The autarky lending rate exceeds the sovereign-ceiling rate.** This is precisely the wedge $\Delta r^{L,\text{aut}} > \gamma \hat s_0$ that makes default-path output losses larger than non-default at impact.

---

## 6. Propositions

**Proposition 1 (Non-Default Channel).** Following a spread crisis with $\Omega_t=1$ and no default, the spread transmits to output through three sequential mechanisms: the **sovereign-ceiling channel** (11)–(12) raises deposit and lending rates contemporaneously, contracting consumption via (3$'$) and output via (8) at $h=0$; the **balance-sheet channel** (14)–(17$'$) progressively erodes net worth in proportion to $b^B_{-1}$, tightening lending (26) and raising the user cost of investment (9) at $h=1$–$4$; the current account adjusts gradually and the fiscal response is mildly expansionary until borrowing capacity is exhausted.

**Proposition 2 (Default-Linked Channel).** Following a spread crisis with default at $t=0$, wholesale funding collapses simultaneously with the debt-service suspension, generating a **two-component sudden stop** (28) and a large, immediate current-account adjustment. The sovereign-ceiling channel is severed, but the autarky deposit rate rises sharply (29), inducing front-loaded consumption contraction. Banks **gamble for resurrection** (30)–(31), accumulating $4$–$6$ pp of sovereign bonds without a conventional private-credit squeeze. The fiscal response is **boom-bust** (32).

**Corollary (The Central Paradox).** The model resolves two findings paradoxical under single-channel frameworks. First, non-default crises generate large output losses *without* default because the sovereign-ceiling channel transmits the spread regardless of repayment. Second, default episodes do *not* exhibit larger private-credit contractions despite the doom loop, because banks redirect capacity toward sovereign bonds rather than contracting. The structural difference: under non-default the sovereign-ceiling channel $\gamma s_t$ persists and deepens; under default it is severed and replaced by the more acute sudden-stop mechanism.

---

## 7. Competitive Equilibrium

A competitive equilibrium is a set of sequences $\{C_t, I_t, K_t, L_t, D_t, F_t, N_t, b^B_t, \ell_t, G_t, B_t, s_t, q_t, R^L_t, R^D_t, R^O_t, w_t\}$ such that: (i) households satisfy (1)–(3); (ii) firms satisfy (7)–(9); (iii) banks optimize subject to (10) and (16) with funding costs (11)–(12) and lending rate (17$'$); (iv) the government satisfies (21) and the default rule (24); (v) foreign investors satisfy (18), yielding (19); (vi) labor, deposit, and wholesale funding markets clear; and (vii) goods-market clearing holds:

$$Y_t = C_t + I_t + G_t + CA_t \tag{39}$$

In the non-crisis steady state, $\beta(R^* + \gamma s^{ss}) = 1$, output is $Y^{ss} = A(K^{ss})^{\alpha}\Psi(R^{L,ss})$ with $R^{L,ss} = R^* + \gamma s^{ss} + \varphi\, b^{B,ss}/N^{ss}$, and $\alpha A(K^{ss})^{\alpha-1} = R^{L,ss}-(1-\delta)$.

---

## 8. Log-Linearization and Analytical Solution

This section log-linearizes the equilibrium around the non-crisis steady state and derives closed-form impulse responses.

### 8.1 Conventions

For variable $X_t$, $\hat x_t \equiv \log(X_t/X^{ss})$. For the spread (which enters additively), $\hat s_t \equiv s_t - s^{ss}$. The spread follows an AR(1):

$$\hat s_t = \rho_s \hat s_{t-1} + \varepsilon_t, \qquad \rho_s \in (0,1) \tag{LL.0}$$

$\rho_s$ is identified from the autocorrelation of EMBIG spreads after onset; $\rho_s \to 1$ corresponds to permanently elevated spreads and the most persistent output losses.

### 8.2 The Log-Linearized System

**Lending rate** (from 17$'$), with $\Omega \equiv \varphi\, b^{B,ss}/N^{ss}$:

$$\Delta r^L_t = \gamma \hat s_t + \Omega(\hat n^{-1}_{t-1} - \acute n_t) \tag{LL.1}$$

In the non-default path $\hat n^{-1}_{t-1}=0$ (no new accumulation):

$$\Delta r^{L,nd}_t = \gamma \hat s_t - \Omega \acute n_t \tag{LL.1a}$$

**Output** (from 8), with $\varepsilon_p \equiv \xi n R^{L,ss}/[1+\xi(R^{L,ss}-1)]$:

$$\hat y_t = \alpha \hat k_t - \varepsilon_p \Delta r^L_t \;\;\Longrightarrow\;\; \hat y^{nd}_t = \alpha \hat k_t - \varepsilon_p(\gamma \hat s_t - \Omega \acute n_t) \tag{LL.2a}$$

**Capital** (from 5 and 9), with $\eta \equiv 1/[(1-\alpha)(R^{L,ss}-(1-\delta))]$:

$$\hat k_{t+1} = (1-\delta)\hat k_t - \eta\,\Delta r^L_t \;\;\Longrightarrow\;\; \hat k^{nd}_{t+1} = (1-\delta)\hat k_t - \eta(\gamma \hat s_t - \Omega \acute n_t) \tag{LL.3a}$$

**Bank net worth** (from 13 and 15), with $\mathrm B \equiv b^{B,ss}/[N^{ss}(1+s^{ss})^2]$ and net-worth multiplier $\Phi^N$:

$$\acute n^{nd}_{t+1} = \Phi^N \acute n_t - \mathrm B\,\hat s_{t+1} - \Phi^O \gamma \hat s_t \tag{LL.5a}$$

driven down by the bond-price loss ($\mathrm B \hat s_{t+1}$, balance-sheet channel) and higher wholesale funding costs ($\Phi^O \gamma \hat s_t$, sovereign-ceiling channel on the liability side). Stability requires $|\Phi^N|<1$.

**Private credit** (from 16):

$$\hat\ell^{nd}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\,\acute n_t \tag{LL.6}$$

**Consumption** (from 3$'$):

$$\hat c_t = \mathbb{E}_t[\hat c_{t+1}] - \frac{\gamma}{\sigma}\frac{\hat s_{t+1}}{R^{D,ss}} \tag{LL.7}$$

### 8.3 Closed-Form Impulse Responses (Non-Default)

Iterating LL.5a forward from $\acute n_0 = -\mathrm B \hat s_0$ with $\hat s_t = \rho_s^t \hat s_0$:

$$\boxed{\;\acute n^{nd}_h = (\Phi^N)^h \acute n_0 - \sum_{h'=1}^{h} (\Phi^N)^{h-h'}\!\left(\mathrm B\,\rho_s^{h'} + \Phi^O \gamma\,\rho_s^{h'-1}\right)\hat s_0\;} \tag{IR.1}$$

> **Correction.** An earlier draft factored this as $-(\mathrm B + \Phi^O\gamma)\sum(\Phi^N)^{h-h'}\rho_s^{h'}\hat s_0$, which is valid only at $\rho_s=1$. The $\mathrm B$ term sums $\rho_s^{h'}$ (the spread at $h'$) while the $\Phi^O\gamma$ term sums $\rho_s^{h'-1}$ (the spread one period earlier). The corrected form above keeps the two summations distinct. The qualitative results are unaffected; only the numerical IRF changes. In the code (`16_model_irf.do`) we iterate LL.5a directly, which is numerically exact and sidesteps the factoring entirely.

Iterating LL.3a from $\hat k_0 = 0$:

$$\hat k^{nd}_h = -\eta\gamma \sum_{h'=1}^{h}(1-\delta)^{h-h'}\rho_s^{h'-1}\hat s_0 + \eta\Omega \sum_{h'=1}^{h}(1-\delta)^{h-h'}\acute n_{h'-1} \tag{IR.2}$$

Substituting into LL.2a:

$$\hat y^{nd}_h = \alpha \hat k^{nd}_h - \varepsilon_p\gamma\rho_s^h \hat s_0 + \varepsilon_p\Omega\,\acute n^{nd}_h \tag{IR.3}$$

Three additive channels: the **capital channel** $\alpha\hat k_h$ (zero at $h=0$, growing thereafter); the **sovereign-ceiling working-capital channel** $-\varepsilon_p\gamma\rho_s^h \hat s_0$ (largest at $h=0$, decaying at $\rho_s$); and the **balance-sheet amplification** $\varepsilon_p\Omega\acute n_h$ (deepening as net worth erodes). The net response is negative at all horizons:

$$\hat y^{nd}_h < 0 \quad \forall h \ge 0 \quad\text{whenever}\quad \hat s_0 > 0 \tag{IR.3a}$$

The trough occurs where the deepening capital/balance-sheet effects outweigh the decaying spread shock — between $h=1$ and $h=3$ for plausible parameters, consistent with the empirical $-3.76$ pp trough at $h=1$. Private credit follows:

$$\hat\ell^{nd}_h = \frac{\lambda N^{ss}}{\ell^{ss}}\,\acute n^{nd}_h \tag{IR.4}$$

insignificant at $h=0$ (when $\acute n_0$ is small), deepening to $-3.48$ pp at $h=4$.

**Lemma (Persistence Without Mean Reversion).** Let $\rho_s,\Phi^N \in (0,1)$. Then (i) $\hat y^{nd}_h<0$ for all $h\ge 0$; (ii) the trough occurs at finite $h^*>0$; (iii) the cumulative output loss diverges as $\rho_s \to 1$, because the spread never returns to steady state and capital never recovers.

### 8.4 The Default Path

Three structural breaks modify the system:

**Break 1 — Wholesale funding collapse.** Setting $F_t=0$ replaces LL.7 with the autarky Euler equation using $R^{D,\text{aut}}_t > R^{D,ss}$ (§5.4), inducing a larger front-loaded consumption contraction.

**Break 2 — Gambling for resurrection.** When (30) holds, $\ell^{\text{def}}_t = 0$ and $b^{B,\text{def}}_t = \lambda N_t$, so $\hat\ell^{\text{def}}_t \ll \hat\ell^{nd}_t$ for all $h$.

> **Proposition 3 (Credit Paradox).** Private lending under default is strictly lower than under non-default ($\hat\ell^{\text{def}}_h < \hat\ell^{nd}_h < 0$), yet the *gross* bank balance sheet does not shrink — banks redirect the full capacity $\lambda N_t$ toward sovereign bonds. The econometrician observes bank credit to government rising while credit to the private sector falls, with the aggregate roughly unchanged. This is why default-linked episodes show **no statistically larger private-credit contraction** despite the doom loop.

**Break 3 — Doom-loop accumulation.** $\Delta b^{B,\text{def}}_0 = \lambda N_0 - b^B_{-1}$; log-linearizing:

$$\hat n^{B,\text{def}}_t = \left(\frac{\lambda N^{ss}}{b^{B,ss}} - 1\right) - \frac{\lambda N^{ss}}{b^{B,ss}}\,\mathrm B\,\hat s_0 \tag{LL.10}$$

positive when the reallocation gain exceeds net-worth erosion — exactly the gambling-for-resurrection condition.

**Output under default.** With $\ell^{\text{def}}_t = 0$ the working-capital constraint collapses; investment is suspended so capital depreciates at $\delta$:

$$\hat y^{\text{def}}_h = h\,\ln(1-\delta)\cdot\alpha - \varepsilon_p\,\Delta r^{L,\text{aut}} \tag{IR.5}$$

> **Correction.** An earlier draft wrote the capital term as $\alpha(1-\delta)^h \hat k_0$ with $\hat k_0=0$, which is identically zero and contradicts the stated deepening. With $I_t=0$ the capital stock evolves $K_{t+h}=(1-\delta)^h K_0$, so in log deviation $\hat k^{\text{def}}_h = h\ln(1-\delta) < 0$ — depreciation accumulates and output deepens. The impact loss is the autarky working-capital wedge $-\varepsilon_p \Delta r^{L,\text{aut}}$; the deepening reflects capital depletion.

---

## 9. Mapping to the Empirical Findings

### Table 1 — Model Sign Predictions by Variable, Path, and Horizon

| Variable | Path | $h=0$ | $h=1\text{–}2$ | $h=3\text{–}4$ | Mechanism |
|---|---|:---:|:---:|:---:|---|
| **Output** | Non-default | $-$ | $-$ | $-/\approx 0$ | working capital + balance sheet |
| | Default | $--$ | $--$ | $--$ | autarky working capital + capital depletion |
| | *Comparison* | nd $>$ d | nd $>$ d | nd $>$ d | default always deeper |
| **Private credit** | Non-default | $-$ (small) | $-$ (growing) | $-$ (deep) | net-worth channel (IR.4) |
| | Default | $\approx 0$ | $\approx 0$ | $\approx 0$ | autarky: banks redirect to $b^B$ |
| | *Comparison* | nd $<$ d | nd $<$ d | nd $<$ d | nd has worse *measured* credit |
| **Investment** | Non-default | $-$ | $-$ | $-/\approx 0$ | user cost rises (eq. 9) |
| | Default | $\approx 0$ | $-$ | $--$ | capital depletes ($I=0$) |
| | *Comparison* | nd $>$ d | nd $>$ d | nd $>$ d | |
| **Claims on govt** | Non-default | $\approx 0$ | $\approx 0$ | $\approx 0$ | no doom-loop incentive |
| | Default | $+$ | $++$ | $+$ | gambling for resurrection (31) |
| | *Comparison* | nd $<$ d | nd $<$ d | nd $<$ d | default: $+4$–$6$ pp accumulation |
| **Current account** | Non-default | $+$ (gradual) | $+$ (fading) | $\approx 0$ | wholesale cost rises, flows continue |
| | Default | $++$ | $++$ | $+$ | two-component sudden stop (28) |
| | *Comparison* | nd $<$ d | nd $<$ d | nd $\approx$ d | default front-loaded |
| **Govt expenditure** | Non-default | $+$ | $+$ | $+$ | retained market access |
| | Default | $++$ | $-$ | $--$ | boom-bust: windfall then conditionality |
| | *Comparison* | nd $<$ d | nd $>$ d | nd $>$ d | paths cross near $h\approx 1$ |

### Table 2 — Empirical Findings, Model Counterparts, and Parameter Conditions

| Empirical finding | Model equation(s) | Parameter restriction |
|---|---|---|
| **(a)** ND output loss ($-1.86$ pp $h{=}0$, $-3.76$ pp $h{=}1$) | LL.2a, IR.3 | $\varepsilon_p>0$ (req. $\xi>0$), $\gamma>0$, $\hat s_0>0$ |
| **(b)** ND trough at $h{=}1$, not $h{=}0$ | IR.3 + IR.2 | $\Phi^N\in(0,1)$, $\Omega>0$: net-worth deepening overtakes spread decay at $h{=}1$ |
| **(c)** Default deeper & front-loaded ($-3.42$ pp $h{=}0$, $-5.50$ pp $h{=}1$) | IR.5, eq. 34 | $\Delta r^{L,\text{aut}}>\gamma\hat s_0$, guaranteed by $K^{\text{aut}}<K^{ss}$ |
| **(d)** Credit contraction concentrated in ND ($-4.78$ pp $h{=}4$; $\approx 0$ default) | IR.4, LL.6D (Prop. 3) | ND: $\acute n<0$; default: GFR condition (30) holds $\Rightarrow \ell=0$ |
| **(e)** Doom loop $+4$–$6$ pp in default only | eqs. 30–31, LL.10 | $\lambda N_0 > b^B_{-1}$; $\mu\theta R^{BDK}+\gamma^B N \ge R^L(1-p^{RDL})$ |
| **(f)** Diabolic loop absorbed by lagged $b^B$ | eqs. 14–15, 17$'$ | amplification $\propto b^B_{t-1}$: controlling $L.b^B$ absorbs $\mathrm B$ in LL.5a |
| **(g)** CA adjustment larger & earlier under default ($+2.61$ pp $h{=}1$) | eq. 28, IR (LL.9) | $F_{-1}>0$, $\delta^B>0$: sudden stop $\gg$ gradual adjustment |
| **(h)** Govexp boom-bust in default ($+2.19$ pp $h{=}0$, $-1.69$ pp $h{=}4$) | eq. 32 | $\delta^B>0$ (windfall); exclusion $>2$ periods (conditionality bites) |
| **(i)** Govexp sustained in ND (positive throughout) | eq. 21 + ceiling (11) | $\gamma<\bar\gamma$: spread rise insufficient to trigger solvency crisis |

**Theorem (Joint Consistency).** Under the conditions in Table 2, the analytical impulse responses of the log-linearized model — non-default IR.1–IR.4 and default IR.5, LL.6D, LL.10 — jointly reproduce the qualitative and directional predictions of all nine empirical findings. No individual result requires a parameter restriction that contradicts another; there exists an open set $\Theta \subset \mathbb{R}^{11}_+$ of parameter values under which all findings hold simultaneously.

---

## 10. Calibration and Numerical Solution

The model is taken to the data in three parameter groups and solved as two blocks. Implementation is in `do/14_calibration.do`, `do/15_solve_default.do`, and `do/16_model_irf.do`.

### 10.1 Parameter Groups

**Group 1 — Externally calibrated (literature).** Standard annual EM values: $\beta=0.96$, $\sigma=2$, $\alpha=0.33$, $\delta=0.10$, $R^*=0.04$; recovery $\theta=0.62$ (Cruces-Trebesch 2013, mean haircut $\approx 38\%$); re-entry $\mu=0.22$ (Gelos et al. 2011, exclusion $\approx 4.5$ years); maturing fraction $\delta^B=0.22$; leverage $\lambda=10$.

**Group 2 — Data-determined steady-state ratios.** Computed from *tranquil-period* means (non-crisis, non-continuation country-years) in the panel: the background spread $s^{ss}$; expenditure shares $s_C, s_I, s_G, s_{CA}$; tax rate $\tau$; bank sovereign exposure $b^{B,ss}/Y$; net worth $N^{ss}/Y$; and hence the amplification $\Omega = \varphi\, b^{B,ss}/N^{ss}$ and the structural elasticities $\varepsilon_p$, $\eta$, $\mathrm B$.

**Group 3 — Estimated processes.** Spread persistence $\rho_s$ from a panel AR(1) on EMBIG spreads; income process $(\rho_y, \sigma_y)$ from a detrended AR(1) on log GDP per capita.

The deep transmission parameters $\{\xi, \varphi, \Phi^N\}$ are calibrated by **Simulated Method of Moments to the empirical non-default IRF only**.

### 10.2 Two-Block Solution

**Block 1 — Nonlinear default block** (`15_solve_default.do`). A canonical Arellano (2008) value-function iteration over states $(B, y)$: Tauchen-discretized income, the default decision (24), and the bond-price fixed point (18)–(19) with recovery $\theta$ and re-entry $\mu$. Simulation of 100,000 periods yields model moments (default frequency, mean/SD spread, debt/GDP) and a **crisis-event spread path** ($h=-2$ to $4$) that serves as the forcing process for Block 2.

**Block 2 — Log-linear transmission block** (`16_model_irf.do`). The model spread path feeds the banking equations LL.5a, LL.3a, LL.2a, and IR.4 (iterated directly), producing model IRFs for output, net worth, capital, and credit.

### 10.3 Identification and Validation

The calibration uses **only the non-default IRF** to pin $\{\xi, \varphi, \Phi^N\}$. The **default path is then validated out of sample**: a single autarky lending-rate wedge $\Delta r^{L,\text{aut}}$ (the one default-specific object, §5.4) is pinned to the $h=0$ default data point, after which $h=1$ to $h=4$ of the default IRF are *zero-free-parameter predictions* generated by IR.5. The degree to which these match the empirical default LP constitutes the model's central test, converting the framework from a fit into a falsifiable prediction.

### 10.4 What the Analytical Solution Cannot Pin Down

Two findings require a full numerical calibration rather than the analytical solution: the *specific magnitudes* of the troughs ($-2.33$ pp non-default vs. $-5.50$ pp default), which need value-function iteration to match levels; and the *insignificance of the frontier interaction*, which reflects low statistical power ($\approx 8$ frontier onset episodes) rather than a model feature.

---

## 11. Scope and Limitations

The sovereign-ceiling pass-through $\gamma$ and the wholesale-funding share $\omega$ are treated as constants, whereas in practice both vary with financial integration and the currency composition of bank liabilities. The default decision abstracts from renegotiation dynamics, long-maturity debt, and partial default. TFP $A$ is constant, abstracting from the endogenous productivity feedback in Mendoza and Yue (2012). These simplifications serve the model's purpose — analytical clarity on the transmission channels identified empirically — rather than quantitative discipline. A full quantitative extension would integrate the dual funding structure with a Bocola (2016)-type balance-sheet block for the non-default channel and a Gennaioli-Martin-Rossi (2014) framework for the default-linked channel.

---

## Notation Summary

| Symbol | Meaning | Symbol | Meaning |
|---|---|---|---|
| $C_t, D_t$ | consumption, deposits | $s_t$ | EMBIG sovereign spread |
| $K_t, I_t$ | capital, investment | $q_t$ | sovereign bond price $=1/(1+s_t)$ |
| $\ell_t, b^B_t$ | bank lending, sovereign bond holdings | $\pi_{t+1}$ | default probability |
| $N_t, F_t$ | bank net worth, wholesale funding | $\theta, \mu$ | recovery rate, re-entry prob. |
| $R^D_t, R^O_t, R^L_t$ | deposit, wholesale, lending rates | $\gamma, \varphi$ | ceiling pass-through, balance-sheet pass-through |
| $G_t, B_t$ | govt expenditure, external debt | $\xi, \lambda$ | working-capital share, bank leverage |
| $\Omega, \mathrm B$ | amplification, balance-sheet sensitivity | $\Phi^N, \Phi^O$ | net-worth multiplier, wholesale income share |
| $\varepsilon_p, \eta$ | output & capital elasticities to $R^L$ | $\rho_s, \rho_y$ | spread, income persistence |

---

*Empirical implementation: `do/14_calibration.do` (calibration), `do/15_solve_default.do` (nonlinear default block, Arellano VFI), `do/16_model_irf.do` (log-linear transmission, SMM, model-vs-data overlay). See `README.md` for the full pipeline.*
