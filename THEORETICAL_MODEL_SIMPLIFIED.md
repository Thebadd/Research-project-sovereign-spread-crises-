# Theoretical Model — Simplified Representation

## A Three-Block Model of Sovereign Spread Crises

This appendix presents a parsimonious model that captures the two transmission mechanisms
identified empirically: a **credit-channel recession** under non-default spread crises, and a
**capital-depletion recession** under default-linked crises. The model consists of three blocks —
firms, banks, and the sovereign-ceiling mechanism — plus a mechanical default rule. Despite its
simplicity, the model jointly reproduces all nine empirical findings reported in the paper.

---

## 1. Firms and Production

A representative firm produces output $Y_t$ using capital $K_t$ and labour $L_t$ according to a
Cobb-Douglas production function:

$$Y_t = A K_t^\alpha L_t^{1-\alpha} \qquad\qquad (1)$$

where $A > 0$ is total factor productivity, $\alpha \in (0,1)$ is the capital share (calibrated to
$\alpha = 0.33$), and $L_t$ is labour input. The firm must borrow from domestic banks to
**pre-finance its wage bill** before production takes place — the working-capital constraint. If
the firm borrows a fraction $\xi \in (0,1)$ of its total wage bill $w_t L_t$ at the domestic
lending rate $R^L_t$, its effective unit cost of labour is $(1 + \xi(R^L_t - 1))$ rather than
$w_t$. Optimisation over labour yields an output supply equation in which output is strictly
decreasing in the lending rate:

$$Y_t = A K_t^\alpha \cdot \Psi(R^L_t), \qquad \Psi(R^L_t) \equiv \left[\frac{1-\alpha}{1 + \xi(R^L_t - 1)}\right]^n, \quad n = \frac{1-\alpha}{\alpha} \qquad\qquad (2)$$

where $\Psi'(R^L_t) < 0$: any increase in the lending rate reduces output by raising the
effective cost of the working-capital loan. Log-linearising around the steady state gives the
**output deviation** $\hat{y}_t \equiv \log(Y_t/Y^{ss})$:

$$\hat{y}_t = \alpha \hat{k}_t - \varepsilon_p (R^L_t - R^{L,ss}) \qquad\qquad (3)$$

where $\hat{k}_t \equiv \log(K_t/K^{ss})$ is the capital stock deviation from steady state,
$\varepsilon_p \equiv \xi n R^{L,ss}/[1 + \xi(R^{L,ss}-1)] > 0$ is the **working-capital output
elasticity** — the percentage-point fall in output per unit rise in the lending rate — and
$R^{L,ss}$ is the steady-state lending rate. Two channels drive $\hat{y}_t$: the **capital
channel** $\alpha\hat{k}_t$ (zero at impact since capital is predetermined, deepening over
time) and the **working-capital channel** $-\varepsilon_p(R^L_t - R^{L,ss})$ (operating
immediately when the lending rate rises).

**Capital accumulation.** Investment $I_t$ adds to the capital stock net of depreciation at
rate $\delta \in (0,1)$:

$$K_{t+1} = (1-\delta)K_t + I_t \qquad\qquad (4)$$

where $\delta = 0.10$ (annual). The investment decision is governed by a standard user-cost
condition: firms invest until the net return to capital equals the lending rate. A rise in
$R^L_t$ therefore reduces investment, and capital accumulates more slowly, deepening the
output loss progressively over horizons $h = 1, 2, 3, 4$.

---

## 2. Banks, Net Worth, and the Leverage Constraint

A representative bank holds two assets: domestic sovereign bonds $b^B_t$ and private loans
$\ell_t$ to firms. Its **net worth** $N_t$ (equity capital) satisfies the balance-sheet identity:

$$b^B_t + \ell_t = \text{Total assets} = \text{Debt} + N_t$$

A regulatory or market-imposed **leverage constraint** limits total assets to a multiple
$\lambda > 1$ of net worth:

$$b^B_t + \ell_t \leq \lambda N_t \qquad\qquad (5)$$

where $\lambda = 10$ (bank assets are ten times equity, standard for emerging-market banks).
This constraint is **always binding** in equilibrium: banks maximise returns by fully deploying
their balance-sheet capacity. Equation (5) therefore holds with equality.

**Net worth dynamics.** Net worth evolves through retained earnings minus losses on the
sovereign bond portfolio. When the sovereign spread $s_t$ rises, the mark-to-market value of
sovereign bonds $b^B_t$ falls, generating a **capital loss** proportional to the bank's
pre-existing sovereign exposure. The log-linearised net-worth equation is:

$$\hat{n}_{t+1} = \Phi^N \hat{n}_t - \mathrm{B}\, \hat{s}_{t+1} \qquad\qquad (6)$$

where $\hat{n}_t \equiv \log(N_t/N^{ss})$ is the net-worth deviation, $\Phi^N \in (0,1)$ is the
**net-worth persistence parameter** governing the speed at which banks rebuild equity through
retained earnings (calibrated by SMM), and $\mathrm{B} \equiv b^{B,ss}/[N^{ss}(1+s^{ss})^2] > 0$
is the **balance-sheet sensitivity** — the fraction of net worth lost per unit rise in the spread,
proportional to the steady-state ratio of sovereign bond holdings to bank equity. A bank with
large pre-existing sovereign exposure (high $b^{B,ss}/N^{ss}$) suffers a larger capital loss
for the same spread shock.

**Private credit.** Since the leverage constraint (5) binds with equality and sovereign bond
holdings $b^B_t$ adjust sluggishly (banks cannot instantaneously sell illiquid bonds), private
lending contracts in proportion to the net-worth loss, amplified by leverage:

$$\hat{\ell}_t = \frac{\lambda N^{ss}}{\ell^{ss}}\, \hat{n}_t \qquad\qquad (7)$$

The amplification factor $\lambda N^{ss}/\ell^{ss} > 1$ means that a 1% fall in net worth
generates a more than 1% fall in private credit. This leverage amplification is the core
balance-sheet mechanism: small sovereign shocks produce disproportionately large credit
contractions.

---

## 3. The Sovereign Ceiling and the Domestic Lending Rate

Domestic banks fund themselves on international wholesale markets at the **sovereign ceiling
rate**. The sovereign ceiling is the empirical regularity that no private borrower in a country
can obtain external funding at a rate below its sovereign. Formally, the international
wholesale funding cost $R^O_t$ faced by banks is:

$$R^O_t = R^* + \gamma s_t \qquad\qquad (8)$$

where $R^* > 0$ is the international risk-free rate (calibrated to 4% annually), $s_t \geq 0$
is the sovereign EMBIG spread (in decimal), and $\gamma \in (0,1]$ is the **sovereign-ceiling
pass-through parameter** — the fraction of the sovereign spread transmitted to banks' wholesale
funding cost. When $\gamma = 1$, the full spread is passed through; when $\gamma < 1$, domestic
financial frictions or partial currency mismatch dampen the transmission.

Banks price private loans at a mark-up over their funding cost, adding a **balance-sheet
premium** that rises when their net worth is depleted relative to their sovereign bond
exposure. The domestic lending rate is:

$$R^L_t = R^* + \gamma s_t + \varphi \frac{b^B_t}{N_t} \qquad\qquad (9)$$

where $\varphi > 0$ is the **balance-sheet pass-through parameter** — the sensitivity of the
lending rate to the bank's leverage in sovereign bonds. The ratio $b^B_t/N_t$ rises both when
the bank accumulates sovereign bonds (increasing $b^B_t$) and when net worth erodes (reducing
$N_t$), so the lending-rate premium is self-amplifying during a sovereign stress episode.

Log-linearising (9) around the steady state yields the **lending-rate deviation**:

$$R^L_t - R^{L,ss} = \gamma \hat{s}_t - \Omega \hat{n}_t \qquad\qquad (10)$$

where $\Omega \equiv \varphi\, b^{B,ss}/N^{ss} > 0$ is the **balance-sheet amplification
coefficient** — the sensitivity of the lending rate to net-worth erosion, proportional to the
steady-state sovereign exposure ratio. The two terms in (10) capture the two routes by which a
sovereign spread shock raises lending rates: the **direct route** $\gamma \hat{s}_t$ (the
ceiling pass-through operating on all banks regardless of their sovereign exposure) and the
**indirect route** $-\Omega \hat{n}_t$ (the balance-sheet amplification operating through the
erosion of net worth, which forces banks to charge a higher premium to maintain profitability).

Substituting (10) into the output equation (3) gives the **combined output equation**:

$$\hat{y}_t = \alpha \hat{k}_t - \varepsilon_p \left(\gamma \hat{s}_t - \Omega \hat{n}_t\right) \qquad\qquad (11)$$

This single equation captures the three-channel transmission of a sovereign spread shock to
output: (i) through accumulated capital stock depletion $\alpha\hat{k}_t$; (ii) through the
direct sovereign-ceiling working-capital cost $-\varepsilon_p \gamma \hat{s}_t$; and (iii)
through the balance-sheet amplification $\varepsilon_p \Omega \hat{n}_t$ (positive since
$\hat{n}_t < 0$ during a crisis, making this term negative and adding to the output loss).

---

## 4. Mechanical Default Rule and the Two Paths

Rather than solving a full sovereign optimisation problem, we classify episodes using the
empirical classification from the crisis database (see Section 2 of the paper). The model
then applies a **mechanical default rule**: when an episode is classified as default-linked,
three structural breaks are imposed simultaneously at the onset date $t = 0$.

**Break 1 — Investment suspension.** Market exclusion cuts the country's access to external
financing. With $F_t = 0$ (wholesale funding gone) and domestic savings insufficient to
replace it, investment falls to zero:

$$I_t = 0 \quad \text{for all } t \text{ during exclusion} \qquad\qquad (12)$$

Capital therefore declines at the depreciation rate with no offsetting investment:

$$K^{\text{def}}_{t+h} = (1-\delta)^h K_0 \implies \hat{k}^{\text{def}}_h = h \ln(1-\delta) < 0 \qquad\qquad (13)$$

This is a **cumulative and worsening** mechanism: at $h = 1$, capital has fallen by 10%;
at $h = 4$, it has fallen by 34%. Since output is proportional to $K^\alpha$, the capital
channel alone drives output down by $\alpha \times 34\% \approx 11$ pp by $h = 4$, even
without any change in the lending rate.

**Break 2 — Private lending collapse.** With no external funding and net worth eroded by
the spread shock, banks redirect their entire balance-sheet capacity toward sovereign bonds
(**gambling for resurrection**: by holding government bonds rather than private loans, the bank
retains a claim on the recovery value $\theta > 0$ of restructured debt):

$$\ell^{\text{def}}_t = 0, \qquad b^{B,\text{def}}_t = \lambda N_t \qquad\qquad (14)$$

Private credit collapses to zero. Firms cannot pre-finance the wage bill, so the
working-capital output loss is **maximal**: $\Psi(R^{L,\text{aut}})$ evaluated at the
autarky lending rate $R^{L,\text{aut}} \gg R^{L,ss}$.

**Break 3 — Sovereign ceiling severed.** With $F_t = 0$, banks no longer borrow on
international markets, so the sovereign-ceiling channel $\gamma s_t$ in equation (8) is
**severed**. The domestic lending rate is no longer pinned to the sovereign spread; instead it
is set by the autarky equilibrium, which equates the net return to capital to the lending rate.
Since capital is now scarce ($K^{\text{def}} < K^{ss}$), the marginal product of capital
rises, and the autarky rate $R^{L,\text{aut}} > R^{L,ss}$ even though the sovereign spread
channel is inactive.

**Output under default.** Combining (13) and Break 2–3, the output equation (11) under
default simplifies to:

$$\hat{y}^{\text{def}}_h = \alpha \cdot h \ln(1-\delta) - \varepsilon_p\, \Delta r^{L,\text{aut}} \qquad\qquad (15)$$

where $\Delta r^{L,\text{aut}} \equiv R^{L,\text{aut}} - R^{L,ss} > 0$ is the autarky
lending-rate wedge. The first term is **zero at $h=0$** (capital is predetermined) and
**strictly more negative at each subsequent horizon** — generating the empirically observed
pattern of deepening output losses under default. The second term is the **impact loss at
$h=0$**, driven entirely by the working-capital channel under autarky.

---

## 5. Impulse Responses

### 5.1 Non-Default Path

In the non-default path, all three structural breaks are absent: investment continues
($I_t > 0$), private lending is positive ($\ell_t > 0$), and the sovereign ceiling channel
is active ($R^O_t = R^* + \gamma s_t$). The spread follows an AR(1) process:

$$\hat{s}_t = \rho_s^t \hat{s}_0, \qquad \rho_s \in (0,1) \qquad\qquad (16)$$

where $\rho_s$ is the spread persistence parameter (estimated from the panel AR(1) in the
data, see Section 3 of the paper) and $\hat{s}_0 > 0$ is the onset-period spread deviation.

**Net worth at horizon $h$** (iterating equation 6):

$$\hat{n}^{nd}_h = (\Phi^N)^h \hat{n}_0 - \sum_{h'=1}^{h}(\Phi^N)^{h-h'}\,\mathrm{B}\,\rho_s^{h'}\hat{s}_0 \qquad\qquad (17)$$

where $\hat{n}_0 = -\mathrm{B}\hat{s}_0$ is the **impact balance-sheet loss** at $h=0$. Net
worth deteriorates progressively as the still-elevated spread generates fresh capital losses
at each horizon. Recovery begins only once the spread has decayed sufficiently for retained
earnings to outweigh new losses ($\Phi^N |\hat{n}_{h-1}| > \mathrm{B}\rho_s^h \hat{s}_0$).

**Private credit at horizon $h$** (from equation 7):

$$\hat{\ell}^{nd}_h = \frac{\lambda N^{ss}}{\ell^{ss}}\,\hat{n}^{nd}_h \qquad\qquad (18)$$

Credit is negligible at $h=0$ (since $\hat{n}_0$ is small when $b^{B,ss}/N^{ss}$ is modest),
deepening progressively as net worth erodes. The empirical estimate is $-4.78$ pp at $h=4$.

**Output at horizon $h$** (from equations 3–4 and 10, with $\hat{k}_0 = 0$):

$$\hat{y}^{nd}_h = \alpha\hat{k}^{nd}_h - \varepsilon_p(\gamma\rho_s^h\hat{s}_0 - \Omega\hat{n}^{nd}_h) \qquad\qquad (19)$$

where $\hat{k}^{nd}_h = -\eta\gamma\sum_{h'=1}^{h}(1-\delta)^{h-h'}\rho_s^{h'-1}\hat{s}_0 + \eta\Omega\sum_{h'=1}^{h}(1-\delta)^{h-h'}\hat{n}_{h'-1}$ accumulates
from the user-cost investment channel, and $\eta \equiv 1/[(1-\alpha)(R^{L,ss}-(1-\delta))] > 0$
is the interest semi-elasticity of capital. All three components of (19) are negative, and the
output trough occurs between $h=1$ and $h=3$ for plausible parameters, consistent with the
empirical trough of $-3.76$ pp at $h=1$.

### 5.2 Default Path

Under the mechanical default rule (equations 12–14), the impulse responses simplify
dramatically:

$$\hat{y}^{\text{def}}_h = \underbrace{\alpha \cdot h \ln(1-\delta)}_{\text{capital depletion}} - \underbrace{\varepsilon_p\, \Delta r^{L,\text{aut}}}_{\text{autarky working capital}} \qquad\qquad (20)$$

The capital-depletion term is **strictly increasing in magnitude** at each horizon, generating
the observed deepening of the recession over $h=0,\ldots,4$. The autarky working-capital term
is **constant across horizons** — it generates the impact loss at $h=0$ and persists as long
as the country remains excluded.

---

## 6. The Credit/GDP Paradox

A central empirical finding of the paper is that **private credit/GDP contracts more in
non-default episodes than in default-linked episodes**, despite private lending collapsing to
zero under default. This apparent paradox is resolved by the **GDP-denominator effect**.

The empirical outcome variable is the ratio $\ell_t / Y_t$, not private credit in levels. Under
default, the output collapse is driven by capital depletion (equation 20), which is cumulative
and worsening: by $h=3$, GDP has fallen by $\alpha \times |3\ln(1-\delta)| \approx 9$ pp from
capital depletion alone, in addition to the working-capital loss. The GDP denominator therefore
shrinks **at least as fast as private credit**, leaving the ratio approximately unchanged.

Under non-default, by contrast, the credit channel is the dominant transmission mechanism:
$\ell_t$ (the numerator) falls much faster than $Y_t$ (the denominator), so the ratio
$\ell_t / Y_t$ falls substantially.

**Proposition (Credit/GDP Paradox).** Let $\Delta\ell / Y$ denote the change in the
credit-to-GDP ratio. Then under the conditions of this model:

$$\left|\frac{\Delta \ell}{Y}\right|^{nd} > \left|\frac{\Delta \ell}{Y}\right|^{def}$$

even though $\ell^{def}_t = 0 < \ell^{nd}_t$ for all $h \geq 0$. The resolution is:

$$\text{Non-default:} \quad |\Delta\ell_t| \gg |\Delta Y_t| \quad \Rightarrow \quad \text{ratio falls (credit channel dominates)}$$
$$\text{Default:} \quad |\Delta Y_t| \gg |\Delta\ell_t| \quad \Rightarrow \quad \text{ratio roughly stable (capital-depletion channel dominates)}$$

**A smaller credit/GDP response in default episodes is not evidence that banks are lending more
— it is evidence that GDP is collapsing faster than credit, driven by a completely different
transmission mechanism.**

---

## 7. Summary of Mechanisms by Path

| Variable | Non-default path | Default path | Dominant mechanism |
|---|---|---|---|
| **Output** | $-3.76$ pp at $h=1$ | $-5.50$ pp at $h=1$, deepening | ND: working capital + balance sheet; DEF: capital depletion |
| **Private credit/GDP** | $-4.78$ pp at $h=4$ | $\approx 0$ measured | ND: $\ell_t$ falls faster than $Y_t$; DEF: $Y_t$ falls faster than $\ell_t$ |
| **Investment** | Reduced ($I_t > 0$) | Suspended ($I_t = 0$) | DEF: market exclusion |
| **Sovereign bonds held by banks** | Unchanged ($\hat{n}^B \approx 0$) | $+4$–$6$ pp (doom loop) | DEF: gambling for resurrection |
| **Current account** | Gradual improvement | Front-loaded surge | DEF: sudden stop ($F_t = 0$ + maturing debt) |

---

## 8. Calibration

The model has **three free parameters** calibrated by Simulated Method of Moments (SMM)
targeting the non-default empirical IRFs for output and private credit:

$$\min_{\gamma,\,\varphi,\,\Phi^N} \sum_{h=0}^{4}\left[(\hat{y}^{nd,\text{model}}_h - \hat{y}^{nd,\text{LP}}_h)^2 + (\hat{\ell}^{nd,\text{model}}_h - \hat{\ell}^{nd,\text{LP}}_h)^2\right] \qquad\qquad (21)$$

All remaining parameters are fixed externally:

| Parameter | Value | Source |
|---|---|---|
| $\alpha$ | 0.33 | Standard macro |
| $\delta$ | 0.10 | Annual depreciation |
| $R^*$ | 0.04 | US Treasury rate |
| $\lambda$ | 10 | EM bank leverage |
| $\xi$ | 0.50 | Working-capital share |
| $\theta$ | 0.62 | Cruces-Trebesch (2013) |

The autarky lending-rate wedge $\Delta r^{L,\text{aut}}$ is pinned to the single $h=0$ default
moment, and horizons $h=1,\ldots,4$ under default are **genuine out-of-sample predictions**
with zero free parameters — providing a clean test of whether the capital-depletion mechanism
in equation (15) is consistent with the data.

---

## 9. Notation Summary

| Symbol | Definition |
|---|---|
| $Y_t$, $K_t$, $L_t$ | Output, capital, labour |
| $A$, $\alpha$, $\delta$ | TFP, capital share, depreciation rate |
| $\xi$ | Working-capital share of wage bill |
| $\varepsilon_p$ | Working-capital output elasticity |
| $\eta$ | Interest semi-elasticity of capital |
| $b^B_t$ | Bank holdings of sovereign bonds |
| $\ell_t$ | Bank private loans to firms |
| $N_t$ | Bank net worth (equity) |
| $\lambda$ | Bank leverage multiplier |
| $\mathrm{B}$ | Balance-sheet sensitivity of net worth to spread |
| $\Phi^N$ | Net-worth persistence parameter |
| $R^*$ | International risk-free rate |
| $s_t$ | Sovereign EMBIG spread (decimal) |
| $\gamma$ | Sovereign-ceiling pass-through |
| $R^O_t$ | International wholesale funding rate (sovereign ceiling) |
| $R^L_t$ | Domestic lending rate |
| $\varphi$ | Balance-sheet pass-through to lending rate |
| $\Omega$ | Balance-sheet amplification coefficient ($= \varphi\, b^{B,ss}/N^{ss}$) |
| $\rho_s$ | Spread persistence (AR(1) coefficient) |
| $\hat{x}_t$ | Log-deviation of $x_t$ from steady state ($\log(x_t/x^{ss})$) |
| $h$ | Horizon in years after crisis onset |
| $\Delta r^{L,\text{aut}}$ | Autarky lending-rate wedge ($R^{L,\text{aut}} - R^{L,ss}$) |
| $\theta$ | Recovery rate on restructured sovereign debt |
