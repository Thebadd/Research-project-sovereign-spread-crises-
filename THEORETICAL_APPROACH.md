# Theoretical Approach: Sovereign Spread Crises and Real Transmission

## 1. The Question the Model Must Answer

Sovereign spread crises come in two varieties: episodes where the country
eventually defaults on its debt, and episodes where it does not. Both involve
a large and persistent rise in the sovereign spread. The question is whether
the real economy — output, credit, investment — responds differently in each
case, and if so, through which channels.

Crucially, the model is built **before** looking at the impulse response
estimates. The structure comes from theory (what mechanisms exist in an
economy with sovereign risk and financial intermediaries); the parameters are
calibrated to moments that are entirely independent of the IRFs being tested
(default frequency, mean and standard deviation of spreads, aggregate
debt/GDP). The comparison of model predictions to the LP estimates then
constitutes a genuine out-of-sample test — the model could fail.

---

## 2. Model Structure: Two Blocks

### Block 1 — Sovereign Default Risk (Nonlinear)

The sovereign faces a standard income-fluctuation problem with default option
(Arellano 2008). Each period the government chooses whether to repay or
default, trading off short-run consumption relief against long-run exclusion
from capital markets. Risk-neutral foreign investors price bonds competitively:
the bond price q(B', y) reflects the probability-weighted payoff, where
default delivers recovery θ and repayment delivers 1.

**What this block generates endogenously:**
- A bond price schedule q(B', y) that embeds default risk
- A sovereign spread: s_t = 1/q(B', y) − 1 − R*
- Default events d_t ∈ {0,1} as an endogenous binary outcome
- After default: exclusion from markets with re-entry probability μ per period

**Critically, the model is agnostic about which crises become defaults.** The
default decision is endogenous: it depends on the realized income draw and the
debt level. From a long simulation, we separately collect two types of crisis
episodes:

- **Non-default crisis:** sovereign spread crosses 1000 bps AND no default
  occurs in the following four years.
- **Pre-default episode:** sovereign spread crosses 1000 bps AND default
  occurs within four years.

These are two ex-post realizations of the same stochastic economy. The model
does not assume which path is worse — it computes the spread path for each
type of episode from the simulation and passes both paths to Block 2.

### Block 2 — Banking Transmission (Log-Linear)

The banking sector intermediates between the sovereign and the real economy.
Banks hold two assets: productive loans to firms (ℓ_t) and sovereign bonds
(b^B_t). A binding leverage constraint limits lending to a multiple of bank
net worth N_t:

```
ℓ_t = λ · N_t
```

The sovereign-ceiling assumption links the bank lending rate to the sovereign
spread: R^L_t = R* + γ · s_t. Banks holding sovereign bonds suffer a
balance-sheet loss when spreads rise:

```
ΔN_t ∝ b^B_{t-1} · Δq_t
```

where q_t is the bond price. From these two mechanisms, the model derives
log-linearized recursions for net worth, capital, output, and credit.

**Net worth (LL.5a):**
```
n_h = Φ_N · n_{h-1} − B_sens · s_h − Φ^O · γ · s_{h-1}
```

**Capital accumulation (LL.3a):**
```
k_h = (1−δ) · k_{h-1} − η · (γ · s_{h-1} − Ω · n_{h-1})
```
where η = δ / [(1−α)(R^L_ss − (1−δ))] captures the investment-rate
sensitivity (the δ factor follows from the capital accumulation identity
k_{t+1} = (1−δ)k_t + δi_t).

**Output (LL.2a):**
```
y_h = α · k_h − ε_p · (γ · s_h − Ω · n_h)
```

**Private credit (IR.4):**
```
ℓ_h = (λ · N_ss / ℓ_ss) · n_h
```

These equations jointly determine output and credit responses given any spread
path s_h as input.

---

## 3. What the Model Predicts A Priori (Before Looking at Data)

The two-block structure generates predictions for both types of episode without
any tuning of the transmission parameters to the outcome IRFs.

### Prediction for non-default crises

The model is run with the simulated spread path for non-default episodes (high
spreads, no actual default). Only the balance-sheet and working-capital
channels operate. The capital channel is mild because investment continues —
there is no market exclusion. Output falls, and so does credit. The key
prediction concerns the *ratio*: because credit is the primary transmission
channel (net worth falls → leverage constraint binds → credit contracts), the
model predicts that the credit-to-GDP ratio will fall, since credit contracts
faster than output.

**The model does not assume this outcome.** If the output channel (ε_p term)
is large relative to the credit channel, credit/GDP could also rise. The
prediction is determined by the parameter calibration.

### Prediction for pre-default episodes

The model is run with the simulated spread path for pre-default episodes. In
addition to the balance-sheet channel, the model enters the autarky regime
once default occurs:

- Investment stops: I_t = 0, so the capital law of motion becomes
  k_{h+1} = (1−δ) · k_h − δ (log-linear: the SS investment term δ goes missing)
- The working-capital wedge is set to dR^L_aut > 0 (autarky rate)
- The country re-enters markets at rate μ per period

The survival-weighted average output is:
```
ȳ^def_h = (1−μ)^h · [α · k^excl_h − ε_p · dR^L_aut] · 100
```
where k^excl_h follows the capital depletion path and the (1−μ)^h weight
shrinks as countries re-enter markets.

**The model does not assume the ordering of the two output paths.** The
relative depth of the default vs. non-default output losses is determined by
calibrated parameters. If the balance-sheet channel is large and the autarky
wedge is small, non-default crises could dominate. The simulation determines
which scenario matches the data.

Regarding credit/GDP under default: the model predicts that GDP falls faster
than credit (capital channel dominates). This implies that the credit-to-GDP
ratio may fall *less* under default than under non-default — the opposite
ordering. Again, this is a theoretical prediction, not an assumption.

---

## 4. Calibration Strategy: Independent Targets Only

The transmission parameters {ξ, φ, Φ_N} are calibrated by SMM to minimize
the distance between the model-implied non-default output path and the
empirical LP estimate. All other parameters are set using moments that are
entirely independent of the LP estimates:

| Parameter | Target | Source |
|-----------|--------|--------|
| β (discount factor) | R* = 4% annual | Standard |
| α (capital share) | 0.35 | National accounts |
| δ (depreciation) | 0.07 | National accounts |
| σ (risk aversion) | 2 | Standard DSGE |
| R* (risk-free rate) | 0.04 | US Treasuries average |
| θ (recovery rate) | 0.65 | Cruces-Trebesch (2013) |
| μ (re-entry probability) | 0.10 | Arellano (2008) |
| ρ_y, σ_ε (income process) | HP-filtered EM output moments | IMF WEO |
| h_def (default cost) | Default frequency 2–4% annual | EM data |
| β' (sovereign discount) | Mean spread 200–400 bps | EMBIG data |
| b^B/GDP (bank sovereign exposure) | Bank claims on govt / GDP | BIS, IMF FSI |
| N/GDP (bank net worth ratio) | Bank equity / GDP | BIS |
| credit/GDP | Private credit / GDP | World Bank |
| λ (leverage multiplier) | credit/GDP ÷ N/GDP | Derived |
| γ (sovereign ceiling pass-through) | 0.80 | Bocola (2016) |

The SMM step (fitting {ξ, φ, Φ_N} to the non-default output IRF) is the only
place the empirical IRF enters calibration. The default path is then fully
out-of-sample: the parameters are fixed and the model predicts the default
output path without further tuning.

---

## 5. The Genuine Empirical Test

The model passes two tests if it is a good description of the data:

**Test 1 (non-default, in-sample):** The SMM-fitted model should produce a
non-default output path that lies close to the LP estimate. A poor fit (large
SSE) indicates the model's channel structure is inconsistent with the data,
even after optimal parameter choice.

**Test 2 (default, out-of-sample):** With parameters fixed from Test 1, the
model generates a default output path. The shape and magnitude of this
prediction — deeper trough, slower recovery — should correspond to the
empirical LP estimate for default-linked episodes. The model could fail this
test: if the autarky channel is too large or too small, or if the recovery
speed is wrong, the model prediction will diverge from the data.

**The credit/GDP test (fully out-of-sample for both types):** Neither credit
IRF is used in calibration. The model predicts private credit responses
through the leverage amplifier (IR.4). The ranking — non-default credit/GDP
falls more than default credit/GDP — must emerge from the model's internal
mechanics, not from parameter tuning.

---

## 6. What the Model Does Not Assume

To be explicit, the following findings are NOT imposed by the model structure:

- That non-default crises cause meaningful output losses (the model could
  generate near-zero effects if transmission parameters are small).
- That default crises are deeper than non-default (the balance-sheet channel
  could dominate the capital-depletion channel at the calibrated parameters).
- That the credit-to-GDP ratio falls more in non-default episodes (this
  follows from the relative size of the credit vs. output channels).
- That the default path is more persistent (recovery depends on μ and δ,
  calibrated independently).

All of these are genuine predictions. The model is a description of mechanisms
that operate in any emerging economy with sovereign risk and financial
intermediaries. Whether those mechanisms generate results consistent with the
LP estimates is an empirical question that the model answers ex ante.

---

## 7. Key Equations Summary

### Block 1 — Sovereign Problem

```
V^R(B, y) = max_{B'} { u(y − B + q(B',y)·B') + β' E[V(B',y')] }
V^D(y)    = u(y^def) + β'{ μ·E[V(0,y')] + (1−μ)·E[V^D(y')] }
V(B, y)   = max{ V^R(B,y), V^D(y) }
q(B',y)   = [1−δ(B',y)·(1−θ)] / (1+R*)
```

where δ(B',y) = E[d(B',y')|y] is the next-period default probability.

### Block 2 — Transmission Recursion

```
s_h  → (balance-sheet shock) → n_h → (leverage constraint) → ℓ_h
s_h  → (working-capital wedge) → y_h
s_h  → (investment cost) → k_h → y_h
```

Full log-linear system (deviation from steady state, h = years after onset):

| Equation | Variable | Formula |
|----------|----------|---------|
| LL.5a | Net worth | n_h = Φ_N·n_{h-1} − B_sens·s_h − Φ^O·γ·s_{h-1} |
| LL.3a | Capital | k_h = (1−δ)·k_{h-1} − η·(γ·s_{h-1} − Ω·n_{h-1}) |
| LL.2a | Output | y_h = α·k_h − ε_p·(γ·s_h − Ω·n_h) |
| IR.4  | Credit | ℓ_h = (λ·N_ss/ℓ_ss)·n_h |

η = δ / [(1−α)·(R^L_ss − (1−δ))]  
ε_p = ξ·(1−α)/α·R^L_ss / [1 + ξ·(R^L_ss − 1)]  
Ω = φ·b^B/N  
B_sens = (b^B/N) / (1 + s_ss)²

### Block 2 — Default Regime (Autarky)

```
k^excl_h = (1−δ)·k^excl_{h-1} − δ      (I = 0: full investment stop)
ȳ^def_h  = (1−μ)^h · [α·k^excl_h − ε_p·dR^L_aut] · 100
```

dR^L_aut is pinned by the h=0 output observation; all h≥1 are out-of-sample.

---

## 8. Implementation in the Pipeline

| Script | Role |
|--------|------|
| `do/14_calibration.do` | Sets all independent calibration parameters |
| `do/15_solve_default.do` | Solves VFI; simulates economy; computes event-study spread paths separately for non-default and pre-default episodes |
| `do/16_model_irf.do` | Runs transmission block on both spread paths; SMM-fits {ξ,φ,Φ_N} to non-default output IRF; generates default path fully out-of-sample; produces model-vs-data figures |

The separation of non-default vs. pre-default spread paths in
`15_solve_default.do` is critical for the test to be genuine: the simulation
must identify the two types of episodes *from the model*, not from the data.

---

## 9. References

- Arellano, C. (2008). Default Risk and Income Fluctuations in Emerging Economies. *AER*.
- Bocola, L. (2016). The Pass-Through of Sovereign Risk. *JPE*.
- Cruces, J. & Trebesch, C. (2013). Sovereign Defaults: The Price of Haircuts. *AEJ Macro*.
- Jordà, O. (2005). Estimation and Inference of Impulse Responses by Local Projections. *AER*.
- Mendoza, E. (2010). Sudden Stops, Financial Crises, and Leverage. *AER*.
