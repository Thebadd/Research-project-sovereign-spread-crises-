# 3. Theoretical Framework

The framework consists of two blocks. The first is a quantitative sovereign default model in the tradition of Arellano (2008), which determines the equilibrium bond price schedule and classifies sovereign episodes as either resolved without default or culminating in debt restructuring, according to the government's optimal default policy $`d(B_t, y_t)`$. The second is a log-linear transmission block, adapted from Bocola (2016), which maps a sovereign spread path into output and credit responses through three channels: a balance-sheet channel operating through bank net worth, a working-capital channel operating through firm financing costs, and a capital-accumulation channel operating through investment. The two blocks are recursive: Block 1 supplies the spread paths; Block 2 is evaluated conditional on those paths, abstracting from feedback from real outcomes to sovereign risk. This formulation is consistent with a local projection identification strategy in which the sovereign spread is the forcing variable and real responses are estimated conditional on its realization.

---

## 3.1 Economic Environment

Time is discrete, indexed by $`t`$. The economy is small and open, facing world risk-free rate $`R^*`$, and is populated by a representative household, a government, a unit measure of competitive banks, and a representative firm. The household supplies labor inelastically and saves through bank deposits. The government issues one-period discount bonds to smooth consumption, subject to the option to default. Banks hold sovereign bonds and extend loans to firms, subject to a leverage constraint. Firms produce using capital and labor and borrow working capital from banks to finance a fraction of their wage bill.

A persistent endowment shock $`y_t`$, following a log-AR(1) process, drives the government's default incentive. The equilibrium sovereign spread $`s_t`$, determined in Block 1, transmits to the real economy through the banking and production system of Block 2.

---

## 3.2 Sovereign Default Problem

The government maximizes household welfare, discounting at rate $`\beta' < \beta`$, where $`\beta`$ is the household's discount factor (Arellano 2008; Aguiar and Gopinath 2006). The government issues one-period discount bonds $`B_t`$ at price $`q_t`$ and decides each period whether to repay or default. The value function satisfies:

$$V(B_t, y_t) = \max\bigl\lbrace V^R(B_t, y_t),\; V^D(y_t)\bigr\rbrace$$

**Repayment.** Given the bond price schedule $`q(B_{t+1}, y_t)`$, the government chooses next-period debt $`B_{t+1}`$ to solve:

$$V^R(B_t, y_t) = \max_{B_{t+1}} \left\lbrace u\left(y_t - B_t + q(B_{t+1}, y_t)\, B_{t+1}\right) + \beta'\,\mathbb{E}_{y'|y}\left[V(B_{t+1}, y')\right] \right\rbrace$$

where $`u(c) = c^{1-\sigma}/(1-\sigma)`$ and the budget constraint is $`c_t = y_t - B_t + q(B_{t+1}, y_t)\,B_{t+1}`$.

**Default.** Following default, the sovereign receives autarky income $`y_t^{def} = \min(y_t, \bar{h}\,\bar{y})`$, where $`\bar{y} = \mathbb{E}[y]`$ and $`\bar{h} \in (0,1)`$, and re-enters international markets with probability $`\mu`$ per period:

$$V^D(y_t) = u(y_t^{def}) + \beta'\left\lbrace\mu\,\mathbb{E}_{y'|y}\left[V(0,y')\right] + (1-\mu)\,\mathbb{E}_{y'|y}\left[V^D(y')\right]\right\rbrace$$

**Default policy.** The optimal rule is:

$$d(B_t, y_t) = \mathbf{1}\left\lbrace V^D(y_t) > V^R(B_t, y_t)\right\rbrace$$

Because $`V^R`$ is strictly increasing in $`y_t`$ while $`V^D`$ is independent of $`B_t`$, the default set $`\mathcal{D}(B_t) = \{y_t : d(B_t, y_t) = 1\}`$ takes the form of a lower interval $`(-\infty,\,\hat{y}(B_t)]`$. The same rule $`d(B_t, y_t)`$ generates both crisis types endogenously: depending on the realizations of $`(B_t, y_t)`$ along the equilibrium path, a high-spread episode may be resolved without default or may culminate in restructuring.

---

## 3.3 Bond Pricing and the Sovereign Spread

Foreign investors are risk-neutral and price bonds at the zero-profit level. A bond pays one unit upon repayment and $`\theta \in (0,1)`$ units upon default ($`\theta`$ is the recovery rate). The equilibrium bond price is:

$$q(B_{t+1}, y_t) = \frac{1 - \delta(B_{t+1}, y_t)\,(1-\theta)}{1 + R^*}$$

where $`\delta(B_{t+1}, y_t) = \mathbb{E}_{y'|y}[d(B_{t+1}, y')]`$ is the next-period default probability. The sovereign spread is:

$$s_t = \frac{1}{q(B_{t+1}, y_t)} - 1 - R^*$$

Higher debt or lower income raises the default probability, compresses $`q`$, and widens the spread. Partial recovery bounds the spread as the default probability approaches one.

**Equilibrium.** An equilibrium consists of value functions $`\{V^R, V^D, V\}`$, a bond price schedule $`q(B', y)`$, a default policy $`d(B, y)`$, and a debt policy $`B'(B, y)`$ such that policies are optimal given $`q`$ and $`q`$ reflects the default probability implied by $`d`$. The equilibrium is computed by value function iteration over a discretized state space.

---

## 3.4 Banking Sector

Each bank holds sovereign bonds $`b^B_t`$ and firm loans $`\ell_t`$, funded by net worth $`N_t`$ and short-term deposits $`d_t`$:

$$b^B_t + \ell_t = N_t + d_t$$

**Net worth dynamics.** A widening sovereign spread reduces the market value of sovereign bond holdings, inflicting a mark-to-market loss on bank equity. For a bank holding $`b^B_{t-1}`$ units purchased at price $`q_{t-1}`$, the capital loss is $`b^B_{t-1}\Delta q_t`$ when the bond price falls by $`\Delta q_t`$. Log-linearizing net worth around the steady state:

$$\hat{n}_t = \Phi_N\,\hat{n}_{t-1} - B_{sens}\,\hat{s}_t - \Phi^O\gamma\,\hat{s}_{t-1} \quad\text{(1)}$$

where $`\hat{n}_t = (N_t - N_{ss})/N_{ss}`$ and $`\hat{s}_t = s_t - s_{ss}`$. The balance-sheet sensitivity $`B_{sens} = (b^B_{ss}/N_{ss})/(1+s_{ss})^2`$ is proportional to the sovereign-bond-to-equity ratio at steady state. The parameter $`\Phi_N \in (0,1)`$ governs the persistence of net worth deviations through retained earnings and recapitalization. The term $`\Phi^O\gamma\hat{s}_{t-1}`$ captures the delayed repricing of wholesale funding costs.

**Leverage constraint.** A binding leverage constraint limits lending to a multiple of net worth (Gertler and Karadi 2011; Bocola 2016):

$$\ell_t = \lambda N_t \quad\text{(2)}$$

where $`\lambda > 1`$ is the leverage multiplier. A spread-induced decline in $`N_t`$ forces a proportional contraction in firm credit.

**Sovereign ceiling.** Following the empirical regularity that domestic banks in emerging markets cannot borrow below the sovereign rate (Borensztein et al. 2006), the lending rate satisfies:

$$R^L_t = R^* + \gamma\,s_t \quad\text{(3)}$$

where $`\gamma \in (0,1]`$ is the pass-through coefficient. Banks maintain a fixed sovereign bond stock $`b^B_{ss}`$, so the sensitivity $`B_{sens}`$ is identical across both crisis regimes. When net worth deviates from steady state, an additional term $`\Omega\hat{n}_t`$, where $`\Omega = \phi\,b^B_{ss}/N_{ss}`$, enters the effective lending rate.

---

## 3.5 Firm Production

The representative firm produces with a Cobb-Douglas technology:

$$Y_t = A\,K_t^\alpha\,L_t^{1-\alpha}$$

where $`\alpha`$ is the capital share and $`A`$ is normalized to one in steady state.

**Working-capital channel.** Following Neumeyer and Perri (2005), a fraction $`\xi`$ of the wage bill is financed in advance at the bank lending rate $`R^L_t`$. The effective unit labor cost is:

$$\tilde{w}_t = w_t\bigl[1 + \xi(R^L_t - 1)\bigr]$$

Log-linearizing, the elasticity of output with respect to the lending rate is:

$$\varepsilon_p = \frac{\xi\,(1-\alpha)/\alpha\cdot R^L_{ss}}{1 + \xi(R^L_{ss} - 1)}$$

This elasticity is increasing in $`\xi`$ and in $`R^L_{ss}`$.

**Capital accumulation.** Investment satisfies the zero-profit condition for capital:

$$\alpha\,\frac{Y_t}{K_t} = R^L_t - (1-\delta)$$

Log-linearizing and combining with the capital accumulation identity $`K_{t+1} = (1-\delta)K_t + I_t`$ yields:

$$\hat{k}_{t+1} = (1-\delta)\,\hat{k}_t - \eta\left(\gamma\hat{s}_t - \Omega\hat{n}_t\right) \quad\text{(4)}$$

where $`\hat{k}_t = (K_t - K_{ss})/K_{ss}`$ and the investment sensitivity is:

$$\eta = \frac{\delta}{(1-\alpha)(R^L_{ss} - (1-\delta))}$$

Any spread-induced investment decline at horizon $`h`$ accumulates in the capital stock at $`h+1`$ and depresses output independently of the contemporaneous spread, generating inertia beyond the spread episode itself.

---

## 3.6 Transmission System

The transmission block maps an exogenous spread path $`\{\hat{s}_h\}_{h=0}^{H}`$ into output, credit, and capital responses, with $`h`$ denoting years after crisis onset. Combining the banking and production blocks, output is determined by the capital stock and the net-worth-adjusted lending wedge:

$$\hat{y}_h = \alpha\hat{k}_h - \varepsilon_p\left(\gamma\hat{s}_h - \Omega\hat{n}_h\right) \quad\text{(5)}$$

while net worth, capital, and credit evolve according to (1), (4), and (2), with the period index $`t`$ replaced by the horizon $`h`$.

At impact the capital stock is predetermined ($`\hat{k}_0 = 0`$), so (1) and (5) reduce to $`\hat{n}_0 = -B_{sens}\,\hat{s}_0`$ and $`\hat{y}_0 = -\varepsilon_p(\gamma\hat{s}_0 - \Omega\hat{n}_0)`$: the spread depresses output on impact purely through the working-capital channel. For $`h \geq 1`$, net worth losses propagate forward at rate $`\Phi_N`$ through (1) and translate into credit contraction through (2), while the capital recursion (4) is strictly inertial — investment depressed at $`h-1`$ lowers the capital stock at $`h`$ and feeds into output through the $`\alpha\hat{k}_h`$ term in (5), independently of the contemporaneous spread. Credit loads solely on $`\hat{n}_h`$, whereas output loads on both $`\hat{k}_h`$ and the rate wedge; this difference in loadings generates divergent credit-to-output trajectories across crisis regimes (Proposition 2). Output responses are reported in percentage points: $`y_h = \hat{y}_h \times 100`$.

---

## 3.7 Default Regime

**Partial investment stop.** Upon default, the sovereign enters financial autarky. Capital inflows cease, eliminating the external financing component of investment. Domestic financing — retained earnings, domestic bank credit, and tax-funded public investment — sustains the remaining fraction. Investment in autarky is:

$$I_t = (1-\chi)\,\delta K_t \quad\text{(6)}$$

where $`\chi \in [0,1]`$ is the external investment share, identified from balance-of-payments data over tranquil periods as the fraction of gross fixed capital formation financed by external flows. Substituting into the capital law of motion:

$$K_{t+1} = (1-\chi\delta)\,K_t$$

The autarky capital stock depreciates at effective rate $`\chi\delta`$. Log-linearizing around the steady state with initial condition $`\hat{k}^{excl}_0 = 0`$, the closed-form capital path is:

$$\hat{k}^{excl}_h = -(1-(1-\chi\delta)^h) \quad\text{(7)}$$

Capital depletion deepens monotonically with $`h`$ at a rate governed by $`\chi\delta`$.

**Autarky lending rate.** The sovereign ceiling implies an elevated domestic lending rate during autarky:

$$R^{L,aut} = R^L_{ss} + \Delta R^L_{aut}, \qquad \Delta R^L_{aut} > 0$$

This wedge generates an immediate output contraction through the working-capital channel, prior to the onset of capital depletion.

**Survival-weighted output path.** The local projection estimates the average output response pooling country-year observations across episodes still in autarky and those that have re-entered markets. The model counterpart weights the output of excluded economies by the autarky survival probability $`(1-\mu)^h`$:

$$\bar{y}^{def}_h = (1-\mu)^h\left[\alpha\hat{k}^{excl}_h - \varepsilon_p\,\Delta R^L_{aut}\right] \times 100 \quad\text{(8)}$$

The capital term $`\alpha\hat{k}^{excl}_h`$ deepens through (7); the survival weight $`(1-\mu)^h`$ governs recovery as excluded economies re-enter markets at rate $`\mu`$.

---

## 3.8 Main Results

The following propositions collect the principal implications of the system (1)–(2), (4)–(5), and (6)–(8). Each holds for any parameterization with $`\chi > 0`$, $`\Delta R^L_{aut} > 0`$, and $`\lambda > 1`$.

**Proposition 1 (Output cost ranking and divergence).** *Under a common sovereign spread path, $`\hat{y}^{def}_h \leq \hat{y}^{nd}_h`$ for all $`h \geq 0`$, with strict inequality for $`h \geq 1`$; the gap $`\hat{y}^{nd}_h - \hat{y}^{def}_h`$ is non-decreasing in $`h`$.*

*Proof.* Both regimes share the balance-sheet (1) and working-capital (5) terms under a common spread path. Relative to the without-default path, the default output (8) includes $`-\varepsilon_p\,\Delta R^L_{aut} \leq 0`$ and $`\alpha\hat{k}^{excl}_h \leq 0`$ from (7). Both terms are non-positive for all $`h`$, giving the ranking. Since $`(1-\chi\delta)^h`$ is strictly decreasing in $`h`$ for $`\chi > 0`$, $`\hat{k}^{excl}_h`$ falls monotonically and the gap is non-decreasing. $`\blacksquare`$

The divergence is the signature of the capital-depletion mechanism: the balance-sheet and working-capital channels, common to both regimes, produce output paths that converge rather than diverge over time.

**Proposition 2 (Credit-to-output ratio across regimes).** *The ratio of credit contraction to output contraction is larger in the without-default regime than in the default regime.*

*Proof.* Credit loads only on net worth via (2): $`\hat{\ell}_h = (\lambda N_{ss}/\ell_{ss})\,\hat{n}_h`$. In both regimes, output loads on net worth through (5). In the default regime, output additionally carries the autonomous depletion term $`\alpha\hat{k}^{excl}_h`$ from (7), which is absent from credit. Output therefore falls more relative to credit in the default regime, reducing the credit-to-output ratio. $`\blacksquare`$

**Proposition 3 (Sources of output persistence).** *Output losses persist beyond the spread episode in both regimes, governed by distinct parameters: $`\delta`$ and $`\eta`$ in the without-default regime; $`\mu`$ in the default regime.*

*Proof.* In the without-default regime, (4) propagates investment declines into the capital stock independently of the contemporaneous spread; output remains below trend through $`\alpha\hat{k}_h`$ in (5) even after $`s_h = s_{ss}`$, with recovery speed set by $`\delta`$ and $`\eta`$. In the default regime, (7) is autonomous of the spread, and the survival-weighted path (8) recovers at rate $`(1-\mu)^h`$, governed by $`\mu`$. $`\blacksquare`$

The two regimes produce persistence through non-overlapping mechanisms: endogenous capital inertia in the without-default regime; exogenous market exclusion in the default regime.

---

## References

Acharya, V. V. and Steffen, S. (2015). The "Greatest" Carry Trade Ever?
Understanding Eurozone Bank Risks. *Journal of Financial Economics*, 115(2),
215–236.

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

Brunnermeier, M. K., Langfield, S., Pagano, M., Reis, R., Van Nieuwerburgh,
S. and Vayanos, D. (2016). ESBies: Safety in the Tranches. *Economic Policy*,
31(85), 175–219.

Cruces, J. J. and Trebesch, C. (2013). Sovereign Defaults: The Price of
Haircuts. *American Economic Journal: Macroeconomics*, 5(3), 85–117.

Crosignani, M., Faria-e-Castro, M. and Fonseca, L. (2021). The (Unintended?)
Consequences of the Largest Liquidity Injection Ever. *Journal of Monetary
Economics*, 112, 97–112.

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
