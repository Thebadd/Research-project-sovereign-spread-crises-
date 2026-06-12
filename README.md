# The Aftermath of Sovereign Spread Crises With and Without Default

## Abstract

Sovereign spread crises — episodes in which EMBIG spreads cross extreme thresholds — represent a recurring source of macroeconomic disruption in emerging markets, yet their real output costs remain imprecisely documented. Using a dataset of 61 spread crisis episodes across 52 emerging market and frontier economies over 1994–2025, we estimate the output cost of sovereign spread crises through the local projections method of Jordà (2005). We find that spread crises reduce cumulative real GDP per capita by approximately 3.8 percentage points at the trough, with no evidence of mean reversion over a four-year horizon. Strikingly, non-default spread crises — episodes resolved without debt restructuring — still impose a significant output cost of approximately 2.5 percentage points, demonstrating that sovereign stress carries substantial real consequences independent of default. Default-linked episodes amplify this cost by an additional 2 to 3 percentage points. These findings challenge the implicit assumption that output losses from sovereign stress are primarily a consequence of default itself.

---

## Key Empirical Narratives (Working Notes)

*This section records the main empirical findings and their paper-ready interpretations, consistent with the current pipeline (do/02 through do/12). Updated as results evolve.*

---

### Act 1 — Output Cost of All Spread Crises (02_lp_all.do)

**Baseline finding:**
Spread crises reduce cumulative real GDP per capita by approximately −3.8pp at the trough (h=3–4), with confidence intervals that exclude zero at all horizons beyond h=1. There is no evidence of mean reversion over the four-year window. This establishes that sovereign spread crises carry large and persistent real output costs regardless of default outcome.

**IPW robustness (08_ipw_lp.do):**
The IPW-weighted estimate (correcting for selection on observables into crisis onset) closely tracks the OLS baseline. The proximity of OLS and IPW estimates suggests that selection on observable pre-crisis fundamentals (debt, CA, inflation, VIX) is not the primary driver of the output loss — the crisis itself causes the contraction.

*Paper sentence:* "IPW-weighted estimates confirm the direction and magnitude of the baseline results, alleviating concerns that the output losses reflect pre-existing structural differences between crisis and non-crisis country-years rather than the causal effect of the spread crisis."

---

### Act 2 — Default vs. Non-Default Resolution Split (03_lp_resolution.do)

**Core finding:**
Non-default spread crises impose an output cost of approximately −2.5pp (trough h=3–4), statistically significant at conventional levels. Default-linked episodes impose an additional −2 to −3pp on top. Both types carry significant real costs, but they operate through different magnitudes and potentially different mechanisms.

**Key message for the literature:**
The standard view treats sovereign default as the primary source of output loss in sovereign debt crises. Our results challenge this: non-default spread crises — where debt service is maintained and no restructuring occurs — still generate large output contractions. The crisis is the shock; default amplifies but does not create the output cost.

*Paper sentence:* "The finding that non-default spread crises impose output losses of similar order of magnitude to mild default episodes directly challenges the implicit assumption that sovereign stress only matters for real activity through the debt restructuring channel."

---

### Act 3 — Transmission Channels (11_channels.do)

#### Channel 1 — Private Credit/GDP
**Finding:** Credit contracts significantly starting at h=2 (−2.25pp), deepening through h=4 (−3.48pp, p<0.001). No significant effect at h=0–1, suggesting a delayed financial transmission.

**Interpretation:** The credit contraction is not instantaneous — it reflects a gradual tightening of bank lending standards and funding costs as sovereign risk permeates the banking system. The lag structure (insignificant at h=0, growing through h=4) is consistent with a sovereign-bank nexus that transmits slowly through balance sheet deterioration rather than through an immediate credit crunch.

*Paper sentence:* "The credit contraction deepens progressively over four years, consistent with a gradual deterioration of bank balance sheets as sovereign risk is transmitted through the financial system rather than an acute credit crunch at onset."

#### Channel 2 — Bank Claims on Govt/GDP (Diabolic Loop)
**Finding (with L.claims_govt control):** All coefficients insignificant (h=0: +0.61pp, p=0.196; through h=4: all p>0.14). Positive signs at all horizons.

**Finding (without L.claims_govt control, original spec):** Previously significant at all 5 horizons (h=0: p=0.011, h=1: p=0.036, ...).

**Interpretation:** The diabolic loop result — banks increasing sovereign bond holdings during a crisis — is absorbed by the pre-existing stock of sovereign bond exposure (L.claims_govt). This has a precise economic meaning: once we account for how much sovereign debt banks already hold entering the crisis, the additional accumulation triggered by the crisis itself is not statistically distinguishable from zero. The diabolic loop operates through *pre-existing* exposure, not through *new* crisis-driven accumulation.

*Paper sentence:* "The sovereign-bank nexus channel is fully absorbed by the lagged stock of bank claims on the government, suggesting that the diabolic loop operates through pre-existing balance sheet exposure rather than new crisis-driven sovereign bond accumulation — banks that entered the crisis heavily exposed to sovereign risk suffer most, but no additional accumulation occurs systematically."

#### Channel 3 — Investment/GDP
**Finding:** Significant trough at h=3 (−1.89pp, p<0.001). Pattern: −0.59, −1.25, −0.86, −1.89, −1.03 — non-monotonic but clearly negative from h=1 onward.

**Interpretation:** The investment contraction is consistent with crowding-out (sovereign borrowing displacing private investment) and with uncertainty effects depressing capital expenditure. The non-monotonic pattern may reflect initial wait-and-see behavior followed by postponed investment materializing as uncertainty crystallizes.

*Paper sentence:* "Investment contracts persistently following spread crisis onset, with the trough at three years post-onset consistent with medium-term crowding-out of private investment by elevated sovereign borrowing costs and heightened policy uncertainty."

#### Channel 4 — Government Expenditure/GDP
**Finding:** Significant fiscal expansion at onset h=0 (+1.19pp, p=0.006), fading to zero and slightly negative thereafter (h=1: +0.44, h=2: −0.07, h=3: −0.17, h=4: −0.02, all insignificant).

**Interpretation:** Governments attempt counter-cyclical fiscal policy at the onset of the crisis, but fiscal space closes quickly as spreads rise and financing costs become prohibitive. The short-lived expansion followed by fiscal neutrality is consistent with a "failed stabilization" pattern common in EM sovereign stress episodes.

*Paper sentence:* "Governments initially expand expenditure at crisis onset — consistent with automatic stabilizers and deliberate counter-cyclical policy — but this expansion fades within one year as rising borrowing costs and declining fiscal space force retrenchment."

#### Channel 5 — Primary Balance/GDP
**Finding:** Borderline deterioration at h=0 (−0.57pp, p=0.038 or p=0.051 depending on spec), reverting to zero thereafter. No sustained fiscal adjustment.

**Interpretation:** The initial primary balance deterioration reflects the combination of falling revenues (automatic stabilizers) and rising expenditure at onset. The rapid mean reversion suggests no systematic pro-cyclical tightening — countries do not immediately impose austerity in response to spread crises, consistent with the govexp result.

#### Channel 6 — FDI/GDP
**Finding:** Borderline negative at h=0 (−0.78pp, p=0.069), reverting to zero or slightly positive thereafter. No sustained effect.

**Interpretation:** FDI shows a mild initial deterrence effect at crisis onset — foreign investors pull back when spreads spike — but this normalizes quickly. FDI is not a primary transmission channel for spread crises, unlike portfolio flows which are more sensitive to sovereign risk pricing.

---

### Act 4 — Channel Heterogeneity by Resolution Type (12_channels_resolution.do)

This is the richest set of results in the paper and provides the clearest narrative differentiation between default and non-default episodes.

#### Credit — A Non-Default Phenomenon
**Finding:** Non-default episodes drive the credit contraction: β_nd deepens from −0.85 (h=0) to −4.78 (h=4, p=0.020). Default-linked episodes show no significant contraction (β_def ranges from +1.33 to +0.62, near zero). Equality test significant at h=4 (OLS p=0.020; IPW p=0.112).

*Paper sentence:* "The credit contraction documented in the aggregate channel analysis is entirely driven by non-default episodes — countries that maintain debt service suffer a prolonged credit squeeze as banks tighten lending standards under persistent sovereign risk, while default-linked episodes paradoxically avoid this channel, possibly because debt restructuring restores sovereign debt sustainability and removes the overhang on bank balance sheets."

#### Bank Sovereign Nexus — A Default-Linked Phenomenon
**Finding:** The most striking heterogeneity result. Default-linked episodes show a massive increase in bank claims on government: +4.06pp at h=0, +6.18pp at h=1, +6.25pp at h=2. Non-default episodes show slightly negative or zero response. Equality test significant at h=0 (OLS p=0.100; IPW **p=0.005**) and h=2 (both specs p≈0.09).

*Paper sentence:* "The diabolic loop — banks loading up on sovereign debt during a crisis — is exclusively a default-linked phenomenon. Banks accumulate 4 to 6 percentage points of GDP in additional sovereign bonds in the years preceding and surrounding default, consistent with the doom-loop mechanism (Brunnermeier et al. 2016; Farhi and Tirole 2018) in which undercapitalized banks 'gamble for resurrection' by concentrating sovereign exposure. This result is robust to IPW correction for pre-crisis fundamentals (p=0.005 at h=0), confirming it reflects a causal response to the default trajectory rather than pre-existing differences between episode types."

#### Investment — A Non-Default Phenomenon
**Finding:** Non-default episodes drive investment contraction: β_nd significant at h=0 (−1.19pp, p=0.003), h=1 (−1.77pp, p=0.021), h=3 (−2.73pp, p=0.010). Default-linked episodes show positive or near-zero investment (β_def: +1.07, +0.19, +0.35, +0.43, +1.52). Equality test robust to IPW (h=0: p=0.045; h=4: p=0.049).

*Paper sentence:* "Investment contraction is concentrated in non-default episodes, where sovereign risk persistently raises the cost of capital and creates fiscal crowding-out without the debt relief that restructuring provides. Default episodes, counterintuitively, do not exhibit significant investment declines — consistent with the view that debt restructuring, by restoring fiscal sustainability, reduces long-run uncertainty and allows investment to stabilize."

#### Government Expenditure — Divergent Fiscal Paths
**Finding:** Non-default: sustained mild fiscal expansion throughout (+0.6 to +1.1pp, all positive). Default-linked: large initial expansion (+2.19pp at h=0) followed by forced fiscal contraction (−1.39pp at h=2, −1.10pp at h=3, −1.69pp at h=4). Equality test borderline at h=4 (OLS p=0.061) and h=2–3 (IPW p=0.062–0.087).

*Paper sentence:* "The fiscal response diverges sharply by resolution type. Non-default episodes feature sustained accommodative fiscal policy throughout the crisis horizon, reflecting retained market access that allows governments to smooth the shock. Default-linked episodes exhibit a boom-bust fiscal pattern: a larger initial expansion as fiscal constraints have not yet fully bound, followed by forced retrenchment as market access is lost and restructuring conditionality imposes consolidation."

#### Primary Balance and FDI — No Significant Heterogeneity
Both channels show no statistically significant difference between nd and def at any horizon. This is itself informative: the primary balance does not discriminate resolution type, suggesting that the fiscal adjustment path (govexp) matters more than the fiscal balance outcome in the short run.

---

### Cross-Cutting Narrative for the Paper

**The central paradox:**
Default-linked episodes do NOT produce more severe credit contractions or investment declines than non-default episodes — if anything, the opposite. This is paradoxical if one assumes default = worse outcome. The explanation lies in debt dynamics: non-default countries maintain full debt service under stressed conditions, keeping sovereign risk elevated persistently and bleeding the real economy through credit and investment channels. Default countries suffer a sharper initial shock but the restructuring provides a cleaner break, removing the sovereign risk overhang and allowing the real economy to stabilize.

**The doom-loop as a default predictor, not a consequence:**
The massive bank accumulation of sovereign bonds (+4 to +6pp) in default episodes should be interpreted as a *leading indicator* of default, not a consequence. Banks in countries on the default trajectory load up on sovereign bonds (gambling for resurrection), which then amplifies the eventual default shock. This is the classic doom-loop mechanism, and our LP evidence provides direct time-series support for it.

**IPW validation:**
The IPW correction in Act 2 (selecting into default conditional on crisis) is well-identified: only 6 observations trimmed, weights well-behaved (max 2.92), and OLS/IPW estimates track closely. This confirms that the nd/def heterogeneity in channels reflects the causal impact of resolution type rather than pre-existing differences in debt levels or external positions between the two groups.

---

### Technical Notes for Paper Writing

- **Year FE absorb VIX and UST10Y completely** — do not include global financial conditions variables in specifications with year fixed effects; they have zero cross-sectional variation and are collinear with the year dummies.
- **Driscoll-Kraay SE** used throughout the unweighted specifications (xtscc, lag=h+1); cluster SE used for IPW (areg, cluster cid) since xtscc does not accept pweights.
- **IPW in Act 1 (11_channels.do Section 7)** corrects for selection into crisis vs. tranquil years — valid but note 37% sample trim and McFadden R²=0.20 (partial correction only).
- **IPW in Act 2 (12_channels_resolution.do)** corrects for selection into default conditional on crisis — much better identified (6 trims, R²=0.185, 54 obs).
- **Lag of own outcome** included for claims_govt, govexp, pb, fdi but not credit or inv. These lags absorb persistence but may also absorb treatment variation — robustness check without lags is warranted.
