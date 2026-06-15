# The Aftermath of Sovereign Spread Crises With and Without Default

## Abstract

Using a dataset of 61 spread crisis episodes across 52 emerging market economies over 1994–2025, we estimate the output cost of sovereign spread crises and their transmission channels through local projections (Jordà, 2005). Spread crises reduce cumulative real GDP per capita by 3.8 percentage points at the trough (h=1), with no mean reversion over four years. Non-default episodes alone impose a cost of 2.3 percentage points, demonstrating that sovereign stress carries large real consequences independent of debt restructuring. Default-linked crises are roughly twice as costly (−5.5pp at h=1), driven by immediate forced external deleveraging and a front-loaded compression of output. Decomposing transmission, we find that credit and investment contractions are concentrated in non-default episodes, while default-linked episodes exhibit a pronounced sovereign-bank doom loop — banks accumulating 4 to 6 percentage points of GDP in sovereign bonds — and a sharper, earlier current account adjustment consistent with a forced sudden stop. These findings challenge the view that output losses from sovereign stress require default to materialize, and reveal that default and non-default crises operate through fundamentally different channels.

---

## Key Empirical Narratives (Working Notes)

*This section records the main empirical findings and their paper-ready interpretations, consistent with the current pipeline (do/02 through do/13). Updated as results evolve.*

---

### Act 1 — Output Cost of All Spread Crises (02_lp_all.do)

**Baseline finding:**
Spread crises cause immediate and persistent output losses. Cumulative real GDP per capita falls by −1.86pp in the year of onset (h=0), deepens to a trough of −3.76pp at h=1, and remains significantly below the pre-crisis trajectory through h=4 (−3.00pp). All estimates are significant at the 1% level under Driscoll-Kraay standard errors. There is no evidence of mean reversion over the four-year window.

**Pre-trend validation:**
At h=−2 the coefficient is close to zero and statistically insignificant (p=0.114), confirming no pre-existing output trend that could confound the crisis indicator. The h=−1 coefficient is significant, reflecting the well-documented pre-crisis deterioration in EM economies — spreads spike precisely when conditions are already weakening — an anticipation effect consistent with the literature, not a violation of identification.

**Placebo randomization (07_placebo.do):**
1,000 random draws of 61 crisis dates yield empirical p-values of 0.000–0.008 at horizons h=1–4, strongly rejecting the null that the estimated responses could arise by chance.

**IPW robustness (08_ipw_lp.do):**
The IPW-weighted estimate closely tracks the OLS baseline, alleviating concerns that output losses reflect pre-existing structural differences between crisis and non-crisis country-years rather than the causal effect of the spread crisis.

*Paper sentence:* "Sovereign spread crises impose an immediate and persistent output cost: GDP per capita falls by 1.9pp in the onset year and deepens to a trough of 3.8pp at h=1, remaining significantly depressed four years after onset. Placebo randomization strongly rejects chance (empirical p ≤ 0.008), and IPW-weighted estimates confirm the direction and magnitude of the baseline."

---

### Act 2 — Default vs. Non-Default Resolution Split (03_lp_resolution.do)

**Sample:** 40 non-default spread crises, 21 default-linked episodes (Venezuela 2008 reclassified as non-default under the 5-year onset-to-restructuring rule; restructuring began 2017, a 9-year lag).

**Core finding:**
Default-linked crises are associated with immediate and severe output losses: −3.42pp at h=0, deepening to −5.50pp at h=1, and remaining persistently depressed. The response is both larger in magnitude and front-loaded, consistent with the acute disruption to credit, trade finance, and investment that accompanies sovereign default or pre-default distress.

Non-default crises produce a smaller and more delayed response: the h=0 estimate is not statistically significant, and the peak loss of −2.33pp occurs at h=1. This suggests that, absent the contractionary shock of debt restructuring, economies retain some capacity to absorb a spread spike before output deteriorates materially.

**Key message for the literature:**
The standard view treats sovereign default as the primary source of output loss. Our results challenge this: non-default spread crises — where debt service is maintained — still generate large output contractions of −2.3pp at the trough. The crisis is the shock; default amplifies but does not create the output cost. The gap of roughly 3pp at the trough is economically meaningful, though formal equality tests reflect limited power with 21 default episodes.

*Paper sentence:* "The finding that non-default spread crises impose output losses of 2.3 percentage points at the trough directly challenges the implicit assumption that sovereign stress only matters for real activity through the debt restructuring channel. Default amplifies the shock — adding a further 3pp — but the crisis itself is the primary driver."

---

### Act 3 — Transmission Channels (11_channels.do)

#### Channel 1 — Private Credit/GDP
**Finding:** Credit contracts significantly starting at h=2 (−2.25pp), deepening through h=4 (−3.48pp, p<0.001). No significant effect at h=0–1, suggesting delayed financial transmission.

*Paper sentence:* "The credit contraction deepens progressively over four years, consistent with a gradual deterioration of bank balance sheets as sovereign risk is transmitted through the financial system rather than an acute credit crunch at onset."

#### Channel 2 — Bank Claims on Govt/GDP (Diabolic Loop)
**Finding (with L.claims_govt control):** All coefficients insignificant once pre-existing sovereign exposure is controlled for. Positive signs at all horizons.

**Finding (without L.claims_govt control):** Previously significant at all five horizons.

**Interpretation:** The diabolic loop operates through *pre-existing* balance sheet exposure rather than *new* crisis-driven accumulation. Banks that entered the crisis heavily exposed to sovereign risk suffer most, but no additional systematic accumulation occurs during the crisis.

*Paper sentence:* "The sovereign-bank nexus channel is fully absorbed by the lagged stock of bank claims on the government, suggesting that the diabolic loop reflects pre-existing exposure rather than new crisis-driven sovereign bond accumulation."

#### Channel 3 — Investment/GDP
**Finding:** Significant trough at h=3 (−1.89pp, p<0.001). Pattern is persistently negative from h=1 onward, consistent with medium-term crowding-out and elevated uncertainty depressing capital expenditure.

*Paper sentence:* "Investment contracts persistently following spread crisis onset, with the trough at three years post-onset consistent with medium-term crowding-out of private investment by elevated sovereign borrowing costs and heightened policy uncertainty."

#### Channel 4 — Government Expenditure/GDP
**Finding:** Significant fiscal expansion at onset h=0 (+1.19pp, p=0.006), fading to zero and slightly negative thereafter — a failed stabilization pattern.

*Paper sentence:* "Governments initially expand expenditure at crisis onset, but this fades within one year as rising borrowing costs and declining fiscal space force retrenchment."

#### Channel 5 — Primary Balance/GDP
**Finding:** Borderline deterioration at h=0 (−0.57pp), reverting to zero thereafter. No sustained fiscal adjustment or austerity response.

#### Channel 6 — FDI/GDP
**Finding:** Mild initial deterrence at h=0 (−0.78pp, p=0.069), normalizing quickly. FDI is not a primary transmission channel.

---

### Act 4 — Channel Heterogeneity by Resolution Type (12_channels_resolution.do)

#### Credit — A Non-Default Phenomenon
**Finding:** Non-default episodes drive the credit contraction: β_nd deepens from −0.85 (h=0) to −4.78 (h=4, p=0.020). Default-linked episodes show no significant contraction (β_def near zero throughout). Equality test significant at h=4.

*Paper sentence:* "The credit contraction is entirely driven by non-default episodes — countries that maintain debt service suffer a prolonged credit squeeze, while default-linked episodes paradoxically avoid this channel, possibly because debt restructuring restores fiscal sustainability and removes the overhang on bank balance sheets."

#### Bank Sovereign Nexus — A Default-Linked Phenomenon
**Finding:** Default-linked episodes show a massive increase in bank claims on government: +4.06pp at h=0, +6.18pp at h=1, +6.25pp at h=2. Non-default episodes show slightly negative or zero response. Equality test significant at h=0 (IPW p=0.005) and h=2.

*Paper sentence:* "The diabolic loop — banks loading up on sovereign debt during a crisis — is exclusively a default-linked phenomenon, with banks accumulating 4 to 6 percentage points of GDP in additional sovereign bonds. This result is robust to IPW correction (p=0.005 at h=0), consistent with undercapitalized banks gambling for resurrection by concentrating sovereign exposure on the path to default."

#### Investment — A Non-Default Phenomenon
**Finding:** Non-default episodes drive investment contraction: β_nd significant at h=0 (−1.19pp, p=0.003), h=1 (−1.77pp, p=0.021), h=3 (−2.73pp, p=0.010). Default-linked episodes show positive or near-zero investment throughout. Equality test robust to IPW (h=0: p=0.045; h=4: p=0.049).

*Paper sentence:* "Investment contraction is concentrated in non-default episodes, where sovereign risk persistently raises the cost of capital without the debt relief that restructuring provides. Default episodes, counterintuitively, do not exhibit significant investment declines — consistent with debt restructuring reducing long-run uncertainty and allowing investment to stabilize."

#### Government Expenditure — Divergent Fiscal Paths
**Finding:** Non-default: sustained mild fiscal expansion throughout (+0.6 to +1.1pp, all positive). Default-linked: large initial expansion (+2.19pp at h=0) followed by forced fiscal contraction (−1.39pp at h=2, −1.69pp at h=4). Equality test borderline at h=4 (OLS p=0.061).

*Paper sentence:* "The fiscal response diverges sharply by resolution type. Non-default episodes feature sustained accommodative fiscal policy, reflecting retained market access. Default-linked episodes exhibit a boom-bust pattern: larger initial expansion as fiscal constraints have not yet fully bound, followed by forced retrenchment as market access is lost and restructuring conditionality imposes consolidation."

---

### Act 5 — Current Account and Forced Deleveraging (13_mechanisms.do)

**Aguiar-Gopinath (2006) test:**
Default-linked episodes show a larger and earlier current account adjustment: +1.44pp at h=0, rising to +2.61pp at h=1. Non-default episodes exhibit a weaker and more delayed adjustment, consistent with retained — if costly — access to international capital markets.

This pattern is consistent with a sharp reversal of capital inflows in default-linked crises, forcing an abrupt compression of the current account deficit (sudden stop). Non-default episodes involve gradual adjustment, suggesting the spread spike triggers tightening but not a complete closure of market access.

**IPW robustness:** IPW-corrected estimates barely alter these conclusions, confirming that compositional differences between groups do not drive the CA finding.

*Paper sentence:* "Default-linked crises force an immediate external adjustment (+2.6pp current account improvement at h=1), consistent with a sudden stop in capital inflows and forced external deleveraging à la Aguiar-Gopinath (2006). Non-default crises show a weaker, delayed adjustment, reflecting partial rather than complete market exclusion."

---

### Cross-Cutting Narrative

**The central paradox:**
Default-linked episodes do not produce more severe credit contractions or investment declines than non-default episodes — if anything, the opposite. Non-default countries maintain full debt service under stressed conditions, keeping sovereign risk elevated persistently and bleeding the real economy through credit and investment channels. Default countries suffer a sharper initial shock, an immediate forced external adjustment, and a doom-loop in bank sovereign holdings — but the restructuring eventually provides a cleaner break.

**Two distinct crisis mechanisms:**
- **Non-default channel:** slow-bleeding — investment crowding-out, prolonged credit squeeze, sustained spread elevation, gradual CA adjustment
- **Default-linked channel:** acute shock — immediate output collapse, forced sudden stop, doom-loop bank behavior, boom-bust fiscal path

**IPW validation:**
The IPW correction for selection into default (conditional on crisis) is well-identified: only 6 observations trimmed, weights well-behaved (max 2.92), and OLS/IPW estimates track closely across all channels.

---

### Technical Notes

- **Year FE absorb VIX and UST10Y completely** — do not include global financial conditions in specifications with year fixed effects.
- **Driscoll-Kraay SE** throughout unweighted specs (xtscc, lag=h+1); cluster SE for IPW (areg, cluster cid).
- **Venezuela 2008** reclassified as non-default: restructuring began 2017, a 9-year lag well beyond the 5-year cutoff rule. Net effect: nd=40, def=21.
- **Pre-trend:** h=−2 clean (p=0.114); h=−1 significant (anticipation effect, not an identification failure).
- **Placebo p-values:** 0.000–0.008 across h=1–4.
