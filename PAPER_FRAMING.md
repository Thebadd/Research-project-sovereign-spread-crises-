# Modeling this paper on Asonuma, Chamon, Erce & Sasahara (JIE 2024)

Reference: Asonuma, T., Chamon, M., Erce, A., Sasahara, A. (2024). *Costs of
sovereign debt crises: Restructuring strategies and bank intermediation.*
Journal of International Economics 152, 104002.

This note records how to frame our paper on that article, which results to
present, the table form to emulate, and the four concrete upgrades to build.
It is a roadmap — no code has been changed yet.

---

## 1. Why it is the right template

Their design is almost one-to-one with ours:

| Asonuma et al. | Our project |
|---|---|
| Local projections (Jordà 2005) | same |
| Resolution dimension: post-default / weakly / strictly preemptive | default-linked vs non-default (2-way) |
| Structural amplifier: bank intermediation (credit/GDP) | exposure heterogeneity (nexus, rollover, credit) |
| Outcomes: GDP, investment, bank credit, capital inflows, lending rate | channels: GDP, credit, investment, FDI, … |
| IPW / AIPW (Jordà–Taylor 2016) with first-stage probit | our IPW (08) |
| Restrict to ever-treated countries (AIPW drops never-treated) | we already do this |

**Their contribution is two-dimensional**: the cost of a crisis depends on
(1) the resolution strategy AND (2) a structural characteristic that amplifies
it. That is exactly our resolution × exposure thesis — the paper certifies the
design.

## 2. Framing to adopt

- **Abstract formula**: "Sovereign spread crises are associated with declines
  in GDP, investment, credit and capital flows. The intensity depends on two
  things: whether the crisis is resolved without default, and the country's
  pre-crisis exposure to [rollover / the sovereign-bank nexus]. Default-linked
  episodes are worse — and much of that difference is driven by high-exposure
  economies."
- **Lead with the interaction, not the average.** Their headline is not
  "default is worse"; it is "default is worse *and the gap concentrates where
  the amplifier is strong*." Our analog: the default vs non-default gap
  concentrates in high-nexus / high-rollover economies.
- **Sample**: restrict to ever-treated countries (already done); report the
  full triple Observations / Countries / Episodes everywhere.

## 3. Results to present

- **Multi-outcome LP IRFs overlaid by type** (GDP, investment, credit, capital
  inflows, lending rate) — our channel set. Present OLS *and* the reweighted
  estimator as robustness (their Fig 3 vs Fig 4).
- **Difference tests done properly** — see upgrade D below.
- **Heterogeneity**: they use an above/below-median subsample split of the
  amplifier (their Fig 6); we use an interaction (13b). Both are valid; consider
  showing the median-split version too, it reads more cleanly in a figure.
- **Pre-trends**: present the pre-onset dynamics figure (their Fig F1 analog).

## 4. Table form — emulate their Table 2 (main results)

The single most important presentation change. Their Table 2:

- **Panels = outcome variables** (Panel A: GDP, B: Investment, C: Bank credit,
  D–F: capital inflows, lending rate). One big table per result-family, not many
  small ones.
- **Columns = horizons** h = 1…5.
- **Rows = the types** (one coefficient + SE in parentheses below, stars).
- Under each coefficient row: an **"Observations / Countries / Episodes"** line
  (e.g. `1534/35/61`).
- A dedicated **"Differences between the estimated coefficients"** block at the
  bottom of each panel: one row per pairwise difference, showing
  `[bootstrap 95% CI]` and `(Clogg et al. z)`.
- Long explicit notes: control set, balancing, bootstrap procedure, star legend.

Their **Table 1** (probit first stage) is our first-stage template: columns =
type; rows grouped into "Predictors" and "Baseline controls"; diagnostics rows
(Observations, χ² for predictors, p-value, Area under ROC) at the bottom.

Our current Table 3 / Table 4 already have the right skeleton (panels by
channel, columns by horizon, rows by type). The upgrades to reach their standard
are (a) add Countries and Episodes next to Observations, (b) add the explicit
pairwise-difference block, (c) consolidate into one multi-panel table.

## 5. The four upgrades to build (roadmap — all approved, not yet implemented)

**A. Table 2-style main table (presentation only). [IN PROGRESS]**
Consolidate the channel results into one multi-panel table: panels = outcomes,
columns = horizons, rows = types; add an Observations/Countries/Episodes line
per coefficient and a pairwise-difference block. Touches the table-export blocks
in `do/11_channels.do`, `do/12_channels_resolution.do` (and the main
`do/02_lp_all.do`, `do/03_lp_resolution.do`). No new estimation.
Done so far: Episodes line added to Tables 1–4; a difference block
(Difference = beta(def) − beta(nd), Clogg et al. 1995 z, p(Clogg z), p(Wald))
added to the resolution tables (Table 2 in `03`, Table 4 in `12`). The Clogg z
is the analytic half of upgrade D; the bootstrap 95% CI is still to come.
Still to do: optionally fold Tables 1–4 into a single stacked multi-panel file.

**B. AIPW estimator (Jordà–Taylor 2016).**
Upgrade the IPW second stage in `do/08_ipw_lp.do` (and reuse in channel files)
to doubly-robust AIPW: combine the IPW term with the regression-adjustment term
(their Eq. 3). More credible than plain IPW; keep IPW as a robustness row.

**C. Probit predictors + ROC.**
Extend the first-stage probit (`08`, `12`) with exclusion-restriction
*predictors* alongside the controls — candidates given our data: US rate
(`ust10y`/`fedfunds`), a regional-contagion measure (share/weighted count of
other-country crises in the region), and count of past crises per country.
Report the area under the ROC curve (Stata `lroc`) as a first-stage diagnostic,
as in their Table 1.

**D. Bootstrap + Clogg difference tests.**
Replace the single Wald `p(diff)` on β_def = β_nd with (i) a bootstrap 95% CI on
the coefficient difference (`bsample`, ~1000 draws, percentile CI) and (ii) the
Clogg et al. (1995) z-statistic `(b1 − b2)/sqrt(se1^2 + se2^2)`. Report both in
the difference block of the main table.

### Suggested build order
A (presentation, immediate payoff) → D (cheap, big credibility gain) →
C (strengthens the first stage) → B (the biggest new-code item).

---

## 6. Notes / open questions for later

- Their resolution split is 3-way (post-default / weakly / strictly preemptive);
  ours is 2-way (default-linked / non-default). Consider whether a finer split
  is feasible with our episode classification, or keep 2-way for power.
- Contagion predictor needs a distance/region weighting — we have `region`
  already; a simple region-year leave-one-out crisis share is a feasible start.
- AIPW first stage drops never-treated by construction — our ever-treated
  restriction is already consistent with this.
