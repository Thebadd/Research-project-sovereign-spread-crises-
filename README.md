# The Aftermath of Sovereign Spread Crises With and Without Default

## Abstract

Using a panel of 52 emerging and frontier-market economies over 1994–2026 and 61 identified sovereign spread-crisis episodes, we estimate the output cost of a spread crisis and its transmission channels using local projections (Jordà, 2005). On average, a spread-crisis onset is followed by a real but front-loaded output shortfall — cumulative real GDP falls by roughly 2 to 3 percentage points at the one- and two-year horizon, fading within about three years. This pooled average, however, obscures the paper's central finding: splitting episodes by resolution type shows the output cost is not a generic consequence of crossing the spread threshold, but is concentrated in, and largely specific to, episodes that culminate in default. Non-default episodes carry no reliable GDP cost at any horizon; default-linked episodes carry a persistent cost of 4 to 7 percentage points at every horizon out to five years, and the gap between the two groups is itself statistically significant. This divergence survives a direct stress test against a pre-existing pre-crisis growth pattern detected in a placebo test, though the longer-horizon persistence is only confirmed at short horizons once selection into resolution type is addressed via propensity weighting. A Gelbach (2016) decomposition shows private credit and government expenditure jointly absorb an increasing share of the default-linked cost, reaching roughly 40 percent by the fourth and fifth years — the strongest quantitative candidates for the transmission mechanism, without establishing causal mediation. The sharpest mechanism evidence comes from splitting default-linked episodes by pre-crisis bank exposure to the sovereign: banks with high exposure significantly increase their government claims for three consecutive years following a default, a genuine interaction-design result, though the associated deeper output cost in that same cell does not independently clear significance given only seven treated observations. These findings indicate that sovereign spread crises are costly mainly insofar as they end in default, and that the transmission runs primarily through credit and fiscal channels rather than through investment, the primary balance, or FDI.

---

## Repository structure

| File / directory | Contents |
|---|---|
| `do/` | Stata pipeline. `10`–`18`: from-source panel build (see `DATA_SOURCES.md`). `02`–`13b`: main empirical analysis. `14`–`16`: structural model (calibration, VFI default solver, model-vs-data IRFs). |
| `data/raw/`, `data/clean/` | Source workbooks and built panels (`panel_build.dta`, `panel_lp.dta`, per-result `irf_*.dta`). |
| `output/tables/`, `output/figures/` | Generated Word/RTF tables and PDF/PNG figures. |
| `RESULTS_SECTION_DRAFT.md` | Paper-ready Results section, continuous prose. |
| `EMPIRICAL_ANALYSIS.md` | Fuller internal reference: full methodology, every coefficient table, explicit limitations. |
| `PAPER_FRAMING.md` | How the paper is positioned relative to Asonuma, Chamon, Erce & Sasahara (2024). |
| `THEORETICAL_APPROACH.md`, `THEORETICAL_MODEL.md`, `THEORETICAL_MODEL_SIMPLIFIED.md` | The structural sovereign-default model backing `14`–`16`. |
| `DATA_SOURCES.md` | Provenance of every macro series: portal, exact indicator code, transform, coverage. |

Run the full pipeline with `do "do/00_master.do"`.

---

## Key empirical findings (working notes)

*Consistent with the current pipeline (`do/02`–`do/13b`) and the verified
run this session's fixes produced. Updated as results evolve; see
`EMPIRICAL_ANALYSIS.md` for full tables and `RESULTS_SECTION_DRAFT.md` for
the paper-ready narrative.*

### The average output cost, and why it is not the headline result

`02_lp_all.do` estimates the pooled response across all 61 episodes:
cumulative real GDP is significantly lower at the one- and two-year
horizon (−2.29pp, $p=0.014$; −2.98pp, $p=0.004$), fading to
insignificance by the fourth and fifth years. We treat this as a
motivating fact rather than the paper's central estimate, because
splitting the sample by resolution type shows the pooled number conceals
far more than it summarizes.

### The resolution split is the central result

`03_lp_resolution.do` estimates non-default and default-linked onsets
jointly against tranquil years (40 non-default, 21 default-linked
episodes; Venezuela 2008 is classified as non-default since its
restructuring began nine years after onset, well past the five-year
rule). Non-default episodes show no reliable GDP cost — significant at
only one horizon out of five. Default-linked episodes show a large,
persistent cost (−4.3 to −6.7pp) significant at every horizon, and the
formal equality test between the two rejects at four of five horizons.
**The output cost of a spread crisis is concentrated in, and largely
specific to, episodes that end in default.**

**Pre-trend stress test.** A placebo test on pre-crisis growth finds a
significant positive coefficient for default-linked episodes ($p=0.003$,
a boom-before-the-bust pattern, not a pre-existing decline) and none for
non-default. Because one lag of growth already sits in the baseline
controls, this window is mechanically absorbed already; adding a second,
independent growth lag moves the default-linked coefficient by no more
than 0.1pp at any horizon, with significance preserved throughout
(`table2pt_pretrend_controlled.rtf`). This is a genuine robustness result,
not a restated caveat.

**Selection.** Reweighting by an estimated propensity score attenuates
the pooled Act-1 cost (`08_ipw_lp.do`), rather than amplifying it — the
unweighted average somewhat overstates the typical cost. For the
resolution split, a two-stage design (each type scored against tranquil
years separately) confirms the extra cost of default at the one- and
two-year horizons under both OLS and IPW, but the IPW-weighted
default-linked coefficient loses significance at longer horizons as the
thin default-linked common support shrinks further under trimming — the
resolution gap is well established short-run, only suggestive long-run.

### Transmission channels: co-movement, not established mediation

`11_channels.do` re-estimates six intermediate outcomes (private credit,
bank claims on government, investment, government expenditure, the
primary balance, FDI) in place of GDP. Private credit shows the clearest
pooled pattern — a delayed, compounding contraction reaching −5 to −7pp
by the third and fourth years — echoing the GDP result's own shape. The
other five channels clear significance at no more than one horizon each.

Splitting by resolution type (`12_channels_resolution.do`) is
considerably noisier (default-linked cells of 13–21 episodes): **of
thirty nd/def equality tests across channels and horizons, only two are
significant** (private credit at Year 2, FDI at Year 4). Several
point-estimate gaps look large but do not clear significance, and are
read as power limitations, not established differences — a materially
more cautious reading than earlier drafts of this project asserted.

**Gelbach decomposition** (`12b_gelbach_decomposition.do`) ties the
channels directly to the default-linked GDP result: adding each channel's
own change as a control to the headline resolution-split specification,
private credit and government expenditure absorb a rising share of the
default-linked coefficient, reaching ~40–44% by Years 4–5; the primary
balance and FDI absorb essentially none of it; investment sits at a
stable ~15–22%. This ranks the channels by their statistical relationship
to the GDP cost — it does not establish causal mediation.

### The sovereign-bank nexus: the sharpest mechanism evidence, and the least secure

Splitting default-linked episodes by pre-crisis bank exposure to the
sovereign (`13d_aipw_nexus_split.do`) shows banks in the high-exposure
cell significantly increase government claims for three consecutive
years post-onset — absent in the low-exposure cell and in non-default
episodes regardless of exposure. The associated output cost is roughly
2.5× deeper in the high-exposure cell, but that gap does not itself clear
significance (only 7 treated country-years). The mirror-image result for
non-default crises is cleaner: high exposure there is associated with a
cost statistically indistinguishable from zero, against a significant
cost under low exposure, and — the one place in this exercise where a
*difference* between subsamples clears significance on its own — the gap
between them is itself significant. Read together: the sovereign-bank
linkage is protective when debt is honored and a channel for reallocation
away from private lending when it is not, though this is the paper's most
interesting and least statistically secure result at once.

### The doubly-robust (AIPW) extension — pending re-verification

`08b_aipw.do` implements a doubly-robust AIPW estimator matching the
reference paper's headline design, with channel (`13c`) and nexus (`13d`)
extensions. The GDP-level AIPW headline has not been re-run inside this
session after the horizon relabeling and control-set fixes below; treat
any AIPW GDP figure predating those fixes as stale until re-confirmed.

---

## Design notes and recent fixes worth knowing about

- **Total real GDP, not per capita.** The headline outcome is total real
  GDP (aligned with Asonuma et al.), built from scratch in `18_transforms.do`;
  a per-capita variant (`dy_pc_*`) is retained as a robustness check only.
- **Horizon labeling matches the reference paper's convention**: the
  crisis year itself is **Year 1**; **Year 0** is an explicit, hard-coded
  pre-crisis baseline (always zero, never estimated); a single genuinely
  estimable pre-trend placebo point is shown at **Year −1**. Internal
  variable names (`dy_0`…`dy_4`) are unchanged — only the printed/plotted
  axis shifted.
- **The panel is rebuilt entirely from official sources** (`10`–`18`
  chain: IMF WEO, World Bank WDI/IDS, IMF MFS/BOP, FRED,
  Laeven–Valencia), replacing the earlier opaque master CSV; only the
  spread-crisis dating/classification itself comes from the project's own
  database. See `DATA_SOURCES.md` for full provenance.
- **Episode duration is not weighted.** Each episode contributes exactly
  one treated observation (its onset year); continuation years are
  excluded from the sample entirely. A crisis lasting five years
  identifies exactly as much as one lasting one year — this is the
  standard local-projection convention, not an oversight, but it does
  mean episode duration itself is not a source of variation in this
  design.
- **Common core controls**: one lag of GDP growth, debt/GDP, current
  account/GDP, years a systemic banking crisis has run (Laeven–Valencia,
  zero-filled where none is recorded), government expenditure/GDP, trade
  openness, bank credit depth, and a hyperinflation flag — all measured
  at $t-1$. Global financial conditions (fed funds, VIX, US 10-year yield)
  are excluded from the outcome equation (absorbed by year FE) and enter
  only the first-stage propensity models.
- **Driscoll-Kraay SEs** throughout the single-stage local projections
  (`xtscc`, lag = $\max(1,h+1)$); country-fixed-effects-only outcome
  models with cluster-robust SEs for the two-stage IPW/AIPW estimators,
  matching the reference paper's design for that stage specifically.

---

## What changed this session (for anyone picking this up)

This session fixed a horizon-indexing bug affecting several files after
the Year-1-relabeling pass (matrices sized for the old indexing, causing
a Stata conformability error in `08_ipw_lp.do`'s Act-2 p-value matrices),
added a pre-trend-controlled robustness specification to
`03_lp_resolution.do`, built the Gelbach decomposition
(`12b_gelbach_decomposition.do`), and rewrote `EMPIRICAL_ANALYSIS.md` and
this README against verified current-code runs — correcting, among other
things, an old claim that IPW reweighting *increases* the estimated
output cost (the current, verified run shows the opposite: it
attenuates it) and an old "credit and investment are a purely
non-default phenomenon" channel narrative that the current, corrected
data no longer supports. See `EMPIRICAL_ANALYSIS.md`'s limitations
section and `RESULTS_SECTION_DRAFT.md` for the fully qualified current
picture.
