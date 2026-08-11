/*===========================================================================
  MASTER DO-FILE
  Sovereign Spread Crises — Output Cost Estimation
  Local Projections (Jordà 2005)

  Structure:
    01_build_panel.do   → import, merge classification, save panel_lp.dta
    02_lp_all.do        → Act 1: LP for ALL 61 spread crisis episodes
    03_lp_resolution.do → Act 2: LP split by default-linked vs. non-default
    04_graphs.do        → publication-quality IRF figures
    (09_lp_imf.do removed: IMF selection model too weak for credible inference)

  Required packages (run once):
    ssc install xtscc
    ssc install coefplot
    ssc install estout   (provides esttab/eststo — used for result tables)
    ssc install boottest
===========================================================================*/

clear all
set more off
set scheme s2mono

* ── Paths ──────────────────────────────────────────────────────────────────
global root  "/home/user/Research-project-sovereign-spread-crises-"
global raw   "$root/data/raw"
global clean "$root/data/clean"
global do    "$root/do"
global figs  "$root/output/figures"
global tabs  "$root/output/tables"

* ── Asonuma-aligned common-core control set (outcome-model controls) ─────────
* Used uniformly across all GDP + channel regressions. Also (re)defined in
* 18_transforms.do; set here so analysis files run standalone after a build.
* The plain lagged columns (l_govexp/l_open/l_credit_bank/l_hyperinfl) are
* created in 18_transforms.do and saved in panel_lp.dta.
global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"

* ══════════════════════════════════════════════════════════════════════════
* DATA BUILD (10-18) — from-scratch, fully-sourced panel. Keeps ONLY the own
* spread-crisis DB for the crisis definition; every macro series is rebuilt from
* official sources (IMF WEO / World Bank WDI+IDS / IMF MFS / FRED / Laeven-
* Valencia) and every transform is done in code. Replaces the old 01x chain
* (01_build_panel .. 01e_predictors), which is retained in do/ but no longer run.
* Provenance documented in DATA_SOURCES.md.
* ══════════════════════════════════════════════════════════════════════════
do "$do/10_skeleton.do"      // crisis DB -> skeleton (iso3, onsets, spreads, class., imf legacy)
do "$do/11_weo.do"           // IMF WEO: GDP/cap, GDP, pop, infl, ca, debt, govexp, revenue, inv, pb
do "$do/12_wdi.do"           // World Bank WDI: credit(by banks), fdi, claims_govt, trade openness, REER
do "$do/13_ifs_nexus.do"     // IMF MFS: sovereign-bank nexus (claims on govt/private / assets)
do "$do/14_ids.do"           // World Bank IDS: stdebt_share, reserves, intpay, debt_service
do "$do/15_rates.do"         // FRED: ust10y, fedfunds, vix (global push predictors)
do "$do/16_banking.do"       // Laeven-Valencia banking-crisis dummy
do "$do/17_predictors.do"    // derived Z: contagion + proneness; numeric cid + xtset
do "$do/18_transforms.do"    // gdpg, lags, dy_h, pre-trends, sample -> panel_lp.dta

* ══════════════════════════════════════════════════════════════════════════
* ANALYSIS (unchanged — all read $clean/panel_lp.dta produced by stage 18)
* ══════════════════════════════════════════════════════════════════════════
do "$do/02_lp_all.do"
do "$do/03_lp_resolution.do"
do "$do/04_graphs.do"
do "$do/05_balance_table.do"
do "$do/06_robustness.do"
do "$do/07_placebo.do"
do "$do/08_ipw_lp.do"
do "$do/08b_aipw.do"                  // IPWRA (the reference paper's headline estimator) + stratified bootstrap CIs
do "$do/08c_first_stage_table.do"     // Table 1-style probit first stage (predictors/controls + chi2/ROC)
* do "$do/09_lp_imf.do"   // removed: IMF selection unpredictable from observables
* do "$do/10_heterogeneity.do"  // removed: frontier variable poorly coded, duration data incomplete
do "$do/11_channels.do"
do "$do/11b_nexus_channels.do"
do "$do/12_channels_resolution.do"
do "$do/13_mechanisms.do"
do "$do/13b_exposure_heterogeneity.do"   // Tier-3: onset x pre-crisis channel exposure
do "$do/13c_aipw_channels.do"            // doubly-robust AIPW for the transmission channels (matches 08b)
do "$do/13d_aipw_nexus_split.do"         // bank-intermediation heterogeneity: AIPW by sovereign-bank nexus (Asonuma Fig 6)
do "$do/13e_nexus_mechanism_diagram.do"  // conceptual 2x2 schematic of the nexus x resolution mechanism

* ── Structural model: calibration + nonlinear default block + transmission ──
do "$do/14_calibration.do"      // calibrate params (literature + data moments)
do "$do/15_solve_default.do"    // Arellano-style VFI: endogenous default & spread
do "$do/16_model_irf.do"        // log-linear transmission; model vs. data IRFs
