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

* ── Run pipeline ────────────────────────────────────────────────────────────
do "$do/01_build_panel.do"
do "$do/01b_merge_new_controls.do"
do "$do/01c_merge_nexus.do"
do "$do/01d_merge_vulnerability.do"   // IDS: rollover, reserves, interest burden
do "$do/01e_predictors.do"            // first-stage exclusion-restriction predictors (Z2, Z3)
do "$do/02_lp_all.do"
do "$do/03_lp_resolution.do"
do "$do/04_graphs.do"
do "$do/05_balance_table.do"
do "$do/06_robustness.do"
do "$do/07_placebo.do"
do "$do/08_ipw_lp.do"
* do "$do/09_lp_imf.do"   // removed: IMF selection unpredictable from observables
* do "$do/10_heterogeneity.do"  // removed: frontier variable poorly coded, duration data incomplete
do "$do/11_channels.do"
do "$do/11b_nexus_channels.do"
do "$do/12_channels_resolution.do"
do "$do/13_mechanisms.do"
do "$do/13b_exposure_heterogeneity.do"   // Tier-3: onset x pre-crisis channel exposure

* ── Structural model: calibration + nonlinear default block + transmission ──
do "$do/14_calibration.do"      // calibrate params (literature + data moments)
do "$do/15_solve_default.do"    // Arellano-style VFI: endogenous default & spread
do "$do/16_model_irf.do"        // log-linear transmission; model vs. data IRFs
