/*===========================================================================
  08_IPW_LP.DO
  Inverse Probability Weighted Local Projections (IPW-LP)
  Following Jordà & Taylor (2016) "The Time for Austerity"

  Steps:
    1. First-stage probit: Pr(onset_all=1 | X_{t-1})
       Controls X: l1_gdpg, l2_gdpg, debt, ca, infl, imf
       Predictors Z: fedfunds (single global push), l_reg_crisis_share, past_onsets
       (fedfunds included in probit only; absorbed by year FE in the LP)
    2. Predict propensity scores p_it
    3. Trim extreme scores (< 0.01 or > 0.99) to avoid explosive weights
    4. Construct stabilized IPW weights:
         treated:  w = Pr(D=1) / p_it
         control:  w = Pr(D=0) / (1 - p_it)
    5. Run LP at h=0..4 with IPW weights (xtreg fe, cluster cid)
    6. Compare IPW-LP vs. baseline LP:
         - If close   → selection on observables not driving baseline
         - If diverge → baseline was biased; IPW estimate preferred

  Saves:
    "$clean/irf_ipw.dta"        — IPW-LP IRF dataset
    "$clean/irf_compare.dta"    — baseline + IPW for overlay plot
    "$tabs/ipw_first_stage.csv" — first-stage probit results
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking l_govexp l_open l_credit hyperinf_dummy"

* LP controls: VIX and ust10y have zero cross-sectional variation and are
* fully absorbed by year fixed effects (i.year) in all LP regressions.
* They are excluded from the LP spec but retained in the probit first stage,
* which has no year FE and where global financial conditions are not absorbed.
local controls_lp  $ctrl_core

* First-stage probit = controls X (same fundamentals as the LP) + excluded
* PREDICTORS Z (Jordà–Taylor 2016). The predictors satisfy exclusion because
* they are omitted from the LP second stage:
*   fedfunds        — pure time-series push factor, absorbed by i.year in the LP
*   l_reg_crisis_share (Z2) — regional contagion, country-varying, omitted from LP
*   past_onsets       (Z3) — own crisis history, country-varying, omitted from LP
* Strict parity with the reference paper: the probit BASELINE = the outcome
* baseline ($ctrl_core), so the first stage and the LP share the same controls;
* the predictors Z below are the only first-stage-specific additions (their
* $convar in both stages, $instrument added only to the probit).
local controls_x   $ctrl_core
* Single global-push predictor (fed funds): fedfunds/ust10y/vix all proxy the
* same global-financial-conditions factor, so only one is kept to avoid splitting
* its explanatory power. fedfunds is pure time-series => absorbed by year FE in
* the LP, so the exclusion restriction still holds.
local predictors_z fedfunds l_reg_crisis_share past_onsets

* ══════════════════════════════════════════════════════════════════════════
* STEP 1 — FIRST-STAGE PROBIT
* Dependent variable: onset_all (=1 in first year of any crisis episode)
* Estimated on estimation sample only (continuation years excluded)
* All RHS variables lagged to ensure predetermination
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== FIRST STAGE: Probit of crisis onset ==="

* POOLED probit (no country FE): country FE in a probit separate/perfectly
* predict never-treated countries and carry incidental-parameters bias with few
* events per group; the pooled model with controls + predictors is the standard
* for rare-event propensity estimation and keeps all countries.

* (a) Controls only — baseline for the ROC comparison (their Table 1 logic)
quietly probit onset_all `controls_x' if sample == 1, vce(cluster cid)
quietly lroc, nograph
local roc_x = r(area)

* (b) Controls + excluded predictors — the first stage used for the propensity
probit onset_all `controls_x' `predictors_z' if sample == 1, vce(cluster cid)
estimates store first_stage
quietly lroc, nograph
local roc_xz = r(area)

* First-stage diagnostics
estat classification
di as result "McFadden Pseudo-R2: " e(r2_p)
di as result "Area under ROC: controls only = " %5.3f `roc_x' ///
             "   +predictors = " %5.3f `roc_xz' ///
             "   (gain = " %5.3f (`roc_xz'-`roc_x') ")"

* Export first-stage table
quietly esttab first_stage using "$tabs/ipw_first_stage.csv", ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    title("First-Stage Probit: Crisis Onset") ///
    replace

* ══════════════════════════════════════════════════════════════════════════
* STEP 2 — PREDICT PROPENSITY SCORES
* ══════════════════════════════════════════════════════════════════════════

predict pscore if sample == 1, pr
label var pscore "Propensity score: Pr(onset=1 | X)"

* Inspect distribution of scores
summarize pscore, detail
histogram pscore if sample == 1, ///
    xtitle("Propensity score") title("Distribution of Propensity Scores") ///
    note("Estimation sample only.") ///
    graphregion(color(white)) fcolor("23 55 94%40") lcolor("23 55 94")
graph export "$figs/fig6a_pscore_dist.pdf", replace

* Common support check: overlap between treated and control
twoway ///
    (histogram pscore if onset_all==0 & sample==1, ///
        fcolor("150 150 150%50") lcolor(none) width(0.02))  ///
    (histogram pscore if onset_all==1 & sample==1, ///
        fcolor("23 55 94%60")    lcolor(none) width(0.02)), ///
    xtitle("Propensity score", size(small)) ///
    title("Common Support Check", size(medium)) ///
    subtitle("Control (grey) vs. Treated (blue)", size(small)) ///
    legend(off) graphregion(color(white)) plotregion(color(white))
graph export "$figs/fig6b_common_support.pdf", replace

* ── Kernel-density overlap (Asonuma et al. 2024 Fig 2 analog) ─────────────
* Cleaner common-support view: densities of the estimated propensity score for
* treated (onset) vs control (tranquil). Overlapping densities => the design has
* common support, so IPW/AIPW reweighting is valid.
capture twoway ///
    (kdensity pscore if onset_all==1 & sample==1, ///
        lcolor("157 36 73") lwidth(medthick)) ///
    (kdensity pscore if onset_all==0 & sample==1, ///
        lcolor("23 55 94") lwidth(medthick) lpattern(dash)), ///
    xtitle("Estimated propensity score  Pr(onset=1 | X, Z)", size(small)) ///
    ytitle("Density", size(small)) ///
    title("Act 1 — Propensity-score overlap (common support)", size(medsmall) color(navy)) ///
    legend(order(1 "Onset (treated)" 2 "Tranquil (control)") ring(0) pos(1) size(small)) ///
    note("Overlapping densities => common support holds. First stage: probit of onset on controls + excluded predictors.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))
if _rc == 0 {
    graph export "$figs/fig6c_overlap_act1.pdf", replace
    di as result "Figure saved: fig6c_overlap_act1.pdf"
}
else di as error "  ** Act 1 overlap kdensity failed (rc=" _rc ")"

* ══════════════════════════════════════════════════════════════════════════
* STEP 3 — TRIM EXTREME SCORES
* Scores very close to 0 or 1 produce explosive weights.
* Standard trim: [0.01, 0.99]. Adjust if needed based on distribution.
* ══════════════════════════════════════════════════════════════════════════

local trim_lo = 0.01
local trim_hi = 0.99

gen trimmed = (pscore < `trim_lo' | pscore > `trim_hi') if !missing(pscore)
quietly count if trimmed == 1
di as result "Observations trimmed (score outside [`trim_lo', `trim_hi']): " r(N)

replace pscore = . if trimmed == 1

* ══════════════════════════════════════════════════════════════════════════
* STEP 4 — CONSTRUCT STABILIZED IPW WEIGHTS
* Stabilized weights use marginal probability of treatment Pr(D=1)
* instead of 1, which reduces variance of the weighted estimator.
* ══════════════════════════════════════════════════════════════════════════

quietly summarize onset_all if sample == 1
local p_treat = r(mean)   // marginal Pr(D=1) in estimation sample
local p_ctrl  = 1 - `p_treat'

di as result "Marginal Pr(onset=1) in estimation sample: " %5.4f `p_treat'

gen ipw = .
replace ipw = `p_treat' / pscore         if onset_all == 1 & !missing(pscore)
replace ipw = `p_ctrl'  / (1 - pscore)   if onset_all == 0 & !missing(pscore)
label var ipw "Stabilized IPW weight"

* Inspect weight distribution
summarize ipw if sample==1, detail
di "Mean weight (treated): "
summarize ipw if onset_all==1 & sample==1

* ══════════════════════════════════════════════════════════════════════════
* STEP 5 — IPW-WEIGHTED LOCAL PROJECTIONS
* xtreg fe with pweights and cluster SE
* Note: xtscc does not accept pweights; xtreg fe vce(cluster cid) used.
* For consistency, also re-run baseline with xtreg (not xtscc) so the
* comparison is apples-to-apples in terms of SE method.
* ══════════════════════════════════════════════════════════════════════════

foreach m in b_ipw lo90_ipw hi90_ipw lo95_ipw hi95_ipw ///
             b_ols lo90_ols hi90_ols lo95_ols hi95_ols {
    matrix `m' = J(5, 1, .)
}

di as result _n "=== IPW-LP vs. UNWEIGHTED LP (xtreg fe, cluster cid) ==="
di "h    beta_OLS   SE_OLS   beta_IPW   SE_IPW   delta"

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    local row = `h' + 1

    * Unweighted (areg absorbs country FE, cluster SE — comparable baseline)
    quietly areg dy_`h' onset_all `controls_lp' i.year ///
        if sample==1 & !missing(ipw), absorb(cid) vce(cluster cid)
    matrix b_ols[`row',1]    = _b[onset_all]
    matrix lo90_ols[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
    matrix hi90_ols[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
    matrix lo95_ols[`row',1] = _b[onset_all] - 1.960*_se[onset_all]
    matrix hi95_ols[`row',1] = _b[onset_all] + 1.960*_se[onset_all]
    local b_u = _b[onset_all]
    local se_u = _se[onset_all]

    * IPW-weighted (areg used because xtreg fe requires weights constant within cid)
    quietly areg dy_`h' onset_all `controls_lp' i.year ///
        [aw=ipw] if sample==1 & !missing(ipw), absorb(cid) vce(cluster cid)
    matrix b_ipw[`row',1]    = _b[onset_all]
    matrix lo90_ipw[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
    matrix hi90_ipw[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
    matrix lo95_ipw[`row',1] = _b[onset_all] - 1.960*_se[onset_all]
    matrix hi95_ipw[`row',1] = _b[onset_all] + 1.960*_se[onset_all]
    local b_w = _b[onset_all]
    local se_w = _se[onset_all]

    di "h=" `h' "   " %7.3f `b_u' "   " %6.3f `se_u' ///
           "   " %7.3f `b_w' "   " %6.3f `se_w' ///
           "   " %7.3f (`b_w' - `b_u')
}

* ══════════════════════════════════════════════════════════════════════════
* STEP 6 — BUILD IRF DATASETS AND OVERLAY FIGURE
* ══════════════════════════════════════════════════════════════════════════

* IPW IRF dataset
clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_ipw, names(`m')
    rename `m'1 `m'
}
gen series = "ipw"
save "$clean/irf_ipw.dta", replace

* Unweighted (OLS-FE) IRF dataset
clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_ols, names(`m')
    rename `m'1 `m'
}
gen series = "ols"
save "$clean/irf_ols_fe.dta", replace

* Combined dataset for comparison plot
use "$clean/irf_ipw.dta", clear
append using "$clean/irf_ols_fe.dta"
save "$clean/irf_compare.dta", replace

* ── Figure 7: Baseline vs. IPW-LP overlay ────────────────────────────────
local c_ols "23 55 94"
local c_ipw "157 36 73"

use "$clean/irf_compare.dta", clear

twoway ///
    (rarea lo90 hi90 horizon if series=="ols",                              ///
        color("`c_ols'%20") lwidth(none))                                   ///
    (connected b horizon if series=="ols",                                  ///
        lcolor("`c_ols'") lwidth(medthick) msymbol(circle) mcolor("`c_ols'")) ///
    (rarea lo90 hi90 horizon if series=="ipw",                              ///
        color("`c_ipw'%20") lwidth(none))                                   ///
    (connected b horizon if series=="ipw",                                  ///
        lcolor("`c_ipw'") lwidth(medthick) lpattern(dash)                   ///
        msymbol(square) mcolor("`c_ipw'")),                                 ///
    yline(0, lpattern(dash) lcolor(gray) lwidth(thin))                      ///
    xlabel(0(1)4, labsize(medsmall))                                        ///
    ylabel(, format(%4.1f) labsize(medsmall))                               ///
    xtitle("Years after crisis onset", size(medsmall))                      ///
    ytitle("Cumulative change in log real GDP (pp)", size(medsmall))   ///
    title("Baseline vs. IPW-Weighted Local Projections", size(medium))      ///
    subtitle("All spread crises (N=61). xtreg FE, cluster SE.", size(small)) ///
    legend(order(2 "Baseline OLS-FE" 4 "IPW-weighted")                     ///
           ring(0) pos(3) size(small))                                      ///
    note("IPW weights from probit of onset on lagged macro fundamentals."   ///
         "Stabilized weights, trimmed at [0.01, 0.99].", size(vsmall))      ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig7_ipw_vs_baseline.pdf", replace
graph export "$figs/fig7_ipw_vs_baseline.png", replace width(1200)
di as result _n "Figure 7 saved: $figs/fig7_ipw_vs_baseline.pdf"
di as result "All IPW results saved."

* ══════════════════════════════════════════════════════════════════════════
* ACT 2 IPW — RESOLUTION SPLIT
* Question: conditional on experiencing a spread crisis, does being
*   default-linked (vs. non-default) cause extra output loss?
*
* Treatment: onset_def = 1 (default-linked episode)
* Control:   onset_nd  = 1 (non-default episode)
* Sample:    crisis onset years only (onset_all == 1)
*
* The balance table showed default-linked countries had worse pre-crisis
* fundamentals (lower growth, higher debt, higher spreads). IPW reweights
* non-default episodes to look like default-linked ones on observables,
* isolating the effect of resolution type from pre-existing differences.
* ══════════════════════════════════════════════════════════════════════════

use "$clean/panel_lp.dta", clear

* Among-onsets default propensity — LEAN, WELL-COVERED baseline (NOT the full
* $ctrl_core). This first stage predicts default-vs-non-default among the ~51
* crisis onsets only; the full 8-var $ctrl_core + 3 predictors (11 params) on that
* thin cell PERFECTLY SEPARATES once its coverage-hole lags (l_credit/l_banking/
* l_open) shrink the cell to ~23 obs (AUROC=1, pscore in {0,1} => IPW dies). We
* therefore use a compact, well-covered predetermined baseline and a further lean
* fallback. Asonuma never runs an among-onsets first stage; this cell is ours.
local controls    l1_gdpg l_debt l_ca
local predictors_z2 fedfunds l_reg_crisis_share past_def_onsets

di as result _n "=== ACT 2 IPW: FIRST STAGE — Probit of default-linked onset ==="
di as result    "    Sample: crisis onset years only (onset_all == 1)"
di as result    "    Lean well-covered baseline (l1_gdpg l_debt l_ca) + predictors Z2"
di as result    "    (fed funds, regional contagion, past DEFAULT onsets)."

quietly probit onset_def `controls' if onset_all == 1, vce(robust)
quietly lroc, nograph
local roc2_x = r(area)

capture probit onset_def `controls' `predictors_z2' if onset_all == 1, vce(robust)
if _rc {
    di as error "  ** Act 2 probit (lean X + Z2) errored (rc=" _rc "); minimal fallback."
    probit onset_def l_debt l_ca `predictors_z2' if onset_all == 1, vce(robust)
}
quietly lroc, nograph
local roc2_xz = r(area)

estimates store fs_act2

di as result "Area under ROC (Act 2): controls only = " %5.3f `roc2_x' ///
             "   +predictors Z2 = " %5.3f `roc2_xz'
di as result "McFadden Pseudo-R2: " e(r2_p)

quietly esttab fs_act2 using "$tabs/ipw_act2_first_stage.csv", ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    title("Act 2 First Stage: Pr(Default-Linked | Crisis Onset) - Parsimonious") replace

* Predict propensity scores (among crisis onsets only)
predict pscore2 if onset_all == 1, pr
label var pscore2 "Pr(default-linked | crisis onset, X)"

* ── Separation guard (the `if _rc' above does NOT catch perfect prediction:
* Stata's probit drops perfectly-predicted obs and finishes with rc=0). If the
* propensity has (almost) no interior mass among onsets, it separated — refit the
* minimal Z2-based spec so pscore2 is usable and the IPW weights survive.
quietly count if inrange(pscore2, 0.02, 0.98) & onset_all == 1
if r(N) < 15 {
    di as error "  ** Act 2 propensity separated (only `r(N)' interior scores) — refitting minimal spec (l_debt l_ca + Z2)."
    capture drop pscore2
    probit onset_def l_debt l_ca `predictors_z2' if onset_all == 1, vce(robust)
    quietly lroc, nograph
    local roc2_xz = r(area)
    estimates store fs_act2
    predict pscore2 if onset_all == 1, pr
    label var pscore2 "Pr(default-linked | crisis onset, X)"
    di as result "  Refit AUROC (Act 2, minimal): " %5.3f `roc2_xz'
}

summarize pscore2, detail
di as result "Range of pscore2: [" r(min) ", " r(max) "]"

* ── Kernel-density overlap (Fig 2 analog) — THE decisive check for the thin
* default sample: densities of Pr(default | crisis) for default vs non-default
* onsets. Plotted BEFORE trimming to show the raw support. Overlap => the
* default propensity is estimable off comparable non-default onsets; little
* overlap => the resolution IPW rests on extrapolation (read with caution). ──
capture twoway ///
    (kdensity pscore2 if onset_def==1, lcolor("157 36 73") lwidth(medthick)) ///
    (kdensity pscore2 if onset_nd==1,  lcolor("34 139 34") lwidth(medthick) lpattern(dash)), ///
    xtitle("Estimated propensity score  Pr(default-linked | crisis, X, Z)", size(small)) ///
    ytitle("Density", size(small)) ///
    title("Act 2 — Propensity overlap: default vs non-default onsets", size(medsmall) color(navy)) ///
    subtitle("Common-support check for the resolution IPW (small default cell)", size(small)) ///
    legend(order(1 "Default-linked (treated)" 2 "Non-default (control)") ring(0) pos(1) size(small)) ///
    note("Overlap => default propensity estimable off comparable non-default onsets. Thin default cell (~21 events): read with the ROC.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))
if _rc == 0 {
    graph export "$figs/fig8b_overlap_act2.pdf", replace
    di as result "Figure saved: fig8b_overlap_act2.pdf"
}
else di as error "  ** Act 2 overlap kdensity failed (rc=" _rc ") — likely too few default obs"

* Trim extreme scores
gen trimmed2 = (pscore2 < 0.05 | pscore2 > 0.95) if !missing(pscore2)
quietly count if trimmed2 == 1
di as result "Observations trimmed (score outside [0.05, 0.95]): " r(N)
replace pscore2 = . if trimmed2 == 1

* Stabilized IPW weights (among crisis episodes)
quietly summarize onset_def if onset_all == 1
local p_def = r(mean)
local p_nd  = 1 - `p_def'
di as result "Share of default-linked among crisis onsets: " %5.4f `p_def'

gen ipw2 = .
replace ipw2 = `p_def' / pscore2         if onset_def == 1 & !missing(pscore2)
replace ipw2 = `p_nd'  / (1 - pscore2)   if onset_nd  == 1 & !missing(pscore2)
label var ipw2 "Act 2 stabilized IPW weight"

summarize ipw2 if onset_all == 1, detail

* ── IPW-weighted LP: non-default vs. default-linked ──────────────────────
* Regression on crisis onset sample only, treatment = onset_def
* Controls absorb remaining differences; IPW reweights for pre-crisis chars

foreach m in b_def_ipw lo90_def_ipw hi90_def_ipw b_nd_ipw lo90_nd_ipw hi90_nd_ipw ///
             b_def_ols lo90_def_ols hi90_def_ols b_nd_ols lo90_nd_ols hi90_nd_ols {
    matrix `m' = J(5, 1, .)
}

di as result _n "=== ACT 2 IPW RESULTS ==="
di "h   b_nd_OLS  b_def_OLS  p_OLS   b_nd_IPW  b_def_IPW  p_IPW   delta_nd  delta_def"

* Store equality-test p-values: H0 beta_nd(h) = beta_def(h), OLS and IPW
matrix pval_act2_ols = J(5, 1, .)
matrix pval_act2_ipw = J(5, 1, .)

* Restore full panel for LP (need all country-years, not just onsets)
* IPW2 weight is non-missing only for onset observations — apply via if condition

forvalues h = 0/4 {
    local row = `h' + 1

    * Unweighted (baseline Act 2, OLS-FE, cluster SE)
    quietly areg dy_`h' onset_nd onset_def `controls_lp' i.year ///
        if sample == 1, absorb(cid) vce(cluster cid)
    matrix b_nd_ols[`row',1]   = _b[onset_nd]
    matrix lo90_nd_ols[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
    matrix hi90_nd_ols[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
    matrix b_def_ols[`row',1]  = _b[onset_def]
    matrix lo90_def_ols[`row',1] = _b[onset_def] - 1.645*_se[onset_def]
    matrix hi90_def_ols[`row',1] = _b[onset_def] + 1.645*_se[onset_def]
    local b_nd_u  = _b[onset_nd]
    local b_def_u = _b[onset_def]
    quietly test onset_nd = onset_def
    matrix pval_act2_ols[`row',1] = r(p)
    local p_u = r(p)

    * IPW-weighted: upweight onset observations by ipw2
    * Control observations (onset_all=0) keep weight = 1 (set ipw2=1 for them)
    quietly replace ipw2 = 1 if onset_all == 0 & sample == 1
    quietly areg dy_`h' onset_nd onset_def `controls_lp' i.year ///
        [aw=ipw2] if sample == 1 & !missing(ipw2), absorb(cid) vce(cluster cid)
    quietly test onset_nd = onset_def
    matrix pval_act2_ipw[`row',1] = r(p)
    local p_w = r(p)
    quietly replace ipw2 = . if onset_all == 0

    matrix b_nd_ipw[`row',1]    = _b[onset_nd]
    matrix lo90_nd_ipw[`row',1]  = _b[onset_nd]  - 1.645*_se[onset_nd]
    matrix hi90_nd_ipw[`row',1]  = _b[onset_nd]  + 1.645*_se[onset_nd]
    matrix b_def_ipw[`row',1]   = _b[onset_def]
    matrix lo90_def_ipw[`row',1] = _b[onset_def] - 1.645*_se[onset_def]
    matrix hi90_def_ipw[`row',1] = _b[onset_def] + 1.645*_se[onset_def]
    local b_nd_w  = _b[onset_nd]
    local b_def_w = _b[onset_def]

    di "h=" `h' "  " %7.3f `b_nd_u'  "    " %7.3f `b_def_u' "   " %5.3f `p_u' ///
           "    " %7.3f `b_nd_w'  "    " %7.3f `b_def_w' "   " %5.3f `p_w' ///
           "    " %7.3f (`b_nd_w'  - `b_nd_u') ///
           "    " %7.3f (`b_def_w' - `b_def_u')
}

* ── Save Act 2 IPW results ────────────────────────────────────────────────
clear
set obs 5
gen horizon = _n - 1

foreach m in b_nd b_def lo90_nd hi90_nd lo90_def hi90_def {
    svmat `m'_ipw, names(`m')
    rename `m'1 `m'
}
gen series = "ipw"
save "$clean/irf_act2_ipw.dta", replace

clear
set obs 5
gen horizon = _n - 1
foreach m in b_nd b_def lo90_nd hi90_nd lo90_def hi90_def {
    svmat `m'_ols, names(`m')
    rename `m'1 `m'
}
gen series = "ols"
save "$clean/irf_act2_ols.dta", replace

* ── Figure 8: Act 2 OLS vs IPW overlay ───────────────────────────────────
use "$clean/irf_act2_ipw.dta", clear
append using "$clean/irf_act2_ols.dta"

local c_nd  "0 84 166"
local c_def "157 36 73"

twoway ///
    (connected b_nd horizon if series=="ols", ///
        lcolor("`c_nd'") lwidth(medthick) msymbol(circle) mcolor("`c_nd'")) ///
    (connected b_def horizon if series=="ols", ///
        lcolor("`c_def'") lwidth(medthick) msymbol(square) mcolor("`c_def'")) ///
    (connected b_nd horizon if series=="ipw", ///
        lcolor("`c_nd'") lwidth(medthick) lpattern(dash) ///
        msymbol(circle_hollow) mcolor("`c_nd'")) ///
    (connected b_def horizon if series=="ipw", ///
        lcolor("`c_def'") lwidth(medthick) lpattern(dash) ///
        msymbol(square_hollow) mcolor("`c_def'")), ///
    yline(0, lpattern(dash) lcolor(gray) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years after crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP (pp)", size(medsmall)) ///
    title("Resolution Split: OLS vs. IPW-Weighted LP", size(medium)) ///
    subtitle("Solid = OLS baseline. Dashed = IPW-weighted.", size(small)) ///
    legend(order(1 "Non-default (OLS)" 2 "Default-linked (OLS)" ///
                 3 "Non-default (IPW)" 4 "Default-linked (IPW)") ///
           cols(2) size(small)) ///
    note("IPW reweights non-default episodes to match default-linked" ///
         "on pre-crisis fundamentals (probit first stage).", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig8_act2_ipw.pdf", replace
graph export "$figs/fig8_act2_ipw.png", replace width(1200)
di as result _n "Figure 8 saved: $figs/fig8_act2_ipw.pdf"
di as result "Act 2 IPW complete."
