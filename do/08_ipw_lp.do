/*===========================================================================
  08_IPW_LP.DO
  Inverse Probability Weighted Local Projections (IPW-LP)
  Following Jordà & Taylor (2016) "The Time for Austerity"

  Steps:
    1. First-stage probit: Pr(onset_all=1 | X_{t-1})  [pooled, no country FE]
       Controls X: the common core ($ctrl_core), all predetermined
       Predictors Z: l_fedfunds (single global push, t-1), l_reg_crisis_share, past_onsets
       (predictors enter the probit only; omitted from the country-FE-only LP)
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
    "$figs/fig6c_overlap_act1.pdf"  — propensity overlap, all onsets vs tranquil
    "$figs/fig8b_overlap_act2.pdf"  — propensity overlap BY resolution type
                                      (Asonuma et al. 2024 Fig. 2 analog)
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"

* LP controls = the common core. Paper-aligned: the two-stage outcome regressions
* carry COUNTRY FE only (no year FE), so global-push predictors are excluded from
* the outcome equation simply by being OMITTED from it (not by year-FE absorption).
local controls_lp  $ctrl_core

* First-stage probit = controls X (same fundamentals as the LP) + excluded
* PREDICTORS Z (Jordà–Taylor 2016). The predictors satisfy exclusion because
* they are omitted from the LP second stage:
*   fedfunds        — pure time-series push factor, omitted from the LP outcome eq.
*   l_reg_crisis_share (Z2) — regional contagion, country-varying, omitted from LP
*   past_onsets       (Z3) — own crisis history, country-varying, omitted from LP
* Strict parity with the reference paper: the probit BASELINE = the outcome
* baseline ($ctrl_core), so the first stage and the LP share the same controls;
* the predictors Z below are the only first-stage-specific additions (their
* $convar in both stages, $instrument added only to the probit).
local controls_x   $ctrl_core
* Single global-push predictor, LAGGED to t-1 (l_fedfunds = L.fedfunds), matching
* the reference paper's federal_funds2 = L.federal_funds. fedfunds/ust10y/vix all
* proxy the same global-financial-conditions factor, so only one is kept to avoid
* splitting its explanatory power. Predetermined and excluded from the outcome eq.
local predictors_z l_fedfunds l_reg_crisis_share past_onsets

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
* for rare-event propensity estimation and keeps all countries. Country FE enter
* the OUTCOME regression (Eq. 1 / the LP) only, matching the paper; no year FE.

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
* Colours follow the paper's Fig 2 convention — treatment solid blue, control dashed
* red — so this and the Act-2 panels (fig8b) read as one family.
capture twoway ///
    (kdensity pscore if onset_all==1 & sample==1, ///
        range(0 1) lcolor("23 55 94") lwidth(medthick)) ///
    (kdensity pscore if onset_all==0 & sample==1, ///
        range(0 1) lcolor("157 36 73") lwidth(medthick) lpattern(dash)), ///
    xlabel(0(.2)1, labsize(small)) ///
    xtitle("Predicted probability  Pr(onset=1 | X, Z)", size(small)) ///
    ytitle("Probability density", size(small)) ///
    title("Act 1 — Propensity-score overlap (common support)", size(medsmall) color(navy)) ///
    legend(order(1 "Treatment group (onset)" 2 "Control group (tranquil)") ring(0) pos(1) size(small)) ///
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

* Row 1 = explicit baseline (displayed h=0, hardcoded 0); rows 2-6 = displayed h=1..5
foreach m in b_ipw lo90_ipw hi90_ipw lo95_ipw hi95_ipw ///
             b_ols lo90_ols hi90_ols lo95_ols hi95_ols {
    matrix `m' = J(6, 1, 0)
}

di as result _n "=== IPW-LP vs. UNWEIGHTED LP (xtreg fe, cluster cid) ==="
di "h    beta_OLS   SE_OLS   beta_IPW   SE_IPW   delta"

forvalues h = 0/4 {
    local hd  = `h' + 1
    local lag = max(1, `h'+1)
    local row = `h' + 2

    * Unweighted (areg absorbs country FE, cluster SE — comparable baseline)
    * Country FE only, no year FE (paper-aligned two-stage outcome regression).
    quietly areg dy_`h' onset_all `controls_lp' ///
        if sample==1 & !missing(ipw), absorb(cid) vce(cluster cid)
    matrix b_ols[`row',1]    = _b[onset_all]
    matrix lo90_ols[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
    matrix hi90_ols[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
    matrix lo95_ols[`row',1] = _b[onset_all] - 1.960*_se[onset_all]
    matrix hi95_ols[`row',1] = _b[onset_all] + 1.960*_se[onset_all]
    local b_u = _b[onset_all]
    local se_u = _se[onset_all]

    * IPW-weighted (areg used because xtreg fe requires weights constant within cid)
    quietly areg dy_`h' onset_all `controls_lp' ///
        [aw=ipw] if sample==1 & !missing(ipw), absorb(cid) vce(cluster cid)
    matrix b_ipw[`row',1]    = _b[onset_all]
    matrix lo90_ipw[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
    matrix hi90_ipw[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
    matrix lo95_ipw[`row',1] = _b[onset_all] - 1.960*_se[onset_all]
    matrix hi95_ipw[`row',1] = _b[onset_all] + 1.960*_se[onset_all]
    local b_w = _b[onset_all]
    local se_w = _se[onset_all]

    di "h=" `hd' "   " %7.3f `b_u' "   " %6.3f `se_u' ///
           "   " %7.3f `b_w' "   " %6.3f `se_w' ///
           "   " %7.3f (`b_w' - `b_u')
}

* ══════════════════════════════════════════════════════════════════════════
* STEP 6 — BUILD IRF DATASETS AND OVERLAY FIGURE
* ══════════════════════════════════════════════════════════════════════════

* IPW IRF dataset
clear
set obs 6
gen horizon = _n - 1     // 0 (baseline), 1, 2, 3, 4, 5
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_ipw, names(`m')
    rename `m'1 `m'
}
gen series = "ipw"
save "$clean/irf_ipw.dta", replace

* Unweighted (OLS-FE) IRF dataset
clear
set obs 6
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
    xlabel(0(1)5, labsize(medsmall))                                        ///
    ylabel(, format(%4.1f) labsize(medsmall))                               ///
    xtitle("Year (Year 1 = crisis year)", size(medsmall))                   ///
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
* ACT 2 IPW — RESOLUTION SPLIT, per-type-vs-tranquil (reference-paper structure)
* Question: what is the output cost of each resolution type, and how much MORE
*   costly is a default-linked crisis than a non-default one?
*
* Design (Asonuma et al.): estimate EACH type vs TRANQUIL years, dropping the
* rival type from the sample so the control group is clean tranquil country-years
* (never the other crisis type):
*   Line 1: onset_nd  vs tranquil   (drop onset_def)
*   Line 2: onset_def vs tranquil   (drop onset_nd)
* IPW reweights tranquil years to match each type's pre-crisis fundamentals.
* Extra cost of default = (default line) - (non-default line), Clogg z. This
* mirrors the AIPW two-line design in 08b and avoids the thin among-onsets cell.
* ══════════════════════════════════════════════════════════════════════════

use "$clean/panel_lp.dta", clear

* PER-TYPE-vs-TRANQUIL (reference-paper structure). Each resolution type is scored
* vs TRANQUIL years with the RIVAL type dropped (Asonuma et al.'s sample_for*:
* control = tranquil, never the other type). Two first stages / two weight sets:
*   non-default vs tranquil  (drop onset_def) ;  default-linked vs tranquil (drop onset_nd)
* Full sample (~500 tranquil controls) => the full $ctrl_core is safe (no separation).
local predictors_z2 l_fedfunds l_reg_crisis_share past_def_onsets

di as result _n "=== ACT 2 IPW: two vs-tranquil first stages (rival type dropped) ==="

foreach s in nd def {
    if "`s'" == "nd"  local rival onset_def
    else              local rival onset_nd

    quietly probit onset_`s' $ctrl_core `predictors_z2' ///
        if sample==1 & `rival'==0, vce(cluster cid)
    quietly lroc, nograph
    di as result "  First stage `s' vs tranquil: AUROC = " %5.3f r(area) "   (N = " e(N) ")"

    capture drop pscore_`s'
    predict pscore_`s' if sample==1 & `rival'==0, pr

    * Untrimmed copy for the overlap density below: the Fig-2 analog exists to
    * DIAGNOSE common support, so it must show the raw score distribution —
    * trimming first would hide the very positivity problem the plot reveals.
    capture drop pscore_u_`s'
    quietly gen double pscore_u_`s' = pscore_`s'
    label var pscore_u_`s' "Predicted Pr(onset_`s'=1 | X,Z), untrimmed"

    quietly replace pscore_`s' = . if (pscore_`s'<0.01 | pscore_`s'>0.99) & !missing(pscore_`s')

    * stabilized weights: treated -> Pr(s)/p ; tranquil controls -> Pr(!s)/(1-p)
    quietly summarize onset_`s' if sample==1 & `rival'==0
    local pmarg = r(mean)
    capture drop ipw_`s'
    gen double ipw_`s' = .
    quietly replace ipw_`s' = `pmarg'     / pscore_`s'       if onset_`s'==1 & !missing(pscore_`s')
    quietly replace ipw_`s' = (1-`pmarg') / (1-pscore_`s')   if onset_all==0 & !missing(pscore_`s')
    label var ipw_`s' "Act 2 stabilized IPW weight (`s' vs tranquil)"
}

* ── Kernel-density overlap BY RESOLUTION TYPE (Asonuma et al. 2024 Fig 2 analog) ─
* Their Fig 2 gives ONE PANEL PER restructuring type, each plotting the density of
* the estimated probability for the TREATMENT group (that type's onsets) and for the
* CONTROL group (tranquil years). This is the Act-2 counterpart of fig6c (Act 1):
* our types are non-default / default-linked instead of their preemptive variants.
* Overlapping densities => common support holds for that type, so the vs-tranquil
* reweighting is valid. A control density crushed against 0 with a treated density
* far to the right would signal a positivity problem for that cell.
* Deviation from their Fig 2, deliberate: they plot only [0.01, 0.60], but we show the
* FULL unit interval. Truncating would hide mass in the upper tail — and a treated
* density piling up near 1 is precisely the positivity failure that has to stay visible
* here (it is what makes a country-FE first stage inestimable on this thin sample).
local ovnames
foreach s in nd def {
    if "`s'" == "nd" {
        local rival onset_def
        local slab  "Non-default"
    }
    else {
        local rival onset_nd
        local slab  "Default-linked"
    }

    * Treatment = this type's onsets; control = tranquil years (rival type dropped),
    * exactly the sample its first-stage probit and its IPW weights were built on.
    capture twoway ///
        (kdensity pscore_u_`s' if onset_`s'==1 & sample==1 & `rival'==0, ///
            range(0 1) lcolor("23 55 94")  lwidth(medthick)) ///
        (kdensity pscore_u_`s' if onset_all==0  & sample==1 & `rival'==0, ///
            range(0 1) lcolor("157 36 73") lwidth(medthick) lpattern(dash)), ///
        xlabel(0(.2)1, labsize(small)) ///
        xtitle("Predicted probability", size(small)) ///
        ytitle("Probability density", size(small)) ///
        title("`slab'", size(medsmall) color(navy)) ///
        legend(order(1 "Treatment group" 2 "Control group") ///
               ring(0) pos(1) size(vsmall) region(lstyle(none)) rows(2)) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ov_`s', replace)

    if _rc == 0 local ovnames `ovnames' ov_`s'
    else di as error "  ** Act 2 overlap kdensity failed for `s' (rc=" _rc ")"
}

if "`ovnames'" != "" {
    capture graph combine `ovnames', rows(1) ///
        graphregion(color(white)) ///
        title("Act 2 — Propensity-score overlap by resolution type", ///
              size(medsmall) color(navy)) ///
        note("Kernel densities of the estimated probability of each crisis type: treatment group (that type's onsets)" ///
             "vs control group (tranquil years), with the rival type dropped. Scores shown untrimmed over the full [0,1]." ///
             "Analog of Asonuma et al. (2024) Fig. 2 (they plot [0.01,0.60]); fig6c is the Act-1 counterpart.", size(vsmall))
    if _rc == 0 {
        graph export "$figs/fig8b_overlap_act2.pdf", replace
        capture graph export "$figs/fig8b_overlap_act2.png", replace width(1200)
        di as result "Figure saved: fig8b_overlap_act2.pdf"
    }
    else di as error "  ** Act 2 overlap combine failed (rc=" _rc ")"
    capture graph drop `ovnames'
}

* ── Estimation: OLS baseline (joint, tranquil = omitted) + two vs-tranquil IPW LPs
* Row 1 = explicit baseline (displayed h=0, hardcoded 0); rows 2-6 = displayed h=1..5
foreach m in b_def_ipw lo90_def_ipw hi90_def_ipw b_nd_ipw lo90_nd_ipw hi90_nd_ipw ///
             b_def_ols lo90_def_ols hi90_def_ols b_nd_ols lo90_nd_ols hi90_nd_ols {
    matrix `m' = J(6, 1, 0)
}
matrix pval_act2_ols = J(6, 1, .)
matrix pval_act2_ipw = J(6, 1, .)

di as result _n "=== ACT 2 RESULTS: each resolution type vs tranquil ==="
di "h    b_nd_OLS  b_def_OLS  p_OLS(Wald)    b_nd_IPW  b_def_IPW  p_IPW(Clogg)"

forvalues h = 0/4 {
    local hd  = `h' + 1
    local row = `h' + 2

    * OLS baseline: JOINT LP with both type dummies on the FULL sample, tranquil
    * the omitted category. This is the reference paper's baseline exactly:
    *     reg g_h dum1 dum2 dum3 g_0 $convar, noconstant
    * all types in one regression, sample_for* NOT applied. The rival-drop belongs
    * to their two-stage design (probit + weighted reg) and is applied to the IPW
    * lines below, not here.
    * Country FE only, no year FE (paper-aligned two-stage outcome regression).
    quietly areg dy_`h' onset_nd onset_def `controls_lp' ///
        if sample == 1, absorb(cid) vce(cluster cid)
    matrix b_nd_ols[`row',1]    = _b[onset_nd]
    matrix lo90_nd_ols[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
    matrix hi90_nd_ols[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
    matrix b_def_ols[`row',1]   = _b[onset_def]
    matrix lo90_def_ols[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
    matrix hi90_def_ols[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
    local b_nd_o  = _b[onset_nd]
    local b_def_o = _b[onset_def]
    quietly test onset_nd = onset_def
    matrix pval_act2_ols[`row',1] = r(p)
    local p_o = r(p)

    * IPW line 1: non-default vs tranquil (default dropped), weight ipw_nd
    quietly areg dy_`h' onset_nd `controls_lp' ///
        [aw=ipw_nd] if sample==1 & onset_def==0 & !missing(ipw_nd), ///
        absorb(cid) vce(cluster cid)
    matrix b_nd_ipw[`row',1]    = _b[onset_nd]
    matrix lo90_nd_ipw[`row',1] = _b[onset_nd] - 1.645*_se[onset_nd]
    matrix hi90_nd_ipw[`row',1] = _b[onset_nd] + 1.645*_se[onset_nd]
    local b_nd_w  = _b[onset_nd]
    local se_nd_w = _se[onset_nd]

    * IPW line 2: default-linked vs tranquil (non-default dropped), weight ipw_def
    quietly areg dy_`h' onset_def `controls_lp' ///
        [aw=ipw_def] if sample==1 & onset_nd==0 & !missing(ipw_def), ///
        absorb(cid) vce(cluster cid)
    matrix b_def_ipw[`row',1]    = _b[onset_def]
    matrix lo90_def_ipw[`row',1] = _b[onset_def] - 1.645*_se[onset_def]
    matrix hi90_def_ipw[`row',1] = _b[onset_def] + 1.645*_se[onset_def]
    local b_def_w  = _b[onset_def]
    local se_def_w = _se[onset_def]

    * extra cost of default (IPW) = def - nd, Clogg z on the two vs-tranquil lines
    local zdiff = (`b_def_w' - `b_nd_w') / sqrt(`se_nd_w'^2 + `se_def_w'^2)
    local p_w   = 2*(1 - normal(abs(`zdiff')))
    matrix pval_act2_ipw[`row',1] = `p_w'

    di "h=" `hd' "  " %7.3f `b_nd_o' "  " %7.3f `b_def_o' "      " %5.3f `p_o' ///
              "     " %7.3f `b_nd_w' "  " %7.3f `b_def_w' "      " %5.3f `p_w'
}

* ── Save Act 2 IPW results ────────────────────────────────────────────────
clear
set obs 6
gen horizon = _n - 1     // 0 (baseline), 1, 2, 3, 4, 5

foreach m in b_nd b_def lo90_nd hi90_nd lo90_def hi90_def {
    svmat `m'_ipw, names(`m')
    rename `m'1 `m'
}
gen series = "ipw"
save "$clean/irf_act2_ipw.dta", replace

clear
set obs 6
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
    xlabel(0(1)5, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Year (Year 1 = crisis year)", size(medsmall)) ///
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
