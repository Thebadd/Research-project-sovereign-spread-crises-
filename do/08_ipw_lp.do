/*===========================================================================
  08_IPW_LP.DO
  Inverse Probability Weighted Local Projections (IPW-LP)
  Following Jordà & Taylor (2016) "The Time for Austerity"

  Steps:
    1. First-stage probit: Pr(onset_all=1 | X_{t-1})
       Controls: l1_gdpg, l2_gdpg, debt, ca, infl, vix, ust10y, imf
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

local controls l1_gdpg l2_gdpg ca debt infl imf vix ust10y

* ══════════════════════════════════════════════════════════════════════════
* STEP 1 — FIRST-STAGE PROBIT
* Dependent variable: onset_all (=1 in first year of any crisis episode)
* Estimated on estimation sample only (continuation years excluded)
* All RHS variables lagged to ensure predetermination
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== FIRST STAGE: Probit of crisis onset ==="

* For the probit, use one-period-lagged controls (t-1 predicts t onset)
* l1_gdpg and l2_gdpg already lagged; debt, ca, infl, imf, vix, ust10y
* entered as-is (measured at t, predetermined relative to future outcome)

probit onset_all l1_gdpg l2_gdpg debt ca infl imf vix ust10y ///
    if sample == 1, vce(cluster cid)

estimates store first_stage

* Pseudo R-squared and hit rate
estat classification
di as result "McFadden Pseudo-R2: " e(r2_p)

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
    quietly areg dy_`h' onset_all `controls' i.year ///
        if sample==1 & !missing(ipw), absorb(cid) vce(cluster cid)
    matrix b_ols[`row',1]    = _b[onset_all]
    matrix lo90_ols[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
    matrix hi90_ols[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
    matrix lo95_ols[`row',1] = _b[onset_all] - 1.960*_se[onset_all]
    matrix hi95_ols[`row',1] = _b[onset_all] + 1.960*_se[onset_all]
    local b_u = _b[onset_all]
    local se_u = _se[onset_all]

    * IPW-weighted (areg used because xtreg fe requires weights constant within cid)
    quietly areg dy_`h' onset_all `controls' i.year ///
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
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall))   ///
    title("Baseline vs. IPW-Weighted Local Projections", size(medium))      ///
    subtitle("All spread crises (N=61). xtreg FE, cluster SE.", size(small))///
    legend(order(2 "Baseline OLS-FE" 4 "IPW-weighted")                     ///
           ring(0) pos(3) size(small))                                      ///
    note("IPW weights from probit of onset on lagged macro fundamentals."   ///
         "Stabilized weights, trimmed at [0.01, 0.99].", size(vsmall))      ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig7_ipw_vs_baseline.pdf", replace
graph export "$figs/fig7_ipw_vs_baseline.png", replace width(1200)
di as result _n "Figure 7 saved: $figs/fig7_ipw_vs_baseline.pdf"
di as result "All IPW results saved."
