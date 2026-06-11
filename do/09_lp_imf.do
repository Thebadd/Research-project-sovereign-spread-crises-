/*===========================================================================
  09_LP_IMF.DO
  ACT 3 — Local Projections: Role of IMF Programs

  Question: Does having an IMF program at crisis onset attenuate the
  output cost of sovereign spread crises?

  Specification A: separate LP for each group (OLS baseline)
  Specification B: joint regression with both dummies + equality test
  Specification C: IPW-corrected LP
    - First stage probit: Pr(imf=1 | crisis onset, X)
    - Reweight non-IMF episodes to look like IMF episodes on observables
    - Corrects for selection: IMF countries are systematically different
      (higher debt, lower growth, worse CA) from non-IMF crisis countries

  Note: imf excluded from controls since it is the treatment modifier.

  Saves:
    "$clean/irf_imf.dta"        OLS IRF — IMF-supported crises
    "$clean/irf_no_imf.dta"     OLS IRF — non-IMF crises
    "$clean/irf_imf_ipw.dta"    IPW IRF — IMF-supported crises
    "$clean/irf_no_imf_ipw.dta" IPW IRF — non-IMF crises
    "$tabs/imf_lp_results.csv"  Summary table
    "$tabs/imf_first_stage.csv" First stage probit
===========================================================================*/

use "$clean/panel_lp.dta", clear

local controls l1_gdpg l2_gdpg ca debt infl vix ust10y
* imf excluded from controls — it is the treatment modifier

* ══════════════════════════════════════════════════════════════════════════
* CONSTRUCT IMF INTERACTION DUMMIES
* ══════════════════════════════════════════════════════════════════════════

gen onset_imf    = (onset_all == 1 & imf == 1)
gen onset_no_imf = (onset_all == 1 & imf == 0)

label var onset_imf    "Onset with IMF program"
label var onset_no_imf "Onset without IMF program"

quietly count if onset_imf == 1
local n_imf = r(N)
quietly count if onset_no_imf == 1
local n_no_imf = r(N)

di as result _n "=== IMF PROGRAM STATUS AT CRISIS ONSET ==="
di as result "Episodes with IMF program:    `n_imf'"
di as result "Episodes without IMF program: `n_no_imf'"
di as result "Total:                        `=`n_imf' + `n_no_imf''"

* ══════════════════════════════════════════════════════════════════════════
* SPEC A: SEPARATE LP — OLS BASELINE
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== OLS BASELINE ==="

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_imf    = J(5, 1, .)
    matrix `m'_no_imf = J(5, 1, .)
}

forvalues h = 0/4 {
    local row = `h' + 1
    local lag = max(1, `h'+1)

    xtscc dy_`h' onset_imf `controls' i.year if sample==1, fe lag(`lag')
    matrix b_imf[`row',1]    = _b[onset_imf]
    matrix lo90_imf[`row',1] = _b[onset_imf] - 1.645*_se[onset_imf]
    matrix hi90_imf[`row',1] = _b[onset_imf] + 1.645*_se[onset_imf]
    matrix lo95_imf[`row',1] = _b[onset_imf] - 1.960*_se[onset_imf]
    matrix hi95_imf[`row',1] = _b[onset_imf] + 1.960*_se[onset_imf]

    xtscc dy_`h' onset_no_imf `controls' i.year if sample==1, fe lag(`lag')
    matrix b_no_imf[`row',1]    = _b[onset_no_imf]
    matrix lo90_no_imf[`row',1] = _b[onset_no_imf] - 1.645*_se[onset_no_imf]
    matrix hi90_no_imf[`row',1] = _b[onset_no_imf] + 1.645*_se[onset_no_imf]
    matrix lo95_no_imf[`row',1] = _b[onset_no_imf] - 1.960*_se[onset_no_imf]
    matrix hi95_no_imf[`row',1] = _b[onset_no_imf] + 1.960*_se[onset_no_imf]
}

* ══════════════════════════════════════════════════════════════════════════
* SPEC B: JOINT REGRESSION — equality test
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== JOINT REGRESSION: test beta_imf = beta_no_imf ==="
di "h   beta_imf   beta_no_imf   diff    p(diff)"

matrix pval_imf = J(5, 1, .)

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_imf onset_no_imf `controls' i.year if sample==1, fe lag(`lag')
    test onset_imf = onset_no_imf
    matrix pval_imf[`h'+1, 1] = r(p)
    di "h=" `h' ":  beta_imf="    %6.3f _b[onset_imf] ///
               "   beta_no_imf=" %6.3f _b[onset_no_imf] ///
               "   diff="        %6.3f (_b[onset_imf] - _b[onset_no_imf]) ///
               "   p(diff)="     %5.3f r(p)
}

* ══════════════════════════════════════════════════════════════════════════
* SPEC C: IPW — correct for selection into IMF program
*
* Countries with IMF programs at crisis onset are systematically different:
* higher debt, lower growth, worse current account. The raw comparison
* conflates the effect of IMF support with pre-existing vulnerability.
* IPW reweights non-IMF episodes to resemble IMF episodes on observables.
* First stage estimated on crisis onsets only (onset_all==1).
* Parsimonious spec: use only significant predictors to avoid overfitting
* with small N (~61 onset observations).
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== ACT 3 IPW: FIRST STAGE — Probit of IMF program at onset ==="

* Full first stage to identify significant predictors
probit onset_imf `controls' if onset_all == 1, vce(robust)
estimates store fs_imf_full
di as result "Full model Pseudo-R2: " e(r2_p)

* Parsimonious first stage using significant predictors only
* (inspect full model output above to confirm which are significant)
quietly probit onset_imf debt l1_gdpg ca if onset_all == 1, vce(robust)
local pr2 = e(r2_p)
di as result "Parsimonious Pseudo-R2 (debt, l1_gdpg, ca): " %5.4f `pr2'

probit onset_imf debt l1_gdpg ca if onset_all == 1, vce(robust)
estimates store fs_imf

quietly esttab fs_imf using "$tabs/imf_first_stage.csv", ///
    se star(* 0.10 ** 0.05 *** 0.01) ///
    title("Act 3 First Stage: Pr(IMF Program | Crisis Onset)") replace

* Predict propensity scores
predict pscore3 if onset_all == 1, pr
label var pscore3 "Pr(IMF program | crisis onset, X)"
summarize pscore3, detail

* Trim extreme scores
gen trimmed3 = (pscore3 < 0.05 | pscore3 > 0.95) if !missing(pscore3)
quietly count if trimmed3 == 1
di as result "Observations trimmed [0.05, 0.95]: " r(N)
replace pscore3 = . if trimmed3 == 1

* Stabilized IPW weights
quietly summarize onset_imf if onset_all == 1
local p_imf    = r(mean)
local p_no_imf = 1 - `p_imf'
di as result "Share with IMF program among crisis onsets: " %5.4f `p_imf'

gen ipw3 = .
replace ipw3 = `p_imf'    / pscore3          if onset_imf    == 1 & !missing(pscore3)
replace ipw3 = `p_no_imf' / (1 - pscore3)    if onset_no_imf == 1 & !missing(pscore3)
label var ipw3 "Act 3 stabilized IPW weight"
summarize ipw3 if onset_all == 1, detail

* IPW-weighted LP
di as result _n "=== ACT 3 IPW RESULTS ==="
di "h   b_imf_OLS  b_no_imf_OLS  b_imf_IPW  b_no_imf_IPW  delta_imf  delta_no_imf"

foreach m in b_imf_ipw lo90_imf_ipw hi90_imf_ipw ///
             b_no_imf_ipw lo90_no_imf_ipw hi90_no_imf_ipw {
    matrix `m' = J(5, 1, .)
}

forvalues h = 0/4 {
    local row = `h' + 1

    * Set tranquil year weights to 1
    quietly replace ipw3 = 1 if onset_all == 0 & sample == 1

    quietly areg dy_`h' onset_imf onset_no_imf `controls' i.year ///
        [aw=ipw3] if sample == 1 & !missing(ipw3), absorb(cid) vce(cluster cid)

    quietly replace ipw3 = . if onset_all == 0

    matrix b_imf_ipw[`row',1]      = _b[onset_imf]
    matrix lo90_imf_ipw[`row',1]   = _b[onset_imf]    - 1.645*_se[onset_imf]
    matrix hi90_imf_ipw[`row',1]   = _b[onset_imf]    + 1.645*_se[onset_imf]
    matrix b_no_imf_ipw[`row',1]   = _b[onset_no_imf]
    matrix lo90_no_imf_ipw[`row',1] = _b[onset_no_imf] - 1.645*_se[onset_no_imf]
    matrix hi90_no_imf_ipw[`row',1] = _b[onset_no_imf] + 1.645*_se[onset_no_imf]

    di "h=" `h' "  " %7.3f b_imf[`row',1]    "    " %7.3f b_no_imf[`row',1] ///
               "    " %7.3f _b[onset_imf]     "    " %7.3f _b[onset_no_imf] ///
               "    " %7.3f (_b[onset_imf]    - b_imf[`row',1]) ///
               "    " %7.3f (_b[onset_no_imf] - b_no_imf[`row',1])
}

* ══════════════════════════════════════════════════════════════════════════
* BUILD IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

* OLS datasets
clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_imf, names(`m')
    rename `m'1 `m'
}
gen series = "imf_ols"
save "$clean/irf_imf.dta", replace

clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_no_imf, names(`m')
    rename `m'1 `m'
}
gen series = "no_imf_ols"
save "$clean/irf_no_imf.dta", replace

* IPW datasets
clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 {
    svmat `m'_imf_ipw, names(`m')
    rename `m'1 `m'
}
gen series = "imf_ipw"
save "$clean/irf_imf_ipw.dta", replace

clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 {
    svmat `m'_no_imf_ipw, names(`m')
    rename `m'1 `m'
}
gen series = "no_imf_ipw"
save "$clean/irf_no_imf_ipw.dta", replace

* ══════════════════════════════════════════════════════════════════════════
* FIGURE 9A: OLS baseline — IMF vs No-IMF with CI bands
* ══════════════════════════════════════════════════════════════════════════

use "$clean/irf_imf.dta", clear
append using "$clean/irf_no_imf.dta"

local c_imf    "0 130 79"
local c_no_imf "23 55 94"
local c_zero   "150 150 150"

local p0 = string(pval_imf[1,1], "%4.3f")
local p1 = string(pval_imf[2,1], "%4.3f")
local p2 = string(pval_imf[3,1], "%4.3f")
local p3 = string(pval_imf[4,1], "%4.3f")
local p4 = string(pval_imf[5,1], "%4.3f")

twoway ///
    (rarea lo90 hi90 horizon if series=="imf_ols", ///
        color("`c_imf'%20") lwidth(none)) ///
    (connected b horizon if series=="imf_ols", ///
        lcolor("`c_imf'") lwidth(medthick) msymbol(circle) mcolor("`c_imf'")) ///
    (rarea lo90 hi90 horizon if series=="no_imf_ols", ///
        color("`c_no_imf'%20") lwidth(none)) ///
    (connected b horizon if series=="no_imf_ols", ///
        lcolor("`c_no_imf'") lwidth(medthick) lpattern(dash) ///
        msymbol(square) mcolor("`c_no_imf'")), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years after crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
    title("Output Cost: IMF-Supported vs. Non-IMF Crises", size(medium)) ///
    subtitle("52 EM economies, 1994-2025. OLS baseline.", size(small)) ///
    legend(order(2 "With IMF program (N=`n_imf')" ///
                 4 "Without IMF program (N=`n_no_imf')") ///
           ring(1) pos(6) cols(2) size(small)) ///
    note("90% CI shown. Driscoll-Kraay SE, lag = h+1. Country + year FE." ///
         "p-values (H0: IMF = No-IMF): h0=`p0' h1=`p1' h2=`p2' h3=`p3' h4=`p4'", ///
         size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig9a_irf_imf_ols.pdf", replace
graph export "$figs/fig9a_irf_imf_ols.png", replace width(1200)
di as result "Figure 9a (OLS) saved."

* ══════════════════════════════════════════════════════════════════════════
* FIGURE 9B: IPW — IMF vs No-IMF with CI bands
* ══════════════════════════════════════════════════════════════════════════

use "$clean/irf_imf_ipw.dta", clear
append using "$clean/irf_no_imf_ipw.dta"

twoway ///
    (rarea lo90 hi90 horizon if series=="imf_ipw", ///
        color("`c_imf'%20") lwidth(none)) ///
    (connected b horizon if series=="imf_ipw", ///
        lcolor("`c_imf'") lwidth(medthick) msymbol(circle) mcolor("`c_imf'")) ///
    (rarea lo90 hi90 horizon if series=="no_imf_ipw", ///
        color("`c_no_imf'%20") lwidth(none)) ///
    (connected b horizon if series=="no_imf_ipw", ///
        lcolor("`c_no_imf'") lwidth(medthick) lpattern(dash) ///
        msymbol(square) mcolor("`c_no_imf'")), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years after crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
    title("Output Cost: IMF-Supported vs. Non-IMF Crises (IPW)", size(medium)) ///
    subtitle("52 EM economies, 1994-2025. IPW-corrected.", size(small)) ///
    legend(order(2 "With IMF program (N=`n_imf')" ///
                 4 "Without IMF program (N=`n_no_imf')") ///
           ring(1) pos(6) cols(2) size(small)) ///
    note("90% CI shown. IPW weights from probit of IMF program on debt, GDP growth, CA." ///
         "Cluster SE. Country + year FE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig9b_irf_imf_ipw.pdf", replace
graph export "$figs/fig9b_irf_imf_ipw.png", replace width(1200)
di as result "Figure 9b (IPW) saved."

* ══════════════════════════════════════════════════════════════════════════
* EXPORT COUNTS AND SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

clear
set obs 2
gen group = ""
replace group = "IMF program"    in 1
replace group = "No IMF program" in 2
gen n_episodes = .
replace n_episodes = `n_imf'    in 1
replace n_episodes = `n_no_imf' in 2
export delimited "$tabs/imf_counts.csv", replace

clear
set obs 5
gen horizon = _n - 1
forvalues h = 0/4 {
    local row = `h' + 1
    gen b_imf_h`h'       = b_imf[`row',1]
    gen b_no_imf_h`h'    = b_no_imf[`row',1]
    gen b_imf_ipw_h`h'   = b_imf_ipw[`row',1]
    gen b_no_imf_ipw_h`h' = b_no_imf_ipw[`row',1]
    gen pval_h`h'         = pval_imf[`row',1]
}
export delimited "$tabs/imf_lp_results.csv", replace
di as result "IMF LP results saved: $tabs/imf_lp_results.csv"
di as result _n "Act 3 (IMF role) complete."
