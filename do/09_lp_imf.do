/*===========================================================================
  09_LP_IMF.DO
  ACT 3 — Local Projections: Role of IMF Programs

  Question: Does having an IMF program at crisis onset attenuate the
  output cost of sovereign spread crises?

  Specification A: separate LP for each group
    - onset_imf    = crisis onset with IMF program (imf==1 at onset)
    - onset_no_imf = crisis onset without IMF program (imf==0 at onset)

  Specification B: joint regression with both dummies — tests equality
    H0: beta_imf(h) = beta_no_imf(h)

  Note: imf variable = 1 if country has an active IMF program in year t.
  We interact with onset_all to get crisis-specific IMF dummies.

  Saves:
    "$clean/irf_imf.dta"     — IMF-supported crisis IRF
    "$clean/irf_no_imf.dta"  — non-IMF crisis IRF
    "$tabs/imf_counts.csv"   — episode counts by IMF status
===========================================================================*/

use "$clean/panel_lp.dta", clear

local controls l1_gdpg l2_gdpg ca debt infl vix ust10y
* Note: imf excluded from controls here since it is the treatment modifier

* ══════════════════════════════════════════════════════════════════════════
* CONSTRUCT IMF INTERACTION DUMMIES
* ══════════════════════════════════════════════════════════════════════════

gen onset_imf    = (onset_all == 1 & imf == 1)
gen onset_no_imf = (onset_all == 1 & imf == 0)

label var onset_imf    "Onset with IMF program"
label var onset_no_imf "Onset without IMF program"

* Count episodes by IMF status
quietly count if onset_imf == 1
local n_imf = r(N)
quietly count if onset_no_imf == 1
local n_no_imf = r(N)

di as result _n "=== IMF PROGRAM STATUS AT CRISIS ONSET ==="
di as result "Episodes with IMF program:    `n_imf'"
di as result "Episodes without IMF program: `n_no_imf'"
di as result "Total:                        `=`n_imf' + `n_no_imf''"

* Export counts
clear
set obs 2
gen group = ""
replace group = "IMF program"    in 1
replace group = "No IMF program" in 2
gen n_episodes = .
replace n_episodes = `n_imf'    in 1
replace n_episodes = `n_no_imf' in 2
export delimited "$tabs/imf_counts.csv", replace
use "$clean/panel_lp.dta", clear
gen onset_imf    = (onset_all == 1 & imf == 1)
gen onset_no_imf = (onset_all == 1 & imf == 0)
local controls l1_gdpg l2_gdpg ca debt infl vix ust10y

* ══════════════════════════════════════════════════════════════════════════
* SPEC A-1: CRISES WITH IMF PROGRAM
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== CRISES WITH IMF PROGRAM (N=`n_imf' onsets) ==="

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_imf = J(5, 1, .)
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
    di "h=" `h' ": beta = " %6.3f _b[onset_imf] "  SE = " %6.3f _se[onset_imf] "  p = " %5.3f (2*(1-normal(abs(_b[onset_imf]/_se[onset_imf]))))
}

* ══════════════════════════════════════════════════════════════════════════
* SPEC A-2: CRISES WITHOUT IMF PROGRAM
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== CRISES WITHOUT IMF PROGRAM (N=`n_no_imf' onsets) ==="

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_no_imf = J(5, 1, .)
}

forvalues h = 0/4 {
    local row = `h' + 1
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_no_imf `controls' i.year if sample==1, fe lag(`lag')
    matrix b_no_imf[`row',1]    = _b[onset_no_imf]
    matrix lo90_no_imf[`row',1] = _b[onset_no_imf] - 1.645*_se[onset_no_imf]
    matrix hi90_no_imf[`row',1] = _b[onset_no_imf] + 1.645*_se[onset_no_imf]
    matrix lo95_no_imf[`row',1] = _b[onset_no_imf] - 1.960*_se[onset_no_imf]
    matrix hi95_no_imf[`row',1] = _b[onset_no_imf] + 1.960*_se[onset_no_imf]
    di "h=" `h' ": beta = " %6.3f _b[onset_no_imf] "  SE = " %6.3f _se[onset_no_imf] "  p = " %5.3f (2*(1-normal(abs(_b[onset_no_imf]/_se[onset_no_imf]))))
}

* ══════════════════════════════════════════════════════════════════════════
* SPEC B: JOINT REGRESSION — tests H0: beta_imf(h) = beta_no_imf(h)
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
* BUILD IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

* IMF crisis IRF dataset
clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_imf, names(`m')
    rename `m'1 `m'
}
gen series = "imf"
save "$clean/irf_imf.dta", replace

* No-IMF crisis IRF dataset
clear
set obs 5
gen horizon = _n - 1
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_no_imf, names(`m')
    rename `m'1 `m'
}
gen series = "no_imf"
save "$clean/irf_no_imf.dta", replace

* ══════════════════════════════════════════════════════════════════════════
* FIGURE: IMF vs No-IMF overlay
* ══════════════════════════════════════════════════════════════════════════

use "$clean/irf_imf.dta", clear
append using "$clean/irf_no_imf.dta"

local c_imf    "0 130 79"
local c_no_imf "23 55 94"
local c_zero   "150 150 150"

* Extract p-values for annotation
local p0 = string(pval_imf[1,1], "%4.3f")
local p1 = string(pval_imf[2,1], "%4.3f")
local p2 = string(pval_imf[3,1], "%4.3f")
local p3 = string(pval_imf[4,1], "%4.3f")
local p4 = string(pval_imf[5,1], "%4.3f")

twoway ///
    (connected b horizon if series=="imf", ///
        lcolor("`c_imf'") lwidth(medthick) msymbol(circle) mcolor("`c_imf'")) ///
    (connected b horizon if series=="no_imf", ///
        lcolor("`c_no_imf'") lwidth(medthick) lpattern(dash) ///
        msymbol(square) mcolor("`c_no_imf'")), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years after crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
    title("Output Cost: IMF-Supported vs. Non-IMF Crises", size(medium)) ///
    subtitle("52 EM economies, 1994-2025", size(small)) ///
    legend(order(1 "With IMF program (N=`n_imf')" ///
                 2 "Without IMF program (N=`n_no_imf')") ///
           ring(1) pos(6) cols(2) size(small)) ///
    note("Point estimates only. Driscoll-Kraay SE, lag = h+1. Country + year FE." ///
         "p-values (H0: IMF = No-IMF): h0=`p0' h1=`p1' h2=`p2' h3=`p3' h4=`p4'", ///
         size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig9_irf_imf.pdf",  replace
graph export "$figs/fig9_irf_imf.png",  replace width(1200)
di as result _n "Figure 9 (IMF role) saved."

* ══════════════════════════════════════════════════════════════════════════
* SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

clear
set obs 5
gen horizon = _n - 1
forvalues h = 0/4 {
    local row = `h' + 1
    gen b_imf_h`h'    = .
    gen b_no_imf_h`h' = .
    gen diff_h`h'     = .
    gen pval_h`h'     = .
    replace b_imf_h`h'    = b_imf[`row',1]    in `row'
    replace b_no_imf_h`h' = b_no_imf[`row',1] in `row'
    replace diff_h`h'     = b_imf[`row',1] - b_no_imf[`row',1] in `row'
    replace pval_h`h'     = pval_imf[`row',1] in `row'
}
keep horizon b_imf_h* b_no_imf_h* diff_h* pval_h*
export delimited "$tabs/imf_lp_results.csv", replace
di as result "IMF LP results saved: $tabs/imf_lp_results.csv"
