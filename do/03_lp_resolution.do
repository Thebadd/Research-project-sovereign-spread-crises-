/*===========================================================================
  03_LP_RESOLUTION.DO
  ACT 2 — Local Projections: Non-Default vs. Default-Linked Episodes

  Specification A: separate LP for each group
  Specification B: joint regression with both dummies — tests equality of IRFs

  Treatment dummies:
    onset_nd  = 1 → 39 non-default episodes
    onset_def = 1 → 22 default-linked episodes

  Saves:
    "$clean/irf_nd.dta"    — non-default IRF
    "$clean/irf_def.dta"   — default-linked IRF
    "$clean/irf_joint.dta" — combined for overlay plot
===========================================================================*/

use "$clean/panel_lp.dta", clear

local controls l1_gdpg l2_gdpg ca debt infl imf vix ust10y

* ══════════════════════════════════════════════════════════════════════════
* SPEC A-1: NON-DEFAULT EPISODES ONLY
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== NON-DEFAULT EPISODES (N=39 onsets) ==="

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_nd = J(7, 1, .)
}

* Pre-trend
foreach h_neg in 2 1 {
    local row = 3 - `h_neg'
    xtscc dy_m`h_neg' onset_nd `controls' i.year if sample==1, fe lag(1)
    matrix b_nd[`row',1]    = _b[onset_nd]
    matrix lo90_nd[`row',1] = _b[onset_nd] - 1.645*_se[onset_nd]
    matrix hi90_nd[`row',1] = _b[onset_nd] + 1.645*_se[onset_nd]
    matrix lo95_nd[`row',1] = _b[onset_nd] - 1.960*_se[onset_nd]
    matrix hi95_nd[`row',1] = _b[onset_nd] + 1.960*_se[onset_nd]
    di "h=-`h_neg': beta = " %6.3f _b[onset_nd] "  p = " %5.3f (2*(1-normal(abs(_b[onset_nd]/_se[onset_nd]))))
}

* Main horizons
forvalues h = 0/4 {
    local row = `h' + 3
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_nd `controls' i.year if sample==1, fe lag(`lag')
    matrix b_nd[`row',1]    = _b[onset_nd]
    matrix lo90_nd[`row',1] = _b[onset_nd] - 1.645*_se[onset_nd]
    matrix hi90_nd[`row',1] = _b[onset_nd] + 1.645*_se[onset_nd]
    matrix lo95_nd[`row',1] = _b[onset_nd] - 1.960*_se[onset_nd]
    matrix hi95_nd[`row',1] = _b[onset_nd] + 1.960*_se[onset_nd]
    di "h=" `h' ": beta = " %6.3f _b[onset_nd] "  SE = " %6.3f _se[onset_nd] "  N = " e(N)
}

* ══════════════════════════════════════════════════════════════════════════
* SPEC A-2: DEFAULT-LINKED EPISODES ONLY
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== DEFAULT-LINKED EPISODES (N=22 onsets) ==="

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_def = J(7, 1, .)
}

foreach h_neg in 2 1 {
    local row = 3 - `h_neg'
    xtscc dy_m`h_neg' onset_def `controls' i.year if sample==1, fe lag(1)
    matrix b_def[`row',1]    = _b[onset_def]
    matrix lo90_def[`row',1] = _b[onset_def] - 1.645*_se[onset_def]
    matrix hi90_def[`row',1] = _b[onset_def] + 1.645*_se[onset_def]
    matrix lo95_def[`row',1] = _b[onset_def] - 1.960*_se[onset_def]
    matrix hi95_def[`row',1] = _b[onset_def] + 1.960*_se[onset_def]
    di "h=-`h_neg': beta = " %6.3f _b[onset_def] "  p = " %5.3f (2*(1-normal(abs(_b[onset_def]/_se[onset_def]))))
}

forvalues h = 0/4 {
    local row = `h' + 3
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_def `controls' i.year if sample==1, fe lag(`lag')
    matrix b_def[`row',1]    = _b[onset_def]
    matrix lo90_def[`row',1] = _b[onset_def] - 1.645*_se[onset_def]
    matrix hi90_def[`row',1] = _b[onset_def] + 1.645*_se[onset_def]
    matrix lo95_def[`row',1] = _b[onset_def] - 1.960*_se[onset_def]
    matrix hi95_def[`row',1] = _b[onset_def] + 1.960*_se[onset_def]
    di "h=" `h' ": beta = " %6.3f _b[onset_def] "  SE = " %6.3f _se[onset_def] "  N = " e(N)
}

* ══════════════════════════════════════════════════════════════════════════
* SPEC B: JOINT REGRESSION — tests H0: beta_nd(h) = beta_def(h)
* Both dummies enter simultaneously; difference tests whether resolution
* outcome matters for the magnitude of the output loss.
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== JOINT REGRESSION: test beta_nd = beta_def ==="

matrix pval_diff = J(5, 1, .)   // h=0..4 only (not pre-trend)

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_nd onset_def `controls' i.year if sample==1, fe lag(`lag')

    * Test equality
    test onset_nd = onset_def
    matrix pval_diff[`h'+1, 1] = r(p)

    di "h=" `h' ":  beta_nd=" %6.3f _b[onset_nd] ///
               "  beta_def=" %6.3f _b[onset_def] ///
               "  p(diff)="  %5.3f r(p)
}

* ══════════════════════════════════════════════════════════════════════════
* BUILD IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

* Non-default IRF dataset
clear
set obs 7
gen horizon = _n - 3
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_nd, names(`m')
    rename `m'1 `m'
}
gen series = "nd"
save "$clean/irf_nd.dta", replace

* Default-linked IRF dataset
clear
set obs 7
gen horizon = _n - 3
foreach m in b lo90 hi90 lo95 hi95 {
    svmat `m'_def, names(`m')
    rename `m'1 `m'
}
gen series = "def"
save "$clean/irf_def.dta", replace

* Joint dataset for overlay plot
use "$clean/irf_nd.dta",  clear
append using "$clean/irf_def.dta"
append using "$clean/irf_all.dta"

label var b       "Point estimate (pp)"
label var horizon "Horizon (years after onset)"
save "$clean/irf_joint.dta", replace

di as result _n "All IRF datasets saved."

* Print p-values for difference test
di as result _n "=== P-VALUES: H0: beta_nd(h) = beta_def(h) ==="
forvalues h = 0/4 {
    di "h=" `h' ":  p = " %5.3f pval_diff[`h'+1, 1]
}
