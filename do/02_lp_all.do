/*===========================================================================
  02_LP_ALL.DO
  ACT 1 — Local Projections: ALL Spread Crisis Episodes (N=61)

  Treatment: onset_all = 1 in first year of any spread crisis (crit1 or crit2)
  Sample:    onset + tranquil years (continuation years excluded)
  Horizons:  h = −2, −1 (pre-trend placebo), 0, 1, 2, 3, 4
  SE:        Driscoll-Kraay with lag = max(1, h+1)  [xtscc, fe]

  Saves:
    - Matrices: b_all, lo90_all, hi90_all, lo95_all, hi95_all (7×1, rows=h)
    - Dataset:  "$clean/irf_all.dta"  (for graphing)
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking l_govexp l_open l_credit hyperinf_dummy"

* ── Controls (pre-determined at t, all lagged relative to outcome) ────────
* VIX and ust10y have zero cross-sectional variation and are fully absorbed by
* year fixed effects (i.year). Including them is redundant; they are omitted.
local controls $ctrl_core

* ── Storage matrices (rows: h = -2, -1, 0, 1, 2, 3, 4) ──────────────────
local nhor = 7
foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_all = J(`nhor', 1, .)
}

* ── Row index: 1=h-2, 2=h-1, 3=h0, 4=h1, 5=h2, 6=h3, 7=h4 ─────────────

* ────────────────────────────────────────────────────────────────────────
* PRE-TREND PLACEBO: h = -2, h = -1
* If identification is valid, these coefficients should be ~0
* ────────────────────────────────────────────────────────────────────────

di as result _n "=== PRE-TREND TEST (placebo horizons) ==="

foreach h_neg in 2 1 {
    local row = (3 - `h_neg')   // h=-2 → row 1, h=-1 → row 2

    xtscc dy_m`h_neg' onset_all `controls' i.year ///
        if sample == 1, fe lag(1)

    matrix b_all[`row', 1]    = _b[onset_all]
    matrix lo90_all[`row', 1] = _b[onset_all] - 1.645 * _se[onset_all]
    matrix hi90_all[`row', 1] = _b[onset_all] + 1.645 * _se[onset_all]
    matrix lo95_all[`row', 1] = _b[onset_all] - 1.960 * _se[onset_all]
    matrix hi95_all[`row', 1] = _b[onset_all] + 1.960 * _se[onset_all]

    di "h=-`h_neg': beta = " %6.3f _b[onset_all] ///
       "  SE = " %6.3f _se[onset_all] ///
       "  p = " %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all]))))
}

* ────────────────────────────────────────────────────────────────────────
* MAIN HORIZONS: h = 0, 1, 2, 3, 4
* ────────────────────────────────────────────────────────────────────────

di as result _n "=== MAIN LP RESULTS (ALL CRISES) ==="

eststo clear   // clear any stored estimates before capturing for Table 1

forvalues h = 0/4 {
    local row = `h' + 3   // h=0 → row 3, ..., h=4 → row 7
    local lag = max(1, `h' + 1)

    xtscc dy_`h' onset_all `controls' i.year ///
        if sample == 1, fe lag(`lag')

    * Capture estimates for the publication table (Table 1)
    quietly count if onset_all == 1 & sample == 1 & !missing(dy_`h')
    local nep = r(N)
    eststo t1_h`h'
    estadd scalar nep = `nep'

    * Store point estimate and confidence intervals
    matrix b_all[`row', 1]    = _b[onset_all]
    matrix lo90_all[`row', 1] = _b[onset_all] - 1.645 * _se[onset_all]
    matrix hi90_all[`row', 1] = _b[onset_all] + 1.645 * _se[onset_all]
    matrix lo95_all[`row', 1] = _b[onset_all] - 1.960 * _se[onset_all]
    matrix hi95_all[`row', 1] = _b[onset_all] + 1.960 * _se[onset_all]

    di "h=" `h' ": beta = " %6.3f _b[onset_all] ///
       "  SE = " %6.3f _se[onset_all] ///
       "  N = " e(N) ///
       "  p = " %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all]))))
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — TABLE 1: Output cost of spread crises (all episodes)
*   Word/RTF. Columns = horizons h=0..4; row = onset coefficient.
*   Coefficient with Driscoll-Kraay SE in parentheses below; stars from e(V).
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

capture esttab t1_h0 t1_h1 t1_h2 t1_h3 t1_h4 using "$tabs/table1_output_all.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(onset_all) coeflabel(onset_all "Spread-crisis onset") ///
    mtitles("h=0" "h=1" "h=2" "h=3" "h=4") nonumber ///
    stats(nep N N_g, labels("Episodes (onsets)" "Observations" "Countries") fmt(0 0 0)) ///
    title("Table 1. Output cost of sovereign spread crises (all episodes)") ///
    addnotes("Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h." ///
             "Jorda (2005) local projections. Country and year fixed effects; continuation years excluded." ///
             "Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table1_output_all.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc di as error "  ** Table 1: esttab failed (rc=" _rc ")"
else di as result "Table 1 saved: $tabs/table1_output_all.rtf"

* ────────────────────────────────────────────────────────────────────────
* BUILD IRF DATASET FOR GRAPHING
* ────────────────────────────────────────────────────────────────────────

clear
local nhor = 7
set obs `nhor'

gen horizon = _n - 3     // -2, -1, 0, 1, 2, 3, 4

svmat b_all,    names(b)
svmat lo90_all, names(lo90)
svmat hi90_all, names(hi90)
svmat lo95_all, names(lo95)
svmat hi95_all, names(hi95)

rename b1    b
rename lo901 lo90
rename hi901 hi90
rename lo951 lo95
rename hi951 hi95

gen series = "all"
label var horizon "Horizon (years after onset)"
label var b       "Point estimate (pp)"
label var lo90    "90% CI lower"
label var hi90    "90% CI upper"
label var lo95    "95% CI lower"
label var hi95    "95% CI upper"

save "$clean/irf_all.dta", replace
di as result _n "IRF data saved: $clean/irf_all.dta"
