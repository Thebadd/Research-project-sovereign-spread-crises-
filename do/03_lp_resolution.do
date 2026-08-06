/*===========================================================================
  03_LP_RESOLUTION.DO
  ACT 2 — Local Projections: Non-Default vs. Default-Linked Episodes

  Specification A: separate LP for each group
  Specification B: joint regression with both dummies — tests equality of IRFs

  Treatment dummies:
    onset_nd  = 1 → 40 non-default episodes
    onset_def = 1 → 21 default-linked episodes

  Saves:
    "$clean/irf_nd.dta"    — non-default IRF
    "$clean/irf_def.dta"   — default-linked IRF
    "$clean/irf_joint.dta" — combined for overlay plot
===========================================================================*/

use "$clean/panel_lp.dta", clear

* VIX and ust10y are pure time-series variables, fully absorbed by i.year.
local controls $ctrl_core

* ══════════════════════════════════════════════════════════════════════════
* SPEC A-1: NON-DEFAULT EPISODES ONLY
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== NON-DEFAULT EPISODES (N=40 onsets) ==="

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

di as result _n "=== DEFAULT-LINKED EPISODES (N=21 onsets) ==="

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

eststo clear   // clear any stored estimates before capturing for Table 2

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_nd onset_def `controls' i.year if sample==1, fe lag(`lag')

    * Coefficients and SEs
    local bnd  = _b[onset_nd]
    local bdef = _b[onset_def]
    local snd  = _se[onset_nd]
    local sdef = _se[onset_def]

    * Wald equality test (kept for continuity)
    test onset_nd = onset_def
    matrix pval_diff[`h'+1, 1] = r(p)
    local pd = r(p)   // store before eststo (which can reset r())

    * Episode counts contributing at this horizon (onset obs, non-missing outcome)
    quietly count if onset_nd  == 1 & sample == 1 & !missing(dy_`h')
    local nepnd = r(N)
    quietly count if onset_def == 1 & sample == 1 & !missing(dy_`h')
    local nepdef = r(N)

    * Clogg et al. (1995) difference: default-linked minus non-default
    * (negative => default-linked loss is deeper). z uses the pooled SE of the
    * two independent-within-regression coefficients.
    local bdiff = `bdef' - `bnd'
    local zdiff = `bdiff' / sqrt(`snd'^2 + `sdef'^2)
    local pz    = 2*(1 - normal(abs(`zdiff')))

    * Capture estimates + difference block for the publication table (Table 2)
    eststo t2_h`h'
    estadd scalar bdiff  = `bdiff'
    estadd scalar zdiff  = `zdiff'
    estadd scalar pzdiff = `pz'
    estadd scalar pdiff  = `pd'
    estadd scalar nepnd  = `nepnd'
    estadd scalar nepdef = `nepdef'

    di "h=" `h' ":  beta_nd=" %6.3f `bnd' ///
               "  beta_def=" %6.3f `bdef' ///
               "  diff(def-nd)=" %6.3f `bdiff' ///
               "  Clogg z=" %5.2f `zdiff' ///
               "  p(Wald)="  %5.3f `pd'
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — TABLE 2: Output cost by resolution type
*   Word/RTF. Columns = horizons h=0..4; rows = non-default and default-linked
*   onset coefficients (entered jointly). Extra row: p-value of H0 beta_nd=beta_def.
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

capture esttab t2_h0 t2_h1 t2_h2 t2_h3 t2_h4 using "$tabs/table2_output_resolution.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(onset_nd onset_def) order(onset_nd onset_def) ///
    coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
    mtitles("h=0" "h=1" "h=2" "h=3" "h=4") nonumber ///
    stats(bdiff zdiff pzdiff pdiff nepnd nepdef N N_g, ///
          labels("Difference (default - non-default)" "  Clogg et al. (1995) z" ///
                 "  p (Clogg z)" "  p (Wald, nd = def)" ///
                 "Episodes (non-default)" "Episodes (default)" ///
                 "Observations" "Countries") ///
          fmt(3 3 3 3 0 0 0 0)) ///
    title("Table 2. Output cost by crisis resolution: non-default vs. default-linked") ///
    addnotes("Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h." ///
             "Both onset dummies enter jointly. Jorda (2005) local projections; country and year fixed effects; continuation years excluded." ///
             "Driscoll-Kraay standard errors in parentheses." ///
             "Difference = beta(default) - beta(non-default); negative means the default-linked loss is deeper." ///
             "Clogg z = difference / sqrt(SE_nd^2 + SE_def^2); p(Wald) is the joint equality test. Bootstrap CIs to be added (upgrade D)." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table2_output_resolution.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc di as error "  ** Table 2: esttab failed (rc=" _rc ")"
else di as result "Table 2 saved: $tabs/table2_output_resolution.rtf"

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
