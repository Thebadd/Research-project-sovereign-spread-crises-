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
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_crisis l_govexp l_open l_credit l_hyperinfl"

* VIX and ust10y are pure time-series variables, fully absorbed by i.year.
local controls $ctrl_core

* Pre-trend controls = the core MINUS l1_gdpg. The placebo outcomes are strictly
* pre-onset growth windows ending at t-1 (dy_m1 = 1yr, dy_m2 = 2yr), and dy_m1 is
* identically l1_gdpg while l1_gdpg is a component of dy_m2 — so keeping it on the
* RHS would regress the placebo on itself. Derived from $ctrl_core to stay in sync.
local cc         $ctrl_core
local dropgdpg   l1_gdpg
local controls_pre : list cc - dropgdpg

* ══════════════════════════════════════════════════════════════════════════
* HELPERS (identical to 02_lp_all.do — see the commentary there)
*   _critvals / _pval : t critical values and p-values on e(df_r), so the CIs
*     match the estimator's own table instead of assuming the normal with ~50
*     clusters. Falls back to the normal when no df is posted.
*   _nepcount         : episode counts inside e(sample), i.e. after listwise
*     deletion on the controls — Table 2 previously reported 35/16 where the
*     regressions identify off 38/21.
* ══════════════════════════════════════════════════════════════════════════

capture program drop _critvals
program define _critvals, rclass
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar df  = .
        return scalar c90 = 1.645
        return scalar c95 = 1.960
    }
    else {
        return scalar df  = `dfr'
        return scalar c90 = invttail(`dfr', 0.05)
        return scalar c95 = invttail(`dfr', 0.025)
    }
end

capture program drop _pval
program define _pval, rclass
    args b se
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar p = 2*(1 - normal(abs(`b'/`se')))
    }
    else {
        return scalar p = 2*ttail(`dfr', abs(`b'/`se'))
    }
end

capture program drop _nepcount
program define _nepcount, rclass
    syntax varname(numeric) , Outcome(varname) Controls(varlist)
    tempvar esmp
    quietly gen byte `esmp' = e(sample)
    quietly count if `esmp' == 1
    if r(N) == 0 {
        quietly replace `esmp' = (sample == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar n = r(N)
end

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
    xtscc dy_m`h_neg' onset_nd `controls_pre' i.year if sample==1, fe lag(1)
    local bb = _b[onset_nd]
    local ss = _se[onset_nd]
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    _pval `bb' `ss'
    local pp = r(p)
    matrix b_nd[`row',1]    = `bb'
    matrix lo90_nd[`row',1] = `bb' - `c90'*`ss'
    matrix hi90_nd[`row',1] = `bb' + `c90'*`ss'
    matrix lo95_nd[`row',1] = `bb' - `c95'*`ss'
    matrix hi95_nd[`row',1] = `bb' + `c95'*`ss'
    di "h=-`h_neg': beta = " %6.3f `bb' "  SE = " %6.3f `ss' "  p = " %5.3f `pp'
}

* Main horizons
forvalues h = 0/4 {
    local row = `h' + 3
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_nd `controls' i.year if sample==1, fe lag(`lag')
    local bb = _b[onset_nd]
    local ss = _se[onset_nd]
    local nn = e(N)
    _nepcount onset_nd, outcome(dy_`h') controls(`controls')
    local nep = r(n)
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    _pval `bb' `ss'
    local pp = r(p)
    matrix b_nd[`row',1]    = `bb'
    matrix lo90_nd[`row',1] = `bb' - `c90'*`ss'
    matrix hi90_nd[`row',1] = `bb' + `c90'*`ss'
    matrix lo95_nd[`row',1] = `bb' - `c95'*`ss'
    matrix hi95_nd[`row',1] = `bb' + `c95'*`ss'
    di "h=" `h' ": beta = " %6.3f `bb' "  SE = " %6.3f `ss' "  N = " `nn' ///
       "  episodes = " `nep' "  p = " %5.3f `pp'
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
    xtscc dy_m`h_neg' onset_def `controls_pre' i.year if sample==1, fe lag(1)
    local bb = _b[onset_def]
    local ss = _se[onset_def]
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    _pval `bb' `ss'
    local pp = r(p)
    matrix b_def[`row',1]    = `bb'
    matrix lo90_def[`row',1] = `bb' - `c90'*`ss'
    matrix hi90_def[`row',1] = `bb' + `c90'*`ss'
    matrix lo95_def[`row',1] = `bb' - `c95'*`ss'
    matrix hi95_def[`row',1] = `bb' + `c95'*`ss'
    di "h=-`h_neg': beta = " %6.3f `bb' "  SE = " %6.3f `ss' "  p = " %5.3f `pp'
}

forvalues h = 0/4 {
    local row = `h' + 3
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_def `controls' i.year if sample==1, fe lag(`lag')
    local bb = _b[onset_def]
    local ss = _se[onset_def]
    local nn = e(N)
    _nepcount onset_def, outcome(dy_`h') controls(`controls')
    local nep = r(n)
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    _pval `bb' `ss'
    local pp = r(p)
    matrix b_def[`row',1]    = `bb'
    matrix lo90_def[`row',1] = `bb' - `c90'*`ss'
    matrix hi90_def[`row',1] = `bb' + `c90'*`ss'
    matrix lo95_def[`row',1] = `bb' - `c95'*`ss'
    matrix hi95_def[`row',1] = `bb' + `c95'*`ss'
    di "h=" `h' ": beta = " %6.3f `bb' "  SE = " %6.3f `ss' "  N = " `nn' ///
       "  episodes = " `nep' "  p = " %5.3f `pp'
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

    * Episode counts contributing at this horizon. These must be counted inside
    * e(sample): an onset with a non-missing outcome but a missing control is
    * dropped by listwise deletion and identifies nothing, which is why the old
    * counts (35/16) disagreed with the treated counts the probit reports.
    _nepcount onset_nd,  outcome(dy_`h') controls(`controls')
    local nepnd = r(n)
    _nepcount onset_def, outcome(dy_`h') controls(`controls')
    local nepdef = r(n)

    * Clogg et al. (1995) difference: default-linked minus non-default
    * (negative => default-linked loss is deeper). z uses the pooled SE of the
    * two independent-within-regression coefficients, referred to the estimator's
    * residual df for consistency with the CIs.
    local bdiff = `bdef' - `bnd'
    local zdiff = `bdiff' / sqrt(`snd'^2 + `sdef'^2)
    _pval `bdiff' `=sqrt(`snd'^2 + `sdef'^2)'
    local pz = r(p)

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
* ROBUSTNESS: OUTTURNS ONLY — is the default premium a forecast artifact?
*
* The premium is largest at h=3/h=4, and those are exactly the horizons where an
* onset from 2022 onward has its outcome built from WEO projections rather than
* actuals. This re-runs the joint spec keeping only observations whose outcome
* year t+h is an outturn for that country. If the def-nd gap survives in sign and
* rough magnitude, the persistence result is real; if it collapses, it was the
* forecast. Either answer is worth having, and the sample shrinks either way.
* ══════════════════════════════════════════════════════════════════════════

capture confirm variable gdp_last_actual
if _rc {
    di as error _n "  ** gdp_last_actual not in the panel — outturn-only robustness skipped."
}
else {
    di as result _n "=== ROBUSTNESS: OUTTURNS ONLY (joint spec) ==="
    di as result "h   beta_nd   beta_def   diff(def-nd)   Clogg z   p(Clogg)   N   nd/def episodes"

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)

        capture xtscc dy_`h' onset_nd onset_def `controls' i.year ///
            if sample==1 & (year + `h') <= gdp_last_actual, fe lag(`lag')

        if _rc {
            di as error "  outturn-only joint regression failed at h=`h' (rc=" _rc ")"
            continue
        }

        local bnd  = _b[onset_nd]
        local bdef = _b[onset_def]
        local snd  = _se[onset_nd]
        local sdef = _se[onset_def]
        local nn   = e(N)

        local bdiff = `bdef' - `bnd'
        local sdiff = sqrt(`snd'^2 + `sdef'^2)
        local zdiff = `bdiff' / `sdiff'
        _pval `bdiff' `sdiff'
        local pz = r(p)

        _nepcount onset_nd,  outcome(dy_`h') controls(`controls')
        local rnepnd = r(n)
        _nepcount onset_def, outcome(dy_`h') controls(`controls')
        local rnepdef = r(n)

        eststo t2r_h`h'
        estadd scalar bdiff  = `bdiff'
        estadd scalar zdiff  = `zdiff'
        estadd scalar pzdiff = `pz'
        estadd scalar nepnd  = `rnepnd'
        estadd scalar nepdef = `rnepdef'

        di "h=" `h' "  " %7.3f `bnd' "  " %8.3f `bdef' "  " %10.3f `bdiff' ///
           "  " %8.2f `zdiff' "  " %8.3f `pz' "  " %5.0f `nn' ///
           "   " %3.0f `rnepnd' "/" %3.0f `rnepdef'
    }

    capture esttab t2r_h0 t2r_h1 t2r_h2 t2r_h3 t2r_h4 ///
        using "$tabs/table2r_output_resolution_outturns.rtf", replace ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_nd onset_def) order(onset_nd onset_def) ///
        coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
        mtitles("h=0" "h=1" "h=2" "h=3" "h=4") nonumber ///
        stats(bdiff zdiff pzdiff nepnd nepdef N N_g, ///
              labels("Difference (default - non-default)" "  Clogg et al. (1995) z" ///
                     "  p (Clogg z)" "Episodes (non-default)" "Episodes (default)" ///
                     "Observations" "Countries") ///
              fmt(3 3 3 0 0 0 0)) ///
        title("Table 2R. Output cost by resolution — outturns only (robustness)") ///
        addnotes("As Table 2, but each horizon drops observations whose outcome year t+h falls after the last WEO actual observation for that country." ///
                 "This removes IMF projections from the dependent variable, at the cost of sample size at longer horizons." ///
                 "Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01.")

    if _rc == 608 di as error "  ** table2r_output_resolution_outturns.rtf is OPEN IN WORD — close it and re-run."
    else if _rc di as error "  ** Table 2R: esttab failed (rc=" _rc ")"
    else di as result "Table 2R saved: $tabs/table2r_output_resolution_outturns.rtf"
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
                 "Episodes (non-default, in estimation sample)" "Episodes (default, in estimation sample)" ///
                 "Observations" "Countries") ///
          fmt(3 3 3 3 0 0 0 0)) ///
    title("Table 2. Output cost by crisis resolution: non-default vs. default-linked") ///
    addnotes("Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h." ///
             "Both onset dummies enter jointly. Jorda (2005) local projections; country and year fixed effects; continuation years excluded." ///
             "Driscoll-Kraay standard errors in parentheses." ///
             "Difference = beta(default) - beta(non-default); negative means the default-linked loss is deeper." ///
             "Clogg z = difference / sqrt(SE_nd^2 + SE_def^2); p(Wald) is the joint equality test. Bootstrap CIs to be added (upgrade D)." ///
             "Episode counts are onsets surviving listwise deletion on the control set, i.e. those actually identifying the coefficients." ///
             "p-values and confidence intervals use t critical values on the estimator's residual degrees of freedom." ///
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
