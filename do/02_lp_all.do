/*===========================================================================
  02_LP_ALL.DO
  ACT 1 — Local Projections: ALL Spread Crisis Episodes (N=61)

  Treatment: onset_all = 1 in first year of any spread crisis (crit1 or crit2)
  Sample:    onset + tranquil years (continuation years excluded)
  Horizons:  h = -1 (pre-trend placebo), 0 (baseline, hardcoded 0), 1, 2, 3, 4, 5
             Year numbering matches Asonuma et al.: Year 1 is the crisis year
             itself, Year 0 is the explicit pre-crisis baseline (always 0).
  SE:        Driscoll-Kraay with lag = max(1, h+1)  [xtscc, fe]

  Saves:
    - Matrices: b_all, lo90_all, hi90_all, lo95_all, hi95_all (7×1, rows=h)
    - Dataset:  "$clean/irf_all.dta"  (for graphing)
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"

* ── Controls (pre-determined at t, all lagged relative to outcome) ────────
* VIX and ust10y have zero cross-sectional variation and are fully absorbed by
* year fixed effects (i.year). Including them is redundant; they are omitted.
local controls $ctrl_core

* ══════════════════════════════════════════════════════════════════════════
* HELPERS — small-sample inference and honest episode counts
*
* _critvals: the CIs below were built with the normal 1.645/1.960. With ~50
*   clusters that understates the interval. Both xtscc and areg post e(df_r),
*   so use invttail() on it and our hand-built CIs/p-values then agree with the
*   command's own printed table. Falls back to the normal if no df is posted.
*
* _nepcount: "Episodes" in Table 1 counted onsets with a non-missing OUTCOME,
*   which ignores listwise deletion on $ctrl_core — the table said 51 while the
*   Act-1 probit reported 38 treated. Count onsets inside e(sample) instead.
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
    if r(N) == 0 {          // estimator posted no e(sample) — rebuild it
        quietly replace `esmp' = (sample == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar n = r(N)
end

* ── Storage matrices (rows: h = -2, -1, 0, 1, 2, 3, 4) ──────────────────
local nhor = 7
foreach m in b se lo90 hi90 lo95 hi95 {
    matrix `m'_all = J(`nhor', 1, .)
}

* ── Row index (Asonuma-style display numbering): 1=h-1, 2=h0(baseline,=0),
*    3=h1, 4=h2, 5=h3, 6=h4, 7=h5 ─────────────────────────────────────────

* ────────────────────────────────────────────────────────────────────────
* PRE-TREND PLACEBO: h = -1 (single horizon, display numbering)
* dy_h = F h.ln_gdp - L.ln_gdp for h=0..4 is always measured against the SAME
* t-1 base. Plugging h=-1 into that identical formula gives L.ln_gdp - L.ln_gdp
* = 0 trivially — the base year equals itself, so h=-1 is not an estimable
* placebo, it IS the baseline (matches the standard event-study convention
* where h=-1 is normalized to zero, not estimated). The first REAL pre-crisis
* horizon on this same base is h=-2 (dy_m2, built in 18_transforms). If
* identification is valid, this coefficient should be ~0 — treated countries
* were not already on a different path before the crisis.
*
* CONTROLS DIFFER HERE, deliberately. dy_m2 is algebraically -l1_gdpg (both are
* +/-100*(L.ln_gdp - L2.ln_gdp)), so leaving l1_gdpg on the RHS would regress
* the placebo outcome on itself and guarantee a null. `controls_pre' is the
* common core minus l1_gdpg, derived from $ctrl_core so it stays in sync if
* the core changes.
* ────────────────────────────────────────────────────────────────────────

local cc         $ctrl_core
local dropgdpg   l1_gdpg
local controls_pre : list cc - dropgdpg

di as result _n "=== PRE-TREND TEST (placebo horizon, displayed as h=-1) ==="
di as result "    controls: `controls_pre'"
di as result "    (l1_gdpg excluded — it IS -1 times the dy_m2 outcome; the true baseline, displayed h=0, is never estimated)"

foreach h_neg in 2 {
    local row = (3 - `h_neg')   // dy_m2 -> displayed h=-1 -> row 1

    xtscc dy_m`h_neg' onset_all `controls_pre' i.year ///
        if sample == 1, fe lag(1)

    local bb = _b[onset_all]
    local ss = _se[onset_all]
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    _pval `bb' `ss'
    local pp = r(p)

    matrix b_all[`row', 1]    = `bb'
    matrix se_all[`row', 1]   = `ss'
    matrix lo90_all[`row', 1] = `bb' - `c90' * `ss'
    matrix hi90_all[`row', 1] = `bb' + `c90' * `ss'
    matrix lo95_all[`row', 1] = `bb' - `c95' * `ss'
    matrix hi95_all[`row', 1] = `bb' + `c95' * `ss'

    di "h=-1: beta = " %6.3f `bb' ///
       "  SE = " %6.3f `ss' ///
       "  p = " %5.3f `pp'
}

* ── Explicit baseline point (displayed h=0), matching Asonuma et al.: hard-
* coded zero, never estimated — the pre-crisis anchor every horizon is
* measured against. Row 2 of the 7-row matrices.
foreach m in b se lo90 hi90 lo95 hi95 {
    matrix `m'_all[2, 1] = 0
}

* ────────────────────────────────────────────────────────────────────────
* MAIN HORIZONS: displayed h = 1, 2, 3, 4, 5 (Asonuma numbering; internal
* variables dy_0..dy_4 are unchanged, only the printed/labeled horizon shifts)
* ────────────────────────────────────────────────────────────────────────

di as result _n "=== MAIN LP RESULTS (ALL CRISES) ==="

eststo clear   // clear any stored estimates before capturing for Table 1

forvalues h = 0/4 {
    local hd  = `h' + 1   // displayed horizon (Asonuma: crisis year = 1)
    local row = `h' + 3   // dy_0 → row 3, ..., dy_4 → row 7
    local lag = max(1, `h' + 1)

    xtscc dy_`h' onset_all `controls' i.year ///
        if sample == 1, fe lag(`lag')

    local bb = _b[onset_all]
    local ss = _se[onset_all]
    local nn = e(N)

    * Onsets that actually survive listwise deletion into this regression
    _nepcount onset_all, outcome(dy_`h') controls(`controls')
    local nep = r(n)

    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    local dfr = r(df)
    _pval `bb' `ss'
    local pp = r(p)

    * Capture estimates for the publication table (Table 1)
    eststo t1_h`h'
    estadd scalar nep = `nep'

    * Store point estimate and confidence intervals
    matrix b_all[`row', 1]    = `bb'
    matrix se_all[`row', 1]   = `ss'
    matrix lo90_all[`row', 1] = `bb' - `c90' * `ss'
    matrix hi90_all[`row', 1] = `bb' + `c90' * `ss'
    matrix lo95_all[`row', 1] = `bb' - `c95' * `ss'
    matrix hi95_all[`row', 1] = `bb' + `c95' * `ss'

    di "h=" `hd' ": beta = " %6.3f `bb' ///
       "  SE = " %6.3f `ss' ///
       "  N = " `nn' ///
       "  episodes = " `nep' ///
       "  df = " `dfr' ///
       "  p = " %5.3f `pp'
}

* ══════════════════════════════════════════════════════════════════════════
* ROBUSTNESS: OUTTURNS ONLY (drop horizons that land in the WEO projection window)
*
* The panel runs to 2026 but WEO Apr-2026's last ACTUAL observation is 2024/2025
* for most countries. So for an onset in 2022 or later, dy_3 and dy_4 are built
* from IMF forecasts, not data — and those are precisely the horizons where the
* estimated cost is largest. This block re-runs each horizon keeping only
* observations whose OUTCOME year (t+h) is an outturn for that country.
*
* Read it as a falsification, not a replacement: if the sign and rough magnitude
* survive, the long-horizon result is not a forecast artifact. The sample shrinks
* at h=3/h=4, so wider CIs there are expected and not themselves informative.
* ══════════════════════════════════════════════════════════════════════════

capture confirm variable gdp_last_actual
if _rc {
    di as error _n "  ** gdp_last_actual not in the panel — outturn-only robustness skipped."
    di as error "     (11_weo.do could not read LATEST_ACTUAL_ANNUAL_DATA from the WEO export.)"
}
else {
    di as result _n "=== ROBUSTNESS: OUTTURNS ONLY (no WEO projections in the outcome) ==="
    di as result "h   beta     SE      N     episodes   p"

    * NOTE: no eststo clear here — the headline t1_* estimates must survive for
    * the Table 1 export below. The robustness set is stored under t1r_* instead.
    forvalues h = 0/4 {
        local hd  = `h' + 1
        local lag = max(1, `h' + 1)

        capture xtscc dy_`h' onset_all `controls' i.year ///
            if sample == 1 & (year + `h') <= gdp_last_actual, fe lag(`lag')

        if _rc {
            di as error "  outturn-only regression failed at h=`hd' (rc=" _rc ")"
            continue
        }

        local bb = _b[onset_all]
        local ss = _se[onset_all]
        local nn = e(N)
        _nepcount onset_all, outcome(dy_`h') controls(`controls')
        local nep = r(n)
        _pval `bb' `ss'
        local pp = r(p)

        eststo t1r_h`h'
        estadd scalar nep = `nep'

        di "h=" `hd' "  " %7.3f `bb' "  " %6.3f `ss' "  " %5.0f `nn' ///
           "     " %3.0f `nep' "     " %5.3f `pp'
    }

    capture esttab t1r_h0 t1r_h1 t1r_h2 t1r_h3 t1r_h4 ///
        using "$tabs/table1r_output_all_outturns.rtf", replace ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_all) coeflabel(onset_all "Spread-crisis onset") ///
        mtitles("h=1" "h=2" "h=3" "h=4" "h=5") nonumber ///
        stats(nep N N_g, labels("Episodes (onsets in estimation sample)" "Observations" "Countries") fmt(0 0 0)) ///
        title("Table 1R. Output cost of sovereign spread crises — outturns only (robustness)") ///
        addnotes("As Table 1, but each horizon drops observations whose outcome year t+h falls after the last WEO actual observation for that country (LATEST_ACTUAL_ANNUAL_DATA)." ///
                 "This removes IMF projections from the dependent variable; the sample therefore shrinks at longer horizons." ///
                 "Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01.")

    if _rc == 608 di as error "  ** table1r_output_all_outturns.rtf is OPEN IN WORD — close it and re-run."
    else if _rc di as error "  ** Table 1R: esttab failed (rc=" _rc ")"
    else di as result "Table 1R saved: $tabs/table1r_output_all_outturns.rtf"
}

* ══════════════════════════════════════════════════════════════════════════
* ROBUSTNESS: CONTINUATION YEARS KEPT AS CONTROLS (matches Asonuma et al.'s
* own onset design; see METHODOLOGY.md section 5)
*
* This project's baseline `sample' requires continuation==0: a country three
* years into an ongoing episode is dropped from the regression entirely, not
* treated and not a control. Reading Asonuma et al.'s own replication script
* (the Figure 3 IRF code) shows their design does NOT do this — dum1/dum2/dum3
* are simply 0 on a continuation-year row, so that row stays in the regression
* as an implicit CONTROL, alongside genuinely tranquil years, no matter how
* deep into an unresolved restructuring it is.
*
* `sample_flow' (built in 18_transforms.do for the flow tier) already IS the
* exact restriction this needs -- `sample' minus the continuation==0
* requirement -- so no new column is built here; it is reused purely as a
* SAMPLE restriction under this file's own onset_all treatment, which is a
* legitimate, minimal reuse, not a repurposing of what it means.
*
* Read this as a direct test of the concern raised discussing Asonuma's own
* design: if in-progress restructuring years are themselves still output-
* depressed, folding them into the control pool (rather than dropping them,
* as the baseline does) should bias onset_all's coefficient TOWARD ZERO
* relative to the baseline above. Reported so that direction is checked
* against what actually happens, not assumed.
* ══════════════════════════════════════════════════════════════════════════
capture confirm variable sample_flow
if _rc {
    di as error _n "  ** sample_flow not in the panel — Asonuma-sample robustness skipped."
    di as error "     (re-run 18_transforms.do first.)"
}
else {
    di as result _n "=== ROBUSTNESS: CONTINUATION YEARS KEPT AS CONTROLS (Asonuma et al.'s own sample) ==="
    di as result "h   beta     SE      N     episodes   p"

    forvalues h = 0/4 {
        local hd  = `h' + 1
        local lag = max(1, `h' + 1)

        capture xtscc dy_`h' onset_all `controls' i.year ///
            if sample_flow == 1, fe lag(`lag')

        if _rc {
            di as error "  Asonuma-sample regression failed at h=`hd' (rc=" _rc ")"
            continue
        }

        local bb = _b[onset_all]
        local ss = _se[onset_all]
        local nn = e(N)
        _nepcount onset_all, outcome(dy_`h') controls(`controls')
        local nep = r(n)
        _pval `bb' `ss'
        local pp = r(p)

        eststo t1a_h`h'
        estadd scalar nep = `nep'

        di "h=" `hd' "  " %7.3f `bb' "  " %6.3f `ss' "  " %5.0f `nn' ///
           "     " %3.0f `nep' "     " %5.3f `pp'
    }

    di as result _n "  Compare against the baseline (sample==1) row printed under"
    di as result "  '=== MAIN LP RESULTS (ALL CRISES) ===' above -- same horizons, same"
    di as result "  controls/FE/DK lag, only the sample restriction differs."

    capture esttab t1a_h0 t1a_h1 t1a_h2 t1a_h3 t1a_h4 ///
        using "$tabs/table1a_output_all_asonumasample.rtf", replace ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_all) coeflabel(onset_all "Spread-crisis onset") ///
        mtitles("h=1" "h=2" "h=3" "h=4" "h=5") nonumber ///
        stats(nep N N_g, labels("Episodes (onsets in estimation sample)" "Observations" "Countries") fmt(0 0 0)) ///
        title("Table 1A. Output cost of sovereign spread crises — continuation years kept as controls (robustness)") ///
        addnotes("As Table 1, but the sample restriction is sample_flow==1 rather than sample==1: continuation years of an ongoing" ///
                 "episode are kept in the regression as controls rather than dropped, matching Asonuma et al.'s own onset design" ///
                 "(their dum1/dum2/dum3 are 0, not excluded, on those rows). Driscoll-Kraay standard errors in parentheses." ///
                 "* p<0.10, ** p<0.05, *** p<0.01.")

    if _rc == 608 di as error "  ** table1a_output_all_asonumasample.rtf is OPEN IN WORD — close it and re-run."
    else if _rc di as error "  ** Table 1A: esttab failed (rc=" _rc ")"
    else di as result "Table 1A saved: $tabs/table1a_output_all_asonumasample.rtf"
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
    mtitles("h=1" "h=2" "h=3" "h=4" "h=5") nonumber ///
    stats(nep N N_g, labels("Episodes (onsets in estimation sample)" "Observations" "Countries") fmt(0 0 0)) ///
    title("Table 1. Output cost of sovereign spread crises (all episodes)") ///
    addnotes("Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h." ///
             "Jorda (2005) local projections. Country and year fixed effects; continuation years excluded." ///
             "Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01." ///
             "Episodes counts onsets surviving listwise deletion on the control set, i.e. those actually identifying the coefficient." ///
             "Confidence intervals in the figures use t critical values on the estimator's residual degrees of freedom, not the normal approximation.")

if _rc == 608 di as error "  ** table1_output_all.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc di as error "  ** Table 1: esttab failed (rc=" _rc ")"
else di as result "Table 1 saved: $tabs/table1_output_all.rtf"

* ────────────────────────────────────────────────────────────────────────
* BUILD IRF DATASET FOR GRAPHING
* ────────────────────────────────────────────────────────────────────────

clear
local nhor = 7
set obs `nhor'

gen horizon = _n - 2     // -1, 0, 1, 2, 3, 4, 5 (Asonuma numbering: crisis year = 1)

svmat b_all,    names(b)
svmat se_all,   names(se)
svmat lo90_all, names(lo90)
svmat hi90_all, names(hi90)
svmat lo95_all, names(lo95)
svmat hi95_all, names(hi95)

rename b1    b
rename se1   se
rename lo901 lo90
rename hi901 hi90
rename lo951 lo95
rename hi951 hi95

gen series = "all"
label var horizon "Horizon (years after onset)"
label var b       "Point estimate (pp)"
label var se      "Driscoll-Kraay standard error"
label var lo90    "90% CI lower"
label var hi90    "90% CI upper"
label var lo95    "95% CI lower"
label var hi95    "95% CI upper"

save "$clean/irf_all.dta", replace
di as result _n "IRF data saved: $clean/irf_all.dta"
