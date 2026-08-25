/*===========================================================================
  21B_FIRST_STAGE_TABLE_FLOW.DO
  First-stage probit table, styled on Asonuma et al. (2024) Table 1
  ("Predicting the start of debt restructurings, probit") — the FLOW file's
  analogue of 08c_first_stage_table.do.

  WHY THIS FILE EXISTS
  --------------------
  08c_first_stage_table.do documents the propensity model behind 08b_aipw.do
  (onset coding). 21_aipw_flow.do's Eq. (2) now fits on the SAME kind of
  restricted sample — tranquil and onset rows only, continuation years
  excluded (`if `touse' & continuation==0' inside _aipw, adopted per the
  estimator change documented in 21_aipw_flow.do's header, "SECTIONS 1d/1e"
  and "ESTIMATOR CHANGE, SECTION 1f ONWARD"). This file is that model, styled
  the same way 08c styles its onset probit: PREDICTING THE START of a spread
  crisis, non-default and default-linked, never continuation years — matching
  the reference IMF working paper's own Table 1, whose main specification
  bounds the probit to "the year of the start of the restructuring and in the
  previous year" (their Section 4.1), with a continuation-inclusive version
  relegated to an appendix robustness check rather than the headline table.

  SAMPLE: sample_flow==1 & continuation==0 -- tranquil country-years plus
  ONSET rows only. This is deliberately NOT sample_flow==1 alone: the whole
  point of the Section 1f estimator change was to keep continuation years out
  of the probit's fitting sample while still scoring them (Eq. 1/Eq. 3 keep
  the full flow-coded treatment) -- this table documents exactly the sample
  the actual estimator fits on, not the broader sample it later extrapolates
  to.

  COLUMNS: non-default and default-linked, each predicted vs TRANQUIL (the
  rival type dropped), exactly as 08c does it and as the reference paper's
  own sample_for* design does.

  Rows, grouped as in the source Table 1:
    PREDICTORS (excluded from the LP/AIPW outcome eq. -- $ctrl_flow/$com):
        Fed funds rate (global push) + distance-weighted contagion
        (l_contagion_dist, CEPII great-circle -- see 17_predictors.do and
        21_aipw_flow.do's "SECOND PREDICTOR CHANGE") + years since the
        country's most recent PRIOR default-linked onset. Both columns share
        these predictors, matching 21_aipw_flow.do's cz_recency (used
        identically for both nd and def there).
    PREDICTOR CHANGE FROM past_def_onsets: 24_aipw_channels_flow.do's Section
        1a diagnostic (credit, h=1) found the def arm's propensity model
        severely weight-concentrated, with past_def_onsets (a running COUNT
        that never resets) as the leading suspect -- for a serial defaulter it
        behaves close to a permanent country identifier. Section 1b there
        tested years_since_def_onset as a replacement: individually
        significant (z=-2.84, p=.005), economically sensible sign, and not
        circular for the same reason the continuation==0 restriction above is
        not. It did NOT meaningfully fix the weight concentration (98.9% ->
        98.6%, essentially unchanged) -- adopted for being a better,
        independently-motivated predictor, not as a variance fix. See
        21_aipw_flow.do's header, "PREDICTOR CHANGE," for the full note.
    BASELINE CONTROLS: $ctrl_core_flowbase -- the SAME row-dated ADOPTED set
        21_aipw_flow.do's `cx_active' uses for Eq. (2) (NOT $ctrl_flow/epc_*,
        which the header of 21_aipw_flow.do -- Section 1(a) -- demonstrates
        cannot be used in the probit: epc_X differs from X only on
        continuation rows, so it would mechanically encode treatment).

  Diagnostic rows (source Table 1's bottom block):
    Chi2 (predictors)      — joint Wald test that the 3 predictors are all zero
    p-value                — of that test
    AUROC, controls only   — lroc on a probit with ONLY $ctrl_core_flowbase, no predictors
    AUROC, with predictors — lroc on the full model
    Observations

  WHY BOTH AUROCs, NOT JUST THE FULL MODEL'S: the source paper's actual
  demonstration that predictors earn their place isn't the chi-squared test
  alone -- it's the WITH-vs-WITHOUT-predictors AUROC comparison they report
  directly in their text (0.87 vs 0.79 post-default, 0.94 vs 0.85 preemptive).
  Chi2(predictors) only says the predictors are jointly non-zero; it says
  nothing about how much discriminatory power they actually add on top of
  what the baseline controls already provide. Reporting only the full model's
  AUROC would assert predictive power without showing the counterfactual.

  Pooled probit, no country FE (matches 21_aipw_flow.do's Eq. (2), which is
  "Pooled, no country FE" per that file's own header). Clustered SEs by
  country, matching 08c's choice.

  Output: $tabs/table_first_stage_flow.rtf. Run AFTER 21_aipw_flow.do (reads
  only $clean/panel_lp.dta, so it does not actually depend on 21 having run,
  but is numbered to sit beside it for the reader).
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"

foreach v in in_crisis_nd in_crisis_def sample_flow continuation years_since_def_onset l_contagion_dist {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel_lp.dta — re-run 18_transforms.do first."
        exit 111
    }
}

* CONTROL SET. Default 0 = the ADOPTED flow-tier baseline, $ctrl_core_flowbase
* (18_transforms.do's "ADOPTED FLOW-TIER CORE CONTROL SET": l_banking_crisis,
* tot_chg, l_lninfl in place of l_banking_duration, l_ca, l_hyperinfl). Set to
* 1 to test the bigger, still-exploratory alternate set instead (18_transforms.do's
* "EXPLORATORY ALTERNATE FLOW CONTROL SET").
local use_flowalt_ctrl 0
local X = cond(`use_flowalt_ctrl', "$ctrl_core_flowalt", "$ctrl_core_flowbase")
* l_contagion_dist (distance-weighted, CEPII great-circle) in place of
* l_reg_crisis_share (flat regional share) -- matches 21_aipw_flow.do's
* cz_recency, "SECOND PREDICTOR CHANGE".
local Z l_fedfunds l_contagion_dist years_since_def_onset

eststo clear

* ── Helper: run one column, add the χ²(predictors)/p/AUROC diagnostics ────────
* Sample is FIXED across both columns: sample_flow==1 & continuation==0, plus
* dropping the rival type -- the flow file's actual Eq. (2) fitting universe.
*
* AUROC WITH vs WITHOUT PREDICTORS -- the source Table 1's actual demonstration
* that predictors add classification power (their 0.87 vs 0.79, 0.94 vs 0.85),
* not just the chi2(predictors) joint-significance test. Controls-ONLY model is
* fit FIRST so its AUROC can be estadd-ed onto the full model's stored estimate
* once that is eststo-d (estadd only writes to the currently active estimates).
capture program drop _fscolflow
program define _fscolflow
    args nm dv ifcond xlist zlist
    quietly probit `dv' `xlist' if `ifcond', vce(cluster cid)
    quietly lroc, nograph
    local aucctrl = r(area)

    probit `dv' `xlist' `zlist' if `ifcond', vce(cluster cid)
    eststo `nm'
    quietly test `zlist'
    estadd scalar chi2p = r(chi2)
    estadd scalar pp    = r(p)
    quietly lroc, nograph
    estadd scalar auroc     = r(area)
    estadd scalar aurocctrl = `aucctrl'
end

_fscolflow ffs_nd  "in_crisis_nd"  "sample_flow==1 & continuation==0 & in_crisis_def==0" "`X'" "`Z'"
_fscolflow ffs_def "in_crisis_def" "sample_flow==1 & continuation==0 & in_crisis_nd==0"  "`X'" "`Z'"

* ── Console echo of the diagnostics ──────────────────────────────────────────
di as result _n "=== FLOW FIRST-STAGE PROBIT DIAGNOSTICS (predictors jointly) ==="
di as result "     sample = tranquil + onset rows only (continuation==0)"
di as result "col            chi2(pred)   p        AUROC(ctrl only)  AUROC(+pred)  delta"
local mindlt = .
foreach c in ffs_nd ffs_def {
    quietly estimates restore `c'
    local dlt = e(auroc) - e(aurocctrl)
    local dltsign = cond(`dlt' >= 0, "+", "")
    di as result %-14s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " ///
                 %8.3f e(aurocctrl) "        " %6.3f e(auroc) "       " "`dltsign'" %6.3f `dlt'
    * Track the SMALLER of the two deltas -- the interpretation below should
    * reflect the weaker column, not be driven by whichever column happens
    * to look best.
    if missing(`mindlt') | `dlt' < `mindlt' local mindlt = `dlt'
}
* Threshold is informal, not a hard rule: 0.03 is well below the reference
* paper's own with/without gaps (0.08-0.09), so anything at or above that is
* read as "the predictors are earning their place," consistent with the
* magnitude their own comparison (0.87 vs 0.79, 0.94 vs 0.85) demonstrates.
di as result _n "      Both columns' AUROC delta: " %5.3f `mindlt' " (smaller of the two) -- " ///
    cond(`mindlt' >= 0.03, ///
        "meaningfully positive, comparable in size to the reference paper's own", ///
        "close to zero, well below the reference paper's own")
di as result "      with/without-predictors gap (their 0.87 vs 0.79, 0.94 vs 0.85). " ///
    cond(`mindlt' >= 0.03, ///
        "The predictors are adding real classification power on top of the baseline", ///
        "The predictors are adding little classification power on top of the baseline")
di as result "      controls, not just asserting a role from theory (Chi-squared above tests joint"
di as result "      significance only, and says nothing about how much discriminatory power is added)."

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — Table 1 style (Predictors / Baseline controls blocks + diags)
* ══════════════════════════════════════════════════════════════════════════
capture esttab ffs_nd ffs_def using "$tabs/table_first_stage_flow.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("Non-default" "Default-linked") ///
    order(l_fedfunds l_contagion_dist years_since_def_onset ///
          l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl tot_chg) ///
    coeflabel(l_fedfunds "US fed funds rate (t-1)" ///
              l_contagion_dist "Distance-weighted contagion (t-1)" ///
              years_since_def_onset "Years since last default onset" ///
              l1_gdpg "GDP growth (t-1)" ///
              l_debt "Public debt / GDP (t-1)" ///
              l_banking_crisis "Banking crisis dummy (t-1)" ///
              l_govexp "Govt expenditure / GDP (t-1)" ///
              l_open "Trade openness (t-1)" ///
              l_credit_bank "Private credit by banks / GDP (t-1)" ///
              l_lninfl "Log gross inflation (t-1)" ///
              tot_chg "Terms-of-trade log-change (t-1)") ///
    refcat(l_fedfunds "Predictors" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp aurocctrl auroc N, ///
          labels("Chi-squared (predictors)" "  p-value" "AUROC, controls only" "AUROC, with predictors" "Observations") ///
          fmt(2 3 3 3 0)) ///
    title("Table 1f. First-stage probit for the flow AIPW: predicting the START of a spread crisis") ///
    addnotes("Dependent variable: dummy = 1 in the onset year of the indicated crisis type; each type predicted vs tranquil years, the rival type dropped." ///
             "Sample restricted to tranquil years and ONSET rows only (continuation==0) -- this is the actual fitting sample 21_aipw_flow.do's Eq. (2) now" ///
             "uses (see that file's header, ESTIMATOR CHANGE, SECTION 1f ONWARD), not the full flow-coded sample. Eq. (1) and Eq. (3) there still use the" ///
             "full flow-coded treatment; only the propensity model's fitting sample is bounded to the start year and tranquil years." ///
             "Pooled probit, no country fixed effects. Robust standard errors clustered by country in parentheses." ///
             "Predictors are excluded from the LP/AIPW outcome equation. Chi-squared (predictors) is the joint Wald test that the three predictors are zero." ///
             "AUROC, controls only is from a probit on the baseline controls alone (no predictors); AUROC, with predictors adds the three predictors -- the" ///
             "difference between the two is the source Table 1's own demonstration that predictors add classification power (their 0.87 vs 0.79, 0.94 vs 0.85)," ///
             "not just the chi-squared joint-significance test. * p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage_flow.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc  di as error "  ** Table (flow first stage): esttab failed (rc=" _rc ")"
else di as result "Flow first-stage table saved: $tabs/table_first_stage_flow.rtf"

di as result _n "21b_first_stage_table_flow.do complete."
