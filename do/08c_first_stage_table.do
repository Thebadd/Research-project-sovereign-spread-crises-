/*===========================================================================
  08C_FIRST_STAGE_TABLE.DO
  First-stage probit table, styled on Asonuma et al. (2024) Table 1
  ("Predicting the start of debt restructurings, probit").

  Columns = onset type BY RESOLUTION, each predicted vs TRANQUIL (the rival
  onset type is dropped so the control group is clean non-crisis years —
  exactly as in the paper, and exactly the first stages that feed the AIPW
  two-line Act 2 in 08b):
    (1) Non-default   : onset_nd  vs tranquil        (sample==1 & onset_def==0)
    (2) Default-linked: onset_def vs tranquil        (sample==1 & onset_nd==0)
  The pooled "All onsets vs tranquil" column is deliberately not reported
  here: every onset is either non-default or default-linked, so the pooled
  column carries no propensity-modeling information beyond the union of
  these two, and this table's whole purpose is to report the first stages
  the by-resolution-type AIPW design actually consumes.

  Rows are grouped, as in their Table 1, into:
    PREDICTORS (excluded from the LP/AIPW outcome eq.):
        Fed funds rate (global push) + regional contagion + past
        DEFAULT-linked onsets (past_def_onsets, a default-proneness
        predictor, not a generic crisis-proneness one).
    BASELINE CONTROLS (the SAME $ctrl_core used in the LP/AIPW outcome eq. —
        strict parity with the reference paper's $convar-in-both design):
        l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2

  Diagnostic rows (their Table 1 bottom block):
    Chi2 (predictors)  — joint Wald test that the 4 predictors are all zero
    p-value            — of that test (small => predictors jointly matter)
    Area under ROC     — lroc on the full model
    Observations

  Pooled probit, no country FE (rare-event propensity: country FE separate/overfit
  with ~20 events; they enter the outcome equation only). No year FE.
  Output: $tabs/table_first_stage.rtf. Run AFTER 01e_predictors.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear

* Baseline controls X (both columns) and predictors Z (proneness = past
* DEFAULT-linked onsets, shared by both resolution-type columns).
* Baseline = outcome baseline ($ctrl_core): the Table-1 probit shares its controls
* with the LP/AIPW outcome equation (their $convar in both stages).
local X    $ctrl_core
local Z2   l_fedfunds l_reg_crisis_share past_def_onsets

eststo clear

* ── Helper: run one column, add the χ²(predictors)/p/AUROC diagnostics ────────
* A program cannot see the caller's locals, so X and Z are passed in explicitly.
* args: estimate-name  DVvar  "if-condition"  "X list"  "Z list"
capture program drop _fscol
program define _fscol
    args nm dv ifcond xlist zlist
    probit `dv' `xlist' `zlist' if `ifcond', vce(cluster cid)
    eststo `nm'
    quietly test `zlist'
    estadd scalar chi2p = r(chi2)
    estadd scalar pp    = r(p)
    quietly lroc, nograph
    estadd scalar auroc = r(area)
end

_fscol fs_nd  "onset_nd"  "sample==1 & onset_def==0"      "`X'" "`Z2'"
_fscol fs_def "onset_def" "sample==1 & onset_nd==0"       "`X'" "`Z2'"

* ── Console echo of the diagnostics ──────────────────────────────────────────
di as result _n "=== FIRST-STAGE PROBIT DIAGNOSTICS (predictors jointly) ==="
di as result "col            chi2(pred)   p        AUROC"
foreach c in fs_nd fs_def {
    quietly estimates restore `c'
    di as result %-14s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " %6.3f e(auroc)
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — Table 1 style (Predictors / Baseline controls blocks + diags)
* ══════════════════════════════════════════════════════════════════════════
capture esttab fs_nd fs_def using "$tabs/table_first_stage.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("Non-default" "Default-linked") ///
    order(l_fedfunds l_reg_crisis_share past_def_onsets ///
          l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2) ///
    coeflabel(l_fedfunds "US fed funds rate (t-1)" ///
              l_reg_crisis_share "Regional contagion (t-1)" ///
              past_def_onsets "Past default-linked onsets" ///
              l1_gdpg "GDP growth (t-1)" ///
              l_debt "Public debt / GDP (t-1)" ///
              l_banking_crisis "Systemic banking-crisis dummy (t-1)" ///
              l_govexp "Govt expenditure / GDP (t-1)" ///
              l_open "Trade openness (t-1)" ///
              l_credit_bank "Private credit by banks / GDP (t-1)" ///
              l_lninfl "Log inflation, continuous (t-1)" ///
              exchange2 "Nominal exchange-rate log-change (t-1)") ///
    refcat(l_fedfunds "Predictors" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp auroc N, ///
          labels("Chi-squared (predictors)" "  p-value" "Area under ROC curve" "Observations") ///
          fmt(2 3 3 0)) ///
    title("Table 1. First-stage probit: predicting the start of a spread crisis") ///
    addnotes("Dependent variable: dummy = 1 in the onset year of the indicated crisis type; each type predicted vs tranquil years (the rival onset type is dropped)." ///
             "Pooled probit, no country fixed effects (country FE separate on the thin event count; they enter the outcome equation only). Robust standard errors clustered by country in parentheses." ///
             "Predictors are excluded from the LP/AIPW outcome equation (they are omitted from the second stage, which carries country FE only, no year FE)." ///
             "Chi-squared (predictors) is the joint Wald test that all three predictors are zero. Area under ROC is for the full model." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc  di as error "  ** Table (first stage): esttab failed (rc=" _rc ")"
else di as result "First-stage table saved: $tabs/table_first_stage.rtf"

di as result _n "08c_first_stage_table.do complete."
