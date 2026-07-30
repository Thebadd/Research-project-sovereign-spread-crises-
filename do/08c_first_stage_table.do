/*===========================================================================
  08C_FIRST_STAGE_TABLE.DO
  First-stage probit table, styled on Asonuma et al. (2024) Table 1
  ("Predicting the start of debt restructurings, probit").

  Columns = type of onset, each predicted vs TRANQUIL (the rival onset type is
  dropped so the control group is clean non-crisis years — exactly as in the
  paper, and exactly the first stages that feed the AIPW two-line Act 2 in 08b):
    (1) All onsets    : onset_all vs tranquil        (sample==1)
    (2) Non-default   : onset_nd  vs tranquil        (sample==1 & onset_def==0)
    (3) Default-linked: onset_def vs tranquil        (sample==1 & onset_nd==0)

  Rows are grouped, as in their Table 1, into:
    PREDICTORS (excluded from the LP/AIPW outcome eq.):
        fed funds rate (global push) + regional contagion + past onsets.
        Col 1 (all onsets) uses past_onsets; the resolution columns use
        past_def_onsets (count of past DEFAULT-linked onsets).
    BASELINE CONTROLS (also in the outcome eq.):
        l1_gdpg  l2_gdpg  debt  ca  infl  imf

  Diagnostic rows (their Table 1 bottom block):
    Chi2 (predictors)  — joint Wald test that the 4 predictors are all zero
    p-value            — of that test (small => predictors jointly matter)
    Area under ROC     — lroc on the full model
    Observations

  Pooled probit, no country FE (rare-event propensity; see PAPER_FRAMING §5.C).
  Output: $tabs/table_first_stage.rtf. Run AFTER 01e_predictors.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear

* Baseline controls X (all columns) and predictors: Z1 for Act 1 (all onsets),
* Z2 for the resolution columns (proneness = past DEFAULT-linked onsets).
local X    l1_gdpg l2_gdpg debt ca infl imf
local Z1   fedfunds l_reg_crisis_share past_onsets
local Z2   fedfunds l_reg_crisis_share past_def_onsets

eststo clear

* ── Helper: run one column, add the χ²(predictors)/p/AUROC diagnostics ────────
* A program cannot see the caller's locals, so X and Z are passed in explicitly.
* args: estimate-name  DVvar  "if-condition"  "X list"  "Z list"
capture program drop _fscol
program define _fscol
    args nm dv ifcond xlist zlist
    quietly probit `dv' `xlist' `zlist' if `ifcond', vce(cluster cid)
    eststo `nm'
    quietly test `zlist'
    estadd scalar chi2p = r(chi2)
    estadd scalar pp    = r(p)
    quietly lroc, nograph
    estadd scalar auroc = r(area)
end

_fscol fs_all "onset_all" "sample==1"                    "`X'" "`Z1'"
_fscol fs_nd  "onset_nd"  "sample==1 & onset_def==0"      "`X'" "`Z2'"
_fscol fs_def "onset_def" "sample==1 & onset_nd==0"       "`X'" "`Z2'"

* ── Console echo of the diagnostics ──────────────────────────────────────────
di as result _n "=== FIRST-STAGE PROBIT DIAGNOSTICS (predictors jointly) ==="
di as result "col            chi2(pred)   p        AUROC"
foreach c in fs_all fs_nd fs_def {
    quietly estimates restore `c'
    di as result %-14s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " %6.3f e(auroc)
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — Table 1 style (Predictors / Baseline controls blocks + diags)
* ══════════════════════════════════════════════════════════════════════════
capture esttab fs_all fs_nd fs_def using "$tabs/table_first_stage.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("All onsets" "Non-default" "Default-linked") ///
    order(fedfunds l_reg_crisis_share past_onsets past_def_onsets ///
          l1_gdpg l2_gdpg debt ca infl imf) ///
    coeflabel(fedfunds "US fed funds rate" ///
              l_reg_crisis_share "Regional contagion (t-1)" ///
              past_onsets "Past onsets (any type)" ///
              past_def_onsets "Past default-linked onsets" ///
              l1_gdpg "GDP growth (t-1)" ///
              l2_gdpg "GDP growth (t-2)" ///
              debt "Public debt / GDP" ///
              ca "Current account / GDP" ///
              infl "CPI inflation" ///
              imf "IMF program") ///
    refcat(fedfunds "Predictors" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp auroc N, ///
          labels("Chi-squared (predictors)" "  p-value" "Area under ROC curve" "Observations") ///
          fmt(2 3 3 0)) ///
    title("Table 1. First-stage probit: predicting the start of a spread crisis") ///
    addnotes("Dependent variable: dummy = 1 in the onset year of the indicated crisis type; each type predicted vs tranquil years (the rival onset type is dropped)." ///
             "Pooled probit, no country fixed effects. Robust standard errors clustered by country in parentheses." ///
             "Predictors are excluded from the LP/AIPW outcome equation (ust10y, vix absorbed by year FE; contagion and past onsets are country-varying and omitted from the second stage)." ///
             "Chi-squared (predictors) is the joint Wald test that all four predictors are zero. Area under ROC is for the full model." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc  di as error "  ** Table (first stage): esttab failed (rc=" _rc ")"
else di as result "First-stage table saved: $tabs/table_first_stage.rtf"

di as result _n "08c_first_stage_table.do complete."
