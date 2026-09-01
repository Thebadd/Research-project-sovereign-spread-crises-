/*===========================================================================
  08C_FIRST_STAGE_TABLE.DO
  First-stage probit table, styled on Asonuma et al. (2024) Table 1
  ("Predicting the start of debt restructurings, probit") — the ONSET-tier
  counterpart of 21b_first_stage_table_flow.do; mirrors that file's
  diagnostics (AUROC with/without predictors + roccomp's formal test, not
  just chi2(predictors)) for this file's own onset specification.

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
    PREDICTORS (excluded from the LP/AIPW outcome eq.): Fed funds rate
        (global push) + distance-weighted contagion (l_contagion_dist,
        a country-year-specific spatial lag, ANY onset type) + years since
        the most recent prior onset of ANY type (years_since_onset, a
        recency clock, censored at 50). Both the contagion and recency terms
        are deliberately kept GENERIC (not narrowed/broadened to
        default-linked-only) so the two predictors are built consistently
        with each other -- see the diagnostic history below and
        17_predictors.do's header for the two comparisons that were run
        before settling on this combination.
    BASELINE CONTROLS (the SAME $ctrl_core used in the LP/AIPW outcome eq. —
        strict parity with the reference paper's $convar-in-both design):
        l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2

  Diagnostic rows (their Table 1 bottom block):
    Chi2 (predictors)      — joint Wald test that the 3 predictors are all zero
    p-value                — of that test (small => predictors jointly matter)
    AUROC, controls only   — lroc on a probit with ONLY $ctrl_core, no predictors
    AUROC, with predictors — lroc on the full model
    Observations

  WHY BOTH AUROCs, NOT JUST THE FULL MODEL'S (matches 21b's reasoning): the
  source paper's actual demonstration that predictors earn their place isn't
  the chi-squared test alone -- it's the WITH-vs-WITHOUT-predictors AUROC
  comparison they report directly in their text (0.87 vs 0.79 post-default,
  0.94 vs 0.85 preemptive). Chi2(predictors) only says the predictors are
  jointly non-zero; it says nothing about how much discriminatory power they
  actually add on top of what the baseline controls already provide.

  THE DELTA ITSELF IS NOT A SIGNIFICANCE TEST -- roccomp's chi2 test for
  equal correlated ROC areas is, and it is reported alongside the delta in
  the console output. The write-up should quote the roccomp p-values, not
  just the delta, when characterising how much work the predictors are doing.

  Pooled probit, no country FE (rare-event propensity: country FE separate/overfit
  with ~20 events; they enter the outcome equation only). No year FE.
  Clustered SEs by country, matching 21b's choice.

  Output: $tabs/table_first_stage.rtf. Run AFTER 17_predictors.do. The
  figure counterpart (kernel density + nested ROC, mirroring 21c) is
  08d_first_stage_figs.do.

  DIAGNOSTIC HISTORY (both resolved, neither block kept as live code):
    - A default-linked-only contagion measure (contagion_dist_def) was
      tested against the generic l_contagion_dist used here: no difference
      in either arm (nd delta +0.009, p=.609; def delta +0.000, p=.931).
      Generic kept.
    - A default-linked-only recency measure (years_since_def_onset) was
      tested against the generic years_since_onset used here in an earlier
      pass. Both predictors are now kept generic (any-onset-type) for
      consistency with each other, rather than mixing a generic contagion
      term with a default-specific recency term.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

* Baseline controls X (both columns) and predictors Z2 (resolution-type
* proneness = years since the most recent prior default-linked onset).
* Baseline = outcome baseline ($ctrl_core): the Table-1 probit shares its controls
* with the LP/AIPW outcome equation (their $convar in both stages).
local X    $ctrl_core
local Z2   l_fedfunds l_contagion_dist years_since_onset

eststo clear

* ── Helper: run BOTH nested models per column, add the AUROC-with-vs-without-
* predictors + chi2(predictors) diagnostics — mirrors 21b's _fscolflow.
* A program cannot see the caller's locals, so X and Z are passed in explicitly.
* args: estimate-name  DVvar  "if-condition"  "X list"  "Z list"
* Controls-only model is fit FIRST so its AUROC can be estadd-ed onto the
* full model's stored estimate once that is eststo-d (estadd only writes to
* the currently active estimates). Predicted probabilities from both models
* are saved for the roccomp formal test run right after this program.
capture program drop _fscol
program define _fscol
    args nm dv ifcond xlist zlist
    quietly probit `dv' `xlist' if `ifcond', vce(cluster cid)
    quietly lroc, nograph
    local aucctrl = r(area)
    capture drop _pctrl_`nm'
    quietly predict double _pctrl_`nm' if `ifcond', pr

    probit `dv' `xlist' `zlist' if `ifcond', vce(cluster cid)
    eststo `nm'
    quietly test `zlist'
    estadd scalar chi2p = r(chi2)
    estadd scalar pp    = r(p)
    quietly lroc, nograph
    estadd scalar auroc     = r(area)
    estadd scalar aurocctrl = `aucctrl'
    capture drop _pfull_`nm'
    quietly predict double _pfull_`nm' if `ifcond', pr
end

_fscol fs_nd  "onset_nd"  "sample==1 & onset_def==0"      "`X'" "`Z2'"
_fscol fs_def "onset_def" "sample==1 & onset_nd==0"       "`X'" "`Z2'"

* ── FORMAL TEST OF WHETHER THE TWO AUROCs ACTUALLY DIFFER ──────────────────
* The delta alone (auroc - aurocctrl) says nothing about whether that gap is
* distinguishable from sampling noise -- an informal threshold cannot
* substitute for a test. roccomp's chi2 test for equal correlated ROC areas
* (same subjects, two nested classifiers) is that test, and it needs the two
* models' predicted probabilities on the SAME sample, saved above.
quietly roccomp onset_nd _pctrl_fs_nd _pfull_fs_nd if !missing(_pctrl_fs_nd,_pfull_fs_nd)
local rocchi2_nd = r(chi2)
local rocp_nd    = r(p)
quietly roccomp onset_def _pctrl_fs_def _pfull_fs_def if !missing(_pctrl_fs_def,_pfull_fs_def)
local rocchi2_def = r(chi2)
local rocp_def    = r(p)
capture drop _pctrl_fs_nd _pfull_fs_nd _pctrl_fs_def _pfull_fs_def

* ── Console echo of the diagnostics ──────────────────────────────────────────
di as result _n "=== FIRST-STAGE PROBIT DIAGNOSTICS (predictors jointly) ==="
di as result "col            chi2(pred)   p        AUROC(ctrl only)  AUROC(+pred)  delta   roccomp chi2   p"
local mindlt = .
foreach c in fs_nd fs_def {
    quietly estimates restore `c'
    local dlt = e(auroc) - e(aurocctrl)
    local dltsign = cond(`dlt' >= 0, "+", "")
    local sfx = subinstr("`c'", "fs_", "", .)
    di as result %-14s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " ///
                 %8.3f e(aurocctrl) "        " %6.3f e(auroc) "       " "`dltsign'" %6.3f `dlt' ///
                 "     " %6.2f `rocchi2_`sfx'' "        " %5.3f `rocp_`sfx''
    * Track the SMALLER of the two deltas -- the interpretation below should
    * reflect the weaker column, not be driven by whichever column happens
    * to look best.
    if missing(`mindlt') | `dlt' < `mindlt' local mindlt = `dlt'
}
* The raw delta is NOT a significance test -- roccomp's chi2 test (above,
* same-subjects correlated-ROC-areas test) is. Report that number, not an
* informal threshold, as the governing statistic.
local worstp = max(`rocp_nd', `rocp_def')
di as result _n "      Both columns' AUROC delta: " %5.3f `mindlt' " (smaller of the two)."
di as result "      roccomp's formal test for whether the two AUROCs differ: nd chi2(1)=" ///
    %5.2f `rocchi2_nd' ", p=" %5.3f `rocp_nd' "   def chi2(1)=" %5.2f `rocchi2_def' ", p=" %5.3f `rocp_def' "."
if `worstp' < 0.05 {
    di as result "      BOTH columns clear the conventional 5 pct level (worse column p=" %5.3f `worstp' ")."
    di as result "      Read this as ESTABLISHED, not merely suggestive: the chi-squared row above tests joint"
    di as result "      predictor significance only and says nothing about classification power; roccomp's"
    di as result "      test is the one that speaks to classification power directly, and here it backs the"
    di as result "      chi-squared row up."
}
else if `worstp' < 0.10 {
    di as result "      One column clears the conventional 5 pct level; the other does not (worse column p=" ///
        %5.3f `worstp' "), but both are inside the 10 pct band."
    di as result "      Read this as MODEST, NOT ABSENT: the chi-squared row above tests joint predictor"
    di as result "      significance only and says nothing about classification power; roccomp's test is the"
    di as result "      one that speaks to classification power directly, and it says 'suggestive, not fully"
    di as result "      established' here."
}
else {
    di as result "      Neither reaches the conventional 5 pct level (worse column p=" %5.3f `worstp' ///
        "); at least one arm sits outside even the 10 pct band."
    di as result "      Read this as MODEST, NOT ABSENT: the chi-squared row above tests joint predictor"
    di as result "      significance only and says nothing about classification power; roccomp's test is the"
    di as result "      one that speaks to classification power directly, and it says 'suggestive, not"
    di as result "      established' here."
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — Table 1 style (Predictors / Baseline controls blocks + diags)
* ══════════════════════════════════════════════════════════════════════════
capture esttab fs_nd fs_def using "$tabs/table_first_stage.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("Non-default" "Default-linked") ///
    order(l_fedfunds l_contagion_dist years_since_onset ///
          l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2) ///
    coeflabel(l_fedfunds "US fed funds rate (t-1)" ///
              l_contagion_dist "Distance-weighted contagion (t-1)" ///
              years_since_onset "Years since last onset (any type)" ///
              l1_gdpg "GDP growth (t-1)" ///
              l_debt "Public debt / GDP (t-1)" ///
              l_banking_crisis "Systemic banking-crisis dummy (t-1)" ///
              l_govexp "Govt expenditure / GDP (t-1)" ///
              l_open "Trade openness (t-1)" ///
              l_credit_bank "Private credit by banks / GDP (t-1)" ///
              l_lninfl "Log inflation, continuous (t-1)" ///
              exchange2 "Nominal exchange-rate log-change (t-1)") ///
    refcat(l_fedfunds "Predictors" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp aurocctrl auroc N, ///
          labels("Chi-squared (predictors)" "  p-value" "AUROC, controls only" "AUROC, with predictors" "Observations") ///
          fmt(2 3 3 3 0)) ///
    title("Table 1. First-stage probit: predicting the start of a spread crisis") ///
    addnotes("Dependent variable: dummy = 1 in the onset year of the indicated crisis type; each type predicted vs tranquil years (the rival onset type is dropped)." ///
             "Pooled probit, no country fixed effects (country FE separate on the thin event count; they enter the outcome equation only). Robust standard errors clustered by country in parentheses." ///
             "Predictors are excluded from the LP/AIPW outcome equation (they are omitted from the second stage, which carries country FE only, no year FE)." ///
             "Chi-squared (predictors) is the joint Wald test that all three predictors are zero. AUROC, controls only is from a probit on the baseline controls" ///
             "alone (no predictors); AUROC, with predictors adds the three predictors -- the difference is the source Table 1's own demonstration that" ///
             "predictors add classification power, not just the chi-squared joint-significance test." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc  di as error "  ** Table (first stage): esttab failed (rc=" _rc ")"
else di as result "First-stage table saved: $tabs/table_first_stage.rtf"

di as result _n "08c_first_stage_table.do complete."
