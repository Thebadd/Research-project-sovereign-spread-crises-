/*===========================================================================
  08D_FIRST_STAGE_TABLE_MINDUR.DO
  First-stage probit table under the MINIMUM-DURATION episode definition
  (18b_mindur_variant.do) -- mirrors 08c_first_stage_table.do's own adopted
  probit EXACTLY (same baseline controls X = $ctrl_core, same predictors
  Z2 = l_fedfunds l_contagion_dist_def years_since_def_onset, same pooled/
  no-country-FE design, same clustered SEs), swapping onset_nd/onset_def
  for onset_nd_mindur/onset_def_mindur -- i.e., re-fitting the two columns
  (Non-default vs tranquil, Default-linked vs tranquil) with the 18
  single-year-spike non-default onsets reclassified as tranquil (folded
  into the control pool) rather than treated. See 08c_first_stage_table.do's
  own header for the full argument behind every design choice here
  (default-linked-specific predictors, pooled/no-FE probit, AUROC diagnostic
  convention) -- unchanged, only the treatment columns differ.

  SCOPE: reports the SAME headline diagnostics as 08c (chi2(predictors),
  AUROC with/without predictors, roccomp's formal test) so the mindur
  columns' classification power can be read directly against the baseline
  table -- does NOT repeat 08c's own one-time separation-source and
  country-FE diagnostics (those investigated project-wide questions
  already resolved there, not something to re-derive per robustness
  variant).

  STANDALONE: reads $clean/panel_lp_mindur.dta (built by
  18b_mindur_variant.do). Does not modify 08c_first_stage_table.do or any
  baseline column.

  Output: $tabs/table_first_stage_mindur.rtf.
  Run AFTER 18b_mindur_variant.do. Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp_mindur.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

local X    $ctrl_core
local Z2   l_fedfunds l_contagion_dist_def years_since_def_onset

eststo clear

* ── Helper: identical to 08c_first_stage_table.do's own _fscol ────────────
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

di as result _n "=== FIRST-STAGE PROBIT (MINIMUM-DURATION ROBUSTNESS) ==="
_fscol fs_nd_md  "onset_nd_mindur"  "sample==1 & onset_def_mindur==0"  "`X'" "`Z2'"
_fscol fs_def_md "onset_def_mindur" "sample==1 & onset_nd_mindur==0"   "`X'" "`Z2'"

* ── FORMAL TEST OF WHETHER THE TWO AUROCs ACTUALLY DIFFER (same construction
* as 08c's own roccomp block) ──────────────────────────────────────────────
quietly roccomp onset_nd_mindur _pctrl_fs_nd_md _pfull_fs_nd_md if !missing(_pctrl_fs_nd_md,_pfull_fs_nd_md)
local rocchi2_nd = r(chi2)
local rocp_nd    = r(p)
quietly roccomp onset_def_mindur _pctrl_fs_def_md _pfull_fs_def_md if !missing(_pctrl_fs_def_md,_pfull_fs_def_md)
local rocchi2_def = r(chi2)
local rocp_def    = r(p)
capture drop _pctrl_fs_nd_md _pfull_fs_nd_md _pctrl_fs_def_md _pfull_fs_def_md

di as result _n "=== FIRST-STAGE PROBIT DIAGNOSTICS (predictors jointly, mindur) ==="
di as result "col              chi2(pred)   p        AUROC(ctrl only)  AUROC(+pred)  delta   roccomp chi2   p"
foreach c in fs_nd_md fs_def_md {
    quietly estimates restore `c'
    local dlt = e(auroc) - e(aurocctrl)
    local dltsign = cond(`dlt' >= 0, "+", "")
    local sfx = subinstr(subinstr("`c'", "fs_", "", .), "_md", "", .)
    di as result %-16s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " ///
                 %8.3f e(aurocctrl) "        " %6.3f e(auroc) "       " "`dltsign'" %6.3f `dlt' ///
                 "     " %6.2f `rocchi2_`sfx'' "        " %5.3f `rocp_`sfx''
}
local worstp = max(`rocp_nd', `rocp_def')
di as result _n "      roccomp's formal test: nd chi2(1)=" %5.2f `rocchi2_nd' ", p=" %5.3f `rocp_nd' ///
    "   def chi2(1)=" %5.2f `rocchi2_def' ", p=" %5.3f `rocp_def' "."
di as result "      Read alongside 08c_first_stage_table.do's own baseline columns:"
di as result "      if the mindur columns' AUROC/roccomp read similarly (predictors still add"
di as result "      power, classification still distinguishable from controls-only), the"
di as result "      propensity model is not sensitive to the single-year-spike reclassification."

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — same Table 1 style as 08c_first_stage_table.do
* ══════════════════════════════════════════════════════════════════════════
capture esttab fs_nd_md fs_def_md using "$tabs/table_first_stage_mindur.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("Non-default (mindur)" "Default-linked (mindur)") ///
    order(l_fedfunds l_contagion_dist_def years_since_def_onset ///
          l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2) ///
    coeflabel(l_fedfunds "US federal funds rate" ///
              l_contagion_dist_def "Contagion, based on default-linked crisis" ///
              years_since_def_onset "Years since last default-linked onset" ///
              l1_gdpg "GDP growth" ///
              l_debt "Public debt-to-GDP ratio" ///
              l_banking_crisis "Banking crisis dummy" ///
              l_govexp "Government expenditure-to-GDP ratio" ///
              l_open "Trade openness" ///
              l_credit_bank "Bank credit-to-GDP ratio" ///
              l_lninfl "Inflation" ///
              exchange2 "Nominal exchange-rate") ///
    refcat(l_fedfunds "Predictors" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp aurocctrl auroc N, ///
          labels("Chi-squared for predictors" " p-value of Chi-squared" "AUROC, controls only" "AUROC, with predictors" "Observations") ///
          fmt(2 3 3 3 0)) ///
    title("Table 1 (minimum-duration robustness). First-stage probit: predicting the start of a spread crisis") ///
    addnotes("Dependent variable: dummy = 1 in the onset year of the indicated crisis type UNDER THE MINIMUM-DURATION DEFINITION (single-year" ///
             "crisis_any spikes reclassified tranquil); each type predicted vs tranquil years (the rival onset type is dropped)." ///
             "Pooled probit, no country fixed effects (matching 08c_first_stage_table.do's own design). Robust standard errors clustered by country." ///
             "Predictors are excluded from the LP/AIPW outcome equation. Compare directly against 08c_first_stage_table.do's own" ///
             "table_first_stage.rtf: same controls/predictors/estimator, only the treatment definition differs." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage_mindur.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc  di as error "  ** Table (first stage, mindur): esttab failed (rc=" _rc ")"
else di as result "First-stage table (mindur) saved: $tabs/table_first_stage_mindur.rtf"

di as result _n "08d_first_stage_table_mindur.do complete."
