/*===========================================================================
  08E_FIRST_STAGE_TABLE_SPREADINTENSITY.DO
  ROBUSTNESS / STANDALONE — does the raw INTENSITY of the spread shock
  (not just pre-crisis risk factors) predict whether a spread crisis
  resolves as default-linked rather than non-default?

  MOTIVATION: 08c_first_stage_table.do's Table 1 probit predicts crisis
  type (onset_nd / onset_def, each vs tranquil years) using only:
    Predictors (Z2): l_fedfunds, l_contagion_dist_atdef, years_since_def_onset
    Baseline controls ($ctrl_core): l1_gdpg l_debt l_banking_crisis l_govexp
        l_open l_credit_bank l_lninfl exchange2
  None of these measure the MAGNITUDE of the spread spike itself -- they
  are all either external push factors, contagion, recency-of-default, or
  slow-moving domestic fundamentals. A supervisor comment raised exactly
  this gap: do more severe spread crises mechanically lead to default,
  and milder ones to non-default resolution, independently of the
  fundamentals already in the model? If so, part of this project's
  default-vs-non-default output-cost contrast could reflect crisis
  INTENSITY rather than the resolution type itself.

  spr_max (peak EMBIG spread level, bps, within the country-year) and
  spr_mean (annual average) already exist in the panel -- built in
  10_skeleton.do from EM_Spread_Crisis_DB_FINAL.xlsx and carried through
  panel_build.dta -> panel_lp.dta -- but are not used as predictors
  anywhere in the project. This file is a first look at whether they
  should be.

  THREE CHECKS, in order of increasing modeling assumptions:
    1. Descriptive: distribution of spr_max (and spr_mean) at onset,
       non-default vs default-linked, with a t-test and a Wilcoxon
       rank-sum test (t-test assumes near-normality with ~20/40 obs per
       group; Wilcoxon does not, so it is reported alongside as the
       less assumption-heavy check on the same question).
    2. Direct test: among onset observations only, does spr_max predict
       onset_def (a probit / marginal check), alone and with $ctrl_core?
    3. Augmented Table 1: re-run 08c's exact two-column onset-vs-tranquil
       design with spr_max added to the predictor set, to see whether it
       is itself significant and whether its inclusion changes the
       existing predictors' significance or the AUROC.

  NOT adopted anywhere -- this file does not touch 08c_first_stage_table.do,
  Table 1 of the paper, or any AIPW estimation file. Standalone, not wired
  into 00_master.do, matching this project's convention for exploratory
  robustness checks (27_png_debt_test.do, 28_aipw_extfin_nexus_split_nd.do).
  Run AFTER 17_predictors.do (same as 08c).

  Output:
    Console: descriptive comparison (means/medians, t-test, Wilcoxon),
             direct-probit result, and the diagnostics for the augmented
             Table 1 columns.
    $figs/fig_spreadintensity_by_crisistype.png — box plot of spr_max at
             onset, by resolution type.
    $tabs/table_first_stage_spreadintensity.rtf — augmented Table 1
             (spr_max added to Z2), a SEPARATE file from
             table_first_stage.rtf so the paper's own Table 1 is
             untouched.
===========================================================================*/

use "$clean/panel_lp.dta", clear
xtset cid year

capture confirm variable spr_max
if _rc {
    di as error "  ** spr_max not found in panel_lp.dta -- it may not survive 18_transforms.do's"
    di as error "     variable selection from panel_build.dta. Check 18_transforms.do's keep/drop"
    di as error "     list, or merge it back in from panel_build.dta before re-running this file."
    exit 111
}

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

* ══════════════════════════════════════════════════════════════════════════
* CHECK 0 — RECONCILE the 61-episode headline count (40 ND + 21 def) against
* the sample used below: onset_nd/onset_def fire on ALL identified onsets,
* but Check 1-3 also require sample==1 (this project's main analysis-sample
* flag) and non-missing spr_max. This block lists exactly which onset
* country-years are dropped by each restriction, and why, rather than
* leaving a silent N mismatch for the reader to puzzle over.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== CHECK 0: reconciling N against the paper's 61 identified onsets ==="
quietly count if onset_nd==1
di as result "      onset_nd==1, all rows (no sample/spr_max restriction): " r(N)
quietly count if onset_def==1
di as result "      onset_def==1, all rows (no sample/spr_max restriction): " r(N)

preserve
    quietly keep if onset_nd==1 | onset_def==1
    quietly count if sample!=1
    if r(N) > 0 {
        di as result _n "      Onsets excluded by sample!=1 (" r(N) "):"
        list country year onset_nd onset_def sample if sample!=1, noobs clean
    }
    quietly count if sample==1 & missing(spr_max)
    if r(N) > 0 {
        di as result _n "      Onsets with sample==1 but missing spr_max (" r(N) "):"
        list country year onset_nd onset_def spr_max if sample==1 & missing(spr_max), noobs clean
    }
    quietly count if onset_nd==1 & sample==1 & !missing(spr_max)
    di as result _n "      onset_nd, surviving both restrictions: " r(N) " (Check 1/2 below use this count)"
    quietly count if onset_def==1 & sample==1 & !missing(spr_max)
    di as result "      onset_def, surviving both restrictions: " r(N) " (Check 1/2 below use this count)"
restore

* ══════════════════════════════════════════════════════════════════════════
* CHECK 1 — DESCRIPTIVE: spread intensity at onset, by resolution type
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== CHECK 1: spr_max / spr_mean at onset, non-default vs default-linked ==="

preserve
    quietly keep if sample==1 & (onset_nd==1 | onset_def==1)
    quietly count if missing(spr_max)
    if r(N) > 0 di as error "  ** " r(N) " onset rows have missing spr_max -- reported stats exclude them."

    foreach v in spr_max spr_mean {
        di as result _n "      `v':"
        quietly summarize `v' if onset_nd==1, detail
        di as result "        Non-default   (N=" r(N) "): mean=" %8.1f r(mean) "  p50=" %8.1f r(p50) "  sd=" %8.1f r(sd)
        quietly summarize `v' if onset_def==1, detail
        di as result "        Default-linked(N=" r(N) "): mean=" %8.1f r(mean) "  p50=" %8.1f r(p50) "  sd=" %8.1f r(sd)

        quietly ttest `v', by(onset_def)
        di as result "        t-test (default-linked - non-default): diff=" %6.1f (r(mu_2)-r(mu_1)) "  t=" %5.2f r(t) "  p=" %5.3f r(p)

        quietly ranksum `v', by(onset_def)
        local rsz = r(z)
        local rsp = 2*(1 - normal(abs(`rsz')))
        di as result "        Wilcoxon rank-sum: z=" %5.2f `rsz' "  p=" %5.3f `rsp'
    }

    capture noisily graph box spr_max, over(onset_def, relabel(1 "Non-default" 2 "Default-linked")) ///
        ytitle("Peak EMBIG spread at onset (bps)") title("Spread intensity at onset, by resolution type") ///
        graphregion(color(white)) bgcolor(white)
    if _rc == 0 {
        graph export "$figs/fig_spreadintensity_by_crisistype.png", replace width(1600)
        di as result "      Figure saved: fig_spreadintensity_by_crisistype.png"
    }
restore

* ══════════════════════════════════════════════════════════════════════════
* CHECK 2 — DIRECT TEST, exported Table-1-style: among onsets only (i.e.
* comparing default-linked onsets directly against non-default onsets, NOT
* onset-vs-tranquil as in Check 3/08c), does spr_max predict onset_def?
* THIS is the test that answers the supervisor's actual question -- Check 3
* below compares onset vs tranquil, which spr_max separates near-mechanically
* since the very definition of a spread crisis requires spr_max to cross a
* threshold; Check 2 has no such circularity, since both onset_nd==1 and
* onset_def==1 rows already satisfy that threshold by construction.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== CHECK 2: probit of onset_def on spr_max, among onsets only (exported) ==="

preserve
    quietly keep if sample==1 & (onset_nd==1 | onset_def==1)

    eststo clear
    quietly probit onset_def spr_max, vce(cluster cid)
    eststo dvsnd_bare
    estadd local ctrls "No"

    quietly probit onset_def spr_max $ctrl_core, vce(cluster cid)
    eststo dvsnd_ctrl
    estadd local ctrls "Yes"

    capture esttab dvsnd_bare dvsnd_ctrl using "$tabs/table_default_vs_nondefault_spreadintensity.rtf", replace ///
        b(4) se(4) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
        mtitles("Default vs.\ non-default" "Default vs.\ non-default") ///
        order(spr_max l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2) ///
        coeflabel(spr_max "Peak EMBIG spread at onset (bps)" ///
                  l1_gdpg "GDP growth" ///
                  l_debt "Public debt-to-GDP ratio" ///
                  l_banking_crisis "Banking crisis dummy" ///
                  l_govexp "Government expenditure-to-GDP ratio" ///
                  l_open "Trade openness" ///
                  l_credit_bank "Bank credit-to-GDP ratio" ///
                  l_lninfl "Inflation" ///
                  exchange2 "Nominal exchange-rate") ///
        stats(ctrls N, labels("Baseline controls" "Observations")) ///
        title("Robustness: does spread intensity at onset predict default vs.\ non-default resolution?") ///
        addnotes("Sample restricted to onset years only (onset_nd==1 | onset_def==1), i.e.\ this compares" ///
                 "default-linked onsets directly against non-default onsets -- NOT onset vs.\ tranquil years" ///
                 "as in the paper's Table 1 or in Check 3/table_first_stage_spreadintensity.rtf below." ///
                 "Dependent variable: onset_def (1 = default-linked, 0 = non-default). Pooled probit, robust" ///
                 "standard errors clustered by country in parentheses." ///
                 "* p<0.10, ** p<0.05, *** p<0.01.")
    if _rc == 608 di as error "  ** table_default_vs_nondefault_spreadintensity.rtf is OPEN IN WORD -- close it and re-run."
    else if _rc  di as error "  ** Table (default vs non-default): esttab failed (rc=" _rc ")"
    else di as result "Default-vs-non-default table saved: $tabs/table_default_vs_nondefault_spreadintensity.rtf"
restore

* ══════════════════════════════════════════════════════════════════════════
* CHECK 3 — AUGMENTED TABLE 1: spr_max added to 08c's onset-vs-tranquil design
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== CHECK 3: augmented Table 1 (spr_max added to Z2) ==="

local X    $ctrl_core
local Z2   l_fedfunds l_contagion_dist_atdef years_since_def_onset
local Z3   spr_max

eststo clear

capture program drop _fscol_si
program define _fscol_si
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

_fscol_si fs_nd_si  "onset_nd"  "sample==1 & onset_def==0" "`X'" "`Z2' `Z3'"
_fscol_si fs_def_si "onset_def" "sample==1 & onset_nd==0"  "`X'" "`Z2' `Z3'"

di as result _n "      col               chi2(pred+spr_max)   p       AUROC(ctrl only)  AUROC(+pred+spr_max)"
foreach c in fs_nd_si fs_def_si {
    quietly estimates restore `c'
    di as result %-16s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " %8.3f e(aurocctrl) "        " %6.3f e(auroc)
}
di as result "      Compare AUROC(+pred+spr_max) above to table_first_stage.rtf's own AUROC(with predictors)"
di as result "      (08c_first_stage_table.do) to see whether adding spr_max improves classification, and"
di as result "      check spr_max's own coefficient/z below for whether it is individually significant."

capture esttab fs_nd_si fs_def_si using "$tabs/table_first_stage_spreadintensity.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("Non-default" "Default-linked") ///
    order(spr_max l_fedfunds l_contagion_dist_atdef years_since_def_onset ///
          l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2) ///
    coeflabel(spr_max "Peak EMBIG spread at onset (bps)" ///
              l_fedfunds "US federal funds rate" ///
              l_contagion_dist_atdef "Contagion, based on default-linked crisis (AT-database-wide donors)" ///
              years_since_def_onset "Years since last default-linked onset" ///
              l1_gdpg "GDP growth" ///
              l_debt "Public debt-to-GDP ratio" ///
              l_banking_crisis "Banking crisis dummy" ///
              l_govexp "Government expenditure-to-GDP ratio" ///
              l_open "Trade openness" ///
              l_credit_bank "Bank credit-to-GDP ratio" ///
              l_lninfl "Inflation" ///
              exchange2 "Nominal exchange-rate") ///
    refcat(spr_max "Predictors (incl. spread intensity)" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp aurocctrl auroc N, ///
          labels("Chi-squared for predictors" " p-value of Chi-squared" "AUROC, controls only" "AUROC, with predictors" "Observations") ///
          fmt(2 3 3 3 0)) ///
    title("Robustness: first-stage probit with spread intensity (spr_max) added") ///
    addnotes("Same design as the paper's Table 1 (08c_first_stage_table.do), with the peak EMBIG spread level at onset (spr_max, bps)" ///
             "added to the predictor set. NOT the paper's adopted Table 1 -- a standalone robustness check." ///
             "Dependent variable: dummy = 1 in the onset year of the indicated crisis type; each type predicted vs tranquil years." ///
             "Pooled probit, no country fixed effects. Robust standard errors clustered by country in parentheses." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage_spreadintensity.rtf is OPEN IN WORD -- close it and re-run to refresh."
else if _rc  di as error "  ** Table (spread intensity): esttab failed (rc=" _rc ")"
else di as result "Augmented first-stage table saved: $tabs/table_first_stage_spreadintensity.rtf"

di as result _n "08e_first_stage_table_spreadintensity.do complete."
