/*===========================================================================
  21C_FIRST_STAGE_FIGS_FLOW.DO
  Figures for the flow first-stage probit — the visual counterpart of
  21b_first_stage_table_flow.do's diagnostics, styled on Asonuma et al.
  (2024)'s own first-stage figures (their companion code to Table 1):

    Fig. A  Kernel density of the predicted probability, treated vs control,
            one panel per crisis type — their "propensity overlap" figure.
    Fig. B  Nested ROC comparison, one panel per crisis type — their AUROC
            curves for successively richer specifications.

  WHY THIS IS A SEPARATE FILE, NOT ADDED TO 21b
  ----------------------------------------------
  21b's job is the table export (esttab); it already computes both AUROC
  numbers used here (aurocctrl/auroc) but only as scalars in a table. These
  figures need their own predicted-probability variables (one per nested
  model, per crisis type) and their own graph combine calls, which is a
  different kind of work from esttab — same reasoning this project already
  applies to keeping estimation and figure-drawing visually distinct within
  a file's own sections (e.g. 21_aipw_flow.do's Section 5). Neither 08c
  (the onset-tier analogue) nor 21b has ever built these figures, so there
  is no existing companion file being extended here.

  WHY THE ROC COMPARISON HAS TWO CURVES, NOT THREE
  ---------------------------------------------------------------------
  Their nested ROC compares THREE specifications: Country FEs alone,
  +Controls, +Controls+Predictors — because their probit's $convar/$cf
  carries country dummies (c1-c74, noconstant) directly as regressors.
  This project's flow probit is deliberately POOLED, no country FE: adding
  i.cid was tested in 21_aipw_flow.do and rejected — the def arm collapsed
  to complete separation on this project's much smaller panel (see that
  file's header, "COUNTRY FE IN THE PROBIT — TESTED AND REJECTED"). There is
  therefore no "country FEs alone" curve available here. This file draws the
  two curves this project actually has and that 21b already reports as
  numbers: CONTROLS ONLY ($ctrl_core_flowbase) vs CONTROLS + PREDICTORS
  (+ cz_recency) — the same aurocctrl/auroc pair from 21b's table, as a
  picture rather than two cells.

  SAMPLE AND MODELS — IDENTICAL TO 21b
  -------------------------------------
  sample_flow==1 & continuation==0, rival type dropped, pooled probit,
  clustered SEs by country. `X' = $ctrl_core_flowbase (or the alternate, if
  toggled — same toggle as 21b, kept in sync manually since this file does
  not source 21b). `Z' = cz_recency (l_fedfunds, l_contagion_dist,
  years_since_def_onset). The kernel density figure uses the FULL model's
  (controls+predictors) predicted probability, matching the reference
  paper's own kdensity figure, which is drawn from their single richest
  probit ($cf $convar $instrument), not a nested comparison.

  DENSITY TRIM: predicted probabilities outside [0.01, 0.6] are excluded
  from the kdensity plot only (not from estimation), matching the reference
  figure's own axis and avoiding a density dominated by the near-zero mass
  every pooled probit on a rare-event outcome produces.

  EACH ROC PANEL CARRIES roccomp's FORMAL TEST, NOT JUST THE TWO AUROCs
  -----------------------------------------------------------------------
  The AUROC delta shown in the legend is a point difference, not a
  significance test. roccomp's chi2 test for equal correlated ROC areas is
  printed as a note under each panel (H0: equal areas), and the combined
  figure's subtitle is set from the actual computed p-values at run time
  (not a fixed claim) -- whether it reads "established" or "modest, not
  absent" depends on whether both arms clear the conventional 5 pct level
  on that run. Under the tot_chg-based core this never happened (nd p=.084,
  def p=.055 as of the run that motivated adding this); under the adopted
  exchange2-based core, on a larger fitting sample (exchange2's fuller
  onset coverage recovers rows tot_chg's missingness was dropping), both
  arms have since cleared 5 pct in at least one run (nd p=.028, def
  p=.033). Matches 21b's console output, which reports the same test.

  Outputs
  -------
    "$figs/fig_kdensity_flow.pdf/.png"  predicted-probability density, nd/def
    "$figs/fig_roc_flow.pdf/.png"       nested ROC (controls vs +predictors), nd/def

  Reads only $clean/panel_lp.dta — self-contained, run after 21b for the
  reader's sake (same numbers), not as a data dependency.
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

* Same toggle and control sets as 21b_first_stage_table_flow.do — kept in
* sync manually, since this file does not source that one. 0=adopted core,
* 1=adopted core + tot_chg (the term exchange2 replaced)/l_imf.
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
local X = cond(`flow_ctrl_variant'==1, "$ctrl_core_flowplus", "$ctrl_core_flowbase")
local Z l_fedfunds l_contagion_dist years_since_def_onset

local c_nd   "0 84 166"
local c_def  "157 36 73"

* ══════════════════════════════════════════════════════════════════════════
* FIT BOTH NESTED MODELS PER TYPE, SAVE PREDICTED PROBABILITIES
* ══════════════════════════════════════════════════════════════════════════
foreach s in nd def {
    if "`s'" == "nd" local ifcond "sample_flow==1 & continuation==0 & in_crisis_def==0"
    else             local ifcond "sample_flow==1 & continuation==0 & in_crisis_nd==0"

    capture drop _p1_`s' _p2_`s'
    quietly probit in_crisis_`s' `X' if `ifcond', vce(cluster cid)
    quietly lroc, nograph
    local auc1_`s' : display %4.2f r(area)
    quietly predict double _p1_`s' if `ifcond', pr

    quietly probit in_crisis_`s' `X' `Z' if `ifcond', vce(cluster cid)
    quietly lroc, nograph
    local auc2_`s' : display %4.2f r(area)
    quietly predict double _p2_`s' if `ifcond', pr
}

* ══════════════════════════════════════════════════════════════════════════
* FIGURE A — PREDICTED-PROBABILITY DENSITY, TREATED VS CONTROL
*
* From the richer (controls+predictors) model, matching the reference
* figure's use of its single richest probit. Trimmed to [0.01,0.6] for the
* plot only, matching their axis.
* ══════════════════════════════════════════════════════════════════════════
local gnames_a
foreach s in nd def {
    local ttl = cond("`s'"=="nd", "Non-default", "Default-linked")
    local clr = cond("`s'"=="nd", "`c_nd'", "`c_def'")

    twoway ///
        (kdensity _p2_`s' if in_crisis_`s'==1 & inrange(_p2_`s', 0.01, 0.6), ///
            lwidth(thick) lcolor("`clr'")) ///
        (kdensity _p2_`s' if in_crisis_`s'==0 & inrange(_p2_`s', 0.01, 0.6), ///
            lwidth(thick) lcolor("`clr'") lpattern(dash)), ///
        graphregion(color(white)) plotregion(color(white)) ///
        legend(order(1 "Treatment group" 2 "Control group") ring(0) pos(1) size(small)) ///
        ytitle("Probability density", size(small)) xtitle("Predicted probability", size(small)) ///
        title("`ttl'", size(medium) color("`clr'")) ///
        name(gk_`s', replace) nodraw
    local gnames_a `gnames_a' gk_`s'
}
graph combine `gnames_a', rows(1) graphregion(color(white)) ///
    title("Predicted Probability of Onset: Treated vs Control", size(medsmall) color(navy)) ///
    subtitle("Flow first-stage probit, controls + predictors. Density trimmed to [0.01, 0.6].", size(small)) ///
    xsize(7) ysize(3.2)
capture graph export "$figs/fig_kdensity_flow.pdf", replace
if _rc di as error "  ** fig_kdensity_flow.pdf export failed (rc=" _rc ") — is it open?"
else {
    capture graph export "$figs/fig_kdensity_flow.png", replace width(1400)
    di as result "Figure saved: fig_kdensity_flow.pdf/.png"
}
foreach nm of local gnames_a {
    capture graph drop `nm'
}

* ══════════════════════════════════════════════════════════════════════════
* FIGURE B — NESTED ROC: CONTROLS ONLY vs CONTROLS + PREDICTORS
*
* Two curves, not three — see header for why a "country FEs alone" curve
* does not exist in this project's flow probit.
* ══════════════════════════════════════════════════════════════════════════
local gnames_b
local worstp_b = 0
foreach s in nd def {
    local ttl = cond("`s'"=="nd", "Non-default", "Default-linked")

    * roccomp's own chi2 test for equal correlated ROC areas -- the formal
    * significance test the raw AUROC delta cannot substitute for. Run once,
    * quietly, to capture r(chi2)/r(p) BEFORE the graphing call (whose own
    * r() is not yet populated at the point its own option string is built),
    * then reused as a literal string in the panel note below.
    quietly roccomp in_crisis_`s' _p1_`s' _p2_`s' if !missing(_p1_`s',_p2_`s')
    local rocnote = "H0: equal areas -- chi2(1)=" + string(r(chi2), "%4.2f") + ", p=" + string(r(p), "%5.3f")
    local rocp_`s' = r(p)
    if r(p) > `worstp_b' local worstp_b = r(p)

    roccomp in_crisis_`s' _p1_`s' _p2_`s' if !missing(_p1_`s',_p2_`s'), ///
        graph summary name(gr_`s', replace) graphregion(color(white)) nodraw ///
        plot1opts(lcolor(red) mcolor(red) msymbol(circle)) ///
        plot2opts(lcolor(green) mcolor(green) msymbol(diamond)) ///
        title("`ttl'", size(medium)) ///
        note("`rocnote'", size(vsmall)) ///
        legend(position(5) region(lwidth(none)) size(vsmall) cols(1) ring(0) ///
            order(1 "Controls: `auc1_`s''" 2 "Controls+Predictors: `auc2_`s''"))
    local gnames_b `gnames_b' gr_`s'
}
local subtxt = cond(`worstp_b' < 0.05, ///
    "ROC area under the curve; 0.50 = no classification power, 1.00 = perfect. No country-FE curve —" + ///
        " see header. Both panels' gaps clear the conventional 5 pct level (roccomp chi2 test, see notes below" + ///
        " each panel) -- an established, not merely suggestive, classification gain on this project's sample.", ///
    "ROC area under the curve; 0.50 = no classification power, 1.00 = perfect. No country-FE curve —" + ///
        " see header. At least one panel's gap does not clear the conventional 5 pct level (roccomp chi2 test," + ///
        " see notes below each panel) -- modest, not absent, on this project's sample size.")
graph combine `gnames_b', graphregion(color(white)) ///
    title("First-Stage Classification: Controls vs Controls + Predictors", size(medsmall) color(navy)) ///
    subtitle("`subtxt'", size(small)) ///
    xsize(6) ysize(3.2)
capture graph export "$figs/fig_roc_flow.pdf", replace
if _rc di as error "  ** fig_roc_flow.pdf export failed (rc=" _rc ") — is it open?"
else {
    capture graph export "$figs/fig_roc_flow.png", replace width(1400)
    di as result "Figure saved: fig_roc_flow.pdf/.png"
}
foreach nm of local gnames_b {
    capture graph drop `nm'
}

capture drop _p1_nd _p2_nd _p1_def _p2_def

di as result _n "21c_first_stage_figs_flow.do complete."
