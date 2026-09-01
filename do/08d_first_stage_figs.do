/*===========================================================================
  08D_FIRST_STAGE_FIGS.DO
  Figures for the onset first-stage probit — the visual counterpart of
  08c_first_stage_table.do's diagnostics, mirroring
  21c_first_stage_figs_flow.do (the flow-tier analogue) for this file's own
  onset specification, styled on Asonuma et al. (2024)'s own first-stage
  figures (their companion code to Table 1):

    Fig. A  Kernel density of the predicted probability, treated vs control,
            one panel per crisis type — their "propensity overlap" figure.
    Fig. B  Nested ROC comparison, one panel per crisis type — their AUROC
            curves for successively richer specifications.

  WHY THIS IS A SEPARATE FILE, NOT ADDED TO 08c
  ----------------------------------------------
  08c's job is the table export (esttab); it already computes both AUROC
  numbers used here (aurocctrl/auroc) but only as scalars in a table. These
  figures need their own predicted-probability variables (one per nested
  model, per crisis type) and their own graph combine calls, which is a
  different kind of work from esttab -- same reasoning 21c gives for being
  separate from 21b. Per 21c's own header, neither 08c nor 21b had ever
  built these figures before this file; this closes that gap for the onset
  tier, matching what the flow tier already has.

  WHY THE ROC COMPARISON HAS TWO CURVES, NOT THREE
  ---------------------------------------------------------------------
  Their nested ROC compares THREE specifications: Country FEs alone,
  +Controls, +Controls+Predictors — because their probit's $convar/$cf
  carries country dummies (c1-c74, noconstant) directly as regressors.
  This project's onset probit is deliberately POOLED, no country FE
  (08c's own header: "rare-event propensity: country FE separate/overfit
  with ~20 events"). There is therefore no "country FEs alone" curve
  available here, same as the flow tier. This file draws the two curves
  this project actually has and that 08c already reports as numbers:
  CONTROLS ONLY ($ctrl_core) vs CONTROLS + PREDICTORS (+ Z2) — the same
  aurocctrl/auroc pair from 08c's table, as a picture rather than two cells.

  SAMPLE AND MODELS — IDENTICAL TO 08c
  -------------------------------------
  sample==1, rival type dropped, pooled probit, clustered SEs by country.
  `X' = $ctrl_core. `Z2' = l_fedfunds, l_contagion_dist, years_since_onset --
  both the contagion and recency terms are generic (any onset type), matching
  08c's adopted predictor set. The kernel
  density figure uses the FULL model's (controls+predictors) predicted
  probability, matching the reference paper's own kdensity figure, which is
  drawn from their single richest probit ($cf $convar $instrument), not a
  nested comparison.

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
  on that run. Matches 08c's console output, which reports the same test.

  Outputs
  -------
    "$figs/fig_kdensity.pdf/.png"  predicted-probability density, nd/def
    "$figs/fig_roc.pdf/.png"       nested ROC (controls vs +predictors), nd/def

  Reads only $clean/panel_lp.dta — self-contained, run after 08c for the
  reader's sake (same numbers), not as a data dependency.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

foreach v in onset_nd onset_def sample years_since_onset l_contagion_dist {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel_lp.dta — re-run 18_transforms.do first."
        exit 111
    }
}

local X  $ctrl_core
local Z2 l_fedfunds l_contagion_dist years_since_onset

local c_nd   "0 84 166"
local c_def  "157 36 73"

* ══════════════════════════════════════════════════════════════════════════
* FIT BOTH NESTED MODELS PER TYPE, SAVE PREDICTED PROBABILITIES
* ══════════════════════════════════════════════════════════════════════════
foreach s in nd def {
    if "`s'" == "nd" local ifcond "sample==1 & onset_def==0"
    else             local ifcond "sample==1 & onset_nd==0"

    capture drop _p1_`s' _p2_`s'
    quietly probit onset_`s' `X' if `ifcond', vce(cluster cid)
    quietly lroc, nograph
    local auc1_`s' : display %4.2f r(area)
    quietly predict double _p1_`s' if `ifcond', pr

    quietly probit onset_`s' `X' `Z2' if `ifcond', vce(cluster cid)
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
        (kdensity _p2_`s' if onset_`s'==1 & inrange(_p2_`s', 0.01, 0.6), ///
            lwidth(thick) lcolor("`clr'")) ///
        (kdensity _p2_`s' if onset_`s'==0 & inrange(_p2_`s', 0.01, 0.6), ///
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
    subtitle("Onset first-stage probit, controls + predictors. Density trimmed to [0.01, 0.6].", size(small)) ///
    xsize(7) ysize(3.2)
capture graph export "$figs/fig_kdensity.pdf", replace
if _rc di as error "  ** fig_kdensity.pdf export failed (rc=" _rc ") — is it open?"
else {
    capture graph export "$figs/fig_kdensity.png", replace width(1400)
    di as result "Figure saved: fig_kdensity.pdf/.png"
}
foreach nm of local gnames_a {
    capture graph drop `nm'
}

* ══════════════════════════════════════════════════════════════════════════
* FIGURE B — NESTED ROC: CONTROLS ONLY vs CONTROLS + PREDICTORS
*
* Two curves, not three — see header for why a "country FEs alone" curve
* does not exist in this project's onset probit.
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
    quietly roccomp onset_`s' _p1_`s' _p2_`s' if !missing(_p1_`s',_p2_`s')
    local rocnote = "H0: equal areas -- chi2(1)=" + string(r(chi2), "%4.2f") + ", p=" + string(r(p), "%5.3f")
    local rocp_`s' = r(p)
    if r(p) > `worstp_b' local worstp_b = r(p)

    roccomp onset_`s' _p1_`s' _p2_`s' if !missing(_p1_`s',_p2_`s'), ///
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
capture graph export "$figs/fig_roc.pdf", replace
if _rc di as error "  ** fig_roc.pdf export failed (rc=" _rc ") — is it open?"
else {
    capture graph export "$figs/fig_roc.png", replace width(1400)
    di as result "Figure saved: fig_roc.pdf/.png"
}
foreach nm of local gnames_b {
    capture graph drop `nm'
}

capture drop _p1_nd _p2_nd _p1_def _p2_def

di as result _n "08d_first_stage_figs.do complete."
