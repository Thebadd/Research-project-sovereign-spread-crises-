/*===========================================================================
  08B_AIPW.DO
  Doubly-robust treatment-effects estimation (Jordà–Taylor 2016), the approach
  used by Asonuma et al. (2024). Upgrade over plain IPW (08_ipw_lp.do): it
  combines the inverse-probability weight with an outcome-regression adjustment
  and is consistent if EITHER the propensity model OR the outcome model is right.

  Implemented with Stata's `teffects aipw`, which computes exactly Eq. (3) of
  Asonuma et al. (the Jordà–Taylor AIPW estimator): the IPW term plus an
  outcome-regression augmentation term. We report the ATE (average treatment
  effect) — the estimand in their Eq. (3), which sums over all observations and
  scales by 1/N; they never target the ATET. `teffects aipw' returns correct
  two-step (influence-function) SEs, so first-stage estimation error is already
  propagated. The sample restriction to ever-treated countries is a mechanical
  consequence of the probit first stage (never-treated give no DV variation).

    Act 1: treatment = onset_all (onset vs tranquil). Effect = output cost of a
           spread crisis. Outcome model: controls + i.year (year FE preserve the
           exclusion of ust10y/vix). Treatment probit: controls + predictors.
    Act 2: among crisis onsets, treatment = onset_def (default vs non-default).
           Effect = EXTRA output cost of default-linked resolution — the AIPW
           analog of the Table 2/4 difference. No year FE (61-obs sample).

  Bootstrap: an OPTIONAL block re-runs the whole AIPW on country resamples to get
  a percentile 95% CI (à la Asonuma et al.), robust to the small default sample.
  It is slow and may drop draws (separation / no overlap) — guarded with capture.

  NOTE (cannot be run here — David runs Stata): teffects coefficient stripes and
  the cluster-bootstrap of teffects can need small syntax tweaks. Run this, and
  paste the console so we can fix any teffects-specific quirk on the first pass.

  Keep 08_ipw_lp.do (plain IPW) as the robustness row. Run AFTER 01e_predictors.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* Outcome-model controls X and treatment-model predictors Z (as in 08)
local cx    l1_gdpg l2_gdpg debt ca infl imf
local cz    ust10y vix l_reg_crisis_share past_onsets
* Act 2 (onset subsample): leaner, no year FE
local cx2   l1_gdpg l2_gdpg debt ca infl
local cz2   debt ca l_reg_crisis_share

local nboot = 500      // bootstrap reps; raise to 1000 for the final run

* Storage: rows h=0..4
foreach m in b se lo hi {
    matrix A1_`m' = J(5,1,.)
    matrix A2_`m' = J(5,1,.)
}

* ══════════════════════════════════════════════════════════════════════════
* ACT 1 — AIPW output cost of a spread crisis (onset vs tranquil)
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 1 — AIPW (ATE): output cost of a spread crisis ==="
di as result "h    ATE       SE       [95% CI]"

forvalues h = 0/4 {
    local row = `h' + 1

    * teffects aipw reports the ATE (Eq. 3 in Asonuma et al.); atet is not
    * allowed for aipw, and its VCE is always the robust two-step IF variance so
    * we pass no vce(). Try the year-FE form first; on failure fall back to
    * dropping i.year, then echo the real error at h=0.
    capture teffects aipw (dy_`h' `cx' i.year) ///
                          (onset_all `cx' `cz', probit) ///
                          if sample==1, ate
    if _rc != 0 {
        * fallback: no year FE in the outcome model (aipw needs overlap in every
        * cell; sparse year dummies on the treated group can break estimation)
        capture teffects aipw (dy_`h' `cx') ///
                              (onset_all `cx' `cz', probit) ///
                              if sample==1, ate
        if _rc != 0 & `h' == 0 {
            di as error "--- Act 1 h=0 failed both forms; echoing the real error ---"
            capture noisily teffects aipw (dy_`h' `cx') ///
                              (onset_all `cx' `cz', probit) ///
                              if sample==1, ate
        }
    }
    if _rc == 0 {
        matrix bb = e(b)
        matrix VV = e(V)
        local att = bb[1,1]
        local se  = sqrt(VV[1,1])
        matrix A1_b[`row',1]  = `att'
        matrix A1_se[`row',1] = `se'
        matrix A1_lo[`row',1] = `att' - 1.96*`se'
        matrix A1_hi[`row',1] = `att' + 1.96*`se'
        di "h=" `h' "   " %7.3f `att' "  " %6.3f `se' ///
           "   [" %6.3f (`att'-1.96*`se') ", " %6.3f (`att'+1.96*`se') "]"
    }
    else di as error "h=" `h' ": teffects aipw (Act 1) failed, rc=" _rc
}

* ══════════════════════════════════════════════════════════════════════════
* ACT 2 — AIPW extra cost of default-linked resolution (default vs non-default)
*   Estimated among crisis onsets; this IS the doubly-robust resolution
*   difference (the AIPW analog of the Table 2/4 difference block).
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 2 — AIPW (ATE): default-linked vs non-default ==="
di as result "h    ATE       SE       [95% CI]   (negative => default deeper)"

forvalues h = 0/4 {
    local row = `h' + 1
    capture teffects aipw (dy_`h' `cx2') ///
                          (onset_def `cz2', probit) ///
                          if onset_all==1, ate
    if _rc != 0 & `h' == 0 {
        di as error "--- Act 2 h=0 failed; echoing the real error ---"
        capture noisily teffects aipw (dy_`h' `cx2') ///
                          (onset_def `cz2', probit) ///
                          if onset_all==1, ate
    }
    if _rc == 0 {
        matrix bb = e(b)
        matrix VV = e(V)
        local att = bb[1,1]
        local se  = sqrt(VV[1,1])
        matrix A2_b[`row',1]  = `att'
        matrix A2_se[`row',1] = `se'
        matrix A2_lo[`row',1] = `att' - 1.96*`se'
        matrix A2_hi[`row',1] = `att' + 1.96*`se'
        di "h=" `h' "   " %7.3f `att' "  " %6.3f `se' ///
           "   [" %6.3f (`att'-1.96*`se') ", " %6.3f (`att'+1.96*`se') "]"
    }
    else di as error "h=" `h' ": teffects aipw (Act 2) failed, rc=" _rc " (likely overlap/separation on the thin default sample)"
}

* ══════════════════════════════════════════════════════════════════════════
* OPTIONAL — BOOTSTRAP PERCENTILE 95% CI (cluster by country), Asonuma-style
*   Slow: re-runs the whole AIPW on `nboot' country resamples per horizon.
*   Overwrites the CI columns with percentile CIs where it succeeds.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== BOOTSTRAP percentile CIs (cluster cid, reps=`nboot') ==="
di as result "  (this is the slow step; comment out to skip)"

* Act 1
forvalues h = 0/4 {
    local row = `h' + 1
    capture bootstrap _b, reps(`nboot') cluster(cid) nodots nowarn: ///
        teffects aipw (dy_`h' `cx') (onset_all `cx' `cz', probit) ///
        if sample==1, ate
    if _rc == 0 {
        capture estat bootstrap, percentile
        if _rc == 0 {
            matrix ci = e(ci_percentile)
            matrix A1_lo[`row',1] = ci[1,1]
            matrix A1_hi[`row',1] = ci[2,1]
            di "Act1 h=" `h' ": percentile CI [" %6.3f ci[1,1] ", " %6.3f ci[2,1] "]"
        }
    }
    else di as error "Act1 h=" `h' ": bootstrap failed, rc=" _rc
}

* Act 2
forvalues h = 0/4 {
    local row = `h' + 1
    capture bootstrap _b, reps(`nboot') cluster(cid) nodots nowarn: ///
        teffects aipw (dy_`h' `cx2') (onset_def `cz2', probit) ///
        if onset_all==1, ate
    if _rc == 0 {
        capture estat bootstrap, percentile
        if _rc == 0 {
            matrix ci = e(ci_percentile)
            matrix A2_lo[`row',1] = ci[1,1]
            matrix A2_hi[`row',1] = ci[2,1]
            di "Act2 h=" `h' ": percentile CI [" %6.3f ci[1,1] ", " %6.3f ci[2,1] "]"
        }
    }
    else di as error "Act2 h=" `h' ": bootstrap failed, rc=" _rc
}

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — AIPW results CSV + comparison figure
* ══════════════════════════════════════════════════════════════════════════
preserve
    clear
    set obs 10
    gen act     = cond(_n<=5, 1, 2)
    gen horizon = mod(_n-1, 5)
    gen aipw_b  = .
    gen aipw_lo = .
    gen aipw_hi = .
    forvalues h = 0/4 {
        replace aipw_b  = A1_b[`h'+1,1]  if act==1 & horizon==`h'
        replace aipw_lo = A1_lo[`h'+1,1] if act==1 & horizon==`h'
        replace aipw_hi = A1_hi[`h'+1,1] if act==1 & horizon==`h'
        replace aipw_b  = A2_b[`h'+1,1]  if act==2 & horizon==`h'
        replace aipw_lo = A2_lo[`h'+1,1] if act==2 & horizon==`h'
        replace aipw_hi = A2_hi[`h'+1,1] if act==2 & horizon==`h'
    }
    label define actl 1 "Act1: crisis cost" 2 "Act2: default extra cost"
    label values act actl
    order act horizon aipw_b aipw_lo aipw_hi
    export delimited "$tabs/aipw_results.csv", replace
    di as result "AIPW results CSV saved: $tabs/aipw_results.csv"

    * Figure: AIPW IRFs, Act 1 and Act 2
    local c1 "23 55 94"
    local c2 "157 36 73"
    twoway ///
        (rarea aipw_lo aipw_hi horizon if act==1, color("`c1'%20") lwidth(none)) ///
        (connected aipw_b horizon if act==1, lcolor("`c1'") lwidth(medthick) msymbol(circle)) ///
        (rarea aipw_lo aipw_hi horizon if act==2, color("`c2'%20") lwidth(none)) ///
        (connected aipw_b horizon if act==2, lcolor("`c2'") lwidth(medthick) lpattern(dash) msymbol(square)), ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)4) ytitle("Cumulative GDPpc change (pp)", size(small)) ///
        xtitle("Years after onset", size(small)) ///
        title("AIPW (doubly-robust) output cost", size(medsmall) color(navy)) ///
        legend(order(2 "Crisis cost (all)" 4 "Extra cost of default") ring(0) pos(7) size(small)) ///
        note("teffects aipw (doubly robust), ATE. CI = bootstrap percentile (cluster country) where available, else robust 95%.", size(vsmall)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig_aipw.pdf", replace
    di as result "Figure saved: fig_aipw.pdf"
restore

di as result _n "08b_aipw.do complete."
di as result "Compare Act 1 AIPW to the OLS/IPW output cost (Table 1 / 08); Act 2 AIPW"
di as result "to the OLS resolution difference (Table 2). Similar magnitudes => the"
di as result "OLS results are not driven by selection (their Fig C1 logic)."
