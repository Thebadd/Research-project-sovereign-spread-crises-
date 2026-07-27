/*===========================================================================
  08B_AIPW.DO
  Augmented Inverse Probability Weighted (AIPW) estimation, coded to match
  Asonuma, Chamon, Erce & Sasahara (2024, JIE) — their Eqs. (1)-(3), which are
  the Jordà & Taylor (2016) AIPW estimator. We implement it BY HAND rather than
  via `teffects aipw` so that (a) it is literally their estimator and (b) the
  bootstrap can be done their way (`bsample` + recompute), avoiding the
  teffects-resample failures (rc=451) on our thin default sample.

  THE THREE STEPS (their §3.2):
    Eq. (1) Outcome regression by OLS with country FE:
              dy_h = a_i + b*D + X*g + u        -> conditional means m1, m0
            (common-slope RA: m1 - m0 = b, the regression-adjustment term).
    Eq. (2) Propensity probit:  Pr(D=1) = Phi(X, Z)   -> phat.
            POOLED (no country FE): country FE cause separation and would drop
            the thin default sample (our locked design choice; the paper uses FE
            but has far more data). Predictors Z enter here ONLY.
    Eq. (3) AIPW estimator (their exact algebraic form):
        Lambda_h = (1/N) sum_i {
            [ D*dy/phat - (1-D)*dy/(1-phat) ]
          - (D - phat)/(phat*(1-phat)) * [ (1-phat)*m1 + phat*m0 ] }
        This is the ATE (sum over ALL i, scaled by 1/N) — the paper's estimand.
        (Algebraically identical to the standard mu1-mu0 + IF-correction form.)

  INFERENCE (their Table 2 note): percentile CIs from `bsample` — resample
  countries (cluster bootstrap), recompute Eq. (3) each draw, take the 2.5/97.5
  percentiles. No analytic SE; the SE column reports the bootstrap SD.

  SAMPLE: restriction to ever-treated countries is mechanical — the probit drops
  countries with no variation in D. Act 1 = onset vs tranquil; Act 2 = default vs
  non-default among onsets (the doubly-robust analog of the Table 2/4 difference).

  Keep 08_ipw_lp.do (plain IPW) as the robustness row. Run AFTER 01e_predictors.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* Outcome-model controls X and treatment-model predictors Z (as in 08)
local cx    l1_gdpg l2_gdpg debt ca infl imf
local cz    ust10y vix l_reg_crisis_share past_onsets
* Act 2 (onset subsample): leaner, no country FE (61 obs)
local cx2   l1_gdpg l2_gdpg debt ca infl
local cz2   debt ca l_reg_crisis_share

local nboot = 500      // bootstrap reps; raise to 1000+ for the final run

* Storage: rows h=0..4
foreach m in b se lo hi {
    matrix A1_`m' = J(5,1,.)
    matrix A2_`m' = J(5,1,.)
}

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — hand-coded AIPW point estimate (Eqs. 1-3), returns r(theta), r(N)
*   _aipw dyvar Dvar [if], omodel(X) pmodel(X Z) [fe(idvar)]
*   `fe(idvar)' adds FE to the OUTCOME regression only (Eq. 1). The probit
*   (Eq. 2) is always pooled. Trims phat to [.01,.99] to bound the weights.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    * omodel = Eq.(1) RHS covariates (X); pmodel = Eq.(2) probit RHS (X and Z).
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'           // common sample for Eq.1 & Eq.2

    tempvar xb m0 m1 ps summ
    * --- Eq. (1): outcome regression -> conditional means m1, m0 ---
    * `fe' names the country-FE variable (cid on real data; a fresh resample id
    * under the cluster bootstrap so duplicated countries are distinct FE groups).
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'   // set D=0
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'   // set D=1

    * --- Eq. (2): pooled propensity probit -> phat (trimmed) ---
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'

    * --- Eq. (3): AIPW summand (Asonuma et al. exact form) ---
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    return scalar theta = r(mean)
    return scalar N     = r(N)
end

* ══════════════════════════════════════════════════════════════════════════
* HELPER — point estimate + bsample percentile CI for one Act/horizon
*   Fills the matrix `stub' (A1_ or A2_) row `row'.
* ══════════════════════════════════════════════════════════════════════════
* (inlined below rather than a second program, to keep matrix names in scope)

* ══════════════════════════════════════════════════════════════════════════
* ACT 1 — AIPW output cost of a spread crisis (onset vs tranquil), ATE
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 1 — AIPW (ATE): output cost of a spread crisis ==="
di as result "h    ATE       SE(boot)  [95% percentile CI]"

forvalues h = 0/4 {
    local row = `h' + 1

    * point estimate on the real data (outcome model with country FE = cid)
    capture _aipw dy_`h' onset_all if sample==1, omodel(`cx') pmodel(`cx' `cz') fe(cid)
    if _rc {
        di as error "h=" `h' ": Act 1 point estimate failed, rc=" _rc
        continue
    }
    local pt = r(theta)
    matrix A1_b[`row',1] = `pt'

    * cluster bootstrap: resample countries (fresh FE ids), recompute Eq. (3)
    tempname pf
    tempfile bf
    quietly postfile `pf' double theta using "`bf'", replace
    forvalues b = 1/`nboot' {
        preserve
            capture drop _bid
            bsample, cluster(cid) idcluster(_bid)
            capture _aipw dy_`h' onset_all if sample==1, omodel(`cx') pmodel(`cx' `cz') fe(_bid)
            if _rc == 0 quietly post `pf' (r(theta))
        restore
    }
    quietly postclose `pf'

    local ndraw = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(theta)
        local ndraw = r(N)
        if `ndraw' >= 50 {
            quietly summarize theta
            matrix A1_se[`row',1] = r(sd)
            _pctile theta, p(2.5 97.5)
            matrix A1_lo[`row',1] = r(r1)
            matrix A1_hi[`row',1] = r(r2)
        }
    restore

    di "h=" `h' "   " %7.3f `pt' "   " ///
       %7.3f A1_se[`row',1] "   [" %7.3f A1_lo[`row',1] ", " %7.3f A1_hi[`row',1] "]" ///
       "   (" `ndraw' "/`nboot' draws)"
}

* ══════════════════════════════════════════════════════════════════════════
* ACT 2 — AIPW extra cost of default-linked resolution (default vs non-default)
*   Among crisis onsets. No country FE in the outcome model (61 obs). ATE.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 2 — AIPW (ATE): default-linked vs non-default ==="
di as result "h    ATE       SE(boot)  [95% CI]   (negative => default deeper)"

forvalues h = 0/4 {
    local row = `h' + 1

    capture _aipw dy_`h' onset_def if onset_all==1, omodel(`cx2') pmodel(`cz2')
    if _rc {
        di as error "h=" `h' ": Act 2 point estimate failed, rc=" _rc ///
            " (likely overlap/separation on the thin default sample)"
        continue
    }
    local pt = r(theta)
    matrix A2_b[`row',1] = `pt'

    tempname pf2
    tempfile bf2
    quietly postfile `pf2' double theta using "`bf2'", replace
    forvalues b = 1/`nboot' {
        preserve
            bsample, cluster(cid)
            capture _aipw dy_`h' onset_def if onset_all==1, omodel(`cx2') pmodel(`cz2')
            if _rc == 0 quietly post `pf2' (r(theta))
        restore
    }
    quietly postclose `pf2'

    local ndraw = 0
    preserve
        quietly use "`bf2'", clear
        quietly count if !missing(theta)
        local ndraw = r(N)
        if `ndraw' >= 50 {
            quietly summarize theta
            matrix A2_se[`row',1] = r(sd)
            _pctile theta, p(2.5 97.5)
            matrix A2_lo[`row',1] = r(r1)
            matrix A2_hi[`row',1] = r(r2)
        }
    restore

    di "h=" `h' "   " %7.3f `pt' "   " ///
       %7.3f A2_se[`row',1] "   [" %7.3f A2_lo[`row',1] ", " %7.3f A2_hi[`row',1] "]" ///
       "   (" `ndraw' "/`nboot' draws)"
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
    gen aipw_se = .
    gen aipw_lo = .
    gen aipw_hi = .
    forvalues h = 0/4 {
        replace aipw_b  = A1_b[`h'+1,1]  if act==1 & horizon==`h'
        replace aipw_se = A1_se[`h'+1,1] if act==1 & horizon==`h'
        replace aipw_lo = A1_lo[`h'+1,1] if act==1 & horizon==`h'
        replace aipw_hi = A1_hi[`h'+1,1] if act==1 & horizon==`h'
        replace aipw_b  = A2_b[`h'+1,1]  if act==2 & horizon==`h'
        replace aipw_se = A2_se[`h'+1,1] if act==2 & horizon==`h'
        replace aipw_lo = A2_lo[`h'+1,1] if act==2 & horizon==`h'
        replace aipw_hi = A2_hi[`h'+1,1] if act==2 & horizon==`h'
    }
    label define actl 1 "Act1: crisis cost" 2 "Act2: default extra cost"
    label values act actl
    order act horizon aipw_b aipw_se aipw_lo aipw_hi
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
        title("AIPW (Asonuma et al. Eq. 3) output cost", size(medsmall) color(navy)) ///
        legend(order(2 "Crisis cost (all)" 4 "Extra cost of default") ring(0) pos(7) size(small)) ///
        note("Hand-coded AIPW (Jordà–Taylor 2016 / Asonuma et al. Eq. 3), ATE. CI = bootstrap percentile (cluster by country).", size(vsmall)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig_aipw.pdf", replace
    di as result "Figure saved: fig_aipw.pdf"
restore

di as result _n "08b_aipw.do complete."
di as result "Compare Act 1 AIPW to the OLS/IPW output cost (Table 1 / 08); Act 2 AIPW"
di as result "to the OLS resolution difference (Table 2). Similar magnitudes => the"
di as result "OLS results are not driven by selection (their Fig C1 logic)."
