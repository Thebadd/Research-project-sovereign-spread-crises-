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
  countries with no variation in D. Act 1 = onset vs tranquil. Act 2 = the
  resolution split shown as TWO level IRFs (non-default vs tranquil, default vs
  tranquil) — the doubly-robust analog of the Table 2 nd/def coefficients; their
  vertical gap is the extra cost of default. Mirrors the IPW Act-2 figure (fig8).

  Keep 08_ipw_lp.do (plain IPW) as the robustness row. Run AFTER 01e_predictors.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* Outcome-model controls X and treatment-model predictors Z (as in 08)
local cx    l1_gdpg l2_gdpg debt ca infl imf
* Act 1 predictors Z1: single global push (fed funds) + contagion + proneness.
local cz    vix l_reg_crisis_share past_onsets
* Act 2 predictors Z2: same, but proneness = past DEFAULT-linked onsets.
local cz_def vix l_reg_crisis_share past_def_onsets

local nboot = 500      // bootstrap reps; raise to 1000+ for the final run

* Storage: rows h=0..4
foreach m in b se lo hi {
    matrix A1_`m'    = J(5,1,.)   // Act 1: crisis cost (all onsets)
    matrix A2nd_`m'  = J(5,1,.)   // Act 2: non-default cost vs tranquil
    matrix A2def_`m' = J(5,1,.)   // Act 2: default-linked cost vs tranquil
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
* ACT 2 — AIPW output cost by resolution, as TWO LEVEL lines vs tranquil
*   (the doubly-robust analog of the Table 2 nd/def coefficients; their vertical
*    gap = the extra default cost). Mirrors the IPW Act-2 figure (fig8): one line
*    for non-default onsets, one for default-linked, each estimated against the
*    tranquil control with country FE + the full (X,Z) probit — exactly like
*    Act 1. The OTHER resolution type is dropped from each estimation so the
*    control group is clean tranquil years, not the rival crisis type.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 2 — AIPW (ATE) output cost by resolution (vs tranquil) ==="
di as result "type  h    ATE       SE(boot)  [95% percentile CI]"

* Loop over the two resolution types. `Dvar' = treatment dummy; `drop' = the
* rival onset dummy to exclude; `stub' = storage-matrix prefix; `lab' = label.
foreach spec in "onset_nd onset_def A2nd non-default" ///
                "onset_def onset_nd A2def default-linked" {
    gettoken Dvar spec  : spec
    gettoken drop spec  : spec
    gettoken stub lab   : spec
    local lab = trim("`lab'")

    forvalues h = 0/4 {
        local row = `h' + 1

        * point estimate: this type's onsets vs tranquil (rival type dropped)
        capture _aipw dy_`h' `Dvar' if sample==1 & `drop'==0, ///
            omodel(`cx') pmodel(`cx' `cz_def') fe(cid)
        if _rc {
            di as error "`lab' h=" `h' ": point estimate failed, rc=" _rc
            continue
        }
        local pt = r(theta)
        matrix `stub'_b[`row',1] = `pt'

        * cluster bootstrap (fresh FE ids), recompute Eq. (3)
        tempname pf2
        tempfile bf2
        quietly postfile `pf2' double theta using "`bf2'", replace
        forvalues b = 1/`nboot' {
            preserve
                capture drop _bid
                bsample, cluster(cid) idcluster(_bid)
                capture _aipw dy_`h' `Dvar' if sample==1 & `drop'==0, ///
                    omodel(`cx') pmodel(`cx' `cz_def') fe(_bid)
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
                matrix `stub'_se[`row',1] = r(sd)
                _pctile theta, p(2.5 97.5)
                matrix `stub'_lo[`row',1] = r(r1)
                matrix `stub'_hi[`row',1] = r(r2)
            }
        restore

        di "`lab'  h=" `h' "  " %7.3f `pt' "   " ///
           %7.3f `stub'_se[`row',1] "   [" %7.3f `stub'_lo[`row',1] ", " ///
           %7.3f `stub'_hi[`row',1] "]   (" `ndraw' "/`nboot' draws)"
    }
}

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — AIPW results CSV + figures (Act 1 single line; Act 2 two-line split)
* ══════════════════════════════════════════════════════════════════════════
preserve
    clear
    * Long dataset: 3 series (all / nd / def) x 5 horizons.
    tempname pfx
    tempfile aipwf
    postfile `pfx' str3 series horizon b se lo hi using "`aipwf'", replace
    foreach map in "all A1" "nd A2nd" "def A2def" {
        gettoken sname stub : map
        forvalues h = 0/4 {
            post `pfx' ("`sname'") (`h') ///
                (`stub'_b[`h'+1,1]) (`stub'_se[`h'+1,1]) ///
                (`stub'_lo[`h'+1,1]) (`stub'_hi[`h'+1,1])
        }
    }
    postclose `pfx'
    use "`aipwf'", clear
    label var b  "AIPW ATE (pp)"
    label var lo "95% percentile CI lower"
    label var hi "95% percentile CI upper"
    order series horizon b se lo hi
    export delimited "$tabs/aipw_results.csv", replace
    di as result "AIPW results CSV saved: $tabs/aipw_results.csv"

    * ── Figure A: Act 1 crisis cost (single line + CI band) ──────────────────
    local c1 "23 55 94"
    twoway ///
        (rarea lo hi horizon if series=="all", color("`c1'%20") lwidth(none)) ///
        (connected b horizon if series=="all", lcolor("`c1'") lwidth(medthick) msymbol(circle)), ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)4) ytitle("Cumulative real GDP change (pp)", size(small)) ///
        xtitle("Years after onset", size(small)) ///
        title("AIPW (Asonuma et al. Eq. 3) output cost — all crises", size(medsmall) color(navy)) ///
        legend(off) ///
        note("Hand-coded AIPW (Jordà–Taylor 2016 / Asonuma et al. Eq. 3), ATE. CI = bootstrap percentile (cluster by country).", size(vsmall)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig_aipw.pdf", replace
    di as result "Figure saved: fig_aipw.pdf"

    * ── Figure B: Act 2 resolution split (two lines + CI bands), fig8 palette ─
    local c_nd  "0 84 166"
    local c_def "157 36 73"
    twoway ///
        (rarea lo hi horizon if series=="nd",  color("`c_nd'%18")  lwidth(none)) ///
        (rarea lo hi horizon if series=="def", color("`c_def'%18") lwidth(none)) ///
        (connected b horizon if series=="nd",  lcolor("`c_nd'")  lwidth(medthick) msymbol(circle)) ///
        (connected b horizon if series=="def", lcolor("`c_def'") lwidth(medthick) lpattern(dash) msymbol(square)), ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)4, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Years after crisis onset", size(medsmall)) ///
        ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
        title("AIPW output cost by resolution", size(medium) color(navy)) ///
        subtitle("Doubly-robust (Asonuma et al. Eq. 3), each vs tranquil.", size(small)) ///
        legend(order(3 "Non-default" 4 "Default-linked") ring(0) pos(7) size(small)) ///
        note("Two AIPW level IRFs; shaded = bootstrap 95% percentile CI (cluster by country). Gap = extra cost of default.", size(vsmall)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig_aipw_act2.pdf", replace
    di as result "Figure saved: fig_aipw_act2.pdf"
restore

di as result _n "08b_aipw.do complete."
di as result "fig_aipw.pdf = Act 1 crisis cost; fig_aipw_act2.pdf = the two-line"
di as result "resolution split (non-default vs default, each vs tranquil) — the AIPW"
di as result "twin of the IPW fig8. Compare the nd/def AIPW levels to fig8's IPW lines"
di as result "and to the Table 2 OLS coefficients (similar ordering => selection is not"
di as result "driving the resolution gap; their Fig C1 logic)."
