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
            POOLED (no country FE): with ~20 events, country dummies in a probit
            separate/overfit and collapse the propensity (our locked design choice;
            the paper uses FE but has 194 events). Country FE enter the OUTCOME
            regression only; no year FE. Predictors Z enter here ONLY.
    Eq. (3) AIPW estimator (their exact algebraic form):
        Lambda_h = (1/N) sum_i {
            [ D*dy/phat - (1-D)*dy/(1-phat) ]
          - (D - phat)/(phat*(1-phat)) * [ (1-phat)*m1 + phat*m0 ] }
        This is the ATE (sum over ALL i, scaled by 1/N) — the paper's estimand.
        (Algebraically identical to the standard mu1-mu0 + IF-correction form.)

  ESTIMATOR: IPWRA, the reference paper's headline. The propensity probit runs
  FIRST; the outcome regression producing the conditional means is then INVERSE-
  PROBABILITY WEIGHTED (their `reg g_h dum g_0 $convar [pweight=invwt]`), and the
  AIPW summand is formed from those weighted means. The summand itself is
  algebraically identical to theirs.

  INFERENCE, ALIGNED WITH THE FLOW TIER'S PRESENTATION (21_aipw_flow.do):
    - Act 1 (pooled, single ATE): percentile CI from a STRATIFIED CLUSTER
      bootstrap (resample countries), as before -- no by-type contrast to
      test here, so the row-level scheme below does not apply.
    - Act 2 (by resolution type) and its difference: each level's CI is now
      1.96*SE around the point estimate, where SE is the ANALYTIC
      influence-function standard error (sum of squared demeaned Eq.-3
      summands / N) -- the same construction behind the reference paper's
      own Table 2/Fig. 4 bands, not a bootstrap SD. A conventional t-test
      (b/se) with stars is reported for each level against zero, exactly as
      their Table 2. The bootstrap is reserved for the BETWEEN-TYPE
      DIFFERENCE (def - nd), using ROW-LEVEL resampling within
      control/non-default/default-linked pools -- the paper's own
      bootstrap device (split into a control pool and one pool per
      treatment type, `bsample` each, stack) -- which is a NATURAL fit
      here (not merely matched for parity): an onset row already IS one
      episode, so row resampling does not face the flow tier's own
      caveat about a single chronic country contributing many
      non-independent rows. Clogg et al. (1995)'s z, built from the same
      analytic SEs as the level bands, is reported alongside the
      bootstrap CI as the PERMISSIVE companion statistic (it assumes the
      two cells are independent, which they are not, since they share the
      tranquil control pool) -- the bootstrap CI governs where the two
      disagree.

  SAMPLE: restriction to ever-treated countries is mechanical — the probit drops
  countries with no variation in D. Act 1 = onset vs tranquil. Act 2 = the
  resolution split shown as TWO level IRFs (non-default vs tranquil, default vs
  tranquil) — the doubly-robust analog of the Table 2 nd/def coefficients; their
  vertical gap is the extra cost of default. Mirrors the IPW Act-2 figure (fig8).

  Run AFTER 01e_predictors. (08_ipw_lp.do, the plain-IPW predecessor to this
  doubly-robust estimator, was removed project-wide once AIPW covered the
  same question with a weaker robustness requirement -- see METHODOLOGY.md.)
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* Probit BASELINE = outcome baseline ($ctrl_core), so pmodel(cx cz) = $ctrl_core + Z
* while omodel = $ctrl_core: same baseline in both stages (their $convar design).
local cx    $ctrl_core
* Act 1 predictors Z1: single global push (fed funds) + contagion + proneness.
local cz    l_fedfunds l_reg_crisis_share past_onsets
* Act 2 predictors Z2 (resolution-type): matches the flow tier's adopted
* cz_recency in spirit (a country-year-specific spatial-lag contagion
* measure in place of a coarse regional share; a genuine recency clock in
* place of a raw past-event count that conflates recency and frequency) --
* see 21_aipw_flow.do's "SECOND PREDICTOR CHANGE" section for that history.
* Both terms are DEFAULT-LINKED-SPECIFIC here (l_contagion_dist_def,
* years_since_def_onset), not the generic versions used in cz (Act 1,
* pooled). The justification is economic, not merely statistical: a
* predictor for "risk of a DEFAULT-linked onset specifically" should
* measure default-linked distress/recency, not spread-crisis distress in
* general -- a generic measure blends illiquidity/market pressure with
* actual default propensity and asks it to predict the narrower one.
* Corroborated (not driven) by 08c_first_stage_table.do's diagnostics: the
* default arm's classification power was materially weaker under the
* generic combination (chi2(pred) p=.101, roccomp p=.525 -- barely
* distinguishable from controls-only) than under default-linked-only
* predictors. See 08c's header for the full argument, including the honest
* caveat that the reference paper's own $instrument is NOT tailored per
* column the way cz_def is here. cz (Act 1, pooled) is unchanged: there is
* no resolution type to be specific about in a pooled spec.
local cz_def l_fedfunds l_contagion_dist_def years_since_def_onset

* ── REPRODUCIBILITY: seed the bootstrap ────────────────────────────────────
* Every CI in this file comes from `bsample', which draws at random. Without a
* seed the intervals move between runs of identical code, and on thin cells
* that is not a rounding issue -- in 13d the non-default nexus gap at Year 1
* returned [0.05, 7.37] in one run and [-0.088, 7.686] in the next, opposite
* verdicts at the 5% line from the same data. Seeding makes the reported
* interval a property of the estimator rather than of the session. The value
* is arbitrary and was not chosen by inspecting results.
set seed 20260819

local nboot = 500      // bootstrap reps; raise to 1000+ for the final run

* Storage: rows h=0..4
foreach m in b se lo hi {
    matrix A1_`m'     = J(5,1,.)   // Act 1: crisis cost (all onsets)
    matrix A2nd_`m'   = J(5,1,.)   // Act 2: non-default cost vs tranquil
    matrix A2def_`m'  = J(5,1,.)   // Act 2: default-linked cost vs tranquil
    matrix A2diff_`m' = J(5,1,.)   // Act 2: extra cost of default (def - nd), row bootstrap
}
matrix A2diff_z = J(5,1,.)   // Act 2 difference: Clogg et al. (1995) z (analytic SEs)
matrix A2diff_p = J(5,1,.)   // Act 2 difference: p-value of the Clogg z

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — hand-coded AIPW point estimate (Eqs. 1-3), returns r(theta), r(N)
*   _aipw dyvar Dvar [if], omodel(X) pmodel(X Z) [fe(idvar)]
*   `fe(idvar)' adds country FE to the OUTCOME regression only (Eq. 1). The probit
*   (Eq. 2) is always pooled. Trims phat to [.01,.99] to bound the weights.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    * omodel = Eq.(1) RHS covariates (X); pmodel = Eq.(2) probit RHS (X and Z).
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'           // common sample for Eq.1 & Eq.2

    tempvar xb m0 m1 ps summ iwt
    * ORDER MATTERS: the propensity comes FIRST, because the outcome regression is
    * WEIGHTED by it. The reference paper's estimator is IPWRA, not plain AIPW —
    * their conditional means come from
    *     reg g_h dum g_0 $convar [pweight=invwt]
    * i.e. the regression that produces mu0/mu1 is itself inverse-probability
    * weighted. Estimating that regression unweighted (as this program used to)
    * gives a different, also doubly-robust, estimator — but not theirs.

    * --- Eq. (2): POOLED propensity probit -> phat (trimmed) ---
    * Pooled (no country FE): with ~20 events, country dummies in a probit
    * separate/overfit (perfectly predict never-crisis countries, AUROC->1),
    * collapsing the propensity — the standard rare-event choice keeps it pooled.
    * Country FE enter the OUTCOME regression (Eq. 1) only, matching the paper's
    * outcome design; no year FE.
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    * Winsorise rather than drop, so the cell keeps its observations. The paper
    * does not trim at all; on 21 default events an untrimmed 1/p can explode.
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'

    * --- IPW weight, the paper's UNSTABILISED form invwt = a/p + (1-a)/(1-p) ---
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'

    * --- Eq. (1): IPW-WEIGHTED outcome regression -> conditional means m1, m0 ---
    * `fe' names the country-FE variable (cid on real data; a fresh resample id
    * under the cluster bootstrap so duplicated countries are distinct FE groups).
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' [pweight=`iwt'] if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' [pweight=`iwt'] if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'   // set D=0
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'   // set D=1

    * --- Eq. (3): AIPW summand (Asonuma et al. exact form) ---
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    local th = r(mean)
    local nn = r(N)

    * Analytic (unclustered) influence-function SE -- the construction behind
    * the reference paper's own Table 2/Fig. 4 bands: SE = sqrt(mean((summand
    * - theta)^2) / N). Matches 21_aipw_flow.do's _aipw exactly.
    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta = `th'
    return scalar N     = `nn'
    return scalar se    = `sean'
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — bootstrap stratum: does this COUNTRY ever carry the treatment?
*
* The paper resamples the control pool and EACH TREATED TYPE separately, so the
* number of treated observations of every type is held fixed across draws:
*     drop if dum1==1 | dum2==1 | dum3==1  /  bsample          (controls)
*     keep if dum`s'==1  /  bsample  /  append                 (each type)
*
* A plain `bsample, cluster(cid)' does not do that: it resamples countries and
* takes whatever treatment mix falls out. With 12 default episodes a draw can
* easily land with three or four of them, or none, and those draws are either
* useless or fail outright — which is why the draw counts came in below nboot.
*
* This builds the stratum so `bsample, cluster(cid) strata()' preserves the count
* of ever-treated COUNTRIES in every draw. Country-level, so it is constant
* within cluster as bsample requires. We stratify on countries rather than on
* observations (their choice) because clustering by country is what keeps the
* within-country dependence in the standard error; theirs ignores it.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _mkstrat
program define _mkstrat
    syntax varlist(min=1 max=2) [if], GENerate(name)
    marksample touse
    capture drop `generate'
    tempvar t1 t2
    local d1 : word 1 of `varlist'
    local d2 : word 2 of `varlist'
    quietly gen byte `t1' = `d1' if `touse'
    bysort cid: egen byte `generate' = max(`t1')
    quietly replace `generate' = 0 if missing(`generate')
    if "`d2'" != "" {
        * two treatments (the paired difference): a 4-level stratum so BOTH
        * treated counts are preserved on the same resample.
        quietly gen byte `t2' = `d2' if `touse'
        tempvar s2
        bysort cid: egen byte `s2' = max(`t2')
        quietly replace `s2' = 0 if missing(`s2')
        quietly replace `generate' = `generate' + 2*`s2'
    }
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — paired bootstrap of the DIFFERENCE between two AIPW cells
*   (cell1 - cell2), ROW-LEVEL resampling within control/type1/type2 pools --
*   THE REFERENCE PAPER'S OWN BOOTSTRAP DEVICE (their file splits the data
*   into a control pool and one pool per treatment type, `bsample`s each
*   independently, and stacks: drop if dum1|dum2|dum3; bsample (controls);
*   keep if dum`s'==1; bsample; append). `bsample, strata()` does exactly
*   that in one command -- draws N_k rows WITH REPLACEMENT inside each
*   stratum, so every pool keeps its original size, which is what stops a
*   draw landing with too few default-linked rows. This is a NATURAL match
*   for the onset tier, not merely matched for parity with the flow tier's
*   own file: an onset row already is one episode (unlike a flow-coded
*   continuation row), so row resampling does not face the flow tier's own
*   caveat about one chronic country contributing many non-independent rows.
*   Also returns each cell's own analytic SE (a1/a2), so the console table
*   and the Clogg z below are built from the SAME estimates as the
*   bootstrap, not a separate run.
*   _aipwpair, y() d1() if1() d2() if2() omod() pz() reps()
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipwpair
program define _aipwpair, rclass
    syntax , Y(string) D1(string) IF1(string) D2(string) IF2(string) ///
             OMOD(string) PZ(string) REPS(integer)
    capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b1 = r(theta)
    local a1 = r(se)
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2 = r(theta)
    local a2 = r(se)
    local dh = `b1' - `b2'

    capture drop _pool
    * Parentheses matter: `if1'/`if2' are compound expressions and Stata gives
    * & and | equal precedence, left to right.
    quietly gen byte _pool = 0 if (`if1') | (`if2')
    quietly replace _pool = 1 if `d1' == 1
    quietly replace _pool = 2 if `d2' == 1

    tempname pf
    tempfile bf
    quietly postfile `pf' double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
            local t1 = cond(_rc==0, r(theta), .)
            capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
            local t2 = cond(_rc==0, r(theta), .)
            if !missing(`t1') & !missing(`t2') quietly post `pf' (`t1' - `t2')
        restore
    }
    quietly postclose `pf'
    capture drop _pool
    local se = .
    local lo = .
    local hi = .
    local nd = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(diff)
        local nd = r(N)
        if `nd' >= 50 {
            quietly summarize diff
            local se = r(sd)
            _pctile diff, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
        }
    restore
    return scalar ok = 1
    return scalar dh = `dh'
    return scalar b1 = `b1'
    return scalar b2 = `b2'
    return scalar a1 = `a1'
    return scalar a2 = `a2'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
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
    capture _aipw dy_`h' onset_all if sample==1, omodel($ctrl_core) pmodel(`cx' `cz') fe(cid)
    if _rc {
        di as error "h=" `h' ": Act 1 point estimate failed, rc=" _rc
        continue
    }
    local pt = r(theta)
    matrix A1_b[`row',1] = `pt'

    * STRATIFIED cluster bootstrap: resample countries within ever-treated and
    * never-treated strata, so every draw carries the same number of treated
    * countries. See _mkstrat above for why an unstratified draw is unusable on
    * a cell this thin.
    _mkstrat onset_all if sample==1, generate(_strat)
    tempname pf
    tempfile bf
    quietly postfile `pf' double theta using "`bf'", replace
    forvalues b = 1/`nboot' {
        preserve
            capture drop _bid
            bsample, cluster(cid) strata(_strat) idcluster(_bid)
            capture _aipw dy_`h' onset_all if sample==1, omodel($ctrl_core) pmodel(`cx' `cz') fe(_bid)
            if _rc == 0 quietly post `pf' (r(theta))
        restore
    }
    quietly postclose `pf'
    capture drop _strat

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

    di "h=" `h'+1 "   " %7.3f `pt' "   " ///
       %7.3f A1_se[`row',1] "   [" %7.3f A1_lo[`row',1] ", " %7.3f A1_hi[`row',1] "]" ///
       "   (" `ndraw' "/`nboot' draws)"
}

* ══════════════════════════════════════════════════════════════════════════
* ACT 2 — AIPW output cost by resolution, levels + the DIFFERENCE (def - nd),
* aligned with the flow tier's presentation (21_aipw_flow.do Section 3):
*   - Each level's CI is 1.96*analytic-SE (the paper's own Table 2/Fig. 4
*     construction), not a bootstrap band; a conventional t-test (b/se_a)
*     with stars is shown for each level against zero, exactly as their
*     Table 2. Mirrors the IPW Act-2 figure (fig8): one line for
*     non-default onsets, one for default-linked, each vs tranquil with
*     country FE + the full (X,Z) probit -- the OTHER resolution type
*     dropped so the control group is clean tranquil years.
*   - The DIFFERENCE (def - nd) is bootstrapped directly (row-level, within
*     control/nd/def pools -- the paper's own device, see _aipwpair's
*     header) so it gets a proper percentile CI, the quantity the paper
*     actually tests rather than eyeballing the two level bands.
*   - Clogg et al. (1995)'s z, built from the SAME analytic SEs as the
*     level bands, is reported alongside the bootstrap CI as the
*     PERMISSIVE companion statistic (assumes independence, which does not
*     hold since both cells share the tranquil control pool) -- the
*     bootstrap CI governs where the two disagree.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 2 — AIPW output cost by resolution, levels + difference ==="
di as result "h   ND (se_a)         DEF (se_a)         def-nd   [95% boot CI]      Clogg z    p     draws"
di as result "    se_a = analytic influence-function SE (the paper's own construction)."
di as result "    ND/DEF stars are the paper's 'conventional t-test' vs zero (b/se_a): * p<.10 ** p<.05 *** p<.01."
di as result "    def-nd's own * marks the bootstrap CI excluding 0 -- the conservative, governing test for the difference."

forvalues h = 0/4 {
    local row = `h' + 1

    _aipwpair, y(dy_`h') ///
        d1(onset_def) if1(sample==1 & onset_nd==0) ///
        d2(onset_nd)  if2(sample==1 & onset_def==0) ///
        omod($ctrl_core) pz(`cx' `cz_def') reps(`nboot')

    if !r(ok) {
        di as error "  h=" `h'+1 ": estimate failed (cell too thin)."
        continue
    }

    local B1 = r(b1)   // default-linked ATE
    local B2 = r(b2)   // non-default ATE
    local A1 = r(a1)   // analytic SE, default-linked
    local A2 = r(a2)   // analytic SE, non-default
    local DH = r(dh)
    local SE = r(se)
    local LO = r(lo)
    local HI = r(hi)
    local ND = r(nd)

    matrix A2def_b[`row',1]  = `B1'
    matrix A2nd_b[`row',1]   = `B2'
    matrix A2def_se[`row',1] = `A1'
    matrix A2nd_se[`row',1]  = `A2'
    matrix A2diff_b[`row',1]  = `DH'
    matrix A2diff_se[`row',1] = `SE'
    matrix A2diff_lo[`row',1] = `LO'
    matrix A2diff_hi[`row',1] = `HI'

    * Level CIs = theta +/- 1.96*analytic SE, matching how the reference
    * paper bands its Fig. 4 (`up = irf + 1.96*se'); the bootstrap is
    * reserved for the difference only, exactly as they reserve it.
    matrix A2nd_lo[`row',1]  = `B2' - 1.96*`A2'
    matrix A2nd_hi[`row',1]  = `B2' + 1.96*`A2'
    matrix A2def_lo[`row',1] = `B1' - 1.96*`A1'
    matrix A2def_hi[`row',1] = `B1' + 1.96*`A1'

    * Clogg et al. (1995) z: (irf1 - irf2)/sqrt(se1^2+se2^2), the same
    * analytic SEs as the level bands. Permissive (assumes independence).
    local zz = .
    local pz = .
    if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
        local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
        local pz = 2*(1 - normal(abs(`zz')))
        matrix A2diff_z[`row',1] = `zz'
        matrix A2diff_p[`row',1] = `pz'
    }

    * Conventional t-test for each level vs zero (b / own analytic SE), as
    * the paper's Table 2 -- contrasted against the bootstrap/Clogg methods
    * used for the difference below, not applied to this test.
    local tnd  = cond(`A2'>0, `B2'/`A2', .)
    local pnd  = cond(!missing(`tnd'), 2*(1-normal(abs(`tnd'))), .)
    local sgnd = cond(missing(`pnd'), "", cond(`pnd'<.01,"***",cond(`pnd'<.05,"**",cond(`pnd'<.10,"*",""))))
    local tdef  = cond(`A1'>0, `B1'/`A1', .)
    local pdef  = cond(!missing(`tdef'), 2*(1-normal(abs(`tdef'))), .)
    local sgdef = cond(missing(`pdef'), "", cond(`pdef'<.01,"***",cond(`pdef'<.05,"**",cond(`pdef'<.10,"*",""))))

    local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
    di "  " %1.0f `h'+1 "  " %8.3f `B2' "`sgnd'" " (" %5.3f `A2' ")  " ///
       %8.3f `B1' "`sgdef'" " (" %5.3f `A1' ")  " %8.3f `DH' ///
       " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
       " " %7.3f `zz' " " %5.3f `pz' " " %4.0f `ND'
}

di as result _n "  * = bootstrap 95% percentile CI for the def-nd gap excludes zero."
di as result "  Clogg z assumes the two cells are independent (they share the tranquil"
di as result "  control pool), so it is the permissive statistic. Where it and the"
di as result "  bootstrap disagree, the bootstrap governs."

* ══════════════════════════════════════════════════════════════════════════
* EXPORT (difference) — extra-cost-of-default CSV (def - nd, paired bootstrap)
* ══════════════════════════════════════════════════════════════════════════
preserve
    clear
    tempname pfd
    tempfile difff
    postfile `pfd' horizon dhl se lo hi double cloggz double cloggp using "`difff'", replace
    post `pfd' (0) (0) (0) (0) (0) (.) (.)   // explicit baseline (displayed h=0), matching Asonuma et al.
    forvalues h = 0/4 {
        post `pfd' (`h'+1) (A2diff_b[`h'+1,1]) (A2diff_se[`h'+1,1]) ///
            (A2diff_lo[`h'+1,1]) (A2diff_hi[`h'+1,1]) ///
            (A2diff_z[`h'+1,1]) (A2diff_p[`h'+1,1])
    }
    postclose `pfd'
    use "`difff'", clear
    label var dhl "AIPW extra cost of default (def - nd, pp)"
    label var lo  "95% percentile CI lower (row bootstrap)"
    label var hi  "95% percentile CI upper (row bootstrap)"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    gen byte sig95 = (!missing(lo) & (lo>0 | hi<0))
    label var sig95 "Bootstrap CI excludes 0 (governing test)"
    order horizon dhl se lo hi sig95 cloggz cloggp
    export delimited "$tabs/aipw_act2_diff.csv", replace
    di as result "AIPW extra-cost-of-default CSV saved: $tabs/aipw_act2_diff.csv"
restore

* representative gap (displayed h=3, near the trough) for the Act-2 figure annotation
local g2  = A2diff_b[3,1]
local g2l = A2diff_lo[3,1]
local g2h = A2diff_hi[3,1]
local gapnote ""
if !missing(`g2') & !missing(`g2l') {
    local g2s  : di %4.1f `g2'
    local g2ls : di %4.1f `g2l'
    local g2hs : di %4.1f `g2h'
    local gapnote "Extra cost of default at h=3: `g2s' pp [`g2ls', `g2hs'] (row bootstrap, paper's device)."
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
        post `pfx' ("`sname'") (0) (0) (0) (0) (0)   // explicit baseline, matching Asonuma et al.
        forvalues h = 0/4 {
            post `pfx' ("`sname'") (`h'+1) ///
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
        xlabel(0(1)5) ytitle("Cumulative real GDP change (pp)", size(small)) ///
        xtitle("Year (Year 1 = crisis year)", size(small)) ///
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
        xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year (Year 1 = crisis year)", size(medsmall)) ///
        ytitle("Cumulative change in log real GDP (pp)", size(medsmall)) ///
        title("AIPW output cost by resolution", size(medium) color(navy)) ///
        subtitle("Doubly-robust (Asonuma et al. Eq. 3), each vs tranquil.", size(small)) ///
        legend(order(3 "Non-default" 4 "Default-linked") size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    * Note text (kept as source comment, no longer rendered on the figure --
    * the legend, previously overlapping the plot at ring(0) pos(7), now
    * takes the bottom position this note used to occupy). `gapnote' is a
    * runtime local (its rendered value can't be embedded in a static
    * comment), so the template is shown instead:
    * "Two AIPW level IRFs; shaded = 1.96*analytic SE band (the paper's own construction). Gap = extra cost of default, bootstrapped directly (row-level). <gapnote>"
    graph export "$figs/fig_aipw_act2.pdf", replace
    di as result "Figure saved: fig_aipw_act2.pdf"
restore

di as result _n "08b_aipw.do complete."
di as result "fig_aipw.pdf = Act 1 crisis cost; fig_aipw_act2.pdf = the two-line"
di as result "resolution split (non-default vs default, each vs tranquil) — the AIPW"
di as result "twin of the IPW fig8. Compare the nd/def AIPW levels to fig8's IPW lines"
di as result "and to the Table 2 OLS coefficients (similar ordering => selection is not"
di as result "driving the resolution gap; their Fig C1 logic)."
di as result "aipw_act2_diff.csv bootstraps the def-nd GAP directly (extra cost of default),"
di as result "matching the paper's inference for the between-type contrast; sig95==1 means the"
di as result "gap's 95% CI excludes 0 at that horizon."
