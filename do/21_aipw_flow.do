/*===========================================================================
  21_AIPW_FLOW.DO
  Asonuma, Chamon, Erce & Sasahara (2024) Eqs. (1)-(3), applied to the FLOW
  treatment definition built in 18_transforms.do and estimated by OLS in
  20_lp_flow.do.

  WHY THIS FILE EXISTS
  --------------------
  20_lp_flow.do finds the def-nd output gap significant at Years 1-4 under flow
  coding (in_crisis = 1 in EVERY year of an episode; 234 treated country-years),
  but that tier is OLS only, so it is not corrected for selection.

  08b_aipw.do runs the reference paper's estimator on ONSET coding (one treated
  row per episode). It is NOT modified here. Running the SAME estimator on the
  flow treatment isolates the effect of the treatment definition while holding
  the estimator fixed: if the two agree, the gap is not an artifact of how the
  treatment is coded.

  THE THREE STEPS, IN THE PAPER'S ORDER
  ------------------------------------
  Eq. (1)  Outcome regression with country FE, INVERSE-PROBABILITY WEIGHTED
           (their `reg g_h dum g_0 $convar [pweight=invwt]'), giving the
           conditional means m1 and m0. Their estimator is IPWRA.
  Eq. (2)  Propensity probit Pr(D=1 | X, Z) -> phat. Pooled, no country FE.
  Eq. (3)  Lambda_h = (1/N) sum { [D*y/p - (1-D)*y/(1-p)]
                                 - (D-p)/(p(1-p)) * [(1-p)*m1 + p*m0] }

  D is in_crisis_nd or in_crisis_def. Two type-level lines against tranquil
  years; the DIFFERENCE (def - nd) is the extra output cost of default and is
  bootstrapped directly on paired cluster resamples, which is the paper's own
  inference for a between-type contrast.

  TWO DEPARTURES FROM 20_lp_flow.do, BOTH TO MATCH THE REFERENCE PAPER
  --------------------------------------------------------------------
  * NO YEAR FIXED EFFECTS. Their $convar carries c1-c74 (country dummies) and
    no year dummies. The project's rule is that two-stage estimators drop year
    FE to stay like-for-like with them (METHODOLOGY.md section 2); 20 keeps them
    because it is single-stage.
  * `reg' RATHER THAN `xtscc'. Inference here comes from the paired cluster
    bootstrap, not from analytic SEs, so the DK correction has nothing to do;
    and Eq. (1) feeds conditional means into Eq. (3) rather than being read as
    a coefficient in its own right. It is also mechanical: xtscc does not accept
    pweight, and Eq. (1) is weighted.

  WHERE THIS DEPARTS FROM THEIR CODE — THREE THINGS, ALL DELIBERATE
  ------------------------------------------------------------------
  The estimator algebra is IDENTICAL to theirs. Their
      iptw   = (2a-1)*g*invwt,  invwt = a/p + (1-a)/(1-p)
  is g/p when treated and -g/(1-p) when control, i.e. the first bracket of
  Eq. (3); and their
      mdiff1 = -(a-p)*mu1/p - (a-p)*mu0/(1-p)
  factors to -(a-p)[mu1/p + mu0/(1-p)], which is the adjustment term as written
  above. The mu0/mu1 shift-by-_b[dum] construction is the same too. What differs
  is three implementation choices:

  1. TRIMMING. They do not trim pihat at all. This file winsorises to
     [0.01, 0.99], following 08b, because an untrimmed 1/p explodes on thin
     cells — and under flow coding continuation rows are precisely the ones at
     risk of scores near 1. Section 1 reports how many rows this touches.
  2. BOOTSTRAP UNIT — NO LONGER A DEPARTURE. The baseline now uses their
     scheme: row-level `bsample' within treatment-type strata, no cluster(),
     country FE on cid so duplicated rows of a country share it. Section 3b
     re-runs the identical estimator resampling whole COUNTRIES and prints the
     two intervals side by side, because the choice matters here in a way it
     does not in their paper: under onset coding a treated row is close to an
     independent episode (194 across 76 countries), while under flow coding a
     treated row is a crisis-YEAR and 234 of them carry 61 episodes in 52
     countries, 29 of them Venezuela's. Point estimates are identical across
     the two modes, so the comparison isolates the inference choice. Expect the
     row intervals to be narrower; where the two verdicts differ, report both.
  3. YEAR FE. Neither has them; see below. (This is agreement, not a departure,
     but it is the thing most likely to be misread as one given 20 has them.)

  WHICH CONTROL SET GOES WHERE — AND WHY THEY DIFFER
  --------------------------------------------------
  Eq. (1) outcome model:  $ctrl_flow          (episode-dated, epc_*)
  Eq. (2) propensity:     $ctrl_core_flowbase (row-dated) via `cx_active' — see
                          "ADOPTED FLOW-TIER CORE CONTROL SET" below.

  The row-dated/episode-dated split described next is not an oversight. epc_X is built as
      epc_X = cond(in_crisis==1, X at episode entry, X at own t-1)
  so its VALUE depends on the treatment status of the row. Tranquil and onset
  rows have epc_X == X by construction; continuation rows generally do not.
  That means "epc_X differs from X" is a direct readout of "this row is a
  continuation row" — a covariate that mechanically encodes the treatment.
  Putting it on the RHS of the probit would separate. Section 1 below counts
  the rows where they differ, so the argument is demonstrated in the log rather
  than asserted here.

  $ctrl_core is the only conditioning set defined identically in both arms, so
  it is the one the probit can use. It is not clean either, and that must be
  said plainly: for a continuation row, last year's debt / current account /
  credit are partly OUTCOMES of the crisis the row is already in. Under onset
  coding every treated row is an entry row and the controls are predetermined
  by construction; under flow coding they are not. Unconfoundedness given X is
  therefore a materially stronger assumption here than in 08b.

  The paper's Eq. (2) has no lagged treatment term and none is added, so the
  weights stay bounded. Whether the propensity nonetheless separates treated
  from control is an empirical question, and Section 1 answers it before any
  coefficient is reported.

  DIAGNOSTIC HISTORY — SECTIONS 1a/1b/1c/1d/1e/1f, ALL REMOVED
  -----------------------------------------------------------------------------
  This file previously carried an extended diagnostic sequence establishing,
  step by step, that the flow AIPW's pooled propensity model near-separates on
  continuation rows (mean p(treated) 0.591, 5 rows p>0.99 in the def arm), that
  entry-dating the covariates barely helped (0.591 -> 0.557), that dropping
  l_ca/l_debt or past_def_onsets individually helped only partially (0.591 ->
  0.523 / 0.321), and that the fix actually adopted — restricting the probit's
  FITTING sample to tranquil+onset rows (continuation==0), then extrapolating
  phat to continuation rows — closed the gap far more (0.591 -> 0.147, p>0.99
  5 -> 0) than any single-variable drop. A follow-on predictor swap
  (past_def_onsets -> years_since_def_onset, motivated by 24_aipw_channels_flow.do's
  own weight-concentration diagnostic) was adopted on similar grounds: a
  better-motivated predictor, not a fix for the def arm's wide SEs, which
  remains open. Both decisions — the continuation==0 fitting-sample
  restriction and the years_since_def_onset predictor — are now settled and
  are the estimator this file always runs, with no toggle back to the
  alternatives. The diagnostic code establishing them is not kept in this
  file; see git history if it is ever needed again. `local cz_def' and the
  frozen `cx' (row-dated on the original $ctrl_core) are gone with it — the
  file now defines only the active predictor/control set, `cz_recency'/`cx_active'.

  SECOND PREDICTOR CHANGE — l_reg_crisis_share REPLACED BY l_contagion_dist
  IN cz_recency
  -----------------------------------------------------------------------------
  l_reg_crisis_share (leave-one-out share of OTHER same-region countries with
  an onset) treats every country in a coarse geographic region as an equally
  weighted neighbor and ignores every country outside it. l_contagion_dist
  (17_predictors.do) relaxes that: a distance-weighted sum of every OTHER
  country's onset that year, using the reference paper's own contagion
  construction (their Section 4.1) applied to onset_all instead of a
  restructuring dummy, with great-circle distance between capital cities
  (CEPII geo_cepii) in place of a fixed regional grouping. It is a
  contemporaneous, country-year-specific spatial lag, not a slow-moving
  characteristic of country i itself -- the same property that made
  years_since_def_onset preferable to past_def_onsets. `cz_recency' below now
  carries l_contagion_dist in place of l_reg_crisis_share.

  COUNTRY FE IN THE PROBIT — TESTED AND REJECTED
  -----------------------------------------------------------------------------
  Their first-stage probit carries country dummies directly ($cf = c1-c74,
  noconstant) alongside $convar and $instrument. Adding i.cid to this file's
  Eq. (2) was tested (on top of cx_active/cz_recency, continuation==0) and
  does not survive on this project's much smaller panel: the nd arm lost
  more than half its identifying countries, and the def arm collapsed to
  complete separation — log-likelihood driven to exactly 0, every
  coefficient exploding with no standard error. 21b_first_stage_table_flow.do's
  pooled, clustered-SE probit (the SAME cx_active/cz_recency, no country FE)
  converges cleanly on the identical sample, confirming the failure is
  specifically about adding country dummies as regressors, not about the
  control/predictor set itself. Eq. (2) therefore stays pooled, exactly as
  it always has been; the diagnostic that established this is not kept in
  this file (see git history if it is ever needed again).

  Outputs
  -------
    "$tabs/table10_aipw_flow.rtf"    the two type lines + the difference
    "$tabs/aipw_flow_diff.csv"       def-nd gap, bootstrap CI, Clogg z
    "$clean/irf_aipw_flow_nd.dta"    non-default AIPW IRF
    "$clean/irf_aipw_flow_def.dta"   default-linked AIPW IRF
    "$figs/fig10_aipw_flow.pdf/.png" the two-line resolution split

  INFERENCE — SPLIT EXACTLY AS THE REFERENCE PAPER SPLITS IT
  ----------------------------------------------------------
  LEVELS (each type vs tranquil): analytic influence-function SE, built inside
  the AIPW loop the way they build it —
      sum dr1 ; gen Isq = (dr1 - mean)^2 ; sum Isq ; se = sqrt(r(mean)/r(N))
  with bands theta +/- 1.96*se. This is what stands behind their Fig. 4, and
  the same se feeds their Clogg z.

  DIFFERENCE (def - nd): 1000-draw stratified cluster bootstrap, percentile CI.
  Their bootstrap file computes each type per draw and the post-processing file
  uses those draws ONLY to form the between-type contrasts and take
  centile(2.5 97.5). Nothing about their level CIs touches the bootstrap.

  ONE THING TO KNOW ABOUT THEIR ANALYTIC SE: it is UNCLUSTERED — sqrt of the
  mean squared influence function over ROWS. Their own code computes a
  clustered version (`reg dr1 ATE_IPWRA, nocons cluster(wdicode)') and then does
  not use it. With repeated country-year observations that understates the
  standard error, and flow coding makes it worse, because 234 treated rows still
  carry only 61 episodes in 52 countries. The per-type BOOTSTRAP SD is therefore
  printed beside the analytic SE at every horizon. It is a diagnostic and is not
  used for any interval: if the two are close the concern is mild; if the
  bootstrap SD is much larger, the paper's construction is optimistic on this
  sample and the write-up must say so rather than quoting the narrower band.

  08b_aipw.do takes the other route — bootstrap for the levels too — as a
  deliberate departure from the paper on thin cells. Two files, two documented
  choices: 08b conservative, this one faithful.

  SELF-CONTAINED: this file reads nothing but $clean/panel_lp.dta. It does not
  depend on 08b_aipw.do or 20_lp_flow.do having been run, and does not compare
  itself to either; any such comparison belongs in the write-up, not the code.

  RUNTIME: the paired bootstrap refits Eq. (1)-(3) twice per draw per horizon,
  and Section 3b repeats the whole loop under the other resampling unit. At
  nboot=1000 that is ~20,000 probit+regression pairs. Expect several minutes.
  If that is painful, drop the Section 3b comparison to 500 draws — never the
  Section 3 baseline.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow years_since_def_onset l_contagion_dist {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel_lp.dta — re-run 18_transforms.do first."
        exit 111
    }
}
if "$ctrl_flow" == "" {
    di as error "  ** \$ctrl_flow not set — re-run 18_transforms.do, then this file."
    exit 111
}

* EXPLORATORY: flow_ctrl_variant -- 0 = adopted flow-tier baseline (default);
* 1 = adopted core PLUS the reference paper's own additional predictors,
* l_imf, since exchange2 is now already inside the adopted core (18_transforms.do's "ALTERNATE FLOW CONTROL SET").
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
* `cx_active' is the row-dated propensity control set: the flow tier's
* adopted core control set (l_banking_crisis, l_lninfl, exchange2 in place
* of l_banking_duration, l_hyperinfl, l_ca/tot_chg), or the alternate (core
* + l_imf) if toggled.
local cx_active = cond(`flow_ctrl_variant'==1, "$ctrl_core_flowplus", "$ctrl_core_flowbase")  // Eq. (2) baseline: row-dated, ADOPTED set
local com    = cond(`flow_ctrl_variant'==1, "$ctrl_flow_flowplus", "$ctrl_flow")    // Eq. (1) controls: episode-dated (already the adopted set via 18_transforms.do)
* cz_recency: the adopted predictor set -- l_contagion_dist in place of
* l_reg_crisis_share, years_since_def_onset in place of past_def_onsets.
local cz_recency l_fedfunds l_contagion_dist years_since_def_onset

set seed 20260819
local nboot = 1000

foreach m in b se lo hi {
    matrix Fnd_`m'   = J(5,1,.)
    matrix Fdef_`m'  = J(5,1,.)
    matrix Fdiff_`m' = J(5,1,.)
}
matrix Fdiff_z = J(5,1,.)
matrix Fdiff_p = J(5,1,.)
foreach m in b lo hi {
    matrix Cdiff_`m' = J(5,1,.)     // country-cluster bootstrap, comparison only
}
matrix Cdiff_n = J(5,1,.)

* ══════════════════════════════════════════════════════════════════════════
* PROGRAMS — copied from 08b_aipw.do so this file is self-contained
* (house style already duplicates _critvals/_pval across files)
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'

    tempvar xb m0 m1 ps summ iwt
    * ORDER FOLLOWS THEIR CODE: Eq. (2) runs FIRST, because Eq. (1) is weighted
    * by the inverse propensity. Their replication file is unambiguous —
    *     reg g_h dum g_0 $convar [pweight=invwt], cluster(wdicode) noconstant
    * so the estimator is IPWRA, not plain AIPW. (Their TEXT says Eq. (1) is
    * estimated "by OLS" and prints it unweighted; the code is what produced
    * Table 2 and Fig. 4, so the code governs here. Both forms are doubly
    * robust; this is a choice between two valid estimators, and this file
    * exists to reproduce theirs.)
    * PROPENSITY FIT SAMPLE EXCLUDES CONTINUATION ROWS -- see header "SECTION
    * 1f / ESTIMATOR CHANGE". Eq. (2) is fit only on tranquil and onset rows
    * (continuation==0); the fitted model is then EXTRAPOLATED to every row of
    * `touse', continuation rows included, because Eq. (1) and Eq. (3) still
    * use the FULL flow-coded `D' (continuation years remain treated). Only the
    * first stage's FITTING sample is restricted -- the treatment definition
    * used everywhere else in this program is untouched.
    quietly probit `D' `pmodel' if `touse' & continuation==0
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01                & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'

    * Eq. (1) — IPW-weighted OLS with country FE -> conditional means m1, m0.
    * Their mu0/mu1 construction (predict, then shift by _b[dum]) is
    * algebraically what is done here.
    if "`fe'" != "" quietly reg `y' `D' `omodel' i.`fe' [pweight=`iwt'] if `touse'
    else            quietly reg `y' `D' `omodel'         [pweight=`iwt'] if `touse'
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']     if `touse'
    * Eq. (3) — the paper's exact algebraic form.
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    * Point estimate = the mean of the summand. Store BOTH returned scalars in
    * locals before anything else touches r(), then build the paper's analytic
    * standard error from the influence function:
    *     sum dr1 ; gen Isq = (dr1-mean)^2 ; sum Isq ; se = sqrt(r(mean)/r(N))
    * This is the SE behind the level CIs in their Fig. 4 and behind their
    * Clogg z. It is UNCLUSTERED — see the header.
    quietly summarize `summ' if `touse', meanonly
    local th = r(mean)
    local nn = r(N)

    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta = `th'
    return scalar N     = `nn'
    return scalar se    = `sean'
end

capture program drop _mkstrat
program define _mkstrat
    syntax varlist(min=1 max=2) [if], GENerate(name)
    marksample touse
    capture drop `generate'
    tempvar t1 t2 s2
    local d1 : word 1 of `varlist'
    local d2 : word 2 of `varlist'
    quietly gen byte `t1' = `d1' if `touse'
    bysort cid: egen byte `generate' = max(`t1')
    quietly replace `generate' = 0 if missing(`generate')
    if "`d2'" != "" {
        quietly gen byte `t2' = `d2' if `touse'
        bysort cid: egen byte `s2' = max(`t2')
        quietly replace `s2' = 0 if missing(`s2')
        quietly replace `generate' = `generate' + 2*`s2'
    }
end

* Paired bootstrap of BOTH cells and their difference on the SAME resample.
* Extended vs 08b: also returns the per-cell bootstrap SDs, so the Clogg z is
* built from the same draws as the percentile CI rather than a separate run.
capture program drop _aipwpairflow
program define _aipwpairflow, rclass
    syntax , Y(string) D1(string) IF1(string) D2(string) IF2(string) ///
             OMOD(string) PZ(string) REPS(integer) [BOOT(string)]
    if "`boot'" == "" local boot row
    if !inlist("`boot'","row","cluster") {
        di as error "  ** boot() must be row or cluster"
        exit 198
    }
    capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b1  = r(theta)
    local a1  = r(se)
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2  = r(theta)
    local a2  = r(se)
    local dh  = `b1' - `b2'

    * ── Build the resampling strata ─────────────────────────────────────────
    * boot(row)     — THE REFERENCE PAPER'S SCHEME. Their bootstrap file splits
    *                 the data into a control pool and one pool per treatment
    *                 type, `bsample's each independently, and stacks:
    *                     drop if dum1|dum2|dum3 ; bsample                (controls)
    *                     keep if dum`s'==1 ; bsample ; append            (each type)
    *                 `bsample, strata()' does exactly that in one command: it
    *                 draws N_k rows WITH REPLACEMENT inside each stratum, so
    *                 every pool keeps its original size. That is the device
    *                 that stops a draw landing with no default-linked rows.
    *                 Row-level strata, no cluster(), no idcluster — as theirs.
    * boot(cluster) — resamples whole COUNTRIES within 4 country-level strata,
    *                 with idcluster so a country drawn twice becomes two
    *                 distinct pseudo-countries with separate fixed effects.
    if "`boot'" == "row" {
        capture drop _pool
        * Parentheses matter: `if1'/`if2' are compound expressions and Stata
        * gives & and | equal precedence, left to right.
        quietly gen byte _pool = 0 if (`if1') | (`if2')
        quietly replace _pool = 1 if `d1' == 1
        quietly replace _pool = 2 if `d2' == 1
    }
    else {
        _mkstrat `d1' `d2', generate(_strat)
    }

    tempname pf
    tempfile bf
    quietly postfile `pf' double t1 double t2 double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            if "`boot'" == "row" {
                * Restrict to the estimation universe BEFORE resampling, which
                * is what their `keep if sampmax == 1' does ahead of the
                * split-bsample-stack. Rows outside it would otherwise be
                * resampled into pool 0 and then dropped by the `if', adding
                * nothing but changing the pool size.
                quietly keep if !missing(_pool)
                * Row resampling within pool; country FE stay on cid, so
                * duplicated rows of a country share its fixed effect exactly
                * as they share a c1-c74 dummy in the reference code.
                quietly bsample, strata(_pool)
                capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
                local v1 = cond(_rc==0, r(theta), .)
                capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
                local v2 = cond(_rc==0, r(theta), .)
            }
            else {
                capture drop _bid
                bsample, cluster(cid) strata(_strat) idcluster(_bid)
                capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(_bid)
                local v1 = cond(_rc==0, r(theta), .)
                capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(_bid)
                local v2 = cond(_rc==0, r(theta), .)
            }
            if !missing(`v1') & !missing(`v2') quietly post `pf' (`v1') (`v2') (`v1'-`v2')
        restore
    }
    quietly postclose `pf'
    capture drop _strat _pool

    local se = .
    local lo = .
    local hi = .
    local s1 = .
    local s2 = .
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
            quietly summarize t1
            local s1 = r(sd)
            quietly summarize t2
            local s2 = r(sd)
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
    return scalar s1 = `s1'
    return scalar s2 = `s2'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* 1. TREATMENT/CONTROL COUNTS
*
* The propensity-model diagnostic sequence that used to run here (separation
* checks, the continuation==0 fitting-sample restriction, the predictor swap)
* is settled and no longer kept in this file -- see the header's "DIAGNOSTIC
* HISTORY" note and git history if it is ever needed again. What remains is
* the basic sample count, so the estimation loop below is read against a
* known N rather than a number pulled from nowhere.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "1. SAMPLE COUNTS"
di as result "════════════════════════════════════════════════════════════"

quietly count if in_crisis_nd==1
di as result "  treated, non-default (flow):     " %4.0f r(N) "   (expect 113)"
quietly count if in_crisis_def==1
di as result "  treated, default-linked (flow):  " %4.0f r(N) "   (expect 121)"
quietly count if sample_flow==1 & in_crisis==0
di as result "  tranquil controls:               " %4.0f r(N)

* ══════════════════════════════════════════════════════════════════════════
* 2. HOW MUCH WORK IS THE AUGMENTATION DOING?
*
* Eq. (3) is the IPW estimator plus an augmentation term:
*     Lambda = [ D*y/p - (1-D)*y/(1-p) ]                        <- IPW piece
*            - (D-p)/(p(1-p)) * [ (1-p)*m1 + p*m0 ]             <- augmentation
* Reporting the two separately is a real specification diagnostic, not
* bookkeeping. If the augmentation is small relative to the IPW piece, the
* outcome model and the propensity model broadly agree and the doubly-robust
* correction is doing little. If it is large, the two models disagree about
* the counterfactual and the estimate is leaning heavily on Eq. (1) being
* right — which is exactly the assumption flow coding strains.
*
* Self-contained: this reads nothing from any other file.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "2. AIPW DECOMPOSITION — IPW piece vs augmentation (Year 1)"
di as result "════════════════════════════════════════════════════════════"

foreach s2 in nd def {
    if "`s2'" == "nd"  local rival in_crisis_def
    else               local rival in_crisis_nd

    capture drop _ps2 _m0 _m1 _xb2 _ipwterm _augterm
    * _aipw's fitting-sample restriction (continuation==0, extrapolated to the
    * full sample) applied directly here, so this decomposition is computed
    * under the SAME propensity construction Section 3 uses.
    quietly probit in_crisis_`s2' `cx_active' `cz_recency' if sample_flow==1 & `rival'==0 & continuation==0
    quietly predict double _ps2 if sample_flow==1 & `rival'==0, pr
    quietly replace _ps2 = .01 if _ps2 < .01 & !missing(_ps2)
    quietly replace _ps2 = .99 if _ps2 > .99 & !missing(_ps2)

    quietly reg dy_0 in_crisis_`s2' `com' i.cid if sample_flow==1 & `rival'==0 & !missing(_ps2)
    quietly predict double _xb2 if e(sample), xb
    quietly gen double _m0 = _xb2 - _b[in_crisis_`s2']*in_crisis_`s2'
    quietly gen double _m1 = _m0 + _b[in_crisis_`s2']

    quietly gen double _ipwterm = in_crisis_`s2'*dy_0/_ps2 - (1-in_crisis_`s2')*dy_0/(1-_ps2)
    quietly gen double _augterm = -((in_crisis_`s2'-_ps2)/(_ps2*(1-_ps2))) * ///
                                   ((1-_ps2)*_m1 + _ps2*_m0)
    quietly summarize _ipwterm if !missing(_ipwterm,_augterm), meanonly
    local ipwm = r(mean)
    quietly summarize _augterm if !missing(_ipwterm,_augterm), meanonly
    local augm = r(mean)
    local tot = `ipwm' + `augm'
    local shr = cond(abs(`tot')>1e-12, 100*abs(`augm')/abs(`tot'), .)
    di as result "  `s2':  IPW piece = " %8.3f `ipwm' ///
                 "   augmentation = " %8.3f `augm' ///
                 "   total = " %8.3f `tot' ///
                 "   |aug|/|total| = " %5.1f `shr' " pct"
}
capture drop _ps2 _m0 _m1 _xb2 _ipwterm _augterm
di as result _n "  A large augmentation share means the propensity and outcome models"
di as result "  disagree about the counterfactual, so the estimate rests on Eq. (1)"
di as result "  being correctly specified rather than on the weighting."

* ══════════════════════════════════════════════════════════════════════════
* 3. Eqs. (1)-(3) BY RESOLUTION TYPE, AND THE DIFFERENCE
*
* Each type is scored against TRANQUIL years with the rival type dropped —
* the reference paper's sample_for* design, so the control group is clean
* tranquil country-years and never the other crisis type.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "3. FLOW AIPW BY RESOLUTION TYPE (Year 1 = the crisis year)"
di as result "════════════════════════════════════════════════════════════"
di as result "Year   ND (se_a)          DEF (se_a)         def-nd   [95% boot CI]      Clogg z    p     draws"
di as result "       se_a = analytic influence-function SE (the paper's own construction)."
di as result "       ND/DEF stars are the paper's 'conventional t-test' vs zero (b/se_a): * p<.10 ** p<.05 *** p<.01."
di as result "       def-nd's own * marks the bootstrap CI excluding 0 -- the conservative, governing test for the difference."

forvalues h = 0/4 {
    local row = `h' + 1
    local hd  = `h' + 1

    _aipwpairflow, y(dy_`h') ///
        d1(in_crisis_def) if1(sample_flow==1 & in_crisis_nd==0) ///
        d2(in_crisis_nd)  if2(sample_flow==1 & in_crisis_def==0) ///
        omod(`com') pz(`cx_active' `cz_recency') reps(`nboot') boot(row)

    if !r(ok) {
        di as error "  Year `hd': estimate failed (cell too thin)."
        continue
    }

    * Pull every returned scalar into a local FIRST. r() survives matrix
    * assignment, but not every command, and reading it late is how the
    * earlier files in this project acquired silent bugs.
    local B1 = r(b1)
    local B2 = r(b2)
    local A1 = r(a1)      // analytic influence-function SE, def  (the paper's)
    local A2 = r(a2)      // analytic influence-function SE, nd
    local DH = r(dh)
    local SE = r(se)
    local LO = r(lo)
    local HI = r(hi)
    local ND = r(nd)

    matrix Fdef_b[`row',1]  = `B1'
    matrix Fnd_b[`row',1]   = `B2'
    matrix Fdef_se[`row',1] = `A1'
    matrix Fnd_se[`row',1]  = `A2'
    matrix Fdiff_b[`row',1]  = `DH'
    matrix Fdiff_se[`row',1] = `SE'
    matrix Fdiff_lo[`row',1] = `LO'
    matrix Fdiff_hi[`row',1] = `HI'

    * Level CIs = theta +/- 1.96 * ANALYTIC SE, which is how the reference paper
    * bands its Fig. 4 (`up = irf + 1.96*se'). The bootstrap is reserved for the
    * DIFFERENCE, exactly as they reserve it: their bootstrap draws are used only
    * to form the between-type contrasts and take centile(2.5 97.5).
    matrix Fnd_lo[`row',1]  = `B2' - 1.96*`A2'
    matrix Fnd_hi[`row',1]  = `B2' + 1.96*`A2'
    matrix Fdef_lo[`row',1] = `B1' - 1.96*`A1'
    matrix Fdef_hi[`row',1] = `B1' + 1.96*`A1'

    * Clogg et al. (1995) z, built from the ANALYTIC SEs — their construction is
    *     clogg = (irf1 - irf2)/(se1^2 + se2^2)^0.5
    * with se the same influence-function SE used for the level bands. It is the
    * PERMISSIVE statistic: it treats the two cells as independent, and they are
    * not, since they share the tranquil control pool. The bootstrap percentile
    * CI is the conservative number and governs where the two disagree.
    local zz = .
    local pz = .
    if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
        local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
        local pz = 2*(1 - normal(abs(`zz')))
        matrix Fdiff_z[`row',1] = `zz'
        matrix Fdiff_p[`row',1] = `pz'
    }

    * Conventional t-test for each level vs zero, exactly as the paper's Table 2
    * reports it: coefficient over its OWN analytic SE, no bootstrap involved --
    * "conventional" there is explicitly contrasted against the bootstrap/Clogg
    * methods used for the DIFFERENCE, not applied to this test.
    local tnd  = cond(`A2'>0, `B2'/`A2', .)
    local pnd  = cond(!missing(`tnd'), 2*(1-normal(abs(`tnd'))), .)
    local sgnd = cond(missing(`pnd'), "", cond(`pnd'<.01,"***",cond(`pnd'<.05,"**",cond(`pnd'<.10,"*",""))))
    local tdef  = cond(`A1'>0, `B1'/`A1', .)
    local pdef  = cond(!missing(`tdef'), 2*(1-normal(abs(`tdef'))), .)
    local sgdef = cond(missing(`pdef'), "", cond(`pdef'<.01,"***",cond(`pdef'<.05,"**",cond(`pdef'<.10,"*",""))))

    local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
    di as result "  " %1.0f `hd' "  " %8.3f `B2' "`sgnd'" " (" %5.3f `A2' ")  " ///
                 %8.3f `B1' "`sgdef'" " (" %5.3f `A1' ")  " %8.3f `DH' ///
                 " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
                 " " %7.3f `zz' " " %5.3f `pz' " " %4.0f `ND'
}

* ══════════════════════════════════════════════════════════════════════════
* 3b. THE SAME DIFFERENCE, COUNTRY-CLUSTER BOOTSTRAP — COMPARISON ONLY
*
* Section 3 uses the reference paper's resampling unit: ROWS, within
* treatment-type strata. Under their onset coding a treated row is close to an
* independent episode — 194 restructurings across 76 countries — so that is a
* reasonable unit. Under flow coding a treated row is a crisis-YEAR: 234 of
* them carry only 61 episodes in 52 countries, and Venezuela alone contributes
* 29. Row resampling draws those 29 as 29 independent observations.
*
* This block re-runs the identical estimator resampling whole COUNTRIES, so the
* cost of that choice is a number rather than an argument. The POINT ESTIMATES
* are identical by construction — a bootstrap does not touch them — so the two
* columns differ only in the interval, which is the point.
*
* Expect the row CIs to be narrower. If they are not, _pool is mis-built (it
* must be row-level, not egen max()-ed to the country).
*
* Default-linked evidence comes from 14 countries: 4 with default episodes only
* (Belize, Brazil, Cote d'Ivoire, Zambia) and 10 with both. Row resampling makes
* that concentration invisible in the interval; cluster resampling does not.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "3b. ROW vs COUNTRY-CLUSTER BOOTSTRAP (same point estimates)"
di as result "════════════════════════════════════════════════════════════"
di as result "Year   def-nd    [row CI: theirs]        [country-cluster CI]      width ratio  draws"

forvalues h = 0/4 {
    local row = `h' + 1
    local hd  = `h' + 1

    _aipwpairflow, y(dy_`h') ///
        d1(in_crisis_def) if1(sample_flow==1 & in_crisis_nd==0) ///
        d2(in_crisis_nd)  if2(sample_flow==1 & in_crisis_def==0) ///
        omod(`com') pz(`cx_active' `cz_recency') reps(`nboot') boot(cluster)

    if !r(ok) {
        di as error "  Year `hd': cluster-bootstrap comparison failed."
        continue
    }
    local CD = r(dh)
    local CL = r(lo)
    local CH = r(hi)
    local CN = r(nd)
    matrix Cdiff_b[`row',1]  = `CD'
    matrix Cdiff_lo[`row',1] = `CL'
    matrix Cdiff_hi[`row',1] = `CH'
    matrix Cdiff_n[`row',1]  = `CN'

    * Point estimates MUST match Section 3 — the check that the resampling has
    * not leaked into the estimate itself.
    local gapchk = abs(`CD' - Fdiff_b[`row',1])
    if `gapchk' > 1e-8 {
        di as error "  ** Year `hd': point estimate differs across bootstrap modes by " ///
                    %9.6f `gapchk' " — resampling has leaked into the estimate."
    }

    local wrow = Fdiff_hi[`row',1] - Fdiff_lo[`row',1]
    local wclu = `CH' - `CL'
    local wrat = cond(`wclu' > 0, `wrow'/`wclu', .)
    local sigr = cond(!missing(Fdiff_lo[`row',1]) & (Fdiff_lo[`row',1]>0 | Fdiff_hi[`row',1]<0), "*", " ")
    local sigc = cond(!missing(`CL') & (`CL'>0 | `CH'<0), "*", " ")

    di as result "  " %1.0f `hd' "   " %8.3f `CD' ///
                 "  [" %7.3f Fdiff_lo[`row',1] ", " %7.3f Fdiff_hi[`row',1] "]`sigr'" ///
                 "  [" %7.3f `CL' ", " %7.3f `CH' "]`sigc'" ///
                 "   " %5.2f `wrat' "   " %4.0f `CN'
}

di as result _n "  A width ratio well below 1 means the row bootstrap is treating"
di as result "  repeated crisis-years of the same country as independent draws."
di as result "  Where the two verdicts (*) differ, the write-up reports both and"
di as result "  says which resampling unit produced which."

di as result _n "  * = bootstrap 95% percentile CI for the gap excludes zero."
di as result "  Clogg z assumes the two cells are independent (they share the"
di as result "  tranquil control pool), so it is the permissive statistic. Where"
di as result "  it and the bootstrap disagree, the bootstrap governs."

* ══════════════════════════════════════════════════════════════════════════
* 4. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
preserve
    clear
    tempname pfd
    tempfile diff2
    postfile `pfd' int horizon double dhl double se double lo double hi ///
        double cloggz double cloggp double clu_lo double clu_hi ///
        using "`diff2'", replace
    post `pfd' (0) (0) (0) (0) (0) (.) (.) (.) (.)
    forvalues h = 1/5 {
        post `pfd' (`h') (Fdiff_b[`h',1]) (Fdiff_se[`h',1]) ///
            (Fdiff_lo[`h',1]) (Fdiff_hi[`h',1]) (Fdiff_z[`h',1]) (Fdiff_p[`h',1]) ///
            (Cdiff_lo[`h',1]) (Cdiff_hi[`h',1])
    }
    postclose `pfd'
    use "`diff2'", clear
    label var dhl     "Flow AIPW extra cost of default (def - nd, pp)"
    label var lo      "95% percentile CI lower (bootstrap)"
    label var hi      "95% percentile CI upper (bootstrap)"
    label var cloggz  "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp  "p-value of the Clogg z"
    label var clu_lo  "95% CI lower, COUNTRY-CLUSTER bootstrap (comparison)"
    label var clu_hi  "95% CI upper, COUNTRY-CLUSTER bootstrap (comparison)"
    gen byte sig95 = (!missing(lo) & (lo>0 | hi<0))
    label var sig95 "Row bootstrap CI excludes 0 (baseline, the paper's scheme)"
    gen byte sig95_clu = (!missing(clu_lo) & (clu_lo>0 | clu_hi<0))
    label var sig95_clu "Country-cluster CI excludes 0 (comparison)"
    order horizon dhl se lo hi sig95 clu_lo clu_hi sig95_clu cloggz cloggp
    export delimited "$tabs/aipw_flow_diff.csv", replace
    di as result "Flow AIPW difference saved: $tabs/aipw_flow_diff.csv"
restore

* ── IRF datasets ────────────────────────────────────────────────────────────
foreach g in nd def {
    preserve
        clear
        set obs 6
        gen horizon = _n - 1
        gen double b  = 0
        gen double lo = 0
        gen double hi = 0
        forvalues h = 1/5 {
            quietly replace b  = F`g'_b[`h',1]  in `=`h'+1'
            quietly replace lo = F`g'_lo[`h',1] in `=`h'+1'
            quietly replace hi = F`g'_hi[`h',1] in `=`h'+1'
        }
        gen series = "aipw_flow_`g'"
        save "$clean/irf_aipw_flow_`g'.dta", replace
    restore
}
di as result "IRF datasets saved: irf_aipw_flow_nd.dta, irf_aipw_flow_def.dta"

* ══════════════════════════════════════════════════════════════════════════
* 5. FIGURE — drawn here, as 20_lp_flow.do does, so the file is self-contained
* ══════════════════════════════════════════════════════════════════════════
local c_nd   "0 84 166"
local c_def  "157 36 73"
local c_zero "150 150 150"

* Common y-scale for the Fig.3 / Fig.4 pair — set this to the SAME string as
* the `yrng' local in 20_lp_flow.do once both have run, so the OLS and AIPW
* panels are literally comparable rather than merely similar.
local yrng ""

preserve
    use "$clean/irf_aipw_flow_nd.dta", clear
    append using "$clean/irf_aipw_flow_def.dta"
    twoway ///
        (rarea lo hi horizon if series=="aipw_flow_nd",  color("`c_nd'%20")  lwidth(none)) ///
        (rarea lo hi horizon if series=="aipw_flow_def", color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if series=="aipw_flow_nd",  lcolor("`c_nd'")  mcolor("`c_nd'") ///
            msymbol(circle) lwidth(medthick)) ///
        (connected b horizon if series=="aipw_flow_def", lcolor("`c_def'") mcolor("`c_def'") ///
            msymbol(square) lpattern(dash) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
        xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year (Year 1 = the crisis year)", size(medsmall)) ///
        ytitle("Change in log real GDP (pp)", size(medsmall)) ///
        title("Output While In a Spread Crisis, by Resolution", size(medium)) ///
        subtitle("AIPW — the reference paper's Eq. (3). Tranquil country-years omitted.", size(small)) ///
        `yrng' ///
        note("95% CIs: theta +/- 1.96 x the analytic influence-function SE, the construction" ///
             "behind the reference paper's Fig. 4. The bootstrap is reserved for the def-nd" ///
             "difference, as it is in their design. Note their SE is unclustered; the per-type" ///
             "bootstrap SD is printed in the log for comparison." ///
             "The AIPW half of a matched pair: this is the analogue of their" ///
             "Fig. 4, and fig9b_irf_flow_resolution.pdf (20_lp_flow.do) is their Fig. 3 — the same" ///
             "object estimated by OLS. Same treatment, same axis, same scale, so the two read side" ///
             "by side and the estimator is the only thing that differs." ///
             "Eq. (1) plain OLS with country FE and episode-dated controls; Eq. (2) pooled probit" ///
             "on row-dated controls plus excluded predictors; Eq. (3) the paper's exact AIPW form.", ///
             size(vsmall)) ///
        legend(order(3 "Non-default" 4 "Default-linked") ring(0) pos(7) cols(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    capture graph export "$figs/fig10_aipw_flow.pdf", replace
    if _rc di as error "  ** fig10_aipw_flow.pdf export failed (rc=" _rc ") — is it open?"
    else {
        capture graph export "$figs/fig10_aipw_flow.png", replace width(1200)
        di as result "Figure 10 saved: fig10_aipw_flow.pdf/.png"
    }
restore

capture drop _ps_nd _ps_def

di as result _n "21_aipw_flow.do complete."
