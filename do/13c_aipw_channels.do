/*===========================================================================
  13C_AIPW_CHANNELS.DO
  Doubly-robust AIPW (Asonuma et al. Eq. 3 / Jordà–Taylor 2016) for the
  TRANSMISSION CHANNELS — the same estimator as 08b_aipw.do (GDP per capita),
  applied to each channel outcome so the whole transmission story is estimated
  doubly-robustly and coherently with the headline GDP result.

  Per channel we report BOTH:
    Act 1 — all-crises AIPW (onset_all vs tranquil), ATE.
    Act 2 — resolution split as TWO level lines (non-default vs tranquil and
            default-linked vs tranquil, rival onset dropped), mirroring 08b's
            fig_aipw_act2.

  Channels (level-difference outcomes, ch_v_h = F h.v - L.v, h=0..4), reusing the
  construction and control sets from 11/11b/12/13. ACTIVE (estimated below):
    credit, inv                                         (11 / 12)
    claims_govt                                         (12)
    claimsgov_assets, claimpriv_assets                  (11b nexus)
  SILENCED for now (not important currently -- see the estimation loop's own
  note; outcome construction still runs, only estimation is skipped):
    govexp, pb, fdi                                     (11 / 12)
    ca (current account)                                (13 Aguiar-Gopinath)

  PROPENSITY model = identical to 08b (selection into treatment is the same
  object regardless of outcome). Only the OUTCOME regression is channel-specific.

  DIAGNOSTIC (console output, before the estimation loop): lists the most
  extreme def-arm channel outcome rows per active channel/horizon. Added
  after a full run showed EVERY level t-test at p<.01 across every channel
  and horizon -- implausible on ~20 default-linked episodes, and consistent
  with a handful of extreme rows dominating the outcome regression's
  treatment coefficient (which inflates the point estimate without
  proportionally inflating the analytic SE). Diagnostic only, does not
  change the estimation.

  INFERENCE, ALIGNED WITH 08b_aipw.do / THE FLOW TIER'S PRESENTATION
  (21_aipw_flow.do): Act 1 (single ATE per channel) keeps its cluster
  bootstrap CI as before -- no by-type contrast to test. Act 2 (by
  resolution type) and its difference are now estimated in ONE pass per
  channel/horizon (_aipwpair, replacing the old separate _aipwci-levels
  loop plus a second _aipwpairdiff-levels-recomputed-again loop): each
  level's CI is 1.96*analytic influence-function SE (the paper's own
  Table 2/Fig. 4 construction), with a conventional t-test (b/se) and
  stars per level; the def-nd DIFFERENCE is bootstrapped directly with
  ROW-LEVEL resampling within control/nd/def pools -- the paper's own
  bootstrap device, a natural fit here since an onset row already is one
  episode; Clogg et al. (1995)'s z (from the same analytic SEs) is
  reported alongside the bootstrap CI as the permissive companion
  statistic. See 08b_aipw.do's header for the full argument.

  CRITICAL: bsample (cluster bootstrap) destroys the panel time order, so nothing
  re-estimated inside the bootstrap may use L./F. operators. The outcome ch_*_h is
  precomputed, and every lagged control is pre-generated as a PLAIN column
  (l_credit = L.credit, ...) before estimation. cx/cz are already plain columns.

  Output: $tabs/aipw_channels.csv ; $figs/fig_aipw_ch_act2.pdf (Act 1 and its
          fig_aipw_ch_act1.pdf are SILENCED, see below). Leaves 11/11b/12/13
          (OLS+IPW) untouched.
  Runtime note: heavy (5 active channels x ~15 fits x nboot; 4 more channels'
  outcome construction runs but their estimation is silenced, see above).
  nboot=300 for a practical
  run; raise to 500 for the final.  Run AFTER 17_predictors.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ── REPRODUCIBILITY: seed the bootstrap ────────────────────────────────────
* Every CI in this file comes from `bsample', which draws at random. Without a
* seed the intervals move between runs of identical code, and on thin cells
* that is not a rounding issue -- in 13d the non-default nexus gap at Year 1
* returned [0.05, 7.37] in one run and [-0.088, 7.686] in the next, opposite
* verdicts at the 5% line from the same data. Seeding makes the reported
* interval a property of the estimator rather than of the session. The value
* is arbitrary and was not chosen by inspecting results.
set seed 20260819

local nboot  = 300      // bootstrap reps per (channel, series, horizon)
local cx     $ctrl_core   // retained for reference; propensity baseline now passes `om' (strict parity)
local cz     l_fedfunds l_reg_crisis_share past_onsets       // Act 1 predictors Z1
* Act 2 predictors Z2 (resolution-type): both terms are DEFAULT-LINKED-
* SPECIFIC (l_contagion_dist_def, years_since_def_onset) on economic
* grounds -- a predictor for default-linked risk should measure
* default-linked distress/recency, not spread-crisis distress in general.
* Corroborated by 08c_first_stage_table.do's diagnostics (the default arm's
* classification power was materially weaker under the generic
* combination); see 08c's header for the full argument and the caveat that
* the reference paper's own instrument is not tailored per column.
local cz_def l_fedfunds l_contagion_dist_def years_since_def_onset

* AIPW outcome-model core = the common core, unchanged. The depth term is
* l_credit_bank (WDI FD.AST.PRVT.GD.ZS, credit by banks).
* This block used to swap it for l_credit (FS.AST.PRVT.GD.ZS, all financial
* corporations) on the belief that the by-banks series had the coverage hole.
* 19_sample_audit.do showed the opposite: l_credit_bank is non-missing on EVERY
* row where l_credit is (overlap 1180 = all of l_credit) plus 141 more, and the
* two correlate 0.950. The swap was costing observations, not saving them, so
* the core now carries l_credit_bank everywhere and this file needs no exception.
* Both are plain saved columns, so the bootstrap stays operator-free either way.
local core_aipw l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2

* ── Build channel outcomes ch_v_h = F h.v - L.v (h=0..4) ─────────────────────
* OUTCOME SCALE. Strictly-positive GDP-ratio channels use the LOG REAL LEVEL
* (ln_r_*, built in 18_transforms), matching the reference paper's var2/var3:
* a change in X/GDP confounds X with a GDP that is collapsing, whereas
* ln(X/GDP * GDP) = ln(X) up to a constant, so the outcome is the cumulative
* percent change in X itself. pb, fdi and ca change sign so they keep the ratio;
* claimsgov_assets and claimpriv_assets are shares of BANK ASSETS, not of GDP, so
* the denominator problem does not arise for them either.
foreach v in credit claims_govt inv govexp pb fdi ///
             claimsgov_assets claimpriv_assets ca {
    local src `v'
    if inlist("`v'","credit","inv","govexp") local src ln_r_`v'
    capture drop `v'_base
    gen `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var);
    * added to each outcome model to absorb the channel's own pre-trend momentum.
    capture drop pre_`v'
    gen pre_`v' = L.`src' - L2.`src'
}

* ── Pre-lag every lagged control to a PLAIN column (bsample-safe) ────────────
foreach v in credit claims_govt pb govexp fdi ca ///
             claimsgov_assets claimpriv_assets {
    capture drop l_`v'
    gen l_`v' = L.`v'
}

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — hand-coded AIPW point estimate (Eqs. 1-3); r(theta), r(N)
*   (identical to 08b_aipw.do)
* ══════════════════════════════════════════════════════════════════════════
* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — bootstrap stratum (see 08b_aipw.do for the full rationale).
* The paper resamples the control pool and each treated type separately, holding
* the treated count fixed. A plain cluster bsample does not, and on cells this
* thin a draw can lose the treated group entirely. This makes the stratum so
* `bsample, cluster(cid) strata()' preserves the ever-treated country count.
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
        quietly gen byte `t2' = `d2' if `touse'
        tempvar s2
        bysort cid: egen byte `s2' = max(`t2')
        quietly replace `s2' = 0 if missing(`s2')
        quietly replace `generate' = `generate' + 2*`s2'
    }
end

capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'
    tempvar xb m0 m1 ps summ iwt
    * IPWRA, matching the reference paper: the propensity is estimated FIRST and
    * the outcome regression that produces mu0/mu1 is IPW-WEIGHTED
    * (their `reg g_h dum g_0 $convar [pweight=invwt]`). Estimating it unweighted
    * gives a different, also doubly-robust, estimator - but not theirs.
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' [pweight=`iwt'] if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' [pweight=`iwt'] if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'   // set D=0
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'   // set D=1
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    local th = r(mean)
    local nn = r(N)

    * Analytic (unclustered) influence-function SE -- matches 08b_aipw.do /
    * 21_aipw_flow.do's _aipw exactly: sqrt(mean((summand - theta)^2) / N).
    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    * DIAGNOSTIC-ONLY: country-clustered version of the same SE -- see
    * 08b_aipw.do's _aipw for the full rationale. Not used for the adopted
    * level bands, returned only for comparison.
    quietly regress `summ' if `touse', vce(cluster cid)
    local seclu = _se[_cons]

    return scalar theta  = `th'
    return scalar N      = `nn'
    return scalar se     = `sean'
    return scalar se_clu = `seclu'
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — point estimate + cluster-bootstrap percentile CI for one cell
*   _aipwci <yvar> <Dvar>, ifc(<cond>) omod(<X>) pz(<X Z>) reps(<n>)
*   returns r(ok) r(b) r(se) r(lo) r(hi) r(nd)
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipwci
program define _aipwci, rclass
    syntax anything, IFC(string) OMOD(string) PZ(string) REPS(integer)
    gettoken yv Dv : anything

    capture _aipw `yv' `Dv' if `ifc', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local pt = r(theta)

    * stratum: ever-treated vs never-treated countries within this cell, so every
    * draw keeps the same number of treated countries (see _mkstrat).
    _mkstrat `Dv' if `ifc', generate(_strat)
    tempname pf
    tempfile bf
    quietly postfile `pf' double theta using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            capture drop _bid
            bsample, cluster(cid) strata(_strat) idcluster(_bid)
            capture _aipw `yv' `Dv' if `ifc', omodel(`omod') pmodel(`pz') fe(_bid)
            if _rc == 0 quietly post `pf' (r(theta))
        restore
    }
    quietly postclose `pf'
    capture drop _strat

    local se = .
    local lo = .
    local hi = .
    local nd = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(theta)
        local nd = r(N)
        if `nd' >= 50 {
            quietly summarize theta
            local se = r(sd)
            _pctile theta, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
        }
    restore

    return scalar ok = 1
    return scalar b  = `pt'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — paired bootstrap of the DIFFERENCE between two AIPW cells,
*   ROW-LEVEL resampling within control/type1/type2 pools -- the reference
*   paper's own bootstrap device (see 08b_aipw.do's _aipwpair header for the
*   full rationale). Here cell1 = default vs tranquil, cell2 = non-default
*   vs tranquil, so r(dh) = def - nd = the extra channel response of
*   default. Also returns each cell's own analytic SE (a1/a2) so the console
*   table and Clogg z are built from the same estimates as the bootstrap.
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
    local c1 = r(se_clu)   // diagnostic-only country-clustered SE
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2 = r(theta)
    local a2 = r(se)
    local c2 = r(se_clu)
    local dh = `b1' - `b2'

    capture drop _pool
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
    return scalar c1 = `c1'
    return scalar c2 = `c2'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* ESTIMATE — loop over channels; Act 1 + Act 2 (nd/def); post to results file
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str24 channel str4 series byte horizon double b se lo hi ///
    using "`resf'", replace

* second results file: the def - nd DIFFERENCE per channel x horizon (row boot)
tempname Rd
tempfile diffresf
postfile `Rd' str24 channel byte horizon double dhl bdef bnd se lo hi nd ///
    double cloggz double cloggp using "`diffresf'", replace

* ══════════════════════════════════════════════════════════════════════════
* DIAGNOSTIC: are the def-arm outcome values themselves extreme?
*
* Every level t-test across every channel/horizon came back p<.01 in the
* last full run -- implausible on ~20 default-linked episodes, and a
* classic symptom of the outcome regression's treatment coefficient being
* dominated by a handful of extreme rows (which inflates the point
* estimate without proportionally inflating the analytic SE, since m0/m1
* are a single shared shift, not independently varying per row -- see the
* conversation this diagnostic responds to). For each active channel, list
* the most extreme onset_def==1 rows at h=1 and h=4 (h=0/3 in code), so a
* dominant outlier is visible directly rather than inferred from the
* aggregate coefficient.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== DIAGNOSTIC: extreme def-arm channel outcomes (candidates for outlier-driven ATEs) ==="
foreach ch in credit inv claims_govt claimsgov_assets claimpriv_assets {
    foreach h in 0 3 {
        capture confirm variable ch_`ch'_`h'
        if !_rc {
            di as result "  `ch', h=" `h'+1 " (onset_def==1 rows, sorted by |value|):"
            preserve
                quietly keep if sample==1 & onset_def==1 & !missing(ch_`ch'_`h')
                quietly count
                if r(N) > 0 {
                    gen double _absval = abs(ch_`ch'_`h')
                    gsort -_absval
                    local ntop = min(5, r(N))
                    list country year ch_`ch'_`h' in 1/`ntop', noobs clean
                }
                else di as result "    (no non-missing rows)"
            restore
        }
    }
}

* govexp, pb, ca, fdi SILENCED below -- not important for now (no reliable
* channel signal so far: govexp/pb/ca/fdi's bootstrap CI never excluded
* zero at any horizon in the last full run). Left out of the active list,
* not deleted -- restore by adding them back: "credit claims_govt inv
* govexp pb fdi claimsgov_assets claimpriv_assets ca". Their outcome
* construction (ch_v_h/pre_v/l_v) above still runs regardless, since it is
* cheap and shared -- only the estimation loop below is skipped for them.
foreach ch in credit claims_govt inv ///
              claimsgov_assets claimpriv_assets {

    * channel-specific OUTCOME-model controls (pre-lagged plain columns)
    * AIPW outcome core ($core_aipw = the common core, depth term l_credit_bank) +
    * channel's own pre_<v>; drop the core term equal to the channel's own lagged
    * level (credit->the depth term, govexp->l_govexp). "ca" needed no special
    * case even before the core control set was unified project-wide (l_ca was
    * a core term then; it no longer is, replaced by exchange2, which is a
    * different variable than the current account and carries no tautology
    * risk for the ca outcome) -- it now falls through to the `else' branch.
    if      "`ch'" == "credit"            local om l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
    else if "`ch'" == "govexp"            local om l1_gdpg l_debt l_banking_crisis l_open l_credit_bank l_lninfl exchange2 pre_govexp
    else                                  local om `core_aipw' pre_`ch'

    di as result _n "=== CHANNEL: `ch' ==="

    * ── Act 1: SILENCED. Only Act 2 (by resolution) is of interest now. Left
    * commented, not deleted -- restoring it means adding back "all A1" to
    * the export foreach below and un-silencing Figure A. _aipwci (the
    * bootstrap-CI helper it used) is left defined, just uncalled.
    /*
    di as result "  Act 1 (all onsets):   h    ATE      SE      [95% CI]   draws"
    post `R' ("`ch'") ("all") (0) (0) (0) (0) (0)   // explicit baseline (h=0), matching Asonuma et al.
    forvalues h = 0/4 {
        * Strict parity: propensity baseline = the OUTCOME model `om' (= $ctrl_core
        * + own pre-trend), plus Z — their g_0+$convar in both stages, $instrument added.
        _aipwci ch_`ch'_`h' onset_all, ifc(sample==1) ///
            omod(`om') pz(`om' `cz') reps(`nboot')
        if r(ok) {
            local b=r(b)
            local se=r(se)
            local lo=r(lo)
            local hi=r(hi)
            local nd=r(nd)
            post `R' ("`ch'") ("all") (`h'+1) (`b') (`se') (`lo') (`hi')
            di "    h=" `h'+1 "  " %8.3f `b' "  " %6.3f `se' ///
               "  [" %7.3f `lo' ", " %7.3f `hi' "]  " `nd' "/`nboot'"
        }
        else di as error "    h=" `h'+1 ": Act 1 estimate failed (rc)."
    }
    */

    * ── Act 2: resolution split (levels, analytic-SE bands + t-test) AND the
    * def-nd difference (row bootstrap + Clogg z), estimated together in ONE
    * pass per horizon via _aipwpair -- see 08b_aipw.do's header for why this
    * replaces two separate estimation loops (levels via bootstrap CI, then
    * the diff re-estimating both cells again) with a single, more efficient
    * and internally consistent one.
    di as result "  Act 2:  h   ND (se_a)        DEF (se_a)        def-nd   [95% boot CI]   Clogg z    p"
    di as result "           ND/DEF stars are the conventional t-test vs zero (b/se_a): * p<.10 ** p<.05 *** p<.01."
    post `R' ("`ch'") ("nd")  (0) (0) (0) (0) (0)   // explicit baseline (h=0)
    post `R' ("`ch'") ("def") (0) (0) (0) (0) (0)
    post `Rd' ("`ch'") (0) (0) (0) (0) (0) (0) (0) (0) (.) (.)   // explicit baseline (h=0)
    forvalues h = 0/4 {
        _aipwpair, y(ch_`ch'_`h') ///
            d1(onset_def) if1(sample==1 & onset_nd==0) ///
            d2(onset_nd)  if2(sample==1 & onset_def==0) ///
            omod(`om') pz(`om' `cz_def') reps(`nboot')
        if r(ok) {
            local B1 = r(b1)   // default-linked ATE
            local B2 = r(b2)   // non-default ATE
            local A1 = r(a1)   // analytic SE, default-linked (unclustered, adopted)
            local A2 = r(a2)   // analytic SE, non-default (unclustered, adopted)
            local C1 = r(c1)   // DIAGNOSTIC-ONLY: country-clustered SE, default-linked
            local C2 = r(c2)   // DIAGNOSTIC-ONLY: country-clustered SE, non-default
            local DH = r(dh)
            local SE = r(se)
            local LO = r(lo)
            local HI = r(hi)
            local ND = r(nd)

            * Level CIs = theta +/- 1.96*analytic SE, matching the paper's
            * own Fig. 4 band construction; the bootstrap is reserved for the
            * difference only.
            post `R' ("`ch'") ("nd")  (`h'+1) (`B2') (`A2') (`B2'-1.96*`A2') (`B2'+1.96*`A2')
            post `R' ("`ch'") ("def") (`h'+1) (`B1') (`A1') (`B1'-1.96*`A1') (`B1'+1.96*`A1')

            local zz = .
            local pz = .
            if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
                local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
                local pz = 2*(1 - normal(abs(`zz')))
            }
            post `Rd' ("`ch'") (`h'+1) (`DH') (`B1') (`B2') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')

            * Conventional t-test for each level vs zero (b / own analytic SE),
            * as the paper's Table 2/Fig. 4 -- same construction as 08b_aipw.do.
            local tnd  = cond(`A2'>0, `B2'/`A2', .)
            local pnd  = cond(!missing(`tnd'), 2*(1-normal(abs(`tnd'))), .)
            local sgnd = cond(missing(`pnd'), "", cond(`pnd'<.01,"***",cond(`pnd'<.05,"**",cond(`pnd'<.10,"*",""))))
            local tdef  = cond(`A1'>0, `B1'/`A1', .)
            local pdef  = cond(!missing(`tdef'), 2*(1-normal(abs(`tdef'))), .)
            local sgdef = cond(missing(`pdef'), "", cond(`pdef'<.01,"***",cond(`pdef'<.05,"**",cond(`pdef'<.10,"*",""))))

            local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
            di "    " %1.0f `h'+1 "  " %8.3f `B2' "`sgnd'" " (" %5.3f `A2' ")  " ///
               %8.3f `B1' "`sgdef'" " (" %5.3f `A1' ")  " %8.3f `DH' ///
               " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
               " " %7.3f `zz' " " %5.3f `pz'

            * DIAGNOSTIC-ONLY: clustered-SE comparison, not adopted -- see
            * _aipw's header (08b_aipw.do) for the full rationale.
            local widnd  = cond(`A2'>0, `C2'/`A2', .)
            local widdef = cond(`A1'>0, `C1'/`A1', .)
            di "         [clustered SE diag: ND se_clu=" %5.3f `C2' " (x" %4.2f `widnd' ")" ///
               "   DEF se_clu=" %5.3f `C1' " (x" %4.2f `widdef' ")]"
        }
        else di as error "    h=" `h'+1 ": Act 2 estimate failed (thin sample)."
    }
}
postclose `R'
postclose `Rd'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — CSV + two small-multiple (by-channel) figures
* ══════════════════════════════════════════════════════════════════════════
* ── Difference CSV: def - nd gap per channel x horizon (row bootstrap) ──────
preserve
    use "`diffresf'", clear
    label var dhl  "AIPW def - nd channel gap (pp)"
    label var bdef "Default-linked ATE"
    label var bnd  "Non-default ATE"
    label var lo   "95% percentile CI lower (row bootstrap)"
    label var hi   "95% percentile CI upper (row bootstrap)"
    label var nd   "Valid bootstrap draws"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
    label var sig95 "Bootstrap CI excludes 0 (governing test)"
    order channel horizon dhl bdef bnd se lo hi nd sig95 cloggz cloggp
    export delimited "$tabs/aipw_channels_diff.csv", replace
    di as result _n "AIPW channel def-nd difference CSV saved: $tabs/aipw_channels_diff.csv"
restore

use "`resf'", clear
label var b  "AIPW ATE (pp of the channel ratio)"
label var se "Analytic influence-function SE (Act 2 levels; Act 1 silenced, see header)"
label var lo "95% CI lower = theta-1.96*se (analytic)"
label var hi "95% CI upper = theta+1.96*se (analytic)"
order channel series horizon b se lo hi
export delimited "$tabs/aipw_channels.csv", replace
di as result _n "AIPW channel results CSV saved: $tabs/aipw_channels.csv"

encode channel, gen(chid)

* ── Figure A: Act 1 all-crises AIPW per channel — SILENCED, see header. ─────
/*
local c1 "23 55 94"
capture twoway ///
    (rarea lo hi horizon if series=="all", color("`c1'%18") lwidth(none)) ///
    (connected b horizon if series=="all", lcolor("`c1'") lwidth(medthick) msymbol(circle)), ///
    by(chid, yrescale legend(off) ///
        note("Doubly-robust AIPW (Asonuma et al. Eq. 3), ATE. Shaded = bootstrap 95% percentile CI.", size(vsmall)) ///
        title("AIPW transmission channels — all crises", size(medsmall) color(navy))) ///
    yline(0, lpattern(dash) lcolor(gs8)) ///
    xlabel(0(1)5) xtitle("Year (Year 1 = crisis year)", size(small)) ///
    ytitle("Cumulative change in channel (pp)", size(small)) ///
    graphregion(color(white)) plotregion(color(white))
if _rc == 0 {
    graph export "$figs/fig_aipw_ch_act1.pdf", replace
    di as result "Figure saved: fig_aipw_ch_act1.pdf"
}
else di as error "  ** fig_aipw_ch_act1 failed (rc=" _rc ")"
*/

* ── Figure B: Act 2 resolution split per channel (nd vs def, fig8 palette) ──
local c_nd  "0 84 166"
local c_def "157 36 73"
capture twoway ///
    (rarea lo hi horizon if series=="nd",  color("`c_nd'%16")  lwidth(none)) ///
    (rarea lo hi horizon if series=="def", color("`c_def'%16") lwidth(none)) ///
    (connected b horizon if series=="nd",  lcolor("`c_nd'")  lwidth(medthick) msymbol(circle)) ///
    (connected b horizon if series=="def", lcolor("`c_def'") lwidth(medthick) lpattern(dash) msymbol(square)), ///
    by(chid, yrescale ///
        note("Two AIPW level IRFs per channel; shaded = 1.96*analytic SE band. Gap = extra cost of default, bootstrapped directly (row-level).", size(vsmall)) ///
        title("AIPW transmission channels by resolution", size(medsmall) color(navy))) ///
    yline(0, lpattern(dash) lcolor(gs8)) ///
    xlabel(0(1)5) xtitle("Year (Year 1 = crisis year)", size(small)) ///
    ytitle("Cumulative change in channel (pp)", size(small)) ///
    legend(order(3 "Non-default" 4 "Default-linked") size(small)) ///
    graphregion(color(white)) plotregion(color(white))
if _rc == 0 {
    graph export "$figs/fig_aipw_ch_act2.pdf", replace
    di as result "Figure saved: fig_aipw_ch_act2.pdf"
}
else di as error "  ** fig_aipw_ch_act2 failed (rc=" _rc ")"

di as result _n "13c_aipw_channels.do complete."
di as result "Compare the AIPW channel IRFs to the OLS/IPW versions in 11/12 (same"
di as result "sign, default deeper). aipw_channels.csv has one row per channel x"
di as result "series(all/nd/def) x horizon."
di as result "aipw_channels_diff.csv bootstraps the def-nd GAP per channel directly"
di as result "(the paper's between-type inference); sig95==1 flags a gap whose 95% CI"
di as result "excludes 0 at that horizon."
