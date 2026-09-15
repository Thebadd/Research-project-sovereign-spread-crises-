/*===========================================================================
  13C_AIPW_CHANNELS_ANALYTICSE.DO — DUPLICATE, PAPER-ALIGNED SE VARIANT

  Standalone duplicate of 13c_aipw_channels.do, NOT wired into 00_master.do.
  The ONE change: every channel's LEVEL (nd/def) confidence intervals and
  t-tests/stars use the paper's own ANALYTIC (unclustered influence-function)
  SE (`A1'/`A2', already returned by _aipwpair) instead of 13c's ADOPTED
  row-bootstrap SE (`BSE1'/`BSE2'). See 08b_aipw_analyticSE.do's header for
  the full rationale -- identical reasoning, applied here per channel.
  Everything else (point estimates, the def-nd difference bootstrap and its
  CI, Clogg z) is UNCHANGED. Outputs write to DIFFERENT filenames so this
  never overwrites 13c_aipw_channels.do's own results.

  SELF-CONTAINED, mirroring 13c_aipw_channels.do's own addition: GDP is now
  estimated IN THIS FILE too ("=== CHANNEL: gdp ===" block, pulling its
  spec from 08b_aipw_analyticSE.do's own Act 2 GDP loop -- NOT the plain
  08b_aipw.do -- so the level bands here stay analytic-SE throughout, this
  file's one consistent convention). The combined 6-panel figure and the
  merged six-variable table/CSV below
  ($tabs/aipw_combined_six_resolution_analyticSE.rtf,
  $tabs/aipw_combined_six_analyticSE.csv) are therefore produced entirely
  from this file's own estimation, with no import of 08b_aipw_analyticSE.do's
  saved aipw_results_analyticSE.csv. GDP is estimated redundantly across
  THREE files project-wide as a result (08b_aipw.do, 08b_aipw_analyticSE.do,
  and here) -- accepted for the same reason as the plain file: standalone
  runnability over saving duplicate computation. 08b_aipw_analyticSE.do
  itself is left completely unchanged.
===========================================================================*/

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
    real_lending                                        (18_transforms; new WDI channel)
    fdi                                                 (11 / 12; added for the combined
                                                          Panel A-F GDP/Investment/Bank credit/
                                                          Claims on government/FDI/Real
                                                          lending-rate figure below)
  SILENCED for now (not important currently -- see the estimation loop's own
  note; outcome construction still runs, only estimation is skipped):
    govexp, pb                                          (11 / 12)
    ca (current account)                                (13 Aguiar-Gopinath)

  BALANCED A-D SAMPLE (common_abcd, built in 18_transforms.do): credit, inv,
  and claims_govt are now ACTUALLY estimated on the common episode set
  shared with 08b_aipw.do's AIPW GDP result and 03_lp_resolution.do's
  Table 2 (Asonuma et al.'s own "balanced panels" convention -- see the
  per-channel loop below for the full argument). claimsgov_assets,
  claimpriv_assets, real_lending, and fdi keep their own best-available
  sample, unrestricted.

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

  INFERENCE, ALIGNED WITH 08b_aipw.do: Act 1 (single ATE per channel) keeps
  its cluster bootstrap CI as before -- no by-type contrast to test. Act 2
  (by resolution type) and its difference are estimated in ONE pass per
  channel/horizon (_aipwpair): each level's CI is 1.96*ROW-BOOTSTRAP SE
  (ADOPTED), with a conventional t-test (b/se_boot) and stars per level.
  THIS IS A DELIBERATE DEPARTURE FROM THE PAPER'S OWN ANALYTIC-SE
  CONSTRUCTION -- a direct diagnostic showed that formula understated the
  def arm's true uncertainty by 3.75-5.5x on ~20-episode default arms (see
  08b_aipw.do's header for the full argument and evidence -- the
  analytic-vs-bootstrap/clustered comparison lines that established this
  are no longer printed each run). The def-nd DIFFERENCE is
  bootstrapped directly with ROW-LEVEL resampling within control/nd/def
  pools -- the paper's own bootstrap device, a natural fit here since an
  onset row already is one episode; Clogg et al. (1995)'s z keeps the
  analytic SEs -- confirmed to match a line literally in their own replication script, not an outside convention -- and is reported as the
  permissive companion statistic on a different SE basis than the level
  display.

  CRITICAL: bsample (cluster bootstrap) destroys the panel time order, so nothing
  re-estimated inside the bootstrap may use L./F. operators. The outcome ch_*_h is
  precomputed, and every lagged control is pre-generated as a PLAIN column
  (l_credit = L.credit, ...) before estimation. cx/cz are already plain columns.

  Output: $tabs/aipw_channels.csv ; $figs/fig_aipw_ch_act2.pdf (Act 1 and its
          fig_aipw_ch_act1.pdf are SILENCED, see below); $figs/fig_aipw_combined.pdf
          (Panel A-F: GDP [imported from 08b_aipw.do's aipw_results.csv],
          Investment, Bank credit, Claims on government, FDI, Real lending
          rate). Leaves 11/11b/12/13 (OLS+IPW) untouched.
  Runtime note: heavy (7 active channels x ~15 fits x nboot; 2 more channels'
  outcome construction runs but their estimation is silenced, see above).
  nboot=1000, matching the reference paper's own G=1000 (their bootstrap
  scripts) and 13d_aipw_nexus_split.do's own setting. Run AFTER 17_predictors.do and 08b_aipw.do
  (the combined figure reads 08b's saved aipw_results.csv for GDP).
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

local nboot  = 1000     // matches the reference paper's own G=1000 (their bootstrap scripts)
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
* THIRD PREDICTOR CHANGE, ADOPTED (matching 13c_aipw_channels.do):
* l_contagion_dist_def -> l_contagion_dist_atdef (AT-database-wide donor
* pool; see 17_predictors.do and 08c_first_stage_table.do's comparison).
local cz_def l_fedfunds l_contagion_dist_atdef years_since_def_onset

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
             claimsgov_assets claimpriv_assets ca real_lending {
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

    * Country count in this cell's regression sample -- ported from 13d_aipw_
    * nexus_split.do's own _aipw (same tag(cid)-under-e(sample) mechanism, not
    * invented fresh) so the merged table can print a full Observations/
    * Countries/Episodes triple per level row, matching the reference paper's
    * own Table 2 layout, instead of the coefficient/SE alone.
    tempvar _tagcty
    quietly egen byte `_tagcty' = tag(cid) if `touse'
    quietly count if `_tagcty'==1
    local nctry = r(N)

    * Treated-episode count in this arm: one treated row is one episode under
    * this project's onset-tier design, so it is simply the treated-dummy
    * count within e(sample) -- no separate episode-counting logic needed.
    quietly count if `touse' & `D'==1
    local ntr = r(N)

    * Analytic (unclustered) influence-function SE -- matches 08b_aipw.do /
    * 21_aipw_flow.do's _aipw exactly: sqrt(mean((summand - theta)^2) / N).
    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta  = `th'
    return scalar N      = `nn'
    return scalar nctry  = `nctry'
    return scalar ntreat = `ntr'
    return scalar se     = `sean'
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
    local n1 = r(N)
    local nctry1 = r(nctry)
    local ntreat1 = r(ntreat)
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2 = r(theta)
    local a2 = r(se)
    local n2 = r(N)
    local nctry2 = r(nctry)
    local ntreat2 = r(ntreat)
    local dh = `b1' - `b2'

    capture drop _pool
    quietly gen byte _pool = 0 if (`if1') | (`if2')
    quietly replace _pool = 1 if `d1' == 1
    quietly replace _pool = 2 if `d2' == 1

    tempname pf
    tempfile bf
    * t1/t2 kept, not just their difference -- DIAGNOSTIC-ONLY: lets the
    * caller compute each level's own bootstrap SE vs. the adopted analytic
    * SE (a1/a2). See 08b_aipw.do's _aipw header for the full rationale.
    quietly postfile `pf' double t1 double t2 double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
            local t1 = cond(_rc==0, r(theta), .)
            capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
            local t2 = cond(_rc==0, r(theta), .)
            if !missing(`t1') & !missing(`t2') quietly post `pf' (`t1') (`t2') (`t1' - `t2')
        restore
    }
    quietly postclose `pf'
    capture drop _pool
    local se = .
    local lo = .
    local hi = .
    local nd = 0
    local bse1 = .
    local bse2 = .
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
            local bse1 = r(sd)
            quietly summarize t2
            local bse2 = r(sd)
        }
    restore
    return scalar ok = 1
    return scalar dh = `dh'
    return scalar b1 = `b1'
    return scalar b2 = `b2'
    return scalar a1 = `a1'
    return scalar a2 = `a2'
    return scalar bse1 = `bse1'
    return scalar bse2 = `bse2'
    return scalar n1 = `n1'
    return scalar n2 = `n2'
    return scalar nctry1 = `nctry1'
    return scalar nctry2 = `nctry2'
    return scalar ntreat1 = `ntreat1'
    return scalar ntreat2 = `ntreat2'
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
    long nobs byte nctry ntreat ///
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
foreach ch in credit inv claims_govt claimsgov_assets claimpriv_assets real_lending fdi {
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

* govexp, pb, ca, claimsgov_assets, claimpriv_assets SILENCED below --
* scoped to exactly the project's six headline dependent variables (GDP
* [08b_aipw_analyticSE.do / the in-file gdp block below] + Investment, Bank
* credit, Claims on government, FDI, Real lending rate here), matching
* 13c_aipw_channels.do's own identical scoping. claimsgov_assets/
* claimpriv_assets are the two sovereign-bank nexus shares -- estimated
* separately in 13d_aipw_nexus_split_analyticSE.do for the exposure-
* heterogeneity design, not part of this file's own six-variable scope.
* Restore any of them by adding back to the active list below: "credit
* claims_govt inv govexp pb fdi claimsgov_assets claimpriv_assets ca".
* Their outcome construction (ch_v_h/pre_v/l_v) above still runs
* regardless, since it is cheap and shared -- only the estimation loop
* below is skipped for them.
foreach ch in credit claims_govt inv real_lending fdi {
* SILENCED (not part of the six headline variables -- uncomment to restore):
*             claimsgov_assets claimpriv_assets

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

    * BALANCED A-D SAMPLE (common_abcd, built in 18_transforms.do): Investment
    * and Bank credit (`inv'/`credit') plus Claims on government
    * (`claims_govt') are now ACTUALLY estimated on the same common episode
    * set as 08b_aipw.do's AIPW GDP result and 03_lp_resolution.do's Table 2
    * -- restricted to onsets where all four outcomes are non-missing at
    * every horizon h=0..4 (Asonuma et al.'s own "balanced panels"
    * convention). This feeds BOTH stages of _aipwpair (propensity and
    * outcome), the same as in 08b_aipw.do. claimsgov_assets,
    * claimpriv_assets, real_lending, and fdi keep their own best-available
    * sample, unrestricted.
    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

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
    di as result "  Act 2 (PAPER-ALIGNED SE):  h   ND (se_analytic)   DEF (se_analytic)   def-nd   [95% boot CI]   Clogg z    p"
    di as result "           se_analytic = the PAPER'S OWN analytic (unclustered influence-function) SE. 13c's own"
    di as result "           headline instead ADOPTS a row-bootstrap SE here (see 08b_aipw.do's header) -- this"
    di as result "           duplicate shows the paper-aligned number directly, side by side."
    di as result "           ND/DEF stars are the conventional t-test vs zero (b/se_analytic): * p<.10 ** p<.05 *** p<.01."
    post `R' ("`ch'") ("nd")  (0) (0) (0) (0) (0) (0) (0) (0)   // explicit baseline (h=0)
    post `R' ("`ch'") ("def") (0) (0) (0) (0) (0) (0) (0) (0)
    post `Rd' ("`ch'") (0) (0) (0) (0) (0) (0) (0) (0) (.) (.)   // explicit baseline (h=0)
    forvalues h = 0/4 {
        _aipwpair, y(ch_`ch'_`h') ///
            d1(onset_def) if1(sample==1 & onset_nd==0`balflag') ///
            d2(onset_nd)  if2(sample==1 & onset_def==0`balflag') ///
            omod(`om') pz(`om' `cz_def') reps(`nboot')
        if r(ok) {
            local B1 = r(b1)   // default-linked ATE
            local B2 = r(b2)   // non-default ATE
            local A1 = r(a1)   // analytic SE, default-linked -- PAPER-ALIGNED: now used for the level CI too
            local A2 = r(a2)   // analytic SE, non-default -- PAPER-ALIGNED: now used for the level CI too
            local DH = r(dh)
            local SE = r(se)
            local LO = r(lo)
            local HI = r(hi)
            local ND = r(nd)

            * nobs/nctry/ntreat: the regression sample (obs, countries,
            * treated episodes) behind THIS arm's own outcome model, from
            * _aipw's own e(sample) via _aipwpair's n1/n2, nctry1/nctry2,
            * ntreat1/ntreat2 -- ported from 13d_aipw_nexus_split.do's
            * identical nctry mechanism, extended with the episode count so
            * the merged table can print a full Observations/Countries/
            * Episodes triple per level row (this file's d1=default-linked,
            * d2=non-default, per the calling convention above).
            local N1 = r(n1)
            local N2 = r(n2)
            local NC1 = r(nctry1)
            local NC2 = r(nctry2)
            local NT1 = r(ntreat1)
            local NT2 = r(ntreat2)

            * PAPER-ALIGNED (this duplicate only): level CIs = theta +/-
            * 1.96*ANALYTIC SE, matching Asonuma et al.'s own construction
            * literally. 13c_aipw_channels.do's own headline departs from
            * this (row-bootstrap SE) after a diagnostic found the analytic
            * formula understated the def arm's true uncertainty by
            * 3.75-5.5x on ~20-episode default arms -- unchanged finding;
            * this file exists only to show the alternative side by side.
            post `R' ("`ch'") ("nd")  (`h'+1) (`B2') (`A2') (`B2'-1.96*`A2') (`B2'+1.96*`A2') (`N2') (`NC2') (`NT2')
            post `R' ("`ch'") ("def") (`h'+1) (`B1') (`A1') (`B1'-1.96*`A1') (`B1'+1.96*`A1') (`N1') (`NC1') (`NT1')

            * Clogg z: analytic SEs -- confirmed to match a line literally in their own replication script (see header); unchanged
            * from the headline -- it was already analytic-SE-based there).
            local zz = .
            local pz = .
            if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
                local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
                local pz = 2*(1 - normal(abs(`zz')))
            }
            post `Rd' ("`ch'") (`h'+1) (`DH') (`B1') (`B2') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')

            * PAPER-ALIGNED: conventional t-test for each level vs zero using
            * the analytic SE (b / own analytic SE), not the adopted bootstrap SE.
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
        }
        else di as error "    h=" `h'+1 ": Act 2 estimate failed (thin sample)."
    }
}

* ══════════════════════════════════════════════════════════════════════════
* GDP — self-contained headline estimate (PAPER-ALIGNED SE), mirroring the
* identical addition in 13c_aipw_channels.do. This is a SEPARATE BLOCK, not
* one more iteration of the `foreach ch of local channels' loop above, for
* the same reason as in the plain file: GDP's outcome is dy_`h' (built in
* 18_transforms.do), not ch_`ch'_`h', so folding it into the loop would mean
* special-casing the outcome-variable name inside every iteration for one
* variable out of six.
*
* Pulled from 08b_aipw_analyticSE.do's OWN Act 2 GDP loop (NOT the plain
* 08b_aipw.do) so this file's GDP block matches ITS OWN SE convention:
* level CIs and t-test stars use the ANALYTIC SEs (`A1'/`A2'), matching the
* channel loop above in this same duplicate file. Same omodel ($ctrl_core),
* pz (`cx' `cz_def'), reps(`nboot'), common_abcd restriction, _aipwpair
* signature, and Clogg z as both 08b_aipw_analyticSE.do and the channel
* loop above. Posted into THIS file's own `R'/`Rd' postfiles under
* channel=="gdp".
*
* REDUNDANCY, ACCEPTED DELIBERATELY (same reasoning as the plain file): GDP
* is now estimated in three places project-wide (08b_aipw.do, 08b_aipw_
* analyticSE.do, and here) -- accepted so this duplicate is runnable
* standalone. 08b_aipw_analyticSE.do itself is untouched by this addition.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== CHANNEL: gdp ==="
di as result "  (Estimated redundantly with 08b_aipw_analyticSE.do's own Act 2 GDP loop, to keep"
di as result "  this file self-contained -- see the block comment above. Same spec/seed/bootstrap"
di as result "  mechanics and the same PAPER-ALIGNED analytic-SE level bands as the channel loop above.)"
di as result "  Act 2 (PAPER-ALIGNED SE):  h   ND (se_analytic)   DEF (se_analytic)   def-nd   [95% boot CI]   Clogg z    p"
post `R' ("gdp") ("nd")  (0) (0) (0) (0) (0) (0) (0) (0)   // explicit baseline (h=0)
post `R' ("gdp") ("def") (0) (0) (0) (0) (0) (0) (0) (0)
post `Rd' ("gdp") (0) (0) (0) (0) (0) (0) (0) (0) (.) (.)   // explicit baseline (h=0)
forvalues h = 0/4 {
    _aipwpair, y(dy_`h') ///
        d1(onset_def) if1(sample==1 & onset_nd==0 & common_abcd==1) ///
        d2(onset_nd)  if2(sample==1 & onset_def==0 & common_abcd==1) ///
        omod($ctrl_core) pz(`cx' `cz_def') reps(`nboot')
    if r(ok) {
        local B1 = r(b1)     // default-linked ATE
        local B2 = r(b2)     // non-default ATE
        local A1 = r(a1)     // analytic SE, default-linked -- PAPER-ALIGNED: used for the level CI too
        local A2 = r(a2)     // analytic SE, non-default -- PAPER-ALIGNED: used for the level CI too
        local DH = r(dh)
        local SE = r(se)
        local LO = r(lo)
        local HI = r(hi)
        local ND = r(nd)

        * nobs/nctry/ntreat per arm -- see the channel loop's identical block
        * above for the full rationale (same _aipwpair mechanism).
        local N1 = r(n1)
        local N2 = r(n2)
        local NC1 = r(nctry1)
        local NC2 = r(nctry2)
        local NT1 = r(ntreat1)
        local NT2 = r(ntreat2)

        post `R' ("gdp") ("nd")  (`h'+1) (`B2') (`A2') (`B2'-1.96*`A2') (`B2'+1.96*`A2') (`N2') (`NC2') (`NT2')
        post `R' ("gdp") ("def") (`h'+1) (`B1') (`A1') (`B1'-1.96*`A1') (`B1'+1.96*`A1') (`N1') (`NC1') (`NT1')

        local zz = .
        local pz = .
        if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
            local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
            local pz = 2*(1 - normal(abs(`zz')))
        }
        post `Rd' ("gdp") (`h'+1) (`DH') (`B1') (`B2') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')

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
    }
    else di as error "    h=" `h'+1 ": Act 2 estimate failed (thin sample)."
}

* ══════════════════════════════════════════════════════════════════════════
* ROBUSTNESS: CREDIT CHANNEL AIPW, BULGARIA EXCLUDED (LEAVE-ONE-OUT)
* SILENCED -- costs its own 5-horizon x nboot bootstrap on top of the
* headline six-variable estimation above, and is diagnostic-only (does not
* feed any exported CSV/table/figure). Uncomment to restore. Matches
* 13c_aipw_channels.do's own identical silencing.
/*
di as result _n "=== ROBUSTNESS: credit channel AIPW, Bulgaria excluded (leave-one-out) ==="
di as result "  h   DEF, ex.Bulgaria (se_analytic)   DEF, full sample (se_analytic)"
local om_credit_lo1 l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
forvalues h = 0/4 {
    _aipwpair, y(ch_credit_`h') ///
        d1(onset_def) if1(sample==1 & onset_nd==0 & common_abcd==1 & country!="Bulgaria") ///
        d2(onset_nd)  if2(sample==1 & onset_def==0 & common_abcd==1 & country!="Bulgaria") ///
        omod(`om_credit_lo1') pz(`om_credit_lo1' `cz_def') reps(`nboot')
    if r(ok) {
        local B1_lo1  = r(b1)
        local A1_lo1  = r(a1)
        di "  h=" `h'+1 "   " %10.3f `B1_lo1' " (" %6.3f `A1_lo1' ")"
    }
    else di as error "  h=" `h'+1 ": leave-one-out estimate failed (thin sample)."
}
di as result "  Compare against credit's own DEF row printed under '=== CHANNEL: credit ==='"
di as result "  above (same B1/A1, full sample, both now analytic SE in this duplicate)."
di as result "  If ex.Bulgaria stays large and significant, the headline result is not"
di as result "  just Bulgaria; if it collapses toward zero, most of the def-arm signal"
di as result "  was one country's crisis."
*/

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
    export delimited "$tabs/aipw_channels_diff_analyticSE.csv", replace
    di as result _n "AIPW channel def-nd difference CSV saved: $tabs/aipw_channels_diff_analyticSE.csv"
    di as result "(the difference bootstrap itself is UNCHANGED from 13c_aipw_channels.do --"
    di as result " only written under a separate name so this file never overwrites the original)"
restore

use "`resf'", clear
label var b  "AIPW ATE (pp of the channel ratio)"
label var se "Analytic (unclustered influence-function) SE, matching the paper's own formula (this duplicate only)"
label var lo "95% CI lower = theta-1.96*se (analytic)"
label var hi "95% CI upper = theta+1.96*se (analytic)"
label var nobs   "Observations in this cell's own AIPW outcome-regression sample"
label var nctry  "Countries in this cell's own AIPW outcome-regression sample"
label var ntreat "Treated onsets in this cell's own AIPW outcome-regression sample (Episodes)"
order channel series horizon b se lo hi nobs nctry ntreat
export delimited "$tabs/aipw_channels_analyticSE.csv", replace
di as result _n "AIPW channel results CSV saved: $tabs/aipw_channels_analyticSE.csv"

* ── Figure A: Act 1 all-crises AIPW per channel — SILENCED, see header. ─────
/*
local c1 "blue"
local channels_ord   credit claims_govt inv claimsgov_assets claimpriv_assets
local titlelabels_ch `" "Bank credit" "Bank claims on government" "Investment" "Bank claims on government / assets" "Bank claims on private sector / assets" "'
local i = 1
foreach ch of local channels_ord {
    local tlab : word `i' of `titlelabels_ch'
    * Y-axis title shown only on the leftmost panel of each row (cols(3)
    * rows(2): panels 1 and 4) -- not repeated on every panel.
    local ytit ""
    if inlist(`i', 1, 4) local ytit "Cumulative percent change"
    capture twoway ///
        (rarea lo hi horizon if series=="all" & channel=="`ch'", color("`c1'%18") lwidth(none)) ///
        (connected b horizon if series=="all" & channel=="`ch'", lcolor("`c1'") lwidth(medthick) msymbol(circle)), ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)5, labsize(medium)) ylabel(, labsize(medium) angle(horizontal)) ///
        xtitle("Year", size(medium)) ///
        ytitle("`ytit'", size(medsmall)) ///
        title("`tlab'", size(medlarge) color(navy)) legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(aipwch_`i', replace)
    local ++i
}
capture graph combine aipwch_1 aipwch_2 aipwch_3 aipwch_4 aipwch_5, ///
    cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
if _rc == 0 {
    graph export "$figs/fig_aipw_ch_act1.pdf", replace
    di as result "Figure saved: fig_aipw_ch_act1.pdf"
}
else di as error "  ** fig_aipw_ch_act1 failed (rc=" _rc ")"
forvalues i = 1/5 {
    capture graph drop aipwch_`i'
}
*/

* ── Figure B: Act 2 resolution split per channel ─────────────────────────
* Same construction pattern as 12_channels_resolution.do's multi-panel
* figure: each channel is its own named twoway, combined via graph combine,
* rather than Stata's by() faceting (which gives less control over per-
* panel spacing/sizing and reads differently from the rest of the project's
* figures). UNIFORM IRF STYLE: non-default = blue, default-linked = red,
* both solid, markers match line color, no legend.
local c_nd  "blue"
local c_def "red"
* channels_ord scoped to the active estimation loop's five channels,
* matching 13c_aipw_channels.do's own identical fix.
local channels_ord   credit claims_govt inv real_lending fdi
local titlelabels_ch `" "Bank credit" "Bank claims on government" "Investment" "Real lending rate" "FDI" "'
local i = 1
foreach ch of local channels_ord {
    local tlab : word `i' of `titlelabels_ch'
    * Y-axis title shown only on the leftmost panel of each row (cols(3)
    * rows(2): panels 1 and 4) -- not repeated on every panel.
    local ytit ""
    if inlist(`i', 1, 4) local ytit "Cumulative percent change"
    capture twoway ///
        (rarea lo hi horizon if series=="nd"  & channel=="`ch'", color("`c_nd'%16")  lwidth(none)) ///
        (rarea lo hi horizon if series=="def" & channel=="`ch'", color("`c_def'%16") lwidth(none)) ///
        (connected b horizon if series=="nd"  & channel=="`ch'", lcolor("`c_nd'")  lwidth(medthick) msymbol(circle)) ///
        (connected b horizon if series=="def" & channel=="`ch'", lcolor("`c_def'") lwidth(medthick) msymbol(square)), ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)5, labsize(medium)) ylabel(, labsize(medium) angle(horizontal)) ///
        xtitle("Year", size(medium)) ///
        ytitle("`ytit'", size(medsmall)) ///
        title("`tlab'", size(medlarge) color(navy)) legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(aipwch2_`i', replace)
    local ++i
}
capture graph combine aipwch2_1 aipwch2_2 aipwch2_3 aipwch2_4 aipwch2_5, ///
    cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
if _rc == 0 {
    graph export "$figs/fig_aipw_ch_act2_analyticSE.pdf", replace
    di as result "Figure saved: fig_aipw_ch_act2_analyticSE.pdf"
}
else di as error "  ** fig_aipw_ch_act2_analyticSE failed (rc=" _rc ")"
forvalues i = 1/5 {
    capture graph drop aipwch2_`i'
}

* ══════════════════════════════════════════════════════════════════════════
* Combined 6-panel figure (PAPER-ALIGNED SE): Panel A GDP, B Investment,
* C Bank credit, D Claims on government, E FDI, F Real lending rate. GDP is
* now estimated IN THIS FILE (the GDP block above, analytic-SE bands), posted
* into `resf' under channel=="gdp" -- no import of 08b_aipw_analyticSE.do's
* aipw_results_analyticSE.csv, so this figure (and this file as a whole) no
* longer depends on 08b_aipw_analyticSE.do having run first. Every panel is
* still analytic-SE throughout, since the GDP block above uses the same
* analytic-SE convention as the channel loop. 08b_aipw_analyticSE.do is
* unchanged and keeps producing its own separate outputs.
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`resf'", clear

    local combo_vars   gdp inv credit claims_govt fdi real_lending
    local combo_labels `" "Panel A: GDP" "Panel B: Investment" "Panel C: Bank credit" "Panel D: Claims on govt" "Panel E: FDI" "Panel F: Real lending rate" "'
    local i = 1
    foreach cv of local combo_vars {
        local clab : word `i' of `combo_labels'
        local ytit ""
        if inlist(`i', 1, 4) local ytit "Cumulative percent change"
        capture twoway ///
            (rarea lo hi horizon if series=="nd"  & channel=="`cv'", color("`c_nd'%16")  lwidth(none)) ///
            (rarea lo hi horizon if series=="def" & channel=="`cv'", color("`c_def'%16") lwidth(none)) ///
            (connected b horizon if series=="nd"  & channel=="`cv'", lcolor("`c_nd'")  lwidth(medthick) msymbol(circle)) ///
            (connected b horizon if series=="def" & channel=="`cv'", lcolor("`c_def'") lwidth(medthick) msymbol(square)), ///
            yline(0, lpattern(dash) lcolor(gs8)) ///
            xlabel(0(1)5, labsize(large)) ylabel(, labsize(large) angle(horizontal)) ///
            xtitle("Year", size(large)) ///
            ytitle("`ytit'", size(large)) ///
            title("`clab'", size(medlarge) color(navy)) legend(off) ///
            graphregion(color(white)) plotregion(color(white)) ///
            name(combA_`i', replace)
        local ++i
    }
    capture graph combine combA_1 combA_2 combA_3 combA_4 combA_5 combA_6, ///
        cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
    if _rc == 0 {
        graph export "$figs/fig_aipw_combined_analyticSE.pdf", replace
        di as result "Figure saved: fig_aipw_combined_analyticSE.pdf (Panel A-F, paper-aligned analytic SE)"
    }
    else di as error "  ** fig_aipw_combined_analyticSE failed (rc=" _rc ")"
    forvalues i = 1/6 {
        capture graph drop combA_`i'
    }
restore

* ══════════════════════════════════════════════════════════════════════════
* MERGED SIX-VARIABLE TABLE (PAPER-ALIGNED SE): GDP + the five channels, ONE
* RTF/CSV, one panel per variable -- mirrors 13c_aipw_channels.do's own
* merged table exactly (same hand-built `file write' mechanism, chosen there
* over `ereturn post'+esttab because the level bands, the bootstrap-CI
* difference row, and the Clogg z sit on different SE bases that do not
* collapse into one synthesized (b,V) pair -- see that file's header for the
* full reasoning), with the ONE substantive change: level rows report the
* ANALYTIC SE (this duplicate's own convention throughout), not the adopted
* row-bootstrap SE. The difference row (bootstrap CI) and the Clogg z are
* unchanged either way, since both were already on their respective bases
* in the plain file too.
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`diffresf'", clear
    tempfile _diffuse
    save `_diffuse'
    use "`resf'", clear
    tempfile _resuse
    save `_resuse'

    capture file close aipwtab
    file open aipwtab using "$tabs/aipw_combined_six_resolution_analyticSE.rtf", write replace
    file write aipwtab "{\rtf1\ansi\deff0" _n
    file write aipwtab "{\b Table: AIPW output/channel cost by crisis resolution, all six headline variables (PAPER-ALIGNED SE)\par}" _n
    file write aipwtab "{\i Doubly-robust AIPW (Asonuma et al. 2024 Eq. 3 / Jorda-Taylor 2016), IPWRA. Each panel is one" _n
    file write aipwtab " outcome (GDP + the five transmission channels), two level rows (non-default vs tranquil," _n
    file write aipwtab " default-linked vs tranquil) and a difference row (def - nd). Level SEs and CIs here use the" _n
    file write aipwtab " PAPER'S OWN ANALYTIC (unclustered influence-function) SE -- this file's one deliberate" _n
    file write aipwtab " departure from 13c_aipw_channels.do's own aipw_combined_six_resolution.rtf, which instead" _n
    file write aipwtab " ADOPTS a row-bootstrap SE for the level bands (see that file's header for the diagnostic that" _n
    file write aipwtab " motivated the switch: the analytic formula understated the def arm's true uncertainty by" _n
    file write aipwtab " 3.75-5.5x on this project's thin default arm). Each level row is followed by its own" _n
    file write aipwtab " Observations/Countries/Episodes line (obs./countries from that arm's own AIPW outcome-" _n
    file write aipwtab " regression sample; episodes = treated onsets). Level stars are the conventional t-test vs" _n
    file write aipwtab " zero (b/se_analytic). The difference row is shown as [bootstrap 95% CI] on one line (its *" _n
    file write aipwtab " marks the CI excluding zero -- the governing test for the difference; the bootstrap itself" _n
    file write aipwtab " is UNCHANGED from the plain file) and the coefficient with stars on the next, with Clogg et" _n
    file write aipwtab " al. (1995)'s z (own analytic SEs, confirmed from their replication script) appended in" _n
    file write aipwtab " parentheses as the permissive companion statistic, unchanged from the plain file since it" _n
    file write aipwtab " was already analytic-SE-based there. GDP, Investment, Bank credit and Claims on government" _n
    file write aipwtab " are estimated on the balanced common_abcd sample; FDI and Real lending rate keep their own" _n
    file write aipwtab " best-available sample.\par}" _n
    file write aipwtab "\par" _n
    file write aipwtab "\tab h = 1\tab h = 2\tab h = 3\tab h = 4\tab h = 5\par" _n
    file write aipwtab "\par" _n

    local combo_vars   gdp inv credit claims_govt fdi real_lending
    local combo_labels `" "GDP" "Investment" "Bank credit" "Bank claims on government" "FDI" "Real lending rate" "'
    local i = 1
    foreach cv of local combo_vars {
        local clab : word `i' of `combo_labels'
        file write aipwtab "{\b `clab'}\par" _n

        local ndline ""
        local ndocline ""
        local defline ""
        local defocline ""
        forvalues h = 1/5 {
            use `_resuse', clear
            quietly summarize b if channel=="`cv'" & series=="nd" & horizon==`h', meanonly
            local bnd = r(mean)
            quietly summarize se if channel=="`cv'" & series=="nd" & horizon==`h', meanonly
            local send = r(mean)
            quietly summarize nobs if channel=="`cv'" & series=="nd" & horizon==`h', meanonly
            local ondn = r(mean)
            quietly summarize nctry if channel=="`cv'" & series=="nd" & horizon==`h', meanonly
            local cndn = r(mean)
            quietly summarize ntreat if channel=="`cv'" & series=="nd" & horizon==`h', meanonly
            local endn = r(mean)
            quietly summarize b if channel=="`cv'" & series=="def" & horizon==`h', meanonly
            local bdef = r(mean)
            quietly summarize se if channel=="`cv'" & series=="def" & horizon==`h', meanonly
            local sedef = r(mean)
            quietly summarize nobs if channel=="`cv'" & series=="def" & horizon==`h', meanonly
            local odef = r(mean)
            quietly summarize nctry if channel=="`cv'" & series=="def" & horizon==`h', meanonly
            local cdef = r(mean)
            quietly summarize ntreat if channel=="`cv'" & series=="def" & horizon==`h', meanonly
            local edef = r(mean)

            local tnd = cond(`send'>0 & !missing(`send'), `bnd'/`send', .)
            local pnd = cond(!missing(`tnd'), 2*(1-normal(abs(`tnd'))), .)
            local sgnd = cond(missing(`pnd'), "", cond(`pnd'<.01,"***",cond(`pnd'<.05,"**",cond(`pnd'<.10,"*",""))))
            local tdef = cond(`sedef'>0 & !missing(`sedef'), `bdef'/`sedef', .)
            local pdef = cond(!missing(`tdef'), 2*(1-normal(abs(`tdef'))), .)
            local sgdef = cond(missing(`pdef'), "", cond(`pdef'<.01,"***",cond(`pdef'<.05,"**",cond(`pdef'<.10,"*",""))))

            local ndl : di %7.3f `bnd'
            local defl : di %7.3f `bdef'
            local ndline "`ndline'\tab `ndl'`sgnd' (`: di %5.3f `send'')"
            local defline "`defline'\tab `defl'`sgdef' (`: di %5.3f `sedef'')"

            local ondns : di %5.0f `ondn'
            local cndns : di %3.0f `cndn'
            local endns : di %4.0f `endn'
            local odefs : di %5.0f `odef'
            local cdefs : di %3.0f `cdef'
            local edefs : di %4.0f `edef'
            local ndocline  "`ndocline'\tab `ondns'/`cndns'/`endns'"
            local defocline "`defocline'\tab `odefs'/`cdefs'/`edefs'"
        }
        file write aipwtab "Non-default`ndline'\par" _n
        file write aipwtab "Observations/Countries/Episodes`ndocline'\par" _n
        file write aipwtab "Default-linked`defline'\par" _n
        file write aipwtab "Observations/Countries/Episodes`defocline'\par" _n

        file write aipwtab "{\i Differences between the estimated coefficients,\line [bootstrap 95% CI] (Clogg et al.'s z)}\par" _n
        local ciline ""
        local zline ""
        forvalues h = 1/5 {
            use `_diffuse', clear
            quietly summarize dhl if channel=="`cv'" & horizon==`h', meanonly
            local dhl = r(mean)
            quietly summarize lo if channel=="`cv'" & horizon==`h', meanonly
            local dlo = r(mean)
            quietly summarize hi if channel=="`cv'" & horizon==`h', meanonly
            local dhi = r(mean)
            quietly summarize nd if channel=="`cv'" & horizon==`h', meanonly
            local dndraws = r(mean)
            quietly summarize cloggz if channel=="`cv'" & horizon==`h', meanonly
            local dzz = r(mean)
            quietly summarize cloggp if channel=="`cv'" & horizon==`h', meanonly
            local dpz = r(mean)

            local dsig = cond(`dndraws'>=50 & !missing(`dlo') & (`dlo'>0 | `dhi'<0), "*", "")
            local dl : di %7.3f `dhl'
            local dlos : di %6.2f `dlo'
            local dhis : di %5.2f `dhi'
            local zl : di %6.2f `dzz'
            local zsg = cond(missing(`dpz'), "", cond(`dpz'<.01,"***",cond(`dpz'<.05,"**",cond(`dpz'<.10,"*",""))))
            local ciline "`ciline'\tab [`dlos', `dhis']`dsig'"
            local zline  "`zline'\tab `dl'`dsig' (`zl'`zsg')"
        }
        file write aipwtab "`ciline'\par" _n
        file write aipwtab "`zline'\par" _n
        file write aipwtab "\par" _n
        local ++i
    }
    file write aipwtab "{\i Level rows: Observations/Countries/Episodes from that arm's own AIPW outcome-regression" _n
    file write aipwtab " sample. * p<0.10, ** p<0.05, *** p<0.01 (levels: t-test vs zero on analytic SE;" _n
    file write aipwtab " difference bracket: row-bootstrap 95% percentile CI excludes zero, the governing test;" _n
    file write aipwtab " Clogg z stars, in parentheses on the coefficient line: permissive companion statistic," _n
    file write aipwtab " assumes independence).\par}" _n
    file write aipwtab "}" _n
    file close aipwtab
    di as result "AIPW combined six-variable table saved (paper-aligned SE): $tabs/aipw_combined_six_resolution_analyticSE.rtf"

    * ── Plain merged dataset: levels (resf) + the matching difference-row
    * columns broadcast onto both the nd/def rows of the same channel x
    * horizon -- both objects were already in memory from the loop above.
    use `_diffuse', clear
    rename se dse
    rename lo dlo
    rename hi dhi
    rename nd ndraws_diff
    tempfile _diffmerge
    save `_diffmerge'
    use `_resuse', clear
    merge m:1 channel horizon using `_diffmerge', nogenerate
    label var b  "AIPW ATE (pp; GDP) or (pp of the channel ratio; channels)"
    label var se "Analytic (unclustered influence-function) SE, matching the paper's own formula (this duplicate only)"
    label var lo "95% CI lower = b - 1.96*se (analytic)"
    label var hi "95% CI upper = b + 1.96*se (analytic)"
    label var dhl "AIPW def - nd difference (pp), row bootstrap"
    label var dse "Bootstrap SE of the difference"
    label var dlo "95% percentile CI lower, difference (row bootstrap)"
    label var dhi "95% percentile CI upper, difference (row bootstrap)"
    label var ndraws_diff "Valid bootstrap draws, difference"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    order channel series horizon b se lo hi dhl dse dlo dhi ndraws_diff cloggz cloggp
    export delimited "$tabs/aipw_combined_six_analyticSE.csv", replace
    di as result "AIPW combined six-variable CSV saved (paper-aligned SE): $tabs/aipw_combined_six_analyticSE.csv"
restore

di as result _n "13c_aipw_channels_analyticSE.do complete (paper-aligned SE duplicate)."
di as result "Compare aipw_channels_analyticSE.csv / fig_aipw_ch_act2_analyticSE.pdf against"
di as result "13c_aipw_channels.do's own aipw_channels.csv / fig_aipw_ch_act2.pdf directly:"
di as result "same point estimates, level SE/CI swapped from row-bootstrap to the paper's"
di as result "own analytic formula. The def-nd difference (aipw_channels_diff_analyticSE.csv)"
di as result "is UNCHANGED from the headline -- that already matched the paper's own"
di as result "bootstrap device."
