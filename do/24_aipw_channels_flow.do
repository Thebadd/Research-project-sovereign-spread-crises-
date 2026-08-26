/*===========================================================================
  24_AIPW_CHANNELS_FLOW.DO
  Doubly-robust AIPW (Asonuma et al. Eq. 3) for the transmission channels,
  FLOW treatment, resolution split ONLY — the flow analogue of 13c's Act 2,
  built on 23_channels_resolution_flow.do's channel construction and
  21_aipw_flow.do's (now-fixed) propensity model.

  WHY THIS FILE EXISTS, AND WHY IT IS SCOPED THE WAY IT IS
  ----------------------------------------------------------
  13c_aipw_channels.do reports two things per channel under ONSET coding:
  Act 1 (pooled, all crises vs tranquil) and Act 2 (the resolution split —
  non-default vs tranquil, default-linked vs tranquil, plus the def-nd
  difference via paired bootstrap). This file is the FLOW-treatment
  counterpart of Act 2 ONLY. No Act 1 (pooled flow AIPW) is built here — flow
  coding never gets a pooled AIPW line anywhere in this project (20/21/22/23
  all stop at nd/def), and this file follows that pattern rather than
  breaking it.

  This file is explicitly the flow AIPW analogue of
  23_channels_resolution_flow.do (channels BY RESOLUTION TYPE), not of
  22_channels_flow.do (pooled in_crisis channels). The channel outcome
  construction (ch_<var>_h, epc_pre_<var>) is identical between 22 and 23 —
  copied here from 23, since 23 is already organised around in_crisis_nd /
  in_crisis_def as two separate arms, the same organising logic this AIPW
  file needs.

  Six channels only: credit, claims_govt, inv, govexp, pb, fdi — the ones
  22/23 already built flow outcomes for. The nexus/ca channels 13c adds
  beyond 11/12 (claimsgov_assets, claimpriv_assets, ca) have no flow-outcome
  construction anywhere in this project yet; adding them here would be new
  scope, not a port of existing work, so they are left out.

  THE PROPENSITY MODEL IS NOT RE-DERIVED HERE — IT IS 21_AIPW_FLOW.DO'S,
  UNCHANGED
  ----------------------------------------------------------------------
  21_aipw_flow.do's header documents (as settled history, the diagnostic code
  itself no longer kept there) why the flow AIPW's propensity model
  near-separates on continuation rows and the fix adopted: Eq. (2) is fit
  ONLY on tranquil and onset rows (continuation==0), and the fitted model is
  then EXTRAPOLATED to every row (continuation included) because Eq. (1)/Eq.
  (3) keep the full flow-coded treatment. This file's `_aipw' program is
  21's, copied as it stands today. See 21_aipw_flow.do's header for the
  justification; it is referenced here, not re-derived.

  PREDICTOR CHANGE: past_def_onsets -> years_since_def_onset. This file
  previously carried a diagnostic sequence (Sections 1a/1b) establishing that
  the def arm's IPW weights were severely concentrated (top 5% of rows =
  98.9% of the AIPW summand's variance) with past_def_onsets (a running COUNT
  that never resets) as the leading suspect, and testing years_since_def_onset
  (years since the country's most recent PRIOR default-linked onset) as a
  replacement: individually significant, sensible sign, but it did NOT fix the
  weight concentration (98.9% -> 98.6%, essentially unchanged) -- adopted for
  being a better, independently-motivated predictor, not as a variance fix,
  which remains open. That decision is settled and this file now always uses
  years_since_def_onset (via cz_recency); the diagnostic code establishing it
  is not kept here (see git history if it is ever needed again).

  Only the OUTCOME model (Eq. 1's `omodel') is channel-specific, exactly as
  in 13c: the propensity model is the SAME object regardless of which channel
  is being explained, because selection into a spread crisis does not depend
  on which downstream variable you are about to measure.

  SPECIFICATION
  -------------
  Eq. (1)  ch_<var>_h on D + episode-dated controls (23's ctrl_<ch>: $ctrl_flow
           plus the channel's own epc_pre_<var>, with the core term equal to
           the channel's own level dropped for credit/govexp), country FE,
           IPW-weighted (IPWRA, matching 21).
  Eq. (2)  Pr(D=1 | $ctrl_core, l_fedfunds, l_reg_crisis_share, past_def_onsets),
           fit on tranquil+onset rows only, extrapolated to continuation rows.
  Eq. (3)  The paper's exact AIPW form, as in 21_aipw_flow.do's `_aipw'.

  Each type scored against TRANQUIL years with the rival type dropped (21's
  sample_for*-style design). The def-nd DIFFERENCE is bootstrapped directly
  on paired resamples (row-level `bsample, strata()', 21's adopted baseline),
  not the Clogg z — same reasoning as 20/23: under flow coding a country can
  contribute treated rows to BOTH arms, so the two cells are not independent.

  INFERENCE, MATCHING 21_AIPW_FLOW.DO'S CONVENTION EXACTLY
  ----------------------------------------------------------
  LEVELS (nd/def vs tranquil): analytic influence-function SE from `_aipw',
  bands theta +/- 1.96*se — NOT their own bootstrap. DIFFERENCE (def-nd):
  1000... here nboot draws, percentile CI. This is 21's "faithful" route
  (08b_aipw.do's own file takes the other, "conservative" route — bootstrap
  for levels too — see 21's header for why the two files differ on this).

  The def-nd difference also carries a Clogg et al. (1995) z as a companion
  statistic (`aipw_channels_flow_diff.csv`'s cloggz/cloggp columns), matching
  21_aipw_flow.do's Section 3 and the reference paper's own Table 3, which
  reports both a bootstrap CI and a Clogg z for every between-cell
  difference. The Clogg z is the PERMISSIVE number: it assumes the two cells
  are independent, which they are not here (they share the tranquil control
  pool, and under flow coding a country can contribute treated rows to both
  arms), so the bootstrap CI governs where the two disagree — the z is
  reported for comparability with the reference paper's table format, not
  because it is the more defensible statistic on this project's thinner
  cells.

  NOT BUILT: the country-cluster bootstrap comparison (21's Section 3b). Can
  be added as a follow-up per channel if wanted; out of scope for this file.

  Outputs
  -------
    "$tabs/aipw_channels_flow.csv"       channel x series(nd/def) x horizon
    "$tabs/aipw_channels_flow_diff.csv"  channel x horizon def-nd gap
    "$figs/fig_aipw_ch_flow_act2.pdf/.png"  2x3 grid, nd/def overlay

  RUNTIME: 6 channels x 5 horizons x (2 levels + 1 paired-bootstrap diff) x
  nboot draws. At nboot=300 this is substantially heavier than 21 alone;
  raise to 500+ only for a final run. Self-contained: reads only
  $clean/panel_lp.dta.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow ep_seq years_since_def_onset l_contagion_dist {
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

* ══════════════════════════════════════════════════════════════════════════
* 1. CHANNEL OUTCOMES — identical construction to 22/23
* ══════════════════════════════════════════════════════════════════════════
local channels credit claims_govt inv govexp pb fdi

foreach var of local channels {
    local src `var'
    if inlist("`var'","credit","inv","govexp") local src ln_r_`var'
    capture drop `var'_base
    gen `var'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`src' - `var'_base
    }
    capture drop pre_`var'
    gen pre_`var' = L.`src' - L2.`src'
    capture drop epc_pre_`var' _ent_pre_`var'
    quietly bysort cid ep_seq: egen double _ent_pre_`var' = ///
        max(cond(onset_all==1, pre_`var', .))
    quietly gen double epc_pre_`var' = cond(in_crisis==1, _ent_pre_`var', pre_`var')
    quietly drop _ent_pre_`var'
    * See 22_channels_flow.do for why this sort is needed inside the loop.
    sort cid year
}

* bysort above re-sorts the physical dataset by its by-list as a side effect
* (independent of what xtset declared), leaving it sorted by cid ep_seq rather
* than cid year. xtscc/probit/bsample do their own sort/panel checks and error
* if this is not restored before the next panel command.
sort cid year

local epc_lc epc_l_credit_bank
local epc_lg epc_l_govexp

* CONTROL SET. flow_ctrl_variant: 0 = the ADOPTED flow-tier baseline,
* $ctrl_flow, built in 18_transforms.do from $ctrl_core_flowbase
* (l_banking_duration -> the l_banking_crisis DUMMY, l_ca -> tot_chg,
* l_hyperinfl -> l_lninfl -- see that file's "ADOPTED FLOW-TIER CORE CONTROL
* SET"). 1 = the adopted core PLUS the reference paper's own additional
* predictors, exchange2 and l_imf (18_transforms.do's "ALTERNATE FLOW
* CONTROL SET"). Default 0. This drives the channel outcome-model controls
* (`ctrl_<ch>') below.
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
if `flow_ctrl_variant'==1 local ctrl_flow_base $ctrl_flow_flowplus
else                       local ctrl_flow_base $ctrl_flow

local ctrl_credit      : list ctrl_flow_base - epc_lc
local ctrl_credit      `ctrl_credit' epc_pre_credit
local ctrl_claims_govt `ctrl_flow_base' epc_pre_claims_govt
local ctrl_inv         `ctrl_flow_base' epc_pre_inv
local ctrl_govexp      : list ctrl_flow_base - epc_lg
local ctrl_govexp      `ctrl_govexp' epc_pre_govexp
local ctrl_pb          `ctrl_flow_base' epc_pre_pb
local ctrl_fdi         `ctrl_flow_base' epc_pre_fdi

* ── Propensity model — SAME as 21_aipw_flow.do ──────────────────────────────
* `cx_active' is the row-dated propensity control set: the flow tier's
* adopted core control set (or the alternate, if toggled). Same
* construction as 21_aipw_flow.do's.
local cx_active = cond(`flow_ctrl_variant'==1, "$ctrl_core_flowplus", "$ctrl_core_flowbase")  // Eq. (2): row-dated, ADOPTED set
* cz_recency: the adopted predictor set -- l_contagion_dist (distance-weighted,
* CEPII great-circle) in place of l_reg_crisis_share, years_since_def_onset in
* place of past_def_onsets. See 21_aipw_flow.do's header for the full note.
local cz_recency l_fedfunds l_contagion_dist years_since_def_onset

set seed 20260819
local nboot = 300

sort cid year   // belt-and-suspenders: guaranteed panel order before estimation

* ══════════════════════════════════════════════════════════════════════════
* PROGRAMS — copied from 21_aipw_flow.do AS IT STANDS TODAY. _aipw's
* propensity-fit-sample restriction (continuation==0, extrapolated) is NOT
* re-derived here; see that file's header for the justification.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'

    tempvar xb m0 m1 ps summ iwt
    * PROPENSITY FIT SAMPLE EXCLUDES CONTINUATION ROWS -- see 21_aipw_flow.do's
    * header, "ESTIMATOR CHANGE, SECTION 1f ONWARD". Fit on tranquil+onset
    * only, then EXTRAPOLATE phat to every row of `touse', continuation rows
    * included, because Eq. (1) and Eq. (3) still use the FULL flow-coded `D'.
    quietly probit `D' `pmodel' if `touse' & continuation==0
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01                & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'

    * Eq. (1) — IPW-weighted OLS with country FE -> conditional means m1, m0.
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

* Paired bootstrap of BOTH cells and their difference on the SAME resample —
* copied from 21_aipw_flow.do's _aipwpairflow, unchanged. omod()/pz() are
* passed in per call, so this is directly reusable across channels without
* modification: the channel varies the outcome and omod(), never the program.
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

    if "`boot'" == "row" {
        capture drop _pool
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
                quietly keep if !missing(_pool)
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

sort cid year   // belt-and-suspenders before Section 2's estimation loop

* ══════════════════════════════════════════════════════════════════════════
* 2. ESTIMATION — loop over channels x horizons; levels + def-nd difference
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str24 channel str4 series byte horizon double b double se double lo double hi ///
    using "`resf'", replace

tempname Rd
tempfile diffresf
postfile `Rd' str24 channel byte horizon double dhl double bdef double bnd ///
    double se double lo double hi long nd double cloggz double cloggp ///
    using "`diffresf'", replace

foreach ch of local channels {
    foreach g in nd def {
        foreach m in b lo hi {
            matrix `m'_`ch'_`g' = J(7, 1, .)
            matrix `m'_`ch'_`g'[2,1] = 0
        }
    }
}

di as result _n "════════════════════════════════════════════════════════════"
di as result "FLOW AIPW CHANNELS BY RESOLUTION TYPE (Year 1 = the crisis year)"
di as result "Levels: analytic influence-function SE. Difference: paired bootstrap."
di as result "════════════════════════════════════════════════════════════"

foreach ch of local channels {
    di as result _n "--- CHANNEL: `ch' ---"
    post `R' ("`ch'") ("nd")  (0) (0) (0) (0) (0)
    post `R' ("`ch'") ("def") (0) (0) (0) (0) (0)
    post `Rd' ("`ch'") (0) (0) (0) (0) (0) (0) (0) (0) (.) (.)

    forvalues h = 0/4 {
        local hd  = `h' + 1
        local row = `h' + 3

        _aipwpairflow, y(ch_`ch'_`h') ///
            d1(in_crisis_def) if1(sample_flow==1 & in_crisis_nd==0) ///
            d2(in_crisis_nd)  if2(sample_flow==1 & in_crisis_def==0) ///
            omod(`ctrl_`ch'') pz(`cx_active' `cz_recency') reps(`nboot') boot(row)

        if !r(ok) {
            di as error "  h=`hd': estimate failed for `ch' (cell too thin)."
            continue
        }

        local B1 = r(b1)     // def level ATE
        local B2 = r(b2)     // nd  level ATE
        local A1 = r(a1)     // analytic SE, def
        local A2 = r(a2)     // analytic SE, nd
        local DH = r(dh)
        local SE = r(se)
        local LO = r(lo)
        local HI = r(hi)
        local ND = r(nd)

        matrix b_`ch'_def[`row',1]  = `B1'
        matrix lo_`ch'_def[`row',1] = `B1' - 1.96*`A1'
        matrix hi_`ch'_def[`row',1] = `B1' + 1.96*`A1'
        matrix b_`ch'_nd[`row',1]   = `B2'
        matrix lo_`ch'_nd[`row',1]  = `B2' - 1.96*`A2'
        matrix hi_`ch'_nd[`row',1]  = `B2' + 1.96*`A2'

        post `R' ("`ch'") ("def") (`hd') (`B1') (`A1') (`B1'-1.96*`A1') (`B1'+1.96*`A1')
        post `R' ("`ch'") ("nd")  (`hd') (`B2') (`A2') (`B2'-1.96*`A2') (`B2'+1.96*`A2')

        * Clogg et al. (1995) z, built from the same analytic influence-function
        * SEs used for the level bands -- identical construction to
        * 21_aipw_flow.do's Section 3. It is the PERMISSIVE statistic: it treats
        * the two cells as independent, and they are not, since they share the
        * tranquil control pool. The bootstrap percentile CI is the conservative
        * number and governs where the two disagree.
        local zz = .
        local pz = .
        if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
            local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
            local pz = 2*(1 - normal(abs(`zz')))
        }

        local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
        di as result "  h=" %1.0f `hd' "  ND=" %8.3f `B2' " (" %5.3f `A2' ")" ///
             "  DEF=" %8.3f `B1' " (" %5.3f `A1' ")" ///
             "  diff=" %8.3f `DH' " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
             "  " %4.0f `ND' "/`nboot'" ///
             "  Clogg z=" %7.3f `zz' " p=" %5.3f `pz'

        post `Rd' ("`ch'") (`hd') (`DH') (`B1') (`B2') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')
    }
}
postclose `R'
postclose `Rd'

* ══════════════════════════════════════════════════════════════════════════
* 3. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`diffresf'", clear
    label var dhl  "AIPW flow def - nd channel gap (pp)"
    label var bdef "Default-linked ATE (analytic SE band)"
    label var bnd  "Non-default ATE (analytic SE band)"
    label var se   "Bootstrap SD of the difference"
    label var lo   "95% percentile CI lower (bootstrap)"
    label var hi   "95% percentile CI upper (bootstrap)"
    label var nd   "Valid bootstrap draws"
    gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
    label var sig95 "Gap CI excludes 0"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    order channel horizon dhl bdef bnd se lo hi nd sig95 cloggz cloggp
    export delimited "$tabs/aipw_channels_flow_diff.csv", replace
    di as result _n "AIPW flow channel def-nd difference CSV saved: $tabs/aipw_channels_flow_diff.csv"
restore

use "`resf'", clear
label var b  "AIPW ATE (channel units, cumulative change)"
label var se "Analytic influence-function SE"
label var lo "Level 95% CI lower (theta +/- 1.96*se)"
label var hi "Level 95% CI upper"
order channel series horizon b se lo hi
export delimited "$tabs/aipw_channels_flow.csv", replace
di as result _n "AIPW flow channel results CSV saved: $tabs/aipw_channels_flow.csv"

* ── IRF datasets + figure ────────────────────────────────────────────────
foreach ch of local channels {
    foreach g in nd def {
        preserve
            clear
            set obs 7
            gen horizon = _n - 2
            foreach m in b lo hi {
                svmat `m'_`ch'_`g', names(`m')
                rename `m'1 `m'
            }
            gen series = "`g'"
            save "$clean/irf_aipw_flow_ch_`ch'_`g'.dta", replace
        restore
    }
}
di as result "IRF datasets saved: irf_aipw_flow_ch_<channel>_nd/_def.dta (12)"

local c_nd  "0 84 166"
local c_def "157 36 73"
local c_zero "150 150 150"
local title_credit "Private credit"
local title_claims_govt "Bank claims on govt"
local title_inv "Investment"
local title_govexp "Govt expenditure"
local title_pb "Primary balance"
local title_fdi "FDI"

local i = 1
local gnames
foreach ch of local channels {
    preserve
        use "$clean/irf_aipw_flow_ch_`ch'_nd.dta", clear
        append using "$clean/irf_aipw_flow_ch_`ch'_def.dta"
        keep if horizon >= 0
        twoway ///
            (rarea lo hi horizon if series=="nd",  color("`c_nd'%16")  lwidth(none)) ///
            (rarea lo hi horizon if series=="def", color("`c_def'%16") lwidth(none)) ///
            (connected b horizon if series=="nd",  lcolor("`c_nd'")  mcolor("`c_nd'")  msymbol(circle) lwidth(medthick)) ///
            (connected b horizon if series=="def", lcolor("`c_def'") mcolor("`c_def'") msymbol(square) lpattern(dash) lwidth(medthick)), ///
            yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
            xlabel(0(1)5, labsize(small)) ylabel(, format(%4.1f) labsize(small)) ///
            xtitle("Year", size(small)) ytitle("", size(small)) ///
            title("`title_`ch''", size(small) color(navy)) legend(off) ///
            graphregion(color(white)) plotregion(color(white)) name(gac`i', replace)
    restore
    local gnames `gnames' gac`i'
    local ++i
}
graph combine `gnames', cols(3) ///
    title("AIPW Transmission Channels by Resolution — Flow Specification", size(medsmall) color(navy)) ///
    subtitle("Navy circles = non-default; brick squares = default-linked. Shaded = 95% analytic-SE band.", size(small)) ///
    graphregion(color(white)) xsize(11) ysize(7)
capture graph export "$figs/fig_aipw_ch_flow_act2.pdf", replace
if _rc di as error "  ** fig_aipw_ch_flow_act2.pdf export failed (rc=" _rc ")"
else {
    capture graph export "$figs/fig_aipw_ch_flow_act2.png", replace width(1200)
    di as result "Figure saved: fig_aipw_ch_flow_act2.pdf/.png"
}
foreach nm of local gnames {
    capture graph drop `nm'
}

di as result _n "24_aipw_channels_flow.do complete."
di as result "  Compare each channel's h=0 def/nd point estimates to"
di as result "  23_channels_resolution_flow.do's OLS numbers for the same"
di as result "  channel/type/horizon: same sign, default deeper, is the expected"
di as result "  relationship (13c_aipw_channels.do's own header states this for"
di as result "  the onset-coded case; the same logic applies here)."
