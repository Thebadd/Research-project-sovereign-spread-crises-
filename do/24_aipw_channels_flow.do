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
  21_aipw_flow.do's header (Sections 1a-1f) documents an extended diagnostic
  investigation into why the flow AIPW's propensity model near-separates on
  continuation rows, and settles on a fix: Eq. (2) is fit ONLY on tranquil and
  onset rows (continuation==0), and the fitted model is then EXTRAPOLATED to
  every row (continuation included) because Eq. (1)/Eq. (3) keep the full
  flow-coded treatment. That investigation is not repeated here — this file's
  `_aipw' program is 21's, copied as it stands today, with every original
  control ($ctrl_core, including l_ca/l_debt) kept exactly as 21 keeps them.
  See 21_aipw_flow.do's header for the justification; it is referenced here,
  not re-derived.

  PREDICTOR CHANGE: past_def_onsets -> years_since_def_onset. Section 1a below
  diagnosed severe weight concentration in the def arm even after the
  continuation==0 fix, with past_def_onsets (a running COUNT that never
  resets) as the leading suspect. Section 1b tested years_since_def_onset
  (years since the country's most recent PRIOR default-linked onset) as a
  replacement: individually significant, sensible sign, but it did NOT fix the
  weight concentration (98.9% -> 98.6%, essentially unchanged) -- adopted for
  being a better, independently-motivated predictor, not as a variance fix.
  Section 2's estimation loop uses years_since_def_onset (via cz_recency,
  built alongside the historical cz_def so Sections 1a/1b's own already-run
  diagnostics stay exactly as documented and reproducible).

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

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow ep_seq years_since_def_onset {
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

* EXPLORATORY: set to 1 to test the alternate flow control set (see
* 18_transforms.do's "EXPLORATORY ALTERNATE FLOW CONTROL SET"). Default 0 =
* current baseline. This drives every control-set reference below, INCLUDING
* Sections 1a/1b's diagnostics -- with it on, their printed comparison
* numbers describe the run that produced years_since_def_onset, not this
* run, so read them as history rather than a live check when testing.
local use_flowalt_ctrl 0
if `use_flowalt_ctrl' local ctrl_flow_base $ctrl_flow_flowalt
else                  local ctrl_flow_base $ctrl_flow

local ctrl_credit      : list ctrl_flow_base - epc_lc
local ctrl_credit      `ctrl_credit' epc_pre_credit
local ctrl_claims_govt `ctrl_flow_base' epc_pre_claims_govt
local ctrl_inv         `ctrl_flow_base' epc_pre_inv
local ctrl_govexp      : list ctrl_flow_base - epc_lg
local ctrl_govexp      `ctrl_govexp' epc_pre_govexp
local ctrl_pb          `ctrl_flow_base' epc_pre_pb
local ctrl_fdi         `ctrl_flow_base' epc_pre_fdi

* ── Propensity model — SAME as 21_aipw_flow.do, unchanged ───────────────────
local cx = cond(`use_flowalt_ctrl', "$ctrl_core_flowalt", "$ctrl_core")   // Eq. (2): row-dated, fit on continuation==0 only
* cz_def: the ORIGINAL predictor set, kept so Sections 1a/1b's already-run
* diagnostics stay exactly as documented. Not used by Section 2's estimation.
local cz_def l_fedfunds l_reg_crisis_share past_def_onsets
* cz_recency: the ACTIVE predictor set Section 2 uses -- see header
* "PREDICTOR CHANGE" and 21_aipw_flow.do's header for the full note.
local cz_recency l_fedfunds l_reg_crisis_share years_since_def_onset

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

* ══════════════════════════════════════════════════════════════════════════
* 1a. WHY IS THE def-ARM SE SO LARGE? PROPENSITY/WEIGHT DIAGNOSTIC (credit, h=1)
*
* Section 2 below is expected to show def-arm point estimates comparable in
* size to the nd arm (or larger), but SEs several times bigger, wide enough
* that def is often not individually significant despite large coefficients.
* Two distinct explanations are possible: (a) a handful of continuation rows
* have EXTRAPOLATED propensities still close to the [0.01,0.99] trim
* boundary, so their IPW weight (1/p or 1/(1-p)) is huge and dominates the
* variance of the AIPW summand; or (b) the def sample is just thin and
* genuinely high-variance (14 countries, 10 of them also contributing nd
* rows), and the wide SE is an honest reflection of that, not an estimator
* artifact. This section checks (a) directly, using credit at h=1 as the
* representative case (the channel/horizon with the most extreme def SE),
* rather than assuming either explanation and picking a fix on that basis.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "1a. def-ARM PROPENSITY/WEIGHT DIAGNOSTIC (credit, h=1)"
di as result "════════════════════════════════════════════════════════════"

capture drop _psdiag _wtdiag
quietly probit in_crisis_def `cx' `cz_def' if sample_flow==1 & in_crisis_nd==0 & continuation==0
quietly predict double _psdiag if sample_flow==1 & in_crisis_nd==0, pr
quietly replace _psdiag = .01 if _psdiag < .01                  & !missing(_psdiag)
quietly replace _psdiag = .99 if _psdiag > .99 & !missing(_psdiag)
* Same touse markout _aipw would apply for this specific (y, omodel, pmodel) call.
* missing() takes comma-separated arguments, not a varlist, so each control is
* checked individually in the loop below rather than passed to missing() as a
* group.
quietly gen byte _touse_diag = !missing(ch_credit_0, in_crisis_def, _psdiag)
foreach v of local ctrl_credit {
    quietly replace _touse_diag = 0 if missing(`v')
}
quietly gen double _wtdiag = in_crisis_def/_psdiag + (1-in_crisis_def)/(1-_psdiag) if _touse_diag==1

quietly count if in_crisis_def==1 & _touse_diag==1
local ntr = r(N)
quietly count if in_crisis_def==1 & _touse_diag==1 & onset_all==1
local nonset = r(N)
quietly count if in_crisis_def==1 & _touse_diag==1 & continuation==1
local ncont = r(N)
di as result "  treated (def) rows in this cell: " %4.0f `ntr' "  (onset " %3.0f `nonset' ", continuation " %3.0f `ncont' ")"

di as result _n "  Extrapolated p(treated=def), among def==1 rows:"
quietly summarize _psdiag if in_crisis_def==1 & _touse_diag==1, detail
di as result "    min=" %5.3f r(min) "  p10=" %5.3f r(p10) "  median=" %5.3f r(p50) ///
             "  p90=" %5.3f r(p90) "  max=" %5.3f r(max)
quietly count if in_crisis_def==1 & _touse_diag==1 & (_psdiag<.05 | _psdiag>.95)
di as result "    rows with p<0.05 or p>0.95: " %3.0f r(N) " / " %3.0f `ntr'
quietly count if in_crisis_def==1 & _touse_diag==1 & (_psdiag<=.011 | _psdiag>=.989)
di as result "    rows AT the trim boundary (<=.011 or >=.989): " %3.0f r(N) " / " %3.0f `ntr'

di as result _n "  IPW weight (1/p or 1/(1-p)), among def==1 rows:"
quietly summarize _wtdiag if in_crisis_def==1 & _touse_diag==1, detail
di as result "    min=" %6.2f r(min) "  median=" %6.2f r(p50) "  p90=" %6.2f r(p90) "  max=" %6.2f r(max)
di as result "    mean=" %6.2f r(mean) "  sd=" %6.2f r(sd)

* Isolate whether a SMALL NUMBER of extreme-weight rows dominate the variance
* of the summand -- the direct mechanism, not just a proxy via the weight.
* Done in a preserve/keep block so gsort never has to contend with the
* thousands of rows outside this estimation cell (which would have missing
* _sqdev_diag and sort to the top of a descending gsort, corrupting the
* "top 5%" selection below).
quietly gen double _summand_diag = ///
    in_crisis_def*ch_credit_0/_psdiag - (1-in_crisis_def)*ch_credit_0/(1-_psdiag) if _touse_diag==1

preserve
    quietly keep if _touse_diag==1
    quietly summarize _summand_diag, meanonly
    local summean = r(mean)
    quietly gen double _sqdev_diag = (_summand_diag - `summean')^2
    quietly summarize _sqdev_diag, meanonly
    local alltopvar = r(sum)
    local ntot = r(N)
    local ntop = max(1, round(0.05*`ntot'))
    gsort -_sqdev_diag
    quietly summarize _sqdev_diag in 1/`ntop', meanonly
    local topvar = r(sum)
restore

di as result _n "  Concentration check: top 5% of rows (" %3.0f `ntop' " of " %4.0f `ntot' ///
             ") by squared deviation from the summand's mean account for " ///
             %5.1f 100*`topvar'/`alltopvar' " pct of its total variance."

di as result _n "  READ THIS AS: if a small handful of rows sit at the trim boundary AND"
di as result "  carry a large share of the summand's variance, the SE is inflated by a"
di as result "  few extreme extrapolated weights (fixable: tighter trim, or overlap"
di as result "  weights). If p(treated) is spread out with none near the boundary, and"
di as result "  the variance is NOT concentrated in a handful of rows, the wide SE is"
di as result "  an honest reflection of a thin, heterogeneous def sample -- not an"
di as result "  estimator artifact to engineer away."

capture drop _psdiag _wtdiag _touse_diag _summand_diag
sort cid year   // belt-and-suspenders before Section 1b / Section 2

* ══════════════════════════════════════════════════════════════════════════
* 1b. ALTERNATIVE PREDICTOR — years_since_def_onset IN PLACE OF past_def_onsets
*
* past_def_onsets (a running COUNT of a country's own default-linked onsets)
* is 1a's prime suspect: it never resets, so for a serial defaulter it behaves
* close to a permanent country identifier rather than a genuine time-varying
* predictor -- exactly the profile that produces the weight concentration 1a
* found (top 5% of rows = 98.9% of the summand's variance).
*
* years_since_def_onset measures RECENCY instead of FREQUENCY: years since
* this country's most recent PRIOR default-linked onset (not the current
* episode's own onset -- see construction below). Unlike a count, it does not
* only grow -- it keeps distinguishing "just had a default" from "had one
* long ago" for the same country, so it should not collapse a serial
* defaulter's whole history into one high, unchanging value.
*
* CONSTRUCTION, WHY IT IS NOT CIRCULAR: _aipw fits Eq. (2) ONLY on tranquil
* and onset rows (continuation==0) -- see 21_aipw_flow.do's Section 1f. For an
* ONSET row, "years since the last onset" necessarily refers to an EARLIER,
* DISTINCT episode (an onset row cannot be its own prior onset), so there is
* no way for this to circularly measure the very episode being predicted. For
* a TRANQUIL row it is unambiguous prior history. Continuation rows are never
* part of the fitting step at all -- they are only scored afterward by
* extrapolation, where this variable reads as "years into the current
* episode," a different but well-defined quantity, not a problem for fitting.
* No entry-dating trick (unlike epc_X) is needed for exactly this reason.
*
* CENSORING: countries with no PRIOR default-linked onset have no defined
* "years since." Censored to 50 (safely beyond this panel's ~35-year span,
* i.e. reads as "effectively never"), following past_def_onsets' own
* convention of replacing missing with a fixed value rather than dropping
* rows.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "1b. ALTERNATIVE PREDICTOR: years_since_def_onset (credit, h=1)"
di as result "════════════════════════════════════════════════════════════"

capture drop _defyear _defyear_lag years_since_def
* Carry forward the most recent default-onset YEAR seen so far (including the
* current row if it is itself an onset) -- same bysort(cid year) logic
* past_def_onsets' own cum_def uses, just tracking a year instead of a count.
gen _defyear = year if onset_def==1
bysort cid (year): replace _defyear = _defyear[_n-1] if missing(_defyear) & _n>1
* Lag by one row/year, exactly as past_def_onsets = L.cum_def does, so the
* current row's OWN onset status (if any) is excluded -- this is what makes
* an onset row's value refer to an earlier, distinct episode, not itself.
bysort cid (year): gen _defyear_lag = _defyear[_n-1]
gen double years_since_def = year - _defyear_lag if !missing(_defyear_lag)
replace years_since_def = 50 if missing(years_since_def)
label var years_since_def "Years since country's most recent PRIOR default-linked onset (censored at 50)"
drop _defyear _defyear_lag
sort cid year

local cz_altrecency l_fedfunds l_reg_crisis_share years_since_def

di as result "  (i) def-arm probit coefficients, years_since_def_onset in place of past_def_onsets:"
probit in_crisis_def `cx' `cz_altrecency' if sample_flow==1 & in_crisis_nd==0 & continuation==0

di as result _n "  (ii) Same weight/concentration check as 1a, this predictor in place of past_def_onsets:"
capture drop _psdiag2 _wtdiag2
quietly probit in_crisis_def `cx' `cz_altrecency' if sample_flow==1 & in_crisis_nd==0 & continuation==0
quietly predict double _psdiag2 if sample_flow==1 & in_crisis_nd==0, pr
quietly replace _psdiag2 = .01 if _psdiag2 < .01                  & !missing(_psdiag2)
quietly replace _psdiag2 = .99 if _psdiag2 > .99 & !missing(_psdiag2)
quietly gen byte _touse_diag2 = !missing(ch_credit_0, in_crisis_def, _psdiag2)
foreach v of local ctrl_credit {
    quietly replace _touse_diag2 = 0 if missing(`v')
}
quietly gen double _wtdiag2 = in_crisis_def/_psdiag2 + (1-in_crisis_def)/(1-_psdiag2) if _touse_diag2==1

quietly summarize _psdiag2 if in_crisis_def==1 & _touse_diag2==1, detail
di as result "    p(treated=def): min=" %5.3f r(min) "  median=" %5.3f r(p50) "  max=" %5.3f r(max)
quietly count if in_crisis_def==1 & _touse_diag2==1 & (_psdiag2<.05 | _psdiag2>.95)
local ntr2 = r(N)
quietly count if in_crisis_def==1 & _touse_diag2==1
di as result "    rows with p<0.05 or p>0.95: " %3.0f `ntr2' " / " %3.0f r(N)

quietly summarize _wtdiag2 if in_crisis_def==1 & _touse_diag2==1, detail
di as result "    IPW weight: median=" %6.2f r(p50) "  p90=" %6.2f r(p90) "  max=" %6.2f r(max)

quietly gen double _summand_diag2 = ///
    in_crisis_def*ch_credit_0/_psdiag2 - (1-in_crisis_def)*ch_credit_0/(1-_psdiag2) if _touse_diag2==1
preserve
    quietly keep if _touse_diag2==1
    quietly summarize _summand_diag2, meanonly
    local summean2 = r(mean)
    quietly gen double _sqdev_diag2 = (_summand_diag2 - `summean2')^2
    quietly summarize _sqdev_diag2, meanonly
    local alltopvar2 = r(sum)
    local ntot2 = r(N)
    local ntop2 = max(1, round(0.05*`ntot2'))
    gsort -_sqdev_diag2
    quietly summarize _sqdev_diag2 in 1/`ntop2', meanonly
    local topvar2 = r(sum)
restore
di as result "    concentration: top 5% of rows (" %3.0f `ntop2' " of " %4.0f `ntot2' ///
             ") account for " %5.1f 100*`topvar2'/`alltopvar2' " pct of the summand's variance."

di as result _n "  Compare directly to 1a (past_def_onsets): p range 0.010-0.938, 15/65 rows"
di as result "  with p<0.05 or p>0.95, top 5% of rows = 98.9 pct of variance, max weight 100.00."
di as result "  A materially narrower p range, fewer extreme rows, and a lower concentration"
di as result "  share here would say years_since_def_onset is both genuinely predictive AND"
di as result "  better-behaved -- worth carrying into cz_def as a real replacement, not just"
di as result "  a weaker stand-in chosen to shrink the numbers."

capture drop _psdiag2 _wtdiag2 _touse_diag2 _summand_diag2
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
            omod(`ctrl_`ch'') pz(`cx' `cz_recency') reps(`nboot') boot(row)

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
