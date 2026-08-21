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
  Eq. (1)  Outcome regression by OLS with country FE, IPW-weighted (their
           IPWRA form: `reg g_h dum g_0 $convar [pweight=invwt]'), giving the
           conditional means m1 and m0.
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
  * `reg' RATHER THAN `xtscc'. Mechanical, not elective: xtscc does not accept
    pweight, and Eq. (1) must be inverse-probability weighted. Same reason 08b
    uses reg. Inference therefore comes from the bootstrap, not from DK SEs.

  WHICH CONTROL SET GOES WHERE — AND WHY THEY DIFFER
  --------------------------------------------------
  Eq. (1) outcome model:  $ctrl_flow  (episode-dated, epc_*)
  Eq. (2) propensity:     $ctrl_core  (row-dated)

  This is not an oversight. epc_X is built as
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

  Outputs
  -------
    "$tabs/table10_aipw_flow.rtf"    the two type lines + the difference
    "$tabs/aipw_flow_diff.csv"       def-nd gap, bootstrap CI, Clogg z
    "$tabs/aipw_flow_compare.csv"    onset AIPW vs flow AIPW vs flow OLS
    "$clean/irf_aipw_flow_nd.dta"    non-default AIPW IRF
    "$clean/irf_aipw_flow_def.dta"   default-linked AIPW IRF
    "$figs/fig10_aipw_flow.pdf/.png" the two-line resolution split

  RUNTIME: the paired bootstrap refits Eq. (1)-(3) twice per draw per horizon.
  At nboot=1000 that is 10,000 probit+regression pairs. Expect several minutes.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow {
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

local cx     $ctrl_core          // Eq. (2) baseline: row-dated, see header
local com    $ctrl_flow          // Eq. (1) controls: episode-dated
local cz_def l_fedfunds l_reg_crisis_share past_def_onsets

set seed 20260819
local nboot = 1000

foreach m in b se lo hi {
    matrix Fnd_`m'   = J(5,1,.)
    matrix Fdef_`m'  = J(5,1,.)
    matrix Fdiff_`m' = J(5,1,.)
}
matrix Fdiff_z = J(5,1,.)
matrix Fdiff_p = J(5,1,.)

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
    * Eq. (2) FIRST — the outcome regression is weighted by it (IPWRA).
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01                & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'

    * Eq. (1) — IPW-weighted OLS with country FE -> conditional means.
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
    return scalar theta = r(mean)
    return scalar N     = r(N)
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
             OMOD(string) PZ(string) REPS(integer)
    capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b1 = r(theta)
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2 = r(theta)
    local dh = `b1' - `b2'

    _mkstrat `d1' `d2', generate(_strat)
    tempname pf
    tempfile bf
    quietly postfile `pf' double t1 double t2 double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            capture drop _bid
            bsample, cluster(cid) strata(_strat) idcluster(_bid)
            capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(_bid)
            local v1 = cond(_rc==0, r(theta), .)
            capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(_bid)
            local v2 = cond(_rc==0, r(theta), .)
            if !missing(`v1') & !missing(`v2') quietly post `pf' (`v1') (`v2') (`v1'-`v2')
        restore
    }
    quietly postclose `pf'
    capture drop _strat

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
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar s1 = `s1'
    return scalar s2 = `s2'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* 1. DIAGNOSTICS — PRINTED BEFORE ANY COEFFICIENT
*
* Two questions decide whether anything below is worth reading.
*
* (a) Is the episode-dated control set usable in the probit? No — and this
*     block shows why rather than asserting it. epc_X differs from X only on
*     continuation rows, so the difference is a readout of treatment status.
*
* (b) Does the propensity separate treated from control WITHOUT collapsing the
*     flow treatment back into onset coding? The [0.01,0.99] winsorisation
*     bounds the weights, but if continuation rows are pushed to scores near 1
*     they contribute almost nothing and this file is silently re-running 08b.
*     The survival counts below are the test.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "1. DIAGNOSTICS"
di as result "════════════════════════════════════════════════════════════"

quietly count if in_crisis_nd==1
di as result "  treated, non-default (flow):     " %4.0f r(N) "   (expect 113)"
quietly count if in_crisis_def==1
di as result "  treated, default-linked (flow):  " %4.0f r(N) "   (expect 121)"
quietly count if sample_flow==1 & in_crisis==0
di as result "  tranquil controls:               " %4.0f r(N)

di as result _n "  (a) Why the probit cannot use the episode-dated controls:"
di as result "      rows where epc_X differs from the row-dated X, by control —"
di as result "      these are continuation rows, so the variable encodes treatment."
foreach X of global ctrl_core {
    capture confirm variable epc_`X', exact
    if _rc continue
    quietly count if !missing(`X') & abs(epc_`X' - `X') > 1e-9
    local ndiff = r(N)
    quietly count if !missing(`X') & abs(epc_`X' - `X') > 1e-9 & continuation==1
    di as result "        " %-22s "`X'" %5.0f `ndiff' " rows differ, " ///
                 %5.0f r(N) " of them continuation"
}

di as result _n "  (b) Propensity overlap and continuation survival:"
foreach s in nd def {
    quietly probit in_crisis_`s' `cx' `cz_def' if sample_flow==1
    capture drop _ps_`s'
    quietly predict double _ps_`s' if e(sample), pr
    quietly count if in_crisis_`s'==1 & !missing(_ps_`s')
    local nin = r(N)
    quietly count if in_crisis_`s'==1 & !missing(_ps_`s') & onset_all==1
    local non = r(N)
    quietly count if in_crisis_`s'==1 & !missing(_ps_`s') & continuation==1
    local ncn = r(N)
    quietly summarize _ps_`s' if in_crisis_`s'==1, detail
    local mt = r(mean)
    quietly summarize _ps_`s' if in_crisis_`s'==0 & sample_flow==1, detail
    local mc = r(mean)
    quietly count if in_crisis_`s'==1 & _ps_`s' > .99 & !missing(_ps_`s')
    local nwin = r(N)
    di as result "      `s': treated rows in probit sample = " %4.0f `nin' ///
                 "  (onset " %3.0f `non' ", continuation " %3.0f `ncn' ")"
    di as result "          mean p(treated) = " %5.3f `mt' ///
                 "   mean p(control) = " %5.3f `mc' ///
                 "   treated rows with p>0.99 = " %3.0f `nwin'
    if `ncn' == 0 {
        di as error "      ** `s': NO continuation rows survive — this has collapsed"
        di as error "         into onset coding and the check below is vacuous."
    }
}
di as result _n "      A large share of treated rows at p>0.99, or continuation rows"
di as result "      absent, means the propensity is reading 'already in crisis'"
di as result "      and the flow AIPW is not adding information over 08b."

* ══════════════════════════════════════════════════════════════════════════
* 2. IDENTITY ANCHOR — must reproduce 08b on the onset sample
*
* Restricted to sample==1 the continuation rows are gone, so in_crisis_nd IS
* onset_nd, in_crisis_def IS onset_def, and epc_* collapses to row-dated. The
* two estimators must therefore agree exactly. Same logic that validated
* 20_lp_flow.do; it catches a mis-specified call before anything is read.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "2. IDENTITY ANCHOR — flow AIPW on sample==1 must equal 08b"
di as result "════════════════════════════════════════════════════════════"
quietly _aipw dy_0 in_crisis_def if sample==1 & in_crisis_nd==0, ///
    omodel(`com') pmodel(`cx' `cz_def') fe(cid)
local anch_flow = r(theta)
quietly _aipw dy_0 onset_def if sample==1 & onset_nd==0, ///
    omodel($ctrl_core) pmodel(`cx' `cz_def') fe(cid)
local anch_ons = r(theta)
local anch_gap = abs(`anch_flow' - `anch_ons')
di as result "  flow  (def, Year 1, sample==1): " %10.6f `anch_flow'
di as result "  onset (def, Year 1, sample==1): " %10.6f `anch_ons'
if `anch_gap' > 1e-6 {
    di as error "  ** ANCHOR FAILED (difference " %10.6f `anch_gap' ")."
    di as error "     The flow call is not the onset call on this sample. Stop and fix."
}
else di as result "  MATCH (difference " %9.2e `anch_gap' ") — the transposition is exact."

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
di as result "Year   ND        DEF       def-nd    [95% boot CI]        Clogg z     p     draws"

forvalues h = 0/4 {
    local row = `h' + 1
    local hd  = `h' + 1

    _aipwpairflow, y(dy_`h') ///
        d1(in_crisis_def) if1(sample_flow==1 & in_crisis_nd==0) ///
        d2(in_crisis_nd)  if2(sample_flow==1 & in_crisis_def==0) ///
        omod(`com') pz(`cx' `cz_def') reps(`nboot')

    if !r(ok) {
        di as error "  Year `hd': estimate failed (cell too thin)."
        continue
    }

    * Pull every returned scalar into a local FIRST. r() survives matrix
    * assignment, but not every command, and reading it late is how the
    * earlier files in this project acquired silent bugs.
    local B1 = r(b1)
    local B2 = r(b2)
    local S1 = r(s1)
    local S2 = r(s2)
    local DH = r(dh)
    local SE = r(se)
    local LO = r(lo)
    local HI = r(hi)
    local ND = r(nd)

    matrix Fdef_b[`row',1]  = `B1'
    matrix Fnd_b[`row',1]   = `B2'
    matrix Fdef_se[`row',1] = `S1'
    matrix Fnd_se[`row',1]  = `S2'
    matrix Fdiff_b[`row',1]  = `DH'
    matrix Fdiff_se[`row',1] = `SE'
    matrix Fdiff_lo[`row',1] = `LO'
    matrix Fdiff_hi[`row',1] = `HI'

    * +/-1.96 bootstrap SD around each level line (the difference gets a proper
    * percentile CI above; these bands are for the figure only).
    matrix Fnd_lo[`row',1]  = `B2' - 1.96*`S2'
    matrix Fnd_hi[`row',1]  = `B2' + 1.96*`S2'
    matrix Fdef_lo[`row',1] = `B1' - 1.96*`S1'
    matrix Fdef_hi[`row',1] = `B1' + 1.96*`S1'

    * Clogg et al. (1995) z — the PERMISSIVE statistic, reported for
    * comparability with the reference paper, which gives both. It treats the
    * two cells as independent; they are not, because they share the tranquil
    * control pool. The bootstrap percentile CI above is the conservative and
    * preferred number, and where the two disagree the write-up says so.
    local zz = .
    local pz = .
    if !missing(`S1') & !missing(`S2') & (`S1'^2 + `S2'^2) > 0 {
        local zz = `DH' / sqrt(`S1'^2 + `S2'^2)
        local pz = 2*(1 - normal(abs(`zz')))
        matrix Fdiff_z[`row',1] = `zz'
        matrix Fdiff_p[`row',1] = `pz'
    }

    local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
    di as result "  " %1.0f `hd' "   " %8.3f `B2' "  " %8.3f `B1' "  " %8.3f `DH' ///
                 "  [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
                 "  " %7.3f `zz' "  " %5.3f `pz' "  " %4.0f `ND'
}

di as result _n "  * = bootstrap 95% percentile CI for the gap excludes zero."
di as result "  Clogg z assumes the two cells are independent (they share the"
di as result "  tranquil control pool), so it is the permissive statistic. Where"
di as result "  it and the bootstrap disagree, the bootstrap governs."

* ══════════════════════════════════════════════════════════════════════════
* 4. THE COMPARISON — what actually answers the question
*
* Three estimates of the same object, laid side by side:
*   (i)   onset AIPW      08b_aipw.do        -> $tabs/aipw_act2_diff.csv
*   (ii)  flow OLS        20_lp_flow.do      -> $tabs/flow_lp.csv
*   (iii) flow AIPW       this file
*
* (i) vs (iii) varies the treatment definition holding the estimator fixed.
* (ii) vs (iii) varies the estimator holding the treatment definition fixed.
* Reading both together separates the two, which neither file could do alone.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "4. COMPARISON — def-nd gap under three designs"
di as result "════════════════════════════════════════════════════════════"

tempname CMP
tempfile cmpf
postfile `CMP' int horizon double onset_aipw double flow_ols double flow_aipw ///
    double flow_aipw_lo double flow_aipw_hi using "`cmpf'", replace

* (i) onset AIPW, written by 08b
matrix ONS = J(5,1,.)
capture confirm file "$tabs/aipw_act2_diff.csv"
if _rc {
    di as error "  ** aipw_act2_diff.csv not found — run 08b_aipw.do for the onset column."
}
else {
    preserve
        quietly import delimited "$tabs/aipw_act2_diff.csv", clear varnames(1)
        forvalues h = 1/5 {
            quietly summarize dhl if horizon==`h', meanonly
            if r(N) > 0 matrix ONS[`h',1] = r(mean)
        }
    restore
}

* (ii) flow OLS, written by 20
matrix FOLS = J(5,1,.)
capture confirm file "$tabs/flow_lp.csv"
if _rc {
    di as error "  ** flow_lp.csv not found — run 20_lp_flow.do for the OLS column."
}
else {
    preserve
        quietly import delimited "$tabs/flow_lp.csv", clear varnames(1)
        forvalues h = 1/5 {
            quietly summarize b if spec=="split" & term=="def_minus_nd" & hdisp==`h', meanonly
            if r(N) > 0 matrix FOLS[`h',1] = r(mean)
        }
    restore
}

di as result "Year   onset AIPW   flow OLS   flow AIPW   [95% boot CI]"
forvalues h = 1/5 {
    di as result "  " %1.0f `h' "   " %10.3f ONS[`h',1] "  " %9.3f FOLS[`h',1] ///
                 "  " %9.3f Fdiff_b[`h',1] ///
                 "   [" %7.3f Fdiff_lo[`h',1] ", " %7.3f Fdiff_hi[`h',1] "]"
    post `CMP' (`h') (ONS[`h',1]) (FOLS[`h',1]) (Fdiff_b[`h',1]) ///
        (Fdiff_lo[`h',1]) (Fdiff_hi[`h',1])
}
postclose `CMP'

di as result _n "  Agreement in SIGN and rough magnitude across the three columns"
di as result "  means the resolution gap is not an artifact of either the"
di as result "  treatment definition or the estimator. Divergence localises which."

* ══════════════════════════════════════════════════════════════════════════
* 5. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`cmpf'", clear
    label var onset_aipw   "def-nd, onset AIPW (08b)"
    label var flow_ols     "def-nd, flow OLS (20)"
    label var flow_aipw    "def-nd, flow AIPW (this file)"
    label var flow_aipw_lo "flow AIPW 95% percentile CI lower"
    label var flow_aipw_hi "flow AIPW 95% percentile CI upper"
    export delimited "$tabs/aipw_flow_compare.csv", replace
    di as result "Comparison saved: $tabs/aipw_flow_compare.csv"
restore

preserve
    clear
    tempname pfd
    tempfile diff2
    postfile `pfd' int horizon double dhl double se double lo double hi ///
        double cloggz double cloggp using "`diff2'", replace
    post `pfd' (0) (0) (0) (0) (0) (.) (.)
    forvalues h = 1/5 {
        post `pfd' (`h') (Fdiff_b[`h',1]) (Fdiff_se[`h',1]) ///
            (Fdiff_lo[`h',1]) (Fdiff_hi[`h',1]) (Fdiff_z[`h',1]) (Fdiff_p[`h',1])
    }
    postclose `pfd'
    use "`diff2'", clear
    label var dhl     "Flow AIPW extra cost of default (def - nd, pp)"
    label var lo      "95% percentile CI lower (bootstrap)"
    label var hi      "95% percentile CI upper (bootstrap)"
    label var cloggz  "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp  "p-value of the Clogg z"
    gen byte sig95 = (!missing(lo) & (lo>0 | hi<0))
    label var sig95 "Bootstrap CI excludes 0"
    order horizon dhl se lo hi sig95 cloggz cloggp
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
* 6. FIGURE — drawn here, as 20_lp_flow.do does, so the file is self-contained
* ══════════════════════════════════════════════════════════════════════════
local c_nd   "0 84 166"
local c_def  "157 36 73"
local c_zero "150 150 150"

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
        title("Flow AIPW: Output While In a Spread Crisis, by Resolution", size(medium)) ///
        subtitle("Asonuma et al. Eqs. (1)-(3) on the flow treatment; each type vs tranquil", size(small)) ///
        note("Treatment = 1 in every year of an episode. Eq. (1) IPW-weighted OLS, country FE," ///
             "episode-dated controls; Eq. (2) pooled probit on row-dated controls + predictors;" ///
             "Eq. (3) the paper's AIPW form. Bands are +/-1.96 bootstrap SD from paired cluster" ///
             "resamples. Not the same object as Figure 9: this corrects for selection into crisis.", ///
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
di as result "  Read Section 1 before Section 3: if continuation rows do not"
di as result "  survive the propensity, this file re-ran 08b and says nothing new."
