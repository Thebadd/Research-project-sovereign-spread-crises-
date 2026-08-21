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
  Eq. (1)  Outcome regression by PLAIN OLS with country FE, giving the
           conditional means m1 and m0. This is 20_lp_flow.do's specification,
           less the year FE.
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
    a coefficient in its own right.

  ONE DELIBERATE DIFFERENCE FROM 08b, WORTH RECORDING
  ---------------------------------------------------
  The reference paper's TEXT says Eq. (1) is estimated "by OLS" and shows it
  unweighted. Their CODE, inside the AIPW loop, weights it:
      reg g_h dum g_0 $convar [pweight=invwt]
  which makes their estimator IPWRA rather than plain AIPW. 08b_aipw.do follows
  the code. This file follows the TEXT, and estimates Eq. (1) unweighted, for
  two reasons: it is what the paper as written specifies, and it keeps Eq. (1)
  identical to the model developed in 20_lp_flow.do apart from the year FE,
  which is the whole point of this file. Both are doubly robust — each is
  consistent if EITHER the outcome model or the propensity model is correct —
  so this is a choice between two valid estimators, not a correction.

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

    tempvar xb m0 m1 ps summ
    * Eq. (1) — PLAIN OLS with country FE, as the paper's text specifies
    * ("we estimate the following regression model by OLS"). No weights: this
    * is 20_lp_flow.do's specification, less the year FE. See the header for
    * why this departs from 08b.
    if "`fe'" != "" quietly reg `y' `D' `omodel' i.`fe' if `touse'
    else            quietly reg `y' `D' `omodel'         if `touse'
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']     if `touse'

    * Eq. (2) — propensity probit, winsorised to bound the weights.
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01                & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
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
             OMOD(string) PZ(string) REPS(integer)
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
    quietly probit in_crisis_`s2' `cx' `cz_def' if sample_flow==1 & `rival'==0
    quietly predict double _ps2 if e(sample), pr
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
di as result "       se_a = analytic influence-function SE (the paper's); sd_b = bootstrap SD (diagnostic)"

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
    local A1 = r(a1)      // analytic influence-function SE, def  (the paper's)
    local A2 = r(a2)      // analytic influence-function SE, nd
    local S1 = r(s1)      // bootstrap SD, def  (diagnostic only)
    local S2 = r(s2)      // bootstrap SD, nd   (diagnostic only)
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

    local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
    di as result "  " %1.0f `hd' "  " %8.3f `B2' " (" %5.3f `A2' ")  " ///
                 %8.3f `B1' " (" %5.3f `A1' ")  " %8.3f `DH' ///
                 " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
                 " " %7.3f `zz' " " %5.3f `pz' " " %4.0f `ND'
    di as result "        bootstrap SD for comparison:  nd " %6.3f `S2' ///
                 "   def " %6.3f `S1'
}

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
di as result "  Read Section 1 before Section 3: if continuation rows do not"
di as result "  survive the propensity, this file re-ran 08b and says nothing new."
