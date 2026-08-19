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
  construction and control sets from 11/11b/12/13:
    credit, claims_govt, inv, govexp, pb, fdi           (11 / 12)
    claimsgov_assets, claimpriv_assets                  (11b nexus)
    ca (current account)                                (13 Aguiar-Gopinath)

  PROPENSITY model = identical to 08b (selection into treatment is the same
  object regardless of outcome). Only the OUTCOME regression is channel-specific.

  CRITICAL: bsample (cluster bootstrap) destroys the panel time order, so nothing
  re-estimated inside the bootstrap may use L./F. operators. The outcome ch_*_h is
  precomputed, and every lagged control is pre-generated as a PLAIN column
  (l_credit = L.credit, ...) before estimation. cx/cz are already plain columns.

  Output: $tabs/aipw_channels.csv ; $figs/fig_aipw_ch_act1.pdf ;
          $figs/fig_aipw_ch_act2.pdf. Leaves 11/11b/12/13 (OLS+IPW) untouched.
  Runtime note: heavy (9 channels x ~15 fits x nboot). nboot=300 for a practical
  run; raise to 500 for the final.  Run AFTER 17_predictors.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
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
local cz_def l_fedfunds l_reg_crisis_share past_def_onsets   // Act 2 predictors Z2

* AIPW outcome-model core = the common core, unchanged. The depth term is
* l_credit_bank (WDI FD.AST.PRVT.GD.ZS, credit by banks).
* This block used to swap it for l_credit (FS.AST.PRVT.GD.ZS, all financial
* corporations) on the belief that the by-banks series had the coverage hole.
* 19_sample_audit.do showed the opposite: l_credit_bank is non-missing on EVERY
* row where l_credit is (overlap 1180 = all of l_credit) plus 141 more, and the
* two correlate 0.950. The swap was costing observations, not saving them, so
* the core now carries l_credit_bank everywhere and this file needs no exception.
* Both are plain saved columns, so the bootstrap stays operator-free either way.
local core_aipw l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl

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
    return scalar theta = r(mean)
    return scalar N     = r(N)
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
* PROGRAM — paired bootstrap of the DIFFERENCE between two AIPW cells
*   (cell1 - cell2) recomputed on the SAME cluster resample each draw. Here
*   cell1 = default vs tranquil, cell2 = non-default vs tranquil, so r(dh) =
*   def - nd = the extra channel response of default. Matches the paper's
*   inference for the between-type contrast (they bootstrap the gap directly
*   rather than eyeballing the distance between two separately-banded lines).
*   _aipwpairdiff, y() d1() if1() d2() if2() omod() pz() reps()
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipwpairdiff
program define _aipwpairdiff, rclass
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
    * 4-level stratum (has-d1 x has-d2) so one resample preserves the treated
    * count of BOTH cells; otherwise the thinner cell is routinely lost and the
    * paired difference draw goes with it.
    _mkstrat `d1' `d2', generate(_strat)
    tempname pf
    tempfile bf
    quietly postfile `pf' double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            capture drop _bid
            bsample, cluster(cid) strata(_strat) idcluster(_bid)
            capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(_bid)
            local t1 = cond(_rc==0, r(theta), .)
            capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(_bid)
            local t2 = cond(_rc==0, r(theta), .)
            if !missing(`t1') & !missing(`t2') quietly post `pf' (`t1' - `t2')
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

* second results file: the def - nd DIFFERENCE per channel x horizon (paired boot)
tempname Rd
tempfile diffresf
postfile `Rd' str24 channel byte horizon double dhl bdef bnd se lo hi nd ///
    using "`diffresf'", replace

foreach ch in credit claims_govt inv govexp pb fdi ///
              claimsgov_assets claimpriv_assets ca {

    * channel-specific OUTCOME-model controls (pre-lagged plain columns)
    * AIPW outcome core ($core_aipw = the common core, depth term l_credit_bank) +
    * channel's own pre_<v>; drop the core term equal to the channel's own lagged
    * level (credit->the depth term, govexp->l_govexp, ca->ca).
    if      "`ch'" == "credit"            local om l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_hyperinfl pre_credit
    else if "`ch'" == "govexp"            local om l1_gdpg l_debt l_ca l_banking_duration l_open l_credit_bank l_hyperinfl pre_govexp
    else if "`ch'" == "ca"                local om l1_gdpg l_debt l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl pre_ca
    else                                  local om `core_aipw' pre_`ch'

    di as result _n "=== CHANNEL: `ch' ==="

    * ── Act 1: all-crises AIPW (onset_all vs tranquil) ──────────────────────
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

    * ── Act 2: resolution split, two level lines vs tranquil ────────────────
    di as result "  Act 2 (vs tranquil):  type  h    ATE      SE      [95% CI]"
    foreach spec in "onset_nd onset_def nd" "onset_def onset_nd def" {
        gettoken Dv rest  : spec
        gettoken drop sn  : rest
        post `R' ("`ch'") ("`sn'") (0) (0) (0) (0) (0)   // explicit baseline (h=0)
        forvalues h = 0/4 {
            _aipwci ch_`ch'_`h' `Dv', ifc(sample==1 & `drop'==0) ///
                omod(`om') pz(`om' `cz_def') reps(`nboot')
            if r(ok) {
                local b=r(b)
                local se=r(se)
                local lo=r(lo)
                local hi=r(hi)
                local nd=r(nd)
                post `R' ("`ch'") ("`sn'") (`h'+1) (`b') (`se') (`lo') (`hi')
                di "    `sn'  h=" `h'+1 "  " %8.3f `b' "  " %6.3f `se' ///
                   "  [" %7.3f `lo' ", " %7.3f `hi' "]  " `nd' "/`nboot'"
            }
            else di as error "    `sn' h=" `h'+1 ": estimate failed (thin sample)."
        }
    }

    * ── Act 2 (difference): extra channel response of default (def - nd) ────
    *   bootstrapped directly on paired cluster resamples, so the between-type
    *   gap gets a proper CI (the paper's inference), not an eyeballed distance.
    di as result "  Act 2 (def - nd gap):  h    def-nd   SE      [95% CI]   draws"
    post `Rd' ("`ch'") (0) (0) (0) (0) (0) (0) (0) (0)   // explicit baseline (h=0)
    forvalues h = 0/4 {
        _aipwpairdiff, y(ch_`ch'_`h') ///
            d1(onset_def) if1(sample==1 & onset_nd==0) ///
            d2(onset_nd)  if2(sample==1 & onset_def==0) ///
            omod(`om') pz(`om' `cz_def') reps(`nboot')
        if r(ok) {
            post `Rd' ("`ch'") (`h'+1) (r(dh)) (r(b1)) (r(b2)) (r(se)) (r(lo)) (r(hi)) (r(nd))
            local sig = cond(r(nd)>=50 & (r(lo)>0 | r(hi)<0), " *", "")
            di "    h=" `h'+1 "  " %8.3f r(dh) "  " %6.3f r(se) ///
               "  [" %7.3f r(lo) ", " %7.3f r(hi) "]  " r(nd) "/`nboot'`sig'"
        }
        else di as error "    h=" `h'+1 ": def-nd gap failed (thin sample)."
    }
}
postclose `R'
postclose `Rd'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — CSV + two small-multiple (by-channel) figures
* ══════════════════════════════════════════════════════════════════════════
* ── Difference CSV: def - nd gap per channel x horizon (paired bootstrap) ────
preserve
    use "`diffresf'", clear
    label var dhl  "AIPW def - nd channel gap (pp)"
    label var bdef "Default-linked ATE"
    label var bnd  "Non-default ATE"
    label var lo   "95% percentile CI lower"
    label var hi   "95% percentile CI upper"
    label var nd   "Valid bootstrap draws"
    gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
    label var sig95 "Gap CI excludes 0"
    order channel horizon dhl bdef bnd se lo hi nd sig95
    export delimited "$tabs/aipw_channels_diff.csv", replace
    di as result _n "AIPW channel def-nd difference CSV saved: $tabs/aipw_channels_diff.csv"
restore

use "`resf'", clear
label var b  "AIPW ATE (pp of the channel ratio)"
label var lo "95% percentile CI lower"
label var hi "95% percentile CI upper"
order channel series horizon b se lo hi
export delimited "$tabs/aipw_channels.csv", replace
di as result _n "AIPW channel results CSV saved: $tabs/aipw_channels.csv"

encode channel, gen(chid)

* ── Figure A: Act 1 all-crises AIPW per channel (line + CI band) ────────────
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

* ── Figure B: Act 2 resolution split per channel (nd vs def, fig8 palette) ──
local c_nd  "0 84 166"
local c_def "157 36 73"
capture twoway ///
    (rarea lo hi horizon if series=="nd",  color("`c_nd'%16")  lwidth(none)) ///
    (rarea lo hi horizon if series=="def", color("`c_def'%16") lwidth(none)) ///
    (connected b horizon if series=="nd",  lcolor("`c_nd'")  lwidth(medthick) msymbol(circle)) ///
    (connected b horizon if series=="def", lcolor("`c_def'") lwidth(medthick) lpattern(dash) msymbol(square)), ///
    by(chid, yrescale ///
        note("Two AIPW level IRFs per channel; shaded = bootstrap 95% CI. Gap = extra cost of default.", size(vsmall)) ///
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
