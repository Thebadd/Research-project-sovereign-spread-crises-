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
if "$ctrl_core"=="" global ctrl_core "l1_gdpg debt ca banking_crisis l_govexp l_open l_credit hyperinf_dummy"
sort cid year
xtset cid year

local nboot  = 300      // bootstrap reps per (channel, series, horizon)
local cx     l1_gdpg l2_gdpg debt ca infl imf
local cz     fedfunds l_reg_crisis_share past_onsets       // Act 1 predictors Z1
local cz_def fedfunds l_reg_crisis_share past_def_onsets   // Act 2 predictors Z2

* AIPW outcome-model core = the common core with the thin by-banks credit-depth
* term (l_credit_bank) swapped for the well-covered TOTAL private-credit lag
* (l_credit). On the restricted channel/nexus cells l_credit_bank's coverage hole
* forces listwise deletion that collapses the AIPW reg/probit sample (only the
* credit channel, whose om already omits l_credit_bank, survived otherwise).
* l_credit is pre-lagged below, so the bootstrap stays operator-free.
local core_aipw l1_gdpg debt ca banking_crisis l_govexp l_open l_credit hyperinf_dummy

* ── Build channel outcomes ch_v_h = F h.v - L.v (h=0..4) ─────────────────────
foreach v in credit claims_govt inv govexp pb fdi ///
             claimsgov_assets claimpriv_assets ca {
    capture drop `v'_base
    gen `v'_base = L.`v'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen ch_`v'_`h' = F`h'.`v' - `v'_base
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var);
    * added to each outcome model to absorb the channel's own pre-trend momentum.
    capture drop pre_`v'
    gen pre_`v' = L.`v' - L2.`v'
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
capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'
    tempvar xb m0 m1 ps summ
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
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

    tempname pf
    tempfile bf
    quietly postfile `pf' double theta using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            capture drop _bid
            bsample, cluster(cid) idcluster(_bid)
            capture _aipw `yv' `Dv' if `ifc', omodel(`omod') pmodel(`pz') fe(_bid)
            if _rc == 0 quietly post `pf' (r(theta))
        restore
    }
    quietly postclose `pf'

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
* ESTIMATE — loop over channels; Act 1 + Act 2 (nd/def); post to results file
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str24 channel str4 series byte horizon double b se lo hi ///
    using "`resf'", replace

foreach ch in credit claims_govt inv govexp pb fdi ///
              claimsgov_assets claimpriv_assets ca {

    * channel-specific OUTCOME-model controls (pre-lagged plain columns)
    * AIPW outcome core ($core_aipw: common core with l_credit->total credit) +
    * channel's own pre_<v>; drop the core term equal to the channel's own lagged
    * level (credit->l_credit, govexp->l_govexp, ca->ca).
    if      "`ch'" == "credit"            local om l1_gdpg debt ca banking_crisis l_govexp l_open hyperinf_dummy pre_credit
    else if "`ch'" == "govexp"            local om l1_gdpg debt ca banking_crisis l_open l_credit hyperinf_dummy pre_govexp
    else if "`ch'" == "ca"                local om l1_gdpg debt banking_crisis l_govexp l_open l_credit hyperinf_dummy pre_ca
    else                                  local om `core_aipw' pre_`ch'

    di as result _n "=== CHANNEL: `ch' ==="

    * ── Act 1: all-crises AIPW (onset_all vs tranquil) ──────────────────────
    di as result "  Act 1 (all onsets):   h    ATE      SE      [95% CI]   draws"
    forvalues h = 0/4 {
        _aipwci ch_`ch'_`h' onset_all, ifc(sample==1) ///
            omod(`om') pz(`cx' `cz') reps(`nboot')
        if r(ok) {
            local b=r(b)
            local se=r(se)
            local lo=r(lo)
            local hi=r(hi)
            local nd=r(nd)
            post `R' ("`ch'") ("all") (`h') (`b') (`se') (`lo') (`hi')
            di "    h=" `h' "  " %8.3f `b' "  " %6.3f `se' ///
               "  [" %7.3f `lo' ", " %7.3f `hi' "]  " `nd' "/`nboot'"
        }
        else di as error "    h=" `h' ": Act 1 estimate failed (rc)."
    }

    * ── Act 2: resolution split, two level lines vs tranquil ────────────────
    di as result "  Act 2 (vs tranquil):  type  h    ATE      SE      [95% CI]"
    foreach spec in "onset_nd onset_def nd" "onset_def onset_nd def" {
        gettoken Dv rest  : spec
        gettoken drop sn  : rest
        forvalues h = 0/4 {
            _aipwci ch_`ch'_`h' `Dv', ifc(sample==1 & `drop'==0) ///
                omod(`om') pz(`cx' `cz_def') reps(`nboot')
            if r(ok) {
                local b=r(b)
                local se=r(se)
                local lo=r(lo)
                local hi=r(hi)
                local nd=r(nd)
                post `R' ("`ch'") ("`sn'") (`h') (`b') (`se') (`lo') (`hi')
                di "    `sn'  h=" `h' "  " %8.3f `b' "  " %6.3f `se' ///
                   "  [" %7.3f `lo' ", " %7.3f `hi' "]  " `nd' "/`nboot'"
            }
            else di as error "    `sn' h=" `h' ": estimate failed (thin sample)."
        }
    }
}
postclose `R'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — CSV + two small-multiple (by-channel) figures
* ══════════════════════════════════════════════════════════════════════════
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
    xlabel(0(1)4) xtitle("Years after onset", size(small)) ///
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
    xlabel(0(1)4) xtitle("Years after onset", size(small)) ///
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
