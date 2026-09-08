/*===========================================================================
  08D_AIPW_REENTRY.DO
  Act 2 AIPW robustness under the RE-ENTRY episode definition
  (18c_reentry_variant.do): combines 08b_aipw.do's own Act 2 GDP loop AND
  13c_aipw_channels.do's own Act 2 channel loop into ONE file, same as
  08c_aipw_mindur.do's own construction -- but here the ONLY change from
  the baseline is the ESTIMATION SAMPLE (`sample==1' -> `sample_reentry==1'
  inside every _aipwpair `if1'/`if2' condition); the treatment variables
  (d1(onset_def), d2(onset_nd)) are UNCHANGED, matching
  12c_channels_resolution_reentry.do's own isolate-the-sample-effect
  design (see that file's and 18c_reentry_variant.do's headers for the
  full argument for why re-entry's real effect on the onset-tier design
  runs through the control-pool composition `sample' selects, not through
  the treatment coding).

  LEVEL SE: analytic (unclustered influence-function) SE, matching this
  session's earlier choice for the mindur AIPW build (only one SE
  convention built per robustness variant, not both).

  The _aipw/_aipwpair programs, propensity predictors (cx/cz_def), and
  outcome-model control sets are copied verbatim from 08b_aipw.do/
  13c_aipw_channels.do; nothing about the ESTIMATOR itself changes here.

  STANDALONE: reads $clean/panel_lp_reentry.dta (built by
  18c_reentry_variant.do). Does not modify 08b_aipw.do, 13c_aipw_channels.do,
  18c_reentry_variant.do, or any baseline column.

  Output: $tabs/aipw_results_reentry.csv (GDP),
          $tabs/aipw_channels_reentry.csv (channels),
          $figs/fig_aipw_combined_reentry.pdf (Panel A-F).
  Run AFTER 18c_reentry_variant.do. Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp_reentry.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

set seed 20260819
local nboot = 300

local cx     $ctrl_core
local cz_def l_fedfunds l_contagion_dist_def years_since_def_onset
local core_aipw l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2

* ── Channel outcomes ch_v_h = F h.v - L.v (h=0..4), same construction as
* 13c_aipw_channels.do -- log real level for credit/inv, ratio for the rest.
foreach v in credit claims_govt inv fdi real_lending {
    local src `v'
    if inlist("`v'","credit","inv") local src ln_r_`v'
    capture drop `v'_base
    gen `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    capture drop pre_`v'
    gen pre_`v' = L.`src' - L2.`src'
}

* ══════════════════════════════════════════════════════════════════════════
* PROGRAMS — identical to 08b_aipw.do / 13c_aipw_channels.do
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'
    tempvar xb m0 m1 ps summ iwt
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
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'
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
    return scalar theta  = `th'
    return scalar N      = `nn'
    return scalar se     = `sean'
end

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
    quietly gen byte _pool = 0 if (`if1') | (`if2')
    quietly replace _pool = 1 if `d1' == 1
    quietly replace _pool = 2 if `d2' == 1

    tempname pf
    tempfile bf
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
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* ACT 2 — GDP (sample swapped, treatment UNCHANGED; analytic SE)
* ══════════════════════════════════════════════════════════════════════════
tempname Rg
tempfile resfg
postfile `Rg' str3 series byte horizon double b se lo hi using "`resfg'", replace

di as result _n "=== ACT 2 (RE-ENTRY ROBUSTNESS: SAMPLE ONLY) — AIPW, GDP ==="
foreach s in nd def {
    post `Rg' ("`s'") (0) (0) (0) (0) (0)
}
forvalues h = 0/4 {
    _aipwpair, y(dy_`h') ///
        d1(onset_def) if1(sample_reentry==1 & onset_nd==0 & common_abcd==1) ///
        d2(onset_nd)  if2(sample_reentry==1 & onset_def==0 & common_abcd==1) ///
        omod(`core_aipw') pz(`cx' `cz_def') reps(`nboot')
    if r(ok) {
        local B1 = r(b1)
        local B2 = r(b2)
        local A1 = r(a1)
        local A2 = r(a2)
        post `Rg' ("nd")  (`h'+1) (`B2') (`A2') (`B2'-1.96*`A2') (`B2'+1.96*`A2')
        post `Rg' ("def") (`h'+1) (`B1') (`A1') (`B1'-1.96*`A1') (`B1'+1.96*`A1')
        di "  h=" `h'+1 "  nd=" %8.3f `B2' " (" %5.3f `A2' ")   def=" %8.3f `B1' " (" %5.3f `A1' ")"
    }
    else di as error "  h=" `h'+1 ": GDP estimate failed (cell too thin)."
}
postclose `Rg'

preserve
    use "`resfg'", clear
    label var b  "AIPW ATE (pp)"
    label var se "Analytic (unclustered influence-function) SE, matching the paper's own formula"
    label var lo "95% CI lower = b - 1.96*se (analytic)"
    label var hi "95% CI upper = b + 1.96*se (analytic)"
    order series horizon b se lo hi
    export delimited "$tabs/aipw_results_reentry.csv", replace
    di as result "Saved: $tabs/aipw_results_reentry.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* ACT 2 — CHANNELS (sample swapped, treatment UNCHANGED; analytic SE)
* ══════════════════════════════════════════════════════════════════════════
tempname Rc
tempfile resfc
postfile `Rc' str24 channel str4 series byte horizon double b se lo hi using "`resfc'", replace

foreach ch in credit claims_govt inv fdi real_lending {
    if      "`ch'" == "credit" local om l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
    else                       local om `core_aipw' pre_`ch'

    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

    di as result _n "=== CHANNEL (RE-ENTRY ROBUSTNESS: SAMPLE ONLY): `ch' ==="
    post `Rc' ("`ch'") ("nd")  (0) (0) (0) (0) (0)
    post `Rc' ("`ch'") ("def") (0) (0) (0) (0) (0)
    forvalues h = 0/4 {
        _aipwpair, y(ch_`ch'_`h') ///
            d1(onset_def) if1(sample_reentry==1 & onset_nd==0`balflag') ///
            d2(onset_nd)  if2(sample_reentry==1 & onset_def==0`balflag') ///
            omod(`om') pz(`om' `cz_def') reps(`nboot')
        if r(ok) {
            local B1 = r(b1)
            local B2 = r(b2)
            local A1 = r(a1)
            local A2 = r(a2)
            post `Rc' ("`ch'") ("nd")  (`h'+1) (`B2') (`A2') (`B2'-1.96*`A2') (`B2'+1.96*`A2')
            post `Rc' ("`ch'") ("def") (`h'+1) (`B1') (`A1') (`B1'-1.96*`A1') (`B1'+1.96*`A1')
            di "  h=" `h'+1 "  nd=" %8.3f `B2' "   def=" %8.3f `B1'
        }
        else di as error "  h=" `h'+1 ": `ch' estimate failed (cell too thin)."
    }
}
postclose `Rc'

use "`resfc'", clear
label var b  "AIPW ATE (pp)"
label var se "Analytic (unclustered influence-function) SE, matching the paper's own formula"
label var lo "95% CI lower = b - 1.96*se (analytic)"
label var hi "95% CI upper = b + 1.96*se (analytic)"
order channel series horizon b se lo hi
export delimited "$tabs/aipw_channels_reentry.csv", replace
di as result _n "Saved: $tabs/aipw_channels_reentry.csv"

tempfile _chanres
save `_chanres'

* ══════════════════════════════════════════════════════════════════════════
* COMBINED 6-PANEL FIGURE
* ══════════════════════════════════════════════════════════════════════════
local c_nd  "blue"
local c_def "red"

import delimited "$tabs/aipw_results_reentry.csv", clear varnames(1) case(preserve)
keep if inlist(series, "nd", "def")
gen str24 channel = "gdp"
keep channel series horizon b se lo hi
append using `_chanres'

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
        name(combre_`i', replace)
    local ++i
}
capture graph combine combre_1 combre_2 combre_3 combre_4 combre_5 combre_6, ///
    cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
if _rc == 0 {
    graph export "$figs/fig_aipw_combined_reentry.pdf", replace
    di as result "Figure saved: fig_aipw_combined_reentry.pdf (re-entry robustness, sample only, Panel A-F)"
}
else di as error "  ** fig_aipw_combined_reentry failed (rc=" _rc ")"
forvalues i = 1/6 {
    capture graph drop combre_`i'
}

di as result _n "08d_aipw_reentry.do complete."
di as result "Treatment variables (onset_nd/onset_def) are UNCHANGED from the baseline;"
di as result "any difference against 08b_aipw_analyticSE.do/13c_aipw_channels_analyticSE.do's"
di as result "own outputs comes purely from the swapped estimation sample (sample_reentry)."
