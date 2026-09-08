/*===========================================================================
  12C_CHANNELS_RESOLUTION_REENTRY.DO
  Act 2 GDP + channel robustness under the RE-ENTRY episode definition
  (18c_reentry_variant.do), ALL IN ONE FILE: GDP's own headline joint
  regression (mirroring 03_lp_resolution.do) runs first, then the five
  channels (mirroring 12_channels_resolution.do), then the combined 6-panel
  figure -- self-contained, no dependency on a separate GDP file.

  UNLIKE the minimum-duration robustness build (12b_channels_resolution_
  mindur.do), the ONLY thing that changes here is the ESTIMATION SAMPLE --
  `if sample==1' -> `if sample_reentry==1' -- NOT the treatment variables
  (onset_nd/onset_def stay exactly as they are, unchanged). This isolates
  the control-pool-composition channel through which the re-entry episode
  definition actually affects the onset-tier design (see
  18c_reentry_variant.do's own header for the full correction/argument):
  sample_reentry keeps every onset row plus every year not still "in
  crisis" under the re-entry rule, dropping any post-onset year re-entry
  still counts as ongoing regardless of the baseline's own `continuation'
  flag. Any difference between this file's coefficients and the baseline
  Table 2/Table 4 numbers is therefore coming purely from the changed
  control pool, not from a different treatment definition.

  SCOPE: GDP + the five channels used throughout this session's other
  robustness builds -- Investment, Bank credit, Bank claims on government,
  FDI, Real lending rate (govexp/pb excluded).

  STANDALONE: reads $clean/panel_lp_reentry.dta (built by
  18c_reentry_variant.do). Does not modify 03_lp_resolution.do,
  12_channels_resolution.do, 18c_reentry_variant.do, or any baseline
  column.

  Output: $figs/fig12_combined_ols_reentry.pdf (Panel A-F: GDP, Investment,
          Bank credit, Claims on government, FDI, Real lending rate).
  Run AFTER 18c_reentry_variant.do. Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp_reentry.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year
local controls $ctrl_core
local cc         $ctrl_core
local dropgdpg   l1_gdpg
local controls_pre : list cc - dropgdpg

capture program drop _critvals
program define _critvals, rclass
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar df  = .
        return scalar c90 = 1.645
        return scalar c95 = 1.960
    }
    else {
        return scalar df  = `dfr'
        return scalar c90 = invttail(`dfr', 0.05)
        return scalar c95 = invttail(`dfr', 0.025)
    }
end

* ══════════════════════════════════════════════════════════════════════════
* GDP — mirrors 03_lp_resolution.do's own headline joint regression,
* sample swapped, treatment UNCHANGED
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== ACT 2 (RE-ENTRY ROBUSTNESS: SAMPLE ONLY) — JOINT REGRESSION, GDP ==="

foreach m in b se lo90 hi90 lo95 hi95 {
    matrix `m'_nd  = J(7, 1, .)
    matrix `m'_def = J(7, 1, .)
}

foreach h_neg in 2 {
    local row = 3 - `h_neg'
    xtreg dy_m`h_neg' onset_nd onset_def `controls_pre' if sample_reentry==1 & common_abcd==1, fe vce(robust)
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    foreach g in nd def {
        local bb = _b[onset_`g']
        local ss = _se[onset_`g']
        matrix b_`g'[`row',1]    = `bb'
        matrix se_`g'[`row',1]   = `ss'
        matrix lo90_`g'[`row',1] = `bb' - `c90'*`ss'
        matrix hi90_`g'[`row',1] = `bb' + `c90'*`ss'
        matrix lo95_`g'[`row',1] = `bb' - `c95'*`ss'
        matrix hi95_`g'[`row',1] = `bb' + `c95'*`ss'
    }
}
foreach g in nd def {
    foreach m in b se lo90 hi90 lo95 hi95 {
        matrix `m'_`g'[2, 1] = 0
    }
}
forvalues h = 0/4 {
    local hd  = `h' + 1
    local row = `h' + 3
    xtreg dy_`h' onset_nd onset_def `controls' if sample_reentry==1 & common_abcd==1, fe vce(robust)
    local bnd  = _b[onset_nd]
    local bdef = _b[onset_def]
    local snd  = _se[onset_nd]
    local sdef = _se[onset_def]
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    matrix b_nd[`row',1]    = `bnd'
    matrix se_nd[`row',1]   = `snd'
    matrix lo90_nd[`row',1] = `bnd' - `c90'*`snd'
    matrix hi90_nd[`row',1] = `bnd' + `c90'*`snd'
    matrix lo95_nd[`row',1] = `bnd' - `c95'*`snd'
    matrix hi95_nd[`row',1] = `bnd' + `c95'*`snd'
    matrix b_def[`row',1]    = `bdef'
    matrix se_def[`row',1]   = `sdef'
    matrix lo90_def[`row',1] = `bdef' - `c90'*`sdef'
    matrix hi90_def[`row',1] = `bdef' + `c90'*`sdef'
    matrix lo95_def[`row',1] = `bdef' - `c95'*`sdef'
    matrix hi95_def[`row',1] = `bdef' + `c95'*`sdef'
    test onset_nd = onset_def
    di "h=" `hd' ":  beta_nd=" %6.3f `bnd' "  beta_def=" %6.3f `bdef' "  F-stat p=" %5.3f r(p)
}

preserve
    clear
    set obs 7
    gen horizon = _n - 2
    foreach m in b se lo90 hi90 lo95 hi95 {
        svmat `m'_nd, names(`m')
        rename `m'1 `m'
    }
    gen series = "nd"
    save "$clean/irf_nd_reentry.dta", replace

    clear
    set obs 7
    gen horizon = _n - 2
    foreach m in b se lo90 hi90 lo95 hi95 {
        svmat `m'_def, names(`m')
        rename `m'1 `m'
    }
    gen series = "def"
    save "$clean/irf_def_reentry.dta", replace
restore

* ══════════════════════════════════════════════════════════════════════════
* CHANNELS — mirrors 12_channels_resolution.do's own headline loop,
* sample swapped, treatment UNCHANGED
* ══════════════════════════════════════════════════════════════════════════
local channels credit claims_govt inv fdi real_lending

* Channel outcomes ch_v_h = F h.v - L.v (h=0..4), same construction as
* 12_channels_resolution.do -- panel_lp_reentry.dta does not carry these
* persistently (see 12b_channels_resolution_mindur.do's own fix for the
* identical gap, independently confirmed there this session).
foreach v of local channels {
    local src `v'
    if inlist("`v'","credit","inv") local src ln_r_`v'
    capture drop `v'_base
    gen double `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen double ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    capture drop pre_`v'
    gen double pre_`v' = L.`src' - L2.`src'
}

local ctrl_credit      l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
local ctrl_claims_govt $ctrl_core pre_claims_govt
local ctrl_inv         $ctrl_core pre_inv
local ctrl_fdi         $ctrl_core pre_fdi
local ctrl_real_lending $ctrl_core pre_real_lending

foreach ch of local channels {
    foreach grp in nd def {
        foreach m in b lo90 hi90 lo95 hi95 {
            matrix `m'_`grp'_`ch' = J(6, 1, 0)
        }
    }
    matrix pval_`ch' = J(6, 1, .)
}

di as result _n "=== ACT 2 (RE-ENTRY ROBUSTNESS: SAMPLE ONLY) — CHANNELS ==="

foreach ch of local channels {

    local ctrl `ctrl_`ch''

    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

    di as result _n "CHANNEL: `ch'"
    di "h   b_nd     b_def    p(nd=def)"

    forvalues h = 0/4 {
        local row = `h' + 2

        capture xtreg ch_`ch'_`h' onset_nd onset_def `ctrl' ///
            if sample_reentry == 1`balflag', fe vce(robust)

        if _rc == 0 {
            matrix b_nd_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix lo95_nd_`ch'[`row',1] = _b[onset_nd]  - 1.960*_se[onset_nd]
            matrix hi95_nd_`ch'[`row',1] = _b[onset_nd]  + 1.960*_se[onset_nd]
            matrix b_def_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix lo95_def_`ch'[`row',1]= _b[onset_def] - 1.960*_se[onset_def]
            matrix hi95_def_`ch'[`row',1]= _b[onset_def] + 1.960*_se[onset_def]
            test onset_nd = onset_def
            matrix pval_`ch'[`row',1] = r(p)
            local b_nd_o  = _b[onset_nd]
            local b_def_o = _b[onset_def]
            local p_o     = r(p)
        }
        else {
            local b_nd_o  = .
            local b_def_o = .
            local p_o     = .
            di as error "regression failed for `ch' h=`=`h'+1'"
        }

        di "h=" `h'+1 "  " %7.3f `b_nd_o'  "  " %7.3f `b_def_o' "  " %5.3f `p_o'
    }
}

* ══════════════════════════════════════════════════════════════════════════
* SAVE IRF DATASETS (same schema as 12_channels_resolution.do's own
* irf_nd_<v>.dta/irf_def_<v>.dta)
* ══════════════════════════════════════════════════════════════════════════
foreach ch of local channels {
    foreach grp in nd def {
        preserve
            clear
            set obs 6
            gen horizon = _n - 1
            foreach m in b lo90 hi90 lo95 hi95 {
                svmat `m'_`grp'_`ch', names(`m')
                rename `m'1 `m'
            }
            gen channel = "`ch'"
            gen group   = "`grp'"
            save "$clean/irf_`grp'_`ch'_reentry.dta", replace
        restore
    }
}

* ══════════════════════════════════════════════════════════════════════════
* COMBINED 6-PANEL FIGURE
* ══════════════════════════════════════════════════════════════════════════
local c_nd  "blue"
local c_def "red"

local combo_vars   gdp inv credit claims_govt fdi real_lending
local combo_labels `" "Panel A: GDP" "Panel B: Investment" "Panel C: Bank credit" "Panel D: Claims on govt" "Panel E: FDI" "Panel F: Real lending rate" "'
local i = 1
foreach cv of local combo_vars {
    local clab : word `i' of `combo_labels'
    local ytit ""
    if inlist(`i', 1, 4) local ytit "Cumulative percent change"

    if "`cv'" == "gdp" {
        use "$clean/irf_nd_reentry.dta", clear
        rename series group
        tempfile _gdpnd
        save `_gdpnd'
        use "$clean/irf_def_reentry.dta", clear
        rename series group
        append using `_gdpnd'
        keep if horizon >= 0
    }
    else {
        use "$clean/irf_nd_`cv'_reentry.dta",  clear
        append using "$clean/irf_def_`cv'_reentry.dta"
    }

    twoway ///
        (rarea lo90 hi90 horizon if group=="nd", ///
            color("`c_nd'%20") lwidth(none)) ///
        (rarea lo90 hi90 horizon if group=="def", ///
            color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if group=="nd", ///
            lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) ///
            lwidth(medthick) msize(small)) ///
        (connected b horizon if group=="def", ///
            lcolor("`c_def'") mcolor("`c_def'") msymbol(square) ///
            lwidth(medthick) msize(small)) ///
        , ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)5, labsize(large)) ///
        ylabel(, format(%9.0f) labsize(large) angle(horizontal)) ///
        xtitle("Year", size(large)) ///
        ytitle("`ytit'", size(large)) ///
        title("`clab'", size(medlarge) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(combr_`i', replace)

    local ++i
}
graph combine combr_1 combr_2 combr_3 combr_4 combr_5 combr_6, ///
    cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
graph export "$figs/fig12_combined_ols_reentry.pdf", replace
di as result _n "Figure saved: fig12_combined_ols_reentry.pdf (re-entry robustness, sample only, Panel A-F)"
forvalues i = 1/6 {
    capture graph drop combr_`i'
}

di as result _n "12c_channels_resolution_reentry.do complete."
