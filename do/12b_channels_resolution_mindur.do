/*===========================================================================
  12B_CHANNELS_RESOLUTION_MINDUR.DO
  Act 2 channel robustness under the MINIMUM-DURATION episode definition
  (18b_mindur_variant.do): mirrors 12_channels_resolution.do's own headline
  joint regression per channel EXACTLY (`xtreg ch_<v>_h onset_nd onset_def
  <ctrl> if sample==1[& common_abcd==1], fe vce(robust)`), swapping
  onset_nd/onset_def for onset_nd_mindur/onset_def_mindur. See
  12_channels_resolution.do's own header for the full spec rationale --
  unchanged here, only the treatment columns differ.

  SCOPE: the five channels the user asked for alongside GDP -- credit, inv,
  claims_govt, fdi, real_lending (govexp/pb excluded, matching the combined
  6-panel figure's own scope, not the full 7-panel figure).

  STANDALONE: reads $clean/panel_lp_mindur.dta (built by
  18b_mindur_variant.do). Does not modify 12_channels_resolution.do or any
  baseline column. Does not build Table 4 or the 7-panel figure -- only the
  per-channel IRF datasets and the combined 6-panel figure (which also
  needs GDP from 03b_lp_resolution_mindur.do, run first).

  Output: $clean/irf_nd_<v>_mindur.dta / irf_def_<v>_mindur.dta (v = credit,
  inv, claims_govt, fdi, real_lending), $figs/fig12_combined_ols_mindur.pdf.
  Run AFTER 18b_mindur_variant.do AND 03b_lp_resolution_mindur.do.
  Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp_mindur.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

local channels credit claims_govt inv fdi real_lending

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

di as result _n "=== ACT 2 (MINIMUM-DURATION ROBUSTNESS) — CHANNELS ==="

foreach ch of local channels {

    local ctrl `ctrl_`ch''

    * Same balanced A-D sample restriction as 12_channels_resolution.do's
    * own headline (credit/inv/claims_govt only; fdi/real_lending keep their
    * own best-available sample).
    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

    di as result _n "CHANNEL: `ch'"
    di "h   b_nd     b_def    p(nd=def)"

    forvalues h = 0/4 {
        local row = `h' + 2

        capture xtreg ch_`ch'_`h' onset_nd_mindur onset_def_mindur `ctrl' ///
            if sample == 1`balflag', fe vce(robust)

        if _rc == 0 {
            matrix b_nd_`ch'[`row',1]    = _b[onset_nd_mindur]
            matrix lo90_nd_`ch'[`row',1] = _b[onset_nd_mindur]  - 1.645*_se[onset_nd_mindur]
            matrix hi90_nd_`ch'[`row',1] = _b[onset_nd_mindur]  + 1.645*_se[onset_nd_mindur]
            matrix lo95_nd_`ch'[`row',1] = _b[onset_nd_mindur]  - 1.960*_se[onset_nd_mindur]
            matrix hi95_nd_`ch'[`row',1] = _b[onset_nd_mindur]  + 1.960*_se[onset_nd_mindur]
            matrix b_def_`ch'[`row',1]   = _b[onset_def_mindur]
            matrix lo90_def_`ch'[`row',1]= _b[onset_def_mindur] - 1.645*_se[onset_def_mindur]
            matrix hi90_def_`ch'[`row',1]= _b[onset_def_mindur] + 1.645*_se[onset_def_mindur]
            matrix lo95_def_`ch'[`row',1]= _b[onset_def_mindur] - 1.960*_se[onset_def_mindur]
            matrix hi95_def_`ch'[`row',1]= _b[onset_def_mindur] + 1.960*_se[onset_def_mindur]
            test onset_nd_mindur = onset_def_mindur
            matrix pval_`ch'[`row',1] = r(p)
            local b_nd_o  = _b[onset_nd_mindur]
            local b_def_o = _b[onset_def_mindur]
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
            gen horizon = _n - 1     // 0 (baseline), 1..5
            foreach m in b lo90 hi90 lo95 hi95 {
                svmat `m'_`grp'_`ch', names(`m')
                rename `m'1 `m'
            }
            gen channel = "`ch'"
            gen group   = "`grp'"
            save "$clean/irf_`grp'_`ch'_mindur.dta", replace
        restore
    }
}

* ══════════════════════════════════════════════════════════════════════════
* COMBINED 6-PANEL FIGURE (mirrors 12_channels_resolution.do Section 8):
* Panel A GDP (from 03b_lp_resolution_mindur.do's saved
* irf_nd_mindur.dta/irf_def_mindur.dta), B Investment, C Bank credit,
* D Claims on government, E FDI, F Real lending rate.
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
        use "$clean/irf_nd_mindur.dta", clear
        rename series group
        tempfile _gdpnd
        save `_gdpnd'
        use "$clean/irf_def_mindur.dta", clear
        rename series group
        append using `_gdpnd'
        keep if horizon >= 0
    }
    else {
        use "$clean/irf_nd_`cv'_mindur.dta",  clear
        append using "$clean/irf_def_`cv'_mindur.dta"
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
        name(combm_`i', replace)

    local ++i
}
graph combine combm_1 combm_2 combm_3 combm_4 combm_5 combm_6, ///
    cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
graph export "$figs/fig12_combined_ols_mindur.pdf", replace
di as result _n "Figure saved: fig12_combined_ols_mindur.pdf (minimum-duration robustness, Panel A-F)"
forvalues i = 1/6 {
    capture graph drop combm_`i'
}

di as result _n "12b_channels_resolution_mindur.do complete."
