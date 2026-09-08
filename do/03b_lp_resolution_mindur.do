/*===========================================================================
  03B_LP_RESOLUTION_MINDUR.DO
  Act 2 GDP robustness under the MINIMUM-DURATION episode definition
  (18b_mindur_variant.do): mirrors 03_lp_resolution.do's own headline joint
  regression (`xtreg dy_h onset_nd onset_def $ctrl_core if sample==1 &
  common_abcd==1, fe vce(robust)`) EXACTLY, swapping onset_nd/onset_def for
  onset_nd_mindur/onset_def_mindur -- the 18 non-default onsets that are
  themselves single-year crisis_any spikes are reclassified as tranquil
  (folded into the control pool) rather than treated; 0 default-linked
  onsets are affected. See 03_lp_resolution.do's own header for the full
  spec rationale (country FE only, no year FE, no VIX/UST10Y, matching the
  reference paper's Table I1 design) -- unchanged here, only the treatment
  columns differ.

  STANDALONE: reads $clean/panel_lp_mindur.dta (built by
  18b_mindur_variant.do, NOT panel_lp.dta) -- does not modify
  03_lp_resolution.do or any baseline column. Does not build the full
  Table 2 layout/robustness blocks (Spec A-1/A-2, pre-trend-controlled,
  outturns-only, Asonuma-sample) -- only the headline joint spec needed to
  feed the Act 2 combined-figure build (12b_channels_resolution_mindur.do).

  Output: $clean/irf_nd_mindur.dta, $clean/irf_def_mindur.dta (same schema
  as 03_lp_resolution.do's own irf_nd.dta/irf_def.dta: horizon, b, se,
  lo90, hi90, lo95, hi95, series).
  Run AFTER 18b_mindur_variant.do. Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp_mindur.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
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

capture program drop _pval
program define _pval, rclass
    args b se
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar p = 2*(1 - normal(abs(`b'/`se')))
    }
    else {
        return scalar p = 2*ttail(`dfr', abs(`b'/`se'))
    }
end

di as result _n "=== ACT 2 (MINIMUM-DURATION ROBUSTNESS) — JOINT REGRESSION, GDP ==="

foreach m in b se lo90 hi90 lo95 hi95 {
    matrix `m'_nd  = J(7, 1, .)
    matrix `m'_def = J(7, 1, .)
}

* Pre-trend placebo, same joint spec, displayed h=-1.
foreach h_neg in 2 {
    local row = 3 - `h_neg'
    xtreg dy_m`h_neg' onset_nd_mindur onset_def_mindur `controls_pre' if sample==1 & common_abcd==1, fe vce(robust)
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    foreach g in nd def {
        local bb = _b[onset_`g'_mindur]
        local ss = _se[onset_`g'_mindur]
        _pval `bb' `ss'
        local pp = r(p)
        matrix b_`g'[`row',1]    = `bb'
        matrix se_`g'[`row',1]   = `ss'
        matrix lo90_`g'[`row',1] = `bb' - `c90'*`ss'
        matrix hi90_`g'[`row',1] = `bb' + `c90'*`ss'
        matrix lo95_`g'[`row',1] = `bb' - `c95'*`ss'
        matrix hi95_`g'[`row',1] = `bb' + `c95'*`ss'
        di "h=-1 (`g'): beta = " %6.3f `bb' "  SE = " %6.3f `ss' "  p = " %5.3f `pp'
    }
}

* Explicit baseline (displayed h=0), hardcoded zero, matching Asonuma et al.
foreach g in nd def {
    foreach m in b se lo90 hi90 lo95 hi95 {
        matrix `m'_`g'[2, 1] = 0
    }
}

* Main horizons (displayed h=1..5)
forvalues h = 0/4 {
    local hd  = `h' + 1
    local row = `h' + 3
    xtreg dy_`h' onset_nd_mindur onset_def_mindur `controls' if sample==1 & common_abcd==1, fe vce(robust)

    local bnd  = _b[onset_nd_mindur]
    local bdef = _b[onset_def_mindur]
    local snd  = _se[onset_nd_mindur]
    local sdef = _se[onset_def_mindur]

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

    test onset_nd_mindur = onset_def_mindur
    local pd = r(p)
    local bdiff = `bdef' - `bnd'

    di "h=" `hd' ":  beta_nd=" %6.3f `bnd' ///
               "  beta_def=" %6.3f `bdef' ///
               "  diff(def-nd)=" %6.3f `bdiff' ///
               "  F-stat p="  %5.3f `pd'
}

* ══════════════════════════════════════════════════════════════════════════
* BUILD IRF DATASETS (same schema as 03_lp_resolution.do's own irf_nd.dta/
* irf_def.dta, so 12b_channels_resolution_mindur.do's combined figure can
* import GDP exactly the way 12_channels_resolution.do imports the baseline).
* ══════════════════════════════════════════════════════════════════════════
preserve
    clear
    set obs 7
    gen horizon = _n - 2
    foreach m in b se lo90 hi90 lo95 hi95 {
        svmat `m'_nd, names(`m')
        rename `m'1 `m'
    }
    gen series = "nd"
    save "$clean/irf_nd_mindur.dta", replace

    clear
    set obs 7
    gen horizon = _n - 2
    foreach m in b se lo90 hi90 lo95 hi95 {
        svmat `m'_def, names(`m')
        rename `m'1 `m'
    }
    gen series = "def"
    save "$clean/irf_def_mindur.dta", replace
restore

di as result _n "03b_lp_resolution_mindur.do complete."
di as result "Saved: $clean/irf_nd_mindur.dta, $clean/irf_def_mindur.dta"
