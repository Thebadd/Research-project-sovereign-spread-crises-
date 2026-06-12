/*===========================================================================
  12_CHANNELS_RESOLUTION.DO
  Transmission Channels by Resolution Type (Non-Default vs. Default-Linked)
  with IPW Correction for Selection into Default

  Research question: do the transmission channels of spread crises differ
  between episodes resolved without default and those linked to default?

  STRATEGY:
  ---------
  For each of the 6 channels from 11_channels.do, we run:

  Spec A (OLS): joint LP with onset_nd and onset_def simultaneously
    ch_var(h) = αi + γt + β_nd(h)·onset_nd + β_def(h)·onset_def
               + X_channel·δ + ε     [DK SE, lag=max(1,h+1)]

  Spec B (IPW): same regression with Act 2 stabilized weights
    IPW2: probit Pr(default-linked | crisis onset) on debt and CA
    Stabilized weights reweight non-default episodes to match
    observable pre-crisis characteristics of default-linked countries.
    areg absorb(cid) [aw=ipw2] vce(cluster cid)

  Equality test at each horizon: H0: β_nd(h) = β_def(h)

  OUTPUTS:
  --------
  - fig12a_channels_ols.pdf   : 2×3 grid OLS (nd=navy, def=brick)
  - fig12b_channels_ipw.pdf   : 2×3 grid IPW (nd=navy, def=brick)
  - channels_resolution.csv   : full coefficient table
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
* ══════════════════════════════════════════════════════════════════════════

foreach var in credit claims_govt inv govexp pb fdi {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 2. COVERAGE CHECK BY RESOLUTION TYPE
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== DATA COVERAGE AT ONSET BY RESOLUTION TYPE ==="
di as result "  Variable        nd (39)   def (22)"
foreach var in credit claims_govt inv govexp pb fdi {
    quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`var'_0)
    local n_nd = r(N)
    quietly count if onset_def == 1 & sample == 1 & !missing(ch_`var'_0)
    local n_def = r(N)
    di as result "  `var'" _col(20) `n_nd' " / 39" _col(32) `n_def' " / 22"
}

* ══════════════════════════════════════════════════════════════════════════
* 3. IPW WEIGHTS — ACT 2 (replicated from 08_ipw_lp.do)
*    Parsimonious first stage: debt + CA only (26 obs per covariate)
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== FIRST STAGE: Pr(default-linked | crisis onset) ==="

capture drop pscore2 trimmed2 ipw2

probit onset_def debt ca if onset_all == 1, vce(robust)
di as result "McFadden Pseudo-R2: " e(r2_p)

predict pscore2 if onset_all == 1, pr

gen trimmed2 = (pscore2 < 0.05 | pscore2 > 0.95) if !missing(pscore2)
quietly count if trimmed2 == 1
di as result "Observations trimmed: " r(N)
replace pscore2 = . if trimmed2 == 1

quietly summarize onset_def if onset_all == 1
local p_def = r(mean)
local p_nd  = 1 - `p_def'

gen ipw2 = .
replace ipw2 = `p_def' / pscore2       if onset_def == 1 & !missing(pscore2)
replace ipw2 = `p_nd'  / (1 - pscore2) if onset_nd  == 1 & !missing(pscore2)

summarize ipw2 if onset_all == 1, detail
di as result "IPW2 weights: min=" r(min) "  max=" r(max) "  mean=" r(mean)

* ══════════════════════════════════════════════════════════════════════════
* 4. LP ESTIMATION BY CHANNEL — OLS AND IPW
* ══════════════════════════════════════════════════════════════════════════

local channels   credit claims_govt inv govexp pb fdi

* Controls per channel (same as 11_channels.do)
local ctrl_credit      l1_gdpg l2_gdpg debt infl ca banking_crisis reer_chg vix
local ctrl_claims_govt L.claims_govt L.credit pb banking_crisis vix
local ctrl_inv         l1_gdpg l2_gdpg debt ca L.credit banking_crisis reer_chg vix
local ctrl_govexp      L.govexp debt revenue_gdp vix
local ctrl_pb          l1_gdpg l2_gdpg debt ca L.pb banking_crisis vix
local ctrl_fdi         l1_gdpg L.fdi infl reer_chg vix ust10y

* Initialize storage matrices
foreach ch of local channels {
    foreach spec in ols ipw {
        foreach grp in nd def {
            foreach m in b lo90 hi90 lo95 hi95 {
                matrix `m'_`grp'_`spec'_`ch' = J(5, 1, .)
            }
        }
        matrix pval_`spec'_`ch' = J(5, 1, .)
    }
}

* ── Loop over channels ───────────────────────────────────────────────────

foreach ch of local channels {

    local ctrl `ctrl_`ch''

    di as result _n "========================================"
    di as result "CHANNEL: `ch'"
    di as result "========================================"
    di "h   b_nd_OLS  b_def_OLS  p_ols  b_nd_IPW  b_def_IPW  p_ipw"

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local row = `h' + 1

        * ── OLS: DK SE joint regression ─────────────────────────────────
        capture xtscc ch_`ch'_`h' onset_nd onset_def `ctrl' i.year ///
            if sample == 1, fe lag(`lag')

        if _rc == 0 {
            matrix b_nd_ols_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_ols_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_ols_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix lo95_nd_ols_`ch'[`row',1] = _b[onset_nd]  - 1.960*_se[onset_nd]
            matrix hi95_nd_ols_`ch'[`row',1] = _b[onset_nd]  + 1.960*_se[onset_nd]
            matrix b_def_ols_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_ols_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_ols_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix lo95_def_ols_`ch'[`row',1]= _b[onset_def] - 1.960*_se[onset_def]
            matrix hi95_def_ols_`ch'[`row',1]= _b[onset_def] + 1.960*_se[onset_def]
            test onset_nd = onset_def
            matrix pval_ols_`ch'[`row',1] = r(p)
            local b_nd_o  = _b[onset_nd]
            local b_def_o = _b[onset_def]
            local p_o     = r(p)
        }
        else {
            local b_nd_o  = .
            local b_def_o = .
            local p_o     = .
            di as error "OLS failed for `ch' h=`h'"
        }

        * ── IPW: areg with stabilized weights ────────────────────────────
        quietly replace ipw2 = 1 if onset_all == 0 & sample == 1
        capture areg ch_`ch'_`h' onset_nd onset_def `ctrl' i.year ///
            [aw=ipw2] if sample == 1 & !missing(ipw2), ///
            absorb(cid) vce(cluster cid)
        quietly replace ipw2 = . if onset_all == 0

        if _rc == 0 {
            matrix b_nd_ipw_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_ipw_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_ipw_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix lo95_nd_ipw_`ch'[`row',1] = _b[onset_nd]  - 1.960*_se[onset_nd]
            matrix hi95_nd_ipw_`ch'[`row',1] = _b[onset_nd]  + 1.960*_se[onset_nd]
            matrix b_def_ipw_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_ipw_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_ipw_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix lo95_def_ipw_`ch'[`row',1]= _b[onset_def] - 1.960*_se[onset_def]
            matrix hi95_def_ipw_`ch'[`row',1]= _b[onset_def] + 1.960*_se[onset_def]
            test onset_nd = onset_def
            matrix pval_ipw_`ch'[`row',1] = r(p)
            local b_nd_w  = _b[onset_nd]
            local b_def_w = _b[onset_def]
            local p_w     = r(p)
        }
        else {
            local b_nd_w  = .
            local b_def_w = .
            local p_w     = .
            di as error "IPW failed for `ch' h=`h'"
        }

        di "h=" `h' ///
           "  " %7.3f `b_nd_o'  "  " %7.3f `b_def_o' "  " %5.3f `p_o' ///
           "  " %7.3f `b_nd_w'  "  " %7.3f `b_def_w' "  " %5.3f `p_w'
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 5. SAVE IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

foreach ch of local channels {
    foreach spec in ols ipw {
        foreach grp in nd def {
            preserve
                clear
                set obs 5
                gen horizon = _n - 1
                foreach m in b lo90 hi90 lo95 hi95 {
                    svmat `m'_`grp'_`spec'_`ch', names(`m')
                    rename `m'1 `m'
                }
                gen channel  = "`ch'"
                gen spec_tag = "`spec'"
                gen group    = "`grp'"
                save "$clean/irf_`grp'_`spec'_`ch'.dta", replace
            restore
        }
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 6. EXPORT SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

preserve
    clear
    local nrows = 5 * 6   // 5 horizons × 6 channels
    set obs `nrows'
    gen channel  = ""
    gen horizon  = .
    foreach v in b_nd_ols b_def_ols p_ols b_nd_ipw b_def_ipw p_ipw {
        gen `v' = .
    }

    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel  = "`ch'"                         in `row'
            replace horizon  = `h'                            in `row'
            replace b_nd_ols  = b_nd_ols_`ch'[`h'+1,1]       in `row'
            replace b_def_ols = b_def_ols_`ch'[`h'+1,1]      in `row'
            replace p_ols     = pval_ols_`ch'[`h'+1,1]       in `row'
            replace b_nd_ipw  = b_nd_ipw_`ch'[`h'+1,1]       in `row'
            replace b_def_ipw = b_def_ipw_`ch'[`h'+1,1]      in `row'
            replace p_ipw     = pval_ipw_`ch'[`h'+1,1]       in `row'
            local ++row
        }
    }

    order channel horizon b_nd_ols b_def_ols p_ols b_nd_ipw b_def_ipw p_ipw
    export delimited "$tabs/channels_resolution.csv", replace
    di as result "Table saved: $tabs/channels_resolution.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 7. FIGURES — 2×3 MULTI-PANEL, OLS THEN IPW
* ══════════════════════════════════════════════════════════════════════════

local c_nd  "23 55 94"    // navy   — non-default
local c_def "180 60 40"   // brick  — default-linked

local titlelabels `" "Private Credit/GDP" "Bank Claims on Govt/GDP" "Investment/GDP" "Govt Expenditure/GDP" "Primary Balance/GDP" "FDI/GDP" "'

* ── Figure A: OLS ────────────────────────────────────────────────────────

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titlelabels'

    use "$clean/irf_nd_ols_`ch'.dta",  clear
    append using "$clean/irf_def_ols_`ch'.dta"

    * Show legend only in last panel (bottom-right)
    if `i' == 6 {
        local legopt legend(order(3 "Non-default" 4 "Default-linked") ///
                     ring(0) pos(7) cols(1) size(vsmall) ///
                     region(lcolor(gs12)))
    }
    else {
        local legopt legend(off)
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
            lwidth(medthick) msize(small) lpattern(dash)) ///
        , ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)4, labsize(small)) ///
        ylabel(, format(%5.1f) labsize(small)) ///
        xtitle("Years after onset", size(vsmall)) ///
        ytitle("Cumulative change (pp)", size(vsmall)) ///
        title(`tlab', size(small) color(navy)) ///
        `legopt' ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ols_`i', replace)

    local ++i
}

graph combine ols_1 ols_2 ols_3 ols_4 ols_5 ols_6, ///
    cols(3) rows(2) ///
    title("Transmission Channels by Resolution Type — OLS", ///
          size(medlarge) color(navy)) ///
    note("90% CI. DK SE. Country & year FE. Non-default: N=39. Default-linked: N=22.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig12a_channels_ols.pdf", replace
di as result "Figure saved: fig12a_channels_ols.pdf"
forvalues i = 1/6 {
    capture graph drop ols_`i'
}

* ── Figure B: IPW ────────────────────────────────────────────────────────

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titlelabels'

    use "$clean/irf_nd_ipw_`ch'.dta",  clear
    append using "$clean/irf_def_ipw_`ch'.dta"

    * Show legend only in last panel (bottom-right)
    if `i' == 6 {
        local legopt legend(order(3 "Non-default (IPW)" 4 "Default-linked (IPW)") ///
                     ring(0) pos(7) cols(1) size(vsmall) ///
                     region(lcolor(gs12)))
    }
    else {
        local legopt legend(off)
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
            lwidth(medthick) msize(small) lpattern(dash)) ///
        , ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)4, labsize(small)) ///
        ylabel(, format(%5.1f) labsize(small)) ///
        xtitle("Years after onset", size(vsmall)) ///
        ytitle("Cumulative change (pp)", size(vsmall)) ///
        title(`tlab', size(small) color(navy)) ///
        `legopt' ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ipw_`i', replace)

    local ++i
}

graph combine ipw_1 ipw_2 ipw_3 ipw_4 ipw_5 ipw_6, ///
    cols(3) rows(2) ///
    title("Transmission Channels by Resolution Type — IPW", ///
          size(medlarge) color(navy)) ///
    note("90% CI. Clustered SE. Country & year FE. IPW: probit(def | crisis) on debt, CA.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig12b_channels_ipw.pdf", replace
di as result "Figure saved: fig12b_channels_ipw.pdf"
forvalues i = 1/6 {
    capture graph drop ipw_`i'
}

di as result _n "12_channels_resolution.do complete."
