/*===========================================================================
  11_CHANNELS.DO
  Transmission Channels of Sovereign Spread Crises

  Research question: through which channels does a spread crisis affect
  real output? We run the same LP framework as Act 1 but replace the
  GDP outcome with six channel variables, each with controls tailored
  to the specific mechanism.

  CHANNELS AND OUTCOMES:
  ----------------------
  1. Credit channel          → credit_gdp  (private credit / GDP)
  2. Sovereign-bank nexus    → claims_govt (bank claims on govt / GDP)
  3. Investment channel      → inv         (gross investment / GDP)
  4. Fiscal channel          → govexp      (government expenditure / GDP)
  5. Fiscal adjustment       → pb          (primary balance / GDP)
  6. External financing      → fdi         (FDI / GDP)

  SPECIFICATION (horizon h = 0,...,4):
  -------------------------------------
  Channel_var(i,t+h) - Channel_var(i,t-1) = αi + γt + β(h)·onset_all
                                            + X_channel · δ + ε(i,t+h)

  All outcomes expressed as cumulative change from t-1 (same anchoring
  as main LP). Controls are pre-determined (lagged relative to t).
  DK SE with lag = max(1, h+1) throughout.

  CONTROLS BY CHANNEL (see discussion in paper):
  -----------------------------------------------
  Credit:          l1_gdpg l2_gdpg debt infl vix ust10y
  Sovereign-bank:  l1_gdpg debt pb vix ust10y
  Investment:      l1_gdpg l2_gdpg debt ca L.credit vix ust10y
  Govt expenditure:l1_gdpg debt ca pb vix ust10y
  Primary balance: l1_gdpg l2_gdpg debt ca vix ust10y
  FDI:             l1_gdpg debt ca vix ust10y

  Saves:
    "$clean/irf_ch_*.dta"        — one IRF dataset per channel
    "$figs/fig11_channels.pdf"   — 2×3 multi-panel figure
    "$tabs/channels_summary.csv" — coefficients and SE for all channels
===========================================================================*/

use "$clean/panel_lp.dta", clear

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
*    Same construction as dy_h: cumulative change relative to t-1
* ══════════════════════════════════════════════════════════════════════════

foreach var in credit claims_govt inv govexp pb fdi {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
        label var ch_`var'_`h' "Cum. change `var': F`h' vs t-1"
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 2. COVERAGE CHECK
*    Report how many onset observations have non-missing channel data
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== DATA COVERAGE AT ONSET (sample==1 & onset_all==1) ==="
foreach var in credit claims_govt inv govexp pb fdi {
    quietly count if onset_all == 1 & sample == 1 & !missing(ch_`var'_0)
    di as result "  `var': " r(N) " / 61 onsets with non-missing data at h=0"
}

* ══════════════════════════════════════════════════════════════════════════
* 3. LP ESTIMATION BY CHANNEL
* ══════════════════════════════════════════════════════════════════════════

* Storage: one 5×5 matrix per channel (rows=horizons, cols=b se lo90 hi90 lo95 hi95)
* We store b, se, lo90, hi90, lo95, hi95 in separate 5×1 matrices

local channels credit claims_govt inv govexp pb fdi

foreach ch of local channels {
    foreach m in b lo90 hi90 lo95 hi95 {
        matrix `m'_`ch' = J(5, 1, .)
    }
}

* ── Channel 1: Credit ────────────────────────────────────────────────────
di as result _n "=== CHANNEL 1: PRIVATE CREDIT / GDP ==="
local ctrl_credit l1_gdpg l2_gdpg debt infl vix ust10y

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_credit_`h' onset_all `ctrl_credit' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_credit[`h'+1,1]    = _b[onset_all]
        matrix lo90_credit[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_credit[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_credit[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_credit[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for credit (rc=" _rc ")"
}

* ── Channel 2: Sovereign-bank nexus ──────────────────────────────────────
di as result _n "=== CHANNEL 2: BANK CLAIMS ON GOVT / GDP ==="
local ctrl_claims l1_gdpg debt pb vix ust10y

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_claims_govt_`h' onset_all `ctrl_claims' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_claims_govt[`h'+1,1]    = _b[onset_all]
        matrix lo90_claims_govt[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_claims_govt[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_claims_govt[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_claims_govt[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for claims_govt (rc=" _rc ")"
}

* ── Channel 3: Investment ────────────────────────────────────────────────
di as result _n "=== CHANNEL 3: GROSS INVESTMENT / GDP ==="
local ctrl_inv l1_gdpg l2_gdpg debt ca L.credit vix ust10y

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_inv_`h' onset_all `ctrl_inv' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_inv[`h'+1,1]    = _b[onset_all]
        matrix lo90_inv[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_inv[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_inv[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_inv[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for inv (rc=" _rc ")"
}

* ── Channel 4: Government expenditure ────────────────────────────────────
di as result _n "=== CHANNEL 4: GOVERNMENT EXPENDITURE / GDP ==="
local ctrl_govexp l1_gdpg debt ca pb vix ust10y

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_govexp_`h' onset_all `ctrl_govexp' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_govexp[`h'+1,1]    = _b[onset_all]
        matrix lo90_govexp[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_govexp[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_govexp[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_govexp[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for govexp (rc=" _rc ")"
}

* ── Channel 5: Primary balance ───────────────────────────────────────────
di as result _n "=== CHANNEL 5: PRIMARY BALANCE / GDP ==="
local ctrl_pb l1_gdpg l2_gdpg debt ca vix ust10y

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_pb_`h' onset_all `ctrl_pb' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_pb[`h'+1,1]    = _b[onset_all]
        matrix lo90_pb[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_pb[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_pb[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_pb[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for pb (rc=" _rc ")"
}

* ── Channel 6: FDI ───────────────────────────────────────────────────────
di as result _n "=== CHANNEL 6: FDI / GDP ==="
local ctrl_fdi l1_gdpg debt ca vix ust10y

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_fdi_`h' onset_all `ctrl_fdi' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_fdi[`h'+1,1]    = _b[onset_all]
        matrix lo90_fdi[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_fdi[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_fdi[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_fdi[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for fdi (rc=" _rc ")"
}

* ══════════════════════════════════════════════════════════════════════════
* 4. SAVE IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

local channels      credit claims_govt inv govexp pb fdi
local chan_labels   "Private credit/GDP" "Bank claims on govt/GDP" ///
                    "Investment/GDP" "Govt expenditure/GDP" ///
                    "Primary balance/GDP" "FDI/GDP"

local i = 1
foreach ch of local channels {
    preserve
        clear
        set obs 5
        gen horizon = _n - 1
        foreach m in b lo90 hi90 lo95 hi95 {
            svmat `m'_`ch', names(`m')
            rename `m'1 `m'
        }
        gen channel = "`ch'"
        save "$clean/irf_ch_`ch'.dta", replace
    restore
    local ++i
}

* ══════════════════════════════════════════════════════════════════════════
* 5. EXPORT SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

preserve
    clear
    local nobs = 5 * 6    // 5 horizons × 6 channels
    set obs `nobs'
    gen channel = ""
    gen horizon = .
    gen b       = .
    gen lo95    = .
    gen hi95    = .
    gen lo90    = .
    gen hi90    = .

    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel = "`ch'"              in `row'
            replace horizon = `h'                 in `row'
            replace b       = b_`ch'[`h'+1, 1]   in `row'
            replace lo95    = lo95_`ch'[`h'+1, 1] in `row'
            replace hi95    = hi95_`ch'[`h'+1, 1] in `row'
            replace lo90    = lo90_`ch'[`h'+1, 1] in `row'
            replace hi90    = hi90_`ch'[`h'+1, 1] in `row'
            local ++row
        }
    }

    order channel horizon b lo90 hi90 lo95 hi95
    export delimited "$tabs/channels_summary.csv", replace
    di as result "Channel summary saved: $tabs/channels_summary.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 6. FIGURE — 2×3 MULTI-PANEL IRF
* ══════════════════════════════════════════════════════════════════════════

local c_main "23 55 94"

local channels    credit claims_govt inv govexp pb fdi
local titlelabels `" "Private Credit/GDP" "Bank Claims on Govt/GDP" "Investment/GDP" "Govt Expenditure/GDP" "Primary Balance/GDP" "FDI/GDP" "'
local fignames    fig11a fig11b fig11c fig11d fig11e fig11f

local i = 1
foreach ch of local channels {

    local tlab : word `i' of `titlelabels'

    use "$clean/irf_ch_`ch'.dta", clear

    twoway ///
        (rarea lo95 hi95 horizon, ///
            color("`c_main'%15") lwidth(none)) ///
        (rarea lo90 hi90 horizon, ///
            color("`c_main'%28") lwidth(none)) ///
        (connected b horizon, ///
            lcolor("`c_main'") mcolor("`c_main'") ///
            msymbol(circle) lwidth(medthick) msize(medium)) ///
        , ///
        yline(0, lcolor(gs8) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)4, labsize(medsmall)) ///
        ylabel(, format(%5.2f) labsize(medsmall)) ///
        xtitle("Years after onset", size(small)) ///
        ytitle("Cumulative change (pp)", size(small)) ///
        title(`tlab', size(medsmall) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(`: word `i' of `fignames'', replace)

    local ++i
}

* Combine into 2×3 grid
graph combine fig11a fig11b fig11c fig11d fig11e fig11f, ///
    cols(3) rows(2) ///
    title("Transmission Channels of Sovereign Spread Crises", ///
          size(medlarge) color(navy)) ///
    note("90% and 95% CI. DK SE. Country & year FE. Treatment: all 61 onset episodes.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig11_channels.pdf", replace
di as result "Figure saved: fig11_channels.pdf"

* Close individual graphs from memory
foreach nm of local fignames {
    graph drop `nm'
}

di as result _n "11_channels.do complete."
