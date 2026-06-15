/*===========================================================================
  04_GRAPHS.DO
  Publication-quality IRF figures

  Figure 1: IRF -- All spread crises (Act 1)
  Figure 2: IRF overlay -- Non-default vs. Default-linked (Act 2)
  Figure 3: Pre-trend validation (placebo horizons h=-2,-1 + main h=0..4)
===========================================================================*/

* Color scheme
local c_all  "23 55 94"
local c_nd   "0 84 166"
local c_def  "157 36 73"
local c_zero "150 150 150"

* ════════════════════════════════════════════════════════════════════════════
* FIGURE 1: ALL SPREAD CRISES
* ════════════════════════════════════════════════════════════════════════════

use "$clean/irf_all.dta", clear
keep if horizon >= 0

twoway ///
    (rarea lo95 hi95 horizon, ///
        color("`c_all'%15") lwidth(none)) ///
    (rarea lo90 hi90 horizon, ///
        color("`c_all'%25") lwidth(none)) ///
    (connected b horizon, ///
        lcolor("`c_all'") lwidth(medthick) ///
        mcolor("`c_all'") msize(medium) msymbol(circle)), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years after crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
    title("Output Cost of Sovereign Spread Crises", size(medium)) ///
    subtitle("All episodes (N = 61), 52 EM economies, 1994-2025", size(small)) ///
    note("Local projections (Jorda 2005). Driscoll-Kraay SE." ///
         "Country and year FE. Controls: lagged GDP growth, CA/GDP," ///
         "debt/GDP, inflation, IMF program, VIX, US 10yr yield.", ///
         size(vsmall)) ///
    legend(order(3 "Point estimate" 2 "90% CI" 1 "95% CI") ///
           ring(0) pos(1) size(small)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig1_irf_all.pdf",  replace
graph export "$figs/fig1_irf_all.png",  replace width(1200)
di as result "Figure 1 saved."

* ════════════════════════════════════════════════════════════════════════════
* FIGURE 2: NON-DEFAULT vs. DEFAULT-LINKED -- OVERLAY
* ════════════════════════════════════════════════════════════════════════════

use "$clean/irf_nd.dta", clear
append using "$clean/irf_def.dta"
keep if horizon >= 0

twoway ///
    (rarea lo90 hi90 horizon if series=="nd", ///
        color("`c_nd'%20") lwidth(none)) ///
    (connected b horizon if series=="nd", ///
        lcolor("`c_nd'") lwidth(medthick) ///
        mcolor("`c_nd'") msize(medium) msymbol(circle)) ///
    (rarea lo90 hi90 horizon if series=="def", ///
        color("`c_def'%20") lwidth(none)) ///
    (connected b horizon if series=="def", ///
        lcolor("`c_def'") lwidth(medthick) lpattern(dash) ///
        mcolor("`c_def'") msize(medium) msymbol(square)), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years after crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
    title("Output Cost: Non-Default vs. Default-Linked Crises", size(medium)) ///
    subtitle("52 EM economies, 1994-2025", size(small)) ///
    legend(order(2 "Non-default (N=40)" 4 "Default-linked (N=21)") ///
           ring(0) pos(3) size(small)) ///
    note("Driscoll-Kraay SE, lag = h+1. Country + year FE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig2_irf_resolution.pdf", replace
graph export "$figs/fig2_irf_resolution.png", replace width(1200)
di as result "Figure 2 saved."

* ════════════════════════════════════════════════════════════════════════════
* FIGURE 3: PRE-TREND + MAIN HORIZONS (ALL CRISES)
* ════════════════════════════════════════════════════════════════════════════

use "$clean/irf_all.dta", clear
keep if series == "all"

twoway ///
    (rarea lo90 hi90 horizon, ///
        color("`c_all'%25") lwidth(none)) ///
    (connected b horizon, ///
        lcolor("`c_all'") lwidth(medthick) ///
        mcolor("`c_all'") msize(medium) msymbol(circle)), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xline(-0.5, lpattern(solid) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(-2(1)4, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Years relative to crisis onset", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP p.c. (pp)", size(medsmall)) ///
    title("Pre-Trend Test + Main Horizons", size(medium)) ///
    subtitle("All spread crises (N = 61). Placebo: h = -2, -1.", size(small)) ///
    text(0 -1.5 "Pre-trend (should be ~0)", size(vsmall) color(gray)) ///
    text(0 2 "Post-onset", size(vsmall) color(gray)) ///
    legend(order(2 "Point estimate" 1 "90% CI") ring(0) pos(1) size(small)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig3_pretrend.pdf", replace
graph export "$figs/fig3_pretrend.png", replace width(1200)
di as result "Figure 3 saved."
di as result _n "All figures saved to: $figs"
