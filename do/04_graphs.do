/*===========================================================================
  04_GRAPHS.DO
  Publication-quality IRF figures

  Figure 1: IRF -- All spread crises (Act 1)
  Figure 2: IRF overlay -- Non-default vs. Default-linked (Act 2)
  Figure 3: Pre-trend validation (placebo h=-1 + explicit baseline h=0 + main h=1..5)

  Horizon numbering matches Asonuma et al.: Year 1 is the crisis year itself;
  Year 0 is the explicit pre-crisis baseline (hardcoded 0, no CI).
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
    xlabel(0(1)5, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Year (Year 1 = crisis year)", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP (pp)", size(medsmall)) ///
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
    xlabel(0(1)5, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Year (Year 1 = crisis year)", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP (pp)", size(medsmall)) ///
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
    xline(0.5, lpattern(solid) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(-1(1)5, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Year (Year 0 = pre-crisis baseline, Year 1 = crisis year)", size(medsmall)) ///
    ytitle("Cumulative change in log real GDP (pp)", size(medsmall)) ///
    title("Pre-Trend Test + Main Horizons", size(medium)) ///
    subtitle("All spread crises (N = 61). Placebo: h = -1.", size(small)) ///
    text(0 -1 "Pre-trend (should be ~0)", size(vsmall) color(gray)) ///
    text(0 3 "Post-onset", size(vsmall) color(gray)) ///
    legend(order(2 "Point estimate" 1 "90% CI") ring(0) pos(1) size(small)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig3_pretrend.pdf", replace
graph export "$figs/fig3_pretrend.png", replace width(1200)
di as result "Figure 3 saved."
* ══════════════════════════════════════════════════════════════════════════
* FIGURE 9 — FLOW SPECIFICATION (20_lp_flow.do)
*
* NOTE THE AXIS. This is not Figure 1 with more data. The treatment is being IN
* a crisis in year t, so every treated row is its own reference and h=0 is
* growth DURING that crisis year. h>0 is a forecast horizon from a treated row,
* not time since onset — the label says so deliberately.
* ══════════════════════════════════════════════════════════════════════════
capture confirm file "$clean/irf_flow.dta"
if _rc {
    di as error "  ** irf_flow.dta not found — run 20_lp_flow.do first; Figure 9 skipped."
}
else {
    local c_flow "23 55 94"
    use "$clean/irf_flow.dta", clear
    twoway ///
        (rarea lo95 hi95 horizon, color("`c_flow'%15") lwidth(none)) ///
        (rarea lo90 hi90 horizon, color("`c_flow'%25") lwidth(none)) ///
        (connected b horizon, lcolor("`c_flow'") lwidth(medthick) ///
            mcolor("`c_flow'") msize(medium) msymbol(circle)), ///
        yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
        xlabel(0(1)4, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Forecast horizon h (h=0 = growth during the crisis year)", size(medsmall)) ///
        ytitle("Change in log real GDP (pp)", size(medsmall)) ///
        title("Output While In a Sovereign Spread Crisis", size(medium)) ///
        subtitle("Treatment = every year of an episode (234 country-years, 61 episodes)", size(small)) ///
        note("Local projections, country and year FE, Driscoll-Kraay SE (lag max(2,h+3))." ///
             "Controls: lagged GDP growth, debt, current account, banking-crisis duration," ///
             "government spending, openness, bank credit, hyperinflation flag — all at t-1." ///
             "h>0 is a horizon from a TREATED row, not time since onset: the estimate averages" ///
             "over how long the country has already been in crisis. Not corrected for selection.", ///
             size(vsmall)) ///
        legend(order(3 "Point estimate" 2 "90% CI" 1 "95% CI") ring(0) pos(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig9_irf_flow.pdf", replace
    graph export "$figs/fig9_irf_flow.png", replace width(1200)
    di as result "Figure 9 saved."

    * ── Figure 9b: by resolution type ───────────────────────────────────────
    capture confirm file "$clean/irf_flow_nd.dta"
    if _rc == 0 {
        use "$clean/irf_flow_nd.dta", clear
        append using "$clean/irf_flow_def.dta"
        twoway ///
            (rarea lo90 hi90 horizon if series=="flow_nd",  color("`c_nd'%20")  lwidth(none)) ///
            (rarea lo90 hi90 horizon if series=="flow_def", color("`c_def'%20") lwidth(none)) ///
            (connected b horizon if series=="flow_nd",  lcolor("`c_nd'")  mcolor("`c_nd'")  ///
                msymbol(circle) lwidth(medthick)) ///
            (connected b horizon if series=="flow_def", lcolor("`c_def'") mcolor("`c_def'") ///
                msymbol(square) lpattern(dash) lwidth(medthick)), ///
            yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
            xlabel(0(1)4, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
            xtitle("Forecast horizon h (h=0 = growth during the crisis year)", size(medsmall)) ///
            ytitle("Change in log real GDP (pp)", size(medsmall)) ///
            title("Output While In a Spread Crisis, by Resolution", size(medium)) ///
            subtitle("Joint specification; tranquil country-years omitted", size(small)) ///
            note("90% CIs. Difference tested by lincom, which uses the covariance between the two" ///
                 "coefficients — countries with episodes of both types contribute to both arms." ///
                 "Non-default carries 113 crisis-years, default-linked 121.", size(vsmall)) ///
            legend(order(3 "Non-default" 4 "Default-linked") ring(0) pos(7) cols(1) size(small)) ///
            graphregion(color(white)) plotregion(color(white))
        graph export "$figs/fig9b_irf_flow_resolution.pdf", replace
        di as result "Figure 9b saved."
    }
}

di as result _n "All figures saved to: $figs"
