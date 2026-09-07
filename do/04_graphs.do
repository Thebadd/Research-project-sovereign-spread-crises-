/*===========================================================================
  04_GRAPHS.DO
  Publication-quality IRF figures

  Figure 1: IRF -- All spread crises (Act 1)
  Figure 1A: IRF overlay -- baseline vs. Asonuma-sample robustness (continuation
             years kept as controls; see 02_lp_all.do's robustness block)
  Figure 2: IRF overlay -- Non-default vs. Default-linked (Act 2)
  Figure 2A: same overlay, Asonuma-sample robustness variant (see
             03_lp_resolution.do's robustness block)
  Figure 3: Pre-trend validation (placebo h=-1 + explicit baseline h=0 + main h=1..5)

  Horizon numbering matches Asonuma et al.: Year 1 is the crisis year itself;
  Year 0 is the explicit pre-crisis baseline (hardcoded 0, no CI).
===========================================================================*/

* UNIFORM IRF STYLE (project-wide onset-tier convention): non-default =
* blue, default-linked = red, both solid lines, markers match their line's
* color. Y-axis always "Cumulative percent change", x-axis always "Year",
* title = the plain variable name only.
local c_all  "blue"
local c_nd   "blue"
local c_def  "red"
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
    xtitle("Year", size(medsmall)) ///
    ytitle("Cumulative percent change", size(medsmall)) ///
    title("GDP", size(medium)) ///
    legend(off) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig1_irf_all.pdf",  replace
graph export "$figs/fig1_irf_all.png",  replace width(1200)
di as result "Figure 1 saved."

* ════════════════════════════════════════════════════════════════════════════
* FIGURE 1A: BASELINE vs. ASONUMA-SAMPLE ROBUSTNESS (continuation years kept
* as controls) -- ALL SPREAD CRISES
* ════════════════════════════════════════════════════════════════════════════

capture confirm file "$clean/irf_all_asonumasample.dta"
if _rc {
    di as error _n "  ** irf_all_asonumasample.dta not found -- Figure 1A skipped."
    di as error "     (re-run 02_lp_all.do's Asonuma-sample robustness block first.)"
}
else {
    use "$clean/irf_all.dta", clear
    replace series = "baseline"
    append using "$clean/irf_all_asonumasample.dta"
    keep if horizon >= 0

    twoway ///
        (rarea lo90 hi90 horizon if series=="baseline", ///
            color("`c_all'%20") lwidth(none)) ///
        (connected b horizon if series=="baseline", ///
            lcolor("`c_all'") lwidth(medthick) ///
            mcolor("`c_all'") msize(medium) msymbol(circle)) ///
        (rarea lo90 hi90 horizon if series=="all_asonumasample", ///
            color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if series=="all_asonumasample", ///
            lcolor("`c_def'") lwidth(medthick) ///
            mcolor("`c_def'") msize(medium) msymbol(square)), ///
        yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
        xlabel(0(1)5, labsize(medsmall)) ///
        ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year", size(medsmall)) ///
        ytitle("Cumulative percent change", size(medsmall)) ///
        title("GDP", size(medium)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white))

    graph export "$figs/fig1a_irf_all_asonumasample.pdf", replace
    graph export "$figs/fig1a_irf_all_asonumasample.png", replace width(1200)
    di as result "Figure 1A saved."
}

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
        lcolor("`c_def'") lwidth(medthick) ///
        mcolor("`c_def'") msize(medium) msymbol(square)), ///
    yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
    xlabel(0(1)5, labsize(medsmall)) ///
    ylabel(, format(%4.1f) labsize(medsmall)) ///
    xtitle("Year", size(medsmall)) ///
    ytitle("Cumulative percent change", size(medsmall)) ///
    title("GDP", size(medium)) ///
    legend(off) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig2_irf_resolution.pdf", replace
graph export "$figs/fig2_irf_resolution.png", replace width(1200)
di as result "Figure 2 saved."

* ════════════════════════════════════════════════════════════════════════════
* FIGURE 2A: NON-DEFAULT vs. DEFAULT-LINKED -- ASONUMA-SAMPLE ROBUSTNESS
* (continuation years kept as controls)
* ════════════════════════════════════════════════════════════════════════════

capture confirm file "$clean/irf_nd_asonumasample.dta"
capture confirm file "$clean/irf_def_asonumasample.dta"
if _rc {
    di as error _n "  ** irf_nd/def_asonumasample.dta not found -- Figure 2A skipped."
    di as error "     (re-run 03_lp_resolution.do's Asonuma-sample robustness block first.)"
}
else {
    use "$clean/irf_nd_asonumasample.dta", clear
    append using "$clean/irf_def_asonumasample.dta"
    keep if horizon >= 0

    twoway ///
        (rarea lo90 hi90 horizon if series=="nd_asonumasample", ///
            color("`c_nd'%20") lwidth(none)) ///
        (connected b horizon if series=="nd_asonumasample", ///
            lcolor("`c_nd'") lwidth(medthick) ///
            mcolor("`c_nd'") msize(medium) msymbol(circle)) ///
        (rarea lo90 hi90 horizon if series=="def_asonumasample", ///
            color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if series=="def_asonumasample", ///
            lcolor("`c_def'") lwidth(medthick) ///
            mcolor("`c_def'") msize(medium) msymbol(square)), ///
        yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
        xlabel(0(1)5, labsize(medsmall)) ///
        ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year", size(medsmall)) ///
        ytitle("Cumulative percent change", size(medsmall)) ///
        title("GDP", size(medium)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white))

    graph export "$figs/fig2a_irf_resolution_asonumasample.pdf", replace
    graph export "$figs/fig2a_irf_resolution_asonumasample.png", replace width(1200)
    di as result "Figure 2A saved."
}

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
    xtitle("Year", size(medsmall)) ///
    ytitle("Cumulative percent change", size(medsmall)) ///
    title("GDP", size(medium)) ///
    legend(off) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig3_pretrend.pdf", replace
graph export "$figs/fig3_pretrend.png", replace width(1200)
di as result "Figure 3 saved."
di as result _n "All figures saved to: $figs"
