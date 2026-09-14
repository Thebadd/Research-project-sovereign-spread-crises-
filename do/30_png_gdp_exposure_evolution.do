/*===========================================================================
  30_PNG_GDP_EXPOSURE_EVOLUTION.DO
  EXPLORATORY / ROBUSTNESS-ONLY TEST -- NOT part of the headline six-variable
  design and NOT wired into 00_master.do, matching this project's own
  convention for standalone exploratory files (27_png_debt_test.do,
  29_png_stock_gdp_exposure.do, etc.).

  PURPOSE. A pre/post-crisis DESCRIPTIVE evolution figure for the exposure
  indicator built in 29_png_stock_gdp_exposure.do (png_gdp = PNG debt stock
  / GDP, current US$ both sides -- see that file's header for the
  denominator-fix history), mirroring EXACTLY the style and construction
  27_png_debt_test.do's Part 1 evolution figure uses for the flow-based
  exposure measure: Years -3..5 relative to onset, country-demeaned,
  rebased so Year 0 = the LAST PRE-CRISIS year (= 0 by construction) and
  Year 1 = the crisis year itself, matching 02a_descriptive_facts.do's own
  convention exactly (not this project's usual h=0-is-the-onset-year
  convention used in the OLS/AIPW tables -- a deliberate one-off exception
  for this figure, same choice already made in 27's own evolution figure).

  This file does NOT re-run any regression (no OLS, no AIPW) -- it is
  purely the descriptive companion, answering "how does pre-crisis
  external-debt exposure (stock/GDP) itself move around a spread crisis,"
  before/independent of any test of whether it predicts GDP/credit/
  investment outcomes (that is 29's job, not this file's).

  EXPOSURE MEASURE, self-contained (re-imports independently, does not
  depend on 29 having run first in the same session): png_gdp = PNG debt
  stock (DT.DOD.DPNG.CD, data/raw/PNGtoGNI.xlsx) / GDP in current US$ (WEO
  NGDPD, data/raw/WEOApr2026all.xlsx) * 100 -- the SAME corrected
  construction 29 uses (NOT gdp_real, which is domestic-currency and was
  found to produce nonsense when tried -- see 29's header for the full
  diagnosis).

  OUTCOME PLOTTED: the LEVEL of png_gdp at each year (not a difference from
  a base year) -- a stock ratio like this is naturally read as a level
  path over time (debt/GDP evolution), matching how this project already
  reads other stock-ratio paths (a_nexus in 13d, debt/GDP conventions
  elsewhere), and matching 27's own final choice (level, not differenced)
  for its flow-based exposure measure after the same design question was
  raised there.

  Output: $figs/fig0_evolution_png_gdp_exposure.pdf (Years -3..5,
          country-demeaned mean level path, all onsets pooled / nd / def).
  Run standalone, any time after 18_transforms.do (needs sample, onset_all,
  onset_nd, onset_def, cid, year, iso3).
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* SETUP + BUILD png_gdp (self-contained, same corrected construction as
* 29_png_stock_gdp_exposure.do)
* ══════════════════════════════════════════════════════════════════════════
import excel "$raw/PNGtoGNI.xlsx", sheet("Data") firstrow allstring clear

capture rename CountryCode iso3
capture rename CounterpartAreaName counterpart
capture rename SeriesCode series_code

keep iso3 counterpart series_code YR*
keep if length(iso3) == 3
keep if counterpart == "World"

foreach v of varlist YR* {
    destring `v', replace force
}

tempfile png_raw
save `png_raw'

use `png_raw', clear
keep if series_code == "DT.DOD.DPNG.CD"
reshape long YR, i(iso3) j(year)
rename YR pngdebt_stock
keep iso3 year pngdebt_stock
tempfile png_cy
save `png_cy'

* -- GDP in current US$ (WEO NGDPD, billions) -- same import mechanics as
* 11_weo.do / 29_png_stock_gdp_exposure.do (Excel letter columns AB..CA =
* years 1980..2031).
import excel "$raw/WEOApr2026all.xlsx", sheet("Countries") firstrow clear
capture rename COUNTRYID   iso3
capture rename INDICATORID indid
local col AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP AQ AR AS AT AU ///
          AV AW AX AY AZ BA BB BC BD BE BF BG BH BI BJ BK BL BM BN BO ///
          BP BQ BR BS BT BU BV BW BX BY BZ CA
local y = 1980
foreach c of local col {
    capture rename `c' yr`y'
    local ++y
}
keep iso3 indid yr*
keep if length(iso3) == 3
keep if indid == "NGDPD"
capture destring yr*, replace force
reshape long yr, i(iso3 indid) j(year)
drop if missing(yr)
rename yr gdp_usd_bn
keep iso3 year gdp_usd_bn
label var gdp_usd_bn "GDP, current prices, US$ billions (IMF WEO NGDPD)"
tempfile weo_gdp
save `weo_gdp'

use "$clean/panel_lp.dta", clear
capture drop pngdebt_stock gdp_usd_bn
merge m:1 iso3 year using `png_cy', keep(master match) nogen
merge m:1 iso3 year using `weo_gdp', keep(master match) nogen
sort cid year
xtset cid year

* gdp_usd_bn is in US$ BILLIONS; pngdebt_stock is in raw current US$ --
* convert gdp_usd_bn to raw US$ (x 1e9) before dividing.
capture drop png_gdp
gen double png_gdp = pngdebt_stock / (gdp_usd_bn * 1e9) * 100 ///
    if pngdebt_stock >= 0 & gdp_usd_bn > 0 & !missing(pngdebt_stock, gdp_usd_bn)
label var png_gdp "Private nonguaranteed external debt stock / GDP, pct"

quietly summarize png_gdp if sample==1
di as result _n "  png_gdp coverage/range check (sample==1): N=" r(N) "  mean=" %6.2f r(mean) ///
    "  min=" %6.2f r(min) "  max=" %6.2f r(max)

* -- Coverage diagnostic at onset (of 61 onsets), by resolution type --
capture drop ch_pngexp_0
gen double ch_pngexp_0 = png_gdp
quietly count if onset_all == 1 & sample == 1 & !missing(ch_pngexp_0)
local n_all = r(N)
quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_pngexp_0)
local n_nd  = r(N)
quietly count if onset_def == 1 & sample == 1 & !missing(ch_pngexp_0)
local n_def = r(N)
di as result "  png_gdp coverage at onset: all=" `n_all' " / 61   non-default=" `n_nd' " / 39   default-linked=" `n_def' " / 22"

* ══════════════════════════════════════════════════════════════════════════
* Pre/post-crisis evolution figure -- EXACTLY mirrors 27_png_debt_test.do's
* Part 1 figure (level at each horizon, country-demeaned, rebased so
* Year 0 = last pre-crisis year = 0, Year 1 = crisis onset year,
* matching 02a_descriptive_facts.do's own convention).
* ══════════════════════════════════════════════════════════════════════════
capture drop pngexp_base
gen double pngexp_base = L.png_gdp

capture drop ch_pngexp_m1
gen double ch_pngexp_m1 = L.png_gdp
forvalues k = 2/4 {
    capture drop ch_pngexp_m`k'
    gen double ch_pngexp_m`k' = L`k'.png_gdp
}
forvalues h = 1/4 {
    capture drop ch_pngexp_`h'
    gen double ch_pngexp_`h' = F`h'.png_gdp
}

foreach h in m4 m3 m2 m1 0 1 2 3 4 {
    capture drop cmean_pngexp_`h' dd_pngexp_`h'
    quietly bysort cid: egen double cmean_pngexp_`h' = mean(ch_pngexp_`h') if sample==1
    quietly gen double dd_pngexp_`h' = ch_pngexp_`h' - cmean_pngexp_`h' if sample==1
}

foreach g in all nd def {
    matrix desc_pngexp_`g' = J(9, 1, .)
    quietly summarize dd_pngexp_m1 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[4,1] = r(mean)
    quietly summarize dd_pngexp_m4 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[1,1] = r(mean)
    quietly summarize dd_pngexp_m3 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[2,1] = r(mean)
    quietly summarize dd_pngexp_m2 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[3,1] = r(mean)
    forvalues h = 0/4 {
        quietly summarize dd_pngexp_`h' if onset_`g'==1 & sample==1, meanonly
        matrix desc_pngexp_`g'[`h'+5,1] = r(mean)
    }
    * Rebase so Year 0 (row 4, the LAST PRE-CRISIS year, t-1) = 0 -- matches
    * 02a's/27's own convention exactly.
    scalar _base_`g' = desc_pngexp_`g'[4,1]
    forvalues r = 1/9 {
        matrix desc_pngexp_`g'[`r',1] = desc_pngexp_`g'[`r',1] - _base_`g'
    }
}

preserve
    clear
    svmat desc_pngexp_all, names(b_all)
    svmat desc_pngexp_nd,  names(b_nd)
    svmat desc_pngexp_def, names(b_def)
    gen horizon = _n - 4     // rows are Year -3..5

    local c_nd  "blue"
    local c_def "red"
    twoway ///
        (line b_all1 horizon, lcolor(gs8) lwidth(medthick) lpattern(dash)) ///
        (connected b_nd1  horizon, lcolor("`c_nd'")  mcolor("`c_nd'")  msymbol(circle) lwidth(medthick)) ///
        (connected b_def1 horizon, lcolor("`c_def'") mcolor("`c_def'") msymbol(square) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xline(1, lpattern(dash) lcolor(gs12) lwidth(thin)) ///
        xlabel(-3(1)5, labsize(medsmall)) ///
        ylabel(, format(%9.1f) labsize(medsmall) angle(horizontal)) ///
        xtitle("Year (0 = last pre-crisis year, 1 = crisis onset)", size(small)) ///
        ytitle("PNG debt / GDP, pct, relative to Year 0", size(small)) ///
        title("Private external-debt exposure (PNG debt stock / GDP)", size(medium) color(navy)) ///
        legend(order(1 "All onsets" 2 "Non-default" 3 "Default-linked") ///
               position(6) rows(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig0_evolution_png_gdp_exposure.pdf", replace
    di as result "Figure saved: fig0_evolution_png_gdp_exposure.pdf (Years -3..5, rebased to Year 0 = last pre-crisis year; all/nd/def)"
restore

di as result _n "30_png_gdp_exposure_evolution.do complete. EXPLORATORY/DESCRIPTIVE ONLY -- not"
di as result "wired into 00_master.do. No regression run here -- see 29_png_stock_gdp_exposure.do"
di as result "for the OLS high/low-exposure split on credit/inv/GDP."
