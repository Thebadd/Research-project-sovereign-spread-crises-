/*===========================================================================
  11_WEO.DO   —  FROM-SCRATCH REBUILD, STAGE 2  (IMF WEO, primary macro source)
  Extracts the national-accounts + fiscal aggregates from the IMF World Economic
  Outlook export already on disk (WEOApr2026all.xlsx, sheet "Countries"), keyed
  on ISO3. Mirrors the proven import pattern in 01b (year columns land as Excel
  letters AB..CA = 1980..2031).

  WEO indicator (INDICATOR.ID) -> panel variable:
    NGDPRPC      -> gdppc_real   Real GDP per capita, constant prices, dom. cur.
    NGDP_R       -> gdp_real     Real GDP, constant prices (level)
    NGDP         -> gdp_nominal  Nominal GDP, current prices (for ratios)
    LP           -> pop          Population (millions)
    PCPIPCH      -> infl         CPI inflation, period avg, % change
    BCA_NGDPD    -> ca           Current account balance, % of GDP
    GGXWDG_NGDP  -> debt         General govt gross debt, % of GDP
    GGX_NGDP     -> govexp       General govt total expenditure, % of GDP
    GGR_NGDP     -> revenue_gdp  General govt revenue, % of GDP
    NID_NGDP     -> inv          Total investment, % of GDP
    GGXONLB_NGDP -> pb           Primary net lending/borrowing, % of GDP

  Growth rates, per-capita logs beyond ln_gdppc, lags and pre_* are built in 18.
  Output: merges the above onto $clean/panel_build.dta (iso3 x year).
===========================================================================*/

* ── 1. Import WEO, standardise keys ─────────────────────────────────────────
import excel "$raw/WEOApr2026all.xlsx", sheet("Countries") firstrow clear
capture rename COUNTRYID   iso3          // "COUNTRY.ID"  -> COUNTRYID
capture rename INDICATORID indid         // "INDICATOR.ID"-> INDICATORID

* Year columns imported as Excel letters AB..CA -> rename to yr1980..yr2031
local col AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP AQ AR AS AT AU ///
          AV AW AX AY AZ BA BB BC BD BE BF BG BH BI BJ BK BL BM BN BO ///
          BP BQ BR BS BT BU BV BW BX BY BZ CA
local y = 1980
foreach c of local col {
    capture rename `c' yr`y'
    local ++y
}

* ── 1b. Last OUTTURN year for the LP outcome series ─────────────────────────
* WEO Apr-2026 runs to 2031, so anything past the last actual observation is an
* IMF PROJECTION, not data. That matters here specifically: an onset in 2022+
* has its h=3/h=4 outcome (GDP at t+3, t+4) built partly from forecasts, and those
* are exactly the horizons carrying the default premium. Column N
* (LATEST_ACTUAL_ANNUAL_DATA) gives the boundary per country and indicator; we
* take it from NGDPRPC, the series the LP outcome is built from.
* Values are mostly plain years but some are fiscal-year strings ("FY2024/25"),
* so pull the leading 4-digit year rather than destringing.
preserve
    capture confirm variable LATEST_ACTUAL_ANNUAL_DATA
    if _rc == 0 {
        keep if indid == "NGDPRPC"
        capture confirm string variable LATEST_ACTUAL_ANNUAL_DATA
        if _rc capture tostring LATEST_ACTUAL_ANNUAL_DATA, replace force
        gen int gdp_last_actual = .
        replace gdp_last_actual = real(regexs(0)) ///
            if regexm(LATEST_ACTUAL_ANNUAL_DATA, "(19|20)[0-9][0-9]")
        keep iso3 gdp_last_actual
        keep if length(iso3) == 3
        drop if missing(gdp_last_actual)
        bysort iso3: keep if _n == 1
        label var gdp_last_actual "Last WEO OUTTURN year for real GDP p.c. (later years are projections)"
        tempfile lastactual
        save `lastactual'
        global _weo_lastactual 1
    }
    else {
        di as error "  ** 11_weo: LATEST_ACTUAL_ANNUAL_DATA column not found —"
        di as error "     forecast-contamination robustness cuts will be skipped."
        global _weo_lastactual 0
    }
restore

* ── 2. Keep the needed indicators + year columns ────────────────────────────
keep iso3 indid yr*
keep if length(iso3) == 3
keep if inlist(indid, "NGDPRPC","NGDP_R","NGDP","LP","PCPIPCH","BCA_NGDPD") ///
      | inlist(indid, "GGXWDG_NGDP","GGX_NGDP","GGR_NGDP","NID_NGDP","GGXONLB_NGDP")

* WEO uses "n/a" / "--" for missing -> force to numeric
capture destring yr*, replace force

* ── 3. Long by year, then wide by indicator (one column per series) ─────────
reshape long yr, i(iso3 indid) j(year)
drop if missing(yr)
reshape wide yr, i(iso3 year) j(indid) string
* columns now yr<CODE>; rename to friendly names
capture rename yrNGDPRPC      gdppc_real
capture rename yrNGDP_R       gdp_real
capture rename yrNGDP         gdp_nominal
capture rename yrLP           pop
capture rename yrPCPIPCH      infl
capture rename yrBCA_NGDPD    ca
capture rename yrGGXWDG_NGDP  debt
capture rename yrGGX_NGDP     govexp
capture rename yrGGR_NGDP     revenue_gdp
capture rename yrNID_NGDP     inv
capture rename yrGGXONLB_NGDP pb

* ── 4. Log real GDP per capita (LP-outcome base; growth/lags built in 18) ────
gen double ln_gdppc = ln(gdppc_real) if gdppc_real > 0 & !missing(gdppc_real)

label var gdppc_real  "Real GDP per capita, const. prices (IMF WEO NGDPRPC)"
label var gdp_real    "Real GDP, constant prices (IMF WEO NGDP_R)"
label var gdp_nominal "Nominal GDP (IMF WEO NGDP)"
label var pop         "Population, millions (IMF WEO LP)"
label var ln_gdppc    "Log real GDP per capita"
label var infl        "CPI inflation, % (IMF WEO PCPIPCH)"
label var ca          "Current account, % GDP (IMF WEO BCA_NGDPD)"
label var debt        "Govt gross debt, % GDP (IMF WEO GGXWDG_NGDP)"
label var govexp      "Govt expenditure, % GDP (IMF WEO GGX_NGDP)"
label var revenue_gdp "Govt revenue, % GDP (IMF WEO GGR_NGDP)"
label var inv         "Total investment, % GDP (IMF WEO NID_NGDP)"
label var pb          "Primary balance, % GDP (IMF WEO GGXONLB_NGDP)"

keep iso3 year gdppc_real gdp_real gdp_nominal pop ln_gdppc ///
     infl ca debt govexp revenue_gdp inv pb
tempfile weo
save `weo'

* ── 5. Merge onto the growing panel (skeleton drives the sample) ────────────
use "$clean/panel_build.dta", clear
merge 1:1 iso3 year using `weo', keep(master match) nogen
if "$_weo_lastactual" == "1" {
    capture drop gdp_last_actual
    merge m:1 iso3 using `lastactual', keep(master match) nogen
}
save "$clean/panel_build.dta", replace

* ── 6. Coverage report (onset observations with non-missing WEO macro) ──────
di as result _n "11_weo.do complete. Coverage at onsets (non-missing):"
foreach v in gdppc_real infl ca debt govexp revenue_gdp inv pb {
    quietly count if onset_all==1 & !missing(`v')
    local na = r(N)
    quietly count if onset_nd==1 & !missing(`v')
    local nn = r(N)
    quietly count if onset_def==1 & !missing(`v')
    local nd = r(N)
    di as result "  `v': all=`na'  nd=`nn'  def=`nd'"
}

* Forecast-contamination scale: how many onsets have their h=4 outcome built from
* projections rather than outturns. This is the caveat that hangs over the long
* horizons, so quantify it every run rather than asserting it.
capture confirm variable gdp_last_actual
if _rc == 0 {
    quietly summarize gdp_last_actual, detail
    di as result _n "  last WEO outturn year: median " r(p50) ", min " r(min) ", max " r(max)
    forvalues h = 0/4 {
        quietly count if onset_all==1 & (year + `h') > gdp_last_actual & !missing(gdp_last_actual)
        di as result "  onsets whose h=`h' outcome falls in the projection window: " r(N)
    }
}
