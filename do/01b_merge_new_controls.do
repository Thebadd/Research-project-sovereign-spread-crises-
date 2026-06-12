/*===========================================================================
  01B_MERGE_NEW_CONTROLS.DO
  Merge three new control variables into panel_lp.dta:
    1. reer_chg      — % change in Real Effective Exchange Rate (WDI, 2010=100)
    2. banking_crisis — systemic banking crisis dummy (Laeven & Valencia via WDI)
    3. revenue_gdp   — general govt revenue excl. grants (% GDP, IMF WEO)

  Source files (all wide format):
    REER_INDEX.xlsx              — WDI, year cols named "YYYY [YRYYYY]"
    Banking_crisis_dummy_data_set.xlsx — WDI, same year col format
    WEOApr2026all.xlsx (Countries sheet) — IMF WEO, year cols named "YYYY"
    revenue_to_GDP.xls           — WDI fallback, year cols named "YYYY",
                                   header on row 4 (skip rows 1–3)

  Merge key: country (string name) × year
  ISO → panel name crosswalk defined below.
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* 0. BUILD ISO → PANEL NAME CROSSWALK
* ══════════════════════════════════════════════════════════════════════════

clear
input str3 iso3 str30 country
"ARG" "Argentina"
"ARM" "Armenia"
"AZE" "Azerbaijan"
"BLZ" "Belize"
"BOL" "Bolivia"
"BRA" "Brazil"
"BGR" "Bulgaria"
"CHL" "Chile"
"CHN" "China"
"COL" "Colombia"
"CRI" "Costa Rica"
"CIV" "Cote d'Ivoire"
"DOM" "Dominican Republic"
"ECU" "Ecuador"
"EGY" "Egypt"
"SLV" "El Salvador"
"GHA" "Ghana"
"GTM" "Guatemala"
"HND" "Honduras"
"HUN" "Hungary"
"IND" "India"
"IDN" "Indonesia"
"JAM" "Jamaica"
"JOR" "Jordan"
"KAZ" "Kazakhstan"
"KEN" "Kenya"
"LBN" "Lebanon"
"MYS" "Malaysia"
"MEX" "Mexico"
"MAR" "Morocco"
"NAM" "Namibia"
"NGA" "Nigeria"
"PAK" "Pakistan"
"PAN" "Panama"
"PRY" "Paraguay"
"PER" "Peru"
"PHL" "Philippines"
"POL" "Poland"
"ROU" "Romania"
"RUS" "Russia"
"SEN" "Senegal"
"SRB" "Serbia"
"ZAF" "South Africa"
"LKA" "Sri Lanka"
"TTO" "Trinidad and Tobago"
"TUN" "Tunisia"
"TUR" "Turkey"
"UKR" "Ukraine"
"URY" "Uruguay"
"VEN" "Venezuela"
"VNM" "Vietnam"
"ZMB" "Zambia"
end

tempfile crosswalk
save `crosswalk'

* ══════════════════════════════════════════════════════════════════════════
* 1. REER INDEX → reer_chg
*    Year columns: "YYYY [YRYYYY]" → Stata renames to v# after firstrow import
*    Use cellrange to import only data rows; rename year cols manually.
* ══════════════════════════════════════════════════════════════════════════

import excel "$raw/REER_INDEX.xlsx", sheet("Data") firstrow allstring clear

* After firstrow import, year columns come in as: A1990YR1990, A1991YR1991...
* or Stata may name them with the bracket notation truncated.
* Rename the four non-year columns first, then reshape on everything else.
capture rename SeriesName    series_name
capture rename SeriesCode    series_code
capture rename CountryName   country_name
capture rename CountryCode   iso3

* Drop if not a 3-letter ISO code (removes regional aggregates)
keep if length(iso3) == 3

* Drop non-data columns
capture drop series_name series_code country_name

* Destring all remaining columns (year data)
foreach v of varlist * {
    if "`v'" != "iso3" {
        destring `v', replace force
    }
}

* Reshape to long — all non-iso3 variables are year columns
reshape long , i(iso3) j(year_str) string
rename _j year_col
rename __000000 reer_index

* year_str will be the column name; extract 4-digit year from it
gen year = real(substr(year_col, 1, 4))
drop year_col year_str
keep if !missing(year) & year >= 1989

* Sort and compute % change
sort iso3 year
by iso3: gen reer_chg = (reer_index / reer_index[_n-1] - 1) * 100 if reer_index[_n-1] != .
drop reer_index

* Merge with crosswalk
merge m:1 iso3 using `crosswalk', keep(match) nogen
keep country year reer_chg

sort country year
tempfile reer
save `reer'
di as result "REER: " _N " country-year observations"

* ══════════════════════════════════════════════════════════════════════════
* 2. BANKING CRISIS DUMMY
*    Same WDI wide format as REER file
* ══════════════════════════════════════════════════════════════════════════

import excel "$raw/Banking_crisis_dummy_data_set.xlsx", sheet("Data") firstrow allstring clear

capture rename SeriesName  series_name
capture rename SeriesCode  series_code
capture rename CountryName country_name
capture rename CountryCode iso3

keep if length(iso3) == 3
capture drop series_name series_code country_name

foreach v of varlist * {
    if "`v'" != "iso3" {
        destring `v', replace force
    }
}

reshape long , i(iso3) j(year_str) string
rename _j year_col
rename __000000 banking_crisis

gen year = real(substr(year_col, 1, 4))
drop year_col year_str
keep if !missing(year) & year >= 1989

keep if !missing(banking_crisis)

merge m:1 iso3 using `crosswalk', keep(match) nogen
keep country year banking_crisis

sort country year
tempfile bcrisis
save `bcrisis'
di as result "Banking crisis: " _N " country-year observations"

* ══════════════════════════════════════════════════════════════════════════
* 3. GOVERNMENT REVENUE / GDP — IMF WEO
*    Sheet "Countries", INDICATOR.ID == GGR_NGDP
*    Year columns named as plain integers: 1980, 1981, ..., 2031
* ══════════════════════════════════════════════════════════════════════════

import excel "$raw/WEOApr2026all.xlsx", sheet("Countries") firstrow clear

* WEO column COUNTRY.ID → iso3; dots in names → Stata uses COUNTRYID
capture rename COUNTRYID iso3
capture rename v1        iso3    // fallback

* Keep only revenue indicator
keep if INDICATORID == "GGR_NGDP"

keep iso3 v*
* Year columns are numeric-named: Stata stores them as v5, v6 ... or as the
* actual numbers. After firstrow they become variable names = column headers.
* Rename to yr#### pattern for reshape.

* Identify and rename year columns (1980–2031)
foreach y of numlist 1980/2031 {
    capture rename `y' yr`y'
}

keep iso3 yr*
keep if length(iso3) == 3

reshape long yr, i(iso3) j(year)
rename yr revenue_gdp

drop if missing(revenue_gdp)
keep if year >= 1989

merge m:1 iso3 using `crosswalk', keep(match) nogen
keep country year revenue_gdp

sort country year
tempfile rev_weo
save `rev_weo'
di as result "WEO Revenue: " _N " country-year observations"

* ══════════════════════════════════════════════════════════════════════════
* 4. MERGE ALL NEW CONTROLS INTO PANEL_LP.DTA
* ══════════════════════════════════════════════════════════════════════════

use "$clean/panel_lp.dta", clear

* -- REER change --
merge m:1 country year using `reer', keep(master match) nogen
quietly count if !missing(reer_chg)
di as result "reer_chg: " r(N) " non-missing obs merged"

* -- Banking crisis --
merge m:1 country year using `bcrisis', keep(master match) nogen
quietly count if !missing(banking_crisis)
di as result "banking_crisis: " r(N) " non-missing obs merged"

* -- Revenue / GDP --
merge m:1 country year using `rev_weo', keep(master match) nogen
quietly count if !missing(revenue_gdp)
di as result "revenue_gdp: " r(N) " non-missing obs merged"

* -- Labels --
label var reer_chg       "% change in REER (WDI, 2010=100 index)"
label var banking_crisis "Systemic banking crisis dummy (Laeven-Valencia/WDI)"
label var revenue_gdp    "Govt revenue excl. grants (% GDP, IMF WEO GGR_NGDP)"

* ── Coverage report at onset years ────────────────────────────────────────
di as result _n "=== NEW VARIABLE COVERAGE AT ONSET (onset_all==1, sample==1) ==="
foreach v in reer_chg banking_crisis revenue_gdp {
    quietly count if onset_all == 1 & sample == 1 & !missing(`v')
    di as result "  `v': " r(N) " / 61 onsets with non-missing data"
}

di as result _n "=== NEW VARIABLE COVERAGE: ALL SAMPLE YEARS (sample==1) ==="
foreach v in reer_chg banking_crisis revenue_gdp {
    quietly count if sample == 1 & !missing(`v')
    di as result "  `v': " r(N) " obs"
}

save "$clean/panel_lp.dta", replace
di as result _n "01b_merge_new_controls.do complete — panel_lp.dta updated."
