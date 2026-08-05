/*===========================================================================
  16_BANKING.DO   —  FROM-SCRATCH REBUILD, STAGE 6  (banking-sector crises)
  Systemic banking-crisis indicator from Laeven & Valencia (via World Bank GFDD,
  series GFDD.OI.19), same WB DataBank wide layout as the other WDI/GFDD files.
  Mirrors the proven 01b import.

    banking_crisis   -> 0/1 systemic banking-crisis dummy (L-V via GFDD)
  (No duration variable — the dummy is kept exactly as provided, by decision.)

  Output: merges banking_crisis onto $clean/panel_build.dta.
===========================================================================*/

import excel "$raw/Banking_crisis_dummy_data_set.xlsx", sheet("Data") firstrow allstring clear
capture rename SeriesName  series_name
capture rename SeriesCode  series_code
capture rename CountryName country_name
capture rename CountryCode iso3
keep if length(iso3) == 3
capture drop series_name series_code country_name

foreach v of varlist YR* {
    destring `v', replace force
}
reshape long YR, i(iso3) j(year)
rename YR banking_crisis
keep iso3 year banking_crisis
drop if missing(banking_crisis)
keep if year >= 1989
label var banking_crisis "Systemic banking-crisis dummy (Laeven-Valencia via GFDD)"

tempfile bank
save `bank'

* ── Merge onto the growing panel ────────────────────────────────────────────
use "$clean/panel_build.dta", clear
capture drop banking_crisis     // in case a legacy copy exists
merge 1:1 iso3 year using `bank', keep(master match) nogen
save "$clean/panel_build.dta", replace

di as result _n "16_banking.do complete. Coverage:"
quietly count if !missing(banking_crisis) & sample_base==1
di as result "  banking_crisis non-missing (sample rows): `r(N)'"
quietly count if banking_crisis==1 & onset_all==1
di as result "  onsets coinciding with a banking crisis: `r(N)' of 61"
