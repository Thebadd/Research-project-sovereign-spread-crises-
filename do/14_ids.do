/*===========================================================================
  14_IDS.DO   —  FROM-SCRATCH REBUILD, STAGE 4b  (external-debt vulnerability)
  World Bank International Debt Statistics (IDS), DataBank "YYYY [YRYYYY]" wide
  layout. Series identified by series_code; keyed on ISO3. Mirrors 01d.

    DT.DOD.DSTC.ZS   -> stdebt_share      Short-term debt, % total ext. debt (rollover)
    FI.RES.TOTL.DT.ZS-> reserves_extdebt  Total reserves, % total ext. debt (buffer)
    DT.INT.DECT.GN.ZS-> intpay_gni        Interest on ext. debt, % GNI (interest burden)
    DT.TDS.DECT.EX.ZS-> debt_service      Total debt service, % of exports (liquidity)

  Files in $raw: International_Debt_Statistics.xlsx (first three series) and
  Total_debt_service.xlsx (debt_service). Used by 13b_exposure_heterogeneity.
===========================================================================*/

tempfile ids
local have = 0

* helper: pull one IDS series from a given file into iso3-year long form
capture program drop _idsseries
program define _idsseries
    syntax, FILE(string) CODE(string) NEWNAME(string) SAVING(string)
    import excel "`file'", sheet("Data") firstrow allstring clear
    capture rename CountryCode iso3
    capture rename SeriesCode  series_code
    keep iso3 series_code YR*
    keep if length(iso3) == 3
    keep if series_code == "`code'"
    foreach v of varlist YR* {
        destring `v', replace force
    }
    reshape long YR, i(iso3) j(year)
    rename YR `newname'
    keep iso3 year `newname'
    drop if missing(`newname')
    save "`saving'", replace
end

foreach spec in ///
    "International_Debt_Statistics DT.DOD.DSTC.ZS stdebt_share" ///
    "International_Debt_Statistics FI.RES.TOTL.DT.ZS reserves_extdebt" ///
    "International_Debt_Statistics DT.INT.DECT.GN.ZS intpay_gni" ///
    "Total_debt_service DT.TDS.DECT.EX.ZS debt_service" {
    gettoken fn spec : spec
    gettoken cd spec : spec
    gettoken nm spec : spec
    capture confirm file "$raw/`fn'.xlsx"
    if _rc {
        di as error "  ** `fn'.xlsx not found in $raw — `nm' will be missing."
        continue
    }
    preserve
        tempfile one
        _idsseries, file("$raw/`fn'.xlsx") code("`cd'") newname("`nm'") saving("`one'")
        if `have' == 0 {
            copy "`one'" "`ids'", replace
            local have = 1
        }
        else {
            use "`one'", clear
            merge 1:1 iso3 year using `ids', nogen
            save `ids', replace
        }
    restore
}

if `have' == 0 {
    di as error "  ** No IDS files found — 14_ids.do made no change."
    exit
}

capture label var stdebt_share     "Short-term debt, % total external debt — rollover (WB IDS)"
capture label var reserves_extdebt "Total reserves, % total external debt — buffer (WB IDS)"
capture label var intpay_gni       "Interest on ext. debt, % GNI — interest burden (WB IDS)"
capture label var debt_service     "Total debt service, % of exports — liquidity (WB IDS)"

use "$clean/panel_build.dta", clear
merge 1:1 iso3 year using `ids', keep(master match) nogen
save "$clean/panel_build.dta", replace

di as result _n "14_ids.do complete. Coverage at onsets:"
foreach v in stdebt_share reserves_extdebt intpay_gni debt_service {
    capture confirm variable `v'
    if !_rc {
        quietly count if onset_all==1 & !missing(`v')
        di as result "  `v': `r(N)' / 61 onsets"
    }
}
