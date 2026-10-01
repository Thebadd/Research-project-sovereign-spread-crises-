/*===========================================================================
  15_RATES.DO   —  FROM-SCRATCH REBUILD, STAGE 5  (global-push predictors, FRED)
  US/global financial-conditions series from FRED. Pure time series (identical
  across countries), merged on YEAR. Absorbed by year FE in the LP (exclusion
  restriction); they enter only the first-stage propensity model.

  Expected raw files in $raw (annual; the loader accepts .csv OR .xlsx):
    GS10      -> ust10y    US 10-year Treasury yield (%)          [.csv provided]
    FEDFUNDS  -> fedfunds  US federal funds effective rate (%)    [.xlsx, "Annual"]
    VIXCLS    -> vix        CBOE VIX, annual average               [.xlsx, "Annual"]

  FRED .csv = "observation_date,VALUE" (string date);
  FRED .xlsx = README + "Annual" sheet with observation_date (datetime) + VALUE.

  Output: merges ust10y / fedfunds / vix onto $clean/panel_build.dta (by year).
===========================================================================*/

tempfile rates
local have = 0

* filename stem  ->  target variable
foreach pair in "GS10 ust10y" "FEDFUNDS fedfunds" "VIXCLS vix" {
    gettoken fn tv : pair

    * locate the file (csv preferred, else xlsx)
    local src ""
    capture confirm file "$raw/`fn'.csv"
    if !_rc local src "csv"
    else {
        capture confirm file "$raw/`fn'.xlsx"
        if !_rc local src "xlsx"
    }
    if "`src'" == "" {
        di as error "  ** `fn'(.csv/.xlsx) not found in $raw — `tv' will be missing."
        continue
    }

    preserve
        if "`src'" == "csv"  import delimited "$raw/`fn'.csv", varnames(1) clear
        else                 import excel "$raw/`fn'.xlsx", sheet("Annual") firstrow clear

        * standardise the date column name
        capture confirm variable observation_date
        if _rc {
            ds
            local first : word 1 of `r(varlist)'
            rename `first' observation_date
        }

        * year: handle string dates (csv) vs Stata date/datetime (xlsx import)
        capture confirm numeric variable observation_date
        if !_rc {
            local fmt : format observation_date
            if strpos("`fmt'", "%tc") > 0  gen int year = year(dofc(observation_date))
            else                            gen int year = year(observation_date)
        }
        else gen int year = real(substr(observation_date, 1, 4))

        * value column = the one non-date variable
        local val
        foreach v of varlist * {
            if !inlist("`v'", "observation_date", "year") local val "`val' `v'"
        }
        local val : word 1 of `val'
        capture destring `val', replace force
        rename `val' `tv'

        keep year `tv'
        drop if missing(year)
        collapse (mean) `tv', by(year)      // safety if any duplicate/sub-annual rows

        if `have' == 0 {
            save `rates', replace
            local have = 1
        }
        else {
            merge 1:1 year using `rates', nogen
            save `rates', replace
        }
    restore
}

if `have' == 0 {
    di as error "  ** No FRED files found — 15_rates.do made no change."
    exit
}

capture label var ust10y   "US 10y Treasury yield, % (FRED GS10)"
capture label var fedfunds "US fed funds effective rate, % (FRED FEDFUNDS)"
capture label var vix      "CBOE VIX, annual avg (FRED VIXCLS)"

use "$clean/panel_build.dta", clear
merge m:1 year using `rates', keep(master match) nogen
save "$clean/panel_build.dta", replace

di as result _n "15_rates.do complete. Coverage over sample years:"
foreach v in ust10y fedfunds vix {
    capture confirm variable `v'
    if !_rc {
        quietly count if !missing(`v') & sample_base==1
        di as result "  `v': `r(N)' non-missing sample rows"
    }
    else di as result "  `v': (not built — file missing)"
}
