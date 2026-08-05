/*===========================================================================
  15_RATES.DO   —  FROM-SCRATCH REBUILD, STAGE 5  (global-push predictors, FRED)
  US/global financial-conditions series from FRED. These are pure time series
  (identical across countries), merged on YEAR. In the LP they are absorbed by
  year FE (exclusion restriction); they enter only the first-stage propensity
  model.

  Expected raw files in $raw (annual, FRED "observation_date,VALUE" layout):
    GS10.csv     -> ust10y    US 10-year Treasury yield (%)          [PROVIDED]
    FEDFUNDS.csv -> fedfunds  US federal funds rate (%)              [NEEDED]
    VIXCLS.csv   -> vix        CBOE VIX, annual average               [NEEDED]

  Output: merges ust10y / fedfunds / vix onto $clean/panel_build.dta (by year).
===========================================================================*/

tempfile rates
local have = 0

* filename  ->  target variable
foreach pair in "GS10 ust10y" "FEDFUNDS fedfunds" "VIXCLS vix" {
    gettoken fn tv : pair
    capture confirm file "$raw/`fn'.csv"
    if _rc {
        di as error "  ** `fn'.csv not found in $raw — `tv' will be missing. Add the FRED annual export."
        continue
    }
    preserve
        import delimited "$raw/`fn'.csv", varnames(1) clear
        * FRED col 1 = observation_date (YYYY-...); col 2 = the series
        capture confirm variable observation_date
        if _rc {
            * fallback: first var is the date
            ds
            local first : word 1 of `r(varlist)'
            rename `first' observation_date
        }
        gen int year = real(substr(observation_date,1,4))
        * the value column is the only non-date variable
        local val
        foreach v of varlist * {
            if "`v'"!="observation_date" & "`v'"!="year" local val "`val' `v'"
        }
        local val : word 1 of `val'
        capture destring `val', replace force
        rename `val' `tv'
        keep year `tv'
        drop if missing(year)
        collapse (mean) `tv', by(year)     // safety if any sub-annual rows
        if `have'==0 {
            save `rates', replace
            local have = 1
        }
        else {
            merge 1:1 year using `rates', nogen
            save `rates', replace
        }
    restore
}

if `have'==0 {
    di as error "  ** No FRED files found — 15_rates.do made no change."
    exit
}

label var ust10y   "US 10y Treasury yield, % (FRED GS10)"
capture label var fedfunds "US fed funds rate, % (FRED FEDFUNDS)"
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
