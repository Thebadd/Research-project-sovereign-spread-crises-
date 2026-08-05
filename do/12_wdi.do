/*===========================================================================
  12_WDI.DO   —  FROM-SCRATCH REBUILD, STAGE 3  (World Bank WDI)
  Private credit, FDI, claims on central government, and trade openness
  (exports + imports) from World Bank WDI. Standard WDI wide layout:
    Country Name | Country Code(ISO3) | Indicator Name | Indicator Code | 1960..2025
  One indicator per file; first sheet (sheet names vary / have trailing spaces,
  so we do NOT name the sheet — import defaults to the first). Year columns are
  renamed POSITIONALLY (metadata first, then years ascending) so we don't depend
  on how Stata names numeric headers.

  Expected files in $raw (name them exactly as below) -> variable [WDI code]:
    domesticcredittoprivatesector.xlsx -> credit   [FS.AST.PRVT.GD.ZS]  (CHANNEL)
    domesticcreditprivatebanks.xlsx    -> credit_bank [FD.AST.PRVT.GD.ZS] (by banks, robustness)
    fdinetinflowsgdp.xlsx              -> fdi      [BX.KLT.DINV.WD.GD.ZS]
    claimsoncentralgovernmentgdp.xlsx  -> claims_govt [FS.AST.CGOV.GD.ZS]
    exportofgoodservicesgdp.xlsx       -> exp_gdp  [NE.EXP.GNFS.ZS]
    importofgoodservicesgdp.xlsx       -> imp_gdp  [NE.IMP.GNFS.ZS]
  Then open = exp_gdp + imp_gdp (trade openness, % GDP; Asonuma control).

  NOTE: the channel variable `credit` = domestic credit to private sector, all
  financial corporations (FS.AST.PRVT.GD.ZS). The by-banks measure
  (FD.AST.PRVT.GD.ZS) is kept as `credit_bank` for robustness / the bank-
  intermediation comparison.

  Output: merges credit, fdi, claims_govt, exp_gdp, imp_gdp, open onto
          $clean/panel_build.dta (iso3 x year).
===========================================================================*/

tempfile wdi
local have = 0

foreach spec in "domesticcredittoprivatesector credit" ///
                "domesticcreditprivatebanks credit_bank" ///
                "fdinetinflowsgdp fdi" ///
                "claimsoncentralgovernmentgdp claims_govt" ///
                "exportofgoodservicesgdp exp_gdp" ///
                "importofgoodservicesgdp imp_gdp" {
    gettoken fn tv : spec
    capture confirm file "$raw/`fn'.xlsx"
    if _rc {
        di as error "  ** `fn'.xlsx not found in $raw — `tv' will be missing."
        continue
    }
    preserve
        import excel "$raw/`fn'.xlsx", firstrow allstring clear
        capture rename CountryCode iso3
        * year columns = everything except the 4 text metadata columns; they are
        * laid out in ascending year order -> rename positionally to yr1960..
        local yrvars
        foreach v of varlist * {
            if !inlist("`v'","CountryName","iso3","IndicatorName","IndicatorCode") ///
                local yrvars `yrvars' `v'
        }
        local y = 1960
        foreach v of local yrvars {
            capture rename `v' yr`y'
            local ++y
        }
        keep iso3 yr*
        keep if length(iso3) == 3
        foreach v of varlist yr* {
            destring `v', replace force      // ".." -> missing
        }
        reshape long yr, i(iso3) j(year)
        rename yr `tv'
        drop if missing(`tv')
        keep if year >= 1989
        if `have' == 0 {
            save `wdi', replace
            local have = 1
        }
        else {
            merge 1:1 iso3 year using `wdi', nogen
            save `wdi', replace
        }
    restore
}

if `have' == 0 {
    di as error "  ** No WDI files found — 12_wdi.do made no change."
    exit
}

* ── REER (World Bank, DataBank "YYYY [YRYYYY]" wide layout — mirrors 01b) ────
*    REER_INDEX.xlsx -> reer_index -> reer_chg (% YoY change). Merged into the
*    same accumulator (iso3 x year). Skipped if the file is absent.
capture confirm file "$raw/REER_INDEX.xlsx"
if !_rc {
    preserve
        import excel "$raw/REER_INDEX.xlsx", sheet("Data") firstrow allstring clear
        capture rename CountryCode iso3
        keep if length(iso3) == 3
        capture drop SeriesName SeriesCode CountryName
        foreach v of varlist YR* {
            destring `v', replace force
        }
        reshape long YR, i(iso3) j(year)
        rename YR reer_index
        keep if !missing(year) & year >= 1989
        sort iso3 year
        by iso3: gen double reer_chg = (reer_index/reer_index[_n-1] - 1)*100 ///
            if !missing(reer_index[_n-1])
        keep iso3 year reer_chg
        drop if missing(reer_chg)
        merge 1:1 iso3 year using `wdi', nogen
        save `wdi', replace
    restore
}

* ── Merge the WDI accumulator onto the growing panel ────────────────────────
use "$clean/panel_build.dta", clear
merge 1:1 iso3 year using `wdi', keep(master match) nogen

* trade openness = exports + imports (% GDP)
capture confirm variable exp_gdp
local hasexp = (_rc==0)
capture confirm variable imp_gdp
local hasimp = (_rc==0)
if `hasexp' & `hasimp' {
    capture drop open
    gen double open = exp_gdp + imp_gdp
    label var open "Trade openness = exports+imports, % GDP (WDI NE.EXP+NE.IMP)"
}
capture label var credit      "Domestic credit to private sector, % GDP (WDI FS.AST.PRVT.GD.ZS)"
capture label var credit_bank "Domestic credit to private sector by banks, % GDP (WDI FD.AST.PRVT.GD.ZS)"
capture label var fdi          "FDI net inflows, % GDP (WDI BX.KLT.DINV.WD.GD.ZS)"
capture label var claims_govt  "Claims on central government, % GDP (WDI FS.AST.CGOV.GD.ZS)"
capture label var exp_gdp      "Exports of goods & services, % GDP (WDI NE.EXP.GNFS.ZS)"
capture label var imp_gdp      "Imports of goods & services, % GDP (WDI NE.IMP.GNFS.ZS)"
capture label var reer_chg     "REER, % YoY change (WDI, 2010=100 index)"

save "$clean/panel_build.dta", replace

di as result _n "12_wdi.do complete. Coverage at onsets (non-missing):"
foreach v in credit fdi claims_govt open {
    capture confirm variable `v'
    if !_rc {
        quietly count if onset_all==1 & !missing(`v')
        local a=r(N)
        quietly count if onset_nd==1 & !missing(`v')
        local n=r(N)
        quietly count if onset_def==1 & !missing(`v')
        local d=r(N)
        di as result "  `v': all=`a'  nd=`n'  def=`d'"
    }
}
