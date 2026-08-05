/*===========================================================================
  13_IFS_NEXUS.DO   —  FROM-SCRATCH REBUILD, STAGE 4  (sovereign-bank nexus)
  Bank balance-sheet composition from IMF Monetary & Financial Statistics
  (Other Depository Corporations), data/raw/Nexus_Sovereign_Bank.xlsx. Same
  proven logic as the old 01c; merged into panel_build on country x year (the
  IMF names are harmonised to the panel's short names, which the skeleton keeps).

    claimsgov_assets    = Total claims on government / Total assets   (%)
    claimpriv_assets    = Claims on private sector  / Total assets    (%)   [BANK]
    netclaimsgov_assets = Net claims on government  / Total assets    (%)

  Source: IMF MFS, aggregate ODC (deposit-taking banks ex central bank) balance
  sheet, local-currency levels, 2001-2024. All three are shares of TOTAL BANK
  ASSETS, so claimpriv_assets is the bank-only private-credit measure that is
  internally consistent with the by-banks `credit` (FD.AST.PRVT.GD.ZS) in 12.

  Coverage: 2001+ only; 6 panel countries absent from the IMF file and stay
  missing: China (mainland), Ecuador, El Salvador, India, Lebanon, Vietnam.
===========================================================================*/

* ── Harmonise IMF country names to the panel's short names ────────────────
capture program drop _nexusnames
program define _nexusnames
    replace country = "Armenia"            if strpos(country, "Armenia")
    replace country = "Azerbaijan"         if strpos(country, "Azerbaijan")
    replace country = "Cote d'Ivoire"      if strpos(country, "Ivoire")
    replace country = "Dominican Republic" if strpos(country, "Dominican")
    replace country = "Egypt"              if strpos(country, "Egypt")
    replace country = "Kazakhstan"         if strpos(country, "Kazakhstan")
    replace country = "Poland"             if strpos(country, "Poland")
    replace country = "Russia"             if strpos(country, "Russian")
    replace country = "Serbia"             if strpos(country, "Serbia")
    replace country = "Turkey"             if strpos(country, "rkiye")
    replace country = "Venezuela"          if strpos(country, "Venezuela")
end

local nexus "$raw/Nexus_Sovereign_Bank.xlsx"

* 1. Total claims on govt / assets (already a ratio in the workbook)
import excel "`nexus'", sheet("Tot_claim_gov_pct") cellrange(A2) firstrow clear
_nexusnames
destring Totclaimsgov_assetpct*, replace force
reshape long Totclaimsgov_assetpct, i(country) j(year)
gen claimsgov_assets = 100 * Totclaimsgov_assetpct
drop Totclaimsgov_assetpct
tempfile tg
save `tg'

* 2. Net claims on govt / assets
import excel "`nexus'", sheet("Net_claim_gov_pct") cellrange(A2) firstrow clear
_nexusnames
destring netclaimsgov_assetpct*, replace force
reshape long netclaimsgov_assetpct, i(country) j(year)
gen netclaimsgov_assets = 100 * netclaimsgov_assetpct
drop netclaimsgov_assetpct
tempfile tn
save `tn'

* 3. Claims on private sector (level, LCU)
import excel "`nexus'", sheet("Claim_priv") cellrange(A2) firstrow clear
_nexusnames
destring claim_priv*, replace force
reshape long claim_priv, i(country) j(year)
tempfile tp
save `tp'

* 4. Total assets (level, LCU) from the stacked Source sheet
import excel "`nexus'", sheet("Source") cellrange(A2) clear
drop if A == "country"
keep if B == "Assets"
rename A country
drop B
ds country, not
local y = 2001
foreach v of varlist `r(varlist)' {
    rename `v' assets`y'
    local ++y
}
_nexusnames
destring assets*, replace force
reshape long assets, i(country) j(year)
tempfile ta
save `ta'

* 5. Assemble nexus country-year file
use `tp', clear
merge 1:1 country year using `ta', keep(match) nogen
gen claimpriv_assets = 100 * claim_priv / assets if assets > 0 & !missing(assets)
drop claim_priv assets
merge 1:1 country year using `tg', nogen
merge 1:1 country year using `tn', nogen
label var claimsgov_assets    "Bank claims on govt / assets, % (IMF MFS/ODC)"
label var claimpriv_assets    "Bank claims on private / assets, % (IMF MFS/ODC)"
label var netclaimsgov_assets "Bank net claims on govt / assets, % (IMF MFS/ODC)"
keep country year claimsgov_assets claimpriv_assets netclaimsgov_assets
tempfile nexus_cy
save `nexus_cy'

* 6. Merge into the growing panel (on country x year — skeleton keeps short names)
use "$clean/panel_build.dta", clear
capture drop claimsgov_assets claimpriv_assets netclaimsgov_assets
merge m:1 country year using `nexus_cy', keep(master match) nogen
save "$clean/panel_build.dta", replace

di as result _n "13_ifs_nexus.do complete. Nexus coverage at onsets:"
foreach v in claimsgov_assets claimpriv_assets netclaimsgov_assets {
    quietly count if onset_all==1 & !missing(`v')
    di as result "  `v': `r(N)' / 61 onsets"
}
