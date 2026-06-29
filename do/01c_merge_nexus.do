/*===========================================================================
  01C_MERGE_NEXUS.DO
  Sovereign-bank nexus variables from IMF Monetary & Financial Statistics
  Source file: data/raw/Nexus_Sovereign_Bank.xlsx
               (aggregate balance sheet of Other Depository Corporations,
                148 reporting economies, 2001-2024, local-currency levels)

  Builds three banking-sector balance-sheet ratios, all scaled by total
  bank assets so they capture how banks allocate the balance sheet between
  the sovereign and the private sector:

    claimsgov_assets    = Total claims on government / Total assets   (%)
    claimpriv_assets    = Claims on private sector  / Total assets    (%)
    netclaimsgov_assets = Net claims on government  / Total assets    (%)

  claimsgov_assets is the headline "doom-loop" intensity measure (the
  portfolio share of bank balance sheets tied to the sovereign).
  claimpriv_assets is its private-sector counterpart; together they show
  reallocation between sovereign and private claims during spread crises.

  Coverage notes:
    - Years 2001-2024 only  -> pre-2001 onsets have no nexus data.
    - 6 panel countries are absent from the IMF file and stay missing:
      China (mainland), Ecuador, El Salvador, India, Lebanon, Vietnam.

  Merge key: country (IMF names harmonised to panel names) x year.
  Run AFTER 01b_merge_new_controls.do.
===========================================================================*/

* ── Harmonise IMF country names to the panel's short names ────────────────
* Applied to the `country' string after each import. Uses strpos() on a
* distinctive substring to be robust to accents/punctuation/suffixes.
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

* ══════════════════════════════════════════════════════════════════════════
* 1. TOTAL CLAIMS ON GOVERNMENT / ASSETS  (already a ratio in the workbook)
* ══════════════════════════════════════════════════════════════════════════
import excel "`nexus'", sheet("Tot_claim_gov_pct") cellrange(A2) firstrow clear
_nexusnames
destring Totclaimsgov_assetpct*, replace force
reshape long Totclaimsgov_assetpct, i(country) j(year)
gen claimsgov_assets = 100 * Totclaimsgov_assetpct     // fraction -> percent
drop Totclaimsgov_assetpct
tempfile tg
save `tg'

* ══════════════════════════════════════════════════════════════════════════
* 2. NET CLAIMS ON GOVERNMENT / ASSETS
* ══════════════════════════════════════════════════════════════════════════
import excel "`nexus'", sheet("Net_claim_gov_pct") cellrange(A2) firstrow clear
_nexusnames
destring netclaimsgov_assetpct*, replace force
reshape long netclaimsgov_assetpct, i(country) j(year)
gen netclaimsgov_assets = 100 * netclaimsgov_assetpct
drop netclaimsgov_assetpct
tempfile tn
save `tn'

* ══════════════════════════════════════════════════════════════════════════
* 3. CLAIMS ON PRIVATE SECTOR (level, LCU)
* ══════════════════════════════════════════════════════════════════════════
import excel "`nexus'", sheet("Claim_priv") cellrange(A2) firstrow clear
_nexusnames
destring claim_priv*, replace force
reshape long claim_priv, i(country) j(year)
tempfile tp
save `tp'

* ══════════════════════════════════════════════════════════════════════════
* 4. TOTAL ASSETS (level, LCU) — from the stacked Source sheet
*    Source layout: A=country, B=var, C.. = years (no clean header names),
*    so import without firstrow (cols named A, B, C, ...) and rename.
* ══════════════════════════════════════════════════════════════════════════
import excel "`nexus'", sheet("Source") cellrange(A2) clear
drop if A == "country"          // the header row read in as data
keep if B == "Assets"
rename A country
drop B
* C, D, ... are the year columns 2001, 2002, ... in order
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

* ══════════════════════════════════════════════════════════════════════════
* 5. ASSEMBLE NEXUS COUNTRY-YEAR FILE
* ══════════════════════════════════════════════════════════════════════════
use `tp', clear
merge 1:1 country year using `ta', keep(match) nogen
gen claimpriv_assets = 100 * claim_priv / assets if assets > 0 & !missing(assets)
drop claim_priv assets

merge 1:1 country year using `tg', nogen
merge 1:1 country year using `tn', nogen

label var claimsgov_assets    "Bank claims on govt / assets (%)"
label var claimpriv_assets    "Bank claims on private / assets (%)"
label var netclaimsgov_assets "Bank net claims on govt / assets (%)"

keep country year claimsgov_assets claimpriv_assets netclaimsgov_assets
tempfile nexus_cy
save `nexus_cy'

* ══════════════════════════════════════════════════════════════════════════
* 6. MERGE INTO PANEL
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp.dta", clear
merge m:1 country year using `nexus_cy', keep(master match) nogen

sort cid year
xtset cid year
save "$clean/panel_lp.dta", replace

* ── Coverage report ───────────────────────────────────────────────────────
di as result _n "=== NEXUS COVERAGE AT ONSET (onset_all==1 & sample==1) ==="
foreach v in claimsgov_assets claimpriv_assets netclaimsgov_assets {
    quietly count if onset_all == 1 & sample == 1 & !missing(`v')
    di as result "  `v': " r(N) " / 61 onsets with non-missing data"
}
di as result _n "=== NEXUS COVERAGE: ALL SAMPLE YEARS (sample==1) ==="
foreach v in claimsgov_assets claimpriv_assets netclaimsgov_assets {
    quietly count if sample == 1 & !missing(`v')
    di as result "  `v': " r(N) " obs"
}

di as result _n "01c_merge_nexus.do complete — panel_lp.dta updated."
