/*===========================================================================
  16_BANKING.DO   —  FROM-SCRATCH REBUILD, STAGE 6  (banking-sector crises)
  Systemic banking-crisis dummy from the ORIGINAL Laeven & Valencia database,
  2026 update: data/raw/SYSTEMIC_BANKING_CRISES_DATABASE_2026.xlsx.

    banking_crisis  -> 1 in every year of a systemic banking crisis, 0 otherwise

  WHY THIS REPLACED THE OLD GFDD EXTRACT (this is not a cosmetic swap):
  The previous source (Banking_crisis_dummy_data_set.xlsx, World Bank GFDD series
  GFDD.OI.19) is a country-year panel whose year columns are 100% MISSING from
  2018 onward. `banking_crisis` therefore existed only 1990-2017, so
  `l_banking_crisis = L.banking_crisis` existed only 1991-2018 — and because l_banking_crisis
  sits in $ctrl_core, that single control truncated the WHOLE pipeline at 2018
  (every xtscc table showed "2019-2026 (omitted)"). Listwise deletion then cut
  the headline LP from 61 onsets to ~20, excluding Argentina 2018, Turkey 2018,
  Lebanon 2019, Zambia 2020, Sri Lanka 2022, Ghana 2022 and COVID.

  The L-V database is an EPISODE LIST (one row per crisis, with Start and End
  years), not a panel. That distinction is the fix: a country-year absent from
  the list is a genuine ZERO (no systemic banking crisis), not an unknown. We
  expand each [Start, End] window to country-years, then zero-fill the rest of
  the panel — so no observation is silently deleted for want of a banking flag.

  Source layout (verified): sheet "Crisis Resolution and Outcomes", header on
  row 1, data from row 2, columns A=Country, B=blank, C=Start, D=End.
  164 episodes, 120 countries, 1976-2023. Footnote markers ("5/") appear inside
  BOTH the country and the year cells and are stripped below.
  NOTE: the "Resolution Details" sheet is NOT used — its date cells are corrupted
  (Argentina's 2001 crisis reads as 2026-11-01); this sheet has clean integer years.

  Coverage vs. our panel: 43 of 52 countries match L-V by exact name; Türkiye is
  the only naming mismatch (handled below). The remaining 8 — Belize, Guatemala,
  Honduras, Namibia, Pakistan, Serbia, South Africa, Trinidad and Tobago — carry
  no L-V systemic banking crisis at all, so they are legitimately 0 throughout.

  ONGOING CRISES (David's decision): a blank End year means the crisis had not
  ended when the database was published, so the country is STILL in crisis and
  stays coded 1 through the last panel year (Vietnam 2022-, Sri Lanka 2023-).

  Output: merges banking_crisis onto $clean/panel_build.dta (country x year).
===========================================================================*/

* ── 0. Panel span (drives the ongoing-crisis fill and the year window) ──────
use "$clean/panel_build.dta", clear
quietly summarize year
local ymin = r(min)
local ymax = r(max)
di as result _n "16_banking: panel spans `ymin'-`ymax'."

* ── 1. Import the crisis episode list ───────────────────────────────────────
* cellrange(A2) skips the header row, so columns arrive as Excel letters A, B, C, D
* (the project's rename-by-position idiom, as in 10_skeleton / 13_ifs_nexus). The
* header is multi-line and footnoted, so `firstrow' would produce unusable names.
import excel "$raw/SYSTEMIC_BANKING_CRISES_DATABASE_2026.xlsx", ///
    sheet("Crisis Resolution and Outcomes") cellrange(A2) allstring clear

foreach v in A C D {
    capture confirm variable `v'
    if _rc {
        di as error "  ** 16_banking: expected column `v' not found in the workbook."
        di as error "     Sheet must be 'Crisis Resolution and Outcomes' with A=Country, C=Start, D=End."
        exit 111
    }
}
keep A C D
rename (A C D) (country start_raw end_raw)

* ── 2. Clean footnote markers ("Argentina 5/", "1994 5/", stray double spaces) ─
replace country   = strtrim(itrim(regexr(country,   "[ ]*[0-9]+/[ ]*$", "")))
replace start_raw = strtrim(regexr(start_raw, "[ ]*[0-9]+/[ ]*$", ""))
replace end_raw   = strtrim(regexr(end_raw,   "[ ]*[0-9]+/[ ]*$", ""))

* Years: pull the 4-digit year out of each cell (robust to leftover text)
gen int start_yr = .
gen int end_yr   = .
replace start_yr = real(regexs(0)) if regexm(start_raw, "(19|20)[0-9][0-9]")
replace end_yr   = real(regexs(0)) if regexm(end_raw,   "(19|20)[0-9][0-9]")

* Drop the trailing footnote/blank rows (they have no country or no start year)
drop if missing(country) | country == "" | missing(start_yr)
di as result "  episodes parsed: " _N

* ── 3. Harmonise the one mismatched country name ────────────────────────────
* strpos on an ASCII substring avoids any UTF-8 issue with "Türkiye" — the same
* trick _nexusnames uses in 13_ifs_nexus.do.
replace country = "Turkey" if strpos(country, "rkiye")

* ── 4. Ongoing crises: blank End = still in crisis -> run to the panel end ───
quietly count if missing(end_yr)
di as result "  ongoing crises (no End year) extended to `ymax': " r(N)
replace end_yr = `ymax' if missing(end_yr)
replace end_yr = start_yr if end_yr < start_yr      // data sanity

* ── 5. Expand each [Start, End] window into country-years ───────────────────
gen long epid = _n                                   // unique episode id

* BANKING-CRISIS DURATION. The reference paper controls for
* banking_duration_lv2018 (Laeven-Valencia 2018), a duration rather than a flag.
* Their replication files import it pre-built and only lag it, so the exact
* definition is not recoverable from their code and we have to choose one.
*
* banking_duration = HOW LONG THE CRISIS HAS ALREADY LASTED at that year:
* a counter running 1, 2, 3, ... from the start year, 0 outside a crisis.
* Ghana's 2015-2019 crisis is 1 in 2015, 2 in 2016, ... 5 in 2019.
*
* This is the right form for a control here, for two reasons the total-length
* alternative fails:
*   (a) it is genuinely PREDETERMINED. Years-so-far uses only information
*       available at t. Total episode length embeds the End year, which is not
*       knowable at t, so lagging it to t-1 would not make it t-1 information the
*       way every other core control is.
*   (b) it needs no censoring assumption. The two ONGOING crises (blank End,
*       extended to `ymax' in step 4) have a well-defined years-so-far at every
*       year; their TOTAL length is unknown and would be understated.
* It also nests the dummy: any positive value means "in a banking crisis", so
* l_banking_crisis is the special case where the counter is collapsed to 0/1.
*
* banking_duration_total (end - start + 1, constant within an episode) is built
* alongside as the robustness alternative, carrying both caveats above.
gen int banking_duration_total = end_yr - start_yr + 1

expand end_yr - start_yr + 1
bysort epid: gen int year = start_yr + _n - 1
gen byte banking_crisis = 1
bysort epid (year): gen int banking_duration = _n

* One row per country-year. (max) on overlapping/repeat episodes: a country-year
* covered by two episodes takes the longer count, consistent with the dummy
* taking 1.
collapse (max) banking_crisis banking_duration banking_duration_total, ///
    by(country year)
keep if inrange(year, `ymin', `ymax')

tempfile bank
save `bank'

* ── 6. Merge onto the panel and ZERO-FILL ───────────────────────────────────
use "$clean/panel_build.dta", clear
capture drop banking_crisis banking_duration banking_duration_total
merge m:1 country year using `bank', keep(master match) nogen

* THE FIX: absence from the L-V episode list means "no systemic banking crisis",
* not "unknown". Without this line every non-crisis country-year stays missing and
* is dropped from every regression by listwise deletion on $ctrl_core.
* The durations zero-fill on the same logic: no crisis => a crisis lasting 0 years.
replace banking_crisis         = 0 if missing(banking_crisis)
replace banking_duration       = 0 if missing(banking_duration)
replace banking_duration_total = 0 if missing(banking_duration_total)

label var banking_crisis   "Systemic banking-crisis dummy (Laeven-Valencia 2026; 1 during each crisis)"
label var banking_duration "Years the banking crisis has lasted so far (1,2,3...); 0 if none - COMMON CORE"
label var banking_duration_total "Total banking-crisis length, years (end-start+1); 0 if none - robustness alt."
save "$clean/panel_build.dta", replace

* ── 7. Coverage report ──────────────────────────────────────────────────────
di as result _n "16_banking.do complete. Coverage:"

* The headline check: banking_crisis must now be non-missing EVERYWHERE. Any
* remaining missing would again delete observations via $ctrl_core listwise deletion.
quietly count if missing(banking_crisis)
di as result "  rows still missing banking_crisis: `r(N)'   (must be 0; the old GFDD"
di as result "     extract left every 2018+ row missing, truncating the LP at 2018)"
quietly count if !missing(banking_crisis) & sample_base==1
di as result "  banking_crisis non-missing (sample rows): `r(N)'"

quietly summarize year if banking_crisis==1
di as result "  crisis years span: " r(min) "-" r(max)
quietly count if banking_crisis==1
di as result "  crisis country-years in panel: `r(N)'"
quietly count if banking_crisis==1 & onset_all==1
di as result "  onsets coinciding with a banking crisis: `r(N)' of 61"

* Duration diagnostics: the core now carries banking_duration, so its spread and
* its coverage at onsets both matter.
quietly summarize banking_duration if banking_crisis==1
di as result "  years-so-far within a crisis: mean " %4.1f r(mean) "  min " r(min) "  max " r(max)
di as result "     (min must be 1 - every crisis has a first year)"
quietly summarize banking_duration_total if banking_crisis==1
di as result "  total crisis length (robustness alt.): mean " %4.1f r(mean) "  min " r(min) "  max " r(max)
quietly count if !missing(banking_duration)
di as result "  rows with a non-missing banking_duration: `r(N)' (must equal the panel)"

* Name-mismatch guard. 44 of our 52 countries appear in L-V (43 exact name matches +
* Türkiye), but only 23 have a crisis inside their OWN panel coverage. Two reasons,
* both legitimate: (a) 11 countries (Chile, Cote d'Ivoire, Egypt, El Salvador, India,
* Jordan, Morocco, Panama, Peru, Senegal, Tunisia) were last hit before 1994; (b) the
* panel is UNBALANCED — it starts when a country enters the EMBIG index — so another
* ~10 (Armenia, Bolivia, Costa Rica, Indonesia, Jamaica, Kenya, Paraguay, Poland,
* Romania, Zambia) have only early/mid-1990s L-V crises that fall before their first
* panel year.
*   EXPECT ROUGHLY 23-26. It was exactly 23 before 10_skeleton began adding carry-in
*   rows; those extend each country's coverage 3 years earlier, so a couple of 1990s
*   crises (Poland, Indonesia) now fall inside it. Anything in that band is fine. A
*   number materially BELOW it means a country name failed to merge and its zeros are
*   FALSE zeros, so the crisis-free list is printed below for eyeballing.
preserve
    quietly keep if banking_crisis==1
    quietly levelsof country, local(bc)
    local nbc : word count `bc'
restore
di as result "  panel countries with >=1 banking crisis in own coverage: `nbc' (expect ~23-26 of 52)"

* Expected crisis-free list (29): the 8 never in L-V at all (Belize, Guatemala,
* Honduras, Namibia, Pakistan, Serbia, South Africa, Trinidad and Tobago), the 11 last
* hit before 1994, and the 10 whose crises predate their own EMBIG entry.
preserve
    collapse (max) banking_crisis, by(country)
    quietly levelsof country if banking_crisis==0, local(nocrisis)
    di as result "  panel countries with no crisis in own coverage (expect ~26-29):"
    foreach c of local nocrisis {
        di as result "      `c'"
    }
restore

* ── 8. Spot-checks that must hold if the episode expansion is right ─────────
*   Ghana        2017 = 1  (L-V crisis 2015-2019)
*   South Africa 2010 = 0  (no L-V crisis ever)
*   Sri Lanka    2024 = 1  (crisis from 2023, still ongoing -> runs to panel end)
*   Turkey       2001 = 1  (verifies the Türkiye -> Turkey name fix worked)
capture program drop _bankcheck
program define _bankcheck
    args cty yr want
    quietly count if country == "`cty'" & year == `yr' & banking_crisis == `want'
    if r(N) > 0  di as result "  check OK:     `cty' `yr' = `want'"
    else         di as error  "  ** check FAILED: `cty' `yr' should be `want'"
end
_bankcheck "Ghana"        2017 1
_bankcheck "South Africa" 2010 0
_bankcheck "Sri Lanka"    2024 1
_bankcheck "Turkey"       2001 1
