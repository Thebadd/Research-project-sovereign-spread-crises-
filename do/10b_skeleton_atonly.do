/*===========================================================================
  10B_SKELETON_ATONLY.DO
  Extends the country-year skeleton (built in 10_skeleton.do) to cover every
  country in the full Asonuma-Trebesch database that is NOT already one of
  this project's 52 spread-crisis-panel countries. Sits between
  10_skeleton.do and 11_weo.do in the active chain — inserted here, not
  later, because stages 11-17 merge macro data onto whatever country-year
  rows already exist in panel_build.dta; a new country needs its skeleton
  row to exist BEFORE those generic merges run, or it never receives any
  macro data at all.

  WHY THIS FILE EXISTS
  ---------------------
  Per Pescatori & Sy (2007) — "we define debt crises as events occurring
  when either a country defaults or its bond spreads are above a critical
  threshold" — AT's own restructuring record is a primary, independent
  authority on default status; it does not require spread confirmation to
  count (METHODOLOGY.md §5.3 states this precisely, including why the
  reverse claim, tried and rejected earlier in this project's own design
  discussion, is wrong). This project's spread-crisis database was only
  ever built for 52 countries, so a country with a genuine AT-recorded
  default but no spread-crisis coverage was, until this file, entirely
  invisible to the panel. Checked directly (data/raw/WEOApr2026all.xlsx:
  197 countries; the WDI extracts: 265; IMF MFS Nexus_Sovereign_Bank.xlsx:
  148; Laeven-Valencia SYSTEMIC_BANKING_CRISES_DATABASE_2026.xlsx: 135) —
  every downstream macro source is GLOBAL, not pre-filtered to the 52, so
  no new raw data needs sourcing; only the skeleton needs extending, and
  11_weo.do through 17_predictors.do's existing generic iso3/country x year
  merges do the rest automatically, unchanged.

  WHAT THESE NEW COUNTRIES CAN AND CANNOT CONTRIBUTE
  -----------------------------------------------------------------------
  They can ONLY ever populate the default-linked side of the broadened
  debt-crisis taxonomy (17b_merge_at_full.do / 26_lp_debtcrisis_flow.do).
  They can NEVER contribute a non-default spread-crisis observation
  (onset_nd stays structurally 0 for every row below), because that
  category is a positive finding under the spread criterion specifically,
  and these countries were never run through it (METHODOLOGY.md §5.3,
  "Consequence for country coverage"). Every spread-side column below
  (spr_max, spr_mean, crit1, crit2, crisis_any, onset_all, onset_nd,
  onset_def, nondefault, continuation, ep_id, ep_status, classification) is
  therefore fixed at its "never tested" value (0 or missing, as
  appropriate) — not estimated, not guessed.

  WINDOW SIZING: WIDE ENOUGH TO IDENTIFY, NOT JUST TO COVER THE CRISIS
  -----------------------------------------------------------------------
  A country contributing ONLY its own crisis year(s), with no surrounding
  untreated years, is fully absorbed by its own fixed effect under this
  project's country+year FE design and identifies nothing (confirmed
  directly for the narrower "only 3 of 39 countries have multiple,
  differently-typed episodes" case explored earlier). Giving every new
  country a real surrounding window — not just its AT episode years —
  fixes this: the buffer years are genuine untreated (dc_default_year==0)
  rows, so even a single-episode country now has real within-country
  contrast to identify against. Per country: `yr_min` = 3 years before its
  EARLIEST AT-recorded episode onset (matching 10_skeleton's own carry-in
  lag-buffer convention); `yr_max` = this panel's own final year (2026),
  uniformly, so every new country gets the same generous forward-looking
  control pool the original 52 already have, regardless of how long ago
  its own AT history ends.

  REGION IS LEFT MISSING, DELIBERATELY — NOT AN OVERSIGHT
  -----------------------------------------------------------------------
  17_predictors.do's `reg_crisis_share` contagion predictor is computed as
  a share WITHIN each region group. Assigning a new country to one of the
  EXISTING 52 countries' region codes would change that share's denominator
  for every existing country in that region, silently altering every AIPW
  file's (21, etc.) already-published first-stage numbers — unacceptable
  under this project's own discipline of not disturbing existing files'
  results. Leaving region missing puts every new-country row in its own
  separate ("missing") bysort group, which changes nothing for the
  existing 52. `reg_crisis_share`/`contagion_dist` are not used anywhere in
  `26_lp_debtcrisis_flow.do` in any case, so nothing here is lost for this
  file's own purpose.

  Outputs: appends new-country rows onto $clean/panel_build.dta with the
  identical column schema 10_skeleton.do produces, so 11-18 need no changes
  at all to pick them up.
===========================================================================*/

use "$clean/panel_build.dta", clear
quietly levelsof iso3, local(existing_iso3)
local n_existing : word count `existing_iso3'
di as result _n "10b_skeleton_atonly: `n_existing' existing panel countries found."

* ══════════════════════════════════════════════════════════════════════════
* 1. IMPORT + COLLAPSE THE FULL AT DATABASE — same positional-column import
*    and 24-month collapse idiom as 17b_merge_at_full.do (this file only
*    needs iso3, a clean country name, and each country's episode-year span
*    for window sizing; the actual at_type/preemptive-post classification
*    is still built later, unchanged, by 17b once these rows exist to merge
*    onto).
* ══════════════════════════════════════════════════════════════════════════
preserve
    import excel "$raw/Asonuma_Trebesch_full_database.xlsx", ///
        sheet("DATASET Defaults & Restruct.") cellrange(B7) clear
    rename B case_nr
    rename C new_case
    rename D case_cruces
    rename E country_case
    rename F iso3_raw
    rename G start_date
    rename H end_date
    rename I alt_end_date
    rename J strictly_preempt
    rename K weakly_preempt
    rename L post_default
    rename M default_date
    rename N announcement_date
    rename O no_exact_start
    capture drop P

    keep if !missing(country_case)
    rename iso3_raw iso3
    replace iso3 = trim(iso3)
    * Same remap 17b_merge_at_full.do applies -- must match here too, or
    * Romania (already an existing panel country) would incorrectly appear
    * as a spurious "new" country under its AT code "ROM".
    replace iso3 = "ROU" if iso3 == "ROM"
    keep if !missing(iso3)

    gen int start_year = year(start_date) if !missing(start_date)
    gen int end_year   = year(end_date) if !missing(end_date)
    replace end_year = start_year if missing(end_year)
    drop if missing(start_year)

    * Restrict to countries NOT already in the panel -- everything else in
    * this file only concerns those.
    gen byte _is_new = 1
    foreach c of local existing_iso3 {
        quietly replace _is_new = 0 if iso3 == "`c'"
    }
    keep if _is_new == 1
    drop _is_new

    quietly count
    if r(N) == 0 {
        di as result "  No AT-only countries found beyond the existing panel -- nothing to add."
        di as result "10b_skeleton_atonly.do complete (no-op)."
        exit 0
    }

    * Clean display name: strip AT's parenthetical creditor-track suffixes
    * ("Dominican Rep. (bank debt)" -> "Dominican Rep."; irrelevant for the
    * genuinely new countries kept here, but harmless to apply uniformly).
    gen str80 country_clean = strtrim(regexr(country_case, "[ ]*\([^)]*\)[ ]*$", ""))

    * ── Collapse to one row per NEW country: earliest onset, latest end ────
    bysort iso3: egen int g_onset = min(start_year)
    bysort iso3: egen int g_end   = max(end_year)
    duplicates drop iso3, force
    keep iso3 country_clean g_onset g_end
    rename country_clean country

    local n_new = _N
    quietly levelsof iso3, local(new_iso3)
    di as result "  AT-only countries beyond the existing `n_existing'-country panel: `n_new'"
    list iso3 country g_onset g_end, noobs

    tempfile atonly_countries
    save `atonly_countries'
restore

* ══════════════════════════════════════════════════════════════════════════
* 2. BUILD THE SKELETON ROWS — one row per new country per year in
*    [earliest AT onset - 3, panel's own final year], matching
*    10_skeleton.do's column schema exactly so 11-18 need no changes.
* ══════════════════════════════════════════════════════════════════════════
quietly summarize year
local panel_ymax = r(max)

preserve
    use `atonly_countries', clear
    gen int yr_min = g_onset - 3
    gen int yr_max = `panel_ymax'
    gen int nyears = yr_max - yr_min + 1
    expand nyears
    bysort iso3: gen int year = yr_min + _n - 1

    * ── Spread-side columns: fixed at "never tested" -- see header ─────────
    gen double spr_max = .
    gen double spr_mean = .
    gen byte crit1 = 0
    gen byte crit2 = 0
    gen byte crisis_any = 0
    gen byte onset_all = 0
    gen byte onset_nd = 0
    gen byte onset_def = 0
    gen byte nondefault = .
    gen byte continuation = 0
    gen str1 ep_id = ""
    gen str1 ep_status = ""
    gen str1 classification = ""
    gen str1 region = ""          // deliberately missing -- see header

    * carry-in: the leading 3 years of each new country's own window are lag
    * scaffolding, identical role and exclusion rule to 10_skeleton.do's own
    * carry-in rows (never in sample_base / sample / sample_flow).
    gen byte carryin = (year < g_onset)
    gen byte sample_base = (continuation == 0) & carryin == 0

    * MUST-EXCLUDE flag: every onset/flow file 02-25 reads $clean/panel_lp.dta's
    * global `sample'/`sample_flow' flags, built once in 18_transforms.do. Without
    * this flag, adding these rows would silently dilute EVERY existing file's
    * control pool with untested-for-spread country-years -- exactly the kind of
    * change METHODOLOGY.md §5.3 promises does not happen. 18_transforms.do
    * excludes atonly_country==1 from `sample'/`sample_flow'; only
    * 26_lp_debtcrisis_flow.do builds its own broadened sample flag that
    * includes these rows.
    gen byte atonly_country = 1
    label var atonly_country "1 = AT-only country (10b), excluded from sample/sample_flow; used only by 26"

    label var iso3 "ISO3 country code (universal merge key)"
    label var carryin "1 = pre-entry scaffolding row (lag support only; never in sample)"
    label var sample_base "Onset + tranquil years (excl. continuation & carry-in); GDP-availability added in 18"

    keep iso3 country region year spr_max spr_mean crit1 crit2 crisis_any ///
         onset_all onset_nd onset_def nondefault continuation ep_id ep_status ///
         classification carryin sample_base atonly_country

    local n_rows = _N
    tempfile new_skeleton
    save `new_skeleton'
restore

* ══════════════════════════════════════════════════════════════════════════
* 3. APPEND onto the existing skeleton and re-save
* ══════════════════════════════════════════════════════════════════════════
append using `new_skeleton'
* Existing 52-country rows never had atonly_country -- append leaves them
* missing, not 0. Fix explicitly rather than relying on missing-as-false.
capture confirm variable atonly_country
if _rc gen byte atonly_country = 0
replace atonly_country = 0 if missing(atonly_country)
label var atonly_country "1 = AT-only country (10b), excluded from sample/sample_flow; used only by 26"

sort iso3 year
save "$clean/panel_build.dta", replace

di as result _n "10b_skeleton_atonly.do complete."
di as result "  Added `n_new' AT-only countries, `n_rows' country-year rows."
quietly count if carryin == 1
di as result "  New rows are carry-in (lag scaffolding): (subset of the total above)"
quietly levelsof iso3, local(all_iso3)
di as result "  Panel now covers " `: word count `all_iso3'' " countries total."
