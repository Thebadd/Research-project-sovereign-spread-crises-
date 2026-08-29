/*===========================================================================
  17B_MERGE_AT_FULL.DO
  Full Asonuma & Trebesch (2016, JEEA) default/restructuring database —
  merged onto panel_build.dta as AT-SIDE columns only, independent of the
  spread-crisis dating built later in 18_transforms.do. Sits between
  17_predictors.do (which produces panel_build.dta) and 18_transforms.do
  (which consumes it and produces panel_lp.dta) in the active 10-18 build
  chain — NOT part of the retired 01_build_panel..01e_predictors chain.

  WHY THIS FILE EXISTS
  --------------------
  The project's existing default-linked classification only checks whether
  an AT restructuring coincides with one of the 61 EMBIG-based spread-crisis
  episodes (Episode_Summary sheet of EM_Spread_Crisis_DB_FINAL.xlsx). That
  means a genuine sovereign default that never produced a detected spread
  crisis is invisible to this project. This file imports the FULL AT
  database (197 case rows, not scoped to this project's existing episodes)
  so that later files can identify default-linked debt crises EITHER via a
  spread crisis that coincides with a restructuring, OR directly via AT when
  it does not.

  Source file: data/raw/Asonuma_Trebesch_full_database.xlsx, sheet
  "DATASET Defaults & Restruct." (the original Asonuma-Trebesch 2016 JEEA
  database, header row 6, columns B:O; column A is blank in the source).

  SCOPE: THIS FILE ONLY BUILDS THE AT-SIDE COLUMNS. It does NOT build the
  combined "debt crisis" classification (that needs onset_all/in_crisis,
  which do not exist until 18_transforms.do runs) — 26_lp_debtcrisis_flow.do
  does that combination, self-contained, the same way 25_aipw_nexus_split_
  flow.do builds its own a_nexus rather than modifying 18_transforms.do.

  TWO DESIGN DECISIONS, STATED EXPLICITLY, NOT LEFT IMPLICIT
  -------------------------------------------------------------------------
  1. SUB-CASE COLLAPSING. 197 raw case rows cover many crises split across
     several creditor tracks (e.g. Ukraine has 6 separate rows: Eurobonds,
     Chase loan, Commercial loans, Global Exchange, ING debt/Merrill Lynch,
     OVDPs non-resid.). Rows for the same country whose start years are
     within 24 months of each other are collapsed into ONE default episode:
     onset = the EARLIEST start year among the group; end = the LATEST end
     year among the group; type = the most informative classification
     present, checked STRICTLY PREEMPTIVE first, then weakly preemptive,
     then post-default (if a country restructures one creditor track
     post-default in the same window, the episode is a post-default crisis
     at the sovereign level regardless of any other track in it).
     THIS ORDER MATTERS AND WAS INITIALLY WRONG: an earlier version of this
     file checked weakly preemptive before strictly preemptive, on the
     reasonable-looking but false assumption that the three flags are
     mutually exclusive severity tiers. A diagnostic (listing every raw row
     with strictly_preempt==1) showed every one of those 27 rows ALSO has
     weakly_preempt==1 set on the identical row — strictly preemptive is
     coded as a stricter SUB-condition of weakly preemptive in this AT
     vintage, not a separate category. Checking weakly first meant
     "Strictly preemptive" could never be assigned to any episode, for any
     panel built from this file — not rare, structurally unreachable. Since
     strictly preemptive is the economically LEAST disruptive outcome (no
     missed payment at all, not even brief technical arrears), it is
     checked first now, so a genuinely strictly-preemptive episode is no
     longer masked by its own co-occurring weakly-preemptive flag.
  2. VINTAGE CUTOFF. The source file's latest case starts Dec 2019 (plus a
     handful "ongoing as of Sept 2020"). Nothing after ~2020 is in it —
     confirmed absent: Zambia 2020-23, Sri Lanka 2022, Ghana 2022, Ethiopia
     2021, Suriname 2020, Belarus 2022. A USER-MAINTAINED supplement block
     below lets these be added by hand, with their own source/date verified
     by the user rather than guessed by this file.

  RANGE-MERGE MECHANICS: Stata has no native range-join, so this uses the
  standard `joinby' (many-to-many on iso3) + `keep if inrange(year,...)'
  idiom: join every panel row to every AT episode for that country, then
  keep only the rows where the panel year actually falls inside that
  episode's [onset_year, end_year] window, then reduce back to one row per
  country-year (flagging, not silently resolving, the rare case where two
  AT episodes for one country somehow overlap in year).

  Outputs added to panel_build.dta (all AT-side, independent of spread coding)
  -------------------------------------------------------------------------
    at_default_year     1 if this country-year falls inside an AT-recorded
                        default/restructuring window (collapsed episode)
    at_episode_onset    the containing AT episode's own onset year
    at_type             "Strictly preemptive" / "Weakly preemptive" /
                        "Post-default" for the containing episode
    at_overlap_flag     1 if two AT episodes for this country overlap in
                        year (should be rare; printed, not silently resolved)
===========================================================================*/

use "$clean/panel_build.dta", clear
capture xtset cid year
sort cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. IMPORT RAW AT DATABASE — positional column names, not firstrow
*
* The source headers are long free-text strings ("Start of default or
* restructuring process: default or announcement") that Stata would auto-
* truncate/sanitize unpredictably under `firstrow'. Importing by POSITION
* from the data row and naming columns explicitly avoids guessing at
* Stata's sanitized names.
* ══════════════════════════════════════════════════════════════════════════
preserve
    import excel "$raw/Asonuma_Trebesch_full_database.xlsx", ///
        sheet("DATASET Defaults & Restruct.") cellrange(B7) clear

    * import excel WITHOUT firstrow names variables after their ACTUAL
    * spreadsheet column letters, not a fresh A-onward sequence starting at
    * the cellrange's first column. Data starts at column B, so the first
    * variable is named B, not A. The sheet's used range also extends one
    * column further than the data (a blank trailing column P), dropped below.
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
    * WDI code is already ISO3 — the whole point of using this column
    * directly rather than fuzzy-matching `Country case' text.
    rename iso3_raw iso3
    replace iso3 = trim(iso3)
    keep if !missing(iso3)

    gen int start_year = year(start_date) if !missing(start_date)
    gen int end_year = year(end_date) if !missing(end_date)
    * `alt_end_date' ("Alternative end date / follow-up restructurings") is
    * FREE-TEXT in the source (e.g. "Two-follow up deals: Sept. 2010 ...
    * April 2016"), not a clean date -- it imports as a string, and
    * attempting year() on it throws a type mismatch (confirmed the hard
    * way). Not used for end_year; the collapse step's own 24-month grouping
    * window already absorbs a documented follow-up as a separate nearby
    * case rather than needing this field parsed.
    replace end_year = start_year if missing(end_year)
    drop if missing(start_year)

    di as result "  AT full database: `=_N' raw case rows imported, " ///
        "`=_N' after dropping missing-start-year rows"

    * ── Collapse sub-case rows into episodes ────────────────────────────
    * Group rows for the same country whose start years are within 24
    * months of each other. Running-counter idiom, same pattern as
    * 18_transforms.do's ep_seq.
    sort iso3 start_year
    by iso3: gen gap = start_year - start_year[_n-1]
    by iso3: gen byte new_group = (_n==1) | (gap > 2)
    gen long grp = sum(new_group)

    * Severity ranking for the collapsed type: post-default > weakly
    * preemptive > strictly preemptive — stated in the header, not implicit.
    bysort grp: egen byte g_post    = max(post_default)
    bysort grp: egen byte g_weakly  = max(weakly_preempt)
    bysort grp: egen byte g_strict  = max(strictly_preempt)
    bysort grp: egen int  g_onset   = min(start_year)
    bysort grp: egen int  g_end     = max(end_year)
    bysort grp: egen byte g_ncases  = count(start_year)

    gen str20 at_type = ""
    replace at_type = "Strictly preemptive" if g_strict==1
    replace at_type = "Weakly preemptive"   if at_type=="" & g_weakly==1
    replace at_type = "Post-default"        if at_type=="" & g_post==1
    quietly count if at_type==""
    if r(N) > 0 di as error "  ** AT full database: " r(N) " collapsed episodes have no preemptive/post-default flag set — check source rows."

    * One row per collapsed episode.
    duplicates drop iso3 grp, force
    keep iso3 grp g_onset g_end at_type g_ncases
    rename g_onset at_episode_onset
    rename g_end   at_episode_end

    local n_episodes = _N
    quietly levelsof iso3, local(atc)
    * levelsof does not reliably set r(N) to the distinct-value count (it can
    * carry over r(N) from whatever command ran before it) -- count the
    * words in the returned local instead. Confirmed the hard way: this
    * printed "127 countries" (the EPISODE count) before the fix.
    local n_countries : word count `atc'
    di as result "  AT full database: collapsed to `n_episodes' default episodes " ///
        "across `n_countries' countries."

    tempfile at_episodes
    save `at_episodes'
restore

* ══════════════════════════════════════════════════════════════════════════
* 2. USER-MAINTAINED SUPPLEMENT — post-~2020 defaults not in the AT vintage
*
* Add one line per event you have identified and verified (source, date,
* preemptive/post-default) — do NOT let this be guessed from memory; AT's
* own rigor (verified announcement/default dates, creditor-track detail) is
* not something to approximate. Columns must match at_episodes exactly:
* iso3, at_episode_onset, at_episode_end, at_type, g_ncases (set g_ncases=1
* for a hand-added single-case entry unless you are merging several
* creditor tracks yourself, in which case set it to that count).
*
* LEFT EMPTY at build time. Filling this in requires no other code change —
* the merge below treats it identically to the AT-derived episode file.
* ══════════════════════════════════════════════════════════════════════════
tempfile at_supplement
preserve
    clear
    gen str3 iso3 = ""
    gen int at_episode_onset = .
    gen int at_episode_end = .
    gen str20 at_type = ""
    gen byte g_ncases = .
    * ── ADD ROWS HERE, e.g.: ──────────────────────────────────────────
    *   set obs `=_N'+1
    *   replace iso3 = "ZMB" in L
    *   replace at_episode_onset = 2020 in L
    *   replace at_episode_end   = 2023 in L
    *   replace at_type = "Post-default" in L
    *   replace g_ncases = 1 in L
    * ────────────────────────────────────────────────────────────────
    save `at_supplement'
restore

preserve
    use `at_episodes', clear
    quietly append using `at_supplement'
    quietly count
    di as result "  AT episodes + user supplement: `=_N' total episodes going into the range-merge."
    save `at_episodes', replace
restore

* ══════════════════════════════════════════════════════════════════════════
* 3. RANGE-MERGE ONTO THE PANEL — joinby (matches only) + merge back, not
*    joinby's own unmatched(master)
*
* CRITICAL FIX: `joinby iso3 using ..., unmatched(master)' explodes every row
* of a MATCHED country into one copy per that country's AT episodes, and
* `unmatched(master)' only rescues countries with ZERO AT episodes at all --
* it does NOTHING for a matched country's non-window years, since every one
* of that year's episode-copies fails the window test and gets dropped
* together. The result: any year of any country that ever appears in AT, but
* that isn't inside one of ITS OWN episode windows, was deleted from the
* panel entirely instead of being kept as an untreated (at_default_year=0)
* row -- confirmed the hard way (in_crisis collapsed from 234 to 81, onsets
* from 61 to 21, xtset reported "with gaps": actual row loss, not reordering).
*
* THE FIX: build the country-year x AT-episode MATCH SET on its own (joinby,
* THEN filter to window matches, THEN reduce to one row per country-year),
* and MERGE that small match-only file back onto the FULL original panel
* (`panel_pre_at', saved before any joinby) via a plain 1:1 merge on
* iso3/year. Every original panel row survives regardless of whether it
* matches; only the AT columns are added where a match exists.
* ══════════════════════════════════════════════════════════════════════════
capture drop at_default_year at_episode_onset at_type at_overlap_flag

tempfile panel_pre_at
save `panel_pre_at'

preserve
    * Country-year x AT-episode candidate matches ONLY -- countries with no
    * AT episode at all simply do not appear here, which is correct: they
    * will get at_default_year=0 for every row once merged back below.
    joinby iso3 using `at_episodes'
    gen byte _in_window = inrange(year, at_episode_onset, at_episode_end)
    quietly count if _in_window==0
    local ndrop = r(N)
    keep if _in_window==1
    drop _in_window
    di as result "  Range-merge: `ndrop' joinby candidate rows fell outside their AT episode's year window (dropped)."

    * Overlap check: two AT episodes covering the same country-year (rare --
    * flag, do not silently resolve). Flagged BEFORE reducing to one row per
    * country-year, so the flag itself survives the reduction below.
    duplicates tag iso3 year, gen(_dup)
    gen byte at_overlap_flag = (_dup > 0)
    quietly count if at_overlap_flag==1
    if r(N) > 0 {
        di as error "  ** AT full database: " r(N) " country-year rows have OVERLAPPING AT episodes — listed below, not silently resolved:"
        list iso3 year at_episode_onset at_episode_end at_type if at_overlap_flag==1, noobs sepby(iso3)
    }
    drop _dup

    * Reduce to exactly one row per country-year for the merge (an overlap's
    * second/third episode is not lost -- it is already flagged above; the
    * merge only needs one row to attach at_episode_onset/at_type to).
    bysort iso3 year: gen byte _first = (_n==1)
    keep if _first==1
    drop _first

    keep iso3 year at_episode_onset at_episode_end at_type at_overlap_flag
    tempfile at_matches
    save `at_matches'
restore

use `panel_pre_at', clear
local n_before = _N
merge 1:1 iso3 year using `at_matches', keep(master match) nogen

gen byte at_default_year = !missing(at_episode_onset)
replace at_overlap_flag = 0 if missing(at_overlap_flag)
label var at_default_year "1 if this country-year falls inside an AT-recorded default/restructuring episode window"
label var at_episode_onset "Onset year of the containing AT episode (collapsed)"
label var at_type "AT classification of the containing episode: Strictly/Weakly preemptive or Post-default"
label var at_overlap_flag "1 if two AT episodes overlap for this country-year (flagged, not resolved)"

* Every row of the original panel must survive this file untouched in count.
local n_after = _N
if `n_after' != `n_before' di as error "  ** AT full database: panel row count changed from `n_before' to `n_after'" ///
    " -- the range-merge lost or duplicated rows, investigate before saving."
else di as result "  Range-merge preserved all `n_after' original panel rows (correct)."

sort cid year
xtset cid year
save "$clean/panel_build.dta", replace

di as result _n "17b_merge_at_full.do complete."
quietly count if at_default_year==1
di as result "  Country-years inside an AT-recorded default episode: `r(N)'"
quietly count if at_default_year==1 & at_type=="Post-default"
di as result "    of which Post-default: `r(N)'"
quietly count if at_default_year==1 & at_type=="Weakly preemptive"
di as result "    of which Weakly preemptive: `r(N)'"
quietly count if at_default_year==1 & at_type=="Strictly preemptive"
di as result "    of which Strictly preemptive: `r(N)'"
