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
     year among the group; type = the MOST SEVERE classification present
     (post-default > weakly preemptive > strictly preemptive) — a country
     that restructures one creditor track post-default and another
     preemptively in the same crisis window is a post-default crisis at the
     sovereign level, not a preemptive one.
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
    gen int end_year_raw = year(end_date) if !missing(end_date)
    * Alternative/follow-up end date, when later, extends the episode window
    * (a case with a documented follow-up restructuring is still "in" that
    * crisis through the follow-up, not just the first exchange).
    gen int alt_end_year = year(alt_end_date) if !missing(alt_end_date)
    gen int end_year = max(end_year_raw, alt_end_year)
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
    by grp: egen byte g_post    = max(post_default)
    by grp: egen byte g_weakly  = max(weakly_preempt)
    by grp: egen byte g_strict  = max(strictly_preempt)
    by grp: egen int  g_onset   = min(start_year)
    by grp: egen int  g_end     = max(end_year)
    by grp: egen byte g_ncases  = count(start_year)

    gen str20 at_type = ""
    replace at_type = "Post-default"        if g_post==1
    replace at_type = "Weakly preemptive"   if at_type=="" & g_weakly==1
    replace at_type = "Strictly preemptive" if at_type=="" & g_strict==1
    quietly count if at_type==""
    if r(N) > 0 di as error "  ** AT full database: " r(N) " collapsed episodes have no preemptive/post-default flag set — check source rows."

    * One row per collapsed episode.
    duplicates drop iso3 grp, force
    keep iso3 grp g_onset g_end at_type g_ncases
    rename g_onset at_episode_onset
    rename g_end   at_episode_end

    local n_episodes = _N
    quietly levelsof iso3, local(atc)
    local n_countries = r(N)
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
* 3. RANGE-MERGE ONTO THE PANEL — joinby + inrange, no native range-join
* ══════════════════════════════════════════════════════════════════════════
capture drop at_default_year at_episode_onset at_type at_overlap_flag

tempfile panel_pre_at
save `panel_pre_at'

joinby iso3 using `at_episodes', unmatched(master)
* _merge-equivalent from joinby: rows with no AT episode for that country
* keep grp/at_type as missing already; only keep true window matches below,
* but a country with NO AT episode at all must still be RETAINED (unmatched
* master rows), so filter on the window condition OR "never had any AT row".
gen byte _in_window = inrange(year, at_episode_onset, at_episode_end)
quietly count if !missing(at_episode_onset) & _in_window==0
local ndrop = r(N)
drop if !missing(at_episode_onset) & _in_window==0
drop _in_window
di as result "  Range-merge: dropped `ndrop' joinby rows outside their AT episode's year window."

* Overlap check: two AT episodes for the same country-year (rare — flag,
* do not silently resolve).
duplicates tag iso3 year, gen(_dup)
gen byte at_overlap_flag = (_dup > 0) & !missing(at_episode_onset)
quietly count if at_overlap_flag==1
if r(N) > 0 {
    di as error "  ** AT full database: " r(N) " country-year rows have OVERLAPPING AT episodes — listed below, not silently resolved:"
    list iso3 year at_episode_onset at_episode_end at_type if at_overlap_flag==1, noobs sepby(iso3)
}
drop _dup

gen byte at_default_year = !missing(at_episode_onset)
label var at_default_year "1 if this country-year falls inside an AT-recorded default/restructuring episode window"
label var at_episode_onset "Onset year of the containing AT episode (collapsed)"
label var at_type "AT classification of the containing episode: Strictly/Weakly preemptive or Post-default"
label var at_overlap_flag "1 if two AT episodes overlap for this country-year (flagged, not resolved)"

* Rows can duplicate if joinby matched multiple episodes before the window
* filter (should not happen after the overlap flag/window filter above, but
* guard rather than assume).
duplicates tag cid year, gen(_dt)
quietly count if _dt > 0
if r(N) > 0 di as error "  ** AT full database: " r(N) " duplicated cid-year rows survived the range-merge — investigate before saving."
drop _dt

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
