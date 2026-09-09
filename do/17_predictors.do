/*===========================================================================
  17_PREDICTORS.DO   —  FROM-SCRATCH REBUILD, STAGE 7  (first-stage predictors Z)
  Derived exclusion-restriction predictors (no raw import). Mirrors the old 01e.
  Also creates the numeric panel id (cid) and xtset — the first stage that needs
  time-series operators.

    cid                 numeric country id (from iso3), + xtset cid year
    reg_crisis_share    leave-one-out share of OTHER same-region countries with
                        an onset in year t
    l_reg_crisis_share  Z2: lagged regional crisis share (contagion)
    contagion_dist      Z2b: distance-weighted SHARE (bounded [0,1], weights
                        normalized to sum to 1) of OTHER countries IN-CRISIS
                        years (onset|continuation, a STOCK of regional distress,
                        not just onset years) at year t, CEPII great-circle
                        distance (GEO_CEPII.xlsx) -- matches the reference
                        paper's own contagion predictor's bounded scale
                        (their Table B3: mean 0.05, range [0,0.88])
    l_contagion_dist    Z2b lagged (predetermined) -- ADOPTED in cz (Act 1,
                        pooled): the GENERIC (any-onset-type) contagion
                        measure, since there is no resolution type to be
                        specific about in a pooled spec
    contagion_dist_def  Z2b(def): same construction, DEFAULT-LINKED in-crisis
                        years only (onset_def|continuation of a default-
                        linked episode)
    l_contagion_dist_def  Z2b(def) lagged -- ADOPTED in cz_def (resolution-
                        type): the default arm's classification power was
                        materially weaker under the generic contagion/recency
                        combination (08c_first_stage_table.do's diagnostics),
                        so cz_def is narrowed to default-linked-only
    contagion_dist_atdef  Z2b(atdef) NEW, ADDITIVE variant: same construction
                        as contagion_dist_def, but the donor pool k is
                        widened from the 52-country spread-crisis panel to
                        that panel UNION every country in the full Asonuma-
                        Trebesch (2016) default/restructuring database
                        (17b_merge_at_full.do's own source, re-read here
                        directly since 17b runs after this file), and the
                        donor-in-crisis flag is "this donor country-year
                        falls inside an AT-recorded default/restructuring
                        window" rather than "this donor has a default-linked
                        spread-crisis onset/continuation". NOT adopted into
                        cz_def or any AIPW/first-stage file -- built here as
                        a candidate for a future robustness test only.
    l_contagion_dist_atdef  Z2b(atdef) lagged (predetermined)
    past_onsets         Z3: cumulative own onsets before year t (proneness)
    past_def_onsets     Z3(def): cumulative own default-linked onsets before t
    years_since_def_onset  Z3(def)-recency: years since most recent prior
                        DEFAULT-linked onset -- ADOPTED in cz_def, same
                        reasoning as l_contagion_dist_def above
    years_since_onset   Z3-recency: generic (any-onset-type) version -- kept
                        built but not currently used by any adopted
                        predictor set (see 08c's diagnostic history)
===========================================================================*/

use "$clean/panel_build.dta", clear

* ── numeric panel id + xtset ────────────────────────────────────────────────
capture drop cid
egen cid = group(iso3)
label var cid "Numeric country id (from iso3)"
xtset cid year
sort cid year

* ── Contagion: leave-one-out regional onset share (Z2) ──────────────────────
* The region aggregates are restricted to carryin==0: carry-in rows are pre-EMBIG
* lag scaffolding with no spread data, so counting them would inflate the member
* denominator and dilute the share for the real observations.
* Dropped SEPARATELY, not as one combined `drop A B C D': reg_n_onset/
* reg_n_members are working variables dropped again at the end of this block
* (line ~37) and so are never saved to panel_build.dta, while reg_crisis_share/
* l_reg_crisis_share ARE saved -- so on any run after the first, this list
* mixes variables that exist with ones that don't. `drop' fails as a WHOLE
* when EITHER is missing, and `capture' then silently swallows that failure,
* leaving reg_crisis_share undropped and the `gen' below erroring "already
* defined" -- exactly what happened on the second run.
capture drop reg_n_onset
capture drop reg_n_members
capture drop reg_crisis_share
capture drop l_reg_crisis_share
bysort region year: egen reg_n_onset   = total(onset_all) if carryin==0
bysort region year: egen reg_n_members = count(cid)       if carryin==0
gen double reg_crisis_share = (reg_n_onset - onset_all) / (reg_n_members - 1) ///
    if reg_n_members > 1 & carryin==0
label var reg_crisis_share "Share of OTHER same-region countries with onset (year t)"
drop reg_n_onset reg_n_members

xtset cid year
gen double l_reg_crisis_share = L.reg_crisis_share
label var l_reg_crisis_share "Z2: lagged regional crisis share (contagion predictor)"

* ── Distance-weighted contagion (Z2b): CEPII geo-based inverse-distance sum ──
* Contagion_it = sum_k [ onset_all_kt / W_ik ],  W_ik = Dist_ik / sum_k' Dist_ik'
* Same construction as the reference paper's restructuring-contagion variable
* (their Section 4.1), applied here to ANY spread-crisis onset rather than a
* restructuring dummy specifically. Dist_ik is the great-circle distance
* between capital cities (Haversine formula on lat/lon), from CEPII's
* geo_cepii reference file (data/raw/GEO_CEPII.xlsx). Dividing by W_ik (a
* SHARE of country i's total distance to every other country) amplifies close
* neighbours and shrinks distant ones -- an inverse-distance weighting, just
* parameterised as a division rather than a multiplication by 1/distance.
* Two panel countries carry non-standard codes in this CEPII vintage:
* Romania is "ROM" (not ISO "ROU") and Serbia is "YUG" ("Serbia and
* Montenegro" -- same capital, Belgrade, so the coordinates are still valid);
* both are remapped to this project's iso3 codes below before merging.
preserve
    import excel "$raw/GEO_CEPII.xlsx", sheet(geo_cepii) firstrow clear
    * GEO_CEPII lists one row per MAJOR CITY, not one row per country: several
    * countries here have 2+ rows (a capital plus a secondary/former capital or
    * largest city), flagged by `cap' (1 = the official capital; 2 = a
    * secondary/other capital; 0 = a major city, not a capital). Filtering to
    * cap==1 BEFORE deduplicating is essential -- without it, `duplicates drop'
    * keeps whichever row the sheet happens to list first, which for this
    * panel is the WRONG city for four countries: Bolivia (La Paz, cap=2,
    * ahead of Sucre, cap=1), Brazil (Sao Paulo, cap=0, ahead of Brasilia),
    * Nigeria (Lagos, cap=2, ahead of Abuja), and Turkey (Istanbul, cap=0,
    * ahead of Ankara) -- confirmed by inspecting the raw file directly.
    keep if cap == 1
    keep iso3 lat lon
    destring lat lon, replace force
    replace iso3 = "ROU" if iso3=="ROM"
    replace iso3 = "SRB" if iso3=="YUG"
    duplicates drop iso3, force
    tempfile geo
    save `geo'
restore

* Distinct panel countries with lat/lon -- carryin==0, matching
* reg_crisis_share's own restriction (carry-in rows are pre-EMBIG scaffolding,
* not real countries to weight contagion against).
preserve
    keep if carryin==0
    contract iso3
    merge 1:1 iso3 using `geo', keep(match master) nogen
    quietly count if missing(lat)
    if r(N) > 0 {
        di as error "  ** GEO_CEPII: `r(N)' panel countries have no lat/lon match -- contagion_dist will be missing for them:"
        list iso3 if missing(lat)
    }
    tempfile ctylist
    save `ctylist'
restore

* All ordered pairs (i,k), i != k, and their Haversine distance + weight.
* Time-invariant (geography doesn't change), so this is built once, not per year.
preserve
    use `ctylist', clear
    rename iso3 iso3_k
    rename lat lat_k
    rename lon lon_k
    tempfile klist
    save `klist'

    use `ctylist', clear
    rename iso3 iso3_i
    rename lat lat_i
    rename lon lon_i
    cross using `klist'
    drop if iso3_i == iso3_k
    drop if missing(lat_i) | missing(lat_k)

    * Great-circle distance via the spherical law of cosines, in km (Earth
    * radius 6371 km). The acos() argument is clipped to [-1,1]: floating-point
    * rounding can push it fractionally outside that range for very close
    * points, which would otherwise return a missing distance.
    gen double _arg = sin(lat_i*c(pi)/180)*sin(lat_k*c(pi)/180) ///
        + cos(lat_i*c(pi)/180)*cos(lat_k*c(pi)/180)*cos((lon_k-lon_i)*c(pi)/180)
    replace _arg = min(1, max(-1, _arg))
    gen double dist_ik = 6371 * acos(_arg)
    drop _arg

    * NORMALIZED INVERSE-DISTANCE WEIGHT, summing to 1 across every donor k
    * for a given i -- so that a weighted sum of a 0/1 donor flag below is a
    * bounded WEIGHTED SHARE (in [0,1]), not an unbounded quantity. w_ik is
    * proportional to 1/dist_ik (a closer country gets more weight), then
    * divided by the sum of every country's inverse distance so the weights
    * sum to exactly 1 for each i, the same way reg_crisis_share's flat
    * regional weights (1/(reg_n_members-1) each) sum to 1 -- this is that
    * same "share of X" construction, with a distance-based weight in place
    * of a flat same-region one, rather than an unbounded weighted sum.
    gen double _invdist_ik = 1 / dist_ik
    bysort iso3_i: egen double _suminvdist_i = total(_invdist_ik)
    gen double w_ik = _invdist_ik / _suminvdist_i
    keep iso3_i iso3_k w_ik
    tempfile weights
    save `weights'
restore

* Join the (i,k) weights to every (k,t) onset, then collapse to (i,t):
* contagion_dist_it = sum over k of in_crisis_kt * w_ik, where in_crisis_kt =
* onset_all|continuation -- a STOCK of regional distress (every year another
* country IS in a crisis, not just the year its crisis STARTED). An earlier
* version of this variable summed onset_all alone (a shock indicator: a
* neighbor's crisis contributed only in its own onset year, then contributed
* nothing in every subsequent year it remained in crisis, even though it was
* still visibly in distress) -- corrected here to measure "how much distress
* currently surrounds country i", the intended meaning of this predictor.
* Since w_ik sums to 1 across every donor k, this is a distance-WEIGHTED
* SHARE bounded in [0,1] -- a country surrounded entirely by in-crisis
* neighbors approaches 1, one with no in-crisis neighbors is 0 -- matching
* the scale of a genuine "share of peers in crisis" measure (the same
* bounded scale as reg_crisis_share, and as the reference paper's own
* contagion predictor: Table B3 reports it as mean 0.05, range [0,0.88]).
preserve
    keep if carryin==0
    gen byte donor_in_crisis = (onset_all==1 | continuation==1)
    keep iso3 year donor_in_crisis
    rename iso3 iso3_k
    tempfile donors
    save `donors'

    use `weights', clear
    joinby iso3_k using `donors'
    gen double _contrib = donor_in_crisis * w_ik
    collapse (sum) contagion_dist = _contrib, by(iso3_i year)
    rename iso3_i iso3
    label var contagion_dist "Z2b: distance-weighted SHARE of OTHER countries in-crisis (onset|continuation, year t), CEPII great-circle, bounded [0,1]"
    tempfile contagion
    save `contagion'
restore

* Dropped before the merge/gen below, not combined -- same reasoning as
* reg_crisis_share/past_onsets above: both persist to panel_build.dta, so a
* second run would otherwise hit "already defined" (merge conflict on
* contagion_dist, then gen conflict on l_contagion_dist).
capture drop contagion_dist
capture drop l_contagion_dist
merge m:1 iso3 year using `contagion', keep(master match) nogen
xtset cid year
gen double l_contagion_dist = L.contagion_dist
label var l_contagion_dist "Z2b: lagged distance-weighted contagion share (in-crisis stock, any type, bounded [0,1]), predetermined"

* ── DEFAULT-LINKED contagion stock (Z2b-def): ADOPTED for cz_def ────────────
* Same "in-crisis stock" construction as contagion_dist above, but the donor
* flag is restricted to a country's DEFAULT-LINKED in-crisis years (its own
* episode's resolution type, not just its onset year). 08c_first_stage_table.do
* found the DEFAULT-linked arm's classification power essentially disappears
* under the generic contagion_dist/years_since_onset combination (chi2(pred)
* p=.101, roccomp p=.525 -- barely distinguishable from the controls-only
* model), unlike the non-default arm (p<.001, roccomp p=.039). Narrowing both
* predictors to default-linked-only is tested here as the fix for that arm
* specifically.
*
* The resolution-type fill (nd_ep in 18_transforms.do) does not exist yet at
* this stage of the build, so it is reproduced locally here, restricted to
* carryin==0 exactly like the donors block above: ep_seq_tmp is a running
* onset counter per country, nd_ep_tmp is nondefault's episode fill over that
* counter -- identical logic to 18_transforms.do's own construction, just
* computed early because this predictor is needed before that file runs.
preserve
    keep if carryin==0
    sort cid year
    bysort cid (year): gen byte ep_seq_tmp = sum(onset_all)
    bysort cid ep_seq_tmp: egen byte nd_ep_tmp = max(nondefault)
    gen byte donor_in_crisis_def = (onset_all==1 | continuation==1) & nd_ep_tmp==0
    keep iso3 year donor_in_crisis_def
    rename iso3 iso3_k
    tempfile donors_def
    save `donors_def'

    use `weights', clear
    joinby iso3_k using `donors_def'
    gen double _contrib = donor_in_crisis_def * w_ik
    collapse (sum) contagion_dist_def = _contrib, by(iso3_i year)
    rename iso3_i iso3
    label var contagion_dist_def "Z2b(def): distance-weighted SHARE of OTHER countries DEFAULT-linked in-crisis, CEPII great-circle, bounded [0,1]"
    tempfile contagion_def
    save `contagion_def'
restore

capture drop contagion_dist_def
capture drop l_contagion_dist_def
merge m:1 iso3 year using `contagion_def', keep(master match) nogen
sort cid year
xtset cid year
gen double l_contagion_dist_def = L.contagion_dist_def
label var l_contagion_dist_def "Z2b(def): lagged distance-weighted contagion share (default-linked in-crisis stock, bounded [0,1]), predetermined"

quietly summarize year
local panel_yr_min = r(min)
local panel_yr_max = r(max)

* ── DEFAULT-LINKED contagion stock, AT-DATABASE-WIDE DONOR POOL (Z2b-atdef
*    variant): NEW, ADDITIVE, NOT adopted into cz_def ───────────────────────
* contagion_dist_atdef / l_contagion_dist_atdef -- same distance-weighted-
* share construction as contagion_dist_def above (same Haversine/CEPII
* weight machinery), but the donor pool k widens from "this project's
* 52-country spread-crisis panel" to that panel UNION every country in the
* full Asonuma-Trebesch (2016) default/restructuring database
* (17b_merge_at_full.do's own source, re-imported here directly because 17b
* runs AFTER this file in the build chain -- see 00_master.do -- so its
* at_default_year column does not exist yet at this point). Recipient
* country i stays restricted to the spread-crisis panel (carryin==0), same
* as contagion_dist_def -- only the donor side widens, since this predictor
* is only ever needed for the 52-country onset-tier propensity model.
*
* DONOR-IN-CRISIS FLAG IS A UNION OF TWO SOURCES, not the AT window alone,
* because the two sources cover different, only PARTIALLY overlapping sets
* of events: a spread-panel donor can be mid-spread-crisis in a year that
* falls outside its own AT restructuring window (e.g. the crisis year
* precedes the eventual restructuring, or the restructuring window as coded
* by AT is narrower than the spread episode), and conversely an AT-recorded
* default/restructuring can fall in a year this project's own spread dating
* does not mark as onset/continuation (e.g. a quiet debt-service suspension
* with no EMBIG spread spike). Using either source alone would silently drop
* donor-in-crisis years the other source captures. The flag is therefore:
*   donor_in_crisis_atdef = 1 if EITHER
*     (a) the donor country-year falls inside an AT-recorded default/
*         restructuring window (any raw AT case row, start year through end
*         year -- see below), OR
*     (b) the donor country-year is a default-linked spread-crisis year per
*         this project's OWN existing dating (onset_all==1 | continuation==1,
*         carryin==0) -- i.e. the identical donor_in_crisis flag the
*         contagion_dist_def block above already uses.
* For a spread-panel donor, both (a) and (b) are available and OR'd. For an
* AT-only donor (no spread-panel presence at all), (b) is structurally
* unavailable (never tested for a spread crisis -- see 10b_skeleton_atonly.do's
* header), so the flag is determined by (a) alone.
*
* Window mechanics: raw AT case rows are used directly here (not 17b's own
* same-country/same-year COLLAPSED episodes) because collapsing changes which
* year-range GROUP similarly-dated cases fall into, not which years are
* covered -- the UNION of windows a country's raw cases span is identical
* either way for a pure in/out-of-window 0/1 flag, so the extra collapsing
* step (and the type/severity logic that goes with it) is unneeded work for
* this specific use; contagion_dist_atdef never reads AT's preemptive/
* post-default type, only whether a window is open at all.
*
* ASSUMPTION FLAGGED: the same user-maintained post-2020 supplement 17b
* carries (Zambia, Ghana, Sri Lanka, Lebanon, Ukraine -- see 17b_merge_at_
* full.do's header, Section 2) is reproduced here so the two files' AT-
* sourced windows agree. If that supplement is ever edited in 17b, this
* block must be updated to match by hand -- it cannot read 17b's in-file
* locals, since 17b has not run yet at this point in the build chain.
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
    replace iso3 = "ROU" if iso3 == "ROM"
    keep if !missing(iso3)

    gen int start_year = year(start_date) if !missing(start_date)
    gen int end_year = year(end_date) if !missing(end_date)
    replace end_year = start_year if missing(end_year)
    drop if missing(start_year)

    keep iso3 start_year end_year
    tempfile at_windows_raw
    save `at_windows_raw'
restore

preserve
    clear
    gen str3 iso3 = ""
    gen int start_year = .
    gen int end_year = .
    set obs 5
    * Same 5 post-2020 events, same dates, as 17b_merge_at_full.do's own
    * user-maintained supplement -- see that file's header for sourcing.
    replace iso3 = "ZMB" in 1
    replace start_year = 2020 in 1
    replace end_year   = 2024 in 1
    replace iso3 = "GHA" in 2
    replace start_year = 2022 in 2
    replace end_year   = 2024 in 2
    replace iso3 = "LKA" in 3
    replace start_year = 2022 in 3
    replace end_year   = 2024 in 3
    replace iso3 = "LBN" in 4
    replace start_year = 2020 in 4
    replace end_year   = 2026 in 4
    replace iso3 = "UKR" in 5
    replace start_year = 2022 in 5
    replace end_year   = 2026 in 5
    append using `at_windows_raw'
    rename iso3 iso3_k
    tempfile at_windows_k
    save `at_windows_k'
restore

* Union donor pool: AT-database countries (distinct iso3 in at_windows_k)
* UNION the existing 52-country spread panel (`ctylist', already built
* above), matched against the FULL (unrestricted) `geo' file for lat/lon.
preserve
    use `at_windows_k', clear
    rename iso3_k iso3
    contract iso3
    drop _freq
    merge 1:1 iso3 using `geo', keep(match master) nogen
    tempfile at_ctylist
    save `at_ctylist'
restore

preserve
    use `ctylist', clear
    quietly count if !missing(lat)
    local n_narrow = r(N)
    append using `at_ctylist'
    duplicates drop iso3, force
    quietly count if missing(lat)
    if r(N) > 0 {
        di as error "  ** GEO_CEPII: `r(N)' AT-database donor countries have no lat/lon match -- excluded from contagion_dist_atdef's donor pool:"
        list iso3 if missing(lat)
    }
    drop if missing(lat)
    tempfile ctylist_wide
    save `ctylist_wide'
    quietly count
    local n_wide = r(N)
    di as result "  contagion_dist_atdef: donor pool widened from `n_narrow' (spread-panel-only) to `n_wide' (spread panel UNION full AT database) countries -- +`=`n_wide'-`n_narrow'' donors."
restore

* Weight matrix for the WIDENED donor pool: i restricted to the 52-country
* spread panel (unchanged recipient side, matching contagion_dist_def
* above), k drawn from `ctylist_wide'. Same Haversine/normalization logic as
* the weights block near the top of this file, just re-run against the
* wider k set -- `weights' above was built i x k over the narrow 52-country
* pool only and cannot be reused here.
preserve
    use `ctylist_wide', clear
    rename iso3 iso3_k
    rename lat lat_k
    rename lon lon_k
    tempfile klist_wide
    save `klist_wide'

    use `ctylist', clear
    rename iso3 iso3_i
    rename lat lat_i
    rename lon lon_i
    cross using `klist_wide'
    drop if iso3_i == iso3_k
    drop if missing(lat_i) | missing(lat_k)

    gen double _arg = sin(lat_i*c(pi)/180)*sin(lat_k*c(pi)/180) ///
        + cos(lat_i*c(pi)/180)*cos(lat_k*c(pi)/180)*cos((lon_k-lon_i)*c(pi)/180)
    replace _arg = min(1, max(-1, _arg))
    gen double dist_ik = 6371 * acos(_arg)
    drop _arg

    gen double _invdist_ik = 1 / dist_ik
    bysort iso3_i: egen double _suminvdist_i = total(_invdist_ik)
    gen double w_ik = _invdist_ik / _suminvdist_i
    keep iso3_i iso3_k w_ik
    tempfile weights_wide
    save `weights_wide'
restore

* Donor-in-crisis flag (a): AT window, over the widened pool x every year in
* the panel's own span. Country-year x AT-window candidate matches via
* joinby (many-to-many on iso3_k), filtered to in-window rows, reduced to a
* 0/1 per country-year (two AT cases overlapping for one country in one year
* both flag the same year -- max, not sum, is the right reduction here).
preserve
    use `ctylist_wide', clear
    keep iso3
    rename iso3 iso3_k
    local nyr = `panel_yr_max' - `panel_yr_min' + 1
    expand `nyr'
    bysort iso3_k: gen int year = `panel_yr_min' + _n - 1
    tempfile donor_years
    save `donor_years'

    joinby iso3_k using `at_windows_k'
    gen byte _in_window = inrange(year, start_year, end_year)
    keep if _in_window==1
    gen byte flag_atwindow = 1
    collapse (max) flag_atwindow, by(iso3_k year)
    tempfile at_flag
    save `at_flag'
restore

* Donor-in-crisis flag (b): this project's own spread-crisis dating
* (onset_all==1 | continuation==1, carryin==0) -- identical source to
* contagion_dist_def's donor_in_crisis above, available only for donors that
* are themselves spread-panel countries.
preserve
    keep if carryin==0
    gen byte flag_spread = (onset_all==1 | continuation==1)
    keep iso3 year flag_spread
    rename iso3 iso3_k
    tempfile spread_flag
    save `spread_flag'
restore

* Union: OR the two flags together. Left-joining flag_spread onto the full
* donor-year skeleton leaves AT-only donors (never in the spread panel) with
* flag_spread missing, treated as 0 -- exactly the "(b) structurally
* unavailable, (a) alone determines the flag" case described above.
preserve
    use `donor_years', clear
    merge 1:1 iso3_k year using `at_flag', keep(master match) nogen
    replace flag_atwindow = 0 if missing(flag_atwindow)
    merge 1:1 iso3_k year using `spread_flag', keep(master match) nogen
    replace flag_spread = 0 if missing(flag_spread)
    gen byte donor_in_crisis_atdef = max(flag_atwindow, flag_spread)
    keep iso3_k year donor_in_crisis_atdef
    tempfile donors_atdef
    save `donors_atdef'
restore

* Join the widened (i,k) weights to the union donor-in-crisis flag, then
* collapse to (i,t) -- same mechanics as contagion_dist/contagion_dist_def
* above: contagion_dist_atdef_it = sum over k of donor_in_crisis_atdef_kt *
* w_ik, a distance-weighted SHARE bounded in [0,1] since w_ik sums to 1
* across every donor k for a given i.
preserve
    use `weights_wide', clear
    joinby iso3_k using `donors_atdef'
    gen double _contrib = donor_in_crisis_atdef * w_ik
    collapse (sum) contagion_dist_atdef = _contrib, by(iso3_i year)
    rename iso3_i iso3
    label var contagion_dist_atdef "Z2b(atdef) NEW/ADDITIVE variant: distance-weighted SHARE of OTHER countries in default/restructuring crisis (year t) -- donor-in-crisis flag is the UNION of an AT-recorded default/restructuring window and this project's own spread-crisis dating; donor pool widened to the full AT database UNION the 52-country spread panel; bounded [0,1]. NOT adopted in cz_def."
    tempfile contagion_atdef
    save `contagion_atdef'
restore

capture drop contagion_dist_atdef
capture drop l_contagion_dist_atdef
merge m:1 iso3 year using `contagion_atdef', keep(master match) nogen
sort cid year
xtset cid year
gen double l_contagion_dist_atdef = L.contagion_dist_atdef
label var l_contagion_dist_atdef "Z2b(atdef) NEW/ADDITIVE: lagged distance-weighted contagion share (AT-database-wide donor pool, union of AT-window and spread-crisis donor-in-crisis flags), predetermined. NOT adopted in any cz/cz_def -- additive robustness variant only."

* Diagnostic comparison: does the widened predictor actually differ from
* contagion_dist_def, restricted to the default-linked arm's onset rows
* (the sample this variant would be tested on in 08c_first_stage_table.do)?
di as result _n "  contagion_dist_atdef diagnostic (default-linked onset rows, carryin==0 & onset_def==1):"
quietly summarize contagion_dist_def if carryin==0 & onset_def==1
di as result "    contagion_dist_def:   mean=" %6.4f r(mean) "  N=" r(N)
quietly summarize contagion_dist_atdef if carryin==0 & onset_def==1
di as result "    contagion_dist_atdef: mean=" %6.4f r(mean) "  N=" r(N)
quietly corr contagion_dist_def contagion_dist_atdef if carryin==0 & onset_def==1
di as result "    corr(contagion_dist_def, contagion_dist_atdef) on those rows: " %6.4f r(rho) "  N=" r(N)


* ── Proneness: cumulative own onsets, lagged (Z3) ───────────────────────────
* Dropped SEPARATELY -- same reasoning as reg_crisis_share above: cum_onset/
* cum_def are working variables dropped again below and never saved, while
* past_onsets/past_def_onsets ARE saved, so a combined drop fails as a whole
* on any run after the first.
capture drop cum_onset
capture drop past_onsets
capture drop cum_def
capture drop past_def_onsets
bysort cid (year): gen cum_onset = sum(onset_all)
gen past_onsets = L.cum_onset
replace past_onsets = 0 if missing(past_onsets)
label var past_onsets "Z3: own onsets before year t (proneness predictor)"
drop cum_onset

bysort cid (year): gen cum_def = sum(onset_def)
gen past_def_onsets = L.cum_def
replace past_def_onsets = 0 if missing(past_def_onsets)
label var past_def_onsets "Z3(def): own default-linked onsets before year t"
drop cum_def

* ── Recency alternative to past_def_onsets: years since the most recent PRIOR
* default-linked onset (Z3(def)-recency). past_def_onsets is a running COUNT
* that never resets, so for a serial defaulter it behaves close to a
* permanent country identifier rather than a genuine time-varying predictor —
* diagnosed in 24_aipw_channels_flow.do's Section 1a as the main driver of
* severe weight concentration in the flow AIPW's def-arm propensity model
* (top 5% of rows = 98.9% of the AIPW summand's variance). Section 1b there
* tested this recency measure as a replacement: individually significant
* (z=-2.84, p=.005, credit h=1 specification) with the economically sensible
* sign (more years since the last default, lower probability of a new one),
* and adopted as the active predictor across the flow files on that basis —
* NOT because it reduces the weight-concentration problem, which Section 1b
* also found it does NOT meaningfully fix (98.9% -> 98.6%, essentially
* unchanged). See 21_aipw_flow.do's header for the full adoption note.
*
* Construction: for an ONSET row, this necessarily refers to an EARLIER,
* distinct episode (an onset row cannot be its own prior onset), so it is not
* circular the way epc_X would be. For a TRANQUIL row it is unambiguous prior
* history. Countries with no PRIOR default-linked onset are censored to 50
* (safely beyond this panel's ~35-year span), following past_def_onsets' own
* convention of replacing missing with a fixed value rather than dropping rows.
* Dropped SEPARATELY, same reasoning: _defyear/_defyear_lag never persist to
* panel_build.dta (dropped below), years_since_def_onset does.
capture drop _defyear
capture drop _defyear_lag
capture drop years_since_def_onset
gen _defyear = year if onset_def==1
bysort cid (year): replace _defyear = _defyear[_n-1] if missing(_defyear) & _n>1
bysort cid (year): gen _defyear_lag = _defyear[_n-1]
gen double years_since_def_onset = year - _defyear_lag if !missing(_defyear_lag)
replace years_since_def_onset = 50 if missing(years_since_def_onset)
label var years_since_def_onset "Z3(def) recency: years since most recent PRIOR default-linked onset (censored at 50). ADOPTED in cz_def."
drop _defyear _defyear_lag

* ── Generic (any-onset-type) recency (Z3-recency) ───────────────────────────
* Same construction as years_since_def_onset, but the clock resets on ANY
* prior onset (onset_all), not just a default-linked one. Built and tested
* against years_since_def_onset in 08c_first_stage_table.do; not currently
* used in any adopted predictor set (cz_def is default-linked-only -- see
* years_since_def_onset above), kept built for reference/future comparison.
capture drop _anyyear
capture drop _anyyear_lag
capture drop years_since_onset
gen _anyyear = year if onset_all==1
bysort cid (year): replace _anyyear = _anyyear[_n-1] if missing(_anyyear) & _n>1
bysort cid (year): gen _anyyear_lag = _anyyear[_n-1]
gen double years_since_onset = year - _anyyear_lag if !missing(_anyyear_lag)
replace years_since_onset = 50 if missing(years_since_onset)
label var years_since_onset "Z3-recency: years since most recent PRIOR onset of any type (censored at 50)"
drop _anyyear _anyyear_lag

save "$clean/panel_build.dta", replace

di as result _n "17_predictors.do complete."
foreach v in l_reg_crisis_share l_contagion_dist l_contagion_dist_def l_contagion_dist_atdef past_onsets past_def_onsets years_since_def_onset years_since_onset {
    quietly count if !missing(`v') & sample_base==1
    di as result "  `v': `r(N)' non-missing sample rows"
}
