/*===========================================================================
  17_PREDICTORS.DO   —  FROM-SCRATCH REBUILD, STAGE 7  (first-stage predictors Z)
  Derived exclusion-restriction predictors (no raw import). Mirrors the old 01e.
  Also creates the numeric panel id (cid) and xtset — the first stage that needs
  time-series operators.

    cid                 numeric country id (from iso3), + xtset cid year
    reg_crisis_share    leave-one-out share of OTHER same-region countries with
                        an onset in year t
    l_reg_crisis_share  Z2: lagged regional crisis share (contagion)
    contagion_dist      Z2b: distance-weighted sum of OTHER countries' onsets
                        (year t), CEPII great-circle distance (GEO_CEPII.xlsx)
    l_contagion_dist    Z2b lagged (predetermined)
    contagion_def_dist  Z2c: same, but donors restricted to default-linked
                        onsets (onset_def) only
    l_contagion_def_dist Z2c lagged (predetermined)
    past_onsets         Z3: cumulative own onsets before year t (proneness)
    past_def_onsets     Z3(def): cumulative own default-linked onsets before t
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

    bysort iso3_i: egen double _sumdist_i = total(dist_ik)
    gen double w_ik = dist_ik / _sumdist_i
    keep iso3_i iso3_k w_ik
    tempfile weights
    save `weights'
restore

* Join the (i,k) weights to every (k,t) onset, then collapse to (i,t):
* contagion_dist_it = sum over k of onset_all_kt / w_ik.
preserve
    keep if carryin==0
    keep iso3 year onset_all
    rename iso3 iso3_k
    tempfile donors
    save `donors'

    use `weights', clear
    joinby iso3_k using `donors'
    gen double _contrib = onset_all / w_ik
    collapse (sum) contagion_dist = _contrib, by(iso3_i year)
    rename iso3_i iso3
    label var contagion_dist "Z2b: distance-weighted sum of OTHER countries' onset_all (year t), CEPII great-circle"
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
label var l_contagion_dist "Z2b: lagged distance-weighted contagion (any onset), predetermined"

* ── Distance-weighted contagion, DEFAULT-LINKED ONSETS ONLY (Z2c) ───────────
* Same construction and weights (`weights', built above) as contagion_dist,
* but the donor numerator is restricted to onset_def rather than onset_all.
* This is closer to the reference paper's own contagion variable, which is
* built from restructuring onsets specifically, not any spread-crisis onset —
* onset_all above is this project's generalisation, not theirs. Built here as
* a SEPARATE candidate predictor, tested against l_contagion_dist in
* 21b_first_stage_table_flow.do's AUROC comparison before any decision about
* replacing cz_recency's l_contagion_dist with it.
preserve
    keep if carryin==0
    keep iso3 year onset_def
    rename iso3 iso3_k
    tempfile donors_def
    save `donors_def'

    use `weights', clear
    joinby iso3_k using `donors_def'
    gen double _contrib_def = onset_def / w_ik
    collapse (sum) contagion_def_dist = _contrib_def, by(iso3_i year)
    rename iso3_i iso3
    label var contagion_def_dist "Z2c: distance-weighted sum of OTHER countries' onset_def (year t), CEPII great-circle"
    tempfile contagion_def
    save `contagion_def'
restore

capture drop contagion_def_dist
capture drop l_contagion_def_dist
merge m:1 iso3 year using `contagion_def', keep(master match) nogen
xtset cid year
gen double l_contagion_def_dist = L.contagion_def_dist
label var l_contagion_def_dist "Z2c: lagged distance-weighted contagion (default-linked onset only), predetermined"

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
label var years_since_def_onset "Z3(def) recency: years since most recent PRIOR default-linked onset (censored at 50)"
drop _defyear _defyear_lag

save "$clean/panel_build.dta", replace

di as result _n "17_predictors.do complete."
foreach v in l_reg_crisis_share l_contagion_dist l_contagion_def_dist past_onsets past_def_onsets years_since_def_onset {
    quietly count if !missing(`v') & sample_base==1
    di as result "  `v': `r(N)' non-missing sample rows"
}
