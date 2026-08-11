/*===========================================================================
  17_PREDICTORS.DO   —  FROM-SCRATCH REBUILD, STAGE 7  (first-stage predictors Z)
  Derived exclusion-restriction predictors (no raw import). Mirrors the old 01e.
  Also creates the numeric panel id (cid) and xtset — the first stage that needs
  time-series operators.

    cid                 numeric country id (from iso3), + xtset cid year
    reg_crisis_share    leave-one-out share of OTHER same-region countries with
                        an onset in year t
    l_reg_crisis_share  Z2: lagged regional crisis share (contagion)
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
capture drop reg_n_onset reg_n_members reg_crisis_share l_reg_crisis_share
bysort region year: egen reg_n_onset   = total(onset_all) if carryin==0
bysort region year: egen reg_n_members = count(cid)       if carryin==0
gen double reg_crisis_share = (reg_n_onset - onset_all) / (reg_n_members - 1) ///
    if reg_n_members > 1 & carryin==0
label var reg_crisis_share "Share of OTHER same-region countries with onset (year t)"
drop reg_n_onset reg_n_members

xtset cid year
gen double l_reg_crisis_share = L.reg_crisis_share
label var l_reg_crisis_share "Z2: lagged regional crisis share (contagion predictor)"

* ── Proneness: cumulative own onsets, lagged (Z3) ───────────────────────────
capture drop cum_onset past_onsets cum_def past_def_onsets
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

save "$clean/panel_build.dta", replace

di as result _n "17_predictors.do complete."
foreach v in l_reg_crisis_share past_onsets past_def_onsets {
    quietly count if !missing(`v') & sample_base==1
    di as result "  `v': `r(N)' non-missing sample rows"
}
