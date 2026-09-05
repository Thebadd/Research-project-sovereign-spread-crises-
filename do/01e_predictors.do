/*===========================================================================
  01E_PREDICTORS.DO
  Exclusion-restriction PREDICTORS for the first-stage probit (Jordà–Taylor
  2016 style), adapted to sovereign spread crises.

  A valid excluded predictor is (i) relevant — predicts crisis onset;
  (ii) excludable — affects output only through the crisis, given controls;
  (iii) predetermined — measured at t-1.

  Z1 — Global push = US 10y (ust10y). Already in the panel. Relevance: push
       factor for EM spreads. Exclusion: enforced by year FE in the LP (a pure
       time-series variable is absorbed by i.year in the second stage), so it
       cannot enter the outcome equation. Used directly in the probit; nothing
       to build here.

  Z2 — Regional contagion = regional crisis share (built here). For country i in
       year t, the share of the OTHER countries in i's region that had a crisis
       onset in year t. Equal-weight within region (no lat/long for distance
       weighting). The one-year lag (l_reg_crisis_share) is the predictor.

  Z3 — Crisis proneness = past-onset count (built here, robustness only). The
       cumulative number of country i's own onsets strictly before year t.
       Weaker exclusion (output hysteresis) → robustness, not baseline.

  Exclusion discipline: Z2 and Z3 are country-varying, so they must enter the
  probit ONLY and be omitted from the LP/AIPW second stage (Z1 is auto-excluded
  by year FE).

  Run AFTER 01d_merge_vulnerability.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* Z2 — REGIONAL CRISIS SHARE (leave-one-out, same region, same year)
* ══════════════════════════════════════════════════════════════════════════
* region grouping: `region' (string) is available; rid is its numeric id.
capture drop reg_n_onset reg_n_members reg_crisis_share l_reg_crisis_share

* Region-year totals over ALL countries with an observation that year
bysort region year: egen reg_n_onset   = total(onset_all)
bysort region year: egen reg_n_members = count(cid)

* Leave-one-out share: exclude country i's own onset and own membership
gen reg_crisis_share = (reg_n_onset - onset_all) / (reg_n_members - 1) ///
    if reg_n_members > 1
label var reg_crisis_share "Share of OTHER same-region countries with onset (year t)"

* Restore the panel sort (bysort region year above re-sorted the data, which
* breaks the L. time-series operator — hence the r(5) "not sorted" error).
sort cid year
xtset cid year

* Predictor = one-year lag (last year's regional crises predict this year's onset)
gen l_reg_crisis_share = L.reg_crisis_share
label var l_reg_crisis_share "Z2: lagged regional crisis share (contagion predictor)"

drop reg_n_onset reg_n_members

* ══════════════════════════════════════════════════════════════════════════
* Z3 — PAST-ONSET COUNT (own crisis history strictly before year t)
* ══════════════════════════════════════════════════════════════════════════
capture drop cum_onset past_onsets

* Running total of onsets including the current year, then lag by one so the
* count refers to onsets strictly before t (predetermined).
bysort cid (year): gen cum_onset = sum(onset_all)
gen past_onsets = L.cum_onset
replace past_onsets = 0 if missing(past_onsets)   // first observed year: 0 prior onsets
label var past_onsets "Z3: number of own onsets before year t (proneness predictor)"

drop cum_onset

* Z3(def) — count of own DEFAULT-LINKED onsets strictly before year t.
* Act 2's proneness predictor: default-repeat predicts default-type resolution.
capture drop cum_def_onset past_def_onsets
bysort cid (year): gen cum_def_onset = sum(onset_def)
gen past_def_onsets = L.cum_def_onset
replace past_def_onsets = 0 if missing(past_def_onsets)
label var past_def_onsets "Z3(def): number of own default-linked onsets before year t"

drop cum_def_onset

* ══════════════════════════════════════════════════════════════════════════
* SAVE + COVERAGE / RELEVANCE PEEK
* ══════════════════════════════════════════════════════════════════════════
sort cid year
xtset cid year
save "$clean/panel_lp.dta", replace

di as result _n "=== PREDICTOR COVERAGE (sample==1) ==="
foreach v in fedfunds l_reg_crisis_share past_onsets past_def_onsets {
    quietly count if sample == 1 & !missing(`v')
    di as result "  `v': " r(N) " non-missing sample obs"
}

di as result _n "=== PREDICTOR MEANS BY ONSET STATUS (sample==1) — relevance peek ==="
di as result "  (predictors should differ between onset and tranquil obs)"
foreach v in fedfunds l_reg_crisis_share past_onsets past_def_onsets {
    quietly summarize `v' if sample==1 & onset_all==1
    local m1 = r(mean)
    quietly summarize `v' if sample==1 & onset_all==0
    local m0 = r(mean)
    di as result "  `v': onset=" %7.3f `m1' "   tranquil=" %7.3f `m0'
}

di as result _n "01e_predictors.do complete — panel_lp.dta updated (Z2, Z3 built; Z1 = ust10y)."
