/*===========================================================================
  01_BUILD_PANEL.DO
  - Import main dataset (CSV)
  - Merge default/non-default classification from Episode_Summary (Excel)
  - Generate LP outcome variables (cumulative log GDP changes)
  - xtset and save as panel_lp.dta
===========================================================================*/

clear all

* ══════════════════════════════════════════════════════════════════════════
* 1. IMPORT MAIN DATASET
* ══════════════════════════════════════════════════════════════════════════

import delimited "$raw/FINAL_DATASET_STATA_V2.csv", ///
    varnames(1) stringcols(2 3 9 10 37 38) clear

* Rename to clean Stata names
* import delimited lowercases all headers; only rename where new name differs
rename country_id        cid_str
rename region_id         rid
rename d_it              onset_all
rename d_it_x_frontier   onset_frontier
rename excl_continuation continuation
rename episode_status    ep_status
rename episode_id        ep_id
rename crit1_flag        crit1
rename crit2_flag        crit2
rename log_gdppc         ln_gdppc
rename gdppc_real_lcu    gdppc_lcu
rename gdp_real_growth   gdpg
rename l1_gdp_growth     l1_gdpg
rename l2_gdp_growth     l2_gdpg
rename spread_max        spr_max
rename spread_mean       spr_mean
rename l_spread_max      l_spr_max
rename l_spread_mean     l_spr_mean
rename ca_gdp            ca
rename debt_gdp          debt
rename primary_balance_gdp pb
rename inflation_cpi     infl
rename investment_gdp    inv
rename tot_index         tot
rename tot_pct           dtot
rename credit_gdp        credit
rename fed_funds         fedfunds
rename treasury_10y      ust10y
rename imf_program       imf
rename onset_month_num   onset_month
rename onset_half        onset_h2
rename fdi_gdp           fdi
rename govexp_gdp        govexp
rename claims_govt_gdp   claims_govt

* Destring numeric variables
destring cid_str rid year onset_all onset_frontier continuation ///
         crit1 crit2 crisis_any frontier ///
         ln_gdppc gdppc_lcu gdpg l1_gdpg l2_gdpg ///
         spr_max spr_mean l_spr_max l_spr_mean ///
         ca debt pb infl inv tot dtot credit ///
         vix fedfunds ust10y imf onset_month onset_h2 ///
         fdi govexp claims_govt, replace force

* Generate numeric country ID from string (for xtset)
encode country, gen(cid)

* Sort
sort cid year

* ══════════════════════════════════════════════════════════════════════════
* 2. MERGE DEFAULT / NON-DEFAULT CLASSIFICATION
*    Source: Episode_Summary sheet in Excel
*    Key:    country × onset_year
* ══════════════════════════════════════════════════════════════════════════

preserve
    * Use column letters (A, C, H) to avoid header-naming ambiguity.
    * Row 2 = headers, Row 3+ = data → cellrange(A3) skips both note and headers.
    * Columns: A=Country, C=Onset Year, H=Classification
    import excel "$raw/EM_Spread_Crisis_DB_FINAL.xlsx", ///
        sheet("Episode_Summary") cellrange(A3) clear

    keep A C H
    rename A country
    rename C year
    rename H classification
    destring year, replace force

    * Drop legend rows at bottom of sheet
    drop if missing(country) | missing(year) | missing(classification)
    drop if !inlist(classification, "NON-DEFAULT Spread Crisis", "Default-Linked")

    * Binary flag: 1 = non-default, 0 = default-linked
    gen     nondefault = (classification == "NON-DEFAULT Spread Crisis")

    * Keep only onset years
    keep country year nondefault classification
    duplicates drop country year, force
    tempfile classif
    save `classif'
restore

merge m:1 country year using `classif', keep(master match) nogen

* ── Manual reclassification ──────────────────────────────────────────────
* Venezuela 2008: restructuring began in 2017 (9-year lag). Episodes where
* onset precedes default by more than 5 years are reclassified as non-default
* because the spread crisis reflected genuine non-default distress at inception.
replace nondefault     = 1  if country == "Venezuela" & year == 2008
replace classification = "NON-DEFAULT Spread Crisis" if country == "Venezuela" & year == 2008

* For non-onset years: set nondefault to missing (not applicable)
replace nondefault      = . if onset_all == 0
replace classification  = "" if onset_all == 0

* Derive two onset dummies for Act 2
gen onset_nd  = (onset_all == 1 & nondefault == 1)   // non-default onsets
gen onset_def = (onset_all == 1 & nondefault == 0)   // default-linked onsets
replace onset_nd  = 0 if missing(nondefault)
replace onset_def = 0 if missing(nondefault)

label var onset_all  "Onset dummy — all 61 crisis episodes"
label var onset_nd   "Onset dummy — 40 non-default episodes"
label var onset_def  "Onset dummy — 21 default-linked episodes"

* ══════════════════════════════════════════════════════════════════════════
* 3. PANEL SETUP
* ══════════════════════════════════════════════════════════════════════════

xtset cid year
xtdescribe

* ══════════════════════════════════════════════════════════════════════════
* 4. GENERATE LP OUTCOME VARIABLES
*    dy_h  = log_gdppc(t+h) − log_gdppc(t−1)   for h = −2, −1, 0, 1, 2, 3, 4
*    Computed on the FULL panel (before dropping continuation years)
*    so that F. operator works correctly across gaps.
* ══════════════════════════════════════════════════════════════════════════

* Baseline: log GDP per capita at t-1
gen ln_gdppc_base = L.ln_gdppc
label var ln_gdppc_base "log real GDPpc at t-1 (LP baseline)"

* Forward differences (cumulative change relative to t-1)
forvalues h = 0/4 {
    gen dy_`h' = (F`h'.ln_gdppc - ln_gdppc_base) * 100
    label var dy_`h' "Cum. % change in log GDPpc: F`h' vs t-1"
}

* Pre-trend placebo horizons (must be ~0 under correct identification)
gen dy_m1 = (L0.ln_gdppc - ln_gdppc_base) * 100   // h = 0 relative to t-1 = 0 by construction; use L2 for h=-1
gen dy_m2 = (L2.ln_gdppc - ln_gdppc_base) * 100

* Re-define correctly:
drop dy_m1 dy_m2
gen dy_m1 = (ln_gdppc    - L2.ln_gdppc) * 100   // t vs t-2 → h=-1 pretrend
gen dy_m2 = (L.ln_gdppc  - L3.ln_gdppc) * 100   // t-1 vs t-3 → h=-2 pretrend
label var dy_m1 "Pre-trend h=-1: GDPpc(t) - GDPpc(t-2)"
label var dy_m2 "Pre-trend h=-2: GDPpc(t-1) - GDPpc(t-3)"

* ══════════════════════════════════════════════════════════════════════════
* 5. SAMPLE RESTRICTION FLAG
*    Estimation sample: excl_continuation == 0 (drops continuation years)
*    Applied via "if sample==1" in regressions — do NOT drop rows here
*    (dropping would break F. and L. operators)
* ══════════════════════════════════════════════════════════════════════════

gen sample = (continuation == 0) & !missing(ln_gdppc_base)
label var sample "=1 in estimation sample (onset + tranquil, not continuation)"

tab onset_all sample, miss
tab onset_nd  sample, miss

* ══════════════════════════════════════════════════════════════════════════
* 6. VARIABLE LABELS AND SAVE
* ══════════════════════════════════════════════════════════════════════════

label var cid        "Country numeric ID"
label var year       "Year"
label var ln_gdppc   "Log real GDP per capita"
label var gdpg       "Real GDP growth (%)"
label var l1_gdpg    "L1 real GDP growth"
label var l2_gdpg    "L2 real GDP growth"
label var ca         "Current account / GDP"
label var debt       "Public debt / GDP"
label var infl       "CPI inflation"
label var spr_mean   "EMBIG mean spread (bps)"
label var l_spr_mean "L1 EMBIG mean spread (bps)"
label var imf        "IMF program dummy"
label var vix        "VIX (annual average)"
label var ust10y     "US 10yr Treasury yield"
label var nondefault "1 = non-default episode onset"

compress
save "$clean/panel_lp.dta", replace

di as result "panel_lp.dta saved. Obs: `=_N', Countries: `=r(max)' (from xtdescribe above)"

* ══════════════════════════════════════════════════════════════════════════
* 7. CRISIS TREATMENT DIAGNOSTIC
*    Identify ever-treated vs. never-treated countries
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== CRISIS TREATMENT DIAGNOSTIC ==="

bysort country: egen ever_crisis = max(onset_all)
bysort country: egen n_onsets   = total(onset_all)

* Ever-treated: number of episodes per country
di as result _n "--- Ever-treated countries (n_onsets >= 1) ---"
preserve
    bysort country: keep if _n == 1
    keep if ever_crisis == 1
    sort n_onsets country
    list country n_onsets, sep(0) noobs abbrev(30)
    quietly count
    di as result "Total ever-treated countries: " r(N)
restore

* Never-treated: pure control group
di as result _n "--- Never-treated countries (pure control group) ---"
preserve
    bysort country: keep if _n == 1
    keep if ever_crisis == 0
    sort country
    list country, sep(0) noobs abbrev(30)
    quietly count
    di as result "Total never-treated countries: " r(N)
restore

drop ever_crisis n_onsets
