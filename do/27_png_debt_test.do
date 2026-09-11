/*===========================================================================
  27_PNG_DEBT_TEST.DO
  EXPLORATORY / ROBUSTNESS-ONLY TEST — NOT part of the headline six-variable
  design (GDP, investment, bank credit, claims on government, FDI, real
  lending rate) and NOT wired into 00_master.do, matching this project's own
  convention for standalone exploratory files (18b_mindur_variant.do,
  13e_nexus_bars.do, etc.).

  MOTIVATION: a candidate new transmission channel — private nonguaranteed
  (PNG) external debt scaled by GNI — as a proxy for "corporate external
  borrowing," the mechanism a related paper on corporate external borrowing
  in sovereign debt restructurings uses. This file tests whether that proxy
  behaves like a transmission channel in THIS project's own design, using
  the project's own established estimators (OLS resolution-split LP, AIPW),
  applied verbatim to one new variable. It does not touch $ctrl_core, any
  existing channel's construction, or any of 12_channels_resolution.do /
  13c_aipw_channels.do / 03c_table_combined_six.do's own active channel lists.

  DATA: data/raw/NetFlowToGNI.xlsx, sheet "Data" — World Bank IDS/WDI long
  (stacked-series) format, NOT this project's usual wide-by-indicator WDI
  layout. Two series stacked per country:
    DT.NFL.DPNG.CD  Net flows on external debt, PNG (NFL, current US$) — a
                    FLOW (can be negative: net repayment years)
    NY.GNP.MKTP.CD  GNI (current US$)
  Country Code is already ISO3, so this merges directly onto the panel's own
  `iso3' key — no name crosswalk needed.

  STOCK -> FLOW SWITCH (this file originally used data/raw/PNGtoGNI.xlsx's
  DT.DOD.DPNG.CD, a debt STOCK; superseded here). A stock's year-over-year
  change conflates new borrowing with valuation/exchange-rate effects on
  existing debt and debt that matures or is written off -- none of which is
  "the private sector substituting toward domestic financing," the actual
  mechanism this file tests. A FLOW (net new external borrowing that year)
  is the economically direct measure of that substitution margin, so the
  channel is rebuilt on DT.NFL.DPNG.CD instead.

  CHANNEL CONSTRUCTION CHOICE: plain ratio (nfl_gni = netflow/gni*100, ppt
  of GNI), NOT a log-real-level transform — GNI as denominator is kept for
  consistency with the file's own established naming/convention (the same
  denominator-procyclicality caveat that applies to any ratio channel in
  this project applies here too, not specific to GNI vs GDP).
  OUTCOME AT EACH HORIZON IS THE LEVEL, NOT A DIFFERENCE: a flow variable's
  own value already is the economically meaningful quantity (how much new
  external borrowing is happening that year) -- differencing a flow against
  a base-year level (as the retired stock version did) would remove it one
  step further from what is actually being measured. ch_pngdebt_h =
  nfl_gni at year t+h directly; pre_pngdebt = nfl_gni at year t-1 (the
  predetermined pre-crisis flow level, not a level-of-a-level difference).

  ══════════════════════════════════════════════════════════════════════════
  HOW TO RUN THIS FILE ONE PART AT A TIME
  ══════════════════════════════════════════════════════════════════════════
  This file is split into FOUR independent parts, each separated by a
  clearly marked "PART n" banner. Each part is fully self-contained: it
  starts with its own `use', builds anything it needs from scratch, and does
  not depend on any other part having just run in the same Stata session.
  PART 0 must be run at least ONCE first (it builds and saves
  $clean/panel_lp_png_test.dta, the checkpoint every later part reads from).
  After that, PARTS 1-3 can each be run independently, in any order, any
  number of times, by highlighting that part's block of code in the Stata
  do-file editor and running the selection (or running this whole file and
  just reading the console output section by section — both work).

    PART 0 — Import + merge + channel construction + coverage diagnostic.
             Run this first. Saves $clean/panel_lp_png_test.dta.
    PART 1 — Descriptive statistics (mirrors 03d_summary_statistics.do).
             Reads $clean/panel_lp_png_test.dta. Exports
             $tabs/png_debt_test_summary.xlsx.
    PART 2 — One-stage OLS (mirrors 12_channels_resolution.do's per-channel
             design). Reads $clean/panel_lp_png_test.dta. Exports
             $tabs/png_debt_test_ols.csv.
    PART 3 — AIPW (mirrors 13c_aipw_channels.do's per-channel design). Reads
             $clean/panel_lp_png_test.dta. Exports
             $tabs/png_debt_test_aipw.csv. Slowest part (1000-draw bootstrap
             x 5 horizons) — run this one last/separately if you just want
             the quick descriptive/OLS read first.
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* PART 0 — IMPORT + MERGE + CHANNEL CONSTRUCTION + COVERAGE DIAGNOSTIC
*   Run this first. Builds nfl_gni and its channel outcome, then SAVES
*   $clean/panel_lp_png_test.dta so Parts 1-3 can each start fresh from it
*   without re-running this import/merge step.
* ══════════════════════════════════════════════════════════════════════════

* -- Import + reshape data/raw/NetFlowToGNI.xlsx (long stacked-series -> wide) --
import excel "$raw/NetFlowToGNI.xlsx", sheet("Data") firstrow allstring clear

* Header cleanup: "Country Code" -> iso3 (already ISO3, no crosswalk needed);
* "1970 [YR1970]" style headers land as YR1970 ... YR2032 in Stata (matches
* 01d_merge_vulnerability.do's REER/IDS import convention exactly).
capture rename CountryCode iso3
capture rename CounterpartAreaName counterpart
capture rename SeriesCode series_code

keep iso3 counterpart series_code YR*
* Drop trailing metadata/footer rows (blank Country Code) and confirm the
* Counterpart-Area filter is the no-op it is documented to be above.
keep if length(iso3) == 3
quietly levelsof counterpart, local(cplevels) clean
di as result "  Counterpart-Area Name values present: `cplevels'"
keep if counterpart == "World"

foreach v of varlist YR* {
    destring `v', replace force   // ".." -> missing
}

tempfile png_raw
save `png_raw'

* -- Net flows on PNG external debt (DT.NFL.DPNG.CD) --
use `png_raw', clear
keep if series_code == "DT.NFL.DPNG.CD"
reshape long YR, i(iso3) j(year)
rename YR netflow
keep iso3 year netflow
tempfile t_png
save `t_png'

* -- GNI (NY.GNP.MKTP.CD) --
use `png_raw', clear
keep if series_code == "NY.GNP.MKTP.CD"
reshape long YR, i(iso3) j(year)
rename YR gni
keep iso3 year gni
tempfile t_gni
save `t_gni'

use `t_png', clear
merge 1:1 iso3 year using `t_gni', nogen
sort iso3 year
tempfile png_cy
save `png_cy'

* -- Merge onto panel_lp.dta + build nfl_gni (ratio, ppt of GNI) --
* CHANNEL CONSTRUCTION CHOICE: plain ratio, NOT a log-real-level transform.
* See header (STOCK -> FLOW SWITCH) for the full reasoning. netflow CAN be
* negative (net repayment years, a real and informative observation), so
* -- unlike the retired stock version -- there is no netflow>=0 filter here.
use "$clean/panel_lp.dta", clear
capture drop netflow gni pngdebt
merge m:1 iso3 year using `png_cy', keep(master match) nogen
sort cid year
xtset cid year

capture drop nfl_gni
gen double nfl_gni = netflow / gni * 100 if gni > 0 & !missing(netflow, gni)
label var nfl_gni "Net flows on PNG external debt / GNI, pct (flow-based external-financing exposure)"

* Channel outcome IS THE LEVEL at each horizon (not a difference from a
* base year, see header) -- ch_pngdebt_h = F h.nfl_gni; pre_pngdebt =
* L.nfl_gni (previous year's net-flow level, predetermined). Variable NAME
* ch_pngdebt_* kept unchanged from the retired stock version so Parts 1-3's
* downstream code needs no renaming -- only what feeds it changed.
forvalues h = 0/4 {
    capture drop ch_pngdebt_`h'
    gen double ch_pngdebt_`h' = F`h'.nfl_gni
}
capture drop pre_pngdebt
gen double pre_pngdebt = L.nfl_gni
label var pre_pngdebt "L1 level of nfl_gni (own pre-crisis flow level, predetermined)"

* -- Coverage diagnostic at onset (of 61 onsets), by resolution type --
di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 0 COMPLETE — COVERAGE: nfl_gni (ch_pngdebt_0) AT ONSET, BY RESOLUTION TYPE (of 61)"
di as result "════════════════════════════════════════════════════════════"
quietly count if onset_all == 1 & sample == 1 & !missing(ch_pngdebt_0)
local n_all = r(N)
quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_pngdebt_0)
local n_nd  = r(N)
quietly count if onset_def == 1 & sample == 1 & !missing(ch_pngdebt_0)
local n_def = r(N)
di as result "  ch_pngdebt_0: all=" `n_all' " / 61   non-default=" `n_nd' " / 39   default-linked=" `n_def' " / 22"

* -- Save the checkpoint Parts 1-3 will each read from --
save "$clean/panel_lp_png_test.dta", replace
di as result _n "Checkpoint saved: $clean/panel_lp_png_test.dta"
di as result "You can now run Part 1, Part 2, and/or Part 3 independently (any order, any number of times)."


* ══════════════════════════════════════════════════════════════════════════
* PART 1 — DESCRIPTIVE STATISTICS (mirrors 03d_summary_statistics.do's style)
*   Self-contained: reads $clean/panel_lp_png_test.dta directly. Requires
*   Part 0 to have been run at least once already.
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp_png_test.dta", clear
xtset cid year

di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 1 — DESCRIPTIVE STATISTICS: ch_pngdebt_0 (net flow / GNI level, Year 1, sample==1)"
di as result "════════════════════════════════════════════════════════════"
summarize ch_pngdebt_0 if sample==1

tempname S
tempfile sumf
postfile `S' str32 variable long obs double mean double sd double min double max using "`sumf'", replace
quietly summarize ch_pngdebt_0 if sample==1
post `S' ("Net flow on PNG debt / GNI (Year 1 level, ppt)") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max))
postclose `S'

preserve
    use "`sumf'", clear
    label var variable "Variable"
    label var obs      "Obs."
    label var mean      "Mean"
    label var sd        "Std. Dev."
    label var min        "Min"
    label var max        "Max"
    export excel "$tabs/png_debt_test_summary.xlsx", replace firstrow(varlabels)
restore
di as result "PART 1 COMPLETE — Descriptive statistics exported: $tabs/png_debt_test_summary.xlsx"
di as result "(native scale, no display rescaling — matches 03d_summary_statistics.do's current convention)"

* ── Pre/post-crisis evolution figure (style mirrors 02a_descriptive_facts.do) ──
* nfl_gni is a FLOW, so this plots its LEVEL at each year relative to onset
* (Year -3..-1 pre-crisis, Year 0..5 post-crisis) -- NOT a cumulative
* change from a base year, unlike the retired stock version and unlike
* 02a's own level-differenced channels. ch_pngdebt_m2/m3/m4 = the flow
* level at t-2/t-3/t-4 respectively (L2/L3/L4.nfl_gni); ch_pngdebt_0..4 =
* the flow level at t..t+4 (already built above). Still country-demeaned
* (dd_pngdebt_h = ch_pngdebt_h - country mean, over sample==1) so the
* figure reads a within-country deviation, matching every other
* descriptive figure's convention, even though the underlying quantity is
* now a level rather than a change.
* Year -1 is NOT zero by construction here (unlike the retired stock-diff
* version, where the base year was defined to be zero relative to itself)
* -- it is the flow level at t-1, computed like every other horizon.
capture drop ch_pngdebt_m1
gen double ch_pngdebt_m1 = L.nfl_gni
forvalues k = 2/4 {
    capture drop ch_pngdebt_m`k'
    gen double ch_pngdebt_m`k' = L`k'.nfl_gni
}

foreach h in m4 m3 m2 m1 0 1 2 3 4 {
    capture drop cmean_pngdebt_`h' dd_pngdebt_`h'
    quietly bysort cid: egen double cmean_pngdebt_`h' = mean(ch_pngdebt_`h') if sample==1
    quietly gen double dd_pngdebt_`h' = ch_pngdebt_`h' - cmean_pngdebt_`h' if sample==1
}

foreach g in all nd def {
    matrix desc_pngdebt_`g' = J(9, 1, .)
    quietly summarize dd_pngdebt_m1 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngdebt_`g'[4,1] = r(mean)
    quietly summarize dd_pngdebt_m4 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngdebt_`g'[1,1] = r(mean)
    quietly summarize dd_pngdebt_m3 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngdebt_`g'[2,1] = r(mean)
    quietly summarize dd_pngdebt_m2 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngdebt_`g'[3,1] = r(mean)
    forvalues h = 0/4 {
        quietly summarize dd_pngdebt_`h' if onset_`g'==1 & sample==1, meanonly
        matrix desc_pngdebt_`g'[`h'+5,1] = r(mean)
    }
    * Rebase so Year 0 (row 5, the crisis year) = 0 for this group's own
    * line -- easier to read as "change since crisis onset" than a level
    * deviation from each country's own long-run average. Each group (all/
    * nd/def) is rebased on ITS OWN Year-0 value, so the three lines stay
    * directly comparable to the un-rebased figure in level terms (only
    * the vertical offset changes, the shape/gap between lines does not).
    scalar _base_`g' = desc_pngdebt_`g'[5,1]
    forvalues r = 1/9 {
        matrix desc_pngdebt_`g'[`r',1] = desc_pngdebt_`g'[`r',1] - _base_`g'
    }
}

preserve
    clear
    svmat desc_pngdebt_all, names(b_all)
    svmat desc_pngdebt_nd,  names(b_nd)
    svmat desc_pngdebt_def, names(b_def)
    gen horizon = _n - 4     // rows are Year -3..5

    local c_nd  "blue"
    local c_def "red"
    twoway ///
        (line b_all1 horizon, lcolor(gs8) lwidth(medthick) lpattern(dash)) ///
        (connected b_nd1  horizon, lcolor("`c_nd'")  mcolor("`c_nd'")  msymbol(circle) lwidth(medthick)) ///
        (connected b_def1 horizon, lcolor("`c_def'") mcolor("`c_def'") msymbol(square) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xline(0, lpattern(dash) lcolor(gs10) lwidth(thin)) ///
        xlabel(-3(1)5, labsize(medsmall)) ///
        ylabel(, format(%9.1f) labsize(medsmall) angle(horizontal)) ///
        xtitle("Year (0 = crisis onset)", size(small)) ///
        ytitle("Country-demeaned net flow / GNI, pct", size(small)) ///
        title("Private external financing (net flows on PNG debt / GNI)", size(medium) color(navy)) ///
        legend(order(1 "All onsets" 2 "Non-default" 3 "Default-linked") ///
               position(6) rows(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig0_descriptive_pngdebt.pdf", replace
    di as result "Figure saved: fig0_descriptive_pngdebt.pdf (Years -3..5, country-demeaned mean, all/nd/def)"
restore


* ══════════════════════════════════════════════════════════════════════════
* PART 2 — ONE-STAGE OLS (mirrors 12_channels_resolution.do's per-channel
*   design EXACTLY: separate xtreg per arm, rival type dropped, country FE,
*   robust SE, $ctrl_core + pre_pngdebt; paired row-bootstrap difference +
*   Clogg z companion. NO common_abcd restriction — own best-available
*   sample.) Self-contained: reads $clean/panel_lp_png_test.dta directly.
*   Requires Part 0 to have been run at least once already.
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp_png_test.dta", clear
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
local ctrl_pngdebt $ctrl_core pre_pngdebt

set seed 20260819
local nboot = 1000     // matches 12_channels_resolution.do / 08b_aipw.do's own G=1000

* PROGRAM — copied verbatim from 12_channels_resolution.do's `_lpdiffboot'
* (see that file's header for the full argument: paired row-bootstrap of the
* def-nd difference, stratified so the tranquil pool and each arm's treated
* count are held fixed; Clogg et al. (1995) z as a companion statistic built
* from the two arms' own analytic SEs, permissive/independence-assuming,
* reported alongside the bootstrap, never in place of it).
capture program drop _lpdiffboot
program define _lpdiffboot, rclass
    syntax , Y(string) DND(string) IFND(string) DDEF(string) IFDEF(string) ///
             CTRLND(string) CTRLDEF(string) REPS(integer)

    capture xtreg `y' `dnd' `ctrlnd' if `ifnd', fe vce(robust)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bnd   = _b[`dnd']
    local send  = _se[`dnd']
    local nnd   = e(N)
    local ngnd  = e(N_g)

    capture xtreg `y' `ddef' `ctrldef' if `ifdef', fe vce(robust)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bdef  = _b[`ddef']
    local sedef = _se[`ddef']
    local ndef  = e(N)
    local ngdef = e(N_g)

    local dh = `bdef' - `bnd'

    capture drop _pool
    quietly gen byte _pool = 0 if (`ifnd') | (`ifdef')
    quietly replace _pool = 1 if `dnd' == 1
    quietly replace _pool = 2 if `ddef' == 1

    tempname pf
    tempfile bf
    quietly postfile `pf' double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture regress `y' `dnd' `ctrlnd' i.cid if `ifnd', vce(robust)
            local t1 = cond(_rc==0, _b[`dnd'], .)
            capture regress `y' `ddef' `ctrldef' i.cid if `ifdef', vce(robust)
            local t2 = cond(_rc==0, _b[`ddef'], .)
            if !missing(`t1') & !missing(`t2') quietly post `pf' (`t2' - `t1')
        restore
    }
    quietly postclose `pf'
    capture drop _pool

    local se = .
    local lo = .
    local hi = .
    local nd = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(diff)
        local nd = r(N)
        if `nd' >= 50 {
            quietly summarize diff
            local se = r(sd)
            _pctile diff, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
        }
    restore

    local pdiff = .
    if !missing(`se') & `se' > 0 local pdiff = 2*(1 - normal(abs(`dh'/`se')))

    local cloggz = .
    local cloggp = .
    if !missing(`send') & !missing(`sedef') & (`send'^2 + `sedef'^2) > 0 {
        local cloggz = `dh' / sqrt(`send'^2 + `sedef'^2)
        local cloggp = 2*(1 - normal(abs(`cloggz')))
    }

    return scalar ok    = 1
    return scalar bnd   = `bnd'
    return scalar send  = `send'
    return scalar nnd   = `nnd'
    return scalar ngnd  = `ngnd'
    return scalar bdef  = `bdef'
    return scalar sedef = `sedef'
    return scalar ndef  = `ndef'
    return scalar ngdef = `ngdef'
    return scalar dh    = `dh'
    return scalar se    = `se'
    return scalar lo    = `lo'
    return scalar hi    = `hi'
    return scalar nboot = `nd'
    return scalar p     = `pdiff'
    return scalar cloggz = `cloggz'
    return scalar cloggp = `cloggp'
end

di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 2 — ONE-STAGE OLS: Net flow on PNG debt / GNI channel, non-default vs default-linked"
di as result "════════════════════════════════════════════════════════════"
di "h   b_nd     b_def    p(nd=def)   Clogg z (p)"

tempname O
tempfile olsf
postfile `O' byte horizon double b_nd se_nd double b_def se_def ///
    double diff se_diff lo95 hi95 pdiff double cloggz cloggp ///
    long n_nd n_def using "`olsf'", replace

forvalues h = 0/4 {
    _lpdiffboot, y(ch_pngdebt_`h') ///
        dnd(onset_nd)  ifnd(sample==1 & onset_def==0) ///
        ddef(onset_def) ifdef(sample==1 & onset_nd==0) ///
        ctrlnd(`ctrl_pngdebt') ctrldef(`ctrl_pngdebt') reps(`nboot')

    if r(ok) {
        local bnd  = r(bnd)
        local bdef = r(bdef)
        local snd  = r(send)
        local sdef = r(sedef)
        local dh_  = r(dh)
        local se_  = r(se)
        local lo_  = r(lo)
        local hi_  = r(hi)
        local p_   = r(p)
        local cz_  = r(cloggz)
        local cp_  = r(cloggp)
        local nnd_ = r(nnd)
        local ndef_= r(ndef)

        post `O' (`h') (`bnd') (`snd') (`bdef') (`sdef') ///
            (`dh_') (`se_') (`lo_') (`hi_') (`p_') (`cz_') (`cp_') (`nnd_') (`ndef_')

        di "h=" `h'+1 "  " %7.3f `bnd' "  " %7.3f `bdef' "  " %5.3f `p_' ///
           "  Clogg z=" %6.3f `cz_' " (p=" %5.3f `cp_' ")"
    }
    else {
        post `O' (`h') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
        di as error "regression failed for nfl_gni h=" `h'+1
    }
}
postclose `O'

preserve
    use "`olsf'", clear
    label var horizon "Horizon h (0..4)"
    label var b_nd    "Non-default onset coefficient"
    label var se_nd   "Robust SE, non-default"
    label var b_def   "Default-linked onset coefficient"
    label var se_def  "Robust SE, default-linked"
    label var diff    "Difference (default - non-default)"
    label var se_diff "Paired row-bootstrap SE"
    label var lo95    "95% bootstrap CI, lower"
    label var hi95    "95% bootstrap CI, upper"
    label var pdiff   "p-value (paired row bootstrap)"
    label var cloggz  "Clogg et al. (1995) z (companion, analytic SEs)"
    label var cloggp  "p-value of Clogg z"
    label var n_nd    "Observations, non-default regression"
    label var n_def   "Observations, default-linked regression"
    export delimited "$tabs/png_debt_test_ols.csv", replace
restore
di as result "PART 2 COMPLETE — OLS results exported: $tabs/png_debt_test_ols.csv"


* ══════════════════════════════════════════════════════════════════════════
* PART 3 — AIPW (mirrors 13c_aipw_channels.do's per-channel design EXACTLY:
*   Act 2 resolution split via _aipwpair, row-bootstrap level SEs (ADOPTED),
*   Clogg z companion, cz_def propensity predictors. NO common_abcd
*   restriction.) Self-contained: reads $clean/panel_lp_png_test.dta
*   directly. Requires Part 0 to have been run at least once already.
*   SLOWEST PART — 1000-draw bootstrap x 5 horizons x 2 arms.
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp_png_test.dta", clear
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
local cz_def l_fedfunds l_contagion_dist_atdef years_since_def_onset

* Outcome-model controls: core_aipw + pre_pngdebt. Judgment call on whether
* any $ctrl_core term should be dropped the way l_credit_bank is dropped for
* the credit channel (that channel's own outcome IS effectively the lagged
* dependent variable, correlation 0.950): l_credit_bank is BANK credit to the
* private sector, a domestic-banking-system claim; nfl_gni is PRIVATE
* NONGUARANTEED EXTERNAL debt, i.e. corporate borrowing from foreign
* creditors, a conceptually distinct balance-sheet object (external vs.
* domestic counterparty) with no accounting identity linking the two, unlike
* the credit/credit_bank case (same "private credit" concept measured two
* ways). No term in $ctrl_core is the channel's own lagged level or a close
* proxy for it, so none is dropped -- the full core_aipw + pre_pngdebt set is
* used, matching claims_govt/inv/fdi/real_lending's own treatment. This is a
* judgment call, not a tested correlation -- worth checking directly
* (correlate l_credit_bank l_debt nfl_gni) once this part is actually run.
local core_aipw l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2
local om_pngdebt `core_aipw' pre_pngdebt

set seed 20260819

* PROGRAMS — copied verbatim from 13c_aipw_channels.do (_mkstrat, _aipw,
* _aipwpair). See that file's header for the full estimator rationale.
capture program drop _mkstrat
program define _mkstrat
    syntax varlist(min=1 max=2) [if], GENerate(name)
    marksample touse
    capture drop `generate'
    tempvar t1 t2
    local d1 : word 1 of `varlist'
    local d2 : word 2 of `varlist'
    quietly gen byte `t1' = `d1' if `touse'
    bysort cid: egen byte `generate' = max(`t1')
    quietly replace `generate' = 0 if missing(`generate')
    if "`d2'" != "" {
        quietly gen byte `t2' = `d2' if `touse'
        tempvar s2
        bysort cid: egen byte `s2' = max(`t2')
        quietly replace `s2' = 0 if missing(`s2')
        quietly replace `generate' = `generate' + 2*`s2'
    }
end

capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'
    tempvar xb m0 m1 ps summ iwt
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' [pweight=`iwt'] if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' [pweight=`iwt'] if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    local th = r(mean)
    local nn = r(N)

    tempvar _tagcty
    quietly egen byte `_tagcty' = tag(cid) if `touse'
    quietly count if `_tagcty'==1
    local nctry = r(N)

    quietly count if `touse' & `D'==1
    local ntr = r(N)

    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta  = `th'
    return scalar N      = `nn'
    return scalar nctry  = `nctry'
    return scalar ntreat = `ntr'
    return scalar se     = `sean'
end

capture program drop _aipwpair
program define _aipwpair, rclass
    syntax , Y(string) D1(string) IF1(string) D2(string) IF2(string) ///
             OMOD(string) PZ(string) REPS(integer)
    capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b1 = r(theta)
    local a1 = r(se)
    local n1 = r(N)
    local nctry1 = r(nctry)
    local ntreat1 = r(ntreat)
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2 = r(theta)
    local a2 = r(se)
    local n2 = r(N)
    local nctry2 = r(nctry)
    local ntreat2 = r(ntreat)
    local dh = `b1' - `b2'

    capture drop _pool
    quietly gen byte _pool = 0 if (`if1') | (`if2')
    quietly replace _pool = 1 if `d1' == 1
    quietly replace _pool = 2 if `d2' == 1

    tempname pf
    tempfile bf
    quietly postfile `pf' double t1 double t2 double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
            local t1 = cond(_rc==0, r(theta), .)
            capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
            local t2 = cond(_rc==0, r(theta), .)
            if !missing(`t1') & !missing(`t2') quietly post `pf' (`t1') (`t2') (`t1' - `t2')
        restore
    }
    quietly postclose `pf'
    capture drop _pool
    local se = .
    local lo = .
    local hi = .
    local nd = 0
    local bse1 = .
    local bse2 = .
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(diff)
        local nd = r(N)
        if `nd' >= 50 {
            quietly summarize diff
            local se = r(sd)
            _pctile diff, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
            quietly summarize t1
            local bse1 = r(sd)
            quietly summarize t2
            local bse2 = r(sd)
        }
    restore
    return scalar ok = 1
    return scalar dh = `dh'
    return scalar b1 = `b1'
    return scalar b2 = `b2'
    return scalar a1 = `a1'
    return scalar a2 = `a2'
    return scalar bse1 = `bse1'
    return scalar bse2 = `bse2'
    return scalar n1 = `n1'
    return scalar n2 = `n2'
    return scalar nctry1 = `nctry1'
    return scalar nctry2 = `nctry2'
    return scalar ntreat1 = `ntreat1'
    return scalar ntreat2 = `ntreat2'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

local nboot_aipw = 1000     // matches 13c_aipw_channels.do's own G=1000

di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 3 — AIPW (Act 2): Net flow on PNG debt / GNI channel, non-default vs default-linked"
di as result "════════════════════════════════════════════════════════════"
di as result "  h   ND (se_boot)     DEF (se_boot)     def-nd   [95% boot CI]   Clogg z    p"
di as result "  se_boot = ROW-BOOTSTRAP SE (ADOPTED, matching 13c_aipw_channels.do's departure from"
di as result "  the analytic-SE formula -- see that file's header for the diagnostic that motivated it)."

tempname A
tempfile aipwf
postfile `A' byte horizon double b_nd se_nd double b_def se_def ///
    double diff se_diff lo95 hi95 double cloggz cloggp ///
    long n_nd n_def nctry_nd nctry_def ntreat_nd ntreat_def using "`aipwf'", replace

forvalues h = 0/4 {
    _aipwpair, y(ch_pngdebt_`h') ///
        d1(onset_def) if1(sample==1 & onset_nd==0) ///
        d2(onset_nd)  if2(sample==1 & onset_def==0) ///
        omod(`om_pngdebt') pz(`om_pngdebt' `cz_def') reps(`nboot_aipw')

    if r(ok) {
        local B1 = r(b1)     // default-linked ATE
        local B2 = r(b2)     // non-default ATE
        local A1 = r(a1)
        local A2 = r(a2)
        local BSE1 = r(bse1)
        local BSE2 = r(bse2)
        local DH = r(dh)
        local SE = r(se)
        local LO = r(lo)
        local HI = r(hi)
        local ND = r(nd)
        local N1 = r(n1)
        local N2 = r(n2)
        local NC1 = r(nctry1)
        local NC2 = r(nctry2)
        local NT1 = r(ntreat1)
        local NT2 = r(ntreat2)

        local zz = .
        local pz = .
        if !missing(`A1') & !missing(`A2') & (`A1'^2 + `A2'^2) > 0 {
            local zz = `DH' / sqrt(`A1'^2 + `A2'^2)
            local pz = 2*(1 - normal(abs(`zz')))
        }

        post `A' (`h') (`B2') (`BSE2') (`B1') (`BSE1') ///
            (`DH') (`SE') (`LO') (`HI') (`zz') (`pz') ///
            (`N2') (`N1') (`NC2') (`NC1') (`NT2') (`NT1')

        local tnd  = cond(`BSE2'>0, `B2'/`BSE2', .)
        local pnd  = cond(!missing(`tnd'), 2*(1-normal(abs(`tnd'))), .)
        local sgnd = cond(missing(`pnd'), "", cond(`pnd'<.01,"***",cond(`pnd'<.05,"**",cond(`pnd'<.10,"*",""))))
        local tdef  = cond(`BSE1'>0, `B1'/`BSE1', .)
        local pdef  = cond(!missing(`tdef'), 2*(1-normal(abs(`tdef'))), .)
        local sgdef = cond(missing(`pdef'), "", cond(`pdef'<.01,"***",cond(`pdef'<.05,"**",cond(`pdef'<.10,"*",""))))
        local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")

        di "    " %1.0f `h'+1 "  " %8.3f `B2' "`sgnd'" " (" %5.3f `BSE2' ")  " ///
           %8.3f `B1' "`sgdef'" " (" %5.3f `BSE1' ")  " %8.3f `DH' ///
           " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
           " " %7.3f `zz' " " %5.3f `pz'
    }
    else {
        post `A' (`h') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.)
        di as error "    h=" `h'+1 ": Act 2 estimate failed (thin sample)."
    }
}
postclose `A'

preserve
    use "`aipwf'", clear
    label var horizon "Horizon h (0..4)"
    label var b_nd    "Non-default ATE"
    label var se_nd   "Row-bootstrap SE, non-default (ADOPTED)"
    label var b_def   "Default-linked ATE"
    label var se_def  "Row-bootstrap SE, default-linked (ADOPTED)"
    label var diff    "def - nd difference"
    label var se_diff "Row-bootstrap SE of the difference"
    label var lo95    "95% bootstrap CI lower"
    label var hi95    "95% bootstrap CI upper"
    label var cloggz  "Clogg et al. (1995) z (companion, analytic SEs)"
    label var cloggp  "p-value of Clogg z"
    label var n_nd    "Observations, non-default outcome-regression sample"
    label var n_def   "Observations, default-linked outcome-regression sample"
    label var nctry_nd  "Countries, non-default sample"
    label var nctry_def "Countries, default-linked sample"
    label var ntreat_nd  "Treated onsets, non-default sample"
    label var ntreat_def "Treated onsets, default-linked sample"
    export delimited "$tabs/png_debt_test_aipw.csv", replace
restore
di as result "PART 3 COMPLETE — AIPW results exported: $tabs/png_debt_test_aipw.csv"

di as result _n "27_png_debt_test.do complete. EXPLORATORY/ROBUSTNESS ONLY -- not wired into"
di as result "00_master.do, $ctrl_core, or any headline channel list."
