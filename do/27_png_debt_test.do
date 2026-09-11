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

  DATA: data/raw/PNGtoGNI.xlsx, sheet "Data" — World Bank IDS/WDI long
  (stacked-series) format, NOT this project's usual wide-by-indicator WDI
  layout. Two series stacked per country:
    DT.DOD.DPNG.CD  External debt stocks, PNG (DOD, current US$) — a STOCK
    NY.GNP.MKTP.CD  GNI (current US$)
  Country Code is already ISO3, so this merges directly onto the panel's own
  `iso3' key — no name crosswalk needed.

  CHANNEL CONSTRUCTION CHOICE: plain ratio (png_gni = pngdebt/gni*100, ppt
  of GNI), NOT a log-real-level transform — see PART 0 below for the full
  reasoning (coverage/zero-value tradeoff, same logic 18_transforms.do
  already applies to claims_govt).

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
*   Run this first. Builds png_gni and its channel outcome, then SAVES
*   $clean/panel_lp_png_test.dta so Parts 1-3 can each start fresh from it
*   without re-running this import/merge step.
* ══════════════════════════════════════════════════════════════════════════

* -- Import + reshape data/raw/PNGtoGNI.xlsx (long stacked-series -> wide) --
import excel "$raw/PNGtoGNI.xlsx", sheet("Data") firstrow allstring clear

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

* -- PNG debt stock (DT.DOD.DPNG.CD) --
use `png_raw', clear
keep if series_code == "DT.DOD.DPNG.CD"
reshape long YR, i(iso3) j(year)
rename YR pngdebt
keep iso3 year pngdebt
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

* -- Merge onto panel_lp.dta + build png_gni (ratio, ppt of GNI) --
* CHANNEL CONSTRUCTION CHOICE: plain ratio, NOT a log-real-level transform.
* 18_transforms.do's own reasoning for keeping claims_govt on a plain ratio
* applies here on similar grounds: PNG debt is IDS-reported and covers only
* debtor-reporting low/middle-income economies, with several of this panel's
* higher-income default cases expected to be missing outright. Multiplying
* by gdp_real and logging would only ever LOSE further observations
* relative to the ratio (any zero PNG-debt year -- a real, informative "no
* corporate external borrowing" observation -- drops out under a log), on a
* variable whose coverage is already the binding constraint. The ratio-to-
* GNI form keeps every non-missing observation, including true zeros, and
* is directly interpretable (ppt of GNI) the same way claims_govt is.
use "$clean/panel_lp.dta", clear
capture drop pngdebt gni
merge m:1 iso3 year using `png_cy', keep(master match) nogen
sort cid year
xtset cid year

capture drop png_gni
gen double png_gni = pngdebt / gni * 100 if pngdebt >= 0 & gni > 0 & !missing(pngdebt, gni)
label var png_gni "Private nonguaranteed external debt / GNI, pct (test channel; plain ratio, see header)"

* Channel outcome: ch_pngdebt_h = F h.png_gni - L.png_gni, h=0..4;
* pre_pngdebt = L.png_gni - L2.png_gni (own pre-trend control, matching
* 12_channels_resolution.do / 13c_aipw_channels.do's identical idiom).
capture drop pngdebt_base
gen double pngdebt_base = L.png_gni
forvalues h = 0/4 {
    capture drop ch_pngdebt_`h'
    gen double ch_pngdebt_`h' = F`h'.png_gni - pngdebt_base
}
capture drop pre_pngdebt
gen double pre_pngdebt = L.png_gni - L2.png_gni
label var pre_pngdebt "L1-L2 change in png_gni (own pre-crisis trend, predetermined)"

* -- Coverage diagnostic at onset (of 61 onsets), by resolution type --
di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 0 COMPLETE — COVERAGE: png_gni (ch_pngdebt_0) AT ONSET, BY RESOLUTION TYPE (of 61)"
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
di as result "PART 1 — DESCRIPTIVE STATISTICS: ch_pngdebt_0 (Year 1 cumulative change, sample==1)"
di as result "════════════════════════════════════════════════════════════"
summarize ch_pngdebt_0 if sample==1

tempname S
tempfile sumf
postfile `S' str32 variable long obs double mean double sd double min double max using "`sumf'", replace
quietly summarize ch_pngdebt_0 if sample==1
post `S' ("PNG debt / GNI (Year 1 chg, ppt)") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max))
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
di as result "PART 2 — ONE-STAGE OLS: PNG debt / GNI channel, non-default vs default-linked"
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
        di as error "regression failed for png_gni h=" `h'+1
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
* private sector, a domestic-banking-system claim; png_gni is PRIVATE
* NONGUARANTEED EXTERNAL debt, i.e. corporate borrowing from foreign
* creditors, a conceptually distinct balance-sheet object (external vs.
* domestic counterparty) with no accounting identity linking the two, unlike
* the credit/credit_bank case (same "private credit" concept measured two
* ways). No term in $ctrl_core is the channel's own lagged level or a close
* proxy for it, so none is dropped -- the full core_aipw + pre_pngdebt set is
* used, matching claims_govt/inv/fdi/real_lending's own treatment. This is a
* judgment call, not a tested correlation -- worth checking directly
* (correlate l_credit_bank l_debt png_gni) once this part is actually run.
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
di as result "PART 3 — AIPW (Act 2): PNG debt / GNI channel, non-default vs default-linked"
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
