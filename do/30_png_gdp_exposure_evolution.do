/*===========================================================================
  30_PNG_GDP_EXPOSURE_EVOLUTION.DO
  EXPLORATORY / ROBUSTNESS-ONLY TEST -- NOT part of the headline six-variable
  design and NOT wired into 00_master.do, matching this project's own
  convention for standalone exploratory files (27_png_debt_test.do,
  29_png_stock_gdp_exposure.do, etc.).

  PURPOSE. The full 27_png_debt_test.do-style battery (descriptive
  evolution figure + one-stage OLS + AIPW), but for the exposure indicator
  built in 29_png_stock_gdp_exposure.do -- png_gdp = PNG debt stock / GDP
  (current US$ both sides, WEO NGDPD denominator -- see 29's header for
  the full denominator-fix history) -- instead of the retired flow-based
  measure. Tests whether the exposure STOCK ITSELF moves differently by
  resolution type (nd vs def), exactly mirroring 27's own question but for
  this stock/GDP indicator. This is DISTINCT from
  29_png_stock_gdp_exposure.do, which tests whether pre-crisis exposure
  LEVEL predicts heterogeneous GDP/credit/investment OUTCOMES (a
  high/low-split amplifier design) -- this file does not touch that
  question at all.

  OUTCOME SCALE: png_gdp's own LEVEL at each horizon (not a difference from
  a base year) -- a stock ratio like debt/GDP is naturally read as a level
  path, matching 27's own final choice for its flow-based measure and this
  file's own Part 1 evolution figure.

  ══════════════════════════════════════════════════════════════════════════
  HOW TO RUN THIS FILE ONE PART AT A TIME
  ══════════════════════════════════════════════════════════════════════════
  Same convention as 27_png_debt_test.do: FOUR independent parts, each
  self-contained. PART 0 must be run at least once first (builds and saves
  $clean/panel_lp_pngexp_test.dta, the checkpoint every later part reads
  from). After that, PARTS 1-3 can each be run independently, any order,
  any number of times, by highlighting that part's block and running the
  selection.

    PART 0 -- Import + merge + construction + coverage diagnostic.
              Run this first. Saves $clean/panel_lp_pngexp_test.dta.
    PART 1 -- Descriptive evolution figure (mirrors 02a's Year-0 convention,
              same as 27's own Part 1). Exports
              $figs/fig0_evolution_png_gdp_exposure.pdf.
    PART 2 -- One-stage OLS (mirrors 27's Part 2 / 12_channels_resolution.do's
              design). Exports $tabs/png_gdp_exposure_ols.csv.
    PART 3 -- AIPW (mirrors 27's Part 3 / 13c_aipw_channels.do's design).
              Exports $tabs/png_gdp_exposure_aipw.csv. Slowest part.
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* PART 0 -- IMPORT + MERGE + CONSTRUCTION + COVERAGE DIAGNOSTIC
*   Run this first. Builds png_gdp and its outcome series, then SAVES
*   $clean/panel_lp_pngexp_test.dta so Parts 1-3 can each start fresh.
* ══════════════════════════════════════════════════════════════════════════
import excel "$raw/PNGtoGNI.xlsx", sheet("Data") firstrow allstring clear

capture rename CountryCode iso3
capture rename CounterpartAreaName counterpart
capture rename SeriesCode series_code

keep iso3 counterpart series_code YR*
keep if length(iso3) == 3
keep if counterpart == "World"

foreach v of varlist YR* {
    destring `v', replace force
}

tempfile png_raw
save `png_raw'

use `png_raw', clear
keep if series_code == "DT.DOD.DPNG.CD"
reshape long YR, i(iso3) j(year)
rename YR pngdebt_stock
keep iso3 year pngdebt_stock
tempfile png_cy
save `png_cy'

* -- GDP in current US$ (WEO NGDPD, billions) -- same mechanics as 11_weo.do
*    / 29_png_stock_gdp_exposure.do (Excel letters AB..CA = years 1980..2031).
import excel "$raw/WEOApr2026all.xlsx", sheet("Countries") firstrow clear
capture rename COUNTRYID   iso3
capture rename INDICATORID indid
local col AB AC AD AE AF AG AH AI AJ AK AL AM AN AO AP AQ AR AS AT AU ///
          AV AW AX AY AZ BA BB BC BD BE BF BG BH BI BJ BK BL BM BN BO ///
          BP BQ BR BS BT BU BV BW BX BY BZ CA
local y = 1980
foreach c of local col {
    capture rename `c' yr`y'
    local ++y
}
keep iso3 indid yr*
keep if length(iso3) == 3
keep if indid == "NGDPD"
capture destring yr*, replace force
reshape long yr, i(iso3 indid) j(year)
drop if missing(yr)
rename yr gdp_usd_bn
keep iso3 year gdp_usd_bn
label var gdp_usd_bn "GDP, current prices, US$ billions (IMF WEO NGDPD)"
tempfile weo_gdp
save `weo_gdp'

use "$clean/panel_lp.dta", clear
capture drop pngdebt_stock gdp_usd_bn
merge m:1 iso3 year using `png_cy', keep(master match) nogen
merge m:1 iso3 year using `weo_gdp', keep(master match) nogen
sort cid year
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

* gdp_usd_bn is in US$ BILLIONS; pngdebt_stock is in raw current US$ --
* convert gdp_usd_bn to raw US$ (x 1e9) before dividing.
capture drop png_gdp
gen double png_gdp = pngdebt_stock / (gdp_usd_bn * 1e9) * 100 ///
    if pngdebt_stock >= 0 & gdp_usd_bn > 0 & !missing(pngdebt_stock, gdp_usd_bn)
label var png_gdp "Private nonguaranteed external debt stock / GDP, pct"

* Outcome IS THE LEVEL at each horizon (not a difference from a base year,
* see header) -- ch_pngexp_h = F h.png_gdp; pre_pngexp = own pre-crisis
* trend control (L.png_gdp - L2.png_gdp), matching every other channel's
* pre_<v> convention in this project even though the outcome itself is a
* level here.
capture drop ch_pngexp_0
gen double ch_pngexp_0 = png_gdp
forvalues h = 1/4 {
    capture drop ch_pngexp_`h'
    gen double ch_pngexp_`h' = F`h'.png_gdp
}
capture drop pre_pngexp
gen double pre_pngexp = L.png_gdp - L2.png_gdp
label var pre_pngexp "L1-L2 change in png_gdp (own pre-crisis trend, predetermined)"

* -- Coverage diagnostic at onset (of 61 onsets), by resolution type --
di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 0 COMPLETE -- COVERAGE: png_gdp (ch_pngexp_0) AT ONSET, BY RESOLUTION TYPE (of 61)"
di as result "════════════════════════════════════════════════════════════"
quietly count if onset_all == 1 & sample == 1 & !missing(ch_pngexp_0)
local n_all = r(N)
quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_pngexp_0)
local n_nd  = r(N)
quietly count if onset_def == 1 & sample == 1 & !missing(ch_pngexp_0)
local n_def = r(N)
di as result "  ch_pngexp_0: all=" `n_all' " / 61   non-default=" `n_nd' " / 39   default-linked=" `n_def' " / 22"

quietly summarize png_gdp if sample==1
di as result _n "  png_gdp range check (sample==1): N=" r(N) "  mean=" %6.2f r(mean) ///
    "  min=" %6.2f r(min) "  max=" %6.2f r(max)

* -- Save the checkpoint Parts 1-3 will each read from --
save "$clean/panel_lp_pngexp_test.dta", replace
di as result _n "Checkpoint saved: $clean/panel_lp_pngexp_test.dta"
di as result "You can now run Part 1, Part 2, and/or Part 3 independently (any order, any number of times)."


* ══════════════════════════════════════════════════════════════════════════
* PART 1 -- DESCRIPTIVE EVOLUTION FIGURE (mirrors 02a_descriptive_facts.do's
*   Year-0 convention EXACTLY, same construction as 27_png_debt_test.do's
*   own Part 1). Self-contained: reads the checkpoint directly. Requires
*   Part 0 to have been run at least once already.
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp_pngexp_test.dta", clear
xtset cid year

capture drop pngexp_base
gen double pngexp_base = L.png_gdp

capture drop ch_pngexp_m1
gen double ch_pngexp_m1 = L.png_gdp
forvalues k = 2/4 {
    capture drop ch_pngexp_m`k'
    gen double ch_pngexp_m`k' = L`k'.png_gdp
}

foreach h in m4 m3 m2 m1 0 1 2 3 4 {
    capture drop cmean_pngexp_`h' dd_pngexp_`h'
    quietly bysort cid: egen double cmean_pngexp_`h' = mean(ch_pngexp_`h') if sample==1
    quietly gen double dd_pngexp_`h' = ch_pngexp_`h' - cmean_pngexp_`h' if sample==1
}

foreach g in all nd def {
    matrix desc_pngexp_`g' = J(9, 1, .)
    quietly summarize dd_pngexp_m1 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[4,1] = r(mean)
    quietly summarize dd_pngexp_m4 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[1,1] = r(mean)
    quietly summarize dd_pngexp_m3 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[2,1] = r(mean)
    quietly summarize dd_pngexp_m2 if onset_`g'==1 & sample==1, meanonly
    matrix desc_pngexp_`g'[3,1] = r(mean)
    forvalues h = 0/4 {
        quietly summarize dd_pngexp_`h' if onset_`g'==1 & sample==1, meanonly
        matrix desc_pngexp_`g'[`h'+5,1] = r(mean)
    }
    * Rebase so Year 0 (row 4, the LAST PRE-CRISIS year, t-1) = 0.
    scalar _base_`g' = desc_pngexp_`g'[4,1]
    forvalues r = 1/9 {
        matrix desc_pngexp_`g'[`r',1] = desc_pngexp_`g'[`r',1] - _base_`g'
    }
}

preserve
    clear
    svmat desc_pngexp_all, names(b_all)
    svmat desc_pngexp_nd,  names(b_nd)
    svmat desc_pngexp_def, names(b_def)
    gen horizon = _n - 4     // rows are Year -3..5

    local c_nd  "blue"
    local c_def "red"
    twoway ///
        (line b_all1 horizon, lcolor(gs8) lwidth(medthick) lpattern(dash)) ///
        (connected b_nd1  horizon, lcolor("`c_nd'")  mcolor("`c_nd'")  msymbol(circle) lwidth(medthick)) ///
        (connected b_def1 horizon, lcolor("`c_def'") mcolor("`c_def'") msymbol(square) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xline(1, lpattern(dash) lcolor(gs12) lwidth(thin)) ///
        xlabel(-3(1)5, labsize(medsmall)) ///
        ylabel(, format(%9.1f) labsize(medsmall) angle(horizontal)) ///
        xtitle("Year (0 = last pre-crisis year, 1 = crisis onset)", size(small)) ///
        ytitle("PNG debt / GDP, pct, relative to Year 0", size(small)) ///
        title("Private external-debt exposure (PNG debt stock / GDP)", size(medium) color(navy)) ///
        legend(order(1 "All onsets" 2 "Non-default" 3 "Default-linked") ///
               position(6) rows(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig0_evolution_png_gdp_exposure.pdf", replace
    di as result "Figure saved: fig0_evolution_png_gdp_exposure.pdf (Years -3..5, rebased to Year 0 = last pre-crisis year; all/nd/def)"
restore
di as result "PART 1 COMPLETE."


* ══════════════════════════════════════════════════════════════════════════
* PART 2 -- ONE-STAGE OLS (mirrors 27_png_debt_test.do's Part 2 /
*   12_channels_resolution.do's per-channel design EXACTLY: separate xtreg
*   per arm, rival type dropped, country FE, robust SE, $ctrl_core +
*   pre_pngexp; paired row-bootstrap difference + Clogg z companion.
*   Self-contained: reads $clean/panel_lp_pngexp_test.dta directly.
*   Requires Part 0 to have been run at least once already.
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp_pngexp_test.dta", clear
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
local ctrl_pngexp $ctrl_core pre_pngexp

set seed 20260819
local nboot = 1000

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
di as result "PART 2 -- ONE-STAGE OLS: PNG debt/GDP exposure, non-default vs default-linked"
di as result "════════════════════════════════════════════════════════════"
di "h   b_nd     b_def    p(nd=def)   Clogg z (p)"

tempname O
tempfile olsf
postfile `O' byte horizon double b_nd se_nd double b_def se_def ///
    double diff se_diff lo95 hi95 pdiff double cloggz cloggp ///
    long n_nd n_def using "`olsf'", replace

forvalues h = 0/4 {
    _lpdiffboot, y(ch_pngexp_`h') ///
        dnd(onset_nd)  ifnd(sample==1 & onset_def==0) ///
        ddef(onset_def) ifdef(sample==1 & onset_nd==0) ///
        ctrlnd(`ctrl_pngexp') ctrldef(`ctrl_pngexp') reps(`nboot')

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
        di as error "regression failed for png_gdp h=" `h'+1
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
    export delimited "$tabs/png_gdp_exposure_ols.csv", replace
restore
di as result "PART 2 COMPLETE -- OLS results exported: $tabs/png_gdp_exposure_ols.csv"


* ══════════════════════════════════════════════════════════════════════════
* PART 3 -- AIPW (mirrors 27_png_debt_test.do's Part 3 / 13c_aipw_channels.do's
*   per-channel design EXACTLY: Act 2 resolution split via _aipwpair,
*   row-bootstrap level SEs, Clogg z companion, cz_def propensity
*   predictors. Self-contained: reads $clean/panel_lp_pngexp_test.dta
*   directly. Requires Part 0 to have been run at least once already.
*   SLOWEST PART -- 1000-draw bootstrap x 5 horizons x 2 arms.
* ══════════════════════════════════════════════════════════════════════════
use "$clean/panel_lp_pngexp_test.dta", clear
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
local cz_def l_fedfunds l_contagion_dist_atdef years_since_def_onset

local core_aipw l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2
local om_pngexp `core_aipw' pre_pngexp

set seed 20260819

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

local nboot_aipw = 1000

di as result _n "════════════════════════════════════════════════════════════"
di as result "PART 3 -- AIPW (Act 2): PNG debt/GDP exposure, non-default vs default-linked"
di as result "════════════════════════════════════════════════════════════"
di as result "  h   ND (se_boot)     DEF (se_boot)     def-nd   [95% boot CI]   Clogg z    p"

tempname A
tempfile aipwf
postfile `A' byte horizon double b_nd se_nd double b_def se_def ///
    double diff se_diff lo95 hi95 double cloggz cloggp ///
    long n_nd n_def nctry_nd nctry_def ntreat_nd ntreat_def using "`aipwf'", replace

forvalues h = 0/4 {
    _aipwpair, y(ch_pngexp_`h') ///
        d1(onset_def) if1(sample==1 & onset_nd==0) ///
        d2(onset_nd)  if2(sample==1 & onset_def==0) ///
        omod(`om_pngexp') pz(`om_pngexp' `cz_def') reps(`nboot_aipw')

    if r(ok) {
        local B1 = r(b1)
        local B2 = r(b2)
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
    label var se_nd   "Row-bootstrap SE, non-default"
    label var b_def   "Default-linked ATE"
    label var se_def  "Row-bootstrap SE, default-linked"
    label var diff    "def - nd difference"
    label var se_diff "Row-bootstrap SE of the difference"
    label var lo95    "95% bootstrap CI lower"
    label var hi95    "95% bootstrap CI upper"
    label var cloggz  "Clogg et al. (1995) z (companion, analytic SEs)"
    label var cloggp  "p-value of the Clogg z"
    label var n_nd    "Observations, non-default outcome-regression sample"
    label var n_def   "Observations, default-linked outcome-regression sample"
    label var nctry_nd  "Countries, non-default sample"
    label var nctry_def "Countries, default-linked sample"
    label var ntreat_nd  "Treated onsets, non-default sample"
    label var ntreat_def "Treated onsets, default-linked sample"
    export delimited "$tabs/png_gdp_exposure_aipw.csv", replace
restore
di as result "PART 3 COMPLETE -- AIPW results exported: $tabs/png_gdp_exposure_aipw.csv"

di as result _n "30_png_gdp_exposure_evolution.do complete. EXPLORATORY/ROBUSTNESS ONLY -- not"
di as result "wired into 00_master.do. Tests whether the EXPOSURE ITSELF moves differently by"
di as result "resolution type -- for whether exposure LEVEL predicts heterogeneous GDP/credit/"
di as result "investment outcomes, see 29_png_stock_gdp_exposure.do instead."
