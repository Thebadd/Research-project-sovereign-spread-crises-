/*===========================================================================
  28_AIPW_EXTFIN_NEXUS_SPLIT_ND.DO
  EXPLORATORY / FIRST LOOK -- NOT wired into 00_master.do, not part of the
  headline six-variable design. Standalone, self-contained (rebuilds png_gni
  from its own raw import, does not require 27_png_debt_test.do to have
  been run first in the same session).

  HYPOTHESIS BEING TESTED. In non-default spread crises, bank credit to the
  private sector does not fall (a robust null result throughout this
  project's channel files) -- one candidate explanation is that firms shift
  from external to domestic financing when a crisis raises the cost/
  availability of foreign borrowing, cushioning domestic credit even as
  output slows. That substitution margin should be narrower for countries
  whose private sector was MORE dependent on external financing pre-crisis
  -- those firms have proportionally more external exposure to unwind and
  less domestic capacity to absorb it. This file tests that directly:
  median-split non-default onsets by pre-crisis private external-financing
  exposure (png_gni = private nonguaranteed external debt / GNI, built in
  27_png_debt_test.do), and compare AIPW-estimated GDP/credit/investment
  paths between the high- and low-exposure halves.

  SCOPE, DELIBERATELY NARROW FOR THIS FIRST PASS:
    - NON-DEFAULT ARM ONLY. No default-linked estimation, no def-vs-nd
      comparison -- "test it only for nd spread crises at first to see
      what is going on" (explicit user instruction). The default arm and
      a possible joint nexus x external-financing design (either a 2x2
      split or a triple interaction with the sovereign-bank nexus
      amplifier from 13d_aipw_nexus_split.do) are natural next steps if
      this shows a signal worth chasing, not built here.
    - SINGLE EXPOSURE. This is the one-amplifier version of the design
      discussed (external financing alone), not the joint nexus x
      external-financing test also discussed -- deliberately simplified so
      a first read doesn't depend on a fragile double median-split or a
      triple interaction with already-thin cells.
    - NO Table-3-style RTF export. CSVs, figures, and console output only
      -- that formatting effort is only worth building once/if this
      becomes a result worth reporting formally.

  COVERAGE CAVEAT, STATED UP FRONT, NOT DISCOVERED LATER: png_gni comes
  from World Bank IDS/WDI private-nonguaranteed-debt reporting, which has
  real gaps relative to claimsgov_assets' coverage (the amplifier used in
  13d_aipw_nexus_split.do) -- several onsets, default and non-default
  alike, may have no reported PNG debt at all. The coverage diagnostic
  below reports the actual usable non-default onset count before any
  estimate is read; if too few onsets survive in either bin, that is the
  finding of this file, not a bug to work around.

  DESIGN, COPIED VERBATIM FROM 13d_aipw_nexus_split.do (not reinvented):
    - _aipw: the same IPWRA engine (Eqs. 1-3) used throughout this
      project's AIPW files (08b/13c/13d/21/24/25) -- propensity-weighted
      outcome regression, doubly-robust, country FE.
    - _aipwdiff: the same HIGH-LOW-within-arm bootstrap program from 13d
      (row-level resampling stratified by {control, high-treated,
      low-treated} pool, G=1000 draws, seeded), reporting each level's own
      row-bootstrap SE (ADOPTED, matching 08b/13c/13d's departure from the
      analytic-SE formula, see those files' headers) plus the Clogg et al.
      (1995) z as an analytic-SE companion statistic on the HIGH-LOW
      difference.
    - cz_def predictor set: l_fedfunds l_contagion_dist_atdef
      years_since_def_onset -- the project-adopted nd/def-arm propensity
      predictors (08c_first_stage_table.do's head-to-head comparison).

  Output: $tabs/aipw_extfin_nexus_split_nd.csv (outcome x bank x horizon,
          levels) ; $tabs/aipw_extfin_nexus_diff_nd.csv (outcome x horizon,
          HIGH-LOW gap + bootstrap CI + Clogg z) ;
          $figs/fig_extfin_nexus_<outcome>.pdf (one per outcome: gdp,
          credit, inv -- two lines, high-extfin orange / low-extfin green,
          reusing 13d's own high/low palette since there is only one
          resolution-type panel here, not a by() grid).
  Run standalone, any time after 17_predictors.do (needs $ctrl_core,
  l_fedfunds, l_contagion_dist_atdef, years_since_def_onset, sample,
  onset_nd, onset_def, ln_r_credit, ln_r_inv).
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* SETUP + BUILD png_gni (self-contained, copied from 27_png_debt_test.do)
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
rename YR pngdebt
keep iso3 year pngdebt
tempfile t_png
save `t_png'

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

use "$clean/panel_lp.dta", clear
capture drop pngdebt gni
merge m:1 iso3 year using `png_cy', keep(master match) nogen
sort cid year
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

capture drop png_gni
gen double png_gni = pngdebt / gni * 100 if pngdebt >= 0 & gni > 0 & !missing(pngdebt, gni)
label var png_gni "Private nonguaranteed external debt / GNI, pct (external-financing exposure)"

* ══════════════════════════════════════════════════════════════════════════
* AMPLIFIER: pre-crisis external-financing exposure + median split OVER ND ONSETS
*   Median computed over onset_nd==1 only (not onset_all) -- this first pass
*   is ND-only, so the split threshold is a property of the population
*   actually being split, matching 13d_aipw_nexus_split.do's own logic.
* ══════════════════════════════════════════════════════════════════════════
capture drop a_extfin a_extfin_cm highextfin
gen a_extfin = L.png_gni                          // predetermined (year before onset)
bysort cid: egen a_extfin_cm = mean(png_gni)
replace a_extfin = a_extfin_cm if missing(a_extfin) // country-mean fill (thin coverage)

quietly summarize a_extfin if sample==1 & onset_nd==1, detail
local med = r(p50)
gen highextfin = (a_extfin >= `med') if !missing(a_extfin)
label define hef 0 "Low external financing" 1 "High external financing"
label values highextfin hef

di as result _n "=== EXTERNAL-FINANCING EXPOSURE MEDIAN SPLIT (NON-DEFAULT ONSETS ONLY) ==="
di as result "  Amplifier = private nonguaranteed external debt / GNI (pre-crisis, country-mean filled)"
di as result "  Median cutoff among NON-DEFAULT crisis onsets = " %6.2f `med'
quietly count if sample==1 & onset_nd==1 & highextfin==1
local nh = r(N)
quietly count if sample==1 & onset_nd==1 & highextfin==0
local nl = r(N)
quietly count if sample==1 & onset_nd==1 & missing(highextfin)
local nm = r(N)
di as result "  onset_nd: high=" `nh' "  low=" `nl' "  unclassified (no png_gni data)=" `nm'
di as result "  (COVERAGE CAVEAT: if `nm' is large relative to `nh'+`nl', png_gni's known"
di as result "   IDS/WDI reporting gaps are binding on this sample -- read the estimates below"
di as result "   with that in mind, not as a fully powered test.)"

* ── Confound check: which countries fall in each bin? ────────────────────
di as result _n "=== EXTFIN-BIN COUNTRY COMPOSITION (non-default crisis onsets) ==="
capture noisily tabulate country highextfin if sample==1 & onset_nd==1, ///
    row nofreq
di as result _n "  Mean pre-crisis external-financing exposure (PNG debt/GNI) by bin, over ND onsets:"
capture noisily tabstat a_extfin if sample==1 & onset_nd==1, ///
    by(highextfin) statistics(mean min max n) format(%6.2f)
di as result "  (Read alongside any income/development ranking of these countries: if the two"
di as result "   columns are not obviously split by development, the result is not merely a"
di as result "   development proxy.)"

* ══════════════════════════════════════════════════════════════════════════
* CHANNEL OUTCOMES — GDP (already dy_h), credit, inv
* ══════════════════════════════════════════════════════════════════════════
foreach v in credit inv {
    local src ln_r_`v'
    capture drop `v'_base
    gen `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    capture drop pre_`v'
    gen pre_`v' = L.`src' - L2.`src'
}

* AIPW outcome-model core, matching 13d_aipw_nexus_split.do's core_aipw.
local core_aipw l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2
local cz_def l_fedfunds l_contagion_dist_atdef years_since_def_onset

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — _aipw (Eqs. 1-3), copied verbatim from 13d_aipw_nexus_split.do
* ══════════════════════════════════════════════════════════════════════════
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

    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta  = `th'
    return scalar N      = `nn'
    return scalar nctry  = `nctry'
    return scalar se     = `sean'
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — HIGH - LOW extfin DIFFERENCE within the nd arm, copied verbatim
*   (structure) from 13d_aipw_nexus_split.do's _aipwdiff.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipwdiff
program define _aipwdiff, rclass
    syntax anything, IFCH(string) IFCL(string) OMOD(string) PZ(string) REPS(integer)
    gettoken yv Dv : anything
    capture _aipw `yv' `Dv' if `ifch', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bh = r(theta)
    local ah = r(se)
    local nh = r(N)
    local nctryh = r(nctry)
    capture _aipw `yv' `Dv' if `ifcl', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bl = r(theta)
    local al = r(se)
    local nl = r(N)
    local nctryl = r(nctry)
    local dh = `bh' - `bl'

    capture drop _pool
    quietly gen byte _pool = 0 if (`ifch') | (`ifcl')
    quietly replace _pool = 1 if `Dv' == 1 & (`ifch')
    quietly replace _pool = 2 if `Dv' == 1 & (`ifcl')

    tempname pf
    tempfile bf
    quietly postfile `pf' double th double tl double diff using "`bf'", replace
    forvalues b = 1/1000 {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture _aipw `yv' `Dv' if `ifch', omodel(`omod') pmodel(`pz') fe(cid)
            local th = cond(_rc==0, r(theta), .)
            capture _aipw `yv' `Dv' if `ifcl', omodel(`omod') pmodel(`pz') fe(cid)
            local tl = cond(_rc==0, r(theta), .)
            if !missing(`th') & !missing(`tl') quietly post `pf' (`th') (`tl') (`th' - `tl')
        restore
    }
    quietly postclose `pf'
    capture drop _pool
    local se = .
    local lo = .
    local hi = .
    local nd = 0
    local bseh = .
    local bsel = .
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
            quietly summarize th
            local bseh = r(sd)
            quietly summarize tl
            local bsel = r(sd)
        }
    restore
    return scalar ok = 1
    return scalar dh = `dh'
    return scalar bh = `bh'
    return scalar bl = `bl'
    return scalar ah = `ah'
    return scalar al = `al'
    return scalar nh = `nh'
    return scalar nl = `nl'
    return scalar nctryh = `nctryh'
    return scalar nctryl = `nctryl'
    return scalar bseh = `bseh'
    return scalar bsel = `bsel'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* REPRODUCIBILITY: same seed convention as 13d_aipw_nexus_split.do
* ══════════════════════════════════════════════════════════════════════════
set seed 20260819

* ══════════════════════════════════════════════════════════════════════════
* ESTIMATE — high/low extfin subsamples, NON-DEFAULT ARM ONLY
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str18 outcome str4 bank byte horizon double b se lo hi ntreat ///
    long nobs byte nctry using "`resf'", replace

tempname D
tempfile diff_resf
postfile `D' str18 outcome byte horizon double dhl bhi blo se lo hi nd ///
    double cloggz double cloggp using "`diff_resf'", replace

foreach oc in "gdp dy" "credit ch_credit" "inv ch_inv" {
    gettoken ocl   oc : oc
    gettoken ystem oc : oc

    if "`ocl'" == "gdp" local om `core_aipw'
    else                local om `core_aipw' pre_`ocl'

    di as result _n "############### OUTCOME: `ocl' (non-default, high vs low extfin) ###############"
    di as result "    h   LOW (se_boot)     HIGH (se_boot)    high-low  [95% boot CI]   Clogg z    p"
    di as result "        se_boot = ROW-BOOTSTRAP SE (ADOPTED, matching 08b/13c/13d's own convention)."
    di as result "        LOW/HIGH stars: conventional t-test vs zero (b/se_boot): * p<.10 ** p<.05 *** p<.01."

    post `R' ("`ocl'") ("low")  (0) (0) (0) (0) (0) (0) (0) (0)
    post `R' ("`ocl'") ("high") (0) (0) (0) (0) (0) (0) (0) (0)
    post `D' ("`ocl'") (0) (0) (0) (0) (0) (0) (0) (.) (.)

    forvalues h = 0/4 {
        local ifch sample==1 & onset_def==0 & (onset_nd==0 | (onset_nd==1 & highextfin==1))
        local ifcl sample==1 & onset_def==0 & (onset_nd==0 | (onset_nd==1 & highextfin==0))

        quietly count if onset_nd==1 & highextfin==1 & sample==1
        local ntrh = r(N)
        quietly count if onset_nd==1 & highextfin==0 & sample==1
        local ntrl = r(N)

        _aipwdiff `ystem'_`h' onset_nd, ifch(`ifch') ifcl(`ifcl') ///
            omod(`om') pz(`om' `cz_def') reps(1000)
        if r(ok) {
            local BH = r(bh)
            local BL = r(bl)
            local AH = r(ah)
            local AL = r(al)
            local BSEH = r(bseh)
            local BSEL = r(bsel)
            local DH = r(dh)
            local SE = r(se)
            local LO = r(lo)
            local HI = r(hi)
            local ND = r(nd)
            local NL = r(nl)
            local NH = r(nh)
            local NCL = r(nctryl)
            local NCH = r(nctryh)

            post `R' ("`ocl'") ("low")  (`h'+1) (`BL') (`BSEL') (`BL'-1.96*`BSEL') (`BL'+1.96*`BSEL') (`ntrl') (`NL') (`NCL')
            post `R' ("`ocl'") ("high") (`h'+1) (`BH') (`BSEH') (`BH'-1.96*`BSEH') (`BH'+1.96*`BSEH') (`ntrh') (`NH') (`NCH')

            local zz = .
            local pz2 = .
            if !missing(`AH') & !missing(`AL') & (`AH'^2 + `AL'^2) > 0 {
                local zz  = `DH' / sqrt(`AH'^2 + `AL'^2)
                local pz2 = 2*(1 - normal(abs(`zz')))
            }
            post `D' ("`ocl'") (`h'+1) (`DH') (`BH') (`BL') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz2')

            local tlo  = cond(`BSEL'>0, `BL'/`BSEL', .)
            local plo  = cond(!missing(`tlo'), 2*(1-normal(abs(`tlo'))), .)
            local sglo = cond(missing(`plo'), "", cond(`plo'<.01,"***",cond(`plo'<.05,"**",cond(`plo'<.10,"*",""))))
            local thi  = cond(`BSEH'>0, `BH'/`BSEH', .)
            local phi  = cond(!missing(`thi'), 2*(1-normal(abs(`thi'))), .)
            local sghi = cond(missing(`phi'), "", cond(`phi'<.01,"***",cond(`phi'<.05,"**",cond(`phi'<.10,"*",""))))

            local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
            di "    " %1.0f `h'+1 "  " %8.3f `BL' "`sglo'" " (" %5.3f `BSEL' ")  " ///
               %8.3f `BH' "`sghi'" " (" %5.3f `BSEH' ")  " %8.3f `DH' ///
               " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
               " " %7.3f `zz' " " %5.3f `pz2'
        }
        else di as error "    h=" `h'+1 ": estimate failed (too thin)."
    }
}
postclose `R'
postclose `D'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — CSVs
* ══════════════════════════════════════════════════════════════════════════
use "`diff_resf'", clear
label var dhl "AIPW (high - low external financing) difference (pp)"
label var bhi "High-extfin ATE"
label var blo "Low-extfin ATE"
label var lo  "95% CI lower (row bootstrap)"
label var hi  "95% CI upper (row bootstrap)"
label var nd  "Valid bootstrap draws"
label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
label var cloggp "p-value of the Clogg z"
gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
label var sig95 "Bootstrap CI excludes 0 (governing test)"
order outcome horizon dhl bhi blo se lo hi nd sig95 cloggz cloggp
export delimited "$tabs/aipw_extfin_nexus_diff_nd.csv", replace
di as result _n "Extfin high-low DIFFERENCE (ND arm) CSV saved: $tabs/aipw_extfin_nexus_diff_nd.csv"

use "`resf'", clear
label var b  "AIPW ATE on outcome (pp)"
label var se "Row-bootstrap SE (ADOPTED)"
label var lo "95% CI lower = b - 1.96*se (bootstrap)"
label var hi "95% CI upper = b + 1.96*se (bootstrap)"
label var ntreat "Treated ND onsets in cell"
label var nobs "Observations in this cell's own AIPW outcome-regression sample"
label var nctry "Countries in this cell's own AIPW outcome-regression sample"
order outcome bank horizon b se lo hi ntreat nobs nctry
export delimited "$tabs/aipw_extfin_nexus_split_nd.csv", replace
di as result _n "Extfin-split AIPW results (ND arm) CSV saved: $tabs/aipw_extfin_nexus_split_nd.csv"

* ══════════════════════════════════════════════════════════════════════════
* FIGURES — one per outcome, high (orange) vs low (green) extfin, ND arm only
* ══════════════════════════════════════════════════════════════════════════
local c_hi "230 126 34"   // high extfin = orange (matches 13d's high/low palette)
local c_lo "34 139 34"    // low  extfin = green
foreach oc in gdp credit inv {
    if "`oc'" == "gdp"    local ptit "GDP"
    if "`oc'" == "credit" local ptit "Bank credit"
    if "`oc'" == "inv"    local ptit "Investment"
    local fnm "fig_extfin_nexus_`oc'"

    capture twoway ///
        (rarea lo hi horizon if bank=="high" & outcome=="`oc'", color("`c_hi'%16") lwidth(none)) ///
        (rarea lo hi horizon if bank=="low"  & outcome=="`oc'", color("`c_lo'%16") lwidth(none)) ///
        (connected b horizon if bank=="high" & outcome=="`oc'", lcolor("`c_hi'") lwidth(medthick) msymbol(square)) ///
        (connected b horizon if bank=="low"  & outcome=="`oc'", lcolor("`c_lo'") lwidth(medthick) msymbol(circle)), ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)5, labsize(medium)) ylabel(, labsize(medium) angle(horizontal)) ///
        xtitle("Year", size(medium)) ///
        ytitle("Cumulative percent change", size(medsmall)) ///
        title("`ptit' (non-default, by external-financing exposure)", size(medlarge) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white))
    if _rc == 0 {
        graph export "$figs/`fnm'.pdf", replace
        di as result "Figure saved: `fnm'.pdf"
    }
    else di as error "  ** `fnm' failed (rc=" _rc ")"
}

di as result _n "28_aipw_extfin_nexus_split_nd.do complete. EXPLORATORY/FIRST LOOK -- ND arm only,"
di as result "single exposure (external financing). Read the credit and investment high-low gaps"
di as result "against the coverage caveat printed above before treating either as a real signal."
