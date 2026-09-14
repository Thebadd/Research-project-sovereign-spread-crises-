/*===========================================================================
  29_PNG_STOCK_GDP_EXPOSURE.DO
  EXPLORATORY / ROBUSTNESS-ONLY TEST -- NOT part of the headline six-variable
  design and NOT wired into 00_master.do, matching this project's own
  convention for standalone exploratory files (27_png_debt_test.do,
  13e_nexus_bars.do, etc.).

  MOTIVATION. Distinct from 27/28's FLOW-based test (net PNG debt flows /
  GNI, abandoned after a robust null on both OLS and AIPW). This file goes
  back to the PNG debt STOCK (as originally used in 27's retired first
  version), but scales it by GDP instead of GNI, and tests a DIFFERENT
  hypothesis shape entirely: not "does the channel itself move by
  resolution type," but "does pre-crisis private-sector external-debt
  EXPOSURE amplify how much a spread crisis costs GDP" -- i.e. an
  amplifier/heterogeneity design, the same TEMPLATE as
  13d_aipw_nexus_split.do's sovereign-bank-nexus split (mirroring Asonuma
  et al.'s own domestic-aggregate-demand split, their eqs. 5-6), but with
  PNG debt/GDP as the amplifier instead of claimsgov_assets, and BOTH
  resolution types in scope (not ND-only, unlike 28's first-look design).

  EXPOSURE MEASURE: png_gdp = PNG debt stock (t-1) / GDP (t-1) * 100.
  A STOCK ratio, matching this project's other pre-crisis exposure
  measures (a_nexus in 13d, debt-to-GDP conventions elsewhere) -- NOT the
  flow measure retired in 27/28. This directly answers "how much of the
  private sector's own balance sheet is externally, unconditionally
  financed," independent of the sovereign's own debt level.

  DENOMINATOR CAVEAT, STATED PLAINLY: this project has no WDI current-US$
  GDP series already imported (only `gdp_real`, inherited from the
  original Excel skeleton import, whose exact price/currency basis is not
  independently re-verified here). `gdp_real` is used because it is
  ALREADY the denominator this project uses to recover levels from WDI
  %-of-GDP ratios elsewhere (ln_r_credit, as_r_claims_govt in
  18_transforms.do multiply WDI ratio variables by gdp_real to get back a
  level) -- so using it here is consistent with that existing convention,
  not a new assumption. But PNG debt stock is reported in CURRENT US$,
  and if `gdp_real` turns out to be on a different basis (constant-price,
  or a different currency convention) than the WDI ratios it is already
  multiplied against elsewhere, this ratio could be off by a scale factor.
  CHECK THE MAGNITUDE of png_gdp once run (should mostly sit in a
  plausible 0-100%-ish range for external private debt/GDP) before
  trusting it -- flagged as an open verification, not resolved here.

  MEDIAN SPLIT: computed ONCE, pooling all onsets (onset_all==1 &
  sample==1), not separately by resolution type -- matching Asonuma et
  al.'s own DAD-split convention (one threshold, comparable across nd/def)
  per their eqs. (5)-(6), and per the explicit recommendation this
  approach was proposed with. Country-mean-filled where the t-1 value is
  missing (same coverage-thinness fallback as 13d's a_nexus).

  DESIGN, mirrors 13d_aipw_nexus_split.do's own high/low-split machinery
  EXACTLY, but OLS not AIPW (this is a first pass -- AIPW is a natural
  extension IF this shows a signal worth chasing, not built here, per the
  same "test first, extend if promising" discipline used for 27/28):
    Four cells: {nd,def} x {high,low png_gdp exposure}, each vs the
    tranquil pool, rival resolution type dropped from each cell's own
    sample -- same onset-tier design as 03_lp_resolution.do/
    12_channels_resolution.do throughout this project.
    Within each resolution type, the HIGH-LOW difference is tested via a
    paired row-bootstrap (G=1000, seeded), stratified over {control,
    high-treated, low-treated} -- the same _aipwdiff-style mechanism 13d
    uses for its own high-low nexus contrast, adapted here to plain OLS
    (`regress`) instead of AIPW's `_aipw`. Clogg et al. (1995) z reported
    alongside as the permissive analytic-SE companion, same convention as
    every other difference test in this project.

  OUTCOME: GDP (dy_h, already built in panel_lp.dta) only, h=0..4 -- not
  the channels, per explicit scope ("not with the flow but with the stock
  and gdp"). A channel extension (credit, investment) is a natural next
  step, not built here.

  INTERPRETIVE CAUTION, stated per this project's own writing standard:
  png_gdp is NOT randomly assigned. Countries with deep private external
  financing access also tend to have deeper financial markets, higher
  income, and different institutional quality -- all things that
  plausibly affect crisis type and severity on their own. This design
  documents a HETEROGENEOUS ASSOCIATION, not an identified causal
  moderating effect, unless png_gdp can be argued predetermined and
  unconfounded conditional on $ctrl_core -- not argued here, stated
  plainly rather than left implicit. The same confound-check diagnostic
  13d prints (country composition, income proxy check) is reproduced
  below for the same reason.

  Output: $tabs/png_stock_gdp_exposure.csv (arm x level x horizon, levels +
          the within-arm high-low difference/CI/Clogg z) ;
          $figs/fig_png_stock_gdp_exposure.pdf (2-panel: nd, def; lines =
          high/low exposure, orange/green, matching 13d's own high/low
          palette).
  Run standalone, any time after 17_predictors.do and 18_transforms.do
  (needs $ctrl_core, gdp_real, dy_h, sample, onset_nd, onset_def).
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* SETUP + BUILD png_gdp (self-contained: re-imports the PNG debt STOCK from
* data/raw/PNGtoGNI.xlsx, the SAME raw file 27_png_debt_test.do's retired
* stock-based version used -- DT.DOD.DPNG.CD -- but scaled by gdp_real
* instead of GNI, so it does not depend on 27/28 having run first)
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

use "$clean/panel_lp.dta", clear
capture drop pngdebt_stock
merge m:1 iso3 year using `png_cy', keep(master match) nogen
sort cid year
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

capture drop png_gdp
gen double png_gdp = pngdebt_stock / gdp_real * 100 ///
    if pngdebt_stock >= 0 & gdp_real > 0 & !missing(pngdebt_stock, gdp_real)
label var png_gdp "Private nonguaranteed external debt stock / GDP, pct (exposure amplifier)"

quietly summarize png_gdp if sample==1
di as result _n "  png_gdp coverage/range check (sample==1): N=" r(N) "  mean=" %6.2f r(mean) ///
    "  min=" %6.2f r(min) "  max=" %6.2f r(max)
di as result "  (SANITY CHECK: should mostly sit in a plausible 0-100%%-ish range for external"
di as result "   private debt/GDP -- see header DENOMINATOR CAVEAT if this looks implausible.)"

* ══════════════════════════════════════════════════════════════════════════
* AMPLIFIER: pre-crisis exposure, entry level (t-1, country-mean filled) +
*   ONE median split pooling ALL onsets (onset_all==1), matching Asonuma et
*   al.'s own DAD-threshold convention -- see header.
* ══════════════════════════════════════════════════════════════════════════
capture drop a_pngexp a_pngexp_cm high_png
gen a_pngexp = L.png_gdp
bysort cid: egen a_pngexp_cm = mean(png_gdp)
replace a_pngexp = a_pngexp_cm if missing(a_pngexp)

quietly summarize a_pngexp if sample==1 & onset_all==1, detail
local med = r(p50)
gen high_png = (a_pngexp >= `med') if !missing(a_pngexp)
label define hpe 0 "Low PNG/GDP exposure" 1 "High PNG/GDP exposure"
label values high_png hpe

di as result _n "=== PNG/GDP EXPOSURE MEDIAN SPLIT (ALL ONSETS, ONE POOLED THRESHOLD) ==="
di as result "  Amplifier = PNG debt stock / GDP (pre-crisis, country-mean filled)"
di as result "  Median cutoff among ALL crisis onsets = " %6.2f `med'
foreach t in onset_all onset_nd onset_def {
    quietly count if sample==1 & `t'==1 & high_png==1
    local nh = r(N)
    quietly count if sample==1 & `t'==1 & high_png==0
    local nl = r(N)
    quietly count if sample==1 & `t'==1 & missing(high_png)
    local nm = r(N)
    di as result "  `t': high=" `nh' "  low=" `nl' "  unclassified=" `nm'
}

* ── Confound check: which countries fall in each bin? ────────────────────
di as result _n "=== PNG-EXPOSURE-BIN COUNTRY COMPOSITION (all crisis onsets) ==="
capture noisily tabulate country high_png if sample==1 & onset_all==1, row nofreq
di as result _n "  Mean pre-crisis PNG/GDP exposure by bin, over all onsets:"
capture noisily tabstat a_pngexp if sample==1 & onset_all==1, ///
    by(high_png) statistics(mean min max n) format(%6.2f)
di as result "  (Read alongside any income/development ranking of these countries -- see header's"
di as result "   interpretive caution: this split is NOT randomly assigned.)"

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM -- HIGH-LOW-within-arm difference, paired row bootstrap + Clogg z.
*   Mirrors 13d_aipw_nexus_split.do's _aipwdiff mechanism exactly, adapted
*   to plain OLS (`regress`) instead of AIPW's `_aipw`.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _lpdiffboot_hl
program define _lpdiffboot_hl, rclass
    syntax , Y(string) D(string) IFHIGH(string) IFLOW(string) ///
             CTRL(string) REPS(integer)

    capture xtreg `y' `d' `ctrl' if `ifhigh', fe vce(robust)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bh   = _b[`d']
    local ah   = _se[`d']
    local nh   = e(N)

    capture xtreg `y' `d' `ctrl' if `iflow', fe vce(robust)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bl   = _b[`d']
    local al   = _se[`d']
    local nl   = e(N)

    local dh = `bh' - `bl'

    capture drop _pool
    quietly gen byte _pool = 0 if (`ifhigh') | (`iflow')
    quietly replace _pool = 1 if `d' == 1 & (`ifhigh')
    quietly replace _pool = 2 if `d' == 1 & (`iflow')

    tempname pf
    tempfile bf
    quietly postfile `pf' double th double tl double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture regress `y' `d' `ctrl' i.cid if `ifhigh', vce(robust)
            local th = cond(_rc==0, _b[`d'], .)
            capture regress `y' `d' `ctrl' i.cid if `iflow', vce(robust)
            local tl = cond(_rc==0, _b[`d'], .)
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

    local cloggz = .
    local cloggp = .
    if !missing(`ah') & !missing(`al') & (`ah'^2 + `al'^2) > 0 {
        local cloggz = `dh' / sqrt(`ah'^2 + `al'^2)
        local cloggp = 2*(1 - normal(abs(`cloggz')))
    }

    return scalar ok    = 1
    return scalar bh    = `bh'
    return scalar bl    = `bl'
    return scalar ah    = `ah'
    return scalar al    = `al'
    return scalar nh    = `nh'
    return scalar nl    = `nl'
    return scalar bseh  = `bseh'
    return scalar bsel  = `bsel'
    return scalar dh    = `dh'
    return scalar se    = `se'
    return scalar lo    = `lo'
    return scalar hi    = `hi'
    return scalar nboot = `nd'
    return scalar cloggz = `cloggz'
    return scalar cloggp = `cloggp'
end

* ══════════════════════════════════════════════════════════════════════════
* ESTIMATE -- GDP, four cells {nd,def} x {high,low}, within-arm difference
* ══════════════════════════════════════════════════════════════════════════
set seed 20260819
local nboot = 1000

tempname R
tempfile resf
postfile `R' str4 part str4 level byte horizon double b se lo hi ntreat long nobs using "`resf'", replace

tempname D
tempfile diff_resf
postfile `D' str4 part byte horizon double dhl bhi blo se lo hi nd double cloggz cloggp using "`diff_resf'", replace

di as result _n "############### OUTCOME: GDP -- PNG/GDP exposure split, high vs low ###############"

foreach cell in "nd onset_nd onset_def" "def onset_def onset_nd" {
    gettoken part cell : cell
    gettoken Dv   cell : cell
    gettoken riv  cell : cell

    di as result _n "--- GDP | Part `part' ---"
    di as result "    h   LOW (se_boot)     HIGH (se_boot)    high-low  [95% boot CI]   Clogg z    p"

    post `R' ("`part'") ("low")  (0) (0) (0) (0) (0) (0) (0)
    post `R' ("`part'") ("high") (0) (0) (0) (0) (0) (0) (0)
    post `D' ("`part'") (0) (0) (0) (0) (0) (0) (0) (.) (.)

    forvalues h = 0/4 {
        local ifhigh sample==1 & `riv'==0 & (`Dv'==0 | (`Dv'==1 & high_png==1))
        local iflow  sample==1 & `riv'==0 & (`Dv'==0 | (`Dv'==1 & high_png==0))

        quietly count if `Dv'==1 & high_png==1 & sample==1 & `riv'==0
        local ntrh = r(N)
        quietly count if `Dv'==1 & high_png==0 & sample==1 & `riv'==0
        local ntrl = r(N)

        _lpdiffboot_hl, y(dy_`h') d(`Dv') ifhigh(`ifhigh') iflow(`iflow') ///
            ctrl(`ctrl_core') reps(`nboot')

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
            local ND = r(nboot)
            local NH = r(nh)
            local NL = r(nl)

            post `R' ("`part'") ("low")  (`h'+1) (`BL') (`BSEL') (`BL'-1.96*`BSEL') (`BL'+1.96*`BSEL') (`ntrl') (`NL')
            post `R' ("`part'") ("high") (`h'+1) (`BH') (`BSEH') (`BH'-1.96*`BSEH') (`BH'+1.96*`BSEH') (`ntrh') (`NH')

            local zz = r(cloggz)
            local pz = r(cloggp)
            post `D' ("`part'") (`h'+1) (`DH') (`BH') (`BL') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')

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
               " " %7.3f `zz' " " %5.3f `pz'
        }
        else di as error "    h=" `h'+1 ": estimate failed (too thin)."
    }
}
postclose `R'
postclose `D'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT -- CSV
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`diff_resf'", clear
    label var dhl "AIPW-style OLS (high - low PNG/GDP exposure) difference, GDP (pp)"
    label var bhi "High-exposure onset coefficient"
    label var blo "Low-exposure onset coefficient"
    label var lo  "95% CI lower (row bootstrap)"
    label var hi  "95% CI upper (row bootstrap)"
    label var nd  "Valid bootstrap draws"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
    label var sig95 "Bootstrap CI excludes 0 (governing test)"
    order part horizon dhl bhi blo se lo hi nd sig95 cloggz cloggp
    save "`diff_resf'", replace
restore

preserve
    use "`resf'", clear
    label var b  "OLS coefficient on onset dummy, GDP (pp)"
    label var se "Row-bootstrap SE"
    label var lo "95% CI lower = b - 1.96*se (bootstrap)"
    label var hi "95% CI upper = b + 1.96*se (bootstrap)"
    label var ntreat "Treated onsets in cell"
    label var nobs   "Observations in this cell's own regression sample"
    order part level horizon b se lo hi ntreat nobs
    tempfile levf
    save `levf'
restore

* Merge levels + difference into one long export for convenience.
preserve
    use `levf', clear
    tempfile levf2
    save `levf2'
    use "`diff_resf'", clear
    gen str4 level = "diff"
    rename dhl b
    append using `levf2'
    export delimited "$tabs/png_stock_gdp_exposure.csv", replace
restore
di as result _n "Results CSV saved: $tabs/png_stock_gdp_exposure.csv"

* ══════════════════════════════════════════════════════════════════════════
* FIGURE -- 2-panel (nd, def), high (orange) vs low (green) exposure
* ══════════════════════════════════════════════════════════════════════════
preserve
    use `levf', clear
    gen byte partid = 1 if part=="nd"
    replace partid = 2 if part=="def"
    label define pl2 1 "Non-default" 2 "Default-linked"
    label values partid pl2

    local c_hi "230 126 34"
    local c_lo "34 139 34"
    capture twoway ///
        (rarea lo hi horizon if level=="high", color("`c_hi'%16") lwidth(none)) ///
        (rarea lo hi horizon if level=="low",  color("`c_lo'%16") lwidth(none)) ///
        (connected b horizon if level=="high", lcolor("`c_hi'") lwidth(medthick) msymbol(square)) ///
        (connected b horizon if level=="low",  lcolor("`c_lo'") lwidth(medthick) msymbol(circle)), ///
        by(partid, yrescale legend(off) note("") graphregion(color(white)) title("GDP", size(medlarge) color(navy))) ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)5, labsize(medium)) ylabel(, labsize(medium) angle(horizontal)) ///
        xtitle("Year", size(medium)) ///
        ytitle("Cumulative percent change", size(medsmall)) ///
        graphregion(color(white)) plotregion(color(white))
    if _rc == 0 {
        graph export "$figs/fig_png_stock_gdp_exposure.pdf", replace
        di as result "Figure saved: fig_png_stock_gdp_exposure.pdf"
    }
    else di as error "  ** figure failed (rc=" _rc ")"
restore

di as result _n "29_png_stock_gdp_exposure.do complete. EXPLORATORY/ROBUSTNESS ONLY -- not wired"
di as result "into 00_master.do or any headline channel list. See header's interpretive caution:"
di as result "this documents a heterogeneous association, not an identified causal effect."
