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

  DENOMINATOR -- RESOLVED, NOT `gdp_real`. The first version of this file
  used `gdp_real` (IMF WEO NGDP_R, already in the panel from 11_weo.do) as
  the denominator, on the reasoning that this project already multiplies
  WDI %-of-GDP ratios by gdp_real elsewhere to recover levels. That
  reasoning was WRONG for this specific use: `gdp_real`/`gdp_nominal` are
  both DOMESTIC-CURRENCY WEO series (confirmed directly from 11_weo.do's
  own indicator-mapping comment: NGDP_R/NGDP, "dom. cur."), whereas PNG
  debt stock is reported in raw CURRENT US$ -- dividing one by the other
  mixes currencies and produced nonsense when actually run (png_gdp came
  back in the billions/tens-of-billions instead of a plausible 0-100%-ish
  ratio, confirmed live). FIXED: this file now imports WEO's own NGDPD
  series (GDP, current prices, US DOLLARS, billions) directly from
  data/raw/WEOApr2026all.xlsx -- confirmed present in that file's own
  indicator list, not previously imported by any other file in this
  project -- and uses it (converted from billions to raw US$) as the
  denominator instead. Both sides of the ratio are now genuinely in the
  same currency.

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
    Cells: {nd [,def]} x {high,low png_gdp exposure}, each vs the tranquil
    pool, rival resolution type dropped from each cell's own sample --
    same onset-tier design as 03_lp_resolution.do/
    12_channels_resolution.do throughout this project.
    Within each resolution type, the HIGH-LOW difference is tested via a
    paired row-bootstrap (G=1000, seeded), stratified over {control,
    high-treated, low-treated} -- the same _aipwdiff-style mechanism 13d
    uses for its own high-low nexus contrast, adapted here to plain OLS
    (`regress`) instead of AIPW's `_aipw`. Clogg et al. (1995) z reported
    alongside as the permissive analytic-SE companion, same convention as
    every other difference test in this project.

  SCOPE, PER EXPLICIT REQUEST (speed + priority): the DEFAULT-LINKED arm
  is MUTED (commented out) in the estimation loop's `foreach cell in ...`
  list -- runs currently estimate the ND arm only. Restore it by
  uncommenting that line once the ND-only read is done.

  OUTCOME: credit (bank credit to the private sector, ch_credit_h) and inv
  (investment, ch_inv_h) -- built fresh here (not persisted in
  panel_lp.dta), same log-real-level differenced construction as
  13c_aipw_channels.do/13d_aipw_nexus_split.do -- PLUS GDP (dy_h, already
  in the panel), all h=0..4. Credit and investment were added, and the
  loop reordered to run them FIRST, per explicit priority: "je souhaite
  voir en priorité comment se comporte le credit au secteur prive et
  l'investment chez les nd spread crises en fonction de si on a une low
  or high exposure."

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

  Output: $tabs/png_stock_gdp_exposure.csv (outcome x arm x level x
          horizon, levels + the within-arm high-low difference/CI/Clogg
          z) ; $figs/fig_png_stock_gdp_<outcome>.pdf, one per outcome
          (credit, inv, gdp), ND arm only while def stays muted -- high
          (orange)/low (green) exposure lines, no by()-panel needed with
          a single resolution type in scope.
  Run standalone, any time after 17_predictors.do and 18_transforms.do
  (needs $ctrl_core, ln_r_credit, ln_r_inv, dy_h, sample, onset_nd,
  onset_def).
===========================================================================*/

* ══════════════════════════════════════════════════════════════════════════
* SETUP + BUILD png_gdp (self-contained: re-imports the PNG debt STOCK from
* data/raw/PNGtoGNI.xlsx -- DT.DOD.DPNG.CD -- AND a GDP-in-current-US$
* series from the IMF WEO raw file already used by 11_weo.do.
*
* DENOMINATOR FIX (supersedes the header's original "use gdp_real" plan,
* which was tested and found WRONG): gdp_real/gdp_nominal (11_weo.do) are
* IMF WEO's NGDP_R/NGDP, both DOMESTIC-CURRENCY series (confirmed from
* 11_weo.do's own indicator-mapping comment) -- dividing a raw current-US$
* PNG debt stock by a domestic-currency GDP figure produced nonsense
* (checked live: png_gdp came back in the billions/tens-of-billions
* instead of a 0-100%-ish ratio). The correct series is WEO's NGDPD (GDP,
* current prices, US DOLLARS, reported in billions) -- confirmed present
* in data/raw/WEOApr2026all.xlsx's own indicator list, not currently
* imported by any file in this project, so it is imported here directly,
* self-contained, matching 11_weo.do's own import mechanics (Excel letter
* columns AB..CA = years 1980..2031).
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

* -- GDP in current US$ (WEO NGDPD, billions) -- same import mechanics as
* 11_weo.do (Excel letter columns AB..CA = years 1980..2031).
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

* gdp_usd_bn is in US$ BILLIONS; pngdebt_stock is in raw current US$ (WDI/IDS
* convention) -- convert gdp_usd_bn to raw US$ (x 1e9) before dividing, so
* both sides of the ratio are on the same unit.
capture drop png_gdp
gen double png_gdp = pngdebt_stock / (gdp_usd_bn * 1e9) * 100 ///
    if pngdebt_stock >= 0 & gdp_usd_bn > 0 & !missing(pngdebt_stock, gdp_usd_bn)
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
* CHANNEL OUTCOMES -- credit, inv (ch_v_h = F h.ln_r_v - L.ln_r_v), built
*   fresh here (not persisted in panel_lp.dta), same construction as
*   13c_aipw_channels.do/13d_aipw_nexus_split.do. Added per explicit
*   request to prioritize credit-to-private-sector and investment under
*   the ND arm before extending to GDP/def.
* ══════════════════════════════════════════════════════════════════════════
foreach v in credit inv {
    local src ln_r_`v'
    capture drop `v'_base
    gen double `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen double ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    capture drop pre_`v'
    gen double pre_`v' = L.`src' - L2.`src'
}

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
* ESTIMATE -- {credit, inv, gdp} x {nd [,def]} x {high,low}, within-arm diff
*
* SCOPE, PER EXPLICIT REQUEST (speed + priority): the DEFAULT-LINKED arm is
* MUTED (commented out) below so a run only estimates the ND arm -- credit
* to the private sector and investment first, GDP kept in the loop too
* since it was already built and costs nothing extra to leave in. To
* restore the def arm, uncomment the "def onset_def onset_nd" line in the
* `foreach cell in ...' list below.
* ══════════════════════════════════════════════════════════════════════════
set seed 20260819
local nboot = 1000

tempname R
tempfile resf
postfile `R' str12 outcome str4 part str4 level byte horizon double b se lo hi ntreat long nobs using "`resf'", replace

tempname D
tempfile diff_resf
postfile `D' str12 outcome str4 part byte horizon double dhl bhi blo se lo hi nd double cloggz cloggp using "`diff_resf'", replace

foreach oc in "credit ch_credit" "inv ch_inv" "gdp dy" {
    gettoken ocl   oc : oc
    gettoken ystem oc : oc

    * Outcome-model controls: GDP uses $ctrl_core as-is; credit drops its
    * own domestic-credit control (l_credit_bank correlates 0.950 with the
    * credit outcome itself, matching 13c/13d's own established exception)
    * and adds pre_credit; inv adds pre_inv.
    if      "`ocl'" == "gdp"    local om $ctrl_core
    else if "`ocl'" == "credit" local om l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
    else                        local om $ctrl_core pre_`ocl'

    di as result _n "############### OUTCOME: `ocl' -- PNG/GDP exposure split, high vs low ###############"

    foreach cell in "nd onset_nd onset_def" /* "def onset_def onset_nd" -- MUTED, see header note above */ {
        gettoken part cell : cell
        gettoken Dv   cell : cell
        gettoken riv  cell : cell

        di as result _n "--- `ocl' | Part `part' ---"
        di as result "    h   LOW (se_boot)     HIGH (se_boot)    high-low  [95% boot CI]   Clogg z    p"

        post `R' ("`ocl'") ("`part'") ("low")  (0) (0) (0) (0) (0) (0) (0)
        post `R' ("`ocl'") ("`part'") ("high") (0) (0) (0) (0) (0) (0) (0)
        post `D' ("`ocl'") ("`part'") (0) (0) (0) (0) (0) (0) (0) (0) (.) (.)

        forvalues h = 0/4 {
            local ifhigh sample==1 & `riv'==0 & (`Dv'==0 | (`Dv'==1 & high_png==1))
            local iflow  sample==1 & `riv'==0 & (`Dv'==0 | (`Dv'==1 & high_png==0))

            quietly count if `Dv'==1 & high_png==1 & sample==1 & `riv'==0
            local ntrh = r(N)
            quietly count if `Dv'==1 & high_png==0 & sample==1 & `riv'==0
            local ntrl = r(N)

            _lpdiffboot_hl, y(`ystem'_`h') d(`Dv') ifhigh(`ifhigh') iflow(`iflow') ///
                ctrl(`om') reps(`nboot')

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

                post `R' ("`ocl'") ("`part'") ("low")  (`h'+1) (`BL') (`BSEL') (`BL'-1.96*`BSEL') (`BL'+1.96*`BSEL') (`ntrl') (`NL')
                post `R' ("`ocl'") ("`part'") ("high") (`h'+1) (`BH') (`BSEH') (`BH'-1.96*`BSEH') (`BH'+1.96*`BSEH') (`ntrh') (`NH')

                local zz = r(cloggz)
                local pz = r(cloggp)
                post `D' ("`ocl'") ("`part'") (`h'+1) (`DH') (`BH') (`BL') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')

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
}
postclose `R'
postclose `D'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT -- CSV (now outcome x part x level x horizon; part=="def" absent
*   from the data entirely while the def cell above stays muted)
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`diff_resf'", clear
    label var dhl "AIPW-style OLS (high - low PNG/GDP exposure) difference (pp)"
    label var bhi "High-exposure onset coefficient"
    label var blo "Low-exposure onset coefficient"
    label var lo  "95% CI lower (row bootstrap)"
    label var hi  "95% CI upper (row bootstrap)"
    label var nd  "Valid bootstrap draws"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
    label var sig95 "Bootstrap CI excludes 0 (governing test)"
    order outcome part horizon dhl bhi blo se lo hi nd sig95 cloggz cloggp
    save "`diff_resf'", replace
restore

preserve
    use "`resf'", clear
    label var b  "OLS coefficient on onset dummy (pp)"
    label var se "Row-bootstrap SE"
    label var lo "95% CI lower = b - 1.96*se (bootstrap)"
    label var hi "95% CI upper = b + 1.96*se (bootstrap)"
    label var ntreat "Treated onsets in cell"
    label var nobs   "Observations in this cell's own regression sample"
    order outcome part level horizon b se lo hi ntreat nobs
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
* FIGURE -- one per outcome (credit, inv, gdp), ND arm only while def stays
*   muted; high (orange) vs low (green) exposure. No by()-panel needed
*   since only one resolution type is estimated right now -- a plain
*   twoway per outcome, not 13d's by(partid) grid (that returns once the
*   def arm is restored).
* ══════════════════════════════════════════════════════════════════════════
preserve
    use `levf', clear
    local c_hi "230 126 34"
    local c_lo "34 139 34"
    foreach oc in credit inv gdp {
        if "`oc'" == "credit" local ptit "Bank credit"
        if "`oc'" == "inv"    local ptit "Investment"
        if "`oc'" == "gdp"    local ptit "GDP"
        local fnm "fig_png_stock_gdp_`oc'"
        capture twoway ///
            (rarea lo hi horizon if level=="high" & outcome=="`oc'" & part=="nd", color("`c_hi'%16") lwidth(none)) ///
            (rarea lo hi horizon if level=="low"  & outcome=="`oc'" & part=="nd", color("`c_lo'%16") lwidth(none)) ///
            (connected b horizon if level=="high" & outcome=="`oc'" & part=="nd", lcolor("`c_hi'") lwidth(medthick) msymbol(square)) ///
            (connected b horizon if level=="low"  & outcome=="`oc'" & part=="nd", lcolor("`c_lo'") lwidth(medthick) msymbol(circle)), ///
            yline(0, lpattern(dash) lcolor(gs8)) ///
            xlabel(0(1)5, labsize(medium)) ylabel(, labsize(medium) angle(horizontal)) ///
            xtitle("Year", size(medium)) ///
            ytitle("Cumulative percent change", size(medsmall)) ///
            title("`ptit' (non-default, PNG/GDP exposure split)", size(medlarge) color(navy)) ///
            legend(order(3 "High exposure" 4 "Low exposure") position(6) rows(1) size(small)) ///
            graphregion(color(white)) plotregion(color(white))
        if _rc == 0 {
            graph export "$figs/`fnm'.pdf", replace
            di as result "Figure saved: `fnm'.pdf"
        }
        else di as error "  ** `fnm' failed (rc=" _rc ")"
    }
restore

di as result _n "29_png_stock_gdp_exposure.do complete. EXPLORATORY/ROBUSTNESS ONLY -- not wired"
di as result "into 00_master.do or any headline channel list. See header's interpretive caution:"
di as result "this documents a heterogeneous association, not an identified causal effect."
