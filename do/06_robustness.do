/*===========================================================================
  06_ROBUSTNESS.DO
  Robustness checks for Act 1 (all crises) and Act 2 (resolution split)

  R1 — Threshold sensitivity
       a. Criterion 1 only (crit1_flag onset, 1000bps)
       b. Criterion 2 only (crit2_flag onset, HRT 378bps)

  R2 — Sample exclusions
       a. Drop Argentina (multiple extreme episodes, dominates LatAM)
       b. Drop Venezuela (extreme outlier, hyperinflation)
       c. Drop both Argentina and Venezuela
       d. Exclude COVID years 2020–2021

  R3 — Global crisis episodes
       a. Drop 2008–2009 (GFC)
       b. Drop 1997–1999 (Asian/EM crisis)

  R4 — Alternative SE: country clustering (xtreg fe cluster) vs DK

  Each robustness run stores IRF at h=0..4 in a summary matrix.
  Final table: robustness_summary.csv  (beta at h=2 across all specs)
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg debt ca banking_crisis l_govexp l_open l_credit hyperinf_dummy"

* VIX & ust10y absorbed by i.year — omitted from all LP specs
    local controls $ctrl_core
local horizons "0 1 2 3 4"

* Master storage: rows = specs, cols = horizons 0..4
* We store point estimates and 90% CI for h=0..4
local specs "baseline crit1 crit2 drop_arg drop_ven drop_both drop_covid drop_gfc drop_afc clust"
local nspecs : word count `specs'

matrix RB_b   = J(`nspecs', 5, .)
matrix RB_lo  = J(`nspecs', 5, .)
matrix RB_hi  = J(`nspecs', 5, .)

* Spec labels for output
local slabels `""Baseline (all, DK)" "Crit1 only (1000bps)" "Crit2 only (HRT)" "Drop Argentina" "Drop Venezuela" "Drop ARG+VEN" "Drop 2020-2021" "Drop 2008-2009" "Drop 1997-1999" "Country cluster""'

* ─────────────────────────────────────────────────────────────────────────
* Helper program: run LP for one spec, store in row `srow' of RB matrices
* Uses a pre-generated binary flag `_smpl' as the estimation sample.
* ─────────────────────────────────────────────────────────────────────────
capture program drop run_rob
program define run_rob
    args srow treatment

    * VIX & ust10y absorbed by i.year — omitted from all LP specs
    local controls $ctrl_core

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local col = `h' + 1

        capture xtscc dy_`h' `treatment' `controls' i.year ///
            if _smpl==1, fe lag(`lag')

        if _rc == 0 {
            matrix RB_b[`srow',  `col'] = _b[`treatment']
            matrix RB_lo[`srow', `col'] = _b[`treatment'] - 1.645*_se[`treatment']
            matrix RB_hi[`srow', `col'] = _b[`treatment'] + 1.645*_se[`treatment']
        }
        else {
            di as error "  run_rob: xtscc failed at h=`h', _rc=" _rc
        }
    }
end

* ══════════════════════════════════════════════════════════════════════════
* R0 — BASELINE (replicates 02_lp_all.do, all crises, DK SE)
* ══════════════════════════════════════════════════════════════════════════
di as result _n "R0: Baseline"
gen _smpl = sample
run_rob 1 onset_all

* ══════════════════════════════════════════════════════════════════════════
* R1 — THRESHOLD SENSITIVITY
* ══════════════════════════════════════════════════════════════════════════

* Need criterion-specific onset dummies
* Crit1 onset: crit1==1 at onset_all==1
gen onset_crit1 = (onset_all==1 & crit1==1)
* Crit2 onset: crit2==1 at onset_all==1
gen onset_crit2 = (onset_all==1 & crit2==1)

di as result _n "R1a: Criterion 1 only (1000bps)"
replace _smpl = sample
run_rob 2 onset_crit1

di as result _n "R1b: Criterion 2 only (HRT 378bps)"
replace _smpl = sample
run_rob 3 onset_crit2

* ══════════════════════════════════════════════════════════════════════════
* R2 — SAMPLE EXCLUSIONS
* ══════════════════════════════════════════════════════════════════════════

di as result _n "R2a: Drop Argentina"
replace _smpl = sample & (country != "Argentina")
run_rob 4 onset_all

di as result _n "R2b: Drop Venezuela"
replace _smpl = sample & (country != "Venezuela")
run_rob 5 onset_all

di as result _n "R2c: Drop Argentina + Venezuela"
replace _smpl = sample & !inlist(country, "Argentina", "Venezuela")
run_rob 6 onset_all

di as result _n "R2d: Drop 2020-2021 (COVID)"
replace _smpl = sample & !inrange(year, 2020, 2021)
run_rob 7 onset_all

* ══════════════════════════════════════════════════════════════════════════
* R3 — DROP GLOBAL CRISIS WINDOWS
* ══════════════════════════════════════════════════════════════════════════

di as result _n "R3a: Drop GFC window (2008-2009)"
replace _smpl = sample & !inrange(year, 2008, 2009)
run_rob 8 onset_all

di as result _n "R3b: Drop Asian/EM crisis window (1997-1999)"
replace _smpl = sample & !inrange(year, 1997, 1999)
run_rob 9 onset_all

* ══════════════════════════════════════════════════════════════════════════
* R4 — ALTERNATIVE SE: country clustering (xtreg fe cluster)
* ══════════════════════════════════════════════════════════════════════════

di as result _n "R4: Country-clustered SE (xtreg fe)"

replace _smpl = sample
forvalues h = 0/4 {
    local col = `h' + 1
    quietly xtreg dy_`h' onset_all `controls' i.year if _smpl==1, ///
        fe vce(cluster cid)
    matrix RB_b[10,  `col'] = _b[onset_all]
    matrix RB_lo[10, `col'] = _b[onset_all] - 1.645*_se[onset_all]
    matrix RB_hi[10, `col'] = _b[onset_all] + 1.645*_se[onset_all]
}
drop _smpl

* ══════════════════════════════════════════════════════════════════════════
* EXPORT ROBUSTNESS SUMMARY TABLE (beta at h=2 across all specs)
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== ROBUSTNESS SUMMARY: beta(h=2) across specifications ==="
di "Spec                          beta     [90% CI lo,  hi]"

local i = 1
foreach s of local specs {
    local label : word `i' of `slabels'
    di %-32s "`label'" ///
       " " %7.3f RB_b[`i',3] ///
       "  [" %7.3f RB_lo[`i',3] ", " %7.3f RB_hi[`i',3] "]"
    local ++i
}

* Save full robustness dataset
clear
local nspecs = 10
set obs `nspecs'
gen spec_num = _n
gen spec = ""
local i = 1
foreach s of local specs {
    replace spec = "`s'" in `i'
    local ++i
}

forvalues h = 0/4 {
    local col = `h' + 1
    gen b_h`h'   = .
    gen lo_h`h'  = .
    gen hi_h`h'  = .
    forvalues r = 1/10 {
        replace b_h`h'  = RB_b[`r',  `col'] in `r'
        replace lo_h`h' = RB_lo[`r', `col'] in `r'
        replace hi_h`h' = RB_hi[`r', `col'] in `r'
    }
}

export delimited "$tabs/robustness_summary.csv", replace
di as result "Robustness table saved: $tabs/robustness_summary.csv"

* ══════════════════════════════════════════════════════════════════════════
* ROBUSTNESS FIGURE: coefficient plot at h=2 across specs
* ══════════════════════════════════════════════════════════════════════════

gen spec_label = ""
replace spec_label = "Baseline (DK SE)"         in 1
replace spec_label = "Criterion 1 only (1000bps)" in 2
replace spec_label = "Criterion 2 only (HRT)"   in 3
replace spec_label = "Drop Argentina"            in 4
replace spec_label = "Drop Venezuela"            in 5
replace spec_label = "Drop ARG + VEN"            in 6
replace spec_label = "Excl. 2020–2021 (COVID)"  in 7
replace spec_label = "Excl. 2008–2009 (GFC)"    in 8
replace spec_label = "Excl. 1997–1999 (AFC)"    in 9
replace spec_label = "Country-clustered SE"      in 10

* Reverse order for reading top-to-bottom in plot
gen ypos = 11 - spec_num

* Define value label for y-axis ticks
label define spec_lbl                    ///
    1  "Country-clustered SE"            ///
    2  "Excl. 1997-1999 (AFC)"           ///
    3  "Excl. 2008-2009 (GFC)"           ///
    4  "Excl. 2020-2021 (COVID)"         ///
    5  "Drop ARG + VEN"                  ///
    6  "Drop Venezuela"                  ///
    7  "Drop Argentina"                  ///
    8  "Criterion 2 only (HRT)"          ///
    9  "Criterion 1 only (1000bps)"      ///
    10 "Baseline (DK SE)"
label values ypos spec_lbl

twoway ///
    (rcap lo_h2 hi_h2 ypos,                                                  ///
        horizontal lcolor("23 55 94") lwidth(medthick))                       ///
    (scatter ypos b_h2,                                                       ///
        msymbol(circle) mcolor("23 55 94") msize(medium)),                    ///
    xline(0, lpattern(dash) lcolor(gray) lwidth(thin))                        ///
    ylabel(1(1)10, valuelabel angle(0) labsize(small) noticks nogrid)         ///
    xlabel(, format(%4.1f) labsize(small))                                    ///
    xtitle("beta at h=2 (pp)", size(small))                                   ///
    ytitle("")                                                                 ///
    title("Robustness: Output Effect at h=2", size(medium))                   ///
    subtitle("All spread crises. 90% confidence intervals.", size(small))     ///
    legend(off)                                                                ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig4_robustness.pdf", replace
graph export "$figs/fig4_robustness.png", replace width(1200)
di as result "Figure 4 (robustness) saved."
