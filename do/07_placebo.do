/*===========================================================================
  07_PLACEBO.DO
  Placebo Randomization Test — All Spread Crises (onset_all)

  Procedure:
    1. Draw 61 random country-year observations from the estimation sample
    2. Assign placebo onset dummy = 1 to those observations
    3. Run LP for h = 0..4 using the placebo dummy
    4. Repeat 1000 times → empirical null distribution of beta(h)
    5. Compare true beta(h) [from 02_lp_all.do] to the distribution
    6. Compute empirical p-values and plot

  Notes:
    - xtreg fe cluster(cid) used inside the loop (faster than xtscc)
      The goal is the null distribution of beta, not precise SEs per draw
    - xtscc with DK SE used for the true beta (already in irf_all.dta)
    - set seed ensures full replication
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking l_govexp l_open l_credit hyperinf_dummy"

* VIX and ust10y absorbed by year FE (i.year) — omitted.
local controls  $ctrl_core
local n_treat   = 61      // number of real onsets
local n_reps    = 1000
local seed      = 20250101

set seed `seed'

* ── Estimation sample index ───────────────────────────────────────────────
* Keep only rows eligible for placebo assignment:
*   - in estimation sample (continuation==0)
*   - non-missing on all controls and outcome dy_2 (h=2 as reference)
quietly {
    gen in_pool = (sample == 1) & !missing(dy_2)
    foreach v of local controls {
        replace in_pool = 0 if missing(`v')
    }
}

quietly count if in_pool == 1
local N_pool = r(N)
di as result "Pool size for placebo draws: `N_pool' obs"
di as result "True onsets to replicate:    `n_treat'"
di as result "Replications:                `n_reps'"

* ── Storage matrix: n_reps × 5 horizons ──────────────────────────────────
matrix PL = J(`n_reps', 5, .)

* ── Placebo loop ─────────────────────────────────────────────────────────
di as result _n "Running placebo loop (this takes a few minutes)..."

forvalues r = 1/`n_reps' {

    * Draw n_treat unique observations from the pool
    quietly {
        gen     _u      = runiform() if in_pool == 1
        egen    _rank   = rank(_u),   unique
        gen     placebo = (_rank <= `n_treat') & !missing(_rank)
        drop    _u _rank
    }

    * Run LP at each horizon using placebo onset
    forvalues h = 0/4 {
        quietly xtreg dy_`h' placebo `controls' i.year ///
            if in_pool == 1, fe vce(cluster cid)
        matrix PL[`r', `h'+1] = _b[placebo]
    }

    drop placebo

    if mod(`r', 100) == 0 di "  `r' / `n_reps' done"
}

di as result "Placebo loop complete."

* ── Save placebo distribution ─────────────────────────────────────────────
clear
svmat PL, names(col)
rename c1 b_h0
rename c2 b_h1
rename c3 b_h2
rename c4 b_h3
rename c5 b_h4
save "$clean/placebo_dist.dta", replace

* ── Load true betas from Act 1 results ───────────────────────────────────
use "$clean/irf_all.dta", clear
keep if horizon >= 0
* True betas: b at horizon 0..4 (rows 1..5 in irf_all, horizon==0..4)
forvalues h = 0/4 {
    quietly summarize b if horizon == `h'
    scalar true_b_h`h' = r(mean)
}

* ── Compute empirical p-values and percentiles ────────────────────────────
use "$clean/placebo_dist.dta", clear

di as result _n "=== PLACEBO TEST RESULTS ==="
di "Horizon  True_beta   p5_null   p50_null  p95_null  Emp_pval"
di "{hline 65}"

matrix PVAL = J(5, 4, .)   // cols: true_b, p5, p50, emp_pval

forvalues h = 0/4 {
    quietly summarize b_h`h', detail
    local p5   = r(p5)
    local p50  = r(p50)
    local p95  = r(p95)

    * Empirical p-value: share of placebo draws <= true beta (left-tail)
    quietly count if b_h`h' <= true_b_h`h'
    local emp_p = r(N) / `n_reps'

    matrix PVAL[`h'+1, 1] = true_b_h`h'
    matrix PVAL[`h'+1, 2] = `p5'
    matrix PVAL[`h'+1, 3] = `p50'
    matrix PVAL[`h'+1, 4] = `emp_p'

    di "h=" `h' "        " %7.3f true_b_h`h' ///
       "     " %7.3f `p5'  ///
       "     " %7.3f `p50' ///
       "     " %7.3f `p95' ///
       "     " %5.3f `emp_p'
}

* ── Figure 5: Placebo distribution vs. true beta ─────────────────────────
* One panel per horizon h=0..4
* Use xline() to mark true beta — simpler and more reliable than pci

local c_main "23 55 94"
local c_null "180 180 180"

forvalues h = 0/4 {
    * Extract true beta, emp p-value and null distribution bounds
    local tb  = PVAL[`h'+1, 1]
    local ep  = PVAL[`h'+1, 4]
    local tb_str  = string(`tb',  "%5.3f")
    local ep_str  = string(`ep',  "%5.3f")

    * Compute x-axis bounds: include true beta with 0.5pp margin on each side
    quietly summarize b_h`h'
    local x_min = min(`tb', r(min)) - 0.5
    local x_max = max(r(max), 0) + 0.5

    * Round to nearest 0.5 for clean axis labels
    local x_min = floor(`x_min' * 2) / 2
    local x_max = ceil(`x_max' * 2) / 2

    twoway ///
        (histogram b_h`h', width(0.3) freq start(`x_min')          ///
            fcolor("`c_null'%50") lcolor("`c_null'") lwidth(vthin)), ///
        xline(`tb', lcolor("`c_main'") lwidth(medthick) lpattern(solid)) ///
        xscale(range(`x_min' `x_max'))                               ///
        xlabel(`x_min'(0.5)`x_max', labsize(small))                  ///
        xtitle("Placebo beta (pp)", size(small))                      ///
        ytitle("Frequency", size(small))                              ///
        title("Placebo Distribution - h=`h'", size(medium))          ///
        subtitle("Vertical line = true beta. `n_reps' random draws.", ///
                 size(small))                                         ///
        note("True beta = `tb_str'. Emp. p-value = `ep_str'",        ///
             size(vsmall))                                            ///
        legend(off)                                                   ///
        graphregion(color(white)) plotregion(color(white))

    graph export "$figs/fig5_placebo_h`h'.pdf", replace
}

di as result _n "Placebo figures saved (fig5_placebo_h0 to h4)."

* ── Export p-value table ──────────────────────────────────────────────────
* Use svmat to convert PVAL matrix directly — avoids repeated gen in loop
clear
svmat PVAL, names(col)
rename c1 true_b
rename c2 p5_null
rename c3 p50_null
rename c4 emp_pval
gen horizon = _n - 1
order horizon true_b p5_null p50_null emp_pval
export delimited "$tabs/placebo_pvalues.csv", replace
di as result "Placebo p-values saved: $tabs/placebo_pvalues.csv"
