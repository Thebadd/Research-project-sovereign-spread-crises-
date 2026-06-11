/*===========================================================================
  10_HETEROGENEITY.DO
  Heterogeneity Analysis — Two Dimensions

  BLOCK A: Financial Development
    Frontier markets (MSCI/FM classification) vs. Established EM
    onset_frontier  = onset in frontier market country  (subset of onset_all)
    onset_em        = onset in non-frontier EM country  (subset of onset_all)

  BLOCK B: Episode Duration as Moderator
    Duration = number of years the episode spans (onset + continuation years)
    Short = duration <= 2 years    (median split, computed from ep_id)
    Long  = duration >= 3 years
    onset_short / onset_long constructed at onset year only

  For each block:
    Spec A — separate LP for each subgroup (DK SE, lag = max(1,h+1))
    Spec B — joint regression + equality test H0: beta_1(h) = beta_2(h)

  Saves IRF datasets to $clean/ and figures to $figs/
===========================================================================*/

use "$clean/panel_lp.dta", clear

local controls l1_gdpg l2_gdpg ca debt infl imf vix ust10y

* ══════════════════════════════════════════════════════════════════════════
* BLOCK A — FRONTIER VS. ESTABLISHED EM
* ══════════════════════════════════════════════════════════════════════════
/*
  onset_frontier already exists in panel_lp (= d_it_x_frontier from CSV):
    = 1 when onset_all==1 AND country is classified as Frontier Market
  onset_em = residual: onset_all==1 AND frontier==0
*/

capture drop onset_em
gen onset_em = (onset_all == 1 & frontier == 0)
label var onset_em      "Onset dummy — non-frontier EM (established EM)"
label var onset_frontier "Onset dummy — frontier market episodes"

di as result _n "=== BLOCK A: EPISODE COUNTS BY MARKET TYPE ==="
count if onset_frontier == 1 & sample == 1
local n_fr = r(N)
di as result "  Frontier onsets in estimation sample: " `n_fr'
count if onset_em == 1 & sample == 1
local n_em = r(N)
di as result "  Established EM onsets in estimation sample: " `n_em'

* Diagnostic: distribution of frontier flag across all countries
di as result _n "  Frontier flag distribution (all years, estimation sample):"
tab frontier if sample == 1, miss

/*
  Guard: if either group has fewer than 15 onset episodes, the LP subgroup
  is not credibly estimable (too few treated units for DK SE with country FE).
  Skip Block A regressions and issue a warning.
*/
if `n_fr' < 15 | `n_em' < 15 {
    di as error _n "WARNING: Block A skipped."
    di as error "  Frontier onsets = `n_fr', Established EM onsets = `n_em'."
    di as error "  At least one group has fewer than 15 onset episodes."
    di as error "  Check the 'frontier' variable in panel_lp.dta — the split"
    di as error "  may be inverted or the variable may not reflect the intended"
    di as error "  classification. Run: tab country frontier if onset_all==1"
    di as error "  to inspect which countries are coded as frontier."

    * Show which countries are coded as frontier at onset
    di as result _n "  Countries coded frontier==1 at onset:"
    tab country if onset_all == 1 & frontier == 1 & sample == 1
    di as result _n "  Countries coded frontier==0 at onset:"
    tab country if onset_all == 1 & frontier == 0 & sample == 1

    * Create empty IRF placeholders so downstream append does not break
    foreach grp in fr em {
        preserve
            clear
            set obs 5
            gen horizon = _n - 1
            foreach v in b lo90 hi90 lo95 hi95 {
                gen `v' = .
            }
            gen series = "`grp'"
            save "$clean/irf_`grp'.dta", replace
        restore
    }
}
else {

* ── Spec A-1: Frontier markets ───────────────────────────────────────────

di as result _n "--- FRONTIER MARKETS ---"

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_fr = J(5, 1, .)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_frontier `controls' i.year if sample==1, fe lag(`lag')
    matrix b_fr[`h'+1,1]    = _b[onset_frontier]
    matrix lo90_fr[`h'+1,1] = _b[onset_frontier] - 1.645*_se[onset_frontier]
    matrix hi90_fr[`h'+1,1] = _b[onset_frontier] + 1.645*_se[onset_frontier]
    matrix lo95_fr[`h'+1,1] = _b[onset_frontier] - 1.960*_se[onset_frontier]
    matrix hi95_fr[`h'+1,1] = _b[onset_frontier] + 1.960*_se[onset_frontier]
    di "h=" `h' ": beta_frontier=" %6.3f _b[onset_frontier] ///
       "  SE=" %6.3f _se[onset_frontier] "  N=" e(N)
}

* ── Spec A-2: Established EM ─────────────────────────────────────────────

di as result _n "--- ESTABLISHED EM ---"

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_em = J(5, 1, .)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_em `controls' i.year if sample==1, fe lag(`lag')
    matrix b_em[`h'+1,1]    = _b[onset_em]
    matrix lo90_em[`h'+1,1] = _b[onset_em] - 1.645*_se[onset_em]
    matrix hi90_em[`h'+1,1] = _b[onset_em] + 1.645*_se[onset_em]
    matrix lo95_em[`h'+1,1] = _b[onset_em] - 1.960*_se[onset_em]
    matrix hi95_em[`h'+1,1] = _b[onset_em] + 1.960*_se[onset_em]
    di "h=" `h' ": beta_em=" %6.3f _b[onset_em] ///
       "  SE=" %6.3f _se[onset_em] "  N=" e(N)
}

* ── Spec B: Joint regression + equality test ─────────────────────────────

di as result _n "=== JOINT REGRESSION: H0: beta_frontier(h) = beta_em(h) ==="

matrix pval_fr = J(5, 1, .)

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_frontier onset_em `controls' i.year if sample==1, fe lag(`lag')
    test onset_frontier = onset_em
    matrix pval_fr[`h'+1, 1] = r(p)
    di "h=" `h' ":  beta_fr=" %6.3f _b[onset_frontier] ///
               "  beta_em=" %6.3f _b[onset_em] ///
               "  p(diff)=" %5.3f r(p)
}

* ── Save IRF datasets (Block A) ──────────────────────────────────────────

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 lo95 hi95 {
        svmat `m'_fr, names(`m')
        rename `m'1 `m'
    }
    gen series = "frontier"
    save "$clean/irf_frontier.dta", replace
restore

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 lo95 hi95 {
        svmat `m'_em, names(`m')
        rename `m'1 `m'
    }
    gen series = "em"
    save "$clean/irf_em.dta", replace
restore

} // end Block A guard

* ══════════════════════════════════════════════════════════════════════════
* BLOCK B — EPISODE DURATION AS MODERATOR
* ══════════════════════════════════════════════════════════════════════════
/*
  Duration = total years of an episode (onset year + continuation years).
  Computed by counting observations per ep_id (non-missing ep_id covers
  all episode years). Merged back onto the onset row only.

  Threshold: short = 1-2 years, long = 3+ years.
  Rationale: median split; a 1-year spike vs. a multi-year retrenchment
  likely has very different transmission through the real economy.
*/

* ── Step 1: compute episode duration from ep_id ─────────────────────────

preserve
    keep if !missing(ep_id)
    bysort ep_id: gen ep_duration = _N
    keep country year ep_id ep_duration
    keep if !missing(ep_id)
    * Keep only the onset row per episode (smallest year within ep_id)
    bysort ep_id (year): keep if _n == 1
    keep country year ep_duration
    tempfile durations
    save `durations'
restore

merge m:1 country year using `durations', keep(master match) nogen

* ep_duration is only defined at onset year; set missing elsewhere
replace ep_duration = . if onset_all == 0

di as result _n "=== EPISODE DURATION DISTRIBUTION ==="
tab ep_duration if onset_all == 1 & sample == 1

quietly summarize ep_duration if onset_all == 1 & sample == 1
di as result "  Median duration: " r(p50) "  Mean: " %4.1f r(mean)

* ── Step 2: create short/long onset dummies ──────────────────────────────

capture drop onset_short onset_long
gen onset_short = (onset_all == 1 & ep_duration <= 2) & !missing(ep_duration)
gen onset_long  = (onset_all == 1 & ep_duration >= 3) & !missing(ep_duration)
replace onset_short = 0 if missing(ep_duration) & onset_all == 0
replace onset_long  = 0 if missing(ep_duration) & onset_all == 0

label var onset_short "Onset dummy — short episodes (duration <= 2 years)"
label var onset_long  "Onset dummy — long episodes  (duration >= 3 years)"

di as result _n "=== BLOCK B: EPISODE COUNTS BY DURATION ==="
count if onset_short == 1 & sample == 1
di as result "  Short onsets (<=2y): " r(N)
count if onset_long == 1 & sample == 1
di as result "  Long onsets  (>=3y): " r(N)

* ── Spec A-1: Short episodes ─────────────────────────────────────────────

di as result _n "--- SHORT EPISODES (<=2 years) ---"

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_sh = J(5, 1, .)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_short `controls' i.year if sample==1, fe lag(`lag')
    matrix b_sh[`h'+1,1]    = _b[onset_short]
    matrix lo90_sh[`h'+1,1] = _b[onset_short] - 1.645*_se[onset_short]
    matrix hi90_sh[`h'+1,1] = _b[onset_short] + 1.645*_se[onset_short]
    matrix lo95_sh[`h'+1,1] = _b[onset_short] - 1.960*_se[onset_short]
    matrix hi95_sh[`h'+1,1] = _b[onset_short] + 1.960*_se[onset_short]
    di "h=" `h' ": beta_short=" %6.3f _b[onset_short] ///
       "  SE=" %6.3f _se[onset_short] "  N=" e(N)
}

* ── Spec A-2: Long episodes ──────────────────────────────────────────────

di as result _n "--- LONG EPISODES (>=3 years) ---"

foreach m in b lo90 hi90 lo95 hi95 {
    matrix `m'_lg = J(5, 1, .)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_long `controls' i.year if sample==1, fe lag(`lag')
    matrix b_lg[`h'+1,1]    = _b[onset_long]
    matrix lo90_lg[`h'+1,1] = _b[onset_long] - 1.645*_se[onset_long]
    matrix hi90_lg[`h'+1,1] = _b[onset_long] + 1.645*_se[onset_long]
    matrix lo95_lg[`h'+1,1] = _b[onset_long] - 1.960*_se[onset_long]
    matrix hi95_lg[`h'+1,1] = _b[onset_long] + 1.960*_se[onset_long]
    di "h=" `h' ": beta_long=" %6.3f _b[onset_long] ///
       "  SE=" %6.3f _se[onset_long] "  N=" e(N)
}

* ── Spec B: Joint regression + equality test ─────────────────────────────

di as result _n "=== JOINT REGRESSION: H0: beta_short(h) = beta_long(h) ==="

matrix pval_dur = J(5, 1, .)

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    xtscc dy_`h' onset_short onset_long `controls' i.year if sample==1, fe lag(`lag')
    test onset_short = onset_long
    matrix pval_dur[`h'+1, 1] = r(p)
    di "h=" `h' ":  beta_sh=" %6.3f _b[onset_short] ///
               "  beta_lg=" %6.3f _b[onset_long] ///
               "  p(diff)=" %5.3f r(p)
}

* ── Save IRF datasets (Block B) ──────────────────────────────────────────

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 lo95 hi95 {
        svmat `m'_sh, names(`m')
        rename `m'1 `m'
    }
    gen series = "short"
    save "$clean/irf_short.dta", replace
restore

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 lo95 hi95 {
        svmat `m'_lg, names(`m')
        rename `m'1 `m'
    }
    gen series = "long"
    save "$clean/irf_long.dta", replace
restore

* ══════════════════════════════════════════════════════════════════════════
* FIGURES
* ══════════════════════════════════════════════════════════════════════════

local c_fr   "23 55 94"     // navy   — frontier
local c_em   "180 60 40"    // brick  — established EM
local c_sh   "23 55 94"     // navy   — short
local c_lg   "60 140 60"    // green  — long

* ── Figure A: Frontier vs. Established EM ────────────────────────────────

use "$clean/irf_frontier.dta", clear
append using "$clean/irf_em.dta"

twoway ///
    (rarea lo95 hi95 horizon if series=="frontier", ///
        color("`c_fr'%15") lwidth(none)) ///
    (rarea lo95 hi95 horizon if series=="em", ///
        color("`c_em'%15") lwidth(none)) ///
    (rarea lo90 hi90 horizon if series=="frontier", ///
        color("`c_fr'%25") lwidth(none)) ///
    (rarea lo90 hi90 horizon if series=="em", ///
        color("`c_em'%25") lwidth(none)) ///
    (connected b horizon if series=="frontier", ///
        lcolor("`c_fr'") mcolor("`c_fr'") msymbol(circle) ///
        lwidth(medthick) msize(medium)) ///
    (connected b horizon if series=="em", ///
        lcolor("`c_em'") mcolor("`c_em'") msymbol(square) ///
        lwidth(medthick) msize(medium) lpattern(dash)) ///
    , ///
    yline(0, lcolor(gs8) lpattern(dash) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medlarge)) ///
    ylabel(, format(%5.1f) labsize(medlarge)) ///
    xtitle("Years after onset", size(medlarge)) ///
    ytitle("Cumulative % change in real GDP per capita", size(medlarge)) ///
    title("Output Cost by Market Classification", size(large) color(navy)) ///
    subtitle("90% and 95% CI | DK SE | Country & year FE", size(small) color(gs6)) ///
    legend(order(5 "Frontier markets" 6 "Established EM") ///
           ring(1) pos(6) cols(2) size(medsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig10a_frontier_vs_em.pdf", replace
di as result "Figure saved: fig10a_frontier_vs_em.pdf"

* ── Figure B: Short vs. Long Episodes ───────────────────────────────────

use "$clean/irf_short.dta", clear
append using "$clean/irf_long.dta"

twoway ///
    (rarea lo95 hi95 horizon if series=="short", ///
        color("`c_sh'%15") lwidth(none)) ///
    (rarea lo95 hi95 horizon if series=="long", ///
        color("`c_lg'%15") lwidth(none)) ///
    (rarea lo90 hi90 horizon if series=="short", ///
        color("`c_sh'%25") lwidth(none)) ///
    (rarea lo90 hi90 horizon if series=="long", ///
        color("`c_lg'%25") lwidth(none)) ///
    (connected b horizon if series=="short", ///
        lcolor("`c_sh'") mcolor("`c_sh'") msymbol(circle) ///
        lwidth(medthick) msize(medium)) ///
    (connected b horizon if series=="long", ///
        lcolor("`c_lg'") mcolor("`c_lg'") msymbol(square) ///
        lwidth(medthick) msize(medium) lpattern(dash)) ///
    , ///
    yline(0, lcolor(gs8) lpattern(dash) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medlarge)) ///
    ylabel(, format(%5.1f) labsize(medlarge)) ///
    xtitle("Years after onset", size(medlarge)) ///
    ytitle("Cumulative % change in real GDP per capita", size(medlarge)) ///
    title("Output Cost by Episode Duration", size(large) color(navy)) ///
    subtitle("Short: <=2 years | Long: >=3 years | 90% and 95% CI", ///
             size(small) color(gs6)) ///
    legend(order(5 "Short episodes (<=2y)" 6 "Long episodes (>=3y)") ///
           ring(1) pos(6) cols(2) size(medsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig10b_short_vs_long.pdf", replace
di as result "Figure saved: fig10b_short_vs_long.pdf"

* ── Export p-value tables ─────────────────────────────────────────────────

preserve
    clear
    svmat pval_fr
    rename pval_fr1 p_diff
    gen horizon = _n - 1
    order horizon p_diff
    export delimited "$tabs/pval_frontier_vs_em.csv", replace
    di as result "P-values saved: pval_frontier_vs_em.csv"
restore

preserve
    clear
    svmat pval_dur
    rename pval_dur1 p_diff
    gen horizon = _n - 1
    order horizon p_diff
    export delimited "$tabs/pval_short_vs_long.csv", replace
    di as result "P-values saved: pval_short_vs_long.csv"
restore

di as result _n "10_heterogeneity.do complete."
