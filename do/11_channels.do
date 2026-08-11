/*===========================================================================
  11_CHANNELS.DO
  Transmission Channels of Sovereign Spread Crises

  Research question: through which channels does a spread crisis affect
  real output? We run the same LP framework as Act 1 but replace the
  GDP outcome with six channel variables, all sharing the UNIFORM common-core
  control set ($ctrl_core) plus each channel's own pre-crisis change (pre_<v>).

  CHANNELS AND OUTCOMES:
  ----------------------
  1. Credit channel          → credit_gdp  (private credit / GDP)
  2. Sovereign-bank nexus    → claims_govt (bank claims on govt / GDP)
  3. Investment channel      → inv         (gross investment / GDP)
  4. Fiscal channel          → govexp      (government expenditure / GDP)
  5. Fiscal adjustment       → pb          (primary balance / GDP)
  6. External financing      → fdi         (FDI / GDP)

  SPECIFICATION (horizon h = 0,...,4):
  -------------------------------------
  Channel_var(i,t+h) - Channel_var(i,t-1) = αi + γt + β(h)·onset_all
                                            + X_core (common core + pre_<v>) · δ + ε(i,t+h)

  All outcomes expressed as cumulative change from t-1 (same anchoring
  as main LP). Controls are pre-determined (lagged relative to t).
  DK SE with lag = max(1, h+1) throughout.

  CONTROLS (uniform common core, Asonuma-aligned):
  -----------------------------------------------
  $ctrl_core = l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank
               l_hyperinfl      (+ each channel's own pre_<v>)
  The core term equal to a channel's own lagged level is dropped from its own
  regression (credit -> the depth term; govexp -> l_govexp). Global time-series
  factors (fed funds / VIX / UST10Y) are omitted — absorbed by year FE.

  Saves:
    "$clean/irf_ch_*.dta"        — one IRF dataset per channel
    "$figs/fig11_channels.pdf"   — 2×3 multi-panel figure
    "$tabs/channels_summary.csv" — coefficients and SE for all channels
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
*    Same construction as dy_h: cumulative change relative to t-1
* ══════════════════════════════════════════════════════════════════════════

foreach var in credit claims_govt inv govexp pb fdi {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
        label var ch_`var'_`h' "Cum. change `var': F`h' vs t-1"
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    * own-outcome pre-trend control, matching the AIPW spec in 13c.
    capture drop pre_`var'
    gen pre_`var' = L.`var' - L2.`var'
}

* ══════════════════════════════════════════════════════════════════════════
* 2. COVERAGE CHECK
*    Report how many onset observations have non-missing channel data
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== DATA COVERAGE AT ONSET (sample==1 & onset_all==1) ==="
foreach var in credit claims_govt inv govexp pb fdi {
    quietly count if onset_all == 1 & sample == 1 & !missing(ch_`var'_0)
    di as result "  `var': " r(N) " / 61 onsets with non-missing data at h=0"
}

* ══════════════════════════════════════════════════════════════════════════
* 3. LP ESTIMATION BY CHANNEL
* ══════════════════════════════════════════════════════════════════════════

/*
  CONTROL SET (uniform across all six channels):
  ----------------------------------------------
  $ctrl_core = l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank
               l_hyperinfl   + each channel's own pre_<v> (Asonuma g_0).
  Rationale: pre-crisis GDP momentum, fiscal solvency (debt) and external
  balance (ca), banking distress, government spending, trade openness, private-
  credit depth, and a hyperinflation flag — the Asonuma $convar adapted to
  spread crises. The core term equal to a channel's own lagged level is dropped
  from that channel's regression (credit -> the depth term; govexp -> l_govexp).
  Global push factors (fed funds / VIX / UST10Y) enter the first-stage
  propensity model only; in the LP they are absorbed by year FE.
*/

local channels credit claims_govt inv govexp pb fdi

* Controls: uniform common core ($ctrl_core) + each channel's own pre_<v>; drop the core term equal to the channel's own lagged level.
* Common-core controls (Asonuma-aligned $ctrl_core) + each channel's own pre_<v>;
* the core term measuring the channel's own lagged level is dropped from its own reg.
* For the credit channel that term is l_credit_bank: not literally the lag of the
* outcome series (credit = all financial corps, l_credit_bank = banks only), but the
* two correlate 0.950, so keeping it would be close to putting the lagged dependent
* variable on the RHS. It stays dropped.
local ctrl_credit      l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_hyperinfl pre_credit
local ctrl_claims_govt $ctrl_core pre_claims_govt
local ctrl_inv         $ctrl_core pre_inv
local ctrl_govexp      l1_gdpg l_debt l_ca l_banking_duration l_open l_credit_bank l_hyperinfl pre_govexp
local ctrl_pb          $ctrl_core pre_pb
local ctrl_fdi         $ctrl_core pre_fdi

foreach ch of local channels {
    foreach m in b lo90 hi90 lo95 hi95 {
        matrix `m'_`ch' = J(5, 1, .)
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 3b. DIAGNOSTIC: FULL REGRESSION TABLE AT h=0 PER CHANNEL
*     Shows coefficient and significance of each control variable.
*     Run without capture so the full xtscc table is displayed.
*     Purpose: verify that controls have the expected signs and are
*     relevant — if a control is systematically insignificant across
*     channels, it may not belong in that specification.
* ══════════════════════════════════════════════════════════════════════════

di as result _n "========================================================"
di as result "DIAGNOSTIC: FULL REGRESSION TABLES AT h=0 (lag=1)"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    di as result _n "--- CHANNEL: `ch' ---"
    di as result "Controls: `ctrl'"
    * Guarded: a channel that fails here must NOT halt the whole file
    * (this block is display-only; estimation happens in the loops below).
    capture xtscc ch_`ch'_0 onset_all `ctrl' i.year if sample==1, fe lag(1)
    if _rc != 0 di as error "  diagnostic xtscc failed for `ch' (rc=" _rc ") — skipped"
}

di as result _n "========================================================"
di as result "END DIAGNOSTIC"
di as result "========================================================"

* ── Channel 1: Credit ────────────────────────────────────────────────────
di as result _n "=== CHANNEL 1: PRIVATE CREDIT / GDP ==="

eststo clear   // capture channel estimates for Table 3 (real loops only, not diagnostic)

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_credit_`h' onset_all `ctrl_credit' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_credit_`h')
        local nep = r(N)
        eststo t3_credit_`h', title("h=`h'")
        estadd scalar nep = `nep'
        local elist_credit `elist_credit' t3_credit_`h'
        matrix b_credit[`h'+1,1]    = _b[onset_all]
        matrix lo90_credit[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_credit[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_credit[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_credit[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for credit (rc=" _rc ")"
}

* ── Channel 2: Sovereign-bank nexus ──────────────────────────────────────
di as result _n "=== CHANNEL 2: BANK CLAIMS ON GOVT / GDP ==="

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_claims_govt_`h' onset_all `ctrl_claims_govt' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_claims_govt_`h')
        local nep = r(N)
        eststo t3_claims_govt_`h', title("h=`h'")
        estadd scalar nep = `nep'
        local elist_claims_govt `elist_claims_govt' t3_claims_govt_`h'
        matrix b_claims_govt[`h'+1,1]    = _b[onset_all]
        matrix lo90_claims_govt[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_claims_govt[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_claims_govt[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_claims_govt[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for claims_govt (rc=" _rc ")"
}

* ── Channel 3: Investment ────────────────────────────────────────────────
di as result _n "=== CHANNEL 3: GROSS INVESTMENT / GDP ==="

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_inv_`h' onset_all `ctrl_inv' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_inv_`h')
        local nep = r(N)
        eststo t3_inv_`h', title("h=`h'")
        estadd scalar nep = `nep'
        local elist_inv `elist_inv' t3_inv_`h'
        matrix b_inv[`h'+1,1]    = _b[onset_all]
        matrix lo90_inv[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_inv[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_inv[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_inv[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for inv (rc=" _rc ")"
}

* ── Channel 4: Government expenditure ────────────────────────────────────
di as result _n "=== CHANNEL 4: GOVERNMENT EXPENDITURE / GDP ==="

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_govexp_`h' onset_all `ctrl_govexp' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_govexp_`h')
        local nep = r(N)
        eststo t3_govexp_`h', title("h=`h'")
        estadd scalar nep = `nep'
        local elist_govexp `elist_govexp' t3_govexp_`h'
        matrix b_govexp[`h'+1,1]    = _b[onset_all]
        matrix lo90_govexp[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_govexp[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_govexp[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_govexp[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for govexp (rc=" _rc ")"
}

* ── Channel 5: Primary balance ───────────────────────────────────────────
di as result _n "=== CHANNEL 5: PRIMARY BALANCE / GDP ==="

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_pb_`h' onset_all `ctrl_pb' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_pb_`h')
        local nep = r(N)
        eststo t3_pb_`h', title("h=`h'")
        estadd scalar nep = `nep'
        local elist_pb `elist_pb' t3_pb_`h'
        matrix b_pb[`h'+1,1]    = _b[onset_all]
        matrix lo90_pb[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_pb[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_pb[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_pb[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for pb (rc=" _rc ")"
}

* ── Channel 6: FDI ───────────────────────────────────────────────────────
di as result _n "=== CHANNEL 6: FDI / GDP ==="

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    capture xtscc ch_fdi_`h' onset_all `ctrl_fdi' i.year ///
        if sample==1, fe lag(`lag')
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_fdi_`h')
        local nep = r(N)
        eststo t3_fdi_`h', title("h=`h'")
        estadd scalar nep = `nep'
        local elist_fdi `elist_fdi' t3_fdi_`h'
        matrix b_fdi[`h'+1,1]    = _b[onset_all]
        matrix lo90_fdi[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_fdi[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_fdi[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_fdi[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h' ": xtscc failed for fdi (rc=" _rc ")"
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — TABLE 3: Transmission channels (pooled, all episodes)
*   Word/RTF, multi-panel: one panel per channel, columns = horizons h=0..4.
*   Each panel reports the onset coefficient (DK SE in parentheses) on the
*   cumulative change in the channel variable. First panel uses replace;
*   remaining panels append to the same file.
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

local t3note "Dependent variable: cumulative change in the channel variable (pp) from t-1 to t+h. Jorda (2005) local projections; country and year fixed effects; common-core controls plus the channel's own pre-crisis change; continuation years excluded. Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01."

* Per-channel panel titles (Panel A carries the overall table caption)
local ptitle_credit      "Table 3. Transmission channels (all episodes) -- Panel A: Private credit/GDP"
local ptitle_claims_govt "Panel B: Bank claims on govt/GDP"
local ptitle_inv         "Panel C: Investment/GDP"
local ptitle_govexp      "Panel D: Govt expenditure/GDP"
local ptitle_pb          "Panel E: Primary balance/GDP"
local ptitle_fdi         "Panel F: FDI/GDP"

* Write one panel per channel to a single RTF. First successful panel uses
* "replace" (creates the file); the rest "append". Each esttab is wrapped in
* capture so a locked file (open in Word) or a missing estimate warns and is
* skipped instead of halting the whole do-file.
local writemode replace
local t3fail 0

foreach ch in credit claims_govt inv govexp pb fdi {

    * Skip channel entirely if no horizon estimate was stored
    if "`elist_`ch''" == "" {
        di as error "  ** Table 3: no estimates for channel `ch' — panel skipped"
        local t3fail 1
        continue
    }

    * Attach the methodology note to the last channel only
    local t3extra
    if "`ch'" == "fdi" local t3extra addnotes("`t3note'")

    capture esttab `elist_`ch'' using "$tabs/table3_channels.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_all) coeflabel(onset_all "Spread-crisis onset") ///
        mtitles nonumber ///
        stats(nep N N_g, labels("Episodes (onsets)" "Observations" "Countries") fmt(0 0 0)) ///
        title("`ptitle_`ch''") `t3extra'

    if _rc == 608 {
        di as error "  ** table3_channels.rtf is OPEN IN WORD — close it and re-run to refresh."
        local t3fail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 3: esttab failed for panel `ch' (rc=" _rc ")"
        local t3fail 1
        continue
    }

    * After the first successful write, switch to append mode
    local writemode append
}

if `t3fail' == 0 di as result "Table 3 saved: $tabs/table3_channels.rtf"
else di as error "Table 3 written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* 4. SAVE IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

local channels      credit claims_govt inv govexp pb fdi
local chan_labels   "Private credit/GDP" "Bank claims on govt/GDP" ///
                    "Investment/GDP" "Govt expenditure/GDP" ///
                    "Primary balance/GDP" "FDI/GDP"

local i = 1
foreach ch of local channels {
    preserve
        clear
        set obs 5
        gen horizon = _n - 1
        foreach m in b lo90 hi90 lo95 hi95 {
            svmat `m'_`ch', names(`m')
            rename `m'1 `m'
        }
        gen channel = "`ch'"
        save "$clean/irf_ch_`ch'.dta", replace
    restore
    local ++i
}

* ══════════════════════════════════════════════════════════════════════════
* 5. EXPORT SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

preserve
    clear
    local nobs = 5 * 6    // 5 horizons × 6 channels
    set obs `nobs'
    gen channel = ""
    gen horizon = .
    gen b       = .
    gen lo95    = .
    gen hi95    = .
    gen lo90    = .
    gen hi90    = .

    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel = "`ch'"              in `row'
            replace horizon = `h'                 in `row'
            replace b       = b_`ch'[`h'+1, 1]   in `row'
            replace lo95    = lo95_`ch'[`h'+1, 1] in `row'
            replace hi95    = hi95_`ch'[`h'+1, 1] in `row'
            replace lo90    = lo90_`ch'[`h'+1, 1] in `row'
            replace hi90    = hi90_`ch'[`h'+1, 1] in `row'
            local ++row
        }
    }

    order channel horizon b lo90 hi90 lo95 hi95
    export delimited "$tabs/channels_summary.csv", replace
    di as result "Channel summary saved: $tabs/channels_summary.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 6. FIGURE — 2×3 MULTI-PANEL IRF
* ══════════════════════════════════════════════════════════════════════════

local c_main "23 55 94"

local channels    credit claims_govt inv govexp pb fdi
local titlelabels `" "Private Credit/GDP" "Bank Claims on Govt/GDP" "Investment/GDP" "Govt Expenditure/GDP" "Primary Balance/GDP" "FDI/GDP" "'
local fignames    fig11a fig11b fig11c fig11d fig11e fig11f

local i = 1
foreach ch of local channels {

    local tlab : word `i' of `titlelabels'

    use "$clean/irf_ch_`ch'.dta", clear

    twoway ///
        (rarea lo95 hi95 horizon, ///
            color("`c_main'%15") lwidth(none)) ///
        (rarea lo90 hi90 horizon, ///
            color("`c_main'%28") lwidth(none)) ///
        (connected b horizon, ///
            lcolor("`c_main'") mcolor("`c_main'") ///
            msymbol(circle) lwidth(medthick) msize(medium)) ///
        , ///
        yline(0, lcolor(gs8) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)4, labsize(medsmall)) ///
        ylabel(, format(%5.2f) labsize(medsmall)) ///
        xtitle("Years after onset", size(small)) ///
        ytitle("Cumulative change (pp)", size(small)) ///
        title(`tlab', size(medsmall) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(`: word `i' of `fignames'', replace)

    local ++i
}

* Combine into 2×3 grid
graph combine fig11a fig11b fig11c fig11d fig11e fig11f, ///
    cols(3) rows(2) ///
    title("Transmission Channels of Sovereign Spread Crises", ///
          size(medlarge) color(navy)) ///
    note("90% and 95% CI. DK SE. Country & year FE. Treatment: all 61 onset episodes.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig11_channels.pdf", replace
di as result "Figure saved: fig11_channels.pdf"

* Close individual graphs from memory
foreach nm of local fignames {
    graph drop `nm'
}

* ══════════════════════════════════════════════════════════════════════════
* 7. IPW-WEIGHTED CHANNEL REGRESSIONS
*    Replicates section 3 with stabilized IPW weights to correct for
*    selection: onset country-years are not randomly drawn from the panel —
*    they have systematically worse pre-crisis fundamentals. IPW reweights
*    tranquil observations to look like onset observations on observables,
*    making the comparison credible.
*
*    Method: areg absorb(cid) [aw=ipw] vce(cluster cid)
*    Note:   xtscc does not accept pweights; areg used here (same as 08_ipw_lp.do).
*            SE are clustered by country rather than Driscoll-Kraay.
*
*    Output: side-by-side comparison of unweighted vs. IPW beta at each horizon.
* ══════════════════════════════════════════════════════════════════════════

* ── 7a. Reload panel and construct IPW weights ───────────────────────────

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* Regenerate channel outcome variables + own pre-trend (same as section 1;
* the reload above dropped them, and ctrl_<ch> now includes pre_<var>)
foreach var in credit claims_govt inv govexp pb fdi {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
    }
    capture drop pre_`var'
    gen pre_`var' = L.`var' - L2.`var'
}

* First-stage probit: same spec as 08_ipw_lp.do
di as result _n "=== IPW FIRST STAGE (channels): Probit of crisis onset ==="
* Controls + excluded predictors (Z), consistent with 08_ipw_lp.do:
* fedfunds (US monetary conditions) + l_reg_crisis_share/past_onsets as excluded predictors (Z).
* Baseline = outcome baseline ($ctrl_core) + excluded predictors Z (strict parity).
local controls_ps $ctrl_core l_fedfunds l_reg_crisis_share past_onsets
* Pooled probit (country FE separate on the thin event count; FE in the outcome only).
probit onset_all `controls_ps' if sample == 1, vce(cluster cid)
di as result "McFadden Pseudo-R2: " e(r2_p)

predict pscore_ch if sample == 1, pr

* Trim [0.01, 0.99]
gen trimmed_ch = (pscore_ch < 0.01 | pscore_ch > 0.99) if !missing(pscore_ch)
quietly count if trimmed_ch == 1
di as result "Trimmed observations: " r(N)
replace pscore_ch = . if trimmed_ch == 1

* Stabilized IPW weights
quietly summarize onset_all if sample == 1
local p_treat = r(mean)
local p_ctrl  = 1 - `p_treat'
di as result "Marginal Pr(onset=1): " %5.4f `p_treat'

gen ipw_ch = .
replace ipw_ch = `p_treat' / pscore_ch        if onset_all == 1 & !missing(pscore_ch)
replace ipw_ch = `p_ctrl'  / (1 - pscore_ch)  if onset_all == 0 & !missing(pscore_ch)
label var ipw_ch "Stabilized IPW weight (channels)"

summarize ipw_ch if sample == 1, detail

* ── 7b. Re-declare control locals (same as section 3) ────────────────────

* Common-core controls (Asonuma-aligned $ctrl_core) + each channel's own pre_<v>;
* the core term measuring the channel's own lagged level is dropped from its own reg.
* For the credit channel that term is l_credit_bank: not literally the lag of the
* outcome series (credit = all financial corps, l_credit_bank = banks only), but the
* two correlate 0.950, so keeping it would be close to putting the lagged dependent
* variable on the RHS. It stays dropped.
local ctrl_credit      l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_hyperinfl pre_credit
local ctrl_claims_govt $ctrl_core pre_claims_govt
local ctrl_inv         $ctrl_core pre_inv
local ctrl_govexp      l1_gdpg l_debt l_ca l_banking_duration l_open l_credit_bank l_hyperinfl pre_govexp
local ctrl_pb          $ctrl_core pre_pb
local ctrl_fdi         $ctrl_core pre_fdi

* ── 7c. IPW-weighted LP: store results ───────────────────────────────────

local channels credit claims_govt inv govexp pb fdi

foreach ch of local channels {
    foreach m in b_ipw lo90_ipw hi90_ipw lo95_ipw hi95_ipw {
        matrix `m'_`ch' = J(5, 1, .)
    }
}

* ── 7d. Run IPW-weighted LP and print side-by-side comparison ─────────────

di as result _n "========================================================"
di as result "IPW-WEIGHTED vs. UNWEIGHTED: CHANNEL COEFFICIENTS"
di as result "  (areg absorb(cid) [aw=ipw] vce(cluster cid))"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''

    di as result _n "--- CHANNEL: `ch' ---"
    di as result "h    beta_OLS   SE_OLS    beta_IPW   SE_IPW    delta"

    forvalues h = 0/4 {
        local row = `h' + 1

        * Unweighted baseline (areg, for apples-to-apples SE comparison)
        * Country FE only, no year FE (paper-aligned two-stage outcome regression).
        * Guarded like the IPW fit below: one failing channel must warn and be
        * skipped, not halt the whole file (a missing control used to throw r(111)
        * here and kill the run mid-loop, after the pooled figure had been written).
        capture areg ch_`ch'_`h' onset_all `ctrl' ///
            if sample == 1 & !missing(ipw_ch), absorb(cid) vce(cluster cid)
        if _rc == 0 {
            local b_ols  = _b[onset_all]
            local se_ols = _se[onset_all]
        }
        else {
            local b_ols  = .
            local se_ols = .
            di as error "h=" `h' ": areg OLS baseline failed for `ch' (rc=" _rc ")"
        }

        * IPW-weighted
        capture areg ch_`ch'_`h' onset_all `ctrl' ///
            [aw=ipw_ch] if sample == 1 & !missing(ipw_ch), absorb(cid) vce(cluster cid)
        if _rc == 0 {
            matrix b_ipw_`ch'[`row',1]    = _b[onset_all]
            matrix lo90_ipw_`ch'[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
            matrix hi90_ipw_`ch'[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
            matrix lo95_ipw_`ch'[`row',1] = _b[onset_all] - 1.960*_se[onset_all]
            matrix hi95_ipw_`ch'[`row',1] = _b[onset_all] + 1.960*_se[onset_all]
            local b_ipw  = _b[onset_all]
            local se_ipw = _se[onset_all]
            di "h=" `h' "   " %7.3f `b_ols'  "   " %6.3f `se_ols' ///
                   "    " %7.3f `b_ipw' "   " %6.3f `se_ipw' ///
                   "    " %7.3f (`b_ipw' - `b_ols')
        }
        else {
            di as error "h=" `h' ": areg IPW failed for `ch' (rc=" _rc ")"
        }
    }
}

di as result _n "========================================================"
di as result "END IPW COMPARISON"
di as result "========================================================"

* ── 7e. Save IPW IRF datasets ─────────────────────────────────────────────

foreach ch of local channels {
    preserve
        clear
        set obs 5
        gen horizon = _n - 1
        foreach m in b lo90 hi90 lo95 hi95 {
            svmat `m'_ipw_`ch', names(`m')
            rename `m'1 `m'
        }
        gen channel = "`ch'"
        save "$clean/irf_ch_`ch'_ipw.dta", replace
    restore
}

* ── 7f. Figure — IPW vs. unweighted overlay per channel ──────────────────

local c_ols "23 55 94"
local c_ipw "157 36 73"

local channels    credit claims_govt inv govexp pb fdi
local titlelabels `" "Private Credit/GDP" "Bank Claims on Govt/GDP" "Investment/GDP" "Govt Expenditure/GDP" "Primary Balance/GDP" "FDI/GDP" "'
local fignames_cmp fig11a_cmp fig11b_cmp fig11c_cmp fig11d_cmp fig11e_cmp fig11f_cmp

local i = 1
foreach ch of local channels {

    local tlab : word `i' of `titlelabels'

    * Load unweighted IRF
    use "$clean/irf_ch_`ch'.dta", clear
    gen series = "ols"
    tempfile ols_`ch'
    save `ols_`ch''

    * Load IPW IRF
    use "$clean/irf_ch_`ch'_ipw.dta", clear
    gen series = "ipw"
    append using `ols_`ch''

    twoway ///
        (rarea lo90 hi90 horizon if series=="ols", ///
            color("`c_ols'%20") lwidth(none)) ///
        (connected b horizon if series=="ols", ///
            lcolor("`c_ols'") lwidth(medthick) msymbol(circle) mcolor("`c_ols'")) ///
        (rarea lo90 hi90 horizon if series=="ipw", ///
            color("`c_ipw'%20") lwidth(none)) ///
        (connected b horizon if series=="ipw", ///
            lcolor("`c_ipw'") lwidth(medthick) lpattern(dash) ///
            msymbol(square) mcolor("`c_ipw'")), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xlabel(0(1)4, labsize(medsmall)) ///
        ylabel(, format(%5.2f) labsize(medsmall)) ///
        xtitle("Years after onset", size(small)) ///
        ytitle("Cumulative change (pp)", size(small)) ///
        title(`tlab', size(medsmall) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        nodraw name(`: word `i' of `fignames_cmp'', replace)

    local ++i
}

graph combine fig11a_cmp fig11b_cmp fig11c_cmp fig11d_cmp fig11e_cmp fig11f_cmp, ///
    cols(3) rows(2) ///
    title("Channels: Unweighted vs. IPW-Weighted", size(medlarge) color(navy)) ///
    note("Blue solid = unweighted (DK SE). Red dashed = IPW-weighted (cluster SE)." ///
         "IPW weights from probit of onset on lagged macro fundamentals.", size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig11_channels_ipw_compare.pdf", replace
di as result "Figure saved: fig11_channels_ipw_compare.pdf"

foreach nm of local fignames_cmp {
    graph drop `nm'
}

di as result _n "11_channels.do complete."
