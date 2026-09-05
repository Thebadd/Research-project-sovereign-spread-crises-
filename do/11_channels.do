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

  ESTIMATOR (matches 02_lp_all.do's own methodology switch): country FE
  only, plain heteroskedasticity-robust SE, no year FE. xtreg ... fe
  vce(robust) is the Stata-idiomatic equivalent of the reference paper's
  reg ..., vce(robust) noconstant with explicit country dummies -- see
  02_lp_all.do's header for the full argument (Figure 3/Table I1's own
  construction script, confirmed directly from Asonuma et al.'s
  replication code).

  IPW REMOVED: this file used to carry an IPW-weighted robustness section
  (Section 7). It was removed project-wide once 08b_aipw.do's doubly-robust
  AIPW estimator superseded plain IPW as the estimator this project reports
  -- see METHODOLOGY.md and 08b_aipw.do's header.

  CONTROLS (uniform common core, Asonuma-aligned):
  -----------------------------------------------
  $ctrl_core = l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2
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
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
*    Same construction as dy_h: cumulative change relative to t-1
* ══════════════════════════════════════════════════════════════════════════

* OUTCOME SCALE. Strictly-positive GDP-ratio channels use the LOG REAL LEVEL
* (ln_r_*, built in 18_transforms) so the outcome is a cumulative percent change in
* the variable itself rather than in its ratio to a GDP that is collapsing — the
* reference paper's var2/var3 construction. pb and fdi change sign, so no log is
* possible and they keep the ratio; for a balance the ratio is the right object.
foreach var in credit claims_govt inv govexp pb fdi {
    local src `var'
    if inlist("`var'","credit","inv","govexp") local src ln_r_`var'
    capture drop `var'_base
    gen `var'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`src' - `var'_base
        label var ch_`var'_`h' "Cum. change `var': F`h' vs t-1"
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    * own-outcome pre-trend control, matching the AIPW spec in 13c.
    capture drop pre_`var'
    gen pre_`var' = L.`src' - L2.`src'
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
  $ctrl_core = l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2
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
local ctrl_credit      l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
local ctrl_claims_govt $ctrl_core pre_claims_govt
local ctrl_inv         $ctrl_core pre_inv
local ctrl_govexp      l1_gdpg l_debt l_banking_crisis l_open l_credit_bank l_lninfl exchange2 pre_govexp
local ctrl_pb          $ctrl_core pre_pb
local ctrl_fdi         $ctrl_core pre_fdi

foreach ch of local channels {
    foreach m in b lo90 hi90 lo95 hi95 {
        matrix `m'_`ch' = J(6, 1, 0)
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 3b. DIAGNOSTIC: FULL REGRESSION TABLE AT h=0 PER CHANNEL
*     Shows coefficient and significance of each control variable.
*     Run without capture so the full xtreg table is displayed.
*     Purpose: verify that controls have the expected signs and are
*     relevant — if a control is systematically insignificant across
*     channels, it may not belong in that specification.
* ══════════════════════════════════════════════════════════════════════════

di as result _n "========================================================"
di as result "DIAGNOSTIC: FULL REGRESSION TABLES AT h=0"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    di as result _n "--- CHANNEL: `ch' ---"
    di as result "Controls: `ctrl'"
    * Guarded: a channel that fails here must NOT halt the whole file
    * (this block is display-only; estimation happens in the loops below).
    capture xtreg ch_`ch'_0 onset_all `ctrl' if sample==1, fe vce(robust)
    if _rc != 0 di as error "  diagnostic xtreg failed for `ch' (rc=" _rc ") — skipped"
}

di as result _n "========================================================"
di as result "END DIAGNOSTIC"
di as result "========================================================"

* ── Channel 1: Credit ────────────────────────────────────────────────────
di as result _n "=== CHANNEL 1: PRIVATE CREDIT / GDP ==="

eststo clear   // capture channel estimates for Table 3 (real loops only, not diagnostic)

forvalues h = 0/4 {
    capture xtreg ch_credit_`h' onset_all `ctrl_credit' ///
        if sample==1, fe vce(robust)
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_credit_`h')
        local nep = r(N)
        eststo t3_credit_`h', title("h=`=`h'+1'")
        estadd scalar nep = `nep'
        local elist_credit `elist_credit' t3_credit_`h'
        matrix b_credit[`h'+2,1]    = _b[onset_all]
        matrix lo90_credit[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_credit[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_credit[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_credit[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h'+1 ": xtreg failed for credit (rc=" _rc ")"
}

* ── Channel 2: Sovereign-bank nexus ──────────────────────────────────────
di as result _n "=== CHANNEL 2: BANK CLAIMS ON GOVT / GDP ==="

forvalues h = 0/4 {
    capture xtreg ch_claims_govt_`h' onset_all `ctrl_claims_govt' ///
        if sample==1, fe vce(robust)
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_claims_govt_`h')
        local nep = r(N)
        eststo t3_claims_govt_`h', title("h=`=`h'+1'")
        estadd scalar nep = `nep'
        local elist_claims_govt `elist_claims_govt' t3_claims_govt_`h'
        matrix b_claims_govt[`h'+2,1]    = _b[onset_all]
        matrix lo90_claims_govt[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_claims_govt[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_claims_govt[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_claims_govt[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h'+1 ": xtreg failed for claims_govt (rc=" _rc ")"
}

* ── Channel 3: Investment ────────────────────────────────────────────────
di as result _n "=== CHANNEL 3: GROSS INVESTMENT / GDP ==="

forvalues h = 0/4 {
    capture xtreg ch_inv_`h' onset_all `ctrl_inv' ///
        if sample==1, fe vce(robust)
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_inv_`h')
        local nep = r(N)
        eststo t3_inv_`h', title("h=`=`h'+1'")
        estadd scalar nep = `nep'
        local elist_inv `elist_inv' t3_inv_`h'
        matrix b_inv[`h'+2,1]    = _b[onset_all]
        matrix lo90_inv[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_inv[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_inv[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_inv[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h'+1 ": xtreg failed for inv (rc=" _rc ")"
}

* ── Channel 4: Government expenditure ────────────────────────────────────
di as result _n "=== CHANNEL 4: GOVERNMENT EXPENDITURE / GDP ==="

forvalues h = 0/4 {
    capture xtreg ch_govexp_`h' onset_all `ctrl_govexp' ///
        if sample==1, fe vce(robust)
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_govexp_`h')
        local nep = r(N)
        eststo t3_govexp_`h', title("h=`=`h'+1'")
        estadd scalar nep = `nep'
        local elist_govexp `elist_govexp' t3_govexp_`h'
        matrix b_govexp[`h'+2,1]    = _b[onset_all]
        matrix lo90_govexp[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_govexp[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_govexp[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_govexp[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h'+1 ": xtreg failed for govexp (rc=" _rc ")"
}

* ── Channel 5: Primary balance ───────────────────────────────────────────
di as result _n "=== CHANNEL 5: PRIMARY BALANCE / GDP ==="

forvalues h = 0/4 {
    capture xtreg ch_pb_`h' onset_all `ctrl_pb' ///
        if sample==1, fe vce(robust)
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_pb_`h')
        local nep = r(N)
        eststo t3_pb_`h', title("h=`=`h'+1'")
        estadd scalar nep = `nep'
        local elist_pb `elist_pb' t3_pb_`h'
        matrix b_pb[`h'+2,1]    = _b[onset_all]
        matrix lo90_pb[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_pb[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_pb[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_pb[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h'+1 ": xtreg failed for pb (rc=" _rc ")"
}

* ── Channel 6: FDI ───────────────────────────────────────────────────────
di as result _n "=== CHANNEL 6: FDI / GDP ==="

forvalues h = 0/4 {
    capture xtreg ch_fdi_`h' onset_all `ctrl_fdi' ///
        if sample==1, fe vce(robust)
    if _rc == 0 {
        quietly count if onset_all==1 & sample==1 & !missing(ch_fdi_`h')
        local nep = r(N)
        eststo t3_fdi_`h', title("h=`=`h'+1'")
        estadd scalar nep = `nep'
        local elist_fdi `elist_fdi' t3_fdi_`h'
        matrix b_fdi[`h'+2,1]    = _b[onset_all]
        matrix lo90_fdi[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_fdi[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
        matrix lo95_fdi[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
        matrix hi95_fdi[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
        di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
           "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
    }
    else di as error "h=" `h'+1 ": xtreg failed for fdi (rc=" _rc ")"
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — TABLE 3: Transmission channels (pooled, all episodes)
*   Word/RTF, multi-panel: one panel per channel, columns = horizons h=0..4.
*   Each panel reports the onset coefficient (robust SE in parentheses) on the
*   cumulative change in the channel variable. First panel uses replace;
*   remaining panels append to the same file.
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

local t3note "Dependent variable: cumulative change in the channel variable (pp) from t-1 to t+h. Jorda (2005) local projections; country fixed effects only (no year FE); common-core controls plus the channel's own pre-crisis change; continuation years excluded. Robust (heteroskedasticity-only) standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01."

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
        set obs 6
        gen horizon = _n - 1     // 0 (baseline), 1..5
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
            replace horizon = `h'+1               in `row'
            replace b       = b_`ch'[`h'+2, 1]   in `row'
            replace lo95    = lo95_`ch'[`h'+2, 1] in `row'
            replace hi95    = hi95_`ch'[`h'+2, 1] in `row'
            replace lo90    = lo90_`ch'[`h'+2, 1] in `row'
            replace hi90    = hi90_`ch'[`h'+2, 1] in `row'
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
        xlabel(0(1)5, labsize(medsmall)) ///
        ylabel(, format(%5.2f) labsize(medsmall)) ///
        xtitle("Year (Year 1 = crisis year)", size(small)) ///
        ytitle("Cumulative change (% or pp — see note)", size(small)) ///
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
    note("90% and 95% CI. Robust (heteroskedasticity-only) SE. Country FE only (no year FE). Treatment: all 61 onset episodes." ///
         "Units differ by channel: private credit, bank claims on govt, investment and govt expenditure are LOG REAL LEVELS, so their scale is cumulative percent change (comparable to the GDP result). Primary balance and FDI change sign, so no log is possible and they remain ratios to GDP, in percentage points.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig11_channels.pdf", replace
di as result "Figure saved: fig11_channels.pdf"

* Close individual graphs from memory
foreach nm of local fignames {
    graph drop `nm'
}


di as result _n "11_channels.do complete."
