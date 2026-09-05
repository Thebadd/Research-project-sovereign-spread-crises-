/*===========================================================================
  12_CHANNELS_RESOLUTION.DO
  Transmission Channels by Resolution Type (Non-Default vs. Default-Linked)

  Research question: do the transmission channels of spread crises differ
  between episodes resolved without default and those linked to default?

  STRATEGY:
  ---------
  For each of the 6 channels from 11_channels.do: JOINT LP with onset_nd and
  onset_def entered simultaneously on the FULL sample, tranquil the omitted
  category (the reference paper's baseline):
    ch_var(h) = αi + β_nd(h)·onset_nd + β_def(h)·onset_def
               + X_core (common core + pre_<v>)·δ + ε
  Country FE only, no year FE, plain robust SE -- the Stata-idiomatic
  equivalent of the reference paper's own reg ..., vce(robust) noconstant
  with explicit country dummies (see 02_lp_all.do's header). Since both
  type dummies are estimated in ONE regression, their covariance is
  exactly estimable and the difference test below uses the Wald
  F-statistic (test onset_nd = onset_def) directly -- the reference
  paper's own difference-test convention (Table I1) -- not Clogg z or a
  bootstrap, which would only approximate what the F-test already gives
  exactly.

  IPW REMOVED: this file used to carry a parallel IPW-weighted (Spec B)
  comparison, dropped project-wide once 08b_aipw.do's doubly-robust AIPW
  estimator superseded plain IPW as the estimator this project reports --
  see METHODOLOGY.md and 08b_aipw.do's header.

  OUTPUTS:
  --------
  - fig12a_channels_ols.pdf   : 2×3 grid (nd=navy, def=brick)
  - channels_resolution.csv   : full coefficient table
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
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
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    * own-outcome pre-trend control, matching the AIPW spec in 13c.
    capture drop pre_`var'
    gen pre_`var' = L.`src' - L2.`src'
}

* ══════════════════════════════════════════════════════════════════════════
* 2. COVERAGE CHECK BY RESOLUTION TYPE
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== DATA COVERAGE AT ONSET BY RESOLUTION TYPE ==="
di as result "  Variable        nd (39)   def (22)"
foreach var in credit claims_govt inv govexp pb fdi {
    quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`var'_0)
    local n_nd = r(N)
    quietly count if onset_def == 1 & sample == 1 & !missing(ch_`var'_0)
    local n_def = r(N)
    di as result "  `var'" _col(20) `n_nd' " / 39" _col(32) `n_def' " / 22"
}

* ══════════════════════════════════════════════════════════════════════════
* 3. LP ESTIMATION BY CHANNEL
* ══════════════════════════════════════════════════════════════════════════

local channels   credit claims_govt inv govexp pb fdi

* Controls: common core ($ctrl_core) + each channel's own pre_<v> (same as 11_channels.do)
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

* Initialize storage matrices
foreach ch of local channels {
    foreach grp in nd def {
        foreach m in b lo90 hi90 lo95 hi95 {
            matrix `m'_`grp'_`ch' = J(6, 1, 0)
        }
    }
    matrix pval_`ch' = J(6, 1, .)
}

* ── Loop over channels ───────────────────────────────────────────────────

eststo clear   // capture OLS resolution-split estimates for Table 4

foreach ch of local channels {

    local ctrl `ctrl_`ch''

    di as result _n "========================================"
    di as result "CHANNEL: `ch'"
    di as result "========================================"
    di "h   b_nd     b_def    p(nd=def)"

    forvalues h = 0/4 {
        local row = `h' + 2

        * JOINT LP, both type dummies, FULL sample. The reference paper's
        * baseline is a single joint regression with all type dummies,
        * country dummies, vce(robust), no year FE, and tranquil as the
        * omitted category (reg g_h dum1 dum2 dum3 g_0 $convar, vce(robust),
        * noconstant).
        capture xtreg ch_`ch'_`h' onset_nd onset_def `ctrl' ///
            if sample == 1, fe vce(robust)

        if _rc == 0 {
            matrix b_nd_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix lo95_nd_`ch'[`row',1] = _b[onset_nd]  - 1.960*_se[onset_nd]
            matrix hi95_nd_`ch'[`row',1] = _b[onset_nd]  + 1.960*_se[onset_nd]
            matrix b_def_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix lo95_def_`ch'[`row',1]= _b[onset_def] - 1.960*_se[onset_def]
            matrix hi95_def_`ch'[`row',1]= _b[onset_def] + 1.960*_se[onset_def]
            test onset_nd = onset_def
            matrix pval_`ch'[`row',1] = r(p)
            local pd_`ch'_`h' = r(p)   // store before eststo (which can reset r())
            local pf_`ch'_`h' = r(F)   // F-statistic itself, Table I1's own convention

            * Difference block: point estimate + the exact, covariance-correct
            * Wald F-test (not Clogg z/bootstrap -- see header) + episode counts.
            local bnd  = _b[onset_nd]
            local bdef = _b[onset_def]
            local bdiff = `bdef' - `bnd'
            quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepnd = r(N)
            quietly count if onset_def == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepdef = r(N)

            eststo t4_`ch'_`h', title("h=`=`h'+1'")
            estadd scalar bdiff  = `bdiff'
            estadd scalar fdiff  = `pf_`ch'_`h''
            estadd scalar pdiff  = `pd_`ch'_`h''
            estadd scalar nepnd  = `nepnd'
            estadd scalar nepdef = `nepdef'
            local elist_`ch' `elist_`ch'' t4_`ch'_`h'
            local b_nd_o  = `bnd'
            local b_def_o = `bdef'
            * NOT r(p): eststo/estadd above clear r(), so r(p) is empty by this
            * point and the console column printed "." for every row.
            local p_o     = `pd_`ch'_`h''
        }
        else {
            local b_nd_o  = .
            local b_def_o = .
            local p_o     = .
            di as error "regression failed for `ch' h=`=`h'+1'"
        }

        di "h=" `h'+1 "  " %7.3f `b_nd_o'  "  " %7.3f `b_def_o' "  " %5.3f `p_o'
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 4. TABLE EXPORT — TABLE 4: Transmission channels by resolution type
*   Word/RTF, multi-panel: one panel per channel, columns = horizons h=0..4.
*   Each panel reports non-default and default-linked onset coefficients from the
*   JOINT regression (the reference paper's baseline), robust SE in parentheses,
*   plus the difference and the Wald F-test of their equality (test onset_nd =
*   onset_def) -- the covariance-correct answer for two coefficients estimated
*   in one regression, and the reference paper's own difference-test convention
*   (Table I1), not Clogg z or a bootstrap.
*   OLS spec (matches Tables 1-3). First panel replaces; the rest append.
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

local t4note "Dependent variable: cumulative change in the channel variable (pp) from t-1 to t+h. Both onset dummies enter jointly with tranquil years as the omitted category, matching the reference paper's baseline. Jorda (2005) local projections; country fixed effects only (no year FE); common-core controls plus the channel's own pre-crisis change; continuation years excluded. Robust (heteroskedasticity-only) standard errors in parentheses. F(nd=def) and p(nd=def) are the Wald equality test. * p<0.10, ** p<0.05, *** p<0.01."

* Per-channel panel titles (Panel A carries the overall table caption)
local ptitle_credit      "Table 4. Channels by resolution (nd vs. def, joint) -- Panel A: Private credit/GDP"
local ptitle_claims_govt "Panel B: Bank claims on govt/GDP"
local ptitle_inv         "Panel C: Investment/GDP"
local ptitle_govexp      "Panel D: Govt expenditure/GDP"
local ptitle_pb          "Panel E: Primary balance/GDP"
local ptitle_fdi         "Panel F: FDI/GDP"

* Write one panel per channel to a single RTF. First successful panel uses
* "replace"; the rest "append". Each esttab is wrapped in capture so a locked
* file (open in Word) or a missing estimate warns and is skipped instead of
* halting the do-file.
local writemode replace
local t4fail 0

foreach ch in credit claims_govt inv govexp pb fdi {

    if "`elist_`ch''" == "" {
        di as error "  ** Table 4: no estimates for channel `ch' — panel skipped"
        local t4fail 1
        continue
    }

    local t4extra
    if "`ch'" == "fdi" local t4extra addnotes("`t4note'")

    capture esttab `elist_`ch'' using "$tabs/table4_channels_resolution.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_nd onset_def) order(onset_nd onset_def) ///
        coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
        mtitles nonumber ///
        stats(bdiff fdiff pdiff nepnd nepdef N N_g, ///
              labels("Difference (default - non-default)" "  F (Wald, nd = def)" ///
                     "  p (Wald, nd = def)" ///
                     "Episodes (non-default)" "Episodes (default)" ///
                     "Observations" "Countries") ///
              fmt(3 2 3 0 0 0 0)) ///
        title("`ptitle_`ch''") `t4extra'

    if _rc == 608 {
        di as error "  ** table4_channels_resolution.rtf is OPEN IN WORD — close it and re-run."
        local t4fail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 4: esttab failed for panel `ch' (rc=" _rc ")"
        local t4fail 1
        continue
    }

    local writemode append
}

if `t4fail' == 0 di as result "Table 4 saved: $tabs/table4_channels_resolution.rtf"
else di as error "Table 4 written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* 5. SAVE IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

foreach ch of local channels {
    foreach grp in nd def {
        preserve
            clear
            set obs 6
            gen horizon = _n - 1     // 0 (baseline), 1..5
            foreach m in b lo90 hi90 lo95 hi95 {
                svmat `m'_`grp'_`ch', names(`m')
                rename `m'1 `m'
            }
            gen channel = "`ch'"
            gen group   = "`grp'"
            save "$clean/irf_`grp'_`ch'.dta", replace
        restore
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 6. EXPORT SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

preserve
    clear
    local nrows = 5 * 6   // 5 horizons × 6 channels
    set obs `nrows'
    gen channel  = ""
    gen horizon  = .
    foreach v in b_nd b_def p {
        gen `v' = .
    }

    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel = "`ch'"              in `row'
            replace horizon = `h'                 in `row'
            replace b_nd    = b_nd_`ch'[`h'+2,1]  in `row'
            replace b_def   = b_def_`ch'[`h'+2,1] in `row'
            replace p       = pval_`ch'[`h'+2,1]  in `row'
            local ++row
        }
    }

    order channel horizon b_nd b_def p
    export delimited "$tabs/channels_resolution.csv", replace
    di as result "Table saved: $tabs/channels_resolution.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 7. FIGURE — 2×3 MULTI-PANEL
* ══════════════════════════════════════════════════════════════════════════

local c_nd  "23 55 94"    // navy   — non-default
local c_def "180 60 40"   // brick  — default-linked

local titlelabels `" "Private Credit/GDP" "Bank Claims on Govt/GDP" "Investment/GDP" "Govt Expenditure/GDP" "Primary Balance/GDP" "FDI/GDP" "'

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titlelabels'

    use "$clean/irf_nd_`ch'.dta",  clear
    append using "$clean/irf_def_`ch'.dta"

    * Show legend only in last panel (bottom-right), bottom-anchored (no
    * ring(0)/pos()) rather than overlapping the plotted lines.
    if `i' == 6 {
        local legopt legend(order(3 "Non-default" 4 "Default-linked") ///
                     cols(2) size(vsmall))
    }
    else {
        local legopt legend(off)
    }

    twoway ///
        (rarea lo90 hi90 horizon if group=="nd", ///
            color("`c_nd'%20") lwidth(none)) ///
        (rarea lo90 hi90 horizon if group=="def", ///
            color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if group=="nd", ///
            lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) ///
            lwidth(medthick) msize(small)) ///
        (connected b horizon if group=="def", ///
            lcolor("`c_def'") mcolor("`c_def'") msymbol(square) ///
            lwidth(medthick) msize(small) lpattern(dash)) ///
        , ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)5, labsize(small)) ///
        ylabel(, format(%5.1f) labsize(small)) ///
        xtitle("Years after onset", size(vsmall)) ///
        ytitle("Cumulative change (% or pp — see note)", size(vsmall)) ///
        title(`tlab', size(small) color(navy)) ///
        `legopt' ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ols_`i', replace)

    local ++i
}

graph combine ols_1 ols_2 ols_3 ols_4 ols_5 ols_6, ///
    cols(3) rows(2) ///
    title("Transmission Channels by Resolution Type", ///
          size(medlarge) color(navy)) ///
    note("90% CI. Robust (heteroskedasticity-only) SE. Country FE only (no year FE). Non-default: N=40. Default-linked: N=21." ///
         "Units differ by channel: private credit, bank claims on govt, investment and govt expenditure are LOG REAL LEVELS, so their scale is cumulative percent change (comparable to the GDP result). Primary balance and FDI change sign, so no log is possible and they remain ratios to GDP, in percentage points.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig12a_channels_ols.pdf", replace
di as result "Figure saved: fig12a_channels_ols.pdf"
forvalues i = 1/6 {
    capture graph drop ols_`i'
}

di as result _n "12_channels_resolution.do complete."
