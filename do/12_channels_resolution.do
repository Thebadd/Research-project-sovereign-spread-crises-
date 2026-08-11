/*===========================================================================
  12_CHANNELS_RESOLUTION.DO
  Transmission Channels by Resolution Type (Non-Default vs. Default-Linked)
  with IPW Correction for Selection into Default

  Research question: do the transmission channels of spread crises differ
  between episodes resolved without default and those linked to default?

  STRATEGY:
  ---------
  For each of the 6 channels from 11_channels.do, we run:

  Spec A (OLS): joint LP with onset_nd and onset_def simultaneously
    ch_var(h) = αi + γt + β_nd(h)·onset_nd + β_def(h)·onset_def
               + X_core (common core + pre_<v>)·δ + ε     [DK SE, lag=max(1,h+1)]

  Spec B (IPW): per-type-vs-tranquil (reference-paper structure). Each type is
    scored vs TRANQUIL years with the rival type dropped — two first stages:
      probit onset_nd  $ctrl_core + Z2 if sample & onset_def==0   -> ipw_nd
      probit onset_def $ctrl_core + Z2 if sample & onset_nd==0    -> ipw_def
    Two SEPARATE weighted LPs: onset_nd vs tranquil [aw=ipw_nd], onset_def vs
    tranquil [aw=ipw_def]. areg absorb(cid) vce(cluster cid).

  Extra cost of default = (default line) - (non-default line), Clogg z per horizon.

  OUTPUTS:
  --------
  - fig12a_channels_ols.pdf   : 2×3 grid OLS (nd=navy, def=brick)
  - fig12b_channels_ipw.pdf   : 2×3 grid IPW (nd=navy, def=brick)
  - channels_resolution.csv   : full coefficient table
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking l_govexp l_open l_credit l_lninfl"
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
* ══════════════════════════════════════════════════════════════════════════

foreach var in credit claims_govt inv govexp pb fdi {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    * own-outcome pre-trend control, matching the AIPW spec in 13c.
    capture drop pre_`var'
    gen pre_`var' = L.`var' - L2.`var'
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
* 3. IPW WEIGHTS — ACT 2 (replicated from 08_ipw_lp.do)
*    Full first stage: X (l1_gdpg l2_gdpg debt ca infl imf) + Z2 (l_fedfunds l_reg_crisis_share past_def_onsets); lean debt-ca fallback
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== FIRST STAGE: Pr(default-linked | crisis onset) ==="

capture drop pscore_nd pscore_def ipw_nd ipw_def

* PER-TYPE-vs-TRANQUIL propensities (reference-paper structure): each resolution
* type is scored vs TRANQUIL years with the RIVAL type dropped, so the control is
* clean tranquil country-years (not the other crisis type). Full sample => the full
* $ctrl_core is safe (no among-onsets separation). Two weight sets: ipw_nd, ipw_def.
foreach s in nd def {
    if "`s'" == "nd"  local rival onset_def
    else              local rival onset_nd

    quietly probit onset_`s' $ctrl_core ///
        l_fedfunds l_reg_crisis_share past_def_onsets if sample==1 & `rival'==0, vce(cluster cid)
    quietly lroc, nograph
    di as result "  First stage `s' vs tranquil: AUROC = " %5.3f r(area) "   (N = " e(N) ")"

    predict pscore_`s' if sample==1 & `rival'==0, pr
    quietly replace pscore_`s' = . if (pscore_`s'<0.01 | pscore_`s'>0.99) & !missing(pscore_`s')

    quietly summarize onset_`s' if sample==1 & `rival'==0
    local pmarg = r(mean)
    gen double ipw_`s' = .
    quietly replace ipw_`s' = `pmarg'     / pscore_`s'     if onset_`s'==1 & !missing(pscore_`s')
    quietly replace ipw_`s' = (1-`pmarg') / (1-pscore_`s') if onset_all==0 & !missing(pscore_`s')
    label var ipw_`s' "Act 2 stabilized IPW weight (`s' vs tranquil)"
}

* ══════════════════════════════════════════════════════════════════════════
* 4. LP ESTIMATION BY CHANNEL — OLS AND IPW
* ══════════════════════════════════════════════════════════════════════════

local channels   credit claims_govt inv govexp pb fdi

* Controls: common core ($ctrl_core) + each channel's own pre_<v> (same as 11_channels.do)
* Common-core controls (Asonuma-aligned $ctrl_core) + each channel's own pre_<v>;
* the core term equal to the channel's own lagged level is dropped from its own reg.
local ctrl_credit      l1_gdpg l_debt l_ca l_banking l_govexp l_open l_lninfl pre_credit
local ctrl_claims_govt $ctrl_core pre_claims_govt
local ctrl_inv         $ctrl_core pre_inv
local ctrl_govexp      l1_gdpg l_debt l_ca l_banking l_open l_credit l_lninfl pre_govexp
local ctrl_pb          $ctrl_core pre_pb
local ctrl_fdi         $ctrl_core pre_fdi

* Initialize storage matrices
foreach ch of local channels {
    foreach spec in ols ipw {
        foreach grp in nd def {
            foreach m in b lo90 hi90 lo95 hi95 {
                matrix `m'_`grp'_`spec'_`ch' = J(5, 1, .)
            }
        }
        matrix pval_`spec'_`ch' = J(5, 1, .)
    }
}

* ── Loop over channels ───────────────────────────────────────────────────

eststo clear   // capture OLS resolution-split estimates for Table 4

foreach ch of local channels {

    local ctrl `ctrl_`ch''

    di as result _n "========================================"
    di as result "CHANNEL: `ch'"
    di as result "========================================"
    di "h   b_nd_OLS  b_def_OLS  p_ols  b_nd_IPW  b_def_IPW  p_ipw"

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local row = `h' + 1

        * ── OLS: DK SE joint regression ─────────────────────────────────
        capture xtscc ch_`ch'_`h' onset_nd onset_def `ctrl' i.year ///
            if sample == 1, fe lag(`lag')

        if _rc == 0 {
            matrix b_nd_ols_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_ols_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_ols_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix lo95_nd_ols_`ch'[`row',1] = _b[onset_nd]  - 1.960*_se[onset_nd]
            matrix hi95_nd_ols_`ch'[`row',1] = _b[onset_nd]  + 1.960*_se[onset_nd]
            matrix b_def_ols_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_ols_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_ols_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix lo95_def_ols_`ch'[`row',1]= _b[onset_def] - 1.960*_se[onset_def]
            matrix hi95_def_ols_`ch'[`row',1]= _b[onset_def] + 1.960*_se[onset_def]
            test onset_nd = onset_def
            matrix pval_ols_`ch'[`row',1] = r(p)
            local pd_`ch'_`h' = r(p)   // store before eststo (which can reset r())

            * Difference block (Clogg et al. 1995) + episode counts
            local bnd  = _b[onset_nd]
            local bdef = _b[onset_def]
            local snd  = _se[onset_nd]
            local sdef = _se[onset_def]
            local bdiff = `bdef' - `bnd'
            local zdiff = `bdiff' / sqrt(`snd'^2 + `sdef'^2)
            local pz    = 2*(1 - normal(abs(`zdiff')))
            quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepnd = r(N)
            quietly count if onset_def == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepdef = r(N)

            * Capture OLS estimates + difference block for the publication table (Table 4)
            eststo t4_`ch'_`h', title("h=`h'")
            estadd scalar bdiff  = `bdiff'
            estadd scalar zdiff  = `zdiff'
            estadd scalar pzdiff = `pz'
            estadd scalar pdiff  = `pd_`ch'_`h''
            estadd scalar nepnd  = `nepnd'
            estadd scalar nepdef = `nepdef'
            local elist_`ch' `elist_`ch'' t4_`ch'_`h'
            local b_nd_o  = _b[onset_nd]
            local b_def_o = _b[onset_def]
            local p_o     = r(p)
        }
        else {
            local b_nd_o  = .
            local b_def_o = .
            local p_o     = .
            di as error "OLS failed for `ch' h=`h'"
        }

        * ── IPW: two SEPARATE vs-tranquil weighted LPs (rival dropped) ────
        * Country FE only, no year FE (paper-aligned two-stage outcome regression).
        * non-default vs tranquil (drop onset_def), weight ipw_nd
        capture areg ch_`ch'_`h' onset_nd `ctrl' ///
            [aw=ipw_nd] if sample==1 & onset_def==0 & !missing(ipw_nd), ///
            absorb(cid) vce(cluster cid)
        if _rc == 0 {
            matrix b_nd_ipw_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_ipw_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_ipw_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix lo95_nd_ipw_`ch'[`row',1] = _b[onset_nd]  - 1.960*_se[onset_nd]
            matrix hi95_nd_ipw_`ch'[`row',1] = _b[onset_nd]  + 1.960*_se[onset_nd]
            local b_nd_w  = _b[onset_nd]
            local se_nd_w = _se[onset_nd]
        }
        else {
            local b_nd_w = .
            local se_nd_w = .
            di as error "IPW nd-vs-tranquil failed for `ch' h=`h'"
        }
        * default-linked vs tranquil (drop onset_nd), weight ipw_def
        capture areg ch_`ch'_`h' onset_def `ctrl' ///
            [aw=ipw_def] if sample==1 & onset_nd==0 & !missing(ipw_def), ///
            absorb(cid) vce(cluster cid)
        if _rc == 0 {
            matrix b_def_ipw_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_ipw_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_ipw_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix lo95_def_ipw_`ch'[`row',1]= _b[onset_def] - 1.960*_se[onset_def]
            matrix hi95_def_ipw_`ch'[`row',1]= _b[onset_def] + 1.960*_se[onset_def]
            local b_def_w  = _b[onset_def]
            local se_def_w = _se[onset_def]
        }
        else {
            local b_def_w = .
            local se_def_w = .
            di as error "IPW def-vs-tranquil failed for `ch' h=`h'"
        }
        * extra cost of default (IPW) = def - nd, Clogg z on the two vs-tranquil lines
        local p_w = .
        if !missing(`b_nd_w') & !missing(`b_def_w') {
            local zdiff = (`b_def_w' - `b_nd_w') / sqrt(`se_nd_w'^2 + `se_def_w'^2)
            local p_w   = 2*(1 - normal(abs(`zdiff')))
        }
        matrix pval_ipw_`ch'[`row',1] = `p_w'

        di "h=" `h' ///
           "  " %7.3f `b_nd_o'  "  " %7.3f `b_def_o' "  " %5.3f `p_o' ///
           "  " %7.3f `b_nd_w'  "  " %7.3f `b_def_w' "  " %5.3f `p_w'
    }
}

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — TABLE 4: Transmission channels by resolution type
*   Word/RTF, multi-panel: one panel per channel, columns = horizons h=0..4.
*   Each panel reports non-default and default-linked onset coefficients
*   (entered jointly, DK SE in parentheses) plus the p-value of their equality.
*   OLS spec (matches Tables 1-3). First panel replaces; the rest append.
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

local t4note "Dependent variable: cumulative change in the channel variable (pp) from t-1 to t+h. Both onset dummies enter jointly. Jorda (2005) local projections; country and year fixed effects; common-core controls plus the channel's own pre-crisis change; continuation years excluded. Driscoll-Kraay standard errors in parentheses. p(nd=def) is the p-value of the equality test. * p<0.10, ** p<0.05, *** p<0.01."

* Per-channel panel titles (Panel A carries the overall table caption)
local ptitle_credit      "Table 4. Channels by resolution (nd vs. def) -- Panel A: Private credit/GDP"
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
        stats(bdiff zdiff pzdiff pdiff nepnd nepdef N N_g, ///
              labels("Difference (default - non-default)" "  Clogg et al. (1995) z" ///
                     "  p (Clogg z)" "  p (Wald, nd = def)" ///
                     "Episodes (non-default)" "Episodes (default)" ///
                     "Observations" "Countries") ///
              fmt(3 3 3 3 0 0 0 0)) ///
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
    foreach spec in ols ipw {
        foreach grp in nd def {
            preserve
                clear
                set obs 5
                gen horizon = _n - 1
                foreach m in b lo90 hi90 lo95 hi95 {
                    svmat `m'_`grp'_`spec'_`ch', names(`m')
                    rename `m'1 `m'
                }
                gen channel  = "`ch'"
                gen spec_tag = "`spec'"
                gen group    = "`grp'"
                save "$clean/irf_`grp'_`spec'_`ch'.dta", replace
            restore
        }
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
    foreach v in b_nd_ols b_def_ols p_ols b_nd_ipw b_def_ipw p_ipw {
        gen `v' = .
    }

    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel  = "`ch'"                         in `row'
            replace horizon  = `h'                            in `row'
            replace b_nd_ols  = b_nd_ols_`ch'[`h'+1,1]       in `row'
            replace b_def_ols = b_def_ols_`ch'[`h'+1,1]      in `row'
            replace p_ols     = pval_ols_`ch'[`h'+1,1]       in `row'
            replace b_nd_ipw  = b_nd_ipw_`ch'[`h'+1,1]       in `row'
            replace b_def_ipw = b_def_ipw_`ch'[`h'+1,1]      in `row'
            replace p_ipw     = pval_ipw_`ch'[`h'+1,1]       in `row'
            local ++row
        }
    }

    order channel horizon b_nd_ols b_def_ols p_ols b_nd_ipw b_def_ipw p_ipw
    export delimited "$tabs/channels_resolution.csv", replace
    di as result "Table saved: $tabs/channels_resolution.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 7. FIGURES — 2×3 MULTI-PANEL, OLS THEN IPW
* ══════════════════════════════════════════════════════════════════════════

local c_nd  "23 55 94"    // navy   — non-default
local c_def "180 60 40"   // brick  — default-linked

local titlelabels `" "Private Credit/GDP" "Bank Claims on Govt/GDP" "Investment/GDP" "Govt Expenditure/GDP" "Primary Balance/GDP" "FDI/GDP" "'

* ── Figure A: OLS ────────────────────────────────────────────────────────

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titlelabels'

    use "$clean/irf_nd_ols_`ch'.dta",  clear
    append using "$clean/irf_def_ols_`ch'.dta"

    * Show legend only in last panel (bottom-right)
    if `i' == 6 {
        local legopt legend(order(3 "Non-default" 4 "Default-linked") ///
                     ring(0) pos(7) cols(1) size(vsmall) ///
                     region(lcolor(gs12)))
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
        xlabel(0(1)4, labsize(small)) ///
        ylabel(, format(%5.1f) labsize(small)) ///
        xtitle("Years after onset", size(vsmall)) ///
        ytitle("Cumulative change (pp)", size(vsmall)) ///
        title(`tlab', size(small) color(navy)) ///
        `legopt' ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ols_`i', replace)

    local ++i
}

graph combine ols_1 ols_2 ols_3 ols_4 ols_5 ols_6, ///
    cols(3) rows(2) ///
    title("Transmission Channels by Resolution Type — OLS", ///
          size(medlarge) color(navy)) ///
    note("90% CI. DK SE. Country & year FE. Non-default: N=40. Default-linked: N=21.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig12a_channels_ols.pdf", replace
di as result "Figure saved: fig12a_channels_ols.pdf"
forvalues i = 1/6 {
    capture graph drop ols_`i'
}

* ── Figure B: IPW ────────────────────────────────────────────────────────

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titlelabels'

    use "$clean/irf_nd_ipw_`ch'.dta",  clear
    append using "$clean/irf_def_ipw_`ch'.dta"

    * Show legend only in last panel (bottom-right)
    if `i' == 6 {
        local legopt legend(order(3 "Non-default (IPW)" 4 "Default-linked (IPW)") ///
                     ring(0) pos(7) cols(1) size(vsmall) ///
                     region(lcolor(gs12)))
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
        xlabel(0(1)4, labsize(small)) ///
        ylabel(, format(%5.1f) labsize(small)) ///
        xtitle("Years after onset", size(vsmall)) ///
        ytitle("Cumulative change (pp)", size(vsmall)) ///
        title(`tlab', size(small) color(navy)) ///
        `legopt' ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ipw_`i', replace)

    local ++i
}

graph combine ipw_1 ipw_2 ipw_3 ipw_4 ipw_5 ipw_6, ///
    cols(3) rows(2) ///
    title("Transmission Channels by Resolution Type — IPW", ///
          size(medlarge) color(navy)) ///
    note("90% CI. Clustered SE. Country & year FE. IPW: probit(def | crisis) on macro controls + fedfunds, contagion, past default onsets.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(7)

graph export "$figs/fig12b_channels_ipw.pdf", replace
di as result "Figure saved: fig12b_channels_ipw.pdf"
forvalues i = 1/6 {
    capture graph drop ipw_`i'
}

di as result _n "12_channels_resolution.do complete."
