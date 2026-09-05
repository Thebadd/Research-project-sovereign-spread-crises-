/*===========================================================================
  11B_NEXUS_CHANNELS.DO
  Sovereign-bank nexus channels (IMF MFS, merged in 01c_merge_nexus.do)

  Outcomes (cumulative change from t-1, in pp of the ratio):
    claimsgov_assets  — bank claims on government / assets   (doom-loop)
    claimpriv_assets  — bank claims on private    / assets   (reallocation)

  Reading the two together shows whether banks shift the balance sheet
  toward the sovereign and away from the private sector after a spread
  crisis. claimsgov_assets is the headline doom-loop channel.

  Estimation mirrors 11_channels.do / 12_channels_resolution.do:
    - Act 1 (pooled, all 61 onsets):  xtreg fe, robust SE, country FE only
    - Act 2 (resolution split):       joint onset_nd + onset_def, xtreg fe,
                                      robust SE, country FE only
  Country FE only, plain robust SE, no year FE -- same switch as 02/03/11
  (see 02_lp_all.do's header). IPW REMOVED: this file used to carry a
  parallel IPW-weighted comparison, dropped project-wide once 08b_aipw.do's
  doubly-robust AIPW estimator superseded plain IPW -- see METHODOLOGY.md.
  Coverage is thin for this channel (IMF data start 2001; 6 panel
  countries absent), so the resolution cells are small — interpret with care.

  Run AFTER 01c_merge_nexus.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ── 1. Channel outcome variables: cumulative change from t-1 ──────────────
* NOT rescaled to log real levels, deliberately. The GDP-ratio channels in 11/12
* are built as ln(gdp_real*x)*100 so a collapsing denominator cannot masquerade as
* a channel response. These two are shares of TOTAL BANK ASSETS, not of GDP, so
* that problem does not arise: the denominator is the banks' own balance sheet and
* the share is exactly the object of interest (how banks reallocate between the
* sovereign and the private sector). They stay in percentage points of assets.
foreach var in claimsgov_assets claimpriv_assets {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
        label var ch_`var'_`h' "Cum. change `var': F`h' vs t-1 (pp)"
    }
    * own-outcome pre-crisis change (Asonuma's g_0), matching 13c/13d. Needs a
    * second lag, so it costs some coverage on the thin 2001+ nexus sample.
    capture drop pre_`var'
    gen pre_`var' = L.`var' - L2.`var'
}

* ── Outcome-model controls: common core ($ctrl_core) + each nexus channel's pre_<v> ──
local channels        claimsgov_assets claimpriv_assets
* Common-core controls ($ctrl_core, Asonuma-aligned) + each nexus channel's pre_<v>.
local ctrl_claimsgov_assets  $ctrl_core pre_claimsgov_assets
local ctrl_claimpriv_assets  $ctrl_core pre_claimpriv_assets

* ── Coverage at onset ─────────────────────────────────────────────────────
di as result _n "=== NEXUS CHANNEL COVERAGE AT ONSET (h=0) ==="
foreach ch of local channels {
    quietly count if onset_all == 1 & sample == 1 & !missing(ch_`ch'_0)
    di as result "  `ch': " r(N) " / 61 onsets"
}

* ══════════════════════════════════════════════════════════════════════════
* 2. ACT 1 — POOLED LP (all onsets)
* ══════════════════════════════════════════════════════════════════════════
foreach ch of local channels {
    foreach m in b lo90 hi90 lo95 hi95 {
        matrix `m'_`ch' = J(6, 1, 0)
    }
}

di as result _n "========================================================"
di as result "ACT 1 — POOLED NEXUS CHANNELS (xtreg fe, robust SE)"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    di as result _n "--- CHANNEL: `ch' ---"
    forvalues h = 0/4 {
        capture xtreg ch_`ch'_`h' onset_all `ctrl' if sample==1, fe vce(robust)
        if _rc == 0 {
            matrix b_`ch'[`h'+2,1]    = _b[onset_all]
            matrix lo90_`ch'[`h'+2,1] = _b[onset_all] - 1.645*_se[onset_all]
            matrix hi90_`ch'[`h'+2,1] = _b[onset_all] + 1.645*_se[onset_all]
            matrix lo95_`ch'[`h'+2,1] = _b[onset_all] - 1.960*_se[onset_all]
            matrix hi95_`ch'[`h'+2,1] = _b[onset_all] + 1.960*_se[onset_all]
            di "h=" `h'+1 ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
               "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
        }
        else di as error "h=" `h'+1 ": xtreg failed for `ch' (rc=" _rc ")"
    }
}

* Save pooled IRF datasets
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
        save "$clean/irf_nx_`ch'.dta", replace
    restore
}

* ── Pooled figure (1x2) ───────────────────────────────────────────────────
local c_main "23 55 94"
local titles `" "Claims on Govt / Assets" "Claims on Private / Assets" "'
local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titles'
    use "$clean/irf_nx_`ch'.dta", clear
    twoway ///
        (rarea lo95 hi95 horizon, color("`c_main'%15") lwidth(none)) ///
        (rarea lo90 hi90 horizon, color("`c_main'%28") lwidth(none)) ///
        (connected b horizon, lcolor("`c_main'") mcolor("`c_main'") ///
            msymbol(circle) lwidth(medthick)), ///
        yline(0, lcolor(gs8) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%5.1f) labsize(medsmall)) ///
        xtitle("Year (Year 1 = crisis year)", size(small)) ///
        ytitle("Cumulative change (pp)", size(small)) ///
        title(`tlab', size(medsmall) color(navy)) legend(off) ///
        graphregion(color(white)) plotregion(color(white)) name(nx_`i', replace)
    local ++i
}
graph combine nx_1 nx_2, cols(2) ///
    title("Sovereign-Bank Nexus Channels (pooled)", size(medlarge) color(navy)) ///
    note("90%/95% CI. Robust SE. Country FE only (no year FE). IMF MFS, 2001-2024.", size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(4)
graph export "$figs/fig11b_nexus_pooled.pdf", replace
forvalues i = 1/2 {
    capture graph drop nx_`i'
}
di as result "Figure saved: fig11b_nexus_pooled.pdf"

* ══════════════════════════════════════════════════════════════════════════
* 3. ACT 2 — RESOLUTION SPLIT (joint nd/def)
* ══════════════════════════════════════════════════════════════════════════
* Reload the panel: the figure loop above left an IRF dataset in memory.
use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year
foreach var in claimsgov_assets claimpriv_assets {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
    }
    capture drop pre_`var'
    gen pre_`var' = L.`var' - L2.`var'
}

foreach ch of local channels {
    foreach grp in nd def {
        foreach m in b lo90 hi90 {
            matrix `m'_`grp'_`ch' = J(6, 1, 0)
        }
    }
    matrix pval_`ch' = J(6, 1, .)
}

di as result _n "========================================================"
di as result "ACT 2 — NEXUS CHANNELS BY RESOLUTION TYPE"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    di as result _n "--- CHANNEL: `ch' ---"
    di "h   b_nd      b_def     p(nd=def)"

    forvalues h = 0/4 {
        local row = `h' + 2

        * JOINT LP, both type dummies, FULL sample, tranquil omitted — the
        * reference paper's baseline (reg g_h dum1 dum2 dum3 g_0 $convar, country
        * dummies, vce(robust), no year FE). The difference test
        * (test onset_nd = onset_def) uses the exact, covariance-correct
        * joint-regression F-statistic -- the paper's own difference-test
        * convention (Table I1) -- not Clogg z/bootstrap.
        capture xtreg ch_`ch'_`h' onset_nd onset_def `ctrl' ///
            if sample == 1, fe vce(robust)
        if _rc == 0 {
            matrix b_nd_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix b_def_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            quietly test onset_nd = onset_def
            matrix pval_`ch'[`row',1] = r(p)
            local b_nd_o  = _b[onset_nd]
            local b_def_o = _b[onset_def]
            local p_o     = r(p)
        }
        else {
            local b_nd_o  = .
            local b_def_o = .
            local p_o     = .
        }

        di "h=" `h'+1 "  " %7.3f `b_nd_o' "  " %7.3f `b_def_o' "  " %5.3f `p_o'
    }
}

* ── Export resolution summary table ───────────────────────────────────────
preserve
    clear
    set obs 10                          // 5 horizons x 2 channels
    gen channel = ""
    gen horizon = .
    foreach v in b_nd b_def p {
        gen `v' = .
    }
    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel = "`ch'"              in `row'
            replace horizon = `h'+1               in `row'
            replace b_nd    = b_nd_`ch'[`h'+2,1]  in `row'
            replace b_def   = b_def_`ch'[`h'+2,1] in `row'
            replace p       = pval_`ch'[`h'+2,1]  in `row'
            local ++row
        }
    }
    order channel horizon b_nd b_def p
    export delimited "$tabs/nexus_channels_resolution.csv", replace
    di as result "Table saved: $tabs/nexus_channels_resolution.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 4. SAVE RESOLUTION IRF DATASETS (one per channel x group)
* ══════════════════════════════════════════════════════════════════════════
foreach ch of local channels {
    foreach grp in nd def {
        preserve
            clear
            set obs 6
            gen horizon = _n - 1     // 0 (baseline), 1..5
            foreach m in b lo90 hi90 {
                svmat `m'_`grp'_`ch', names(`m')
                rename `m'1 `m'
            }
            gen channel = "`ch'"
            gen group   = "`grp'"
            save "$clean/irf_nx_`grp'_`ch'.dta", replace
        restore
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 5. RESOLUTION FIGURE (non-default vs default-linked)
* ══════════════════════════════════════════════════════════════════════════
local c_nd  "23 55 94"     // navy  — non-default
local c_def "180 60 40"    // brick — default-linked
local titles `" "Claims on Govt / Assets" "Claims on Private / Assets" "'

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titles'
    use "$clean/irf_nx_nd_`ch'.dta", clear
    append using "$clean/irf_nx_def_`ch'.dta"
    if `i' == 2 {
        * Bottom-anchored (no ring(0)/pos()), matching the fix applied
        * elsewhere in the project for nd/def legends that used to
        * overlap the plotted lines.
        local legopt legend(order(3 "Non-default" 4 "Default-linked") ///
                     cols(2) size(vsmall))
    }
    else {
        local legopt legend(off)
    }
    twoway ///
        (rarea lo90 hi90 horizon if group=="nd",  color("`c_nd'%20")  lwidth(none)) ///
        (rarea lo90 hi90 horizon if group=="def", color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if group=="nd",  lcolor("`c_nd'")  mcolor("`c_nd'") ///
            msymbol(circle) lwidth(medthick)) ///
        (connected b horizon if group=="def", lcolor("`c_def'") mcolor("`c_def'") ///
            msymbol(square) lwidth(medthick) lpattern(dash)), ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)5, labsize(small)) ylabel(, format(%5.1f) labsize(small)) ///
        xtitle("Years after onset", size(vsmall)) ///
        ytitle("Cumulative change (pp)", size(vsmall)) ///
        title(`tlab', size(small) color(navy)) ///
        `legopt' graphregion(color(white)) plotregion(color(white)) ///
        name(nxr_`i', replace)
    local ++i
}
graph combine nxr_1 nxr_2, cols(2) ///
    title("Nexus Channels by Resolution", size(medlarge) color(navy)) ///
    note("90% CI. Robust SE. Country FE only (no year FE). Non-default vs. default-linked. IMF MFS 2001-2024.", ///
         size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(4)
graph export "$figs/fig11b_nexus_resolution.pdf", replace
forvalues i = 1/2 {
    capture graph drop nxr_`i'
}
di as result "Figure saved: fig11b_nexus_resolution.pdf"

di as result _n "11b_nexus_channels.do complete."
