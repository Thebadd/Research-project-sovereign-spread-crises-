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
    - Act 1 (pooled, all 61 onsets):  xtscc fe, DK SE, lag = max(1,h+1)
    - Act 2 (resolution split):       joint onset_nd + onset_def
                                      OLS (xtscc, DK) and IPW (areg, cluster)
  Coverage is thin for this channel (IMF data start 2001; 6 panel
  countries absent), so the resolution cells are small — interpret with care.

  Run AFTER 01c_merge_nexus.do (and after 08_ipw_lp.do for context).
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ── 1. Channel outcome variables: cumulative change from t-1 ──────────────
foreach var in claimsgov_assets claimpriv_assets {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
        label var ch_`var'_`h' "Cum. change `var': F`h' vs t-1 (pp)"
    }
}

* ── Channel-specific controls (lean, to preserve the short 2001+ sample) ──
local channels        claimsgov_assets claimpriv_assets
local ctrl_claimsgov_assets  L.claimsgov_assets l1_gdpg debt pb banking_crisis
local ctrl_claimpriv_assets  L.claimpriv_assets l1_gdpg l2_gdpg debt ca banking_crisis

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
        matrix `m'_`ch' = J(5, 1, .)
    }
}

di as result _n "========================================================"
di as result "ACT 1 — POOLED NEXUS CHANNELS (xtscc fe, DK SE)"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    di as result _n "--- CHANNEL: `ch' ---"
    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        capture xtscc ch_`ch'_`h' onset_all `ctrl' i.year if sample==1, fe lag(`lag')
        if _rc == 0 {
            matrix b_`ch'[`h'+1,1]    = _b[onset_all]
            matrix lo90_`ch'[`h'+1,1] = _b[onset_all] - 1.645*_se[onset_all]
            matrix hi90_`ch'[`h'+1,1] = _b[onset_all] + 1.645*_se[onset_all]
            matrix lo95_`ch'[`h'+1,1] = _b[onset_all] - 1.960*_se[onset_all]
            matrix hi95_`ch'[`h'+1,1] = _b[onset_all] + 1.960*_se[onset_all]
            di "h=" `h' ": beta=" %7.3f _b[onset_all] "  SE=" %6.3f _se[onset_all] ///
               "  p=" %5.3f (2*(1-normal(abs(_b[onset_all]/_se[onset_all])))) "  N=" e(N)
        }
        else di as error "h=" `h' ": xtscc failed for `ch' (rc=" _rc ")"
    }
}

* Save pooled IRF datasets
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
        xlabel(0(1)4, labsize(medsmall)) ylabel(, format(%5.1f) labsize(medsmall)) ///
        xtitle("Years after onset", size(small)) ///
        ytitle("Cumulative change (pp)", size(small)) ///
        title(`tlab', size(medsmall) color(navy)) legend(off) ///
        graphregion(color(white)) plotregion(color(white)) name(nx_`i', replace)
    local ++i
}
graph combine nx_1 nx_2, cols(2) ///
    title("Sovereign-Bank Nexus Channels (pooled)", size(medlarge) color(navy)) ///
    note("90%/95% CI. DK SE. Country & year FE. IMF MFS, 2001-2024.", size(vsmall)) ///
    graphregion(color(white)) xsize(10) ysize(4)
graph export "$figs/fig11b_nexus_pooled.pdf", replace
forvalues i = 1/2 {
    capture graph drop nx_`i'
}
di as result "Figure saved: fig11b_nexus_pooled.pdf"

* ══════════════════════════════════════════════════════════════════════════
* 3. ACT 2 IPW WEIGHTS (rebuilt as in 08_ipw_lp.do / 12_channels_resolution)
*    Parsimonious first stage: Pr(default-linked | crisis onset) on debt, CA
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
}

capture drop pscore2 trimmed2 ipw2
quietly probit onset_def debt ca if onset_all == 1, vce(robust)
predict pscore2 if onset_all == 1, pr
gen trimmed2 = (pscore2 < 0.05 | pscore2 > 0.95) if !missing(pscore2)
replace pscore2 = . if trimmed2 == 1
quietly summarize onset_def if onset_all == 1
local p_def = r(mean)
local p_nd  = 1 - `p_def'
gen ipw2 = .
replace ipw2 = `p_def' / pscore2       if onset_def == 1 & !missing(pscore2)
replace ipw2 = `p_nd'  / (1 - pscore2) if onset_nd  == 1 & !missing(pscore2)

* ══════════════════════════════════════════════════════════════════════════
* 4. ACT 2 — RESOLUTION SPLIT (joint nd/def), OLS and IPW
* ══════════════════════════════════════════════════════════════════════════
foreach ch of local channels {
    foreach spec in ols ipw {
        foreach grp in nd def {
            foreach m in b lo90 hi90 {
                matrix `m'_`grp'_`spec'_`ch' = J(5, 1, .)
            }
        }
        matrix pval_`spec'_`ch' = J(5, 1, .)
    }
}

di as result _n "========================================================"
di as result "ACT 2 — NEXUS CHANNELS BY RESOLUTION TYPE"
di as result "========================================================"

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    di as result _n "--- CHANNEL: `ch' ---"
    di "h   b_nd_OLS  b_def_OLS  p_OLS   b_nd_IPW  b_def_IPW  p_IPW"

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local row = `h' + 1

        * OLS (xtscc, DK)
        capture xtscc ch_`ch'_`h' onset_nd onset_def `ctrl' i.year ///
            if sample == 1, fe lag(`lag')
        if _rc == 0 {
            matrix b_nd_ols_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_ols_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_ols_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix b_def_ols_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_ols_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_ols_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            quietly test onset_nd = onset_def
            matrix pval_ols_`ch'[`row',1] = r(p)
            local b_nd_o = _b[onset_nd]
            local b_def_o = _b[onset_def]
            local p_o = r(p)
        }
        else {
            local b_nd_o = .
            local b_def_o = .
            local p_o = .
        }

        * IPW (areg, cluster) — control obs get weight 1
        quietly replace ipw2 = 1 if onset_all == 0 & sample == 1
        capture areg ch_`ch'_`h' onset_nd onset_def `ctrl' i.year ///
            [aw=ipw2] if sample == 1 & !missing(ipw2), absorb(cid) vce(cluster cid)
        if _rc == 0 {
            quietly test onset_nd = onset_def
            local p_w = r(p)
            matrix b_nd_ipw_`ch'[`row',1]    = _b[onset_nd]
            matrix lo90_nd_ipw_`ch'[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
            matrix hi90_nd_ipw_`ch'[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
            matrix b_def_ipw_`ch'[`row',1]   = _b[onset_def]
            matrix lo90_def_ipw_`ch'[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
            matrix hi90_def_ipw_`ch'[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
            matrix pval_ipw_`ch'[`row',1] = `p_w'
            local b_nd_w = _b[onset_nd]
            local b_def_w = _b[onset_def]
        }
        else {
            local b_nd_w = .
            local b_def_w = .
            local p_w = .
        }
        quietly replace ipw2 = . if onset_all == 0

        di "h=" `h' "  " %7.3f `b_nd_o' "  " %7.3f `b_def_o' "  " %5.3f `p_o' ///
                  "  " %7.3f `b_nd_w' "  " %7.3f `b_def_w' "  " %5.3f `p_w'
    }
}

* ── Export resolution summary table ───────────────────────────────────────
preserve
    clear
    set obs 10                          // 5 horizons x 2 channels
    gen channel = ""
    gen horizon = .
    foreach v in b_nd_ols b_def_ols p_ols b_nd_ipw b_def_ipw p_ipw {
        gen `v' = .
    }
    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel   = "`ch'"               in `row'
            replace horizon   = `h'                  in `row'
            replace b_nd_ols  = b_nd_ols_`ch'[`h'+1,1]  in `row'
            replace b_def_ols = b_def_ols_`ch'[`h'+1,1] in `row'
            replace p_ols     = pval_ols_`ch'[`h'+1,1]  in `row'
            replace b_nd_ipw  = b_nd_ipw_`ch'[`h'+1,1]  in `row'
            replace b_def_ipw = b_def_ipw_`ch'[`h'+1,1] in `row'
            replace p_ipw     = pval_ipw_`ch'[`h'+1,1]  in `row'
            local ++row
        }
    }
    order channel horizon b_nd_ols b_def_ols p_ols b_nd_ipw b_def_ipw p_ipw
    export delimited "$tabs/nexus_channels_resolution.csv", replace
    di as result "Table saved: $tabs/nexus_channels_resolution.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 5. SAVE RESOLUTION IRF DATASETS (one per channel x spec x group)
* ══════════════════════════════════════════════════════════════════════════
foreach ch of local channels {
    foreach spec in ols ipw {
        foreach grp in nd def {
            preserve
                clear
                set obs 5
                gen horizon = _n - 1
                foreach m in b lo90 hi90 {
                    svmat `m'_`grp'_`spec'_`ch', names(`m')
                    rename `m'1 `m'
                }
                gen channel  = "`ch'"
                gen spec_tag = "`spec'"
                gen group    = "`grp'"
                save "$clean/irf_nx_`grp'_`spec'_`ch'.dta", replace
            restore
        }
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 6. RESOLUTION FIGURES (non-default vs default-linked), OLS then IPW
* ══════════════════════════════════════════════════════════════════════════
local c_nd  "23 55 94"     // navy  — non-default
local c_def "180 60 40"    // brick — default-linked
local titles `" "Claims on Govt / Assets" "Claims on Private / Assets" "'

foreach spec in ols ipw {
    local i = 1
    foreach ch of local channels {
        local tlab : word `i' of `titles'
        use "$clean/irf_nx_nd_`spec'_`ch'.dta", clear
        append using "$clean/irf_nx_def_`spec'_`ch'.dta"
        if `i' == 2 {
            local legopt legend(order(3 "Non-default" 4 "Default-linked") ///
                         ring(0) pos(7) cols(1) size(vsmall) region(lcolor(gs12)))
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
            xlabel(0(1)4, labsize(small)) ylabel(, format(%5.1f) labsize(small)) ///
            xtitle("Years after onset", size(vsmall)) ///
            ytitle("Cumulative change (pp)", size(vsmall)) ///
            title(`tlab', size(small) color(navy)) ///
            `legopt' graphregion(color(white)) plotregion(color(white)) ///
            name(nxr_`spec'_`i', replace)
        local ++i
    }
    if "`spec'" == "ols" local stitle "OLS (DK SE)"
    else                 local stitle "IPW (cluster SE)"
    graph combine nxr_`spec'_1 nxr_`spec'_2, cols(2) ///
        title("Nexus Channels by Resolution — `stitle'", size(medlarge) color(navy)) ///
        note("90% CI. Country & year FE. Non-default vs. default-linked. IMF MFS 2001-2024.", ///
             size(vsmall)) ///
        graphregion(color(white)) xsize(10) ysize(4)
    graph export "$figs/fig11b_nexus_resolution_`spec'.pdf", replace
    forvalues i = 1/2 {
        capture graph drop nxr_`spec'_`i'
    }
    di as result "Figure saved: fig11b_nexus_resolution_`spec'.pdf"
}

di as result _n "11b_nexus_channels.do complete."
