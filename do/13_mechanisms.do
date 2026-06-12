/*===========================================================================
  13_MECHANISMS.DO
  Mechanistic Tests for Transmission Channels

  TEST 1 — Supply vs. Demand in the Credit Channel
  -------------------------------------------------
  Question: is the credit contraction supply-driven (banks substitute
  sovereign bonds for private loans) or demand-driven (firms stop
  borrowing)?

  Method: add L.claims_govt to the credit LP. If β_credit shrinks
  significantly, sovereign bond accumulation absorbs the credit
  contraction → supply-side (portfolio substitution). If β_credit is
  unchanged → demand-side or independent channel.

  Specs compared:
    Baseline : ch_credit_h = αi + γt + β·onset_all
               + l1_gdpg l2_gdpg debt infl ca banking_crisis + ε
    With clms: same + L.claims_govt

  TEST 2 — Credit as Mediator of Investment Contraction
  ------------------------------------------------------
  Question: does the credit contraction mechanically cause the investment
  decline? If crisis → credit contraction → investment decline, then
  controlling for credit should absorb the investment LP coefficient.

  Method: run investment LP with and without L.credit.
    Without credit: β_inv captures the TOTAL effect of the crisis on
                    investment (direct + through credit)
    With credit:    β_inv captures only the DIRECT effect net of credit

  Mediation share = (β_without - β_with) / β_without

  If mediation share > 0 and meaningful → credit channel mediates
  investment contraction.

  Specs compared:
    Without : ch_inv_h = αi + γt + β·onset_all
              + l1_gdpg l2_gdpg debt ca banking_crisis + ε
    With    : same + L.credit   [current 11_channels.do spec]

  Outputs:
    Printed comparison tables at each horizon
    fig13a_credit_supply_demand.pdf  — credit β with vs without claims_govt
    fig13b_inv_mediation.pdf         — investment β with vs without credit
    "$clean/irf_mech_*.dta"          — IRF datasets for both tests
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ── Generate outcome variables ───────────────────────────────────────────

foreach var in credit claims_govt inv {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
    }
}

* ══════════════════════════════════════════════════════════════════════════
* TEST 1 — SUPPLY VS. DEMAND: CREDIT CHANNEL
*   Baseline credit spec vs. credit spec + L.claims_govt
* ══════════════════════════════════════════════════════════════════════════

di as result _n "========================================================"
di as result "TEST 1: SUPPLY VS. DEMAND IN THE CREDIT CHANNEL"
di as result "  Baseline credit spec vs. adding L.claims_govt"
di as result "========================================================"
di as result "h    β_baseline  SE_base   β_+clms_govt  SE_clms   absorbed(%)"

* Storage matrices
foreach m in b_base lo90_base hi90_base b_clms lo90_clms hi90_clms {
    matrix `m' = J(5,1,.)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    local row = `h' + 1

    * Baseline: current 11_channels spec
    capture xtscc ch_credit_`h' onset_all ///
        l1_gdpg l2_gdpg debt infl ca banking_crisis ///
        i.year if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_base[`row',1]    = _b[onset_all]
        matrix lo90_base[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_base[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
        local b0  = _b[onset_all]
        local se0 = _se[onset_all]
    }

    * With L.claims_govt added
    capture xtscc ch_credit_`h' onset_all ///
        l1_gdpg l2_gdpg debt infl ca banking_crisis L.claims_govt ///
        i.year if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_clms[`row',1]    = _b[onset_all]
        matrix lo90_clms[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_clms[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
        local b1  = _b[onset_all]
        local se1 = _se[onset_all]

        * Absorption share
        if `b0' != 0 {
            local absorbed = (`b0' - `b1') / `b0' * 100
        }
        else {
            local absorbed = .
        }

        di "h=" `h' "   " %8.3f `b0' "   " %6.3f `se0' ///
               "      " %8.3f `b1' "   " %6.3f `se1' ///
               "    " %6.1f `absorbed' "%"
    }
}

di as result _n "Interpretation:"
di as result "  If absorbed > 30% → supply-side (portfolio substitution)"
di as result "  If absorbed < 10% → demand-side or independent channel"

* ── Save IRF datasets for Test 1 ─────────────────────────────────────────

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 {
        svmat `m'_base, names(`m')
        rename `m'1 `m'
    }
    gen spec = "baseline"
    save "$clean/irf_mech_credit_base.dta", replace
restore

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 {
        svmat `m'_clms, names(`m')
        rename `m'1 `m'
    }
    gen spec = "with_clms"
    save "$clean/irf_mech_credit_clms.dta", replace
restore

* ── Figure: Test 1 ───────────────────────────────────────────────────────

local c_base "23 55 94"
local c_clms "157 36 73"

use "$clean/irf_mech_credit_base.dta", clear
append using "$clean/irf_mech_credit_clms.dta"

twoway ///
    (rarea lo90 hi90 horizon if spec=="baseline", ///
        color("`c_base'%20") lwidth(none)) ///
    (connected b horizon if spec=="baseline", ///
        lcolor("`c_base'") lwidth(medthick) msymbol(circle) mcolor("`c_base'")) ///
    (rarea lo90 hi90 horizon if spec=="with_clms", ///
        color("`c_clms'%20") lwidth(none)) ///
    (connected b horizon if spec=="with_clms", ///
        lcolor("`c_clms'") lwidth(medthick) lpattern(dash) ///
        msymbol(square) mcolor("`c_clms'")), ///
    yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%5.2f) labsize(medsmall)) ///
    xtitle("Years after onset", size(small)) ///
    ytitle("Cumulative change in credit/GDP (pp)", size(small)) ///
    title("Credit Channel: Supply vs. Demand Test", size(medium) color(navy)) ///
    subtitle("Does adding bank sovereign exposure absorb the credit effect?", size(small)) ///
    legend(order(2 "Baseline (no claims_govt)" 4 "+ L.claims_govt control") ///
           ring(0) pos(7) cols(1) size(small) region(lcolor(none) fcolor(none))) ///
    note("Blue = baseline credit spec. Red = adding lagged bank sovereign bond holdings." ///
         "If red line closer to zero → supply-side (portfolio substitution) channel." ///
         "DK SE. Country & year FE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig13a_credit_supply_demand.pdf", replace
di as result "Figure saved: fig13a_credit_supply_demand.pdf"

* ══════════════════════════════════════════════════════════════════════════
* TEST 2 — CREDIT AS MEDIATOR OF INVESTMENT CONTRACTION
*   Investment LP without credit vs. with L.credit (current spec)
* ══════════════════════════════════════════════════════════════════════════

di as result _n "========================================================"
di as result "TEST 2: CREDIT AS MEDIATOR OF INVESTMENT CONTRACTION"
di as result "  Investment LP without vs. with L.credit"
di as result "========================================================"
di as result "h    β_no_credit  SE_no    β_with_credit  SE_with  mediated(%)"

* Storage matrices
foreach m in b_noc lo90_noc hi90_noc b_wic lo90_wic hi90_wic {
    matrix `m' = J(5,1,.)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    local row = `h' + 1

    * Without L.credit (total effect)
    capture xtscc ch_inv_`h' onset_all ///
        l1_gdpg l2_gdpg debt ca banking_crisis ///
        i.year if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_noc[`row',1]    = _b[onset_all]
        matrix lo90_noc[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_noc[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
        local b0  = _b[onset_all]
        local se0 = _se[onset_all]
    }

    * With L.credit (direct effect net of credit)
    capture xtscc ch_inv_`h' onset_all ///
        l1_gdpg l2_gdpg debt ca banking_crisis L.credit ///
        i.year if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_wic[`row',1]    = _b[onset_all]
        matrix lo90_wic[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_wic[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
        local b1  = _b[onset_all]
        local se1 = _se[onset_all]

        * Mediation share = (total - direct) / total
        if `b0' != 0 {
            local mediated = (`b0' - `b1') / `b0' * 100
        }
        else {
            local mediated = .
        }

        di "h=" `h' "   " %9.3f `b0' "   " %6.3f `se0' ///
               "      " %9.3f `b1' "   " %6.3f `se1' ///
               "    " %6.1f `mediated' "%"
    }
}

di as result _n "Interpretation:"
di as result "  mediated > 30% → credit is a significant mediator of investment"
di as result "  mediated < 10% → investment contraction is direct (not through credit)"
di as result "  Causal chain if mediated: crisis → credit ↓ → investment ↓"

* ── Save IRF datasets for Test 2 ─────────────────────────────────────────

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 {
        svmat `m'_noc, names(`m')
        rename `m'1 `m'
    }
    gen spec = "no_credit"
    save "$clean/irf_mech_inv_noc.dta", replace
restore

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b lo90 hi90 {
        svmat `m'_wic, names(`m')
        rename `m'1 `m'
    }
    gen spec = "with_credit"
    save "$clean/irf_mech_inv_wic.dta", replace
restore

* ── Figure: Test 2 ───────────────────────────────────────────────────────

local c_noc "23 55 94"
local c_wic "157 36 73"

use "$clean/irf_mech_inv_noc.dta", clear
append using "$clean/irf_mech_inv_wic.dta"

twoway ///
    (rarea lo90 hi90 horizon if spec=="no_credit", ///
        color("`c_noc'%20") lwidth(none)) ///
    (connected b horizon if spec=="no_credit", ///
        lcolor("`c_noc'") lwidth(medthick) msymbol(circle) mcolor("`c_noc'")) ///
    (rarea lo90 hi90 horizon if spec=="with_credit", ///
        color("`c_wic'%20") lwidth(none)) ///
    (connected b horizon if spec=="with_credit", ///
        lcolor("`c_wic'") lwidth(medthick) lpattern(dash) ///
        msymbol(square) mcolor("`c_wic'")), ///
    yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%5.2f) labsize(medsmall)) ///
    xtitle("Years after onset", size(small)) ///
    ytitle("Cumulative change in investment/GDP (pp)", size(small)) ///
    title("Investment Channel: Credit Mediation Test", size(medium) color(navy)) ///
    subtitle("Does the credit contraction explain the investment decline?", size(small)) ///
    legend(order(2 "Without L.credit (total effect)" 4 "With L.credit (direct effect)") ///
           ring(0) pos(7) cols(1) size(small) region(lcolor(none) fcolor(none))) ///
    note("Blue = total effect of crisis on investment." ///
         "Red = direct effect after absorbing credit channel." ///
         "Gap = share of investment contraction mediated by credit. DK SE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig13b_inv_mediation.pdf", replace
di as result "Figure saved: fig13b_inv_mediation.pdf"

di as result _n "13_mechanisms.do complete."
