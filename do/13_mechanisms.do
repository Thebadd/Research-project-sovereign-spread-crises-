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
               + l1_gdpg l2_gdpg debt infl ca + ε
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
              + l1_gdpg l2_gdpg debt ca + ε
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

* ── Reload panel after Test 1 figure section changed the dataset ─────────
use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

foreach var in inv {
    capture drop `var'_base
    gen `var'_base = L.`var'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`var' - `var'_base
    }
}

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

* ══════════════════════════════════════════════════════════════════════════
* TEST 3 — CURRENT ACCOUNT LP
*   Theory: Aguiar-Gopinath (2006) predict spread crises force current
*   account adjustment toward surplus (forced deleveraging). Test whether
*   the CA moves toward surplus following onset, overall and by episode type.
*
*   Note: ca is a CONTROL in channel regressions; here it is the OUTCOME
*   so it is excluded from the right-hand side.
*
*   Specs:
*     Aggregate : ch_ca_h = αi + γt + β·onset_all + l1_gdpg l2_gdpg debt + ε
*     By type   : same with onset_nd and onset_def separately
*
*   Prediction: β > 0 (CA moves toward surplus = forced deleveraging)
*               stronger for default episodes (harder market exclusion)
*
*   Output: fig13c_ca_lp.pdf, irf_mech_ca_*.dta
* ══════════════════════════════════════════════════════════════════════════

* ── Reload panel ─────────────────────────────────────────────────────────
use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* Generate CA outcome variables
capture drop ca_base
gen ca_base = L.ca
forvalues h = 0/4 {
    capture drop ch_ca_`h'
    gen ch_ca_`h' = F`h'.ca - ca_base
}

di as result _n "========================================================"
di as result "TEST 3: CURRENT ACCOUNT LP (Aguiar-Gopinath mechanism)"
di as result "  Prediction: CA moves toward surplus after onset (β > 0)"
di as result "  Controls include L.ca to absorb CA persistence"
di as result "========================================================"
di as result "h    β_all    SE      R2_all   β_nd     SE       β_def    SE      R2_split"

* Storage matrices
foreach m in b_all lo90_all hi90_all b_nd lo90_nd hi90_nd b_def lo90_def hi90_def {
    matrix `m' = J(5,1,.)
}

forvalues h = 0/4 {
    local lag = max(1, `h'+1)
    local row = `h' + 1

    * Aggregate — with lagged CA for persistence
    capture xtscc ch_ca_`h' onset_all ///
        l1_gdpg l2_gdpg debt L.ca ///
        i.year if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_all[`row',1]    = _b[onset_all]
        matrix lo90_all[`row',1] = _b[onset_all] - 1.645*_se[onset_all]
        matrix hi90_all[`row',1] = _b[onset_all] + 1.645*_se[onset_all]
        local b0   = _b[onset_all]
        local se0  = _se[onset_all]
        local r2_0 = e(r2_w)
    }

    * Split by episode type — with lagged CA
    capture xtscc ch_ca_`h' onset_nd onset_def ///
        l1_gdpg l2_gdpg debt L.ca ///
        i.year if sample==1, fe lag(`lag')
    if _rc == 0 {
        matrix b_nd[`row',1]    = _b[onset_nd]
        matrix lo90_nd[`row',1] = _b[onset_nd] - 1.645*_se[onset_nd]
        matrix hi90_nd[`row',1] = _b[onset_nd] + 1.645*_se[onset_nd]
        matrix b_def[`row',1]    = _b[onset_def]
        matrix lo90_def[`row',1] = _b[onset_def] - 1.645*_se[onset_def]
        matrix hi90_def[`row',1] = _b[onset_def] + 1.645*_se[onset_def]
        local b1   = _b[onset_nd]
        local se1  = _se[onset_nd]
        local b2   = _b[onset_def]
        local se2  = _se[onset_def]
        local r2_1 = e(r2_w)

        di "h=" `h' "  " %6.3f `b0' "  " %5.3f `se0' "  " %5.3f `r2_0' ///
               "  " %6.3f `b1' "  " %5.3f `se1' ///
               "  " %6.3f `b2' "  " %5.3f `se2' "  " %5.3f `r2_1'
    }
}

di as result _n "Interpretation:"
di as result "  β > 0 → CA moves toward surplus (forced deleveraging, Aguiar-Gopinath)"
di as result "  β_def > β_nd → harder adjustment in default episodes"

* ── Save IRF datasets ────────────────────────────────────────────────────

foreach spec in all nd def {
    preserve
        clear
        set obs 5
        gen horizon = _n - 1
        foreach m in b lo90 hi90 {
            svmat `m'_`spec', names(`m')
            rename `m'1 `m'
        }
        gen spec = "`spec'"
        save "$clean/irf_mech_ca_`spec'.dta", replace
    restore
}

* ── Figure 13c: All episodes ─────────────────────────────────────────────

local c_all "23 55 94"

use "$clean/irf_mech_ca_all.dta", clear

twoway ///
    (rarea lo90 hi90 horizon, ///
        color("`c_all'%20") lwidth(none)) ///
    (connected b horizon, ///
        lcolor("`c_all'") lwidth(medthick) msymbol(circle) mcolor("`c_all'")), ///
    yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%5.2f) labsize(medsmall)) ///
    xtitle("Years after onset", size(small)) ///
    ytitle("Cumulative change in current account/GDP (pp)", size(small)) ///
    title("Current Account Response — All Episodes", size(medium) color(navy)) ///
    subtitle("Forced deleveraging test (Aguiar-Gopinath 2006)", size(small)) ///
    legend(off) ///
    note("All 61 spread crisis episodes. Positive = CA moves toward surplus." ///
         "Controls: l1_gdpg l2_gdpg debt. DK SE. Country & year FE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig13c_ca_all.pdf", replace
di as result "Figure saved: fig13c_ca_all.pdf"

* ── Figure 13d: Non-default vs. Default split ────────────────────────────

local c_nd  "34 139 34"
local c_def "157 36 73"

use "$clean/irf_mech_ca_nd.dta", clear
append using "$clean/irf_mech_ca_def.dta"

twoway ///
    (rarea lo90 hi90 horizon if spec=="nd", ///
        color("`c_nd'%20") lwidth(none)) ///
    (connected b horizon if spec=="nd", ///
        lcolor("`c_nd'") lwidth(medthick) msymbol(triangle) mcolor("`c_nd'")) ///
    (rarea lo90 hi90 horizon if spec=="def", ///
        color("`c_def'%20") lwidth(none)) ///
    (connected b horizon if spec=="def", ///
        lcolor("`c_def'") lwidth(medthick) lpattern(dash) ///
        msymbol(square) mcolor("`c_def'")), ///
    yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%5.2f) labsize(medsmall)) ///
    xtitle("Years after onset", size(small)) ///
    ytitle("Cumulative change in current account/GDP (pp)", size(small)) ///
    title("Current Account Response — By Resolution Type", size(medium) color(navy)) ///
    subtitle("Forced deleveraging test (Aguiar-Gopinath 2006)", size(small)) ///
    legend(order(2 "Non-default (39 episodes)" 4 "Default-linked (22 episodes)") ///
           ring(0) pos(7) cols(1) size(small) region(lcolor(none) fcolor(none))) ///
    note("Green = non-default. Red = default-linked." ///
         "Controls: l1_gdpg l2_gdpg debt. DK SE. Country & year FE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig13d_ca_split.pdf", replace
di as result "Figure saved: fig13d_ca_split.pdf"

* ══════════════════════════════════════════════════════════════════════════
* TEST 3b — IPW-CORRECTED CURRENT ACCOUNT LP (nd vs. default)
*   Reweights non-default episodes to match default-linked on pre-crisis
*   fundamentals. If IPW and OLS results are similar → CA finding is not
*   driven by compositional differences in pre-crisis conditions.
*
*   Weights: probit of onset_def ~ debt ca (among onset_all==1)
*            same parsimonious first stage as in 08_ipw_lp.do Act 2
*   Estimator: areg [aw=ipw2] (xtscc does not accept pweights)
*
*   Output: fig13e_ca_ipw.pdf
* ══════════════════════════════════════════════════════════════════════════

* ── Reload panel ─────────────────────────────────────────────────────────
use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* Regenerate CA outcome variables
capture drop ca_base
gen ca_base = L.ca
forvalues h = 0/4 {
    capture drop ch_ca_`h'
    gen ch_ca_`h' = F`h'.ca - ca_base
}

* ── Reconstruct Act 2 IPW weights ────────────────────────────────────────
di as result _n "=== TEST 3b: Reconstructing Act 2 IPW weights (probit) ==="

probit onset_def debt ca if onset_all == 1, vce(robust)
di as result "McFadden Pseudo-R2 (first stage): " e(r2_p)

predict pscore2_ca if onset_all == 1, pr
label var pscore2_ca "Pr(default-linked | crisis onset, X) — for CA test"

* Trim [0.05, 0.95] consistent with 08_ipw_lp.do
gen trimmed2_ca = (pscore2_ca < 0.05 | pscore2_ca > 0.95) if !missing(pscore2_ca)
quietly count if trimmed2_ca == 1
di as result "Observations trimmed: " r(N)
replace pscore2_ca = . if trimmed2_ca == 1

quietly summarize onset_def if onset_all == 1
local p_def = r(mean)
local p_nd  = 1 - `p_def'

gen ipw2_ca = .
replace ipw2_ca = `p_def' / pscore2_ca         if onset_def == 1 & !missing(pscore2_ca)
replace ipw2_ca = `p_nd'  / (1 - pscore2_ca)   if onset_nd  == 1 & !missing(pscore2_ca)
label var ipw2_ca "Act 2 stabilized IPW weight (CA test)"

summarize ipw2_ca if onset_all == 1, detail

* Set weight = 1 for non-onset years (control observations in LP)
replace ipw2_ca = 1 if onset_all == 0 & sample == 1

* ── Run CA LP: OLS and IPW for nd vs. def split ──────────────────────────
di as result _n "========================================================"
di as result "TEST 3b: IPW-CORRECTED CA LP — Non-default vs. Default"
di as result "========================================================"
di as result "h   b_nd_OLS  b_def_OLS  b_nd_IPW  b_def_IPW"

foreach m in b_nd_ols lo90_nd_ols hi90_nd_ols b_def_ols lo90_def_ols hi90_def_ols ///
             b_nd_ipw lo90_nd_ipw hi90_nd_ipw b_def_ipw lo90_def_ipw hi90_def_ipw {
    matrix `m' = J(5,1,.)
}

forvalues h = 0/4 {
    local row = `h' + 1

    * OLS baseline (areg for comparability with IPW)
    capture areg ch_ca_`h' onset_nd onset_def ///
        l1_gdpg l2_gdpg debt L.ca ///
        i.year if sample == 1, absorb(cid) vce(cluster cid)
    if _rc == 0 {
        matrix b_nd_ols[`row',1]    = _b[onset_nd]
        matrix lo90_nd_ols[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
        matrix hi90_nd_ols[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
        matrix b_def_ols[`row',1]   = _b[onset_def]
        matrix lo90_def_ols[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
        matrix hi90_def_ols[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
        local b_nd_u  = _b[onset_nd]
        local b_def_u = _b[onset_def]
    }

    * IPW-weighted
    capture areg ch_ca_`h' onset_nd onset_def ///
        l1_gdpg l2_gdpg debt L.ca ///
        [aw=ipw2_ca] if sample == 1 & !missing(ipw2_ca), absorb(cid) vce(cluster cid)
    if _rc == 0 {
        matrix b_nd_ipw[`row',1]    = _b[onset_nd]
        matrix lo90_nd_ipw[`row',1] = _b[onset_nd]  - 1.645*_se[onset_nd]
        matrix hi90_nd_ipw[`row',1] = _b[onset_nd]  + 1.645*_se[onset_nd]
        matrix b_def_ipw[`row',1]   = _b[onset_def]
        matrix lo90_def_ipw[`row',1]= _b[onset_def] - 1.645*_se[onset_def]
        matrix hi90_def_ipw[`row',1]= _b[onset_def] + 1.645*_se[onset_def]
        local b_nd_w  = _b[onset_nd]
        local b_def_w = _b[onset_def]

        di "h=" `h' "  " %7.3f `b_nd_u' "    " %7.3f `b_def_u' ///
               "    " %7.3f `b_nd_w' "    " %7.3f `b_def_w'
    }
}

* ── Save IRF datasets ────────────────────────────────────────────────────

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b_nd b_def lo90_nd hi90_nd lo90_def hi90_def {
        svmat `m'_ols, names(`m')
        rename `m'1 `m'
    }
    gen series = "ols"
    save "$clean/irf_ca_act2_ols.dta", replace
restore

preserve
    clear
    set obs 5
    gen horizon = _n - 1
    foreach m in b_nd b_def lo90_nd hi90_nd lo90_def hi90_def {
        svmat `m'_ipw, names(`m')
        rename `m'1 `m'
    }
    gen series = "ipw"
    save "$clean/irf_ca_act2_ipw.dta", replace
restore

* ── Figure 13e: CA LP — OLS vs. IPW, nd vs. def ─────────────────────────

local c_nd  "34 139 34"
local c_def "157 36 73"

use "$clean/irf_ca_act2_ols.dta", clear
append using "$clean/irf_ca_act2_ipw.dta"

twoway ///
    (rarea lo90_nd hi90_nd horizon if series=="ols", ///
        color("`c_nd'%15") lwidth(none)) ///
    (connected b_nd horizon if series=="ols", ///
        lcolor("`c_nd'") lwidth(medthick) msymbol(triangle) mcolor("`c_nd'")) ///
    (rarea lo90_def hi90_def horizon if series=="ols", ///
        color("`c_def'%15") lwidth(none)) ///
    (connected b_def horizon if series=="ols", ///
        lcolor("`c_def'") lwidth(medthick) msymbol(square) mcolor("`c_def'")) ///
    (connected b_nd horizon if series=="ipw", ///
        lcolor("`c_nd'") lwidth(medthick) lpattern(dash) ///
        msymbol(triangle_hollow) mcolor("`c_nd'")) ///
    (connected b_def horizon if series=="ipw", ///
        lcolor("`c_def'") lwidth(medthick) lpattern(dash) ///
        msymbol(square_hollow) mcolor("`c_def'")), ///
    yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
    xlabel(0(1)4, labsize(medsmall)) ///
    ylabel(, format(%5.2f) labsize(medsmall)) ///
    xtitle("Years after onset", size(small)) ///
    ytitle("Cumulative change in current account/GDP (pp)", size(small)) ///
    title("Current Account Response — IPW Robustness Check", size(medium) color(navy)) ///
    subtitle("Solid = OLS baseline. Dashed = IPW-weighted.", size(small)) ///
    legend(order(2 "Non-default (OLS)" 4 "Default-linked (OLS)" ///
                 5 "Non-default (IPW)" 6 "Default-linked (IPW)") ///
           cols(2) size(small) region(lcolor(none) fcolor(none))) ///
    note("IPW reweights non-default episodes to match default-linked on pre-crisis" ///
         "fundamentals (probit: debt & CA). areg, clustered SE.", size(vsmall)) ///
    graphregion(color(white)) plotregion(color(white))

graph export "$figs/fig13e_ca_ipw.pdf", replace
di as result "Figure saved: fig13e_ca_ipw.pdf"

di as result _n "13_mechanisms.do complete."
