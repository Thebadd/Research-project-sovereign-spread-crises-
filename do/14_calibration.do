/*===========================================================================
  14_CALIBRATION.DO
  Calibrate the structural parameters of the theoretical model
  (Theoretical Appendix v3) to the empirical sample.

  Three parameter groups:
    GROUP 1 — externally calibrated from the EM-DSGE literature
    GROUP 2 — data-determined steady-state ratios (tranquil-period means)
    GROUP 3 — stochastic processes estimated from the panel
              (spread persistence rho_s; income process rho_y, sigma_y)

  Output: scalars stored as globals + a calibration table (cal_params.dta)
          consumed by 15_solve_default.do and 16_model_irf.do.

  NOTE: "tranquil" = non-crisis, non-continuation country-years
        (onset_all==0 & continuation==0). This makes the steady state
        reflect normal times, not crisis dynamics.
===========================================================================*/

use "$clean/panel_lp.dta", clear

* ══════════════════════════════════════════════════════════════════════════
* GROUP 1 — EXTERNALLY CALIBRATED (annual frequency, EM literature)
* ══════════════════════════════════════════════════════════════════════════
* Sources: Aguiar-Gopinath (2006), Arellano (2008), Neumeyer-Perri (2005),
*          Cruces-Trebesch (2013), Gelos et al. (2011).

scalar beta   = 0.96      // discount factor (household)
scalar betap  = 0.95      // government discount factor (impatience for default)
scalar sigma  = 2.0       // inverse EIS / CRRA
scalar alpha  = 0.33      // capital share
scalar delta  = 0.10      // depreciation (annual)
scalar Rstar  = 0.04      // world risk-free rate
scalar theta  = 0.62      // recovery rate = 1 - haircut (Cruces-Trebesch ~38%)
scalar mu     = 0.22      // re-entry probability (avg exclusion ~4.5 yrs)
scalar deltaB = 0.22      // fraction of external debt maturing per year
scalar lambda = 10.0      // bank leverage (assets/equity)
scalar xi     = 0.50      // working-capital share of wage bill (Group 3 refines)
scalar gamma  = 0.80      // sovereign-ceiling pass-through (Group 3 refines)
scalar phi    = 0.05      // balance-sheet pass-through to lending rate (Group 3 refines)

* ══════════════════════════════════════════════════════════════════════════
* GROUP 2 — DATA-DETERMINED STEADY-STATE RATIOS (tranquil means)
* ══════════════════════════════════════════════════════════════════════════

gen byte tranquil = (onset_all==0 & continuation==0) & !missing(ln_gdppc)

* Steady-state EMBIG spread (median tranquil spread, bps -> decimal)
quietly summarize spr_mean if tranquil, detail
scalar s_ss = r(p50) / 10000
di as result "s_ss (steady-state spread, decimal) = " s_ss

* Expenditure shares of GDP (tranquil means, percent -> share)
quietly summarize inv    if tranquil
scalar sI = r(mean)/100
quietly summarize govexp if tranquil
scalar sG = r(mean)/100
quietly summarize ca     if tranquil
scalar sCA = r(mean)/100
scalar sC = 1 - sI - sG - sCA           // consumption share (residual)
di as result "Shares: sC=" sC "  sI=" sI "  sG=" sG "  sCA=" sCA

* Tax rate (general govt revenue / GDP)
capture confirm variable revenue_gdp
if _rc == 0 {
    quietly summarize revenue_gdp if tranquil
    scalar tau = r(mean)/100
}
else scalar tau = 0.25
di as result "tau (tax rate) = " tau

* Bank sovereign-bond exposure proxy: claims on govt / GDP (tranquil mean)
* Maps to b^B_ss in the model (banks' steady-state sovereign holdings).
capture confirm variable claims_govt
if _rc == 0 {
    quietly summarize claims_govt if tranquil
    scalar bB_gdp = r(mean)/100
}
else scalar bB_gdp = 0.15
di as result "bB_gdp (bank claims on govt / GDP) = " bB_gdp

* Steady-state external debt / GDP (tranquil mean of debt)
quietly summarize debt if tranquil
scalar B_gdp = r(mean)/100
di as result "B_gdp (public debt / GDP) = " B_gdp

* Implied steady-state balance-sheet amplification Omega = phi * (bB_ss/N_ss).
* Net worth N_ss = bank assets / lambda; bank assets ~ (bB + credit).
quietly summarize credit if tranquil
scalar credit_gdp = r(mean)/100
scalar bankassets_gdp = bB_gdp + credit_gdp
scalar N_gdp = bankassets_gdp / lambda          // N_ss as share of GDP
scalar Omega = phi * (bB_gdp / N_gdp)
di as result "Omega (= phi * bB/N) = " Omega

* ══════════════════════════════════════════════════════════════════════════
* GROUP 3 — STOCHASTIC PROCESSES
* ══════════════════════════════════════════════════════════════════════════

* ── 3a. Spread persistence rho_s : panel AR(1) of EMBIG spread ────────────
* Use crisis-relevant variation: all country-years with valid spreads.
xtset cid year
quietly xtreg spr_mean L.spr_mean, fe
scalar rho_s = _b[L.spr_mean]
di as result "rho_s (spread persistence) = " rho_s

* ── 3b. Income process rho_y, sigma_y : AR(1) on detrended log GDPpc ───────
* Remove country-specific linear trends, then AR(1) on the cyclical residual.
preserve
    keep cid year ln_gdppc
    drop if missing(ln_gdppc)
    xtset cid year
    * country-specific linear detrend
    gen t = year
    quietly reg ln_gdppc i.cid i.cid#c.t
    predict ycyc, resid
    xtset cid year
    quietly xtreg ycyc L.ycyc, fe
    scalar rho_y = _b[L.ycyc]
    scalar sigma_eps = e(rmse)            // SD of AR(1) innovation
    scalar sigma_y   = sigma_eps / sqrt(1 - rho_y^2)   // unconditional SD
restore
di as result "rho_y = " rho_y "   sigma_y = " sigma_y "   sigma_eps = " sigma_eps

* ══════════════════════════════════════════════════════════════════════════
* DERIVED STEADY-STATE OBJECTS (for the log-linear transmission block)
* ══════════════════════════════════════════════════════════════════════════

* Steady-state rates from SS.1, SS.2, SS.4
scalar RD_ss = 1/beta                          // SS.1
scalar RL_ss = 1/beta + phi*(bB_gdp/N_gdp)     // SS.4
* note: R* + gamma*s_ss should ~ 1/beta - 1 in gross terms; we keep RL_ss as gross-1 premium form

* Working-capital output elasticity eps_p (LL.2 definition)
scalar nexp  = (1-alpha)/alpha
scalar eps_p = xi*nexp*RL_ss / (1 + xi*(RL_ss-1))
di as result "eps_p (working-capital output elasticity) = " eps_p

* Interest semi-elasticity of capital (LL.3 definition)
scalar eta   = 1/((1-alpha)*(RL_ss-(1-delta)))
di as result "eta (capital semi-elasticity) = " eta

* Balance-sheet sensitivity to spread (LL.5a definition)
scalar Bsens = (bB_gdp/N_gdp) / (1+s_ss)^2
di as result "Bsens (B in LL.5a) = " Bsens

* Net-worth multiplier Phi_N — calibrated in Group 3 matching (16_model_irf).
* Placeholder default (refined by SMM in 16): governs speed of credit deepening.
scalar Phi_N = 0.70
scalar Phi_O = 0.30                            // wholesale-funding income share

* ══════════════════════════════════════════════════════════════════════════
* EXPORT CALIBRATION → globals + dataset
* ══════════════════════════════════════════════════════════════════════════

foreach p in beta betap sigma alpha delta Rstar theta mu deltaB lambda ///
             xi gamma phi s_ss sC sI sG sCA tau bB_gdp B_gdp Omega ///
             rho_s rho_y sigma_y sigma_eps RD_ss RL_ss eps_p eta Bsens ///
             Phi_N Phi_O N_gdp credit_gdp {
    global `p' = `p'
}

clear
set obs 1
foreach p in beta betap sigma alpha delta Rstar theta mu deltaB lambda ///
             xi gamma phi s_ss sC sI sG sCA tau bB_gdp B_gdp Omega ///
             rho_s rho_y sigma_y sigma_eps RD_ss RL_ss eps_p eta Bsens ///
             Phi_N Phi_O N_gdp credit_gdp {
    gen `p' = ${`p'}
}
save "$clean/cal_params.dta", replace

di as result _n "════════════════════════════════════════════════════════"
di as result "CALIBRATION COMPLETE — saved to cal_params.dta"
di as result "Group 1 (literature): beta sigma alpha delta R* theta mu deltaB lambda"
di as result "Group 2 (data SS):    s_ss shares tau bB/GDP B/GDP Omega"
di as result "Group 3 (processes):  rho_s rho_y sigma_y"
di as result "Group 3 to be SMM-matched in 16: xi gamma phi Phi_N"
di as result "════════════════════════════════════════════════════════"
