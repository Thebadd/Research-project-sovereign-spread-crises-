/*===========================================================================
  18_TRANSFORMS.DO   —  FROM-SCRATCH REBUILD, STAGE 8  (finalise panel_lp.dta)
  Build all econometric transforms IN CODE from the sourced levels, then save
  the analysis file $clean/panel_lp.dta consumed by 02..16 unchanged.

    gdpg              real GDP growth, % = 100*(ln gdp_real - ln L.gdp_real)
    l1_gdpg, l2_gdpg  first/second lag of gdpg (control momentum, = cx)
    ln_gdppc_base     L.ln_gdppc (LP baseline at t-1)
    dy_0..dy_4        cumulative % change in log real GDPpc, F h vs t-1 (LP outcome)
    dy_m1, dy_m2      pre-trend placebos (single-year growth fully before onset)
    l_spr_mean/max    lagged EMBIG spread (balance table)
    sample            onset + tranquil years, excl. continuation, GDP base present
===========================================================================*/

use "$clean/panel_build.dta", clear
capture xtset cid year
sort cid year

* ── Real GDP growth + its lags (Asonuma gdpg2 analog; controls cx) ──────────
capture drop gdpg l1_gdpg l2_gdpg
gen double gdpg = 100*(ln(gdp_real) - ln(L.gdp_real)) ///
    if gdp_real > 0 & L.gdp_real > 0 & !missing(gdp_real, L.gdp_real)
label var gdpg "Real GDP growth, % (WEO NGDP_R, log-difference)"
gen double l1_gdpg = L.gdpg
gen double l2_gdpg = L2.gdpg
label var l1_gdpg "L1 real GDP growth"
label var l2_gdpg "L2 real GDP growth"

* ── HEADLINE LP outcome: cumulative % change in log TOTAL real GDP ──────────
*   Aligned with Asonuma et al. (their g_h = F h.ln(gdp_real) - L.ln(gdp_real)),
*   i.e. cumulative real-GDP growth relative to t-1 (NOT per capita).
capture drop ln_gdp ln_gdp_base dy_0 dy_1 dy_2 dy_3 dy_4 dy_m1 dy_m2
gen double ln_gdp = ln(gdp_real) if gdp_real > 0 & !missing(gdp_real)
gen double ln_gdp_base = L.ln_gdp
label var ln_gdp_base "log real GDP at t-1 (LP baseline)"
forvalues h = 0/4 {
    gen double dy_`h' = (F`h'.ln_gdp - ln_gdp_base) * 100
    label var dy_`h' "Cum. % change in log real GDP (total): F`h' vs t-1"
}
* pre-trend placebos on the same (total real GDP) series
gen double dy_m1 = (ln_gdp - L2.ln_gdp) * 100
gen double dy_m2 = (L.ln_gdp - L3.ln_gdp) * 100
label var dy_m1 "Pre-trend h=-1: GDP(t) - GDP(t-2)"
label var dy_m2 "Pre-trend h=-2: GDP(t-1) - GDP(t-3)"

* ── ROBUSTNESS outcome: per-capita version (kept as dy_pc_*) ─────────────────
capture drop ln_gdppc_base dy_pc_0 dy_pc_1 dy_pc_2 dy_pc_3 dy_pc_4
gen double ln_gdppc_base = L.ln_gdppc
label var ln_gdppc_base "log real GDPpc at t-1 (per-capita LP baseline)"
forvalues h = 0/4 {
    gen double dy_pc_`h' = (F`h'.ln_gdppc - ln_gdppc_base) * 100
    label var dy_pc_`h' "Cum. % change in log real GDPpc: F`h' vs t-1 (robustness)"
}

* ── Lagged spreads (used in the balance table) ──────────────────────────────
capture drop l_spr_mean l_spr_max
gen double l_spr_mean = L.spr_mean
gen double l_spr_max  = L.spr_max
label var l_spr_mean "L1 EMBIG mean spread (bps)"
label var l_spr_max  "L1 EMBIG max spread (bps)"

* ── Asonuma-aligned COMMON-CORE control set (plain predetermined columns) ───
*   One control set used as the OUTCOME-model controls in every GDP + channel
*   regression, mirroring Asonuma's $convar: GDP momentum, gov spending, openness,
*   bank-credit depth, hyperinflation dummy, banking-crisis dummy + our debt/ca.
*   Plain lagged columns so the SAME set works in xtscc LP and bsample AIPW.
capture drop hyperinf_dummy l_govexp l_open l_credit_bank
gen byte   hyperinf_dummy = (L.infl > 50) if !missing(L.infl)
gen double l_govexp       = L.govexp
gen double l_open         = L.open
gen double l_credit_bank  = L.credit_bank
label var hyperinf_dummy "Hyperinflation dummy (L.infl > 50) — Asonuma"
label var l_govexp       "L1 govt expenditure, % GDP"
label var l_open         "L1 trade openness, % GDP"
label var l_credit_bank  "L1 bank credit to private / GDP (financial depth)"

global ctrl_core "l1_gdpg l2_gdpg debt ca banking_crisis l_govexp l_open l_credit_bank hyperinf_dummy"

* ── Estimation sample ───────────────────────────────────────────────────────
capture drop sample
gen byte sample = (continuation==0) & !missing(ln_gdp_base)
label var sample "Estimation sample (onset + tranquil, excl. continuation, GDP base present)"

* ── Save the analysis file (drop-in replacement consumed by 02..16) ─────────
compress
sort cid year
xtset cid year
save "$clean/panel_lp.dta", replace

* ══════════════════════════════════════════════════════════════════════════
* VALIDATION SUMMARY
* ══════════════════════════════════════════════════════════════════════════
di as result _n "18_transforms.do complete — panel_lp.dta rebuilt from source."
quietly count if sample==1
di as result "  sample rows: `r(N)'"
quietly count if onset_all==1 & sample==1
di as result "  onsets in sample: `r(N)'"
di as result _n "Coverage at onsets (non-missing) for the key analysis variables:"
foreach v in dy_0 l1_gdpg debt ca infl imf credit fdi claims_govt inv govexp pb ///
             banking_crisis reer_chg revenue_gdp open claimsgov_assets ///
             claimpriv_assets vix fedfunds ust10y {
    capture confirm variable `v'
    if !_rc {
        quietly count if onset_all==1 & !missing(`v')
        di as result "    `v': `r(N)' / 61"
    }
    else di as error "    `v': MISSING VARIABLE"
}
