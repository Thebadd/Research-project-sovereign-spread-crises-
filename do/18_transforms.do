/*===========================================================================
  18_TRANSFORMS.DO   —  FROM-SCRATCH REBUILD, STAGE 8  (finalise panel_lp.dta)
  Build all econometric transforms IN CODE from the sourced levels, then save
  the analysis file $clean/panel_lp.dta consumed by 02..16 unchanged.

    gdpg              real GDP growth, % = 100*(ln gdp_real - ln L.gdp_real)
    l1_gdpg, l2_gdpg  first/second lag of gdpg (control momentum, = cx)
    ln_gdp_base       L.ln_gdp (headline LP baseline at t-1; ln_gdppc_base = per-capita robustness)
    dy_0..dy_4        cumulative % change in log TOTAL real GDP (ln gdp_real), F h vs t-1 (headline outcome)
    dy_pc_0..dy_pc_4  per-capita version (robustness)
    dy_m1, dy_m2      pre-trend placebos (two-year growth spans ending before onset)
    hyperinf_dummy (L.infl>50, core), l_lninfl (= ln(1+L.infl/100)) and l_infl (raw)
    built as continuous robustness alternatives, not in the core,
    l_govexp, l_open, l_credit (+l_credit_bank robustness), $ctrl_core  common-core ingredients
    tot_chg, exchange2, ex_dum1-5  robustness-tier controls (Asonuma additional)
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
* Guard: these sources are built upstream — infl/govexp in stage 11 (WEO),
* open/credit_bank in stage 12 (WDI). If any is absent, panel_build.dta was not
* built by the full 10->18 chain (e.g. 18 was run alone on a stale file, or a
* WDI raw file was missing). Stop with a clear message instead of a cryptic halt.
foreach req in infl govexp open credit_bank {
    capture confirm variable `req'
    if _rc {
        di as error "  ** 18_transforms: required variable '`req'' is not in panel_build.dta."
        di as error "     Built upstream: infl/govexp -> stage 11 (11_weo); open/credit_bank -> stage 12 (12_wdi)."
        di as error "     Fix: run the FULL chain in order (do 00_master.do), with all data/raw files present."
        exit 111
    }
}
capture drop hyperinf_dummy l_infl l_lninfl l_govexp l_open l_credit_bank l_credit l_debt l_ca l_banking_crisis
* INFLATION IN THE CORE = hyperinf_dummy (David's decision), matching Asonuma et al.
* Raw inflation cannot be used directly here: across the 52 panel countries 1994-2026 the
* median is 9.0% while Venezuela 2018 is 65,374% (9 obs above 1000%), so a raw term would
* let a handful of country-years set the coefficient. The dummy sidesteps that entirely by
* flagging the hyperinflation regime instead of its magnitude — 63 country-years exceed the
* 50% threshold (Venezuela 16, Turkey 11, Bulgaria/Argentina/Lebanon 4 each, 17 countries
* in all), so the flag is well populated in the raw data.
*   NB it collapsed to a single estimation-sample observation in the pre-2026 runs, but
*   that was the OLD banking control truncating the panel at 2018, which killed the entire
*   recent cluster (Argentina 2019-24, Lebanon 2020-23, Turkey 2022-24, Venezuela 2019-26).
*   With the Laeven-Valencia rebuild extending the sample past 2018 those years return.
*   Worth re-checking the probit for "predicts failure perfectly" on the new sample.
gen byte   hyperinf_dummy = (L.infl > 50) if !missing(L.infl)
* Continuous alternatives, built but NOT in the core — available for a dummy-vs-continuous
* robustness spec, and a one-token swap in the 14 $ctrl_core definitions if wanted. Use
* l_lninfl rather than l_infl for that: the log tames the tail (Venezuela ~75x the median
* instead of ~7,000x). Safe because the minimum is -8.53% (Azerbaijan 1999), well above
* the -100% where ln(1+x) breaks.
gen double l_infl         = L.infl
gen double l_lninfl       = ln(1 + L.infl/100) if !missing(L.infl) & L.infl > -100
gen double l_govexp       = L.govexp
gen double l_open         = L.open
gen double l_credit_bank  = L.credit_bank
gen double l_credit       = L.credit
* debt, current account, and the banking-crisis flag lagged to t-1 so the ENTIRE
* common core is predetermined (t-1) — internally coherent and matching the
* reference paper's pre-crisis controls (avoids onset-year simultaneity/bad-control).
gen double l_debt         = L.debt
gen double l_ca           = L.ca
gen byte   l_banking_crisis = L.banking_crisis
* Global-push predictors lagged to t-1 (predetermined), matching the reference
* paper's federal_funds2 = L.federal_funds. Plain saved columns => bootstrap-safe
* in the AIPW (no L. operator at estimation time). Rates are pure time series, so
* L. under xtset cid year gives the year-(t-1) value for every country. vix/ust10y
* built too so a future push-variable swap is one line. Raw fedfunds/vix/ust10y
* columns are kept (coverage report below still uses them).
gen double l_fedfunds     = L.fedfunds
gen double l_vix          = L.vix
gen double l_ust10y       = L.ust10y
label var l_fedfunds "L1 US fed funds rate (predetermined; = Asonuma federal_funds2)"
label var l_vix      "L1 CBOE VIX (predetermined)"
label var l_ust10y   "L1 US 10y Treasury yield (predetermined)"
label var hyperinf_dummy "Hyperinflation dummy (L.infl > 50) — common core"
label var l_infl         "L1 CPI inflation, % (raw; robustness alternative, not in core)"
label var l_lninfl       "L1 log gross inflation = ln(1+L.infl/100) (robustness alt., not in core)"
label var l_govexp       "L1 govt expenditure, % GDP"
label var l_open         "L1 trade openness, % GDP"
label var l_credit_bank  "L1 bank credit to private / GDP (financial depth, by-banks; robustness)"
label var l_credit       "L1 private credit / GDP (financial depth, total; common core)"
label var l_debt         "L1 public debt, % GDP (predetermined)"
label var l_ca           "L1 current account, % GDP (predetermined)"
label var l_banking_crisis      "L1 systemic banking-crisis dummy (predetermined)"

global ctrl_core "l1_gdpg l_debt l_ca l_banking_crisis l_govexp l_open l_credit hyperinf_dummy"

* ── ROBUSTNESS-tier controls (Asonuma additional controls; NOT in the core) ──
*   terms-of-trade change and nominal-FX-change quantile dummies, built the
*   Asonuma way. Predetermined. Used only in robustness specs.
capture drop tot_chg exchange2 ex_dum1 ex_dum2 ex_dum3 ex_dum4 ex_dum5
* terms of trade (optional robustness source) — skip if not present
capture confirm variable tot
if !_rc {
    gen double tot_chg = 100*(ln(L.tot) - ln(L2.tot)) if L.tot>0 & L2.tot>0 & !missing(L.tot,L2.tot)
    label var tot_chg "L1 terms-of-trade log-change, % (WDI; robustness)"
}
else di as error "  ** tot not found — skipping tot_chg (robustness only; add termsoftrade.xlsx to build it)."
* official FX (optional robustness source) — build exchange2 + ex_dum bins, else skip
capture confirm variable exch
if !_rc {
    gen double exchange2 = ln(1+L.exch) - ln(1+L2.exch) if !missing(L.exch,L2.exch)
    label var exchange2 "L1 nominal exchange-rate log-change (robustness)"
    quietly summarize exchange2 if continuation==0 & !missing(ln_gdp_base), detail
    foreach k in 1 2 3 4 5 {
        gen byte ex_dum`k' = 0 if !missing(exchange2)
    }
    replace ex_dum1 = 1 if exchange2 < r(p5)
    replace ex_dum2 = 1 if exchange2 < r(p25)
    replace ex_dum2 = 0 if exchange2 < r(p5)
    replace ex_dum3 = 1 if exchange2 < r(p50)
    replace ex_dum3 = 0 if exchange2 < r(p25)
    replace ex_dum4 = 1 if exchange2 < r(p75)
    replace ex_dum4 = 0 if exchange2 < r(p50)
    replace ex_dum5 = 1 if exchange2 < r(p95)
    replace ex_dum5 = 0 if exchange2 < r(p75)
    label var ex_dum1 "FX-change bin < p5 (Asonuma ex_dum; robustness)"
}
else di as error "  ** exch not found — skipping exchange2/ex_dum (robustness only; add officialexchangerate.xlsx)."

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
foreach v in dy_0 l1_gdpg debt ca infl hyperinf_dummy l_lninfl imf credit fdi claims_govt inv govexp pb ///
             banking_crisis l_banking_crisis reer_chg revenue_gdp open claimsgov_assets ///
             claimpriv_assets vix fedfunds ust10y {
    capture confirm variable `v'
    if !_rc {
        quietly count if onset_all==1 & !missing(`v')
        di as result "    `v': `r(N)' / 61"
    }
    else di as error "    `v': MISSING VARIABLE"
}

* ── Inflation control: confirm the log transform tamed the tail ─────────────
* The core carries hyperinf_dummy. These lines report how many observations trip the
* 50% threshold and what the continuous alternatives look like, so the choice stays
* visible in the run log.
capture confirm variable l_lninfl
if !_rc {
    di as result _n "Inflation control (core uses hyperinf_dummy = L.infl > 50):"
    quietly count if hyperinf_dummy==1 & sample==1
    di as result "    hyperinf_dummy==1 in estimation sample: `r(N)'"
    di as result "       (was 1 pre-rebuild, when the panel was truncated at 2018)"
    quietly summarize l_infl if sample==1, detail
    di as result "    raw  l_infl   median " %8.2f r(p50) "   max " %11.1f r(max)
    quietly summarize l_lninfl if sample==1, detail
    di as result "    log  l_lninfl median " %8.3f r(p50) "   max " %11.3f r(max)   "  (alt., not in core)"
}
