/*===========================================================================
  03D_SUMMARY_STATISTICS.DO
  Summary statistics table, styled on Asonuma et al. (2024) Table B3 --
  adapted to this project's own outcome/control set (no fabricated rows for
  variables this project does not build: Freedom House indices, Paris Club
  dummy, decade dummies, oil prices, and fed-funds-rise/decline dummies are
  NOT in this project's data and are omitted rather than invented -- see the
  NOT INCLUDED note at the bottom of this header).

  SAMPLE: matches the AIPW propensity-weighted regression's own restriction
  in spirit -- Asonuma's script keeps only rows used in AT LEAST ONE of the
  three per-type IPW regressions for GDP (samp1+samp2+samp3>0). This project
  has no analogous per-type sample union to reuse (its onset-tier estimation
  universe is a single `sample' flag, not per-type IPW samples), so summary
  stats here are computed on `sample==1' directly -- the onset + tranquil-
  year estimation universe every LP/AIPW file in this project already uses.
  This is stated explicitly rather than left to look like a literal match to
  their sample-construction step.

  HORIZON: Asonuma's table reports h=1 (the crisis/restructuring year itself,
  their g_v_1 = var - L.var). This project's own Year-1/h=0 convention is the
  direct analog (dy_0, ch_v_0 = F0.src - L.src = src - L.src) -- same
  horizon, this project's own naming.

  CATEGORIES, mapped from Table B3 to what this project actually has:

    Dependent variables       GDP, investment, bank credit (all log-real,
                               *100), claims on govt/GDP and the two nexus
                               shares (level-change, *100 where a ratio),
                               real lending rate (log level, *100).
    Baseline control variables $ctrl_core (8 terms) -- reported on the SAME
                               "b"-suffix x100 scale as 05b_balance_
                               regressions.do, for readability and for direct
                               comparability with Asonuma's own baseline
                               block (their gdpg2/gov_exp2/open2/credit_bank2
                               are also x100-style percent scales, not raw
                               decimals).
    Additional control variables debt/GDP (l_debt, x100), terms of trade
                               (tot_chg), IMF-supported program dummy (l_imf).
                               Freedom House indices and the Paris Club dummy
                               are NOT built in this project (no source data)
                               and are omitted, not fabricated.
    Predictors                 Federal funds rate, contagion (default-linked,
                               distance-weighted), years since the most
                               recent prior default-linked onset (this
                               project's recency-clock analog of Asonuma's
                               "number of past preemptive cases" count).
    Additional dependent vars  Nominal lending rate change, inflation rate
                               change (both *100) -- the two series real
                               lending rate is built from. This project has
                               no separate "total credit" series distinct
                               from bank credit and does not report that row.
    Additional macro controls  NOT built in this project (oil prices, fed
                               funds rise/decline dummies, decade dummies) --
                               omitted, not fabricated.

  Output: $tabs/table_summary_statistics.rtf (via `estpost summarize`+esttab)
          and console echo of every `summarize' call.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

* ── Dependent variables: channel outcomes at h=0 (Year 1 = crisis year),
* identical construction to 02a_descriptive_facts.do / 03c_table_combined_
* six.do (ch_v_h = F h.src - L.src; log-real-level src for credit/inv, raw
* level for claims_govt/claimsgov_assets/claimpriv_assets/real_lending). Not
* persisted in panel_lp.dta, built fresh here. ───────────────────────────
foreach v in credit inv claims_govt claimsgov_assets claimpriv_assets real_lending {
    local src `v'
    if inlist("`v'","credit","inv") local src ln_r_`v'
    capture drop `v'_sbase
    quietly gen double `v'_sbase = L.`src'
    capture drop ch0_`v'
    quietly gen double ch0_`v' = F0.`src' - `v'_sbase
}

* GDP: dy_0 already exists on panel_lp.dta (18_transforms.do), no rebuild.

* ── Rescale to the Asonuma-comparable x100 display scale ──────────────────
* dy_0/ch0_credit/ch0_inv/ch0_real_lending are already log-differenced (need
* no further *100 -- they are already on that scale from their own
* construction). claims_govt/claimsgov_assets/claimpriv_assets are ratio-
* point differences already in percentage-point units (no *100 needed
* either, matching how 02a's own dd_* outputs are read directly). The
* $ctrl_core "b"-suffix rescale below mirrors 05b_balance_regressions.do
* exactly.
local rescale100 l1_gdpg l_debt l_govexp l_open l_credit_bank l_lninfl exchange2
foreach v of local rescale100 {
    capture confirm variable `v'
    if !_rc {
        capture drop `v'_b
        quietly gen double `v'_b = `v' * 100
    }
}

* Nominal lending rate / inflation change at h=0 (same construction as
* 02a_descriptive_facts.do's lending/infl_defl block), *100 not needed since
* these are already plain percentage-point series (lending/infl_defl are
* WDI rates in percent already).
capture confirm variable lending
capture confirm variable infl_defl
if !_rc {
    foreach v in lending infl_defl {
        capture drop `v'_sbase
        quietly gen double `v'_sbase = L.`v'
        capture drop ch0_`v'
        quietly gen double ch0_`v' = F0.`v' - `v'_sbase
    }
}

di as result _n "════════════════════════════════════════════════════════════"
di as result "TABLE [X]: SUMMARY STATISTICS (sample==1, onset + tranquil years)"
di as result "════════════════════════════════════════════════════════════"

di as result _n "-- Dependent variables --"
summarize dy_0 ch0_inv ch0_credit ch0_claims_govt ch0_claimsgov_assets ch0_claimpriv_assets ch0_real_lending if sample==1

di as result _n "-- Baseline control variables ($ctrl_core, x100 scale) --"
summarize l1_gdpg_b l_govexp_b l_open_b l_banking_crisis l_credit_bank_b l_lninfl_b exchange2_b if sample==1

di as result _n "-- Additional control variables --"
summarize l_debt_b tot_chg l_imf if sample==1

di as result _n "-- Predictors --"
summarize l_fedfunds l_contagion_dist_def years_since_def_onset if sample==1

di as result _n "-- Additional dependent variables --"
capture confirm variable ch0_lending
if !_rc summarize ch0_lending ch0_infl_defl if sample==1
else di as error "  ** lending/infl_defl not built -- add lendinginterestrate.xlsx/inflationgdpdeflator.xlsx to \$raw."

di as result _n "NOT INCLUDED (no source data in this project, not fabricated):"
di as result "  Freedom House civil-liberties/political-rights indices, Paris Club dummy,"
di as result "  oil prices, fed-funds rise/decline dummies, decade dummies, a separate"
di as result "  'total credit' series distinct from bank credit."

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT (estpost summarize -> esttab, one block at a time so section
* headers can be preserved via separate `title()' rows rather than one flat
* table with no grouping)
* ══════════════════════════════════════════════════════════════════════════
estpost summarize dy_0 ch0_inv ch0_credit ch0_claims_govt ch0_claimsgov_assets ch0_claimpriv_assets ch0_real_lending if sample==1
esttab using "$tabs/table_summary_statistics.rtf", replace cells("count(fmt(0) label(Obs.)) mean(fmt(2) label(Mean)) sd(fmt(2) label(Std. Dev.)) min(fmt(2) label(Min)) max(fmt(2) label(Max))") ///
    noobs nonumber title("Table [X]. Summary statistics -- Dependent variables") ///
    coeflabel(dy_0 "GDP (log real, x100 change, h=1)" ch0_inv "Investment (log real, x100 change, h=1)" ///
              ch0_credit "Bank credit (log real, x100 change, h=1)" ch0_claims_govt "Claims on government / GDP (pp change, h=1)" ///
              ch0_claimsgov_assets "Bank claims on govt / assets (pp change, h=1)" ch0_claimpriv_assets "Bank claims on private / assets (pp change, h=1)" ///
              ch0_real_lending "Real lending interest rate (log level, x100 change, h=1)")

estpost summarize l1_gdpg_b l_govexp_b l_open_b l_banking_crisis l_credit_bank_b l_lninfl_b exchange2_b if sample==1
esttab using "$tabs/table_summary_statistics.rtf", append cells("count(fmt(0) label(Obs.)) mean(fmt(2) label(Mean)) sd(fmt(2) label(Std. Dev.)) min(fmt(2) label(Min)) max(fmt(2) label(Max))") ///
    noobs nonumber title("Baseline control variables") ///
    coeflabel(l1_gdpg_b "GDP growth rate" l_govexp_b "Government expenditure-to-GDP ratio" l_open_b "Openness" ///
              l_banking_crisis "Banking crisis dummy" l_credit_bank_b "Bank credit-to-GDP ratio" ///
              l_lninfl_b "Log inflation" exchange2_b "Nominal exchange rate change")

estpost summarize l_debt_b tot_chg l_imf if sample==1
esttab using "$tabs/table_summary_statistics.rtf", append cells("count(fmt(0) label(Obs.)) mean(fmt(2) label(Mean)) sd(fmt(2) label(Std. Dev.)) min(fmt(2) label(Min)) max(fmt(2) label(Max))") ///
    noobs nonumber title("Additional control variables") ///
    coeflabel(l_debt_b "Debt-to-GDP ratio" tot_chg "Terms of trade (rate of change)" l_imf "IMF-supported program dummy")

estpost summarize l_fedfunds l_contagion_dist_def years_since_def_onset if sample==1
esttab using "$tabs/table_summary_statistics.rtf", append cells("count(fmt(0) label(Obs.)) mean(fmt(2) label(Mean)) sd(fmt(2) label(Std. Dev.)) min(fmt(2) label(Min)) max(fmt(2) label(Max))") ///
    noobs nonumber title("Predictors") ///
    coeflabel(l_fedfunds "Federal funds rate" l_contagion_dist_def "Contagion, based on default-linked crises" ///
              years_since_def_onset "Years since last default-linked onset")

capture confirm variable ch0_lending
if !_rc {
    estpost summarize ch0_lending ch0_infl_defl if sample==1
    esttab using "$tabs/table_summary_statistics.rtf", append cells("count(fmt(0) label(Obs.)) mean(fmt(2) label(Mean)) sd(fmt(2) label(Std. Dev.)) min(fmt(2) label(Min)) max(fmt(2) label(Max))") ///
        noobs nonumber title("Additional dependent variables") ///
        addnotes("Freedom House indices, Paris Club dummy, oil prices, fed-funds rise/decline dummies, decade dummies, and a separate" ///
                 "total-credit series are not built in this project (no source data) and are not reported." ///
                 "Sample: onset + tranquil years (sample==1). GDP/investment/bank credit/real lending rate are log-differenced, x100." ///
                 "Claims-on-government and the two bank-asset shares are percentage-point changes. Baseline/additional controls marked" ///
                 "'_b' or listed as debt/govexp/openness/credit_bank/log inflation/exchange rate are shown x100, matching Asonuma et al." ///
                 "(2024)'s own Table B3 scale; the banking crisis and IMF-program dummies are unscaled (0/1).") ///
        coeflabel(ch0_lending "Nominal lending rate (pp change, h=1)" ch0_infl_defl "Inflation rate, GDP deflator (pp change, h=1)")
}

if _rc == 608 di as error "  ** table_summary_statistics.rtf is OPEN IN WORD -- close it and re-run to refresh."
else if _rc  di as error "  ** Table (summary statistics): esttab failed (rc=" _rc ")"
else di as result _n "Summary statistics table saved: $tabs/table_summary_statistics.rtf"

di as result _n "03d_summary_statistics.do complete."
