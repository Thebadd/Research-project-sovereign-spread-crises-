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

  SCOPE, PER THE PROJECT OWNER'S EXPLICIT INSTRUCTION: only variables this
  project actually defines AND uses in its own estimation are reported --
  the six headline outcomes (GDP + the five channels the LP/AIPW files
  estimate) and the eight $ctrl_core controls (every regression's own
  control set). Asonuma's Table B3 also reports an "Additional control
  variables" block (debt, terms of trade, Freedom House indices, Paris
  Club/IMF dummies), a "Predictors" block (fed funds, contagion, past
  preemptive cases), an "Additional dependent variables" block (nominal
  lending rate, inflation rate, total credit), and "Additional macro
  controls" (oil prices, fed-funds direction dummies, decade dummies) --
  NONE of those are reproduced here: debt is already inside $ctrl_core (no
  separate row needed), and the rest (tot_chg, l_imf, the first-stage
  probit's own predictors, the raw lending/inflation series real_lending is
  built from, and every "Additional macro controls" item) are either
  robustness-only variables never used in a headline regression, or simply
  not built in this project at all -- reporting them here would describe
  objects this project's own analysis does not actually use.

  CATEGORIES:

    Dependent variables       GDP, investment, bank credit (all log-real,
                               *100), claims on govt/GDP and the two nexus
                               shares (level-change, *100 where a ratio),
                               real lending rate (log level, *100) -- the
                               six outcomes 03_lp_resolution.do/12_channels_
                               resolution.do/08b_aipw.do/13c_aipw_channels.do
                               actually estimate.
    Baseline control variables $ctrl_core (8 terms), reported on their own
                               NATIVE scale exactly as used in every
                               regression project-wide -- no display
                               rescaling (unlike 05b_balance_regressions.do's
                               own "b"-suffix x100 table, which exists
                               specifically to compare against Asonuma's
                               Table H1; this table reports this project's
                               own numbers as they actually are in the code).

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
* NO rescaling here: every variable below is reported on its OWN native
* scale exactly as it sits in this project's code -- dy_0/ch0_* exactly as
* 18_transforms.do/02a_descriptive_facts.do/03c_table_combined_six.do build
* them, and every $ctrl_core term on its own native (decimal) scale exactly
* as used in every regression project-wide. No display-only x100 rescale
* (unlike 05b_balance_regressions.do's own "b"-suffix table, which exists
* specifically to be compared against Asonuma's Table H1 on their scale --
* this table is not that, and reports this project's own numbers as-is).

di as result _n "════════════════════════════════════════════════════════════"
di as result "TABLE [X]: SUMMARY STATISTICS (sample==1, onset + tranquil years, native scale)"
di as result "════════════════════════════════════════════════════════════"

di as result _n "-- Dependent variables --"
summarize dy_0 ch0_inv ch0_credit ch0_claims_govt ch0_claimsgov_assets ch0_claimpriv_assets ch0_real_lending if sample==1

di as result _n "-- Baseline control variables (\$ctrl_core, native scale) --"
summarize l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2 if sample==1

di as result _n "NOT INCLUDED: only the outcomes actually estimated by the LP/AIPW files"
di as result "and the controls in \$ctrl_core (every regression's own control set) are"
di as result "reported -- no additional-control, predictor, or auxiliary-outcome rows"
di as result "for variables not part of the headline design, and no display rescaling."

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT (estpost summarize -> esttab, one block at a time so section
* headers can be preserved via separate title() rows rather than one flat
* table with no grouping)
* ══════════════════════════════════════════════════════════════════════════
local cellspec "count(fmt(0) label(Obs.)) mean(fmt(3) label(Mean)) sd(fmt(3) label(Std. Dev.)) min(fmt(3) label(Min)) max(fmt(3) label(Max))"

estpost summarize dy_0 ch0_inv ch0_credit ch0_claims_govt ch0_claimsgov_assets ch0_claimpriv_assets ch0_real_lending if sample==1
esttab using "$tabs/table_summary_statistics.rtf", replace cells("`cellspec'") ///
    noobs nonumber varwidth(28) msign(none) label ///
    title("Table [X]. Summary statistics -- Dependent variables (h=1, native scale)") ///
    coeflabel(dy_0 "GDP" ch0_inv "Investment" ///
              ch0_credit "Bank credit" ch0_claims_govt "Claims on govt / GDP" ///
              ch0_claimsgov_assets "Bank claims on govt / assets" ch0_claimpriv_assets "Bank claims on private / assets" ///
              ch0_real_lending "Real lending interest rate")

estpost summarize l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2 if sample==1
esttab using "$tabs/table_summary_statistics.rtf", append cells("`cellspec'") ///
    noobs nonumber varwidth(28) msign(none) label ///
    title("Baseline control variables (\$ctrl_core, native scale)") ///
    addnotes("Sample: onset + tranquil years. Every variable is reported on its own native scale, exactly as used in this project's regressions.") ///
    coeflabel(l1_gdpg "GDP growth rate" l_debt "Debt-to-GDP ratio" l_banking_crisis "Banking crisis dummy" ///
              l_govexp "Govt. expenditure-to-GDP ratio" l_open "Openness" ///
              l_credit_bank "Bank credit-to-GDP ratio" l_lninfl "Log inflation" exchange2 "Nominal exchange rate change")

* Plain CSV alongside the RTF -- easier to reformat cleanly in Excel/Word
* than relying on esttab's own RTF rendering.
preserve
    tempname C
    tempfile csvf
    postfile `C' str32 variable long n double mean double sd double min double max using "`csvf'", replace
    foreach v in dy_0 ch0_inv ch0_credit ch0_claims_govt ch0_claimsgov_assets ch0_claimpriv_assets ch0_real_lending ///
                 l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2 {
        quietly summarize `v' if sample==1
        post `C' ("`v'") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max))
    }
    postclose `C'
    use "`csvf'", clear
    export delimited "$tabs/summary_statistics.csv", replace
restore
di as result "Plain CSV also saved: $tabs/summary_statistics.csv"

if _rc == 608 di as error "  ** table_summary_statistics.rtf is OPEN IN WORD -- close it and re-run to refresh."
else if _rc  di as error "  ** Table (summary statistics): esttab failed (rc=" _rc ")"
else di as result _n "Summary statistics table saved: $tabs/table_summary_statistics.rtf"

di as result _n "03d_summary_statistics.do complete."
