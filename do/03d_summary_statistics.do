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
    Predictors                 The pooled (Act 1) propensity-model predictor
                               set `cz' actually used project-wide (l_fedfunds
                               l_reg_crisis_share past_onsets) -- ADDED per
                               the project owner's follow-up request, which
                               revises the earlier deliberate omission above:
                               these are first-stage predictors, not outcome-
                               model variables, but they are genuinely built
                               and used (08b_aipw.do/13c/13d's own `cz'
                               local), so reporting them is accurate, not
                               invented. The default-linked-specific
                               predictor set (`cz_def': l_fedfunds
                               l_contagion_dist_atdef years_since_def_onset)
                               is NOT added here -- it differs from `cz' only
                               in the contagion and recency terms, both
                               already documented in 17_predictors.do, and a
                               second near-duplicate block would add little.

  Output: $tabs/table_summary_statistics.xlsx (native Excel, via `export
          excel') AND $tabs/table_summary_statistics.rtf (a real bordered
          Word table, matching the reference paper's Table B3 visual layout:
          one block per category, a bold block-header row, then one row per
          variable with Obs./Mean/Std. Dev./Min/Max), plus a console echo of
          every `summarize' call.
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

di as result _n "-- Predictors (pooled Act 1 propensity set \`cz', native scale) --"
summarize l_fedfunds l_reg_crisis_share past_onsets if sample==1

di as result _n "NOT INCLUDED: only the outcomes actually estimated by the LP/AIPW files,"
di as result "the controls in \$ctrl_core, and the pooled first-stage predictor set \`cz'"
di as result "are reported -- no additional-control, additional-dependent-variable, or"
di as result "auxiliary-macro-control rows for variables not part of this project's own"
di as result "headline design, and no display rescaling."

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT -- native Excel (.xlsx, one sheet per block) AND a real
* bordered RTF/Word table (one file, three blocks stacked, matching the
* reference paper's Table B3 visual layout: bold block-header row, then
* Variable | Obs. | Mean | Std. Dev. | Min | Max per row).
* ══════════════════════════════════════════════════════════════════════════
local depvars    dy_0 ch0_inv ch0_credit ch0_claims_govt ch0_claimsgov_assets ch0_claimpriv_assets ch0_real_lending
local depvarlab  `" "GDP" "Investment" "Bank credit" "Claims on govt / GDP" "Bank claims on govt / assets" "Bank claims on private / assets" "Real lending interest rate" "'
local ctrlvars   l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2
local ctrlvarlab `" "GDP growth rate" "Debt-to-GDP ratio" "Banking crisis dummy" "Govt. expenditure-to-GDP ratio" "Openness" "Bank credit-to-GDP ratio" "Log inflation" "Nominal exchange rate change" "'
local predvars   l_fedfunds l_reg_crisis_share past_onsets
local predvarlab `" "Federal funds rate" "Regional crisis share (contagion)" "Number of past onsets" "'

tempname S
tempfile sumf
postfile `S' str48 variable long obs double mean double sd double min double max byte blockn using "`sumf'", replace

local i = 1
foreach v of local depvars {
    local lab : word `i' of `depvarlab'
    quietly summarize `v' if sample==1
    post `S' ("`lab'") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max)) (1)
    local ++i
}
local i = 1
foreach v of local ctrlvars {
    local lab : word `i' of `ctrlvarlab'
    quietly summarize `v' if sample==1
    post `S' ("`lab'") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max)) (2)
    local ++i
}
local i = 1
foreach v of local predvars {
    local lab : word `i' of `predvarlab'
    quietly summarize `v' if sample==1
    post `S' ("`lab'") (r(N)) (r(mean)) (r(sd)) (r(min)) (r(max)) (3)
    local ++i
}
postclose `S'

preserve
    use "`sumf'", clear
    label var variable "Variable"
    label var obs      "Obs."
    label var mean      "Mean"
    label var sd        "Std. Dev."
    label var min        "Min"
    label var max        "Max"

    keep if blockn==1
    drop blockn
    export excel "$tabs/table_summary_statistics.xlsx", replace firstrow(varlabels) sheet("Dependent variables")
restore
preserve
    use "`sumf'", clear
    keep if blockn==2
    drop blockn
    label var variable "Variable"
    label var obs      "Obs."
    label var mean      "Mean"
    label var sd        "Std. Dev."
    label var min        "Min"
    label var max        "Max"
    export excel "$tabs/table_summary_statistics.xlsx", sheetreplace firstrow(varlabels) sheet("Baseline controls")
restore
preserve
    use "`sumf'", clear
    keep if blockn==3
    drop blockn
    label var variable "Variable"
    label var obs      "Obs."
    label var mean      "Mean"
    label var sd        "Std. Dev."
    label var min        "Min"
    label var max        "Max"
    export excel "$tabs/table_summary_statistics.xlsx", sheetreplace firstrow(varlabels) sheet("Predictors")
restore

* ── RTF/Word table, real bordered rows (6 cells: label + Obs/Mean/SD/Min/Max) ──
capture program drop _rtfrow6b
program define _rtfrow6b
    syntax , [C1(string) C2(string) C3(string) C4(string) C5(string) C6(string) BOLD SHADE]
    local bd "\clbrdrt\brdrs\brdrw10\clbrdrl\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\clbrdrr\brdrs\brdrw10"
    local sh = cond("`shade'"=="shade", "\clshdng10000\clcbpat22", "")
    local rowdef "\trowd\trgaph80\trleft-108`sh'`bd'\cellx3600`sh'`bd'\cellx4700`sh'`bd'\cellx5800`sh'`bd'\cellx6900`sh'`bd'\cellx8000`sh'`bd'\cellx9100"
    local os = cond("`bold'"=="bold", "\b ", "")
    local oe = cond("`bold'"=="bold", "\b0 ", "")
    file write t3s "`rowdef'" _n
    file write t3s "\pard\intbl\ql {`os'`c1'`oe'}\cell \qc {`os'`c2'`oe'}\cell {`os'`c3'`oe'}\cell {`os'`c4'`oe'}\cell {`os'`c5'`oe'}\cell {`os'`c6'`oe'}\cell " _n
    file write t3s "\row" _n
end

capture file close t3s
file open t3s using "$tabs/table_summary_statistics.rtf", write replace
file write t3s "{\rtf1\ansi\deff0" _n
file write t3s "{\b Table [X]. Summary Statistics\par}" _n
file write t3s "\par" _n
_rtfrow6b, c1("Variable") c2("Obs.") c3("Mean") c4("Std. Dev.") c5("Min") c6("Max") bold shade

foreach blk in 1 2 3 {
    if `blk'==1 local blktit "Dependent variables"
    if `blk'==2 local blktit "Baseline control variables"
    if `blk'==3 local blktit "Predictors"
    _rtfrow6b, c1(`"`blktit'"') bold

    preserve
        use "`sumf'", clear
        keep if blockn==`blk'
        local nr = _N
        forvalues r = 1/`nr' {
            local vv  = variable[`r']
            local oo : display %6.0f obs[`r']
            local mm : display %8.3f mean[`r']
            local ss : display %8.3f sd[`r']
            local mn : display %8.2f min[`r']
            local mx : display %8.2f max[`r']
            _rtfrow6b, c1(`"`vv'"') c2(`"`oo'"') c3(`"`mm'"') c4(`"`ss'"') c5(`"`mn'"') c6(`"`mx'"')
        }
    restore
}

file write t3s "\pard\par" _n
file write t3s "{\i Notes: The table shows summary statistics of variables for the observations included in this project's estimation sample" _n
file write t3s " (sample==1: onset + tranquil years). All variables are on their own native scale, exactly as used in the regressions --" _n
file write t3s " no display rescaling. Dependent variables are reported at h=1 (the crisis onset year, ch0\_v = v - L.v). Baseline control" _n
file write t3s " variables are the eight terms in the project's common core (\$ctrl\_core), used as outcome-model controls in every LP/AIPW" _n
file write t3s " regression. Predictors are the pooled (Act 1) propensity-model predictor set (\`cz').\par}" _n
file write t3s "}" _n
file close t3s

di as result _n "Summary statistics tables saved:"
di as result "  $tabs/table_summary_statistics.xlsx (3 sheets: Dependent variables / Baseline controls / Predictors)"
di as result "  $tabs/table_summary_statistics.rtf  (real bordered Word table, 3 blocks stacked)"
di as result "  Sample: onset + tranquil years (sample==1). Every variable is on its own native scale, exactly as used in this project's regressions."

di as result _n "03d_summary_statistics.do complete."
