/*===========================================================================
  03C_TABLE_COMBINED_SIX.DO
  Single merged onset-tier OLS table, all SIX of the project's headline
  outcome variables (GDP, Bank credit, Investment, Bank claims on
  government/GDP, FDI, Real lending rate) as panels in ONE RTF -- Act 2
  (nd/def joint regression), matching 03_lp_resolution.do's own headline
  GDP spec (table2_output_resolution.rtf) and 12_channels_resolution.do's
  own channel spec (table4_channels_resolution.rtf) EXACTLY, panel for
  panel, with no change to either specification.

  WHY A SEPARATE FILE, NOT JUST ESTTAB-ING ACROSS 03 AND 12: confirmed
  this session that 12_channels_resolution.do's own `eststo clear' (its
  own line 121) wipes out 03_lp_resolution.do's t2_h* estimates when both
  run in the same Stata session (00_master.do runs 03 then 12) -- there is
  also no existing project precedent for `estimates save'/`estimates use'
  across do-files (grepped, none found). Every existing multi-panel table
  in this project re-estimates fresh within ONE do-file's own eststo/
  esttab block, so this file follows that same convention: it re-runs
  BOTH the GDP regression and the five channel regressions itself, in one
  file, so the eststo store is never cleared between them.

  SPEC (identical to both source files, unchanged here):
    xtreg <outcome> onset_nd onset_def <controls> if sample==1[&
    common_abcd==1], fe vce(robust) -- country FE only, no year FE, plain
    robust SE, tranquil years the omitted category, matching Asonuma et
    al.'s own Table I1 design. Difference tested via the Wald F-statistic
    (test onset_nd = onset_def), the covariance-correct test for two
    coefficients from the same joint regression -- not a Clogg z or
    bootstrap (see 03_lp_resolution.do's own header for the full argument).
    GDP + credit/inv/claims_govt run on the balanced A-D sample
    (common_abcd==1); FDI/real_lending keep their own best-available
    sample, unrestricted -- exactly as in both source files.

  STATS ROW STANDARDIZED ACROSS ALL SIX PANELS to Table 4's own set
  (bdiff, fdiff, pdiff, nepnd, nepdef, N, N_g) -- GDP's own source table
  additionally reports per-arm country counts (nctynd/nctydef); those are
  dropped here for a single consistent stats block across the whole
  table, not because they are wrong, just to standardize on the set five
  of the six panels already use. GDP's own fdiff (Wald F-statistic) is
  captured the same way Table 4 already does (r(F) right after `test'),
  which 03_lp_resolution.do's own table2_output_resolution.rtf does NOT
  currently report explicitly (it shows only the p-value) -- added here
  for consistency across panels, not present in the original GDP table.

  Output: $tabs/table_combined_six_resolution.rtf.
  Run AFTER 18_transforms.do. Not wired into 00_master.do. Does not modify
  03_lp_resolution.do, 12_channels_resolution.do, or either of their own
  RTF outputs -- both keep producing their existing tables unchanged.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
local controls $ctrl_core

* ── Helper: episode/country counts inside e(sample), matching
* 03_lp_resolution.do's own _nepcount (used for the GDP panel only, to
* stay faithful to that file's own methodology; the channel panels use
* the simpler in-sample count 12_channels_resolution.do itself uses --
* preserved as-is so each panel's numbers match its own source table
* exactly, not homogenized across the two methodologies).
capture program drop _nepcount
program define _nepcount, rclass
    syntax varname(numeric) , Outcome(varname) Controls(varlist)
    tempvar esmp
    quietly gen byte `esmp' = e(sample)
    quietly count if `esmp' == 1
    if r(N) == 0 {
        quietly replace `esmp' = (sample == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar n = r(N)
end

eststo clear

* ══════════════════════════════════════════════════════════════════════════
* PANEL A — GDP (mirrors 03_lp_resolution.do's own headline joint block)
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== PANEL A: GDP ==="
forvalues h = 0/4 {
    quietly xtreg dy_`h' onset_nd onset_def `controls' if sample==1 & common_abcd==1, fe vce(robust)

    local bnd  = _b[onset_nd]
    local bdef = _b[onset_def]
    local bdiff = `bdef' - `bnd'

    quietly test onset_nd = onset_def
    local pd = r(p)
    local pf = r(F)

    _nepcount onset_nd,  outcome(dy_`h') controls(`controls')
    local nepnd = r(n)
    _nepcount onset_def, outcome(dy_`h') controls(`controls')
    local nepdef = r(n)

    eststo t2_h`h', title("h=`=`h'+1'")
    estadd scalar bdiff  = `bdiff'
    estadd scalar fdiff  = `pf'
    estadd scalar pdiff  = `pd'
    estadd scalar nepnd  = `nepnd'
    estadd scalar nepdef = `nepdef'

    di "h=" `h'+1 ":  beta_nd=" %6.3f `bnd' "  beta_def=" %6.3f `bdef' ///
       "  diff=" %6.3f `bdiff' "  F=" %6.2f `pf' "  p=" %5.3f `pd'
}

* ══════════════════════════════════════════════════════════════════════════
* PANELS B-F — channels (mirrors 12_channels_resolution.do's own loop)
* ══════════════════════════════════════════════════════════════════════════
local channels credit inv claims_govt fdi real_lending

local ctrl_credit       l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
local ctrl_inv          $ctrl_core pre_inv
local ctrl_claims_govt  $ctrl_core pre_claims_govt
local ctrl_fdi          $ctrl_core pre_fdi
local ctrl_real_lending $ctrl_core pre_real_lending

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

    di as result _n "=== PANEL: `ch' ==="
    forvalues h = 0/4 {
        capture xtreg ch_`ch'_`h' onset_nd onset_def `ctrl' ///
            if sample == 1`balflag', fe vce(robust)

        if _rc == 0 {
            local bnd  = _b[onset_nd]
            local bdef = _b[onset_def]
            local bdiff = `bdef' - `bnd'

            test onset_nd = onset_def
            local pd = r(p)
            local pf = r(F)

            quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepnd = r(N)
            quietly count if onset_def == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepdef = r(N)

            eststo t4_`ch'_`h', title("h=`=`h'+1'")
            estadd scalar bdiff  = `bdiff'
            estadd scalar fdiff  = `pf'
            estadd scalar pdiff  = `pd'
            estadd scalar nepnd  = `nepnd'
            estadd scalar nepdef = `nepdef'
            local elist_`ch' `elist_`ch'' t4_`ch'_`h'

            di "h=" `h'+1 ":  beta_nd=" %6.3f `bnd' "  beta_def=" %6.3f `bdef' ///
               "  diff=" %6.3f `bdiff' "  F=" %6.2f `pf' "  p=" %5.3f `pd'
        }
        else di as error "  ** regression failed for `ch' h=" `h'+1 " (rc=" _rc ")"
    }
}

* ══════════════════════════════════════════════════════════════════════════
* SINGLE MERGED ESTTAB — GDP panel first (replace), each channel appended
* ══════════════════════════════════════════════════════════════════════════
local statslabels labels("Difference (default - non-default)" "  F (Wald, nd = def)" ///
                          "  p (Wald, nd = def)" "Episodes (non-default)" ///
                          "Episodes (default)" "Observations" "Countries")
local statsfmt fmt(3 2 3 0 0 0 0)

capture esttab t2_h0 t2_h1 t2_h2 t2_h3 t2_h4 using "$tabs/table_combined_six_resolution.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(onset_nd onset_def) order(onset_nd onset_def) ///
    coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
    mtitles nonumber ///
    stats(bdiff fdiff pdiff nepnd nepdef N N_g, `statslabels' `statsfmt') ///
    title("Table X. Output cost by crisis resolution -- Panel A: GDP") ///
    addnotes("Dependent variable: cumulative change in the indicated outcome (pp) from t-1 to t+h. Both onset dummies enter jointly with tranquil years as the omitted category, matching the reference paper's baseline. Jorda (2005) local projections; country fixed effects only (no year FE); common-core controls plus (for channel panels) the channel's own pre-crisis change. Robust (heteroskedasticity-only) standard errors in parentheses. F(nd=def) and p(nd=def) are the Wald equality test. GDP, Bank credit, and Investment are estimated on the balanced A-D sample (common_abcd); Claims on government/GDP is also balanced; FDI and Real lending rate use their own best-available sample. * p<0.10, ** p<0.05, *** p<0.01.")

local ptitle_credit       "Panel B: Bank credit"
local ptitle_inv          "Panel C: Investment"
local ptitle_claims_govt  "Panel D: Bank claims on government / GDP"
local ptitle_fdi          "Panel E: FDI"
local ptitle_real_lending "Panel F: Real lending rate"

foreach ch of local channels {
    if "`elist_`ch''" == "" {
        di as error "  ** table_combined_six_resolution.rtf: no estimates for `ch' -- panel skipped"
        continue
    }
    capture esttab `elist_`ch'' using "$tabs/table_combined_six_resolution.rtf", append ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_nd onset_def) order(onset_nd onset_def) ///
        coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
        mtitles nonumber ///
        stats(bdiff fdiff pdiff nepnd nepdef N N_g, `statslabels' `statsfmt') ///
        title("`ptitle_`ch''")
}

if _rc == 608 di as error "  ** table_combined_six_resolution.rtf is OPEN IN WORD -- close it and re-run."
else if _rc di as error "  ** esttab failed on the last panel (rc=" _rc ")"
else di as result "Combined six-variable table saved: $tabs/table_combined_six_resolution.rtf"

di as result _n "03c_table_combined_six.do complete."
