/*===========================================================================
  RETIRED — superseded by 02a_descriptive_facts.do
  ---------------------------------------------------------------------------
  This file is no longer called from 00_master.do and its output is not cited
  anywhere in the paper.

  Why: it compares the two resolution groups on the legacy control list
  `gdpg l1_gdpg debt ca infl l_spr_mean imf fdi govexp` -- contemporaneous
  levels of debt, the current account and inflation, plus the retired IMF-
  programme dummy. The paper conditions on $ctrl_core, all measured at t-1, so
  this table described variables the regressions no longer use. That is the
  same defect that had 13_mechanisms.do running on a retired specification.

  The replacement, in 02a_descriptive_facts.do, reports the same comparison on
  the actual common core with tranquil years as a reference column and a t-test
  on the nd-def difference, and exports descriptive_summary.csv.

  Kept rather than deleted: the t-test loop and the export idiom are a working
  template if a balance table on some other variable set is ever wanted. To
  revive it, update the `balvars' local to $ctrl_core and uncomment the call in
  00_master.do.
===========================================================================*/

/*===========================================================================
  05_BALANCE_TABLE.DO
  Pre-crisis balance: Non-Default vs. Default-Linked episodes

  Compares observable characteristics averaged over t-3 to t-1 (3 years
  before onset) between the two groups. A significant difference signals
  selection bias that must be acknowledged.

  Output: balance_table.csv  (import into LaTeX via booktabs)
===========================================================================*/

use "$clean/panel_lp.dta", clear

* ── Keep only onset rows (one obs per episode) ───────────────────────────
keep if onset_all == 1 & !missing(nondefault)

* Variables to compare (pre-crisis averages already in the panel at t)
* We use the controls at the onset year as proxies for pre-crisis conditions
* (debt, ca, gdpg are measured at t; lagged spread is already L1)

local balvars gdpg l1_gdpg debt ca infl l_spr_mean imf fdi govexp

* ── Summary stats by group ───────────────────────────────────────────────
di as result _n "=== BALANCE TABLE: Non-Default vs. Default-Linked ==="
di as result "Comparing pre-crisis observables at onset year"
di ""

* Store results in matrix
local nvars : word count `balvars'
matrix BT = J(`nvars', 5, .)   // cols: mean_nd, sd_nd, mean_def, sd_def, p_diff
local rnames ""

local i = 1
foreach v of local balvars {

    * Non-default group
    quietly summarize `v' if nondefault == 1
    local m_nd = r(mean)
    local s_nd = r(sd)
    local n_nd = r(N)

    * Default-linked group
    quietly summarize `v' if nondefault == 0
    local m_def = r(mean)
    local s_def = r(sd)
    local n_def = r(N)

    * Two-sample t-test (unequal variances)
    quietly ttest `v', by(nondefault) unequal
    local pval = r(p)

    matrix BT[`i', 1] = `m_nd'
    matrix BT[`i', 2] = `s_nd'
    matrix BT[`i', 3] = `m_def'
    matrix BT[`i', 4] = `s_def'
    matrix BT[`i', 5] = `pval'

    di "  `v':  ND = " %6.2f `m_nd' " (" %5.2f `s_nd' ")" ///
       "   DEF = " %6.2f `m_def' " (" %5.2f `s_def' ")" ///
       "   p = " %5.3f `pval' _c
    if `pval' < 0.05 di "  *"
    else if `pval' < 0.10 di "  +"
    else di ""

    local rnames "`rnames' `v'"
    local ++i
}

di ""
di "N non-default: `n_nd'   N default-linked: `n_def'"

* ── Export to CSV ────────────────────────────────────────────────────────
matrix rownames BT = `rnames'
matrix colnames BT = mean_nd sd_nd mean_def sd_def p_diff

* Convert matrix to dataset and save
clear
svmat BT, names(col)
gen varname = ""
local i = 1
foreach v of local balvars {
    replace varname = "`v'" in `i'
    local ++i
}
order varname mean_nd sd_nd mean_def sd_def p_diff

* Add significance stars
gen stars = ""
replace stars = "***" if p_diff < 0.01
replace stars = "**"  if p_diff >= 0.01 & p_diff < 0.05
replace stars = "*"   if p_diff >= 0.05 & p_diff < 0.10

export delimited "$tabs/balance_table.csv", replace
di as result "Balance table saved: $tabs/balance_table.csv"

* ── Regional breakdown of episodes ──────────────────────────────────────
use "$clean/panel_lp.dta", clear
keep if onset_all == 1 & !missing(nondefault)

di as result _n "=== EPISODE COUNT BY REGION ==="
tab region nondefault, row
