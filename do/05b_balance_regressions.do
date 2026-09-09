/*===========================================================================
  05B_BALANCE_REGRESSIONS.DO
  Regression of control variables on onset dummies, styled on Asonuma et al.
  (2024) Table H1 ("Regression of Restructuring Strategies on Control
  Variables") -- the balance/justification table for this project's own
  $ctrl_core, in place of a simple t-test.

  WHAT IT DOES, MATCHING THEIR REPLICATION SCRIPT EXACTLY IN SPIRIT: for each
  control variable used in $ctrl_core, regress that variable ON the onset
  dummies (their dum1/dum2/dum3 -> this project's onset_nd/onset_def) with
  country fixed effects and standard errors clustered by country:

      xtreg <control> onset_nd onset_def if sample==1, fe vce(cluster cid)

  This is a WITHIN-COUNTRY comparison: does a country's own control variable
  move, in the crisis year, relative to that country's own tranquil-year
  average, differently depending on which type of crisis it is? A significant
  coefficient on onset_nd or onset_def signals that the treatment and control
  groups are NOT balanced on that dimension going into the crisis year --
  exactly the selection the first-stage propensity model and the LP/AIPW
  controls are then asked to address, so this table is read as a
  justification for including that control, not as a balance test that has
  to fail to pass.

  WHY THIS IS A DIFFERENT TABLE FROM 08c_first_stage_table.do: 08c asks
  "does the probability of a crisis TYPE depend on these controls" (probit,
  controls as regressors, crisis type as outcome). This file asks the
  reverse question Table H1 asks: "does each control variable ITSELF differ
  in the crisis year vs tranquil years" (linear FE regression, crisis-type
  dummies as regressors, one control at a time as the outcome). Both are
  legitimate balance/justification exercises; this file is the one that
  matches Table H1's specific design.

  SAMPLE: if sample==1, the same onset-tier estimation universe every other
  LP/AIPW file in this project uses (onset + tranquil years, continuation and
  carry-in rows excluded) -- not the full raw panel, so this table describes
  exactly the comparison the outcome regressions are run on, not a different
  population.

  VARIABLES: the 8 terms of $ctrl_core (18_transforms.do), each already at
  its regression-ready construction (t-1, and on whatever scale -- decimal or
  log -- the outcome regressions themselves use):
    l1_gdpg          GDP growth
    l_debt           Debt-to-GDP ratio
    l_banking_crisis Banking crisis dummy
    l_govexp         Government expenditure-to-GDP ratio
    l_open           Trade openness (exports + imports / GDP)
    l_credit_bank    Bank credit-to-GDP ratio
    l_lninfl         Log inflation
    exchange2        Log change in the nominal exchange rate
  No new variables are built here -- all eight already exist on panel_lp.dta.

  Output: $tabs/table_balance_regressions.rtf
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
local balvars $ctrl_core

eststo clear

di as result _n "════════════════════════════════════════════════════════════"
di as result "BALANCE REGRESSIONS -- each control on onset_nd/onset_def, country FE, cluster(cid)"
di as result "════════════════════════════════════════════════════════════"

foreach v of local balvars {
    capture confirm variable `v'
    if _rc {
        di as error "  ** `v' not found -- skipping."
        continue
    }
    quietly xtreg `v' onset_nd onset_def if sample==1, fe vce(cluster cid)
    eststo col_`v'
    di as result %-18s "`v'" "  coef(onset_nd)=" %8.3f _b[onset_nd] "  coef(onset_def)=" %8.3f _b[onset_def] ///
        "  N=" %5.0f e(N) "  N_g=" %3.0f e(N_g) "  r2=" %5.3f e(r2)
}

esttab col_l1_gdpg col_l_debt col_l_banking_crisis col_l_govexp col_l_open ///
       col_l_credit_bank col_l_lninfl col_exchange2 ///
    using "$tabs/table_balance_regressions.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(onset_nd onset_def) order(onset_nd onset_def) ///
    coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
    mtitles("GDP growth" "Debt/GDP" "Banking crisis dummy" "Govt. exp./GDP" ///
             "Openness" "Bank credit/GDP" "Log inflation" "Nominal exchange rate") ///
    nonumber ///
    stats(N N_g r2, labels("Observations" "Countries" "R-squared") fmt(0 0 3)) ///
    title("Table [X]. Regression of onset type on control variables") ///
    addnotes("Dependent variable: the indicated control (measured at t-1, matching \$ctrl_core), regressed on the onset dummies with country fixed" ///
             "effects. Sample: onset + tranquil years (sample==1), the same estimation universe as every LP/AIPW outcome regression in this project." ///
             "Standard errors clustered by country in parentheses. * p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_balance_regressions.rtf is OPEN IN WORD -- close it and re-run to refresh."
else if _rc  di as error "  ** Table (balance regressions): esttab failed (rc=" _rc ")"
else di as result _n "Balance regressions table saved: $tabs/table_balance_regressions.rtf"

di as result _n "05b_balance_regressions.do complete."
