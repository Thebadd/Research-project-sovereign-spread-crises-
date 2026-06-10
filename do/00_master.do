/*===========================================================================
  MASTER DO-FILE
  Sovereign Spread Crises — Output Cost Estimation
  Local Projections (Jordà 2005)

  Structure:
    01_build_panel.do   → import, merge classification, save panel_lp.dta
    02_lp_all.do        → Act 1: LP for ALL 61 spread crisis episodes
    03_lp_resolution.do → Act 2: LP split by default-linked vs. non-default
    04_graphs.do        → publication-quality IRF figures

  Required packages (run once):
    ssc install xtscc
    ssc install coefplot
    ssc install esttab   (part of estout)
    ssc install boottest
===========================================================================*/

clear all
set more off
set scheme s2mono

* ── Paths ──────────────────────────────────────────────────────────────────
global root  "/home/user/Research-project-sovereign-spread-crises-"
global raw   "$root/data/raw"
global clean "$root/data/clean"
global do    "$root/do"
global figs  "$root/output/figures"
global tabs  "$root/output/tables"

* ── Run pipeline ────────────────────────────────────────────────────────────
do "$do/01_build_panel.do"
do "$do/02_lp_all.do"
do "$do/03_lp_resolution.do"
do "$do/04_graphs.do"
do "$do/05_balance_table.do"
do "$do/06_robustness.do"
