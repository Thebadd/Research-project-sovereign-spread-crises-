/*===========================================================================
  20B_EXPOSURE_HETEROGENEITY_FLOW.DO
  Pre-crisis exposure x FLOW-treatment interactions — the flow counterpart
  of 13b_exposure_heterogeneity.do.

  WHY THIS FILE EXISTS
  ---------------------
  13b tests whether a spread crisis hits output MORE where a country was
  MORE EXPOSED, pre-crisis, to a given channel -- by interacting standardized
  pre-crisis exposure with crisis ONSET (one row per episode). This file
  re-runs the identical idea on the FLOW treatment (in_crisis/in_crisis_nd/
  in_crisis_def -- every year of an episode, onset and continuation alike),
  the same relationship every other onset/flow file pair in this project
  already covers (02->20, 11->22, 12->23, 13c->24, 13d->25). 13b itself is
  the one onset-tier design that has not had a flow counterpart until now.

  THE SAME 9 EXPOSURES, THE SAME SIGN-AMBIGUITY CLASSIFICATION -- CARRIED
  OVER VERBATIM, NOT RE-DERIVED
  -----------------------------------------------------------------------
  13b's header states a hard-won lesson: an interaction's sign is evidence of
  channel TRANSMISSION only for an exposure with a single economically
  sensible sign. Most of these 9 do not have one:

    AMBIGUOUS (ratio proxies financial development as much as vulnerability)
      credit, claims_govt, inv, fdi, debt_service
      -- more bank-dependent / capital-intensive  => amplifies   (d < 0)
      -- deeper, more developed financial system  => cushions    (d > 0)
      13b found these come back POSITIVE (the development reading dominates)
      and reports them as UNINFORMATIVE about transmission, not as evidence
      the channel protects.

    UNAMBIGUOUS (only one sign is economically sensible)
      stdebt_share      short-term debt / total external debt -- rollover
                        risk. No development reading: a country is not
                        "more developed" for having more debt maturing next
                        year. THE mechanism the user asked about directly
                        ("harder to refinance").
      claimsgov_assets  bank claims on govt / bank assets -- a PORTFOLIO
                        SHARE, not a size/depth measure, so it does not scale
                        with development the way a GDP ratio does. The
                        sovereign-bank nexus parameter (also tested, more
                        thoroughly, in 13d/25's median-split design).

    intpay_gni (interest payments on external debt / GNI -- the OTHER
    mechanism the user asked about, "higher debt payments") is collinear
    with stdebt_share and debt_service -- 13b's own note: "stdebt_share,
    debt_service and intpay_gni all proxy the same external-debt
    refinancing-burden construct... read them as ALTERNATIVE measures
    (robustness), not independent channels." Tested here on the same terms
    13b tests it, not singled out as a separate finding.

    reserves_extdebt is merged by 01d_merge_vulnerability.do but EXCLUDED
    here, same as 13b ("gave a noisy, wrong-signed default estimate on a
    tiny default cell").

  This classification is about what the EXPOSURE VARIABLE measures, not
  about how the treatment is coded, so flow coding does not change it. Read
  the sign only where the classification above says it is interpretable.

  WHAT IS NEW: EVERY EXPOSURE MUST BE ENTRY-DATED, NOT JUST ROW-DATED
  ---------------------------------------------------------------------------
  13b measures exposure at t-1 (`gen exp_e = L.e`), predetermined relative to
  a single onset row. Under flow coding that is not enough: a CONTINUATION
  row's own L.exposure is partly a crisis OUTCOME by the time an episode is
  several years in (bank credit ratios, claims shares, debt-service ratios
  all move during a crisis) -- exactly the contamination problem
  18_transforms.do's epc_* construction and 25_aipw_nexus_split_flow.do's
  a_nexus both exist to solve. So every exposure here is entry-dated: the
  episode's own onset-year value, held fixed across every continuation row
  of that episode; tranquil rows keep their own row-dated value. Identical
  mechanism to a_nexus, applied to all 9 exposures instead of one.

  STANDARDIZATION, CONFOUND CHECK, CONTROLS, INFERENCE
  ------------------------------------------------------
  Standardized to mean-0/SD-1 over `sample_flow==1` (the flow-sample analogue
  of 13b's own `sample==1`) -- same population 13b standardizes over, just
  swapped to the flow estimation sample. The income-sorting confound check
  (does a median split on the exposure just sort onsets by development?)
  stays at `onset_all==1 & sample==1`, unchanged from 13b: it is a one-time
  cross-sectional check on the onset population and is not something later
  crisis-years' coding affects (numerically identical to 13b's own check,
  since the entry-dated exposure equals the row-dated value at the onset row
  itself).

  Controls: $ctrl_flow (episode-dated core; flow_ctrl_variant toggle, same
  as 20/22/23/24/25), country + year FE (drop_year_fe toggle, same
  convention), Driscoll-Kraay SE at lag(max(2, h+3)) -- the FLOW rule (not
  13b's onset rule max(1,h+1)): the regressor here is itself flow-coded and
  serially correlated within an episode, same reasoning as every other
  flow-tier LP in this project. Sample: sample_flow==1 (tranquil + every
  crisis-year, onset and continuation).

  PART A (pooled) IS KEPT, NOT DROPPED -- UNLIKE 21/24/25
  ---------------------------------------------------------------------------
  The AIPW files (21/24/25) dropped their pooled treatment-vs-tranquil line
  because of a problem specific to the propensity model's near-tautology on
  continuation rows. This file is plain OLS, not AIPW -- 20_lp_flow.do (this
  file's direct sibling) keeps a full pooled section ("1. POOLED FLOW LP")
  alongside its by-type section, and so does this file, mirroring 13b's own
  two-part structure (Part A pooled, Part B by resolution type) exactly.

  PART B'S DIFFERENCE TEST NEEDS NO NEW MACHINERY
  ---------------------------------------------------------------------------
  `test Dnd_e = Ddef_e` is a plain within-regression Wald test on two
  interaction terms from the SAME joint regression -- the covariance between
  them is already correctly accounted for. This is nothing like the AIPW
  files' bootstrap-vs-Clogg-z tension (two SEPARATELY estimated AIPW cells);
  13b already uses exactly this test and it needs no adaptation here.

  MULTIPLE TESTING -- SAME CAVEAT, SAME ARITHMETIC
  ---------------------------------------------------------------------------
  9 exposures x 5 horizons x 2 panels = 90 tests; at 5% about 4-5 significant
  results are expected by chance. No single panel here should be read as
  decisive. Use this file for the two UNAMBIGUOUS exposures above, and read
  any nexus-flavoured result as one leg of the same three-design pattern
  13b's own header describes (13d/25 median split, 13b/this file's continuous
  interaction, 13c/24 AIPW channel gap).

  Outputs
  -------
    "$tabs/table5f_exposure_interactions_flow.rtf"    Word: pooled (Part A)
    "$tabs/table6f_exposure_by_resolution_flow.rtf"    Word: by type (Part B)
    "$tabs/exposure_interactions_flow.csv"             raw pooled coefficients
    "$tabs/exposure_by_resolution_flow.csv"            raw by-type coefficients
    "$figs/fig13f_exposure_interactions_flow.pdf"      pooled d_h, grid
    "$figs/fig13g_exposure_by_resolution_flow.pdf"     d_nd vs d_def, grid
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow sample ///
             onset_all onset_nd onset_def ep_seq gdppc_real {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel_lp.dta — re-run 18_transforms.do first."
        exit 111
    }
}
if "$ctrl_flow" == "" {
    di as error "  ** \$ctrl_flow not set — re-run 18_transforms.do, then this file."
    exit 111
}

* ── Exposure variables and readable labels -- IDENTICAL to 13b's own list ──
* Default-mechanism exposures (access loss / balance-sheet / capital):
*   credit, claims_govt, claimsgov_assets, inv, fdi
* Non-default-mechanism exposures (cost / rollover) -- what the user asked
* about directly:
*   claimpriv_assets (bank private-lending share -> lending-rate pass-through)
*   stdebt_share     (short-term debt / ext. debt -> rollover risk)      [01d]
*   debt_service     (debt service / exports      -> liquidity burden)   [01d]
*   intpay_gni       (interest payments / GNI      -> avg interest cost)  [01d]
* stdebt_share, debt_service and intpay_gni all proxy the same external-debt
* refinancing-burden construct and are collinear -- read them as ALTERNATIVE
* measures (robustness), not independent channels. Same note as 13b.
* (reserves_extdebt is merged by 01d but excluded here, same as 13b: it gave
*  a noisy, wrong-signed default estimate on a tiny default cell there.)
* The 01d (IDS) variables are skipped automatically if not merged yet.
local expvars    credit claims_govt claimsgov_assets claimpriv_assets inv fdi ///
                 stdebt_share debt_service intpay_gni
local lbl_credit           "Private credit/GDP"
local lbl_claims_govt      "Bank claims on govt/GDP"
local lbl_claimsgov_assets "Bank sovereign exposure (nexus)"
local lbl_claimpriv_assets "Bank private lending share"
local lbl_inv              "Investment/GDP"
local lbl_fdi              "FDI/GDP"
local lbl_stdebt_share     "Short-term debt share (rollover)"
local lbl_debt_service     "Debt service (% exports)"
local lbl_intpay_gni       "Interest burden avg (% GNI)"

* Keep only exposures that actually exist in the panel; skip (with a warning)
* any that have not been merged yet, so this file always runs -- same
* graceful-skip pattern as 13b.
local expok
foreach e of local expvars {
    capture confirm variable `e', exact
    if _rc == 0 local expok `expok' `e'
    else di as error "  ** exposure `e' not in panel yet — skipped (run 01d to merge it)"
}
local expvars `expok'

* CONTROL SET. flow_ctrl_variant: 0 = the ADOPTED flow-tier baseline,
* $ctrl_flow, built in 18_transforms.do from $ctrl_core
* (l_banking_duration -> the l_banking_crisis DUMMY, l_ca -> tot_chg -> exchange2,
* l_hyperinfl -> l_lninfl -- see that file's "ADOPTED FLOW-TIER CORE CONTROL
* SET"). 1 = the adopted core PLUS the reference paper's own additional
* predictors, tot_chg (the term exchange2 replaced) and l_imf (18_transforms.do's "ALTERNATE FLOW
* CONTROL SET"). Default 0.
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
if `flow_ctrl_variant'==1 local ctrl_flow_base $ctrl_flow_flowplus
else                       local ctrl_flow_base $ctrl_flow

* EXPLORATORY: set to 1 to drop year FE and match the reference paper's
* single-stage rule (country FE only), matching 20/22/23/24/25's toggle of
* the same name. Default 0 = current baseline.
local drop_year_fe 0
local yearfe = cond(`drop_year_fe', "", "i.year")

local controls `ctrl_flow_base'

* ══════════════════════════════════════════════════════════════════════════
* 1. ENTRY-DATED, STANDARDIZED EXPOSURES AND INTERACTIONS
*
* Entry-dating: the episode's own onset-year value, held fixed across every
* continuation row of that episode. Tranquil rows keep their own row-dated
* value. Identical construction to 18_transforms.do's epc_* terms and
* 25_aipw_nexus_split_flow.do's a_nexus, applied to all 9 exposures here.
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== PRE-CRISIS EXPOSURE COVERAGE (entry-dated, flow sample) ==="
foreach e of local expvars {
    capture drop exp_`e'_row _ent_exp_`e' exp_`e' z_`e' Dz_`e' Dnd_`e' Ddef_`e'

    * Row-dated pre-crisis level, as 13b builds it.
    gen exp_`e'_row = L.`e'

    * Entry-dated: freeze at the episode's onset-year value for every row of
    * that episode (onset + continuation); tranquil rows keep their own
    * row-dated value.
    quietly bysort cid ep_seq: egen double _ent_exp_`e' = ///
        max(cond(onset_all==1, exp_`e'_row, .))
    gen double exp_`e' = cond(in_crisis==1, _ent_exp_`e', exp_`e'_row)
    drop _ent_exp_`e'
    * bysort...egen re-sorts the physical dataset by its by-list, leaving it
    * sorted by cid ep_seq rather than cid year -- xtscc needs cid year.
    sort cid year

    * Standardize over the FLOW estimation sample (mean 0, SD 1) -- the flow-
    * sample analogue of 13b's own `if sample==1'.
    quietly summarize exp_`e' if sample_flow == 1
    gen z_`e' = (exp_`e' - r(mean)) / r(sd)
    label var z_`e' "Std. pre-crisis `lbl_`e'' (entry-dated)"

    * Interactions: flow treatment (pooled and by-type) x standardized exposure
    gen Dz_`e'   = in_crisis     * z_`e'
    gen Dnd_`e'  = in_crisis_nd  * z_`e'
    gen Ddef_`e' = in_crisis_def * z_`e'
    label var Dz_`e'   "Flow-crisis x `lbl_`e''"
    label var Dnd_`e'  "Flow-crisis(non-default) x `lbl_`e''"
    label var Ddef_`e' "Flow-crisis(default) x `lbl_`e''"

    quietly count if sample_flow==1 & in_crisis==1     & !missing(z_`e')
    local n_all = r(N)
    quietly count if sample_flow==1 & in_crisis_nd==1  & !missing(z_`e')
    local n_nd  = r(N)
    quietly count if sample_flow==1 & in_crisis_def==1 & !missing(z_`e')
    local n_def = r(N)
    di as result "  `e': all=" `n_all' "  non-default=" `n_nd' "  default=" `n_def'
}

* ══════════════════════════════════════════════════════════════════════════
* CONFOUND CHECK — same as 13b's, unchanged: does each exposure sort onsets
* on DEVELOPMENT rather than on vulnerability to the channel? Computed at
* the ONSET population only (onset_all==1 & sample==1) -- a one-time cross-
* sectional check, unaffected by how later crisis-years are coded, and
* numerically identical to 13b's own check since the entry-dated exposure
* equals the row-dated value at the onset row itself.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== CONFOUND CHECK: does each exposure sort onsets by income? ==="
di as result "  exposure            mean GDPpc low-Z   high-Z     ratio  n_lo/n_hi"
foreach e of local expvars {
    quietly summarize z_`e' if sample==1 & onset_all==1, detail
    local zmed = r(p50)
    capture drop _hiZ
    quietly gen byte _hiZ = (z_`e' >= `zmed') if !missing(z_`e')
    quietly summarize gdppc_real if sample==1 & onset_all==1 & _hiZ==0
    local glo = r(mean)
    local nlo = r(N)
    quietly summarize gdppc_real if sample==1 & onset_all==1 & _hiZ==1
    local ghi = r(mean)
    local nhi = r(N)
    local rat = .
    if `glo' > 0 & !missing(`glo') & !missing(`ghi') local rat = `ghi'/`glo'
    local flag ""
    if !missing(`rat') & (`rat' > 1.5 | `rat' < 0.667) local flag "   ** sorts on income"
    di as result "  " %-18s "`e'" "  " %12.0f `glo' "  " %10.0f `ghi' ///
       "  " %6.2f `rat' "   " %2.0f `nlo' "/" %2.0f `nhi' "`flag'"
}
capture drop _hiZ
di as result "  (ratio = mean real GDP per capita, high-Z bin / low-Z bin, over onsets.)"
di as result "  Flagged rows sort onsets by income: for those the interaction d cannot be"
di as result "  read as channel exposure. Same check, same population, as 13b's own."

sort cid year   // belt-and-suspenders before the estimation loops

* ══════════════════════════════════════════════════════════════════════════
* PART A — POOLED FLOW INTERACTION (in_crisis x exposure)
*   Coefficient of interest: Dz_`e' (= d_h). Negative => channel amplifies
*   the output loss where pre-crisis exposure was higher. Read the sign only
*   where the classification in the header says it is interpretable.
* ══════════════════════════════════════════════════════════════════════════

foreach e of local expvars {
    foreach m in dcoef dlo90 dhi90 {
        matrix `m'_`e' = J(6, 1, 0)
    }
}

eststo clear

foreach e of local expvars {

    di as result _n "========================================================"
    di as result "PART A — EXPOSURE: `lbl_`e'' (pooled, flow-coded)"
    di as result "========================================================"
    di as result "h    b_crisis  SE       d_interact  SE       p(d)"

    local elistA_`e'

    forvalues h = 0/4 {
        local lag = max(2, `h'+3)
        local row = `h' + 2

        capture xtscc dy_`h' in_crisis z_`e' Dz_`e' `controls' `yearfe' ///
            if sample_flow == 1, fe lag(`lag')

        if _rc == 0 {
            eststo a_`e'_`h', title("h=`=`h'+1'")
            local elistA_`e' `elistA_`e'' a_`e'_`h'

            matrix dcoef_`e'[`row',1] = _b[Dz_`e']
            matrix dlo90_`e'[`row',1] = _b[Dz_`e'] - 1.645*_se[Dz_`e']
            matrix dhi90_`e'[`row',1] = _b[Dz_`e'] + 1.645*_se[Dz_`e']

            local pd = 2*(1 - normal(abs(_b[Dz_`e']/_se[Dz_`e'])))
            di "h=" `h'+1 "   " %7.3f _b[in_crisis] "  " %6.3f _se[in_crisis] ///
               "   " %7.3f _b[Dz_`e'] "  " %6.3f _se[Dz_`e'] ///
               "   " %5.3f `pd'
        }
        else di as error "h=" `h'+1 ": xtscc failed for exposure `e' (rc=" _rc ")"
    }
}

di as result _n "Interpretation (Part A) -- READ THE SIGN ONLY WHERE IT IS IDENTIFIED:"
di as result "  d_h < 0 => output falls MORE where pre-crisis exposure was higher."
di as result "  credit/inv/fdi/claims_govt/debt_service are GDP ratios that proxy"
di as result "  financial DEPTH as much as vulnerability -- a positive d there is NOT"
di as result "  evidence the channel protects, and NOT evidence against the channel."
di as result "  stdebt_share and claimsgov_assets are the two exposures with a single"
di as result "  economically sensible sign; intpay_gni/debt_service are collinear"
di as result "  ALTERNATIVE measures of the same refinancing-burden construct."

* ── TABLE 5f: pooled flow exposure interactions (Word/RTF) ─────────────────
local t5fnote "Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h, FLOW coding (in_crisis = 1 in every year of an episode, onset and continuation alike). Each column adds one channel's standardized, ENTRY-DATED pre-crisis exposure and its interaction with the flow crisis dummy. 'Flow-crisis x exposure' is the effect per 1 SD of pre-crisis exposure; a negative value means the output loss is deeper where exposure was higher. Exposure entry-dated (frozen at the episode's onset-year value across continuation rows) to avoid contamination by the crisis itself. Country and year fixed effects; Driscoll-Kraay SE at lag(max(2,h+3)), the flow-tier rule. * p<0.10, ** p<0.05, *** p<0.01. A negative coefficient is evidence of channel transmission ONLY for an exposure with a single economically sensible sign -- see file header. This file runs 9 exposures x 5 horizons x 2 panels = 90 tests; at the 5% level roughly 4-5 significant cells are expected by chance. No single panel should be read as decisive."

local panel A
local writemode replace
local t5ffail 0

foreach e of local expvars {
    if "`elistA_`e''" == "" {
        di as error "  ** Table 5f: no estimates for exposure `e' — panel skipped"
        local t5ffail 1
        continue
    }
    local t5fextra
    if "`e'" == "intpay_gni" local t5fextra addnotes("`t5fnote'")

    capture esttab `elistA_`e'' using "$tabs/table5f_exposure_interactions_flow.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(in_crisis Dz_`e') order(in_crisis Dz_`e') ///
        coeflabel(in_crisis "Flow spread crisis (at mean exposure)" ///
                  Dz_`e' "Flow-crisis x pre-crisis exposure (per SD)") ///
        mtitles nonumber ///
        stats(N N_g, labels("Observations" "Countries") fmt(0 0)) ///
        title("Table 5f. Exposure heterogeneity, FLOW coding (pooled) -- Panel `panel': `lbl_`e''") ///
        `t5fextra'

    if _rc == 608 {
        di as error "  ** table5f_exposure_interactions_flow.rtf is OPEN IN WORD — close it and re-run."
        local t5ffail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 5f: esttab failed for exposure `e' (rc=" _rc ")"
        local t5ffail 1
        continue
    }
    local writemode append
    if      "`panel'"=="A" local panel B
    else if "`panel'"=="B" local panel C
    else if "`panel'"=="C" local panel D
    else if "`panel'"=="D" local panel E
    else if "`panel'"=="E" local panel F
    else if "`panel'"=="F" local panel G
    else if "`panel'"=="G" local panel H
    else if "`panel'"=="H" local panel I
}
if `t5ffail' == 0 di as result "Table 5f saved: $tabs/table5f_exposure_interactions_flow.rtf"
else di as error "Table 5f written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* PART B — FLOW INTERACTION BY RESOLUTION TYPE (non-default vs default-linked)
*   d_nd = in_crisis_nd x exposure   d_def = in_crisis_def x exposure
*   Test H0: d_nd = d_def, a plain within-regression Wald test -- no new
*   inference machinery needed (unlike the AIPW files' bootstrap/Clogg-z).
* ══════════════════════════════════════════════════════════════════════════

foreach e of local expvars {
    foreach m in bnd blo_nd bhi_nd bdef blo_def bhi_def pdiff {
        matrix `m'_`e' = J(6, 1, 0)
    }
}

eststo clear

foreach e of local expvars {

    di as result _n "========================================================"
    di as result "PART B — EXPOSURE: `lbl_`e'' (by resolution type, flow-coded)"
    di as result "========================================================"
    di as result "h    d_nd      SE       d_def     SE       p(d_nd=d_def)"

    local elistB_`e'

    forvalues h = 0/4 {
        local lag = max(2, `h'+3)
        local row = `h' + 2

        capture xtscc dy_`h' in_crisis_nd in_crisis_def z_`e' Dnd_`e' Ddef_`e' ///
            `controls' `yearfe' if sample_flow == 1, fe lag(`lag')

        if _rc == 0 {
            capture test Dnd_`e' = Ddef_`e'
            if _rc == 0 local pd = r(p)
            else        local pd = .

            eststo b_`e'_`h', title("h=`=`h'+1'")
            estadd scalar pdiff = `pd'
            local elistB_`e' `elistB_`e'' b_`e'_`h'

            matrix bnd_`e'[`row',1]     = _b[Dnd_`e']
            matrix blo_nd_`e'[`row',1]  = _b[Dnd_`e']  - 1.645*_se[Dnd_`e']
            matrix bhi_nd_`e'[`row',1]  = _b[Dnd_`e']  + 1.645*_se[Dnd_`e']
            matrix bdef_`e'[`row',1]    = _b[Ddef_`e']
            matrix blo_def_`e'[`row',1] = _b[Ddef_`e'] - 1.645*_se[Ddef_`e']
            matrix bhi_def_`e'[`row',1] = _b[Ddef_`e'] + 1.645*_se[Ddef_`e']
            matrix pdiff_`e'[`row',1]   = `pd'

            di "h=" `h'+1 "   " %7.3f _b[Dnd_`e'] "  " %6.3f _se[Dnd_`e'] ///
               "   " %7.3f _b[Ddef_`e'] "  " %6.3f _se[Ddef_`e'] ///
               "   " %5.3f `pd'
        }
        else di as error "h=" `h'+1 ": xtscc failed for `e' by type (rc=" _rc ")"
    }
}

di as result _n "Interpretation (Part B):"
di as result "  Compare d_nd vs d_def: a more-negative d in one crisis type means"
di as result "  the channel amplifies the output loss more in that type."
di as result "  p(d_nd=d_def) < 0.10 => the amplification differs by resolution."
di as result "  Same sign caveat as Part A. With 90 tests in this file, ~4-5"
di as result "  significant cells are expected by chance: treat any single panel"
di as result "  as suggestive, not decisive."

* ── TABLE 6f: flow exposure interactions by resolution type (Word/RTF) ─────
local t6fnote "Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h, FLOW coding. Both flow-crisis dummies and both exposure interactions enter jointly. 'Flow-crisis(type) x exposure' is the extra output effect per 1 SD of pre-crisis exposure in that crisis type. p(nd=def) is the p-value of the equality of the two interaction terms (within-regression Wald test). Exposure entry-dated (see file header) and standardized over the flow sample. Country and year fixed effects; Driscoll-Kraay SE at lag(max(2,h+3)). * p<0.10, ** p<0.05, *** p<0.01. A negative coefficient is evidence of channel transmission ONLY for an exposure with a single economically sensible sign -- see file header. This file runs 9 exposures x 5 horizons x 2 panels = 90 tests; at the 5% level roughly 4-5 significant cells are expected by chance. No single panel should be read as decisive."

local panel A
local writemode replace
local t6ffail 0

foreach e of local expvars {
    if "`elistB_`e''" == "" {
        di as error "  ** Table 6f: no estimates for exposure `e' — panel skipped"
        local t6ffail 1
        continue
    }
    local t6fextra
    if "`e'" == "intpay_gni" local t6fextra addnotes("`t6fnote'")

    capture esttab `elistB_`e'' using "$tabs/table6f_exposure_by_resolution_flow.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(Dnd_`e' Ddef_`e') order(Dnd_`e' Ddef_`e') ///
        coeflabel(Dnd_`e' "Flow-crisis(non-default) x exposure (per SD)" ///
                  Ddef_`e' "Flow-crisis(default) x exposure (per SD)") ///
        mtitles nonumber ///
        stats(pdiff N N_g, labels("p (nd = def)" "Observations" "Countries") fmt(3 0 0)) ///
        title("Table 6f. Exposure amplification by resolution, FLOW coding -- Panel `panel': `lbl_`e''") ///
        `t6fextra'

    if _rc == 608 {
        di as error "  ** table6f_exposure_by_resolution_flow.rtf is OPEN IN WORD — close it and re-run."
        local t6ffail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 6f: esttab failed for exposure `e' (rc=" _rc ")"
        local t6ffail 1
        continue
    }
    local writemode append
    if      "`panel'"=="A" local panel B
    else if "`panel'"=="B" local panel C
    else if "`panel'"=="C" local panel D
    else if "`panel'"=="D" local panel E
    else if "`panel'"=="E" local panel F
    else if "`panel'"=="F" local panel G
    else if "`panel'"=="G" local panel H
    else if "`panel'"=="H" local panel I
}
if `t6ffail' == 0 di as result "Table 6f saved: $tabs/table6f_exposure_by_resolution_flow.rtf"
else di as error "Table 6f written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* RAW CSVs
* ══════════════════════════════════════════════════════════════════════════

* Pooled (Part A)
preserve
    clear
    local nexp : word count `expvars'
    local nobs = 5 * `nexp'
    set obs `nobs'
    gen exposure = ""
    gen horizon  = .
    gen d_coef   = .
    gen d_lo90   = .
    gen d_hi90   = .
    local row = 1
    foreach e of local expvars {
        forvalues h = 0/4 {
            replace exposure = "`e'"              in `row'
            replace horizon  = `h'+1              in `row'
            replace d_coef   = dcoef_`e'[`h'+2,1] in `row'
            replace d_lo90   = dlo90_`e'[`h'+2,1] in `row'
            replace d_hi90   = dhi90_`e'[`h'+2,1] in `row'
            local ++row
        }
    }
    order exposure horizon d_coef d_lo90 d_hi90
    export delimited "$tabs/exposure_interactions_flow.csv", replace
    di as result "Pooled flow exposure CSV saved: $tabs/exposure_interactions_flow.csv"
restore

* By type (Part B)
preserve
    clear
    local nexp : word count `expvars'
    local nobs = 5 * `nexp'
    set obs `nobs'
    gen exposure = ""
    gen horizon  = .
    gen d_nd     = .
    gen d_def    = .
    gen p_diff   = .
    local row = 1
    foreach e of local expvars {
        forvalues h = 0/4 {
            replace exposure = "`e'"               in `row'
            replace horizon  = `h'+1               in `row'
            replace d_nd     = bnd_`e'[`h'+2,1]    in `row'
            replace d_def    = bdef_`e'[`h'+2,1]   in `row'
            replace p_diff   = pdiff_`e'[`h'+2,1]  in `row'
            local ++row
        }
    }
    order exposure horizon d_nd d_def p_diff
    export delimited "$tabs/exposure_by_resolution_flow.csv", replace
    di as result "By-type flow exposure CSV saved: $tabs/exposure_by_resolution_flow.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* FIGURE A — pooled interaction d_h across horizons (grid)
* ══════════════════════════════════════════════════════════════════════════

local c_amp "157 36 73"
local fignamesA
local i = 1
foreach e of local expvars {
    preserve
        clear
        set obs 6
        gen horizon = _n - 1     // 0 (baseline), 1..5
        foreach m in dcoef dlo90 dhi90 {
            svmat `m'_`e', names(`m')
            rename `m'1 `m'
        }
        twoway ///
            (rarea dlo90 dhi90 horizon, color("`c_amp'%20") lwidth(none)) ///
            (connected dcoef horizon, ///
                lcolor("`c_amp'") lwidth(medthick) msymbol(circle) mcolor("`c_amp'")), ///
            yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
            xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
            xtitle("Year (Year 1 = crisis year)", size(small)) ///
            ytitle("Flow-crisis x exposure (pp/SD)", size(small)) ///
            title("`lbl_`e''", size(small) color(navy)) ///
            legend(off) graphregion(color(white)) plotregion(color(white)) ///
            name(gAf_`i', replace)
    restore
    local fignamesA `fignamesA' gAf_`i'
    local ++i
}
graph combine `fignamesA', cols(3) ///
    title("Exposure Heterogeneity, Flow Coding (Pooled): Output Loss by Pre-Crisis Exposure", ///
          size(medsmall) color(navy)) ///
    note("Interaction d_h (flow-crisis x standardized, entry-dated pre-crisis exposure), 90% CI. Below zero => output falls more where exposure was higher." ///
         "Driscoll-Kraay SE, lag(max(2,h+3)). Country & year FE. Same GDP outcome as the main flow LP (20_lp_flow.do).", size(vsmall)) ///
    graphregion(color(white)) xsize(11) ysize(7)
graph export "$figs/fig13f_exposure_interactions_flow.pdf", replace
di as result "Figure saved: fig13f_exposure_interactions_flow.pdf"
foreach nm of local fignamesA {
    capture graph drop `nm'
}

* ══════════════════════════════════════════════════════════════════════════
* FIGURE B — by-type: d_nd vs d_def across horizons (grid)
* ══════════════════════════════════════════════════════════════════════════

local c_nd  "34 139 34"
local c_def "157 36 73"
local fignamesB
local i = 1
foreach e of local expvars {
    preserve
        clear
        set obs 6
        gen horizon = _n - 1     // 0 (baseline), 1..5
        foreach m in bnd blo_nd bhi_nd bdef blo_def bhi_def {
            svmat `m'_`e', names(`m')
            rename `m'1 `m'
        }
        twoway ///
            (rarea blo_nd bhi_nd horizon, color("`c_nd'%15") lwidth(none)) ///
            (connected bnd horizon, ///
                lcolor("`c_nd'") lwidth(medthick) msymbol(triangle) mcolor("`c_nd'")) ///
            (rarea blo_def bhi_def horizon, color("`c_def'%15") lwidth(none)) ///
            (connected bdef horizon, ///
                lcolor("`c_def'") lwidth(medthick) lpattern(dash) ///
                msymbol(square) mcolor("`c_def'")), ///
            yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
            xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
            xtitle("Year (Year 1 = crisis year)", size(small)) ///
            ytitle("Flow-crisis x exposure (pp/SD)", size(small)) ///
            title("`lbl_`e''", size(small) color(navy)) ///
            legend(off) graphregion(color(white)) plotregion(color(white)) ///
            name(gBf_`i', replace)
    restore
    local fignamesB `fignamesB' gBf_`i'
    local ++i
}
graph combine `fignamesB', cols(3) ///
    title("Exposure Amplification by Resolution Type, Flow Coding", size(medsmall) color(navy)) ///
    subtitle("Green triangles = non-default; red squares = default-linked", size(small)) ///
    note("Interaction d (flow-crisis-by-type x standardized, entry-dated pre-crisis exposure), 90% CI. More-negative in a type => stronger amplification there." ///
         "Driscoll-Kraay SE, lag(max(2,h+3)). Country & year FE. Same GDP outcome as the main flow LP (20_lp_flow.do).", size(vsmall)) ///
    graphregion(color(white)) xsize(11) ysize(7)
graph export "$figs/fig13g_exposure_by_resolution_flow.pdf", replace
di as result "Figure saved: fig13g_exposure_by_resolution_flow.pdf"
foreach nm of local fignamesB {
    capture graph drop `nm'
}

* ══════════════════════════════════════════════════════════════════════════
* ONSET/FLOW SANITY CHECK -- the two headline UNAMBIGUOUS exposures only,
* against 13b's own onset-coded numbers. Same comparison pattern used
* between every other onset/flow file pair in this project (e.g. 25's
* closing message against 13d).
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "ONSET vs FLOW SANITY CHECK — the two UNAMBIGUOUS exposures"
di as result "════════════════════════════════════════════════════════════"
di as result "  13b (onset-coded) stdebt_share:      d_def = -6.20 (h3, p=.053), -5.77 (h4, p=.033)"
foreach e in stdebt_share {
    capture confirm matrix bdef_`e'
    if _rc == 0 {
        di as result "  This file (flow-coded) stdebt_share:  d_def = " ///
            %6.2f bdef_`e'[5,1] " (h3, p=" %5.3f pdiff_`e'[5,1] "), " ///
            %6.2f bdef_`e'[6,1] " (h4, p=" %5.3f pdiff_`e'[6,1] ")"
    }
}
di as result "  13b (onset-coded) claimsgov_assets:  d_def = -4.43 (h3, p=.034), -4.27 (h4, p=.029)"
foreach e in claimsgov_assets {
    capture confirm matrix bdef_`e'
    if _rc == 0 {
        di as result "  This file (flow-coded) claimsgov_assets: d_def = " ///
            %6.2f bdef_`e'[5,1] " (h3, p=" %5.3f pdiff_`e'[5,1] "), " ///
            %6.2f bdef_`e'[6,1] " (h4, p=" %5.3f pdiff_`e'[6,1] ")"
    }
}
di as result "  Sign agreement (not magnitude) is the primary sanity check here,"
di as result "  matching the pattern used between every other onset/flow file pair."

di as result _n "20b_exposure_heterogeneity_flow.do complete."
