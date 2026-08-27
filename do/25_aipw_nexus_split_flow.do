/*===========================================================================
  25_AIPW_NEXUS_SPLIT_FLOW.DO
  Sovereign-bank nexus heterogeneity split — FLOW treatment

  WHY THIS FILE EXISTS
  ---------------------
  13d_aipw_nexus_split.do (onset-coded) median-splits crisis onsets on
  pre-crisis bank exposure to the sovereign (claimsgov_assets, "the doom-loop
  amplifier") into high/low subsamples, then runs AIPW separately on each —
  for GDP plus four channels (credit, inv, claimpriv_assets, claims_govt) —
  and bootstraps the HIGH-LOW difference within each resolution-type cell.
  Headline claim: default-linked crises are costliest where the sovereign-bank
  nexus is tight pre-crisis.

  This file is the FLOW-treatment counterpart: the same heterogeneity split,
  estimated on in_crisis_nd/in_crisis_def (every year of an episode, onset and
  continuation alike — 18_transforms.do), not just the onset row.

  WHAT IS REUSED, UNCHANGED
  --------------------------
  * _aipw, _mkstrat, _aipwpairflow: copied verbatim from 24_aipw_channels_flow.do,
    which itself copies them verbatim from 21_aipw_flow.do. The propensity
    model's fix for flow-coding's near-separation problem (Eq. 2 fit on
    tranquil+onset rows only, continuation==0, then EXTRAPOLATED to every row)
    is not re-derived here — see 21's header ("DIAGNOSTIC HISTORY") for the
    settled justification; the diagnostic code itself is no longer kept there.
  * Channel outcome construction (ch_<var>_h, pre_<var>, epc_pre_<var>): copied
    from 22/23/24's identical loop, restricted to the four channels 13d uses
    (credit, inv, claims_govt, claimpriv_assets — govexp/pb/fdi were never
    part of 13d and are not added here).
  * cz_recency (years_since_def_onset in place of past_def_onsets,
    l_contagion_dist in place of l_reg_crisis_share): the adopted predictor
    set — see 24's header "PREDICTOR CHANGE" and 21_aipw_flow.do's header
    "SECOND PREDICTOR CHANGE".

  WHAT IS NEW
  -----------
  THE AMPLIFIER MUST BE ENTRY-DATED, NOT ROW-DATED. a_nexus is a PRE-CRISIS
  exposure level. Under onset coding this is unambiguous (L.claimsgov_assets at
  the single onset row). Under flow coding, a continuation row's own
  contemporaneous L.claimsgov_assets is partly a crisis OUTCOME — bank asset
  composition shifts as a crisis runs on — exactly the contamination problem
  18_transforms.do's epc_* construction exists to solve for every other
  control, and this file's own epc_pre_<var> channel controls solve for each
  channel's pre-trend. So a_nexus is built the same way: computed at the
  episode's ENTRY year (onset_all==1 row) and held FIXED across every row of
  that episode, continuation included; tranquil rows keep their own row-dated
  value (they are never part of a "high/low" episode).

  THE MEDIAN CUTOFF IS COMPUTED OVER ONSETS ONLY (sample==1 & onset_all==1),
  matching 13d exactly. The split threshold is a property of the onset-level
  population; it does not move just because flow coding adds continuation
  rows to the estimation sample. Every row of an episode inherits its one
  episode's high/low label from the onset row's a_nexus value.

  THE HIGH-LOW DIFFERENCE BOOTSTRAP reuses _aipwpairflow directly rather than
  writing a new estimator: for a given resolution type (nd or def), two
  0/1 indicator variables are built — "treated AND high-nexus" and "treated
  AND low-nexus" — and passed to _aipwpairflow as d1/d2 over the SAME sample
  condition (tranquil + this type's treated rows, rival type dropped). Because
  the indicators are disjoint (an episode is high-nexus xor low-nexus),
  _aipwpairflow's existing _pool row-bootstrap (0 = tranquil, 1 = d1==1,
  2 = d2==1) is already exactly the right stratification — no new program is
  needed.

  SCOPE: PART B (resolution x exposure) ONLY, NO PART A (pooled high/low).
  13d builds both. Every flow file in this project (20, 21, 22, 23, 24) has
  already dropped the pooled/Act-1 AIPW line and kept only the resolution
  split, because flow coding's near-tautological continuation persistence
  makes the pooled treatment-vs-tranquil comparison less informative once
  continuation rows are added. This file follows that same established
  boundary rather than reopening it: 4 cells only (nd-high, nd-low, def-high,
  def-low), plus the high-low difference within nd and within def.

  INFERENCE, MATCHING 21/24's CONVENTION EXACTLY
  ---------------------------------------------
  LEVELS (each of the 4 cells vs tranquil): analytic influence-function SE
  from _aipw, band theta +/- 1.96*se, PLUS a conventional t-test vs zero
  (b/se_a, normal-based, no df available from an influence-function SE),
  printed with stars next to HIGH/LOW exactly as 21's Section 3 and (now)
  24's Section 2 do.

  DIFFERENCE (high-low, within nd or def): paired row bootstrap via
  _aipwpairflow, percentile CI — the GOVERNING test, reported only when >=50
  valid draws (13d's own threshold, carried forward) — below that the point
  estimate stands with a printed caveat, not a fabricated interval. Also
  carries a Clogg et al. (1995) z as a companion statistic (`cloggz'/`cloggp'
  columns), matching 21's Section 3 and 24's Section 2 exactly. The Clogg z
  is the PERMISSIVE number: it assumes the two cells (high-nexus treated,
  low-nexus treated) are independent, which they are not here (they share
  the same tranquil control pool), so the bootstrap CI governs where the two
  disagree.

  SECTION 4b, PORTED FROM 21/24 -- OFF BY DEFAULT: the same high-low
  difference can be re-estimated resampling whole COUNTRIES instead of rows,
  per outcome/type/horizon, and compared against the row-bootstrap CI. Set
  `run_cluster_boot' to 1 in Section 4b to run it -- it is a second full
  bootstrap pass (roughly doubling this file's already-heavier runtime), so
  it is left off by default. `aipw_nexus_split_flow_diff.csv' always carries
  the clu_lo/clu_hi/sig95_clu columns; they are all missing/0 when the
  toggle is off.

  Outputs
  -------
    "$tabs/aipw_nexus_split_flow.csv"       outcome x cell(nd_high/nd_low/def_high/def_low) x horizon
    "$tabs/aipw_nexus_split_flow_diff.csv"  outcome x type(nd/def) x horizon, high-low gap + CI,
                                             Clogg z, and (if run_cluster_boot=1) country-cluster CI
    "$figs/fig_aipw_nexus_split_flow.pdf"   GDP, nd/def panels, high vs low overlay
    "$figs/fig_nexus_flow_<channel>.pdf"    one per channel (credit, inv, claims_govt, claimpriv_assets)

  RUNTIME: 5 outcomes x 2 types x 5 horizons x nboot draws for the difference
  alone, on top of the 4 point estimates per outcome/horizon. Heavier than 24
  (5 outcomes here vs 6 channels there, but each outcome carries 2 types x 2
  bank-exposure cells vs 24's 2 resolution cells). Section 4b (OFF by
  default) repeats it all again in CLUSTER mode when turned on. Raise nboot
  only for a final run.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow ep_seq ///
             years_since_def_onset claimsgov_assets l_contagion_dist {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel_lp.dta — re-run 18_transforms.do / 01c_merge_nexus.do first."
        exit 111
    }
}
if "$ctrl_flow" == "" {
    di as error "  ** \$ctrl_flow not set — re-run 18_transforms.do, then this file."
    exit 111
}

* ══════════════════════════════════════════════════════════════════════════
* 1. CHANNEL OUTCOMES — identical construction to 22/23/24, restricted to
*    13d's four channels (credit, inv, claims_govt, claimpriv_assets)
* ══════════════════════════════════════════════════════════════════════════
local channels credit inv claims_govt claimpriv_assets

foreach var of local channels {
    local src `var'
    if inlist("`var'","credit","inv") local src ln_r_`var'
    capture drop `var'_base
    gen `var'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`src' - `var'_base
    }
    capture drop pre_`var'
    gen pre_`var' = L.`src' - L2.`src'
    capture drop epc_pre_`var' _ent_pre_`var'
    quietly bysort cid ep_seq: egen double _ent_pre_`var' = ///
        max(cond(onset_all==1, pre_`var', .))
    quietly gen double epc_pre_`var' = cond(in_crisis==1, _ent_pre_`var', pre_`var')
    quietly drop _ent_pre_`var'
    * See 22_channels_flow.do for why this sort is needed inside the loop.
    sort cid year
}
sort cid year

* ══════════════════════════════════════════════════════════════════════════
* 2. THE AMPLIFIER — pre-crisis sovereign-bank nexus, ENTRY-DATED
* ══════════════════════════════════════════════════════════════════════════
capture drop a_nexus_row a_nexus_cm a_nexus _ent_a_nexus highbank_flow
gen a_nexus_row = L.claimsgov_assets
bysort cid: egen a_nexus_cm = mean(claimsgov_assets)
replace a_nexus_row = a_nexus_cm if missing(a_nexus_row)

* Entry-dated: the episode's own onset-year value, held fixed across every
* continuation row of that episode. Tranquil rows keep their own row-dated
* value. Same construction as $ctrl_flow's epc_* terms and this file's own
* epc_pre_<var> channel controls above.
quietly bysort cid ep_seq: egen double _ent_a_nexus = ///
    max(cond(onset_all==1, a_nexus_row, .))
gen double a_nexus = cond(in_crisis==1, _ent_a_nexus, a_nexus_row)
drop _ent_a_nexus
sort cid year

* Median cutoff over ONSETS ONLY, matching 13d exactly — the threshold is a
* property of the onset population, not shifted by adding continuation rows.
quietly summarize a_nexus_row if sample==1 & onset_all==1, detail
local med = r(p50)
gen highbank_flow = (a_nexus >= `med') if !missing(a_nexus)
label define hbf 0 "Low nexus" 1 "High nexus"
label values highbank_flow hbf

di as result _n "════════════════════════════════════════════════════════════"
di as result "SOVEREIGN-BANK NEXUS MEDIAN SPLIT — FLOW (entry-dated amplifier)"
di as result "════════════════════════════════════════════════════════════"
di as result "  Amplifier = claims on govt / bank assets at episode entry (country-mean filled)"
di as result "  Median cutoff among crisis onsets = " %6.2f `med'
foreach t in in_crisis in_crisis_nd in_crisis_def {
    quietly count if sample_flow==1 & `t'==1 & highbank_flow==1
    local nh = r(N)
    quietly count if sample_flow==1 & `t'==1 & highbank_flow==0
    local nl = r(N)
    quietly count if sample_flow==1 & `t'==1 & missing(highbank_flow)
    local nm = r(N)
    di as result "  `t' (flow rows): high=" `nh' "  low=" `nl' "  unclassified=" `nm'
}

di as result _n "=== NEXUS-BIN COUNTRY COMPOSITION (crisis onsets) ==="
capture noisily tabulate country highbank_flow if sample==1 & onset_all==1, row nofreq
di as result "  (If the two columns are not obviously split by development, the nexus"
di as result "   result is not merely a development proxy — same check as 13d's header.)"

* ══════════════════════════════════════════════════════════════════════════
* 3. CONTROL SETS — same convention as 22/23/24
* ══════════════════════════════════════════════════════════════════════════
local epc_lc epc_l_credit_bank

* CONTROL SET. flow_ctrl_variant: 0 = the ADOPTED flow-tier baseline,
* $ctrl_flow, built in 18_transforms.do from $ctrl_core_flowbase
* (l_banking_duration -> the l_banking_crisis DUMMY, l_ca -> tot_chg -> exchange2,
* l_hyperinfl -> l_lninfl -- see that file's "ADOPTED FLOW-TIER CORE CONTROL
* SET"). 1 = the adopted core PLUS the reference paper's own additional
* predictors, tot_chg (the term exchange2 replaced) and l_imf (18_transforms.do's "ALTERNATE FLOW
* CONTROL SET"). Default 0, matching 20/22/23/24's toggle of the same name.
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
if `flow_ctrl_variant'==1 local ctrl_flow_base $ctrl_flow_flowplus
else                       local ctrl_flow_base $ctrl_flow

local ctrl_credit           : list ctrl_flow_base - epc_lc
local ctrl_credit           `ctrl_credit' epc_pre_credit
local ctrl_inv               `ctrl_flow_base' epc_pre_inv
local ctrl_claims_govt       `ctrl_flow_base' epc_pre_claims_govt
local ctrl_claimpriv_assets  `ctrl_flow_base' epc_pre_claimpriv_assets
local ctrl_gdp                `ctrl_flow_base'

* ── Propensity model — SAME as 21/24's ACTIVE (`cx_active') baseline ───────
local cx = cond(`flow_ctrl_variant'==1, "$ctrl_core_flowplus", "$ctrl_core_flowbase")
* l_contagion_dist (distance-weighted, CEPII great-circle) in place of
* l_reg_crisis_share -- matches 21_aipw_flow.do's cz_recency, "SECOND
* PREDICTOR CHANGE".
local cz_recency l_fedfunds l_contagion_dist years_since_def_onset

set seed 20260819
local nboot = 300

sort cid year   // belt-and-suspenders before estimation

* ══════════════════════════════════════════════════════════════════════════
* PROGRAMS — copied verbatim from 24_aipw_channels_flow.do (itself copied
* from 21_aipw_flow.do). See 21's header for the propensity-fit-sample
* justification; it is not re-derived here.
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'

    tempvar xb m0 m1 ps summ iwt
    quietly probit `D' `pmodel' if `touse' & continuation==0
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01                & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'

    if "`fe'" != "" quietly reg `y' `D' `omodel' i.`fe' [pweight=`iwt'] if `touse'
    else            quietly reg `y' `D' `omodel'         [pweight=`iwt'] if `touse'
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']     if `touse'
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    local th = r(mean)
    local nn = r(N)

    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta = `th'
    return scalar N     = `nn'
    return scalar se    = `sean'
end

capture program drop _mkstrat
program define _mkstrat
    syntax varlist(min=1 max=2) [if], GENerate(name)
    marksample touse
    capture drop `generate'
    tempvar t1 t2 s2
    local d1 : word 1 of `varlist'
    local d2 : word 2 of `varlist'
    quietly gen byte `t1' = `d1' if `touse'
    bysort cid: egen byte `generate' = max(`t1')
    quietly replace `generate' = 0 if missing(`generate')
    if "`d2'" != "" {
        quietly gen byte `t2' = `d2' if `touse'
        bysort cid: egen byte `s2' = max(`t2')
        quietly replace `s2' = 0 if missing(`s2')
        quietly replace `generate' = `generate' + 2*`s2'
    }
end

capture program drop _aipwpairflow
program define _aipwpairflow, rclass
    syntax , Y(string) D1(string) IF1(string) D2(string) IF2(string) ///
             OMOD(string) PZ(string) REPS(integer) [BOOT(string)]
    if "`boot'" == "" local boot row
    if !inlist("`boot'","row","cluster") {
        di as error "  ** boot() must be row or cluster"
        exit 198
    }
    capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b1  = r(theta)
    local a1  = r(se)
    capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local b2  = r(theta)
    local a2  = r(se)
    local dh  = `b1' - `b2'

    if "`boot'" == "row" {
        capture drop _pool
        quietly gen byte _pool = 0 if (`if1') | (`if2')
        quietly replace _pool = 1 if `d1' == 1
        quietly replace _pool = 2 if `d2' == 1
    }
    else {
        _mkstrat `d1' `d2', generate(_strat)
    }

    tempname pf
    tempfile bf
    quietly postfile `pf' double t1 double t2 double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            if "`boot'" == "row" {
                quietly keep if !missing(_pool)
                quietly bsample, strata(_pool)
                capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(cid)
                local v1 = cond(_rc==0, r(theta), .)
                capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(cid)
                local v2 = cond(_rc==0, r(theta), .)
            }
            else {
                capture drop _bid
                bsample, cluster(cid) strata(_strat) idcluster(_bid)
                capture _aipw `y' `d1' if `if1', omodel(`omod') pmodel(`pz') fe(_bid)
                local v1 = cond(_rc==0, r(theta), .)
                capture _aipw `y' `d2' if `if2', omodel(`omod') pmodel(`pz') fe(_bid)
                local v2 = cond(_rc==0, r(theta), .)
            }
            if !missing(`v1') & !missing(`v2') quietly post `pf' (`v1') (`v2') (`v1'-`v2')
        restore
    }
    quietly postclose `pf'
    capture drop _strat _pool

    local se = .
    local lo = .
    local hi = .
    local nd = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(diff)
        local nd = r(N)
        if `nd' >= 50 {
            quietly summarize diff
            local se = r(sd)
            _pctile diff, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
        }
    restore
    return scalar ok = 1
    return scalar dh = `dh'
    return scalar b1 = `b1'
    return scalar b2 = `b2'
    return scalar a1 = `a1'
    return scalar a2 = `a2'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* 4. ESTIMATION — per outcome x resolution type: high-nexus cell, low-nexus
*    cell (each vs tranquil, analytic SE), plus the HIGH-LOW difference.
*
*    The high-low difference reuses _aipwpairflow directly rather than a new
*    program: d1/d2 are 0/1 indicators for "treated AND high-nexus" /
*    "treated AND low-nexus" (disjoint by construction, since an episode is
*    high-nexus xor low-nexus), passed over the SAME sample condition
*    (tranquil + this type's treated rows, rival type dropped). _aipwpairflow's
*    existing _pool row-bootstrap (0=tranquil, 1=d1==1, 2=d2==1) is then
*    exactly the right stratification for this comparison — no new estimator.
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str24 outcome str8 cell byte horizon double b double se double lo double hi ///
    using "`resf'", replace

tempname Rd
tempfile diffresf
postfile `Rd' str24 outcome str4 restype byte horizon double dhl double bhigh double blow ///
    double se double lo double hi long nd double cloggz double cloggp using "`diffresf'", replace

local outcomes gdp credit inv claims_govt claimpriv_assets

foreach oc in `outcomes' {
    local ystem = cond("`oc'"=="gdp", "dy", "ch_`oc'")
    local om `ctrl_`oc''

    foreach g in nd def {
        foreach h2 in high low {
            matrix b_`oc'_`g'_`h2' = J(7, 1, .)
            matrix lo_`oc'_`g'_`h2' = J(7, 1, .)
            matrix hi_`oc'_`g'_`h2' = J(7, 1, .)
            matrix b_`oc'_`g'_`h2'[2,1] = 0
            matrix lo_`oc'_`g'_`h2'[2,1] = 0
            matrix hi_`oc'_`g'_`h2'[2,1] = 0
        }
        * Row-bootstrap diff CI, saved per outcome/type/horizon so Section 4b
        * (country-cluster comparison) can print the width ratio against it.
        matrix Ddlo_`oc'_`g' = J(5,1,.)
        matrix Ddhi_`oc'_`g' = J(5,1,.)
    }
}

di as result _n "════════════════════════════════════════════════════════════"
di as result "FLOW AIPW NEXUS SPLIT (Year 1 = the crisis year)"
di as result "Levels: analytic influence-function SE. High-low diff: paired bootstrap."
di as result "HIGH/LOW stars are the conventional t-test vs zero (b/se_a): * p<.10 ** p<.05 *** p<.01."
di as result "hi-lo's own * marks the bootstrap CI excluding 0 -- the conservative, governing test."
di as result "════════════════════════════════════════════════════════════"

foreach oc in `outcomes' {
    local ystem = cond("`oc'"=="gdp", "dy", "ch_`oc'")
    local om `ctrl_`oc''
    di as result _n "############### OUTCOME: `oc' ###############"

    foreach g in nd def {
        local Dv   in_crisis_`g'
        local rivg = cond("`g'"=="nd", "def", "nd")
        local riv  in_crisis_`rivg'

        capture drop _dhigh_`oc'_`g' _dlow_`oc'_`g'
        quietly gen byte _dhigh_`oc'_`g' = (`Dv'==1 & highbank_flow==1)
        quietly gen byte _dlow_`oc'_`g'  = (`Dv'==1 & highbank_flow==0)
        local ifc "sample_flow==1 & `riv'==0"

        post `Rd' ("`oc'") ("`g'") (0) (0) (0) (0) (0) (0) (0) (0)

        forvalues h = 0/4 {
            local hd  = `h' + 1
            local row = `h' + 3

            _aipwpairflow, y(`ystem'_`h') ///
                d1(_dhigh_`oc'_`g') if1(`ifc') ///
                d2(_dlow_`oc'_`g')  if2(`ifc') ///
                omod(`om') pz(`cx' `cz_recency') reps(`nboot') boot(row)

            if !r(ok) {
                di as error "  h=`hd' `g': estimate failed for `oc' (cell too thin)."
                continue
            }

            local BH = r(b1)     // high-nexus level ATE
            local BL = r(b2)     // low-nexus  level ATE
            local AH = r(a1)     // analytic SE, high
            local AL = r(a2)     // analytic SE, low
            local DH = r(dh)
            local SE = r(se)
            local LO = r(lo)
            local HI = r(hi)
            local ND = r(nd)

            matrix b_`oc'_`g'_high[`row',1]  = `BH'
            matrix lo_`oc'_`g'_high[`row',1] = `BH' - 1.96*`AH'
            matrix hi_`oc'_`g'_high[`row',1] = `BH' + 1.96*`AH'
            matrix b_`oc'_`g'_low[`row',1]   = `BL'
            matrix lo_`oc'_`g'_low[`row',1]  = `BL' - 1.96*`AL'
            matrix hi_`oc'_`g'_low[`row',1]  = `BL' + 1.96*`AL'

            post `R' ("`oc'") ("`g'_high") (`hd') (`BH') (`AH') (`BH'-1.96*`AH') (`BH'+1.96*`AH')
            post `R' ("`oc'") ("`g'_low")  (`hd') (`BL') (`AL') (`BL'-1.96*`AL') (`BL'+1.96*`AL')

            * Clogg et al. (1995) z, built from the same analytic influence-
            * function SEs used for the level bands -- identical construction to
            * 21_aipw_flow.do's Section 3 and 24's Section 2. PERMISSIVE: treats
            * the two cells as independent, which they are not (shared tranquil
            * pool), so the bootstrap CI governs where the two disagree.
            local zz = .
            local pz = .
            if !missing(`AH') & !missing(`AL') & (`AH'^2 + `AL'^2) > 0 {
                local zz = `DH' / sqrt(`AH'^2 + `AL'^2)
                local pz = 2*(1 - normal(abs(`zz')))
            }

            * Conventional t-test for each level vs zero, identical construction
            * to 21/24's Section 3/2 (b/se_a, normal-based -- an influence-
            * function SE carries no regression df to build a t distribution from).
            local thi  = cond(`AH'>0, `BH'/`AH', .)
            local phi  = cond(!missing(`thi'), 2*(1-normal(abs(`thi'))), .)
            local sghi = cond(missing(`phi'), "", cond(`phi'<.01,"***",cond(`phi'<.05,"**",cond(`phi'<.10,"*",""))))
            local tlo  = cond(`AL'>0, `BL'/`AL', .)
            local plo  = cond(!missing(`tlo'), 2*(1-normal(abs(`tlo'))), .)
            local sglo = cond(missing(`plo'), "", cond(`plo'<.01,"***",cond(`plo'<.05,"**",cond(`plo'<.10,"*",""))))

            local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
            di as result "  h=" %1.0f `hd' " `g'  HIGH=" %8.3f `BH' "`sghi'" " (" %5.3f `AH' ")" ///
                 "  LOW=" %8.3f `BL' "`sglo'" " (" %5.3f `AL' ")" ///
                 "  hi-lo=" %8.3f `DH' " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
                 "  " %4.0f `ND' "/`nboot'" ///
                 "  Clogg z=" %7.3f `zz' " p=" %5.3f `pz'

            post `Rd' ("`oc'") ("`g'") (`hd') (`DH') (`BH') (`BL') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz')
            matrix Ddlo_`oc'_`g'[`hd',1] = `LO'
            matrix Ddhi_`oc'_`g'[`hd',1] = `HI'
        }
        capture drop _dhigh_`oc'_`g' _dlow_`oc'_`g'
        sort cid year
    }
}
postclose `R'
postclose `Rd'

* ══════════════════════════════════════════════════════════════════════════
* 4b. THE SAME HIGH-LOW DIFFERENCE, COUNTRY-CLUSTER BOOTSTRAP — COMPARISON
*     ONLY, OFF BY DEFAULT (set run_cluster_boot to 1 to run it)
*
* Ported from 21_aipw_flow.do's Section 3b / 24_aipw_channels_flow.do's
* Section 2b, per outcome/type. Section 4 resamples ROWS within treatment-
* type strata; under flow coding a treated row is a crisis-YEAR, not an
* independent episode, and a handful of countries contribute many rows to
* one bank-exposure cell. This block re-runs the identical estimator
* resampling whole COUNTRIES instead. Point estimates are identical by
* construction (a bootstrap does not touch them) -- only the interval
* differs.
*
* OFF BY DEFAULT: this is a second full bootstrap pass (5 outcomes x 2 types
* x 5 horizons x nboot draws again), roughly doubling this file's already
* heavier runtime. The export block always merges against `clusterf' -- when
* this is skipped that file is an empty shell, so clu_lo/clu_hi/clu_nd/
* sig95_clu still appear in the CSV but are all missing/0.
* ══════════════════════════════════════════════════════════════════════════
local run_cluster_boot 0

tempname Rc
tempfile clusterf
postfile `Rc' str24 outcome str4 restype byte horizon double clu_lo double clu_hi long clu_nd ///
    using "`clusterf'", replace
postclose `Rc'   // empty shell so Section 5's merge has something to read when the block below is skipped

if `run_cluster_boot' {
    tempname Rc2
    postfile `Rc2' str24 outcome str4 restype byte horizon double clu_lo double clu_hi long clu_nd ///
        using "`clusterf'", replace

    di as result _n "════════════════════════════════════════════════════════════"
    di as result "4b. ROW vs COUNTRY-CLUSTER BOOTSTRAP, PER OUTCOME/TYPE (same point estimates)"
    di as result "════════════════════════════════════════════════════════════"

    foreach oc in `outcomes' {
        local ystem = cond("`oc'"=="gdp", "dy", "ch_`oc'")
        local om `ctrl_`oc''
        di as result _n "############### OUTCOME: `oc' ###############"

        foreach g in nd def {
            local Dv   in_crisis_`g'
            local rivg = cond("`g'"=="nd", "def", "nd")
            local riv  in_crisis_`rivg'

            capture drop _dhigh_`oc'_`g' _dlow_`oc'_`g'
            quietly gen byte _dhigh_`oc'_`g' = (`Dv'==1 & highbank_flow==1)
            quietly gen byte _dlow_`oc'_`g'  = (`Dv'==1 & highbank_flow==0)
            local ifc "sample_flow==1 & `riv'==0"

            post `Rc2' ("`oc'") ("`g'") (0) (.) (.) (0)

            forvalues h = 0/4 {
                local hd = `h' + 1

                _aipwpairflow, y(`ystem'_`h') ///
                    d1(_dhigh_`oc'_`g') if1(`ifc') ///
                    d2(_dlow_`oc'_`g')  if2(`ifc') ///
                    omod(`om') pz(`cx' `cz_recency') reps(`nboot') boot(cluster)

                if !r(ok) {
                    di as error "  h=`hd' `g': cluster-bootstrap comparison failed for `oc'."
                    continue
                }
                local CD = r(dh)
                local CL = r(lo)
                local CH = r(hi)
                local CN = r(nd)

                local rowlo = Ddlo_`oc'_`g'[`hd',1]
                local rowhi = Ddhi_`oc'_`g'[`hd',1]
                local wrow = `rowhi' - `rowlo'
                local wclu = `CH' - `CL'
                local wrat = cond(`wclu' > 0, `wrow'/`wclu', .)
                local sigr = cond(!missing(`rowlo') & (`rowlo'>0 | `rowhi'<0), "*", " ")
                local sigc = cond(!missing(`CL') & (`CL'>0 | `CH'<0), "*", " ")

                di as result "  h=" %1.0f `hd' " `g'  hi-lo=" %8.3f `CD' ///
                     "  [row: " %7.3f `rowlo' ", " %7.3f `rowhi' "]`sigr'" ///
                     "  [cluster: " %7.3f `CL' ", " %7.3f `CH' "]`sigc'" ///
                     "  width ratio=" %5.2f `wrat' "  " %4.0f `CN' "/`nboot'"

                post `Rc2' ("`oc'") ("`g'") (`hd') (`CL') (`CH') (`CN')
            }
            capture drop _dhigh_`oc'_`g' _dlow_`oc'_`g'
            sort cid year
        }
    }
    postclose `Rc2'

    di as result _n "  A width ratio well below 1 means the row bootstrap is treating"
    di as result "  repeated crisis-years of the same country as independent draws."
    di as result "  Where the two verdicts (*) differ, the write-up reports both and"
    di as result "  says which resampling unit produced which."
}
else {
    di as result _n "  4b. Country-cluster bootstrap comparison SKIPPED (run_cluster_boot=0)."
    di as result "      Set run_cluster_boot to 1 above to run it -- see file header."
}

* ══════════════════════════════════════════════════════════════════════════
* 5. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
preserve
    use "`diffresf'", clear
    quietly merge 1:1 outcome restype horizon using "`clusterf'", nogenerate
    label var dhl    "AIPW flow high - low nexus gap (pp)"
    label var bhigh  "High-nexus ATE (analytic SE band)"
    label var blow   "Low-nexus ATE (analytic SE band)"
    label var se     "Bootstrap SD of the difference"
    label var lo     "95% percentile CI lower (row bootstrap, baseline)"
    label var hi     "95% percentile CI upper (row bootstrap, baseline)"
    label var nd     "Valid row-bootstrap draws"
    gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
    label var sig95 "Row bootstrap CI excludes 0 (baseline, the paper's scheme)"
    label var clu_lo "95% CI lower, COUNTRY-CLUSTER bootstrap (comparison)"
    label var clu_hi "95% CI upper, COUNTRY-CLUSTER bootstrap (comparison)"
    label var clu_nd "Valid country-cluster bootstrap draws"
    gen byte sig95_clu = (!missing(clu_lo) & (clu_lo>0 | clu_hi<0))
    label var sig95_clu "Country-cluster CI excludes 0 (comparison)"
    label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
    label var cloggp "p-value of the Clogg z"
    order outcome restype horizon dhl bhigh blow se lo hi nd sig95 clu_lo clu_hi clu_nd sig95_clu cloggz cloggp
    export delimited "$tabs/aipw_nexus_split_flow_diff.csv", replace
    di as result _n "AIPW flow nexus high-low difference CSV saved: $tabs/aipw_nexus_split_flow_diff.csv"
restore

use "`resf'", clear
label var b  "AIPW ATE (outcome units, cumulative change)"
label var se "Analytic influence-function SE"
label var lo "Level 95% CI lower (theta +/- 1.96*se)"
label var hi "Level 95% CI upper"
order outcome cell horizon b se lo hi
export delimited "$tabs/aipw_nexus_split_flow.csv", replace
di as result _n "AIPW flow nexus split results CSV saved: $tabs/aipw_nexus_split_flow.csv"

* ── IRF datasets + figures ───────────────────────────────────────────────
foreach oc in `outcomes' {
    foreach g in nd def {
        foreach h2 in high low {
            preserve
                clear
                set obs 7
                gen horizon = _n - 2
                foreach m in b lo hi {
                    svmat `m'_`oc'_`g'_`h2', names(`m')
                    rename `m'1 `m'
                }
                gen bank = "`h2'"
                save "$clean/irf_nexus_flow_`oc'_`g'_`h2'.dta", replace
            restore
        }
    }
}
di as result "IRF datasets saved: irf_nexus_flow_<outcome>_<nd|def>_<high|low>.dta (20)"

local c_hi "157 36 73"
local c_lo "0 84 166"
local c_zero "150 150 150"
local ptit_gdp "GDP"
local ptit_credit "Private credit"
local ptit_inv "Investment"
local ptit_claims_govt "Bank claims on govt"
local ptit_claimpriv_assets "Nexus: claims on private"

foreach oc in `outcomes' {
    local fnm = cond("`oc'"=="gdp", "fig_aipw_nexus_split_flow", "fig_nexus_flow_`oc'")
    local i = 1
    local gnames
    foreach g in nd def {
        preserve
            use "$clean/irf_nexus_flow_`oc'_`g'_high.dta", clear
            append using "$clean/irf_nexus_flow_`oc'_`g'_low.dta"
            keep if horizon >= 0
            local gt = cond("`g'"=="nd", "Non-default", "Default-linked")
            twoway ///
                (rarea lo hi horizon if bank=="high", color("`c_hi'%16") lwidth(none)) ///
                (rarea lo hi horizon if bank=="low",  color("`c_lo'%16") lwidth(none)) ///
                (connected b horizon if bank=="high", lcolor("`c_hi'") mcolor("`c_hi'") msymbol(square) lpattern(dash) lwidth(medthick)) ///
                (connected b horizon if bank=="low",  lcolor("`c_lo'") mcolor("`c_lo'") msymbol(circle) lwidth(medthick)), ///
                yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
                xlabel(0(1)5, labsize(small)) ylabel(, format(%4.1f) labsize(small)) ///
                xtitle("Year", size(small)) ytitle("", size(small)) ///
                title("`gt'", size(small) color(navy)) legend(off) ///
                graphregion(color(white)) plotregion(color(white)) name(gnx`i', replace)
        restore
        local gnames `gnames' gnx`i'
        local ++i
    }
    graph combine `gnames', cols(2) ///
        title("`ptit_`oc'' by Sovereign-Bank Nexus — Flow Specification", size(medsmall) color(navy)) ///
        subtitle("Red squares = high nexus (doom-loop); blue circles = low nexus. Shaded = 95% analytic-SE band.", size(small)) ///
        graphregion(color(white)) xsize(9) ysize(5)
    capture graph export "$figs/`fnm'.pdf", replace
    if _rc di as error "  ** `fnm' export failed (rc=" _rc ")"
    else {
        capture graph export "$figs/`fnm'.png", replace width(1200)
        di as result "Figure saved: `fnm'.pdf/.png"
    }
    foreach nm of local gnames {
        capture graph drop `nm'
    }
}

di as result _n "25_aipw_nexus_split_flow.do complete."
di as result "  Compare GDP's def-high vs def-low sign to 13d_aipw_nexus_split.do's own"
di as result "  onset-coded headline (high-nexus default crises deeper) as the same"
di as result "  onset/flow sanity-check relationship used throughout this project."
