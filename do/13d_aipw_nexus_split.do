/*===========================================================================
  13D_AIPW_NEXUS_SPLIT.DO
  Bank-intermediation heterogeneity — DOOM-LOOP version (Asonuma et al. §4 / Fig 6
  analog). Their headline two-dimensional result: the output cost of a crisis
  depends on the resolution strategy AND on bank intermediation. They median-split
  by bank-credit/GDP ("large" vs "small" banking sector) and estimate AIPW
  separately on each subsample. We go one better: split by the SOVEREIGN-BANK
  NEXUS — claimsgov_assets (bank claims on government / total assets, the doom-loop
  intensity Asonuma didn't have) — in a spread-crisis setting.

  Headline claim: default-linked spread crises are costliest where the
  sovereign-bank nexus is tight (banks heavily exposed to the sovereign).

  Outcomes: GDP (dy_h, coherent with 08b) PLUS the transmission channels
  credit, inv, claimpriv_assets, claims_govt — so we can watch each channel evolve
  differently under high vs low nexus (which channel carries the non-default
  cushion vs the default-linked doom-loop loss).

  Design (per outcome, coherent with 08b/13c):
    Amplifier a_nexus = pre-crisis claimsgov_assets (L.claimsgov_assets), filled
      with the country mean where the t-1 value is missing (coverage is thin, 2001+).
    Median split over crisis onsets: highbank = a_nexus >= median(onsets).
    AIPW cells (control = tranquil years only; treated = the specific cell; rival
    onsets dropped), estimated with the same _aipw as 08b/13c (levels from the
    analytic SE, only the high-low difference bootstrapped -- see below):
      Part A (robust headline, all onsets x {high, low}) is SILENCED -- only
        the resolution split below is of interest now.
      Part B (two-dimensional): {nd, def} x {high, low}          -> 4 lines

  SMALL-SAMPLE CAVEAT: the default x {high,low} cells hold only a handful of events.
  Level CIs no longer depend on the bootstrap at all (analytic SE, always reported);
  only the HIGH-LOW DIFFERENCE is bootstrapped, and its CI is shown only when >=50
  valid draws, else the point gap stands with a printed caveat. That threshold is
  ABSOLUTE, so with nboot now at 1000 read the printed nd/nboot RATE, not just the
  count: 292 of 300 is a healthy cell, 292 of 1000 would mean seven draws in ten
  failed to estimate and the surviving subset is not representative. This thinness
  is the honest limit vs Asonuma's 194 restructurings -- read the def x {high,low}
  cell as suggestive.

  Coherence with the Asonuma replication (their Fig 6 / IPWRA engine): (a) each
  channel outcome model controls for the channel's own pre-crisis change pre_<v>
  (their g_0); (b) besides the two level IRFs we bootstrap the HIGH - LOW nexus
  DIFFERENCE within each cell (their within-type high-vs-low contrast) so the
  sign-flip is formally tested, not just eyeballed from two bands.

  INFERENCE, ALIGNED WITH 08b_aipw.do / 13c_aipw_channels.do / THE FLOW
  TIER'S PRESENTATION (21_aipw_flow.do): each (part, bank) level's CI is
  1.96*analytic influence-function SE (the paper's own Table 2/Fig. 4
  construction), with a conventional t-test (b/se) and stars against zero.
  The HIGH-LOW nexus difference within each part is bootstrapped directly
  with ROW-LEVEL resampling within control/high/low pools -- the paper's
  own bootstrap device, a natural fit since an onset row already is one
  episode. Clogg et al. (1995)'s z (from the same analytic SEs as the
  level bands) is reported alongside the bootstrap CI as the permissive
  companion statistic (the two cells share the same tranquil control
  pool, so they are not independent). See 08b_aipw.do's header for the
  full argument.

  Output: $tabs/aipw_nexus_split.csv (outcome x part x bank x horizon, levels) ;
          $tabs/aipw_nexus_diff.csv  (outcome x part x horizon, high-low gap + CI) ;
          $figs/fig_aipw_nexus_split.pdf (GDP) + $figs/fig_nexus_<channel>.pdf.
  Run AFTER 17_predictors.do (needs fedfunds, past_onsets, l_contagion_dist_def,
  years_since_def_onset, nexus vars).
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ── REPRODUCIBILITY: seed the bootstrap ────────────────────────────────────
* Every CI in this file comes from `bsample', which draws at random. Without a
* seed the intervals move between runs of identical code, and on cells this
* thin that is not a rounding issue: the non-default high-minus-low nexus gap
* at Year 1 came back [0.05, 7.37] in one run and [-0.088, 7.686] in the next
* -- opposite verdicts at the 5% line from the same data. Whichever number
* reached the write-up would then depend on which run happened to be cited.
* Seeding fixes the draws so the reported intervals are a property of the
* estimator rather than of the session. The value is arbitrary and was not
* chosen by inspecting results.
set seed 20260819

* ── DRAW COUNT: why 1000 and not 300 ───────────────────────────────────────
* Seeding made the intervals reproducible; it did not make them RELIABLE, and
* the seeded run showed why. Comparing it with the unseeded run -- identical
* data, identical point estimates -- three cells crossed the 5% line purely on
* the change of draws:
*     GDP,         nd,  high-low, h=1 : [-0.088, 7.686] -> [ 0.167, 7.553]
*     GDP,         def x high,    h=3 : [-20.07, -0.42] -> [-20.03, 0.057]
*     claims_govt, def, high-low, h=3 : [-0.461, 12.758] -> [0.972, 12.764]
* Two gained significance, one lost it. A verdict that depends on the seed is
* not a verdict the data has delivered.
*
* The cause is Monte Carlo error in the percentile endpoint, not sampling
* error in the estimate. At 300 draws the 2.5th percentile is about the 7th
* smallest value, and in these cells only 235-274 draws survive estimation, so
* it is roughly the 6th of 235 -- an order statistic estimated with enough
* noise to move the endpoint by more than its distance from zero.
*
* Raising to 1000 cuts that noise by roughly a factor of sqrt(1000/300) ~ 1.8
* and makes the endpoints stable across seeds. It does NOT narrow the
* intervals: the width is driven by 7-10 treated observations per cell, which
* no amount of resampling can improve. The gain is a reliable answer, not a
* tighter one, and the small-cell caveat in Section 8.2 stands either way.
*
* Cost: runtime scales linearly, so ~3.3x the 300-draw run. This file is the
* only one where cells are thin enough for the endpoint noise to decide
* significance; 08b (500 draws, ~40 treated) and 13c (300, ~18) are not close
* to the boundary and are left as they are.
local nboot  = 1000
local cx     $ctrl_core   // retained for reference; propensity baseline now passes `om' (strict parity)
local cz     l_fedfunds l_reg_crisis_share past_onsets       // Act 1 predictors
* Resolution predictors: both terms are DEFAULT-LINKED-SPECIFIC
* (l_contagion_dist_def, years_since_def_onset) on economic grounds -- a
* predictor for default-linked risk should measure default-linked
* distress/recency, not spread-crisis distress in general. Corroborated by
* 08c_first_stage_table.do's diagnostics; see 08c's header for the full
* argument and the caveat that the reference paper's own instrument is not
* tailored per column.
local cz_def l_fedfunds l_contagion_dist_def years_since_def_onset

* AIPW outcome-model core = the common core, unchanged. The depth term is
* l_credit_bank (WDI FD.AST.PRVT.GD.ZS, credit by banks).
* This block used to swap it for l_credit (FS.AST.PRVT.GD.ZS, all financial
* corporations) on the belief that the by-banks series had the coverage hole.
* 19_sample_audit.do showed the opposite: l_credit_bank is non-missing on EVERY
* row where l_credit is (overlap 1180 = all of l_credit) plus 141 more, and the
* two correlate 0.950. The swap was costing observations, not saving them, so
* the core now carries l_credit_bank everywhere and this file needs no exception.
* Both are plain saved columns, so the bootstrap stays operator-free either way.
local core_aipw l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2

* ══════════════════════════════════════════════════════════════════════════
* AMPLIFIER: pre-crisis sovereign-bank nexus + median split over onsets
* ══════════════════════════════════════════════════════════════════════════
capture drop a_nexus a_nexus_cm highbank
gen a_nexus = L.claimsgov_assets                 // predetermined (year before onset)
bysort cid: egen a_nexus_cm = mean(claimsgov_assets)
replace a_nexus = a_nexus_cm if missing(a_nexus) // country-mean fill (thin coverage)

quietly summarize a_nexus if sample==1 & onset_all==1, detail
local med = r(p50)
gen highbank = (a_nexus >= `med') if !missing(a_nexus)
label define hb 0 "Low nexus" 1 "High nexus"
label values highbank hb

di as result _n "=== SOVEREIGN-BANK NEXUS MEDIAN SPLIT ==="
di as result "  Amplifier = claims on govt / bank assets (pre-crisis, country-mean filled)"
di as result "  Median cutoff among crisis onsets = " %6.2f `med'
foreach t in onset_all onset_nd onset_def {
    quietly count if sample==1 & `t'==1 & highbank==1
    local nh = r(N)
    quietly count if sample==1 & `t'==1 & highbank==0
    local nl = r(N)
    quietly count if sample==1 & `t'==1 & missing(highbank)
    local nm = r(N)
    di as result "  `t': high=" `nh' "  low=" `nl' "  unclassified=" `nm'
}

* ── Confound check: which countries fall in each nexus bin? ──────────────────
*   If the high/low split just proxies income / financial development, the two
*   bins will be dominated by systematically different countries. Print the
*   country composition of the bins (over onsets) and the mean nexus by bin so
*   the split can be inspected, not taken on faith.
di as result _n "=== NEXUS-BIN COUNTRY COMPOSITION (crisis onsets) ==="
capture noisily tabulate country highbank if sample==1 & onset_all==1, ///
    row nofreq
di as result _n "  Mean pre-crisis nexus (claims-on-govt/assets) by bin, over onsets:"
capture noisily tabstat a_nexus if sample==1 & onset_all==1, ///
    by(highbank) statistics(mean min max n) format(%6.2f)
di as result "  (Read alongside any income/development ranking of these countries:"
di as result "   if the two columns are not obviously split by development, the nexus"
di as result "   result is not merely a development proxy.)"

* ══════════════════════════════════════════════════════════════════════════
* CHANNEL OUTCOMES — evolve each channel by high/low nexus (as in 13c)
*   Outcomes: GDP (dy_h, already in panel) + credit, inv, claimpriv_assets,
*   claims_govt (ch_v_h = F h.v - L.v). Every lagged control pre-generated as a
*   PLAIN column so the cluster bootstrap (bsample destroys time order) is valid.
* ══════════════════════════════════════════════════════════════════════════
* OUTCOME SCALE. Strictly-positive GDP-ratio channels use the LOG REAL LEVEL
* (ln_r_*, built in 18_transforms), matching the reference paper's var2/var3:
* a change in X/GDP confounds X with a GDP that is collapsing, whereas
* ln(X/GDP * GDP) = ln(X) up to a constant, so the outcome is the cumulative
* percent change in X itself. pb, fdi and ca change sign so they keep the ratio;
* claimsgov_assets and claimpriv_assets are shares of BANK ASSETS, not of GDP, so
* the denominator problem does not arise for them either.
foreach v in credit inv claimpriv_assets claims_govt {
    local src `v'
    if inlist("`v'","credit","inv") local src ln_r_`v'
    capture drop `v'_base
    gen `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    * (a) pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    *     controls for the channel's own pre-trend momentum in each outcome model.
    capture drop pre_`v'
    gen pre_`v' = L.`src' - L2.`src'
}
foreach v in credit claimpriv_assets claims_govt pb {
    capture drop l_`v'
    gen l_`v' = L.`v'
}

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — _aipw (Eqs. 1-3), returns theta/N/analytic SE
* ══════════════════════════════════════════════════════════════════════════

capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'
    tempvar xb m0 m1 ps summ iwt
    * IPWRA, matching the reference paper: the propensity is estimated FIRST and
    * the outcome regression that produces mu0/mu1 is IPW-WEIGHTED
    * (their `reg g_h dum g_0 $convar [pweight=invwt]`). Estimating it unweighted
    * gives a different, also doubly-robust, estimator - but not theirs.
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `iwt' = `D'/`ps' + (1-`D')/(1-`ps') if `touse'
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' [pweight=`iwt'] if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' [pweight=`iwt'] if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'   // set D=0
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'   // set D=1
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    local th = r(mean)
    local nn = r(N)

    * Analytic (unclustered) influence-function SE -- matches 08b_aipw.do /
    * 13c_aipw_channels.do / 21_aipw_flow.do's _aipw exactly: this is what
    * bands the LEVEL estimates below (no bootstrap on levels).
    tempvar isq
    quietly gen double `isq' = (`summ' - `th')^2 if `touse'
    quietly summarize `isq' if `touse', meanonly
    local sean = sqrt(r(mean)/r(N))

    return scalar theta = `th'
    return scalar N     = `nn'
    return scalar se    = `sean'
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — HIGH - LOW nexus DIFFERENCE within one cell (Asonuma's within-type
*   high-vs-low contrast), aligned with 08b/13c's _aipwpair: levels come from
*   the analytic SE (no bootstrap), and ONLY the difference is bootstrapped,
*   with ROW-LEVEL resampling within control/high/low pools -- the reference
*   paper's own bootstrap device (split into a control pool and one pool per
*   treatment type, `bsample` each, stack), a natural fit here since an onset
*   row already is one episode.
*   _aipwdiff <yvar> <Dvar>, ifch(<high cond>) ifcl(<low cond>) omod() pz() reps()
*   returns r(ok) r(dh) [point high-low] r(bh) r(bl) r(ah) r(al) [analytic SEs]
*           r(se) r(lo) r(hi) r(nd) [bootstrap, difference only]
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipwdiff
program define _aipwdiff, rclass
    syntax anything, IFCH(string) IFCL(string) OMOD(string) PZ(string) REPS(integer)
    gettoken yv Dv : anything
    capture _aipw `yv' `Dv' if `ifch', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bh = r(theta)
    local ah = r(se)
    capture _aipw `yv' `Dv' if `ifcl', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bl = r(theta)
    local al = r(se)
    local dh = `bh' - `bl'

    * Row-level pools: 0 = control (either cell's tranquil rows), 1 = treated
    * in the HIGH-nexus cell, 2 = treated in the LOW-nexus cell.
    capture drop _pool
    quietly gen byte _pool = 0 if (`ifch') | (`ifcl')
    quietly replace _pool = 1 if `Dv' == 1 & (`ifch')
    quietly replace _pool = 2 if `Dv' == 1 & (`ifcl')

    tempname pf
    tempfile bf
    quietly postfile `pf' double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture _aipw `yv' `Dv' if `ifch', omodel(`omod') pmodel(`pz') fe(cid)
            local th = cond(_rc==0, r(theta), .)
            capture _aipw `yv' `Dv' if `ifcl', omodel(`omod') pmodel(`pz') fe(cid)
            local tl = cond(_rc==0, r(theta), .)
            if !missing(`th') & !missing(`tl') quietly post `pf' (`th' - `tl')
        restore
    }
    quietly postclose `pf'
    capture drop _pool
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
    return scalar bh = `bh'
    return scalar bl = `bl'
    return scalar ah = `ah'
    return scalar al = `al'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* ESTIMATE — high/low nexus subsamples, all + resolution split; post results
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str18 outcome str4 part str4 bank byte horizon double b se lo hi ntreat nd ///
    using "`resf'", replace

* (b) second results file: the HIGH - LOW nexus difference per outcome x part x h
tempname D
tempfile diff_resf
postfile `D' str18 outcome str4 part byte horizon double dhl bhi blo se lo hi nd ///
    double cloggz double cloggp using "`diff_resf'", replace

* Outcomes: label | outcome-variable stem (dy or ch_<v>) | channel-specific
*   outcome-model controls (om), same specs as 13c. GDP uses cx (unchanged).
foreach oc in "gdp dy" "credit ch_credit" "inv ch_inv" ///
              "claimpriv_assets ch_claimpriv_assets" "claims_govt ch_claims_govt" {
    gettoken ocl   oc : oc
    gettoken ystem oc : oc

    * outcome-model controls; each channel gets its own pre-crisis change pre_<v>
    * (Asonuma g_0). GDP already carries lagged growth (l1/l2_gdpg) in cx.
    * AIPW outcome core ($core_aipw: common core with total-credit lag); GDP uses
    * the core as-is, channels add their own pre_<v> (credit drops the depth term as its
    * own-level term).
    if      "`ocl'" == "gdp"              local om `core_aipw'
    else if "`ocl'" == "credit"           local om l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
    else                                  local om `core_aipw' pre_`ocl'

    di as result _n "############### OUTCOME: `ocl' ###############"

    * Rows: part-label | treatment dummy | rival dummy to drop | predictors
    *   Part A ("all", onset_all, rival = "" -> none) is SILENCED -- only the
    *   resolution split (nd/def) is of interest now. Restoring it means
    *   adding "all onset_all . cz" back as the first entry below:
    *     foreach cell in "all onset_all . cz" ///
    *                     "nd  onset_nd onset_def cz_def" ///
    *                     "def onset_def onset_nd cz_def" {
    foreach cell in "nd  onset_nd onset_def cz_def" ///
                    "def onset_def onset_nd cz_def" {
        gettoken part cell : cell
        gettoken Dv   cell : cell
        gettoken riv  cell : cell
        gettoken pzn  cell : cell
        * predictor set (indirect: pzn is the local NAME cz or cz_def)
        local pz "``pzn''"

        * Levels (analytic SE, no bootstrap) AND the high-low difference (row
        * bootstrap + Clogg z) come from ONE call to _aipwdiff per horizon --
        * see 08b/13c's header for why this replaces two separate estimation
        * passes with a single, more efficient and internally consistent one.
        di as result _n "--- `ocl' | Part `part' ---"
        di as result "    h   LOW (se_a)        HIGH (se_a)       high-low  [95% boot CI]   Clogg z    p"
        di as result "        LOW/HIGH stars are the conventional t-test vs zero (b/se_a): * p<.10 ** p<.05 *** p<.01."
        post `R' ("`ocl'") ("`part'") ("low")  (0) (0) (0) (0) (0) (0) (0)   // explicit baseline (h=0)
        post `R' ("`ocl'") ("`part'") ("high") (0) (0) (0) (0) (0) (0) (0)
        post `D' ("`ocl'") ("`part'") (0) (0) (0) (0) (0) (0) (0) (0) (.) (.)   // explicit baseline (h=0)

        forvalues h = 0/4 {
            if "`riv'" == "." {
                local ifch sample==1 & (onset_all==0 | (`Dv'==1 & highbank==1))
                local ifcl sample==1 & (onset_all==0 | (`Dv'==1 & highbank==0))
            }
            else {
                local ifch sample==1 & `riv'==0 & (onset_all==0 | (`Dv'==1 & highbank==1))
                local ifcl sample==1 & `riv'==0 & (onset_all==0 | (`Dv'==1 & highbank==0))
            }
            quietly count if `Dv'==1 & highbank==1 & sample==1
            local ntrh = r(N)
            quietly count if `Dv'==1 & highbank==0 & sample==1
            local ntrl = r(N)

            _aipwdiff `ystem'_`h' `Dv', ifch(`ifch') ifcl(`ifcl') ///
                omod(`om') pz(`om' `pz') reps(`nboot')
            if r(ok) {
                local BH = r(bh)   // high-nexus ATE
                local BL = r(bl)   // low-nexus ATE
                local AH = r(ah)   // analytic SE, high
                local AL = r(al)   // analytic SE, low
                local DH = r(dh)
                local SE = r(se)
                local LO = r(lo)
                local HI = r(hi)
                local ND = r(nd)

                * Level CIs = theta +/- 1.96*analytic SE; the bootstrap is
                * reserved for the difference only.
                post `R' ("`ocl'") ("`part'") ("low")  (`h'+1) (`BL') (`AL') (`BL'-1.96*`AL') (`BL'+1.96*`AL') (`ntrl') (.)
                post `R' ("`ocl'") ("`part'") ("high") (`h'+1) (`BH') (`AH') (`BH'-1.96*`AH') (`BH'+1.96*`AH') (`ntrh') (.)

                local zz = .
                local pz2 = .
                if !missing(`AH') & !missing(`AL') & (`AH'^2 + `AL'^2) > 0 {
                    local zz  = `DH' / sqrt(`AH'^2 + `AL'^2)
                    local pz2 = 2*(1 - normal(abs(`zz')))
                }
                post `D' ("`ocl'") ("`part'") (`h'+1) (`DH') (`BH') (`BL') (`SE') (`LO') (`HI') (`ND') (`zz') (`pz2')

                * Conventional t-test for each level vs zero (b / own analytic
                * SE), as the paper's Table 2/Fig. 4 -- same construction as
                * 08b_aipw.do / 13c_aipw_channels.do.
                local tlo  = cond(`AL'>0, `BL'/`AL', .)
                local plo  = cond(!missing(`tlo'), 2*(1-normal(abs(`tlo'))), .)
                local sglo = cond(missing(`plo'), "", cond(`plo'<.01,"***",cond(`plo'<.05,"**",cond(`plo'<.10,"*",""))))
                local thi  = cond(`AH'>0, `BH'/`AH', .)
                local phi  = cond(!missing(`thi'), 2*(1-normal(abs(`thi'))), .)
                local sghi = cond(missing(`phi'), "", cond(`phi'<.01,"***",cond(`phi'<.05,"**",cond(`phi'<.10,"*",""))))

                local sig = cond(`ND'>=50 & !missing(`LO') & (`LO'>0 | `HI'<0), " *", "  ")
                di "    " %1.0f `h'+1 "  " %8.3f `BL' "`sglo'" " (" %5.3f `AL' ")  " ///
                   %8.3f `BH' "`sghi'" " (" %5.3f `AH' ")  " %8.3f `DH' ///
                   " [" %7.3f `LO' ", " %7.3f `HI' "]`sig'" ///
                   " " %7.3f `zz' " " %5.3f `pz2'
            }
            else di as error "    h=" `h'+1 ": estimate failed (too thin)."
        }
    }
}
postclose `R'
postclose `D'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — CSV + figure (panels by resolution part; high vs low nexus lines)
* ══════════════════════════════════════════════════════════════════════════
* (b) HIGH - LOW difference table: dhl = point gap, [lo,hi] = bootstrap 95% CI;
*     a CI excluding 0 means the nexus effect differs significantly by cell.
use "`diff_resf'", clear
label var dhl "AIPW (high - low nexus) difference (pp)"
label var bhi "High-nexus ATE"
label var blo "Low-nexus ATE"
label var lo  "95% CI lower (row bootstrap)"
label var hi  "95% CI upper (row bootstrap)"
label var nd  "Valid bootstrap draws"
label var cloggz "Clogg et al. (1995) z (permissive; assumes independence)"
label var cloggp "p-value of the Clogg z"
gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
label var sig95 "Bootstrap CI excludes 0 (governing test)"
order outcome part horizon dhl bhi blo se lo hi nd sig95 cloggz cloggp
export delimited "$tabs/aipw_nexus_diff.csv", replace
di as result _n "Nexus high-low DIFFERENCE CSV saved: $tabs/aipw_nexus_diff.csv"

use "`resf'", clear
label var b  "AIPW ATE on outcome (pp)"
label var se "Analytic influence-function SE (no bootstrap on levels)"
label var lo "95% CI lower = b - 1.96*se (analytic)"
label var hi "95% CI upper = b + 1.96*se (analytic)"
label var ntreat "Treated onsets in cell"
drop nd
order outcome part bank horizon b se lo hi ntreat
export delimited "$tabs/aipw_nexus_split.csv", replace
di as result _n "Nexus-split AIPW results CSV saved: $tabs/aipw_nexus_split.csv"

* numeric part id for by() panels. Part "all" is SILENCED above (see the
* "all onset_all . cz" note), so `part' now only ever takes nd/def -- no
* partid==1/"All onsets" panel exists any more.
gen byte partid = 1 if part=="nd"
replace partid = 2 if part=="def"
label define pl 1 "Non-default" 2 "Default-linked"
label values partid pl

* ── One high-vs-low figure per outcome (GDP keeps its historical look) ───────
local c_hi "157 36 73"    // high nexus = red (the doom-loop)
local c_lo "0 84 166"     // low  nexus = blue
foreach oc in gdp credit inv claimpriv_assets claims_govt {
    if "`oc'" == "gdp" {
        local ptit "Output cost by sovereign-bank nexus"
        local ytit "Cumulative real GDP change (pp)"
        local fnm  "fig_aipw_nexus_split"       // unchanged filename for GDP
    }
    else {
        local ptit "`oc' channel by sovereign-bank nexus"
        local ytit "Cumulative change in `oc' (pp)"
        local fnm  "fig_nexus_`oc'"
    }
    capture twoway ///
        (rarea lo hi horizon if bank=="high" & outcome=="`oc'", color("`c_hi'%16") lwidth(none)) ///
        (rarea lo hi horizon if bank=="low"  & outcome=="`oc'", color("`c_lo'%16") lwidth(none)) ///
        (connected b horizon if bank=="high" & outcome=="`oc'", lcolor("`c_hi'") lwidth(medthick) msymbol(square) lpattern(dash)) ///
        (connected b horizon if bank=="low"  & outcome=="`oc'", lcolor("`c_lo'") lwidth(medthick) msymbol(circle)), ///
        by(partid, yrescale ///
            note("AIPW (Asonuma Eq. 3), ATE. Split by pre-crisis bank claims-on-govt / assets (median). Shaded = 1.96*analytic SE band. High-low gap bootstrapped directly (row-level).", size(vsmall)) ///
            title("`ptit'", size(medsmall) color(navy))) ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)5) xtitle("Year (Year 1 = crisis year)", size(small)) ///
        ytitle("`ytit'", size(small)) ///
        legend(order(3 "High nexus" 4 "Low nexus") size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    if _rc == 0 {
        graph export "$figs/`fnm'.pdf", replace
        di as result "Figure saved: `fnm'.pdf"
    }
    else di as error "  ** `fnm' failed (rc=" _rc ")"
}

di as result _n "13d_aipw_nexus_split.do complete."
di as result "GDP headline (Part def): high-nexus default crises are the deepest (doom-loop"
di as result "amplifies the default cost). The channel outcomes trace the mechanism: expect"
di as result "claims_govt to RISE under def x high (banks absorbing the sovereign = the loss"
di as result "channel). Read Part B with its small-cell / draw-count caveat."
di as result "aipw_nexus_diff.csv formally tests the sign-flip: a high-low gap whose 95% CI"
di as result "excludes 0 (sig95==1) means the nexus effect differs significantly in that cell."
