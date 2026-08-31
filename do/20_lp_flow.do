/*===========================================================================
  20_LP_FLOW.DO
  The output cost of BEING IN a sovereign spread crisis

  WHY THIS FILE EXISTS
  --------------------
  The headline design (02, 03) treats an episode as a point event: onset_all is
  1 in the first year only, and the estimation sample drops every continuation
  year. That answers "what follows the START of an episode". It does not answer
  the question this project is actually asking — "what is the output cost of
  being in a spread crisis" — because the 173 country-years during which
  countries were actually in their crises are discarded before estimation.

  Here the treatment is the STATE, not the event:

      dy_{i,t+h} = a_i + lambda_t + beta_h * Crisis_it + X_{i,t-1}'d + e

  with Crisis_it = 1 for every year that independently clears the annual
  crisis criterion (in_crisis, built in 18_transforms.do; redefined from
  episode membership to the annual criterion — see that file's header). More
  treated rows than the 61 onsets, though fewer than the old episode-
  membership count, since the years that used to be bridged across a gap
  below the criterion are now ordinary control years — see 18's own printed
  redefinition diagnostic for the exact counts.

  HOW TO READ beta_h — AND WHAT IT IS NOT
  ---------------------------------------
  Horizons follow the paper's convention: Year 1 is the crisis year (dy_0),
  Year 0 is the baseline t-1 and is zero by construction. At Year 1 the outcome
  is 100*(ln_gdp(t) - ln_gdp(t-1)), i.e. growth DURING the crisis year in
  question. The onset row contributes growth in crisis year
  1, the next row growth in crisis year 2, and so on. This is the cleanest
  reading in the file and it is exactly the object the research question names.

  Beyond Year 1 the interpretation changes, and the change must not be papered
  over. beta_h is NOT "the response h years after onset", because the row
  generating it is itself already treated: a row three years into an episode
  contributes a window whose baseline is itself a crisis year. It is the
  cumulative change over the following years associated with being in a crisis
  year, averaged over how long the country has ALREADY been in crisis. Since
  elapsed duration is itself an outcome of the crisis, that average has no fixed
  economic referent. The axis is shared with Tables 1 and 2 so the two designs
  can be compared directly, but the two Year-1 coefficients are NOT the same
  object: the onset one is the first year of every episode, the flow one pools
  the first year of one episode with the fourth year of another.

  WHAT THIS DESIGN DOES NOT IDENTIFY
  ----------------------------------
  The two-stage tier (08, 08b, 13c, 13d) is NOT extended to flow coding, and
  the reason is not computational. For a time-varying treatment the strongest
  predictor of being in crisis at t is being in crisis at t-1 — a post-treatment
  outcome of the same episode — so unconditional unconfoundedness given X is
  false whether or not that term is included. Include it and the first stage
  becomes a persistence model whose scores approach 1 on continuation rows,
  which the [0.01,0.99] trim in 08_ipw_lp.do then deletes: precisely the rows
  this design exists to add. The correct estimator for a time-varying treatment
  is a marginal structural model with inverse-probability-of-treatment weights
  over the treatment HISTORY (Robins, Hernan & Brumback 2000), which is out of
  scope here.

  So onset coding keeps the identification claim, and this file reports
  MAGNITUDE AND STATE CONDITIONAL ON SELECTION INTO CRISIS. Read it as the
  descriptive rung of the ladder, not as a second identified effect.

  NO PRE-TREND PLACEBO IS RUN HERE, deliberately. dy_m2 is
  100*(ln_gdp(t-2) - ln_gdp(t-1)); for a row three or more years into an
  episode BOTH endpoints are crisis years, so regressing it on in_crisis is not
  a placebo but the effect re-estimated with the sign flipped, and it would
  "reject" by construction. The valid restriction is ep_year==1, which IS the
  onset sample — that placebo is already reported in 02 and 03.

  INFERENCE
  ---------
  234 treated rows but still 61 episodes in 52 countries, and a handful of
  chronic cases supply a large share of them. Rows are not information: every
  table below prints episodes and countries alongside N.

  Driscoll-Kraay lag: max(2, h+3). The onset design uses max(1, h+1), which
  covers the h+1 overlap of the outcome windows. Flow coding adds a second
  source of dependence the onset design does not have — the REGRESSOR is
  serially correlated within an episode — so the lag carries +2 for persistence
  at the median episode duration of 2 years. The onset lag rule is reported as
  a second row so that a flow-vs-onset comparison is not confounded by a
  simultaneous change of inference.

  Outputs
  -------
    "$tabs/table9_flow_lp.rtf"        pooled + by resolution type
    "$tabs/flow_lp.csv"               raw coefficients, all variants
    "$clean/irf_flow.dta"             pooled IRF
    "$clean/irf_flow_nd.dta"          non-default IRF
    "$clean/irf_flow_def.dta"         default-linked IRF
    "$figs/fig9_irf_flow.pdf/.png"    pooled figure (built HERE, not in 04)
    "$figs/fig9b_irf_flow_resolution.pdf/.png"   by resolution type
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

capture confirm variable in_crisis, exact
if _rc {
    di as error "  ** in_crisis not in panel_lp.dta — re-run 18_transforms.do first."
    exit 111
}

* ══════════════════════════════════════════════════════════════════════════
* HELPERS (identical to 02_lp_all.do, plus an episode counter)
* ══════════════════════════════════════════════════════════════════════════
capture program drop _critvals
program define _critvals, rclass
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar df  = .
        return scalar c90 = 1.645
        return scalar c95 = 1.960
    }
    else {
        return scalar df  = `dfr'
        return scalar c90 = invttail(`dfr', 0.05)
        return scalar c95 = invttail(`dfr', 0.025)
    }
end

capture program drop _pval
program define _pval, rclass
    args b se
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar p = 2*(1 - normal(abs(`b'/`se')))
    }
    else {
        return scalar p = 2*ttail(`dfr', abs(`b'/`se'))
    }
end

* Rows are not episodes under flow coding. This counts TREATED ROWS, DISTINCT
* EPISODES and COUNTRIES inside e(sample), because a table reporting only N
* would imply far more independent information than the design contains.
capture program drop _nflowcount
program define _nflowcount, rclass
    syntax varname(numeric) , Outcome(varname) Controls(varlist) Samp(varname)
    tempvar esmp
    quietly gen byte `esmp' = e(sample)
    quietly count if `esmp' == 1
    if r(N) == 0 {
        quietly replace `esmp' = (`samp' == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar nrow = r(N)

    tempvar tagep tagcty
    quietly egen byte `tagep' = tag(cid ep_seq) if `varlist'==1 & `esmp'==1
    quietly count if `tagep'==1
    return scalar nep = r(N)

    quietly egen byte `tagcty' = tag(cid) if `varlist'==1 & `esmp'==1
    quietly count if `tagcty'==1
    return scalar ncty = r(N)
end

* ══════════════════════════════════════════════════════════════════════════
* 0. WHAT IS TREATED
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "FLOW TREATMENT — what is treated"
di as result "════════════════════════════════════════════════════════════"
quietly count if onset_all==1 & carryin==0
di as result "  onset years:              " %4.0f r(N)
quietly count if continuation==1 & carryin==0
di as result "  continuation years:       " %4.0f r(N)
quietly count if in_crisis==1
di as result "  TREATED (in_crisis):      " %4.0f r(N)
quietly count if sample_flow==1 & in_crisis==0
di as result "  control (tranquil) years: " %4.0f r(N)
quietly count if sample_flow==1
di as result "  flow sample rows:         " %4.0f r(N)
quietly count if sample==1
di as result "  (onset-design sample:     " %4.0f r(N) ")"

di as result _n "  Treatment is EPISODE MEMBERSHIP: every year from onset to episode end."
di as result "  It is not the annual criterion flag — see the note in 18_transforms.do."

* ══════════════════════════════════════════════════════════════════════════
* 1. POOLED FLOW LP
* ══════════════════════════════════════════════════════════════════════════
* BASELINE CONTROLS ARE EPISODE-DATED ($ctrl_flow, built in 18_transforms.do):
* $ctrl_core evaluated at the episode's entry year and held fixed across its
* years, with tranquil rows keeping their own t-1. Row-dated controls would be
* outcomes of the crisis for any continuation row — mediators on the RHS. This
* dating reproduces the reference paper's timing, whose treated rows are always
* onsets and whose L./L2. controls are therefore predetermined by construction.
* The row-dated set runs as a robustness column in section 4.
if "$ctrl_flow" == "" {
    di as error "  ** $ctrl_flow not set — re-run 18_transforms.do, then this file."
    exit 111
}

* YEAR FIXED EFFECTS ARE KEPT, and 21_aipw_flow.do does not have them.
* METHODOLOGY.md section 2: single-stage LPs carry country AND year FE (02, 03);
* two-stage estimators drop them to stay like-for-like with the reference paper,
* whose $convar carries country dummies (c1-c74) and no year dummies. This file
* is single-stage, so it follows the single-stage rule. Year FE absorb the
* global cycle directly rather than through a proxy, which matters here: the
* year dummies are large and economically real (2020 at -8.95, 2019 at -9.92,
* 2009 at -4.01) and default-linked episodes cluster in 2019-2022, so without
* them COVID and the GFC can load on the treatment.
*
* THE CONSEQUENCE FOR THE FIGURE PAIR, WHICH MUST BE STATED. Figure 9b here and
* Figure 10 in 21 are meant to read as the reference paper's Fig. 3 (OLS) and
* Fig. 4 (AIPW) — the same specification estimated two ways. They are not quite
* that: they differ in the ESTIMATOR and in the year FE simultaneously. The
* no-year-FE robustness row (r_noyearfe, section 4) is what separates the two,
* and the write-up should quote it when setting the figures side by side rather
* than attributing the whole difference to the estimator.
* CONTROL SET. flow_ctrl_variant: 0 = the ADOPTED flow-tier baseline,
* $ctrl_flow, built in 18_transforms.do from $ctrl_core_flowbase
* (l_banking_duration swapped for the l_banking_crisis DUMMY, l_ca swapped
* for tot_chg, tot_chg swapped for exchange2 -- see that file's "ADOPTED
* FLOW-TIER CORE CONTROL SET"). 1 = the adopted core PLUS two further
* candidates, tot_chg (the term exchange2 replaced) and an IMF-program
* dummy l_imf (18_transforms.do's "ALTERNATE FLOW CONTROL SET"). Default 0.
* Because the ADOPTED set is no longer term-for-term identical to
* $ctrl_core, the identity check below (section 2) no longer reproduces
* 02_lp_all.do's published onset-coded coefficient under the default -- it
* is a pure self-consistency check (flow collapses to onset coding under
* whichever control set is active), not a literal match to the published
* figure. `ctrl_row' (and r_rowdated, section 4) also now report the ACTIVE
* variant's row-dated form, not the ORIGINAL $ctrl_core's -- there is no
* longer a built-in bridge back to 02_lp_all.do's exact published number
* anywhere in this file; that comparison belongs in 02_lp_all.do itself.
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
local controls  = cond(`flow_ctrl_variant'==1, "$ctrl_flow_flowplus", "$ctrl_flow")
local ctrl_row  = cond(`flow_ctrl_variant'==1, "$ctrl_core_flowplus", "$ctrl_core_flowbase")

* EXPLORATORY: set to 1 to drop year FE and match the reference paper's
* single-stage rule (country FE only, via c1-c74) instead of this project's
* documented single-stage rule (country AND year FE, METHODOLOGY.md section
* 2 -- see this file's own header for why: year dummies are large and
* economically real, and default-linked episodes cluster in 2019-2022, so
* without them COVID/GFC can load on the treatment). Default 0 = current
* baseline. Drives Sections 1, 3 and the robustness table (section 4).
* Does NOT touch the identity check (section 2), which deliberately runs
* BOTH variants side by side regardless of this setting, or r_noyearfe
* (already the no-year-FE row on its own terms).
local drop_year_fe 0
local yearfe = cond(`drop_year_fe', "", "i.year")

* HORIZON CONVENTION — IDENTICAL TO 02/03 AND TO THE REFERENCE PAPER.
* dy_h is differenced against the row's own t-1, so plugging h=-1 into the same
* formula gives y(t-1) - y(t-1) = 0 by construction, for continuation rows
* exactly as for onset rows. The Year-0 zero is therefore a real normalisation
* here, not a fiction: what changes under flow coding is only what that baseline
* year IS — "the year before this crisis year" rather than "the year before the
* crisis started". So the axis matches Tables 1 and 2: Year 0 = 0, Year 1 = the
* crisis year in question (dy_0), out to Year 5 (dy_4). Figure 9 can be laid
* directly over Figure 1.
*
* Row 1 (displayed Year -1) stays MISSING: the pre-trend placebo is not
* estimable under flow coding (see header), so nothing is written there.
local nhor = 7
foreach m in b se lo90 hi90 lo95 hi95 {
    matrix `m'_flow = J(`nhor', 1, .)
    matrix `m'_flow[2, 1] = 0
}

tempname F
tempfile flowf
postfile `F' str16 spec str14 term int hdisp double b double se double p ///
    long nrow long nep long ncty long N using "`flowf'", replace

di as result _n "════════════════════════════════════════════════════════════"
di as result "1. POOLED — output while in a spread crisis (Year 1 = the crisis year)"
di as result "════════════════════════════════════════════════════════════"
di as result "  Year 1 is growth DURING the crisis year. Later years are the"
di as result "  cumulative change from t-1 associated with a crisis year — NOT"
di as result "  'years after onset': the row generating it is itself treated."
di as result ""

forvalues h = 0/4 {
    local hd  = `h' + 1
    local row = `h' + 3
    local lag = max(2, `h' + 3)

    capture noisily xtscc dy_`h' in_crisis `controls' `yearfe' ///
        if sample_flow == 1, fe lag(`lag')
    if _rc {
        di as error "  ** h=`hd' failed (rc=" _rc ")"
        continue
    }

    local bb = _b[in_crisis]
    local ss = _se[in_crisis]
    _pval `bb' `ss'
    local pp = r(p)
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)

    _nflowcount in_crisis, outcome(dy_`h') controls(`controls') samp(sample_flow)
    local nr = r(nrow)
    local ne = r(nep)
    local nc = r(ncty)

    eststo f1_h`h'
    estadd scalar nrow = `nr'
    estadd scalar nep  = `ne'
    estadd scalar ncty = `nc'

    matrix b_flow[`row',1]    = `bb'
    matrix se_flow[`row',1]   = `ss'
    matrix lo90_flow[`row',1] = `bb' - `c90'*`ss'
    matrix hi90_flow[`row',1] = `bb' + `c90'*`ss'
    matrix lo95_flow[`row',1] = `bb' - `c95'*`ss'
    matrix hi95_flow[`row',1] = `bb' + `c95'*`ss'

    di as result "  h=" %1.0f `hd' "   b = " %7.3f `bb' "   SE = " %6.3f `ss' ///
                 "   p = " %5.3f `pp' "   rows = " %4.0f `nr' ///
                 "   episodes = " %3.0f `ne' "   countries = " %3.0f `nc'

    post `F' ("pooled") ("in_crisis") (`hd') (`bb') (`ss') (`pp') (`nr') (`ne') (`nc') (e(N))
}

* ══════════════════════════════════════════════════════════════════════════
* 2. THE IDENTITY CHECK
*
* Restricted to sample==1 — which drops every continuation row — in_crisis IS
* onset_all, and epc_* collapses to the row-dated controls. So the flow
* specification and the onset specification must return the SAME coefficient on
* that sample. This is the check that catches a mis-built treatment or sample
* flag. It must run on sample==1, NOT sample_flow==1: on the flow sample the two
* legitimately differ, because the control group is no longer the same.
*
* TWO VERSIONS ARE RUN, both testing the SAME thing: that flow-coded in_crisis
* on sample==1 (which drops every continuation row) reproduces onset_all's
* coefficient under the ACTIVE control set (`controls'/`ctrl_row', currently
* the adopted $ctrl_flow/$ctrl_core_flowbase by default). This is now a pure
* internal self-consistency check, NOT a literal match to 02_lp_all.do's
* published Table 1 coefficient: since the flow tier's adopted core control
* set (18_transforms.do's $ctrl_core_flowbase) swaps l_banking_duration for
* l_banking_crisis and l_ca for tot_chg, it is no longer term-for-term
* identical to 02's $ctrl_core. Both variants are still reported so a
* failure localises to either the treatment/sample construction (variant i)
* or something FE-specific (variant ii).
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "2. IDENTITY CHECK — flow Year 1 on sample==1 must equal onset Year 1"
di as result "════════════════════════════════════════════════════════════"
quietly count if sample==1 & in_crisis==1 & onset_all==0
if r(N) != 0 {
    di as error "  ** `r(N)' rows have in_crisis==1 & onset_all==0 inside sample==1."
    di as error "     sample is supposed to exclude every continuation year — check 18."
}
* On sample==1 the episode-dated controls must be IDENTICAL to the row-dated
* ones: continuation rows are excluded, so in_crisis==1 only on onset rows,
* where the entry year IS t-1. Checking that first isolates a failure to the
* right place if the coefficient test below fails.
local nredate = 0
foreach X of global ctrl_core_flowbase {
    capture confirm variable epc_`X', exact
    if _rc continue
    quietly count if sample==1 & !missing(`X') & abs(epc_`X' - `X') > 1e-9
    local nredate = `nredate' + r(N)
}
if `nredate' != 0 di as error "  ** `nredate' rows where epc_ differs from row-dated INSIDE sample==1 — re-dating leaked"
else              di as result "  Episode-dated controls collapse to row-dated on sample==1 (correct)"

* (i) no year FE — the baseline specification's own internal equivalence
quietly xtscc dy_0 in_crisis `controls' if sample==1, fe lag(1)
local b_flow0 = _b[in_crisis]
quietly xtscc dy_0 onset_all `ctrl_row' if sample==1, fe lag(1)
local b_ons0 = _b[onset_all]
local gap0 = abs(`b_flow0' - `b_ons0')
di as result "  no year FE   flow: " %9.6f `b_flow0' "   onset: " %9.6f `b_ons0'

* (ii) with year FE — flow vs onset under the SAME active control set. No
* longer comparable to 02_lp_all.do's published number under the default
* (adopted) control set — see the section header.
quietly xtscc dy_0 in_crisis `controls' i.year if sample==1, fe lag(1)
local b_flow0y = _b[in_crisis]
quietly xtscc dy_0 onset_all `ctrl_row' i.year if sample==1, fe lag(1)
local b_ons0y = _b[onset_all]
local gap0y = abs(`b_flow0y' - `b_ons0y')
di as result "  with year FE flow: " %9.6f `b_flow0y' "   onset: " %9.6f `b_ons0y'
di as result "               (self-consistency check under the ACTIVE control set;"
di as result "                not a literal match to 02_lp_all.do's published Table 1"
di as result "                since none of the flow_ctrl_variant options point ctrl_row at \$ctrl_core exactly)"
if `gap0y' > 1e-6 di as error "  ** with-year-FE anchor FAILED (difference " %9.6f `gap0y' ")."
if `gap0' > 1e-6 {
    di as error "  ** IDENTITY CHECK FAILED (difference " %9.6f `gap0' ")."
    di as error "     in_crisis or sample_flow is mis-built. Stop and fix before reading anything below."
}
else di as result "  MATCH (difference " %9.2e `gap0' ") — treatment, sample flag and control dating are all correct."

* ══════════════════════════════════════════════════════════════════════════
* 3. BY RESOLUTION TYPE — joint specification
*
* Both type dummies entered together with tranquil years as the omitted
* category, mirroring Spec B in 03_lp_resolution.do (the headline design there).
*
* The difference is tested with LINCOM, not with the Clogg et al. (1995) z used
* at 03_lp_resolution.do:365-372. That z treats the two coefficients as
* independent and adds their variances; inside a single joint regression they
* are not independent, and under flow coding they are less independent still,
* because a country with both a non-default and a default-linked episode
* contributes treated rows to BOTH arms. lincom uses the estimated covariance
* and is the correct standard error for the difference. (03's Wald test two
* lines above its Clogg z is the covariance-correct number there; the z is a
* redundant second statistic.)
* ══════════════════════════════════════════════════════════════════════════
foreach g in nd def {
    foreach m in b se lo90 hi90 lo95 hi95 {
        matrix `m'_flow_`g' = J(`nhor', 1, .)
        matrix `m'_flow_`g'[2, 1] = 0
    }
}
matrix pdiff_flow = J(5, 1, .)

di as result _n "════════════════════════════════════════════════════════════"
di as result "3. BY RESOLUTION TYPE — joint, tranquil omitted"
di as result "════════════════════════════════════════════════════════════"

forvalues h = 0/4 {
    local hd  = `h' + 1
    local row = `h' + 3
    local lag = max(2, `h' + 3)

    capture noisily xtscc dy_`h' in_crisis_nd in_crisis_def `controls' `yearfe' ///
        if sample_flow == 1, fe lag(`lag')
    if _rc {
        di as error "  ** h=`hd' split failed (rc=" _rc ")"
        continue
    }

    local bnd  = _b[in_crisis_nd]
    local snd  = _se[in_crisis_nd]
    local bdef = _b[in_crisis_def]
    local sdef = _se[in_crisis_def]
    _critvals
    local c90 = r(c90)
    local c95 = r(c95)
    _pval `bnd' `snd'
    local pnd = r(p)
    _pval `bdef' `sdef'
    local pdef = r(p)

    _nflowcount in_crisis_nd, outcome(dy_`h') controls(`controls') samp(sample_flow)
    local nrnd = r(nrow)
    local nend = r(nep)
    _nflowcount in_crisis_def, outcome(dy_`h') controls(`controls') samp(sample_flow)
    local nrdef = r(nrow)
    local nedef = r(nep)

    * Covariance-correct difference.
    quietly lincom in_crisis_def - in_crisis_nd
    local bdiff = r(estimate)
    local sdiff = r(se)
    local dfd   = r(df)
    if !missing(`dfd') & `dfd' > 0 local pdi = 2*ttail(`dfd', abs(`bdiff'/`sdiff'))
    else                           local pdi = 2*(1 - normal(abs(`bdiff'/`sdiff')))
    matrix pdiff_flow[`h'+1, 1] = `pdi'

    eststo f2_h`h'
    estadd scalar nepnd  = `nend'
    estadd scalar nepdef = `nedef'
    estadd scalar nrownd = `nrnd'
    estadd scalar nrowdf = `nrdef'
    estadd scalar bdiff  = `bdiff'
    estadd scalar pdiff  = `pdi'

    foreach g in nd def {
        matrix b_flow_`g'[`row',1]    = `b`g''
        matrix se_flow_`g'[`row',1]   = `s`g''
        matrix lo90_flow_`g'[`row',1] = `b`g'' - `c90'*`s`g''
        matrix hi90_flow_`g'[`row',1] = `b`g'' + `c90'*`s`g''
        matrix lo95_flow_`g'[`row',1] = `b`g'' - `c95'*`s`g''
        matrix hi95_flow_`g'[`row',1] = `b`g'' + `c95'*`s`g''
    }

    di as result "  h=" %1.0f `hd' ///
        "   ND = "  %7.3f `bnd'  " (" %5.3f `pnd'  ", " %3.0f `nend'  " ep)" ///
        "   DEF = " %7.3f `bdef' " (" %5.3f `pdef' ", " %3.0f `nedef' " ep)" ///
        "   diff = " %7.3f `bdiff' "  p = " %5.3f `pdi'

    post `F' ("split") ("in_crisis_nd")  (`hd') (`bnd')  (`snd')  (`pnd')  (`nrnd')  (`nend')  (.) (e(N))
    post `F' ("split") ("in_crisis_def") (`hd') (`bdef') (`sdef') (`pdef') (`nrdef') (`nedef') (.) (e(N))
    post `F' ("split") ("def_minus_nd")  (`hd') (`bdiff') (`sdiff') (`pdi') (.) (.) (.) (e(N))
}

* ══════════════════════════════════════════════════════════════════════════
* 4. ROBUSTNESS
*
* Each variant states the concern it addresses, so that a reader knows what
* would have counted as failure.
*
* (0) YEAR FE. The baseline KEEPS them, following the project's rule for a
*     single-stage LP (METHODOLOGY.md section 2). 21_aipw_flow.do cannot carry
*     them, so the no-year-FE version runs here for both the pooled coefficient
*     and the def-nd difference, and is the row to read before setting Figure 9b
*     beside Figure 10.
*
* (a) INFERENCE. The onset design's DK lag, max(1,h+1), covers the overlap of
*     the outcome windows but not the serial correlation of the REGRESSOR
*     within an episode, which flow coding introduces and onset coding does not
*     have. Reporting the old rule alongside the new one means a flow-vs-onset
*     comparison is not confounded by a simultaneous change of inference.
*     Also reported: country-clustered SEs, on 52 clusters.
*
* (b) TREATMENT DEFINITION. The baseline `in_crisis' is now the ANNUAL
*     CRITERION (crisis_any==1 that year; redefined in 18_transforms.do). This
*     row reports in_crisis_episode, the OLD episode-membership definition
*     (onset_all|continuation), which bridges the 13 mid-episode years that fail
*     the annual criterion back into the treated set instead of leaving them as
*     controls. Reported so the choice — and its cost — stays visible; see
*     18_transforms.do's header for the full redefinition rationale and the
*     exact row-count delta it prints.
*     Also: the gap years dropped entirely (neither treated nor control), to
*     test whether using them as ordinary controls (the new baseline's choice)
*     matters versus excluding them from the sample outright.
*
* (c) CONTROL DATING. The baseline dates $ctrl_core at the episode's entry year.
*     The row-dated version conditions every continuation row on covariates the
*     crisis itself produced — lagged growth most obviously, but debt, the
*     current account, bank credit, government spending and banking duration
*     equally — which will absorb part of the effect being measured and bias
*     beta toward zero. If the two agree, one sentence closes the question; if
*     they disagree, the gap measures how much the row-dated controls absorb.
*     A row-dated-minus-l1_gdpg variant is reported as well, since lagged growth
*     is the single worst offender.
*
* (d) CONCENTRATION. Flow coding weights episodes by their length, so a few
*     chronic cases carry a large share of the treatment. Venezuela alone is in
*     crisis in most of its panel years. If the result rests on it, that must
*     be visible rather than discovered by a referee.
*
* (e) FORECAST DATA. The panel runs to 2026 while realised national accounts end
*     in 2024-2025. This matters MORE under flow coding than for onsets, because
*     long episodes run into the projection window at every horizon.
*
* (f) EXCHANGE-RATE SWAP — SETTLED, DIAGNOSTIC CODE REMOVED. tot_chg
*     (terms-of-trade log-change) has been REPLACED by exchange2 (log FX
*     change) inside $ctrl_core_flowbase itself (18_transforms.do, "THE
*     FOURTH SWAP") -- no longer a robustness row here, because it is now
*     the baseline every row in this section already runs under. The
*     decision was made on LITERATURE grounds (the reference paper's own
*     $convar carries exchange-rate depreciation as a baseline control;
*     terms-of-trade was never in their $convar -- METHODOLOGY.md section 4),
*     not on the empirical test this section used to run: that test found
*     epc_exchange2 individually significant at only 1 of 5 horizons (h=5,
*     p=.037, with year FE) where epc_tot_chg had been significant at 4 of 5
*     (p<.05 at h=2,3,4,5), and found dropping tot_chg shifted in_crisis
*     systematically more negative at every horizon -- evidence AGAINST the
*     swap on pure outcome-equation grounds, overridden here by the
*     literature-fidelity argument. Both sides of that tradeoff are recorded
*     in 18_transforms.do's header, not just the side that was acted on. See
*     git history if the diagnostic code itself is ever needed again.
* ══════════════════════════════════════════════════════════════════════════
local dropgdpg l1_gdpg
local cflow : list ctrl_row - dropgdpg

di as result _n "════════════════════════════════════════════════════════════"
di as result "4. ROBUSTNESS (all horizons written to the CSV)"
di as result "════════════════════════════════════════════════════════════"

forvalues h = 0/4 {
    local hd = `h' + 1

    * (0) NO YEAR FIXED EFFECTS — the specification 21_aipw_flow.do runs.
    *     The baseline here keeps year FE, which is the project's standard for a
    *     single-stage LP (02, 03) and absorbs the global cycle directly. 21
    *     cannot carry them: it is two-stage, and the reference paper's $convar
    *     has country dummies only. So Figure 9b and Figure 10 differ in the
    *     ESTIMATOR and in the year FE at the same time, and this row is what
    *     lets a reader separate the two. Read it against the baseline before
    *     comparing the two figures.
    *
    *     Expect a large gap. The year dummies are economically real (2020 at
    *     -8.95, 2019 at -9.92, 2009 at -4.01 in the year-FE run) and
    *     default-linked episodes cluster in 2019-2022, so dropping them lets
    *     global shocks load on the treatment. METHODOLOGY.md records the onset
    *     analogue: default-linked moves from -5.53/-6.75 to -7.58/-9.48.
    capture quietly xtscc dy_`h' in_crisis `controls' if sample_flow==1, fe lag(`=max(2,`h'+3)')
    if _rc == 0 post `F' ("r_noyearfe") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))

    capture quietly xtscc dy_`h' in_crisis_nd in_crisis_def `controls' if sample_flow==1, fe lag(`=max(2,`h'+3)')
    if _rc == 0 {
        quietly lincom in_crisis_def - in_crisis_nd
        post `F' ("r_noyearfe") ("def_minus_nd") (`hd') (r(estimate)) (r(se)) ///
            (cond(!missing(r(df)) & r(df)>0, 2*ttail(r(df), abs(r(estimate)/r(se))), ///
                  2*(1-normal(abs(r(estimate)/r(se)))))) (.) (.) (.) (e(N))
    }

    * (a) onset-design lag, then country clustering
    capture quietly xtscc dy_`h' in_crisis `controls' `yearfe' if sample_flow==1, fe lag(`=max(1,`h'+1)')
    if _rc == 0 post `F' ("r_oldlag") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))

    capture quietly areg dy_`h' in_crisis `controls' `yearfe' if sample_flow==1, absorb(cid) vce(cluster cid)
    if _rc == 0 post `F' ("r_cluster") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*ttail(e(df_r), abs(_b[in_crisis]/_se[in_crisis]))) (.) (.) (e(N_clust)) (e(N))

    * (b) treatment definition — OLD episode-membership definition, robustness/
    *     historical comparison against the new annual-criterion baseline.
    capture quietly xtscc dy_`h' in_crisis_episode `controls' `yearfe' if sample_flow==1, fe lag(`=max(2,`h'+3)')
    if _rc == 0 post `F' ("r_episodedef") ("in_crisis_episode") (`hd') (_b[in_crisis_episode]) (_se[in_crisis_episode]) ///
        (2*(1-normal(abs(_b[in_crisis_episode]/_se[in_crisis_episode])))) (.) (.) (.) (e(N))

    capture quietly xtscc dy_`h' in_crisis `controls' `yearfe' if sample_flow==1 & gap_year_dropped==0, fe lag(`=max(2,`h'+3)')
    if _rc == 0 post `F' ("r_nogap") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))

    * (c) row-dated controls, and row-dated minus lagged growth
    capture quietly xtscc dy_`h' in_crisis `ctrl_row' `yearfe' if sample_flow==1, fe lag(`=max(2,`h'+3)')
    if _rc == 0 post `F' ("r_rowdated") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))

    capture quietly xtscc dy_`h' in_crisis `cflow' `yearfe' if sample_flow==1, fe lag(`=max(2,`h'+3)')
    if _rc == 0 post `F' ("r_rowdated_nog") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))

    * (d) concentration
    capture quietly xtscc dy_`h' in_crisis `controls' `yearfe' if sample_flow==1 & country!="Venezuela", fe lag(`=max(2,`h'+3)')
    if _rc == 0 post `F' ("r_noven") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
        (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))

    * (e) realised outcomes only
    capture confirm variable gdp_last_actual, exact
    if !_rc {
        capture quietly xtscc dy_`h' in_crisis `controls' `yearfe' ///
            if sample_flow==1 & (year + `h') <= gdp_last_actual, fe lag(`=max(2,`h'+3)')
        if _rc == 0 post `F' ("r_outturn") ("in_crisis") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
            (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))
    }

}
di as result "  (variants estimated; see $tabs/flow_lp.csv)"

postclose `F'

* ══════════════════════════════════════════════════════════════════════════
* 5. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
capture esttab f1_h0 f1_h1 f1_h2 f1_h3 f1_h4 using "$tabs/table9_flow_lp.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(in_crisis) coeflabel(in_crisis "In a spread crisis") ///
    mtitles("Year 1" "Year 2" "Year 3" "Year 4" "Year 5") nonumber ///
    stats(nep ncty nrow N, labels("Episodes" "Countries" "Treated country-years" "Observations") fmt(0 0 0 0)) ///
    title("Table 9. Output while in a sovereign spread crisis (flow specification)") ///
    addnotes("Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h." ///
             "Treatment = 1 in EVERY year of an episode, onset and continuation alike (234 country-years)." ///
             "Horizons follow Tables 1-2: Year 0 is the baseline t-1 (zero by construction), Year 1 is the" ///
             "crisis year. Year 1 here is growth DURING a crisis year, but it pools the first year of one" ///
             "episode with the fourth year of another, so it is NOT the same object as Year 1 in Table 1." ///
             "Later years are cumulative changes from t-1, not responses h years after onset." ///
             "Controls are dated at the EPISODE's entry year and held fixed across its years; tranquil" ///
             "rows keep their own t-1. Row-dated controls would be outcomes of the crisis for any" ///
             "continuation row. Onset rows are unaffected, so this matches Table 1's control set at Year 1." ///
             "Country AND year fixed effects, the project's standard for a single-stage LP." ///
             "21_aipw_flow.do is two-stage and carries country FE only, matching the reference" ///
             "paper, so Figure 9b and Figure 10 differ in the estimator AND in the year FE." ///
             "Robustness row r_noyearfe drops them here, which is what separates the two." ///
             "Driscoll-Kraay SE, lag max(2,h+3): h+1 for the overlap of" ///
             "outcome windows plus 2 for serial correlation of the treatment within an episode." ///
             "Report episodes and countries, not rows: 234 treated rows carry 61 episodes in 52 countries." ///
             "This specification does not correct for selection into crisis — see METHODOLOGY.md." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")
if _rc == 608 di as error "  ** table9_flow_lp.rtf is OPEN IN WORD — close it and re-run."
else if _rc   di as error "  ** Table 9: esttab failed (rc=" _rc ")"
else          di as result "Table 9 saved: $tabs/table9_flow_lp.rtf"

capture esttab f2_h0 f2_h1 f2_h2 f2_h3 f2_h4 using "$tabs/table9_flow_lp.rtf", append ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(in_crisis_nd in_crisis_def) ///
    coeflabel(in_crisis_nd "In a non-default crisis" in_crisis_def "In a default-linked crisis") ///
    mtitles("Year 1" "Year 2" "Year 3" "Year 4" "Year 5") nonumber ///
    stats(bdiff pdiff nepnd nepdef N, ///
          labels("Difference (def - nd)" "p-value of difference (lincom)" ///
                 "Episodes, non-default" "Episodes, default-linked" "Observations") fmt(3 3 0 0 0)) ///
    title("Table 9b. By resolution type, joint specification (tranquil omitted)") ///
    addnotes("Both treatment dummies entered jointly; tranquil country-years are the omitted category." ///
             "The difference and its p-value come from lincom, which uses the estimated covariance" ///
             "between the two coefficients. A Clogg et al. (1995) z would treat them as independent," ///
             "which they are not: countries with episodes of both types contribute to both arms.")
if _rc == 608 di as error "  ** table9_flow_lp.rtf is OPEN IN WORD — close it and re-run."
else if _rc   di as error "  ** Table 9b: esttab failed (rc=" _rc ")"
else          di as result "Table 9b appended: $tabs/table9_flow_lp.rtf"

preserve
    use "`flowf'", clear
    label var spec "Specification variant"
    label var hdisp "Horizon h (0 = growth during the crisis year)"
    label var nep  "Distinct episodes contributing"
    label var ncty "Distinct countries contributing"
    label var nrow "Treated country-years"
    export delimited "$tabs/flow_lp.csv", replace
    di as result "Raw coefficients saved: $tabs/flow_lp.csv"
restore

* ── IRF datasets for 04_graphs.do ───────────────────────────────────────────
foreach g in flow flow_nd flow_def {
    preserve
        clear
        set obs `nhor'
        gen horizon = _n - 2
        foreach m in b se lo90 hi90 lo95 hi95 {
            svmat `m'_`g', names(`m')
            rename `m'1 `m'
        }
        gen series = "`g'"
        save "$clean/irf_`g'.dta", replace
    restore
}
di as result "IRF datasets saved: irf_flow.dta, irf_flow_nd.dta, irf_flow_def.dta"

* ══════════════════════════════════════════════════════════════════════════
* 6. FIGURES
*
* Drawn here rather than in 04_graphs.do so this file is self-contained: run it
* alone and you get the tables, the CSV, the IRF datasets and the figures.
* (02/03 hand their IRFs to 04 instead; the flow tier does not, so there is
* only one place where Figures 9 and 9b are defined.)
*
* SAME AXIS AS FIGURE 1 (Year 0 = baseline t-1, Year 1 = the crisis year), so
* the two can be laid over each other — but the Year-1 points are NOT the same
* object: Figure 1's is the first year of every episode; this one pools the
* first year of one episode with the fourth year of another, because the
* treatment is being IN a crisis, not entering one.
*
* Colours repeat 04_graphs.do:15-17 so the flow figures sit visually alongside
* Figures 1-2. Everything runs inside preserve/restore and every export is
* wrapped in capture: a locked file or a missing $figs must not discard
* estimates that are already computed.
* ══════════════════════════════════════════════════════════════════════════
local c_flow "23 55 94"
local c_nd   "0 84 166"
local c_def  "157 36 73"
local c_zero "150 150 150"

preserve
    use "$clean/irf_flow.dta", clear
    keep if horizon >= 0
    twoway ///
        (rarea lo95 hi95 horizon, color("`c_flow'%15") lwidth(none)) ///
        (rarea lo90 hi90 horizon, color("`c_flow'%25") lwidth(none)) ///
        (connected b horizon, lcolor("`c_flow'") lwidth(medthick) ///
            mcolor("`c_flow'") msize(medium) msymbol(circle)), ///
        yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
        xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year (Year 1 = the crisis year)", size(medsmall)) ///
        ytitle("Change in log real GDP (pp)", size(medsmall)) ///
        title("Output While In a Sovereign Spread Crisis", size(medium)) ///
        subtitle("Treatment = every year of an episode (234 country-years, 61 episodes)", size(small)) ///
        note("Local projections, country and year FE, Driscoll-Kraay SE (lag max(2,h+3))." ///
             "Controls are the common core dated at the EPISODE's entry year and held fixed" ///
             "across its years; tranquil rows keep their own t-1." ///
             "Same axis as Figure 1, but Year 1 is not the same object: it pools the first year of" ///
             "one episode with the fourth year of another. Not corrected for selection.", ///
             size(vsmall)) ///
        legend(order(3 "Point estimate" 2 "90% CI" 1 "95% CI") ring(0) pos(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    capture graph export "$figs/fig9_irf_flow.pdf", replace
    if _rc di as error "  ** fig9_irf_flow.pdf export failed (rc=" _rc ") — is it open?"
    else {
        capture graph export "$figs/fig9_irf_flow.png", replace width(1200)
        di as result "Figure 9 saved: fig9_irf_flow.pdf/.png"
    }
restore

* ── Figure 9b — the OLS half of the Fig.3 / Fig.4 pair ──────────────────────
* `yrng' sets a COMMON y-scale so this and fig10_aipw_flow.pdf can be shown
* together. Left empty until both files have run and the AIPW magnitudes are
* known; then set the SAME string here and in 21_aipw_flow.do, e.g.
*     local yrng "yscale(range(-12 6)) ylabel(-12(3)6, format(%4.1f) labsize(medsmall))"
local yrng ""

capture confirm file "$clean/irf_flow_nd.dta"
if _rc {
    di as error "  ** irf_flow_nd.dta not found — Figure 9b skipped."
}
else {
    preserve
        use "$clean/irf_flow_nd.dta", clear
        append using "$clean/irf_flow_def.dta"
        keep if horizon >= 0
        twoway ///
            (rarea lo95 hi95 horizon if series=="flow_nd",  color("`c_nd'%20")  lwidth(none)) ///
            (rarea lo95 hi95 horizon if series=="flow_def", color("`c_def'%20") lwidth(none)) ///
            (connected b horizon if series=="flow_nd",  lcolor("`c_nd'")  mcolor("`c_nd'")  ///
                msymbol(circle) lwidth(medthick)) ///
            (connected b horizon if series=="flow_def", lcolor("`c_def'") mcolor("`c_def'") ///
                msymbol(square) lpattern(dash) lwidth(medthick)), ///
            yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
            xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
            xtitle("Year (Year 1 = the crisis year)", size(medsmall)) ///
            ytitle("Change in log real GDP (pp)", size(medsmall)) ///
            title("Output While In a Spread Crisis, by Resolution", size(medium)) ///
            subtitle("OLS — the reference paper's Eq. (1). Tranquil country-years omitted.", size(small)) ///
            `yrng' ///
            note("95% CIs. Country and year FE. The OLS half of a pair with fig10_aipw_flow.pdf" ///
                 "(21_aipw_flow.do), the analogues of the reference paper's Fig. 3 and Fig. 4." ///
                 "NOTE: 21 is two-stage and carries country FE only, so the two panels differ in the" ///
                 "estimator AND in the year FE. Robustness row r_noyearfe separates them." ///
                 "Difference tested by lincom, using the covariance between the two coefficients." ///
                 "Non-default carries 113 crisis-years, default-linked 121.", size(vsmall)) ///
            legend(order(3 "Non-default" 4 "Default-linked") ring(0) pos(7) cols(1) size(small)) ///
            graphregion(color(white)) plotregion(color(white))
        capture graph export "$figs/fig9b_irf_flow_resolution.pdf", replace
        if _rc di as error "  ** fig9b export failed (rc=" _rc ") — is it open?"
        else {
            capture graph export "$figs/fig9b_irf_flow_resolution.png", replace width(1200)
            di as result "Figure 9b saved: fig9b_irf_flow_resolution.pdf/.png"
        }
    restore
}

di as result _n "20_lp_flow.do complete."
