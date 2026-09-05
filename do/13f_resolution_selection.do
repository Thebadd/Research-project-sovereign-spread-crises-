/*===========================================================================
  13F_RESOLUTION_SELECTION.DO
  Does pre-crisis bank exposure to the sovereign predict HOW a spread crisis
  is resolved?

  THE QUESTION, AND WHY NOTHING ELSE IN THE PROJECT ANSWERS IT
  ------------------------------------------------------------
  Every propensity model in this project estimates

        Pr(onset_s = 1 | X_{t-1})       s in {all, nd, def}

  against TRANQUIL years, with the rival type dropped from the sample
  (08b_aipw.do, 13d). That is selection INTO A CRISIS. This
  file estimates a different object entirely — selection INTO A RESOLUTION,

        Pr(non-default | spread-crisis onset, nexus_{t-1}, X_{t-1})

  i.e. conditional on a crisis having started, does the pre-crisis share of
  the banking system's balance sheet committed to the sovereign predict
  whether the sovereign services its debt or defaults?

  WHAT WE ALREADY KNOW GOING IN
  -----------------------------
  02a_descriptive_facts.do reports claimsgov_assets at t-1 as 15.53% of bank
  assets for non-default onsets and 12.56% for default-linked, p = 0.554.
  The unconditional difference points toward high exposure and non-default
  resolution, and is nowhere near significant on 61 episodes. This file asks
  whether conditioning helps. The prior is therefore a NULL, not a result to
  be recovered.

  THE TWO READINGS, STATED BEFORE ESTIMATION
  ------------------------------------------
  Per CLAUDE.md, both signs are defined in advance and the data pick.

    b > 0  (high exposure -> NON-default)
      Domestic bank holdings of government debt as a COMMITMENT DEVICE:
      defaulting destroys the domestic banking system, so a sovereign whose
      banks are heavily loaded raises its own cost of default and muddles
      through instead. Gennaioli, Martin & Rossi (2014, JF); Broner, Erce,
      Martin & Ventura (2014, JME).

    b < 0  (high exposure -> DEFAULT)
      Exposure as FRAGILITY: a banking system saturated with sovereign paper
      is a weak banking system, the sovereign cannot lean on it to roll over,
      and default follows anyway. The doom-loop reading.

  Both are published; neither is refuted by the other's sign. A null is the
  third admissible outcome and, given the counts below, the most likely one.

  THE BINDING CONSTRAINT: COUNTS
  ------------------------------
  61 onsets, less the IMF depository-corporations coverage limit (2001+, and
  six countries absent from the ODC file entirely: China, Ecuador,
  El Salvador, India, Lebanon, Vietnam), less listwise deletion on controls.
  The realised estimation sample is printed at every column, together with
  the number of DEFAULT EVENTS, because that — not N — is what binds a binary
  model. At the conventional ten-events-per-parameter rule, roughly fifteen
  default events support one regressor plus one or two controls and no more.
  The specification ladder below is therefore short by design and stops
  before it becomes arithmetic rather than estimation.

  Because the design cannot be made larger, this file reports the MINIMUM
  DETECTABLE EFFECT alongside every estimate: the smallest marginal effect
  that would have cleared 5% given the standard error actually obtained. A
  null accompanied by an MDE of half a probability point is informative; a
  null accompanied by an MDE of thirty points is a statement about the
  sample, not about the world. The distinction is the whole content of this
  file if the coefficient does not clear significance.

  WHAT THIS CANNOT SHOW, UNDER ANY OUTCOME
  ----------------------------------------
  This is an association between a predetermined variable and a resolution
  outcome across roughly forty-five episodes. It cannot establish that the
  nexus CAUSES governments to avoid default. Exposure measured at t-1 is
  predetermined with respect to the onset, which rules out the crisis causing
  the exposure, but not a third factor causing both — and specifically not
  the possibility that banks were already accumulating sovereign paper in
  anticipation of the distress that determines the resolution. Section 4
  below probes that directly by re-running the ladder on exposure measured at
  t-2 and t-3.

  THE COMPLICATION FOR SECTION 8 OF THE WRITE-UP
  ----------------------------------------------
  13b finds that WITHIN default-linked episodes, high pre-crisis exposure
  deepens the output cost (d_def = -4.43 at h3, -4.27 at h4). If this file
  finds that high exposure also selects sovereigns AWAY from default, the two
  results are not additive: the high-exposure episodes that defaulted anyway
  would be a selected and unusual group, and the amplification estimate is
  conditional on that selection. That is a genuine interpretive cost and it
  belongs in the write-up whichever way the coefficient falls. It is a reason
  to run this test, not a reason to avoid it.

  Outputs
  -------
    "$tabs/table13f_resolution_selection.rtf"   Word: the specification ladder
    "$tabs/resolution_selection.csv"            raw coefficients + MDE
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 0. BUILD THE EXPOSURE MEASURES
*
*    Standardization is over the ONSET sample, not over the full panel as in
*    13b. The two files ask different questions and need different scales:
*    13b asks how the output response varies across the panel's exposure
*    distribution, so it standardizes over the panel; here the comparison is
*    strictly between episodes, so the relevant unit is one standard
*    deviation ACROSS EPISODES. The coefficient reads as the change in the
*    probability of non-default per 1 SD of cross-episode exposure.
* ══════════════════════════════════════════════════════════════════════════
local nexvars claimsgov_assets netclaimsgov_assets

foreach v of local nexvars {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel — skipped (run 13_ifs_nexus.do)"
        continue
    }
    capture drop nx1_`v' nx2_`v' nx3_`v' z_`v' z2_`v' z3_`v'

    * predetermined: exposure one, two and three years before onset
    quietly gen double nx1_`v' = L.`v'
    quietly gen double nx2_`v' = L2.`v'
    quietly gen double nx3_`v' = L3.`v'

    * z_  is the t-1 measure (the one used throughout); z2_/z3_ are the
    * deeper lags used by the run-up check in Section 2(c).
    foreach k in 1 2 3 {
        local zn = cond(`k'==1, "z_`v'", "z`k'_`v'")
        quietly summarize nx`k'_`v' if onset_all==1 & sample==1
        if r(N) > 1 & r(sd) > 0 {
            quietly gen double `zn' = (nx`k'_`v' - r(mean)) / r(sd)
        }
        else {
            quietly gen double `zn' = .
            di as error "  ** `v' at lag `k': too few onsets to standardize"
        }
    }
}
capture label var z_claimsgov_assets    "Std. pre-crisis bank sovereign exposure (t-1)"
capture label var z_netclaimsgov_assets "Std. pre-crisis NET bank sovereign exposure (t-1)"

* ── Coverage: the number that determines what this file can and cannot say ──
di as result _n "════════════════════════════════════════════════════════════"
di as result "EPISODE COUNTS — the binding constraint on everything below"
di as result "════════════════════════════════════════════════════════════"
quietly count if onset_all==1 & sample==1
di as result "  Onsets in the estimation sample:                " %4.0f r(N)
quietly count if onset_all==1 & sample==1 & nondefault==1
local n_nd0 = r(N)
quietly count if onset_all==1 & sample==1 & nondefault==0
local n_def0 = r(N)
di as result "    of which non-default:                         " %4.0f `n_nd0'
di as result "    of which default-linked:                      " %4.0f `n_def0'
quietly count if onset_all==1 & sample==1 & !missing(z_claimsgov_assets)
di as result "  Onsets with nexus coverage at t-1:              " %4.0f r(N)
quietly count if onset_all==1 & sample==1 & !missing(z_claimsgov_assets) & nondefault==0
local n_ev = r(N)
di as result "    of which default-linked (the EVENT count):    " %4.0f `n_ev'
di as result "  => at 10 events per parameter this design supports about " ///
             %3.1f (`n_ev'/10) " regressors."

* ══════════════════════════════════════════════════════════════════════════
* 1. THE RESOLUTION PROBIT LADDER
*
*    Outcome is `nondefault` (1 = the sovereign kept servicing), so a POSITIVE
*    coefficient is the commitment-device reading and a NEGATIVE one is the
*    fragility reading.
*
*    The ladder adds controls in the order the DATA say they matter, not in
*    the order convention would suggest. Table 1 of the Data section found
*    exactly two pre-crisis characteristics that separate the two resolution
*    types: lagged GDP growth (4.62 vs 1.21, p=0.001) and the mean EMBIG
*    spread (403.7 vs 592.0, p=0.001). Those are the two confounds with
*    demonstrated purchase in this sample and they enter first. Public debt
*    enters third despite failing that test, because a resolution regression
*    that does not condition on the debt burden would not be believed
*    whatever the p-value on the balance table said.
*
*    The spread control deserves its own note. It is the severity of the
*    crisis as the market priced it going in, so conditioning on it turns the
*    question from "do high-nexus countries have milder crises that happen to
*    end without default" into "conditional on how badly the market had
*    already priced the sovereign, does exposure still predict resolution".
*    The second is the question worth asking, and it is the harder test.
*
*    NO country fixed effects. With 45-odd episodes spread over 52 countries,
*    most countries contribute a single episode and a country effect would
*    perfectly predict its own outcome — the incidental-parameters problem in
*    its most acute form. SEs are clustered on country instead, which handles
*    the repeat-episode countries without discarding the single-episode ones.
* ══════════════════════════════════════════════════════════════════════════
local L1 ""
local L2 "l1_gdpg"
local L3 "l1_gdpg l_spr_mean"
local L4 "l1_gdpg l_spr_mean l_debt"

local nm1 "Exposure only"
local nm2 "+ pre-crisis growth"
local nm3 "+ spread level"
local nm4 "+ public debt"

tempname R
tempfile resf
postfile `R' str12 model str24 xvar double b double se double z double p ///
    double mfx double mfx_se double mde long N long n_events double auroc ///
    using "`resf'", replace

di as result _n "════════════════════════════════════════════════════════════"
di as result "1. RESOLUTION PROBIT:  Pr(non-default | onset, exposure, X)"
di as result "   Positive coefficient = commitment device; negative = fragility"
di as result "════════════════════════════════════════════════════════════"

forvalues c = 1/4 {
    local ctrls "`L`c''"

    capture noisily probit nondefault z_claimsgov_assets `ctrls' ///
        if onset_all==1 & sample==1, vce(cluster cid)

    if _rc {
        di as error "  ** column `c' (`nm`c'') failed to converge (rc=" _rc ") — skipped"
        di as error "     With this few events that is a result about the sample, not a bug."
        continue
    }

    local N   = e(N)
    quietly count if e(sample) & nondefault==0
    local nev = r(N)

    * coefficient on the exposure term
    local b  = _b[z_claimsgov_assets]
    local se = _se[z_claimsgov_assets]
    local zz = `b' / `se'
    local pp = 2*(1 - normal(abs(`zz')))

    * AUROC: how well the model separates the two resolution types at all.
    * Reported because a probit can carry a significant coefficient and still
    * classify no better than a coin, and the reader should see both.
    local auc = .
    capture quietly lroc, nograph
    if _rc == 0 local auc = r(area)

    * Average marginal effect — the quantity that is actually interpretable.
    * Deliberately NOT `margins, post`: posting would overwrite the probit
    * estimates that the AUROC and the table below still need.
    local mfx = .
    local mse = .
    capture quietly margins, dydx(z_claimsgov_assets)
    if _rc == 0 {
        matrix M = r(table)
        local mfx = M[1,1]
        local mse = M[2,1]
    }

    * MINIMUM DETECTABLE EFFECT at 5%, in probability points per SD of
    * exposure. This is 1.96 x the standard error actually obtained: any true
    * effect smaller than this could not have been distinguished from zero by
    * this design no matter what the world looks like. Report it beside every
    * null so the null can be read as evidence rather than as silence.
    local mde = .
    if !missing(`mse') local mde = 1.96 * `mse' * 100

    di as result _n "  Column `c' — `nm`c''"
    di as result "    b(exposure) = " %7.3f `b' "   SE = " %6.3f `se' ///
                 "   z = " %6.2f `zz' "   p = " %5.3f `pp'
    if !missing(`mfx') ///
    di as result "    Avg marginal effect = " %6.1f (100*`mfx') " pp of Pr(non-default) per 1 SD"
    if !missing(`mde') ///
    di as result "    Min detectable eff. = " %6.1f `mde' " pp at 5 pct — smaller true effects were undetectable here"
    if !missing(`auc') ///
    di as result "    AUROC = " %5.3f `auc'
    di as result "    N = " %3.0f `N' "   default events = " %3.0f `nev' ///
                 "   parameters = " %2.0f (wordcount("`ctrls'") + 2)

    post `R' ("col`c'") ("z_claimsgov_assets") (`b') (`se') (`zz') (`pp') ///
        (`mfx') (`mse') (`mde') (`N') (`nev') (`auc')
}

* Re-estimate for the table (margins ... post above destroys the estimates)
forvalues c = 1/4 {
    capture quietly probit nondefault z_claimsgov_assets `L`c'' ///
        if onset_all==1 & sample==1, vce(cluster cid)
    if _rc == 0 {
        quietly estimates store E`c'
        local elist `elist' E`c'
    }
}

capture esttab `elist' using "$tabs/table13f_resolution_selection.rtf", ///
    replace b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    scalars("N Episodes" "r2_p Pseudo R2") ///
    mtitles("`nm1'" "`nm2'" "`nm3'" "`nm4'") ///
    title("Table 13f. Does pre-crisis bank exposure to the sovereign predict resolution?") ///
    addnotes("Probit. Sample = spread-crisis onsets only. Outcome = 1 if the episode was resolved WITHOUT default," ///
             "so a positive coefficient means higher pre-crisis exposure predicts NON-default (the commitment-device" ///
             "reading of Gennaioli, Martin & Rossi 2014); a negative coefficient is the fragility reading. Exposure is" ///
             "bank claims on general government as a share of total bank assets at t-1, standardized across episodes." ///
             "SEs clustered by country. No country fixed effects: most countries contribute one episode. Read every" ///
             "column against the default-event count, not N — see the minimum detectable effects in the log.")
if _rc di as error "  ** esttab failed (rc=" _rc ") — CSV still written"

* ══════════════════════════════════════════════════════════════════════════
* 2. FUNCTIONAL-FORM AND MEASUREMENT ROBUSTNESS
*
*    Three checks, each with a stated failure condition.
*
*    (a) LOGIT and a LINEAR PROBABILITY MODEL on the preferred column. Probit
*        and logit will agree; the LPM is the check that matters, because with
*        this few events a maximum-likelihood binary model can be driven by
*        its functional form in the tails. If the LPM disagrees in sign, the
*        probit result is an artifact of the link function.
*
*    (b) NET exposure (claims on government minus government deposits, over
*        assets). The gross measure counts a bank that holds government paper
*        and government deposits in equal measure as exposed; the net measure
*        does not. If the two disagree, what is being measured is the size of
*        the government's banking relationship rather than its exposure.
*
*    (c) EXPOSURE AT t-2 AND t-3. If banks were accumulating sovereign paper
*        in anticipation of the distress that determines the resolution, the
*        t-1 measure is contaminated by exactly the thing it is meant to
*        predict. A coefficient that holds at deeper lags is measuring a
*        structural feature of the banking system; one that appears only at
*        t-1 is measuring the run-up.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "2. ROBUSTNESS"
di as result "════════════════════════════════════════════════════════════"

di as result _n "  (a) Link function — probit vs logit vs LPM, on column 3"
foreach est in probit logit {
    capture quietly `est' nondefault z_claimsgov_assets `L3' ///
        if onset_all==1 & sample==1, vce(cluster cid)
    if _rc == 0 {
        di as result "    `est': b = " %7.3f _b[z_claimsgov_assets] ///
                     "   p = " %5.3f (2*(1-normal(abs(_b[z_claimsgov_assets]/_se[z_claimsgov_assets]))))
        post `R' ("`est'") ("z_claimsgov_assets") (_b[z_claimsgov_assets]) ///
            (_se[z_claimsgov_assets]) (_b[z_claimsgov_assets]/_se[z_claimsgov_assets]) ///
            (2*(1-normal(abs(_b[z_claimsgov_assets]/_se[z_claimsgov_assets])))) ///
            (.) (.) (.) (e(N)) (.) (.)
    }
    else di as error "    `est' failed (rc=" _rc ")"
}
capture quietly regress nondefault z_claimsgov_assets `L3' ///
    if onset_all==1 & sample==1, vce(cluster cid)
if _rc == 0 {
    di as result "    LPM:    b = " %7.3f _b[z_claimsgov_assets] ///
                 "   p = " %5.3f (2*ttail(e(df_r), abs(_b[z_claimsgov_assets]/_se[z_claimsgov_assets]))) ///
                 "   (directly in probability points per SD)"
    post `R' ("lpm") ("z_claimsgov_assets") (_b[z_claimsgov_assets]) (_se[z_claimsgov_assets]) ///
        (_b[z_claimsgov_assets]/_se[z_claimsgov_assets]) ///
        (2*ttail(e(df_r), abs(_b[z_claimsgov_assets]/_se[z_claimsgov_assets]))) ///
        (_b[z_claimsgov_assets]) (_se[z_claimsgov_assets]) (1.96*_se[z_claimsgov_assets]*100) ///
        (e(N)) (.) (.)
}

di as result _n "  (b) Gross vs net sovereign exposure, on column 3"
foreach v in claimsgov_assets netclaimsgov_assets {
    capture confirm variable z_`v', exact
    if _rc continue
    capture quietly probit nondefault z_`v' `L3' if onset_all==1 & sample==1, vce(cluster cid)
    if _rc == 0 {
        di as result "    `v': b = " %7.3f _b[z_`v'] ///
                     "   p = " %5.3f (2*(1-normal(abs(_b[z_`v']/_se[z_`v'])))) ///
                     "   N = " %3.0f e(N)
        post `R' ("gross_net") ("z_`v'") (_b[z_`v']) (_se[z_`v']) ///
            (_b[z_`v']/_se[z_`v']) (2*(1-normal(abs(_b[z_`v']/_se[z_`v'])))) ///
            (.) (.) (.) (e(N)) (.) (.)
    }
    else di as error "    `v' failed (rc=" _rc ")"
}

di as result _n "  (c) Depth of the lag — is this structure or the run-up?"
foreach k in 1 2 3 {
    local zn = cond(`k'==1, "z_claimsgov_assets", "z`k'_claimsgov_assets")
    capture confirm variable `zn', exact
    if _rc continue
    capture quietly probit nondefault `zn' `L3' ///
        if onset_all==1 & sample==1, vce(cluster cid)
    if _rc == 0 {
        di as result "    exposure at t-`k': b = " %7.3f _b[`zn'] ///
                     "   p = " %5.3f (2*(1-normal(abs(_b[`zn']/_se[`zn'])))) ///
                     "   N = " %3.0f e(N)
        post `R' ("lag`k'") ("`zn'") (_b[`zn']) (_se[`zn']) ///
            (_b[`zn']/_se[`zn']) (2*(1-normal(abs(_b[`zn']/_se[`zn'])))) ///
            (.) (.) (.) (e(N)) (.) (.)
    }
    else di as error "    exposure at t-`k' failed (rc=" _rc ")"
}

* ══════════════════════════════════════════════════════════════════════════
* 3. THE DEVELOPMENT CONFOUND
*
*    13b established that GDP-scaled exposure ratios sort on financial
*    development rather than on vulnerability, which is why this project uses
*    the assets-scaled measure. A portfolio share should not carry that
*    confound — but "should not" is an argument, and the check is cheap.
*    Split the onsets at the median of exposure and compare mean pre-crisis
*    real GDP per capita across the two bins. A ratio far from one means the
*    exposure variable is partly an income variable and the coefficient above
*    cannot be read as being about the nexus.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "3. IS THE EXPOSURE MEASURE SORTING ON DEVELOPMENT?"
di as result "════════════════════════════════════════════════════════════"
capture confirm variable gdppc_real, exact
if _rc == 0 {
    quietly summarize z_claimsgov_assets if onset_all==1 & sample==1, detail
    local med = r(p50)
    quietly summarize gdppc_real if onset_all==1 & sample==1 & z_claimsgov_assets <= `med'
    local inc_lo = r(mean)
    quietly summarize gdppc_real if onset_all==1 & sample==1 & z_claimsgov_assets >  `med'
    local inc_hi = r(mean)
    local ratio = `inc_hi' / `inc_lo'
    di as result "    Mean pre-crisis real GDP per capita, LOW-exposure onsets:  " %10.0f `inc_lo'
    di as result "    Mean pre-crisis real GDP per capita, HIGH-exposure onsets: " %10.0f `inc_hi'
    di as result "    Ratio (high/low) = " %5.2f `ratio'
    if `ratio' > 1.5 | `ratio' < 0.667 {
        di as error "    ** FLAG: the exposure bins differ markedly in income. The"
        di as error "       coefficient in Section 1 is partly an income coefficient"
        di as error "       and must not be reported as a clean nexus result."
    }
    else {
        di as result "    Bins are comparable in income — the portfolio-share construction"
        di as result "    does what it was chosen to do, and the coefficient can be read"
        di as result "    as being about exposure rather than about development."
    }
}
else di as error "  ** gdppc_real not in panel — development check skipped"

* ══════════════════════════════════════════════════════════════════════════
* 4. THE SAME QUESTION FROM THE PANEL SIDE
*
*    Section 1 throws away every tranquil country-year, which is why it has
*    45-odd observations. The same contrast can be recovered from the full
*    panel with a multinomial logit over {tranquil, non-default onset,
*    default-linked onset}, tranquil as the base category. The controls are
*    then pinned by roughly five hundred tranquil years instead of by the
*    episodes themselves, and the difference between the two exposure
*    coefficients is testable INSIDE one model rather than across two.
*
*    Note what this is and is not. The [def]-vs-[nd] contrast in this model is
*    THE SAME OBJECT as the probit in Section 1, estimated a second way — not
*    a second finding. Section 1 is the direct estimate and does the
*    identifying work; this is the efficient one, bought at the price of the
*    independence-of-irrelevant-alternatives assumption, which is not
*    innocuous when the two crisis types are close substitutes. Where they
*    disagree, believe Section 1 and say so.
*
*    The two separate vs-tranquil probits that follow mirror 08b_aipw.do's
*    own propensity design, and are reported for comparability with the
*    paper's existing first stages. THEY ARE A DIAGNOSTIC AND NOTHING ELSE:
*    the weighting models in 08b/13d are deliberately NOT modified, because
*    adding a variable to $ctrl_core there would change every IPW and AIPW
*    estimate in the paper. A Clogg z across those two probits would also be
*    wrong, since both samples share the same tranquil control years and the
*    two coefficients are therefore not independent — which is precisely the
*    problem the multinomial logit solves.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "4. PANEL-SIDE VERSION: multinomial logit, tranquil = base"
di as result "════════════════════════════════════════════════════════════"

capture drop resol3
quietly gen byte resol3 = .
quietly replace resol3 = 0 if sample==1 & onset_all==0
quietly replace resol3 = 1 if sample==1 & onset_nd==1
quietly replace resol3 = 2 if sample==1 & onset_def==1
label define resol3L 0 "Tranquil" 1 "Non-default onset" 2 "Default-linked onset", replace
label values resol3 resol3L

capture noisily mlogit resol3 z_claimsgov_assets $ctrl_core ///
    if !missing(resol3), baseoutcome(0) vce(cluster cid)

if _rc {
    di as error "  ** mlogit failed (rc=" _rc ") — with cells this thin that is"
    di as error "     an expected outcome, not a coding error. Section 1 stands alone."
}
else {
    local b_nd  = _b[1:z_claimsgov_assets]
    local b_def = _b[2:z_claimsgov_assets]
    di as result _n "    Exposure coefficient, non-default vs tranquil:   " %7.3f `b_nd'
    di as result "    Exposure coefficient, default-linked vs tranquil: " %7.3f `b_def'
    di as result "    Difference (nd - def):                           " %7.3f (`b_nd' - `b_def')

    * The test that matters: do the two crisis types load differently on
    * exposure? If they do not, exposure predicts CRISIS but not RESOLUTION,
    * and the answer to the question this file was written to ask is no.
    capture quietly test [1]z_claimsgov_assets = [2]z_claimsgov_assets
    if _rc == 0 {
        di as result "    H0: equal loading on exposure —  chi2 = " %6.2f r(chi2) ///
                     "   p = " %5.3f r(p)
        if r(p) < 0.10 {
            di as result "    => the two resolution types DO load differently on pre-crisis"
            di as result "       exposure. Read the sign of the difference above against the"
            di as result "       two hypotheses stated in the header."
        }
        else {
            di as result "    => no detectable difference in loading. On this evidence"
            di as result "       exposure does not predict how a crisis is resolved, only"
            di as result "       (at most) whether one occurs. Report as a null, with the"
            di as result "       MDE from Section 1 attached."
        }
        post `R' ("mlogit") ("nd_minus_def") (`b_nd' - `b_def') (.) (.) (r(p)) ///
            (.) (.) (.) (e(N)) (.) (.)
    }
}

di as result _n "  Mirror of 08b_aipw.do's first stages (DIAGNOSTIC ONLY —"
di as result "  the weighting models in 08b/13d are not modified):"
* Matches 08b_aipw.do's own cz_def (see its header for the swap rationale,
* citing 21_aipw_flow.do's "SECOND PREDICTOR CHANGE"; both terms are
* DEFAULT-LINKED-SPECIFIC) -- kept in sync since this block exists only to
* mirror 08b's first stage for comparability.
local predictors_z2 l_fedfunds l_contagion_dist_def years_since_def_onset
foreach s in nd def {
    if "`s'" == "nd"  local rival onset_def
    else              local rival onset_nd
    capture quietly probit onset_`s' z_claimsgov_assets $ctrl_core `predictors_z2' ///
        if sample==1 & `rival'==0, vce(cluster cid)
    if _rc == 0 {
        di as result "    onset_`s' vs tranquil: b(exposure) = " %7.3f _b[z_claimsgov_assets] ///
                     "   p = " %5.3f (2*(1-normal(abs(_b[z_claimsgov_assets]/_se[z_claimsgov_assets])))) ///
                     "   N = " %5.0f e(N)
    }
    else di as error "    onset_`s' first stage failed (rc=" _rc ")"
}

* ══════════════════════════════════════════════════════════════════════════
* 5. SAVE
* ══════════════════════════════════════════════════════════════════════════
postclose `R'
preserve
    use "`resf'", clear
    label var b        "Coefficient on standardized pre-crisis exposure"
    label var mfx      "Average marginal effect (probability units per SD)"
    label var mde      "Minimum detectable effect at 5% (pp per SD)"
    label var n_events "Default-linked episodes in the estimation sample"
    label var auroc    "Area under the ROC curve"
    export delimited "$tabs/resolution_selection.csv", replace
    di as result _n "Saved: $tabs/resolution_selection.csv"
restore

di as result _n "════════════════════════════════════════════════════════════"
di as result "HOW TO READ THIS FILE"
di as result "════════════════════════════════════════════════════════════"
di as result "  Report whatever comes out, including a null, and report it with"
di as result "  the default-event count and the minimum detectable effect next to"
di as result "  it. Section 2(c) decides whether any surviving coefficient is a"
di as result "  structural feature of the banking system or an artifact of the"
di as result "  pre-crisis run-up, and Section 3 decides whether it is about the"
di as result "  nexus at all or about income. If the coefficient clears"
di as result "  significance and survives both, the write-up gains a selection"
di as result "  result — and Section 8 acquires the qualification set out in the"
di as result "  header, because the amplification estimated there is then"
di as result "  conditional on a selected set of default-linked episodes."

di as result _n "13f_resolution_selection.do complete."
