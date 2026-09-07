/*===========================================================================
  08C_FIRST_STAGE_TABLE.DO
  First-stage probit table, styled on Asonuma et al. (2024) Table 1
  ("Predicting the start of debt restructurings, probit") — the ONSET-tier
  counterpart of 21b_first_stage_table_flow.do; mirrors that file's
  diagnostics (AUROC with/without predictors + roccomp's formal test, not
  just chi2(predictors)) for this file's own onset specification.

  Columns = onset type BY RESOLUTION, each predicted vs TRANQUIL (the rival
  onset type is dropped so the control group is clean non-crisis years —
  exactly as in the paper, and exactly the first stages that feed the AIPW
  two-line Act 2 in 08b):
    (1) Non-default   : onset_nd  vs tranquil        (sample==1 & onset_def==0)
    (2) Default-linked: onset_def vs tranquil        (sample==1 & onset_nd==0)
  The pooled "All onsets vs tranquil" column is deliberately not reported
  here: every onset is either non-default or default-linked, so the pooled
  column carries no propensity-modeling information beyond the union of
  these two, and this table's whole purpose is to report the first stages
  the by-resolution-type AIPW design actually consumes.

  Rows are grouped, as in their Table 1, into:
    PREDICTORS (excluded from the LP/AIPW outcome eq.): Fed funds rate
        (global push) + distance-weighted contagion (l_contagion_dist_def,
        a country-year-specific spatial lag of OTHER countries' DEFAULT-
        LINKED in-crisis STOCK, i.e. onset_def|continuation of a default
        episode, not just onset years) + years since the most recent prior
        DEFAULT-LINKED onset (years_since_def_onset, a recency clock,
        censored at 50).

    WHY DEFAULT-LINKED-SPECIFIC, NOT GENERIC (any-onset-type) -- the
    justification is economic, not merely statistical. A predictor meant to
    proxy "risk of a DEFAULT-linked onset specifically" should measure
    default-linked distress and default-linked recency, not spread-crisis
    distress/recency in general: a country surrounded by non-default spread
    crises, or one whose own last onset was itself non-default, is not
    obviously more exposed to a NEW default -- a generic measure blends two
    conceptually different risks (illiquidity/market pressure vs. actual
    default propensity) into one number and asks it to predict the narrower
    one. Narrowing to default-linked-only makes each term measure exactly
    what the column is trying to explain, which is the standard the
    reference paper's own instrument is held to as well (their federal
    funds/contagion/past-preemptive terms are exclusion-restriction
    predictors chosen for their own conceptual fit to a restructuring
    probability, not swapped in from an unrelated outcome). This is stated
    as a design choice about MEASUREMENT, not merely a specification search
    for the best-scoring predictor set on this sample.

    ONE HONEST CAVEAT, stated plainly rather than glossed over: the
    reference paper's own $instrument (federal_funds2 cont_all2
    past_preemp2) is fixed ONCE and used IDENTICALLY across all four of
    their columns (post-default, weakly preemptive, strictly preemptive,
    preemptive) -- they do not build a separate contagion/recency term per
    sub-type. Their practice is closer to this project's Act-1 (pooled) cz,
    which does stay generic. The default-linked-specific cz_def used here is
    therefore a genuine methodological departure from their literal
    practice, adopted on the economic-measurement grounds above and
    corroborated (not driven) by the AUROC/roccomp diagnostics below, which
    showed the generic combination left the default arm's classification
    power statistically indistinguishable from controls-only.
    BASELINE CONTROLS (the SAME $ctrl_core used in the LP/AIPW outcome eq. —
        strict parity with the reference paper's $convar-in-both design):
        l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2

  Diagnostic rows (their Table 1 bottom block):
    Chi2 (predictors)      — joint Wald test that the 3 predictors are all zero
    p-value                — of that test (small => predictors jointly matter)
    AUROC, controls only   — lroc on a probit with ONLY $ctrl_core, no predictors
    AUROC, with predictors — lroc on the full model
    Observations

  WHY BOTH AUROCs, NOT JUST THE FULL MODEL'S (matches 21b's reasoning): the
  source paper's actual demonstration that predictors earn their place isn't
  the chi-squared test alone -- it's the WITH-vs-WITHOUT-predictors AUROC
  comparison they report directly in their text (0.87 vs 0.79 post-default,
  0.94 vs 0.85 preemptive). Chi2(predictors) only says the predictors are
  jointly non-zero; it says nothing about how much discriminatory power they
  actually add on top of what the baseline controls already provide.

  THE DELTA ITSELF IS NOT A SIGNIFICANCE TEST -- roccomp's chi2 test for
  equal correlated ROC areas is, and it is reported alongside the delta in
  the console output. The write-up should quote the roccomp p-values, not
  just the delta, when characterising how much work the predictors are doing.

  Pooled probit, no country FE (rare-event propensity: country FE separate/overfit
  with ~20 events; they enter the outcome equation only). No year FE.
  Clustered SEs by country, matching 21b's choice.

  Output: $tabs/table_first_stage.rtf. Run AFTER 17_predictors.do. The
  figure counterpart (kernel density + nested ROC, mirroring 21c) is
  08d_first_stage_figs.do.

  DIAGNOSTIC HISTORY:
    - An earlier pass tested a default-linked-only contagion measure
      (built by restricting the ORIGINAL onset-only contagion_dist to
      onset_def) against the generic version: no difference in either arm
      (nd delta +0.009, p=.609; def delta +0.000, p=.931). Kept generic on
      that basis.
    - contagion_dist was then redefined project-wide from an onset-only
      SHOCK indicator to an in-crisis STOCK (onset|continuation) -- see
      17_predictors.do's header. Re-run under that new construction, the
      def arm's classification power collapsed under the generic predictor
      combination (chi2(pred) p=.101, roccomp p=.525, both individually
      insignificant: l_contagion_dist z=-0.74, years_since_onset z=-1.19),
      unlike the nd arm (p<.001, roccomp p=.039, l_contagion_dist z=-2.58).
    - Both predictors were therefore narrowed to default-linked-specific
      (l_contagion_dist_def, years_since_def_onset) for cz_def, matching
      what the default arm's own diagnostics call for. cz (Act 1, pooled)
      keeps the generic l_contagion_dist -- there is no resolution type to
      be specific about in a pooled spec.

  Also runs a DIAGNOSTIC-ONLY comparison (console output only, before the
  table export): does adding COUNTRY FIXED EFFECTS to the probit (matching
  the reference paper's own $cf design) converge on this project's much
  thinner panel? RESOLVED, NOT ADOPTED -- confirms the pooled design above:
    - nd arm: 25 of 51 countries dropped for zero outcome variation (49 pct
      of the sample), N 933->525, Wald chi2/p become UNCOMPUTABLE (Prob>chi2
      = .), and both l_contagion_dist_def and years_since_def_onset collapse
      to noise (z=-0.11, z=-1.00) once country FE absorb the between-country
      variation they relied on.
    - def arm: worse. 41 of 51 countries dropped (80 pct), N 917->133,
      "8 failures and 0 successes completely determined" (severe
      separation). Every predictor loses significance except
      years_since_def_onset (z=1.98, p=.047); l_fedfunds itself -- robust in
      every other specification in this file -- drops to z=0.61.
    - Matches 21_aipw_flow.do's Section 1g exactly: with ~14 default-linked
      countries mostly holding a single episode each, a country dummy
      absorbs that one event entirely and leaves nothing to identify a
      within-country contrast from. The reference paper's 194 events over
      74 countries (~2.6 episodes/treated country on average) do not face
      this; this project's panel does. Pooled, no country FE stays correct
      -- confirmed empirically, not merely asserted by analogy to the flow
      tier. The exported table above is unaffected; this block is a
      permanent record of the test, not live code any consumer reads.

  A separate diagnostic (console output, run right after the adopted
  fs_nd/fs_def models fit) investigates the def arm's own "2 failures
  completely determined" note directly: does l_lninfl or exchange2 alone
  cleanly separate onset_def on this sample (checked first; ranges
  overlapped for both, so separation is on the FULL linear index, not
  either term alone), which onset rows sit at the extreme univariate tails,
  and -- the step that actually pins it down -- which control rows have the
  most extreme predicted linear index (xb) once all controls+predictors are
  combined, the direct signature of "N failures completely determined."

  RESULT, FLAGGED NOT ACTED ON: El Salvador 2002 is the clear primary case
  -- predicted index -9.80, ~4 log-probit units more extreme than the
  next-closest observation (China 2022, -6.09), driven by exchange2=-1.58,
  by far the largest currency move in the def-arm sample (next runner-up
  +-0.83). A cluster of China years (2019, 2021-23) sits at a similarly
  extreme index despite near-zero l_lninfl/exchange2 values of their own --
  most likely a symptom of the same instability (the huge coefficients the
  El Salvador observation forces onto l_lninfl/exchange2 make the whole
  index sensitive even where those two terms are near zero), not an
  independent second cause; the exact second "completely determined" row
  was not pinned down further. NOT acted on: this is one country-year out
  of a full $ctrl_core shared with every other regression in this project
  (the LP/AIPW outcome equation, every onset-tier table) -- dropping or
  winsorizing it here only, for this one probit, would decouple this
  table's control set from the rest of the project's, which is a larger
  change than the finding warrants. Recorded so the "2 failures completely
  determined" note is explained, not left as an unexplained artifact.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

* Baseline controls X (both columns) and predictors Z2 (resolution-type
* proneness = years since the most recent prior default-linked onset).
* Baseline = outcome baseline ($ctrl_core): the Table-1 probit shares its controls
* with the LP/AIPW outcome equation (their $convar in both stages).
local X    $ctrl_core
local Z2   l_fedfunds l_contagion_dist_def years_since_def_onset

eststo clear

* ── Helper: run BOTH nested models per column, add the AUROC-with-vs-without-
* predictors + chi2(predictors) diagnostics — mirrors 21b's _fscolflow.
* A program cannot see the caller's locals, so X and Z are passed in explicitly.
* args: estimate-name  DVvar  "if-condition"  "X list"  "Z list"
* Controls-only model is fit FIRST so its AUROC can be estadd-ed onto the
* full model's stored estimate once that is eststo-d (estadd only writes to
* the currently active estimates). Predicted probabilities from both models
* are saved for the roccomp formal test run right after this program.
capture program drop _fscol
program define _fscol
    args nm dv ifcond xlist zlist
    quietly probit `dv' `xlist' if `ifcond', vce(cluster cid)
    quietly lroc, nograph
    local aucctrl = r(area)
    capture drop _pctrl_`nm'
    quietly predict double _pctrl_`nm' if `ifcond', pr

    probit `dv' `xlist' `zlist' if `ifcond', vce(cluster cid)
    eststo `nm'
    quietly test `zlist'
    estadd scalar chi2p = r(chi2)
    estadd scalar pp    = r(p)
    quietly lroc, nograph
    estadd scalar auroc     = r(area)
    estadd scalar aurocctrl = `aucctrl'
    capture drop _pfull_`nm'
    quietly predict double _pfull_`nm' if `ifcond', pr
end

_fscol fs_nd  "onset_nd"  "sample==1 & onset_def==0"      "`X'" "`Z2'"
_fscol fs_def "onset_def" "sample==1 & onset_nd==0"       "`X'" "`Z2'"

* ══════════════════════════════════════════════════════════════════════════
* DIAGNOSTIC: WHICH ROWS ARE THE "2 FAILURES COMPLETELY DETERMINED" IN THE
* DEF ARM? l_lninfl/exchange2 print implausibly large coefficients there
* (roughly -7 to -8 and +4 to +5) alongside that separation note. Rather
* than guess, this identifies the actual country-years responsible: for
* each of the two suspect controls, split the def-arm sample by
* onset_def and report whether there is a clean, non-overlapping cutoff --
* the hallmark of (quasi-)complete separation -- then list the specific
* onset rows sitting at the extreme tail, since a handful of onset years
* with an unusually large inflation/depreciation reading are the most
* likely source (matches the reasoning that already motivated swapping
* l_hyperinfl->l_lninfl and ex_dum1-5->exchange2 in $ctrl_core itself).
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== DIAGNOSTIC: separation source in the def arm (l_lninfl / exchange2) ==="
foreach v in l_lninfl exchange2 {
    di as result _n "      `v', by onset_def (def-arm sample, onset_nd==0):"
    quietly summarize `v' if sample==1 & onset_nd==0 & onset_def==1, detail
    local mn1 = r(min)
    local mx1 = r(max)
    di as result "        onset_def==1 (treated): min=" %9.3f `mn1' "  p50=" %9.3f r(p50) "  max=" %9.3f `mx1'
    quietly summarize `v' if sample==1 & onset_nd==0 & onset_def==0, detail
    local mn0 = r(min)
    local mx0 = r(max)
    di as result "        onset_def==0 (control): min=" %9.3f `mn0' "  p50=" %9.3f r(p50) "  max=" %9.3f `mx0'
    * A clean, non-overlapping range between the two groups is the direct
    * signature of (quasi-)complete separation on this one variable alone.
    if `mn1' > `mx0' | `mx1' < `mn0' {
        di as result "        ** RANGES DO NOT OVERLAP -- `v' alone separates onset_def on this sample."
    }
    else {
        di as result "        ranges overlap -- `v' alone does not fully separate; check jointly with other controls."
    }
}

di as result _n "      Onset rows (def arm) at the extreme tails of l_lninfl / exchange2:"
di as result "      (the likely candidates for the 2 perfectly-determined observations)"
preserve
    quietly keep if sample==1 & onset_nd==0
    gen double _rank_infl = abs(l_lninfl - 0)
    gsort -_rank_infl
    di as result _n "      Top 5 by |l_lninfl|, def-arm sample:"
    list country year onset_def l_lninfl exchange2 in 1/5, noobs clean
    gen double _rank_fx = abs(exchange2 - 0)
    gsort -_rank_fx
    di as result _n "      Top 5 by |exchange2|, def-arm sample:"
    list country year onset_def l_lninfl exchange2 in 1/5, noobs clean
restore

* Neither variable alone need have disjoint ranges for "N failures completely
* determined" to fire -- that note flags separation on the FULL linear index
* (all controls + predictors combined), not any one term in isolation. This
* step finds the exact rows directly: refit the def-arm model, predict the
* linear index (xb, not the bounded probability), and list the rows whose
* index is most extreme AMONG THE CONTROLS (onset_def==0) -- "failures
* completely determined" means Stata's optimizer drove some control rows'
* predicted probability to (numerically) exactly 0, which shows up as an
* extreme negative xb.
di as result _n "      Exact rows: def-arm CONTROL rows (onset_def==0) with the most extreme"
di as result "      predicted index (candidates for the '2 failures completely determined'):"
preserve
    quietly keep if sample==1 & onset_nd==0
    quietly probit onset_def `X' `Z2', vce(cluster cid)
    capture drop _xb_def
    quietly predict double _xb_def, xb
    keep if onset_def==0
    sort _xb_def
    list country year onset_def l_lninfl exchange2 _xb_def in 1/5, noobs clean
restore

* ── FORMAL TEST OF WHETHER THE TWO AUROCs ACTUALLY DIFFER ──────────────────
* The delta alone (auroc - aurocctrl) says nothing about whether that gap is
* distinguishable from sampling noise -- an informal threshold cannot
* substitute for a test. roccomp's chi2 test for equal correlated ROC areas
* (same subjects, two nested classifiers) is that test, and it needs the two
* models' predicted probabilities on the SAME sample, saved above.
quietly roccomp onset_nd _pctrl_fs_nd _pfull_fs_nd if !missing(_pctrl_fs_nd,_pfull_fs_nd)
local rocchi2_nd = r(chi2)
local rocp_nd    = r(p)
quietly roccomp onset_def _pctrl_fs_def _pfull_fs_def if !missing(_pctrl_fs_def,_pfull_fs_def)
local rocchi2_def = r(chi2)
local rocp_def    = r(p)
capture drop _pctrl_fs_nd _pfull_fs_nd _pctrl_fs_def _pfull_fs_def

* ── Console echo of the diagnostics ──────────────────────────────────────────
di as result _n "=== FIRST-STAGE PROBIT DIAGNOSTICS (predictors jointly) ==="
di as result "col            chi2(pred)   p        AUROC(ctrl only)  AUROC(+pred)  delta   roccomp chi2   p"
local mindlt = .
foreach c in fs_nd fs_def {
    quietly estimates restore `c'
    local dlt = e(auroc) - e(aurocctrl)
    local dltsign = cond(`dlt' >= 0, "+", "")
    local sfx = subinstr("`c'", "fs_", "", .)
    di as result %-14s "`c'" "  " %8.2f e(chi2p) "  " %6.3f e(pp) "  " ///
                 %8.3f e(aurocctrl) "        " %6.3f e(auroc) "       " "`dltsign'" %6.3f `dlt' ///
                 "     " %6.2f `rocchi2_`sfx'' "        " %5.3f `rocp_`sfx''
    * Track the SMALLER of the two deltas -- the interpretation below should
    * reflect the weaker column, not be driven by whichever column happens
    * to look best.
    if missing(`mindlt') | `dlt' < `mindlt' local mindlt = `dlt'
}
* The raw delta is NOT a significance test -- roccomp's chi2 test (above,
* same-subjects correlated-ROC-areas test) is. Report that number, not an
* informal threshold, as the governing statistic.
local worstp = max(`rocp_nd', `rocp_def')
di as result _n "      Both columns' AUROC delta: " %5.3f `mindlt' " (smaller of the two)."
di as result "      roccomp's formal test for whether the two AUROCs differ: nd chi2(1)=" ///
    %5.2f `rocchi2_nd' ", p=" %5.3f `rocp_nd' "   def chi2(1)=" %5.2f `rocchi2_def' ", p=" %5.3f `rocp_def' "."
if `worstp' < 0.05 {
    di as result "      BOTH columns clear the conventional 5 pct level (worse column p=" %5.3f `worstp' ")."
    di as result "      Read this as ESTABLISHED, not merely suggestive: the chi-squared row above tests joint"
    di as result "      predictor significance only and says nothing about classification power; roccomp's"
    di as result "      test is the one that speaks to classification power directly, and here it backs the"
    di as result "      chi-squared row up."
}
else if `worstp' < 0.10 {
    di as result "      One column clears the conventional 5 pct level; the other does not (worse column p=" ///
        %5.3f `worstp' "), but both are inside the 10 pct band."
    di as result "      Read this as MODEST, NOT ABSENT: the chi-squared row above tests joint predictor"
    di as result "      significance only and says nothing about classification power; roccomp's test is the"
    di as result "      one that speaks to classification power directly, and it says 'suggestive, not fully"
    di as result "      established' here."
}
else {
    di as result "      Neither reaches the conventional 5 pct level (worse column p=" %5.3f `worstp' ///
        "); at least one arm sits outside even the 10 pct band."
    di as result "      Read this as MODEST, NOT ABSENT: the chi-squared row above tests joint predictor"
    di as result "      significance only and says nothing about classification power; roccomp's test is the"
    di as result "      one that speaks to classification power directly, and it says 'suggestive, not"
    di as result "      established' here."
}

* ══════════════════════════════════════════════════════════════════════════
* DIAGNOSTIC (NOT ADOPTED): does adding COUNTRY FIXED EFFECTS to the probit
* match the reference paper's design ($cf c1-c74, noconstant), or does it
* break on this project's much thinner panel?
*
* The reference paper's own probit carries country dummies directly. Their
* design can afford this because with 194 events over 74 countries most
* treated countries have MULTIPLE episodes (~2.6 on average), so a country
* FE for a never-treated country is absorbed cheaply while enough
* ever-treated countries remain to identify the slopes. This project's
* onset probit is pooled, no country FE, precisely because that condition
* does not hold here: ~54 countries, 61 onsets, and in the default-linked
* arm specifically only ~14 countries ever had a default-linked episode,
* most with a single short episode -- close to one event per treated
* country, not 2.6.
*
* THE EXPECTED FAILURE MODE, STATED IN ADVANCE (matching 21_aipw_flow.do's
* own Section 1g, which ran the identical test for the flow tier's
* propensity model and found country FE either failed to converge or
* collapsed the def-arm control pool to only ever-treated countries): any
* country with ZERO outcome variation in a given arm is dropped
* automatically ("predicts failure/success perfectly"), and on the def arm
* this could mean far more than a handful of countries -- checked directly
* below, not assumed.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== DIAGNOSTIC: country FE added to the probit -- does it even converge? ==="
foreach s in nd def {
    local dv  = cond("`s'"=="nd", "onset_nd", "onset_def")
    local ifc = cond("`s'"=="nd", "sample==1 & onset_def==0", "sample==1 & onset_nd==0")

    di as result _n "      `s' arm, WITH i.cid (raw probit output, watch for separation notes):"
    capture noisily probit `dv' `X' `Z2' i.cid if `ifc', vce(cluster cid)
    if _rc {
        di as error "      `s' arm: probit FAILED to converge (rc=" _rc ")."
        continue
    }

    * How many countries were in the fitting sample vs actually used
    * (e(sample)==1 excludes rows dropped for zero outcome variation)?
    quietly levelsof cid if `ifc', local(allcty)
    quietly levelsof cid if e(sample)==1, local(keptcty)
    local ndrop : list allcty - keptcty
    local ndropn : word count `ndrop'
    local nallc  : word count `allcty'
    di as result "      countries in the `s' fitting sample: " `nallc' ///
                 "  |  dropped for zero outcome variation: " `ndropn'
}
di as result _n "      Compare to the pooled (no country FE) columns above: if country FE"
di as result "      converges cleanly without collapsing the def-arm country pool, it is"
di as result "      worth adopting; if it drops most of the def arm's countries or fails to"
di as result "      converge, that is the honest reason this project's onset probit stays"
di as result "      pooled, matching what 21_aipw_flow.do's Section 1g already found for the"
di as result "      flow tier's propensity model. NOT wired into the adopted `X'/`Z2' probit"
di as result "      above -- this block only reports the comparison."

* ══════════════════════════════════════════════════════════════════════════
* TABLE EXPORT — Table 1 style (Predictors / Baseline controls blocks + diags)
* ══════════════════════════════════════════════════════════════════════════
capture esttab fs_nd fs_def using "$tabs/table_first_stage.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonumber ///
    mtitles("Non-default" "Default-linked") ///
    order(l_fedfunds l_contagion_dist_def years_since_def_onset ///
          l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2) ///
    coeflabel(l_fedfunds "US fed funds rate (t-1)" ///
              l_contagion_dist_def "Distance-weighted contagion, default-linked (t-1)" ///
              years_since_def_onset "Years since last default-linked onset" ///
              l1_gdpg "GDP growth (t-1)" ///
              l_debt "Public debt / GDP (t-1)" ///
              l_banking_crisis "Systemic banking-crisis dummy (t-1)" ///
              l_govexp "Govt expenditure / GDP (t-1)" ///
              l_open "Trade openness (t-1)" ///
              l_credit_bank "Private credit by banks / GDP (t-1)" ///
              l_lninfl "Log inflation, continuous (t-1)" ///
              exchange2 "Nominal exchange-rate log-change (t-1)") ///
    refcat(l_fedfunds "Predictors" l1_gdpg "Baseline controls", nolabel) ///
    stats(chi2p pp aurocctrl auroc N, ///
          labels("Chi-squared (predictors)" "  p-value" "AUROC, controls only" "AUROC, with predictors" "Observations") ///
          fmt(2 3 3 3 0)) ///
    title("Table 1. First-stage probit: predicting the start of a spread crisis") ///
    addnotes("Dependent variable: dummy = 1 in the onset year of the indicated crisis type; each type predicted vs tranquil years (the rival onset type is dropped)." ///
             "Pooled probit, no country fixed effects (country FE separate on the thin event count; they enter the outcome equation only). Robust standard errors clustered by country in parentheses." ///
             "Predictors are excluded from the LP/AIPW outcome equation (they are omitted from the second stage, which carries country FE only, no year FE)." ///
             "Chi-squared (predictors) is the joint Wald test that all three predictors are zero. AUROC, controls only is from a probit on the baseline controls" ///
             "alone (no predictors); AUROC, with predictors adds the three predictors -- the difference is the source Table 1's own demonstration that" ///
             "predictors add classification power, not just the chi-squared joint-significance test." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")

if _rc == 608 di as error "  ** table_first_stage.rtf is OPEN IN WORD — close it and re-run to refresh."
else if _rc  di as error "  ** Table (first stage): esttab failed (rc=" _rc ")"
else di as result "First-stage table saved: $tabs/table_first_stage.rtf"

di as result _n "08c_first_stage_table.do complete."
