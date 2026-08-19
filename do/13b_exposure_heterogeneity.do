/*===========================================================================
  13B_EXPOSURE_HETEROGENEITY.DO
  Tier-3 channel evidence: exposure-interaction tests

  IDEA (model-free)
  -----------------
  Showing "crisis -> credit falls" and "crisis -> GDP falls" separately does
  NOT prove credit is the channel through which GDP falls: both could be
  falling for a common reason. A channel leaves a sharper empirical footprint:
  a crisis should depress output MORE in countries that were MORE EXPOSED to
  that channel before the crisis. We test this by interacting crisis onset
  with a pre-determined (pre-crisis) measure of exposure to each channel.

  PART A — POOLED (all 61 onsets)
      dy_{i,t+h} = a_i + l_t + b_h*onset_all
                   + d_h*(onset_all x Z) + phi*Z + controls + e

  PART B — BY RESOLUTION TYPE (non-default vs default-linked)
      dy_{i,t+h} = a_i + l_t + b_nd*onset_nd + b_def*onset_def
                   + d_nd*(onset_nd x Z) + d_def*(onset_def x Z)
                   + phi*Z + controls + e
      with the test  H0: d_nd = d_def  at each horizon.

  The OUTCOME is always the SAME cumulative real-GDP change as the main
  LP (dy_h). Z is the standardized pre-crisis exposure to a channel.

  READING THE SIGN OF d -- AND THE CONFOUND THAT MAKES IT INTERPRETABLE ONLY
  FOR SOME EXPOSURES (revised against the estimates, not assumed in advance)
  ---------------------------------------------------------------------------
  This file used to state a flat prior: "d < 0 means the channel transmits."
  That was never estimated anywhere in the project -- it is an assumption about
  what the exposure VARIABLE measures, and the run rejected it for most of them.
  What the project HAS established is that credit contracts after a crisis
  (11_channels: -7.2pp at Year 3) and that credit absorbs part of the GDP cost
  (12b Gelbach: 43% of the default-linked coefficient by Year 4). Neither
  implies that a higher credit/GDP RATIO means greater vulnerability, and that
  missing step is where the sign prediction came from.

  The interaction d is interpretable ONLY IF the exposure has a single
  economically sensible sign. Most of these do not:

    AMBIGUOUS (ratio proxies financial development as much as vulnerability)
      credit, inv, fdi, claims_govt, debt_service
      - more bank-dependent / capital-intensive  => amplifies      (d < 0)
      - deeper, more developed financial system  => cushions       (d > 0)
      The literature is genuinely two-sided (financial accelerator vs financial
      development vs "too much finance"), so neither sign refutes anything.
      ESTIMATES: credit d = +2.33 (h3, p=.001), +3.87 (h5, p=.001);
      inv +2.06 (h5, p=.033); fdi +1.33 (h3, p=.026) -- all POSITIVE, i.e. the
      development reading dominates. Corroborated inside our own results:
      l_credit_bank enters the Act-1 outcome equation NEGATIVELY and growing
      (-0.026 at h1 to -0.171 at h5, t=-3.98), so high-credit countries grow
      more slowly on trend AND lose less in a crisis -- the signature of an
      income/development proxy, not of a credit channel.
      => Report these as UNINFORMATIVE about transmission. Do not read the
         positive coefficient as "credit protects"; read it as "this measure
         cannot separate exposure from development."

    UNAMBIGUOUS (only one sign is economically sensible)
      stdebt_share   short-term share of external debt. No development reading:
                     a country is not "more developed" for having more debt
                     maturing next year. Pure rollover vulnerability.
                     ESTIMATE: d = -4.70 (h2, p=.000), -3.19 (h3, p=.005).
                     => the one clean dose-response result in the file.
      claimsgov_assets  share of the BANK BALANCE SHEET committed to the
                     sovereign. A portfolio share, not a size/depth measure, so
                     it does not scale with development the way GDP ratios do.
                     Null pooled; the result is in Part B (see below).

  PART B adds the resolution split. Two findings survived the control-set
  alignment and one did not:
      claimsgov_assets  d_def = -4.43 (h3, p(diff)=.034), -4.27 (h4, p=.029),
                        d_nd ~ 0 -- the sovereign-bank nexus amplifies the
                        output cost specifically in default-linked episodes.
      stdebt_share      d_def = -6.20 (h3, p(diff)=.053), -5.77 (h4, p=.033).
      claims_govt/GDP   previously d_def = -7.56 (h1, p(diff)=.006) on the OLD
                        legacy control set; on $ctrl_core the SE triples
                        (2.58 -> 7.33) and p goes to .246. That result was an
                        artifact of an underspecified model and is NOT used.

  MULTIPLE TESTING. 9 exposures x 5 horizons x 2 panels = 90 tests; at 5% about
  4-5 significant results are expected by chance and roughly that many appear.
  With 11-20 treated episodes per default cell no single panel here can carry a
  headline claim. Use this file for the two unambiguous exposures above, and
  read the nexus result as one leg of a three-design pattern (13d median split,
  13b continuous interaction, 13c AIPW channel gap) rather than as standalone
  proof.

  EXPOSURE MEASURES (one per channel mechanism)
  ---------------------------------------------
    credit           Private credit / GDP            AMBIGUOUS (depth confound)
    claims_govt      Bank claims on govt / GDP       AMBIGUOUS (GDP-scaled)
    claimsgov_assets Bank claims on govt / assets    UNAMBIGUOUS -- the
                                                     sovereign-bank NEXUS
                                                     PARAMETER (portfolio share
                                                     of the bank balance sheet)
    inv              Investment / GDP                AMBIGUOUS (depth confound)
    fdi              FDI / GDP                       AMBIGUOUS (depth confound)

  DESIGN CHOICES
  --------------
  * Exposure is PRE-DETERMINED: measured at t-1 (L.var), so it cannot be an
    outcome of the crisis it is supposed to amplify.
  * Exposure is STANDARDIZED (mean 0, SD 1) over the estimation sample, so d
    is "effect per 1 SD of exposure" and b is the onset effect at AVERAGE
    exposure (directly comparable to the Table 1 / Table 2 baselines).
  * Same controls, country + year FE, Driscoll-Kraay SE, continuation years
    excluded -- identical to the main LP, so only the interaction is new.
  * NOTE ON COVERAGE: claimsgov_assets (nexus) is IMF-based, 2001-2024 and
    missing 6 countries, so its onset sample is smaller; the by-type split is
    especially thin there. Coverage is printed below and low-N panels are
    handled gracefully (skipped with a warning, never a halt).

  OUTPUTS
  -------
    "$tabs/table5_exposure_interactions.rtf"       Word: pooled (Part A)
    "$tabs/table6_exposure_by_resolution.rtf"      Word: by type + p(diff) (Part B)
    "$tabs/exposure_interactions.csv"              raw pooled coefficients
    "$tabs/exposure_by_resolution.csv"             raw by-type coefficients
    "$figs/fig13f_exposure_interactions.pdf"       pooled d_h, grid
    "$figs/fig13g_exposure_by_resolution.pdf"      d_nd vs d_def, grid
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

* ── Exposure variables and readable labels ───────────────────────────────
* Default-mechanism exposures (access loss / balance-sheet / capital):
*   credit, claims_govt, claimsgov_assets, inv, fdi
* Non-default-mechanism exposures (cost / rollover):
*   claimpriv_assets (bank private-lending share -> lending-rate pass-through)
*   stdebt_share     (short-term debt / ext. debt -> rollover risk)      [01d]
*   debt_service     (debt service / exports      -> liquidity burden)   [01d]
*   intpay_gni       (interest payments / GNI      -> avg interest cost)  [01d]
* Note: stdebt_share, debt_service and intpay_gni all proxy the same external-
* debt refinancing-burden construct and are collinear — read them as ALTERNATIVE
* measures (robustness), not independent channels.
* (reserves_extdebt is merged by 01d but excluded here: the reserves/ext-debt
*  measure gave a noisy, wrong-signed default estimate on a tiny default cell.)
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
* any that have not been merged yet, so this file always runs.
local expok
foreach e of local expvars {
    capture confirm variable `e', exact
    if _rc == 0 local expok `expok' `e'
    else di as error "  ** exposure `e' not in panel yet — skipped (run 01d to merge it)"
}
local expvars `expok'

* ── Baseline controls (identical to the main output LP, 02_lp_all.do) ─────
* VIX and ust10y are absorbed by year FE; excluded as elsewhere.
local controls $ctrl_core

* ══════════════════════════════════════════════════════════════════════════
* 1. BUILD PRE-DETERMINED, STANDARDIZED EXPOSURES AND INTERACTIONS
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== PRE-CRISIS EXPOSURE COVERAGE AT ONSET ==="
foreach e of local expvars {
    capture drop exp_`e' z_`e' Dz_`e' Dnd_`e' Ddef_`e'

    * Pre-crisis level: value at t-1 (predetermined w.r.t. onset at t)
    gen exp_`e' = L.`e'

    * Standardize over the estimation sample (mean 0, SD 1)
    quietly summarize exp_`e' if sample == 1
    gen z_`e' = (exp_`e' - r(mean)) / r(sd)
    label var z_`e' "Std. pre-crisis `lbl_`e''"

    * Interactions: onset (pooled) and onset-by-type x standardized exposure
    gen Dz_`e'   = onset_all * z_`e'
    gen Dnd_`e'  = onset_nd  * z_`e'
    gen Ddef_`e' = onset_def * z_`e'
    label var Dz_`e'   "Onset x `lbl_`e''"
    label var Dnd_`e'  "Onset(non-default) x `lbl_`e''"
    label var Ddef_`e' "Onset(default) x `lbl_`e''"

    quietly count if onset_all == 1 & sample == 1 & !missing(z_`e')
    local n_all = r(N)
    quietly count if onset_nd  == 1 & sample == 1 & !missing(z_`e')
    local n_nd  = r(N)
    quietly count if onset_def == 1 & sample == 1 & !missing(z_`e')
    local n_def = r(N)
    di as result "  `e': all=" `n_all' "  non-default=" `n_nd' "  default=" `n_def'
}

* ══════════════════════════════════════════════════════════════════════════
* CONFOUND CHECK — is the exposure sorting on DEVELOPMENT rather than on
* vulnerability to the channel?
*
* This is the diagnostic 13d already runs on its nexus median split and that
* this file was missing. It is the check that would have caught the positive
* credit/inv/fdi interactions before they were read as a channel result.
*
* For each exposure: split onsets at the median of Z and report mean pre-crisis
* real GDP per capita in each bin. If the high-exposure bin is systematically
* richer, the interaction is picking up income/development, not exposure to the
* channel, and its sign carries no information about transmission.
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
di as result "  read as channel exposure. Compare with 13d, where the nexus split passes"
di as result "  this check (high bin holds Brazil/Mexico/Turkey AND Ghana/Kenya/Zambia)."

* ══════════════════════════════════════════════════════════════════════════
* PART A — POOLED INTERACTION (onset_all x exposure)
*   Coefficient of interest: Dz_`e' (= d_h). Negative => channel amplifies
*   the output loss where pre-crisis exposure was higher.
* ══════════════════════════════════════════════════════════════════════════

* Storage for the pooled figure: d_h with 90% CI
foreach e of local expvars {
    foreach m in dcoef dlo90 dhi90 {
        matrix `m'_`e' = J(6, 1, 0)
    }
}

eststo clear

foreach e of local expvars {

    di as result _n "========================================================"
    di as result "PART A — EXPOSURE: `lbl_`e'' (pooled)"
    di as result "========================================================"
    di as result "h    b_onset   SE       d_interact  SE       p(d)"

    local elistA_`e'

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local row = `h' + 2

        capture xtscc dy_`h' onset_all z_`e' Dz_`e' `controls' i.year ///
            if sample == 1, fe lag(`lag')

        if _rc == 0 {
            eststo a_`e'_`h', title("h=`=`h'+1'")
            local elistA_`e' `elistA_`e'' a_`e'_`h'

            matrix dcoef_`e'[`row',1] = _b[Dz_`e']
            matrix dlo90_`e'[`row',1] = _b[Dz_`e'] - 1.645*_se[Dz_`e']
            matrix dhi90_`e'[`row',1] = _b[Dz_`e'] + 1.645*_se[Dz_`e']

            local pd = 2*(1 - normal(abs(_b[Dz_`e']/_se[Dz_`e'])))
            di "h=" `h'+1 "   " %7.3f _b[onset_all] "  " %6.3f _se[onset_all] ///
               "   " %7.3f _b[Dz_`e'] "  " %6.3f _se[Dz_`e'] ///
               "   " %5.3f `pd'
        }
        else di as error "h=" `h'+1 ": xtscc failed for exposure `e' (rc=" _rc ")"
    }
}

di as result _n "Interpretation (Part A) -- READ THE SIGN ONLY WHERE IT IS IDENTIFIED:"
di as result "  d_h < 0 => output falls MORE where pre-crisis exposure was higher."
di as result "  BUT this is evidence of transmission only for an exposure with a single"
di as result "  economically sensible sign. credit/inv/fdi/claims_govt/debt_service are"
di as result "  GDP ratios that proxy financial DEPTH as much as vulnerability, so a"
di as result "  positive d there means the development reading dominates -- it is NOT"
di as result "  evidence the channel protects, and NOT evidence against the channel."
di as result "  stdebt_share is the clean case (no development reading): d<0 = rollover risk."

* ── TABLE 5: pooled exposure interactions (Word/RTF) ─────────────────────
local t5note "Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h (same as Table 1). Each column adds one channel's standardized pre-crisis exposure and its interaction with crisis onset. 'Onset x exposure' is the effect per 1 SD of pre-crisis exposure; a negative value means the output loss is deeper where exposure was higher. Exposure measured at t-1. Country and year fixed effects; continuation years excluded. Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01. A negative coefficient is evidence of channel transmission ONLY for an exposure with a single economically sensible sign. credit, investment, FDI, bank claims on government and debt service are GDP ratios that proxy financial development as much as vulnerability, so their positive coefficients indicate the development reading dominates and are uninformative about transmission. The interpretable measures are short-term debt share (pure rollover risk, no development reading) and bank claims on government / bank assets (a portfolio share rather than a size measure). This file runs 9 exposures x 5 horizons x 2 panels = 90 tests; at the 5% level roughly 4-5 significant cells are expected by chance. No single panel should be read as decisive. "

local panel A
local writemode replace
local t5fail 0

foreach e of local expvars {
    if "`elistA_`e''" == "" {
        di as error "  ** Table 5: no estimates for exposure `e' — panel skipped"
        local t5fail 1
        continue
    }
    local t5extra
    if "`e'" == "fdi" local t5extra addnotes("`t5note'")

    capture esttab `elistA_`e'' using "$tabs/table5_exposure_interactions.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_all Dz_`e') order(onset_all Dz_`e') ///
        coeflabel(onset_all "Spread-crisis onset (at mean exposure)" ///
                  Dz_`e' "Onset x pre-crisis exposure (per SD)") ///
        mtitles nonumber ///
        stats(N N_g, labels("Observations" "Countries") fmt(0 0)) ///
        title("Table 5. Exposure heterogeneity (pooled) -- Panel `panel': `lbl_`e''") ///
        `t5extra'

    if _rc == 608 {
        di as error "  ** table5_exposure_interactions.rtf is OPEN IN WORD — close it and re-run."
        local t5fail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 5: esttab failed for exposure `e' (rc=" _rc ")"
        local t5fail 1
        continue
    }
    local writemode append
    if      "`panel'"=="A" local panel B
    else if "`panel'"=="B" local panel C
    else if "`panel'"=="C" local panel D
    else if "`panel'"=="D" local panel E
}
if `t5fail' == 0 di as result "Table 5 saved: $tabs/table5_exposure_interactions.rtf"
else di as error "Table 5 written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* PART B — INTERACTION BY RESOLUTION TYPE (non-default vs default-linked)
*   d_nd  = onset_nd  x exposure   d_def = onset_def x exposure
*   Test H0: d_nd = d_def  -> is the channel's amplification stronger in one
*   crisis type? (e.g., is the nexus exposure the default-specific margin?)
* ══════════════════════════════════════════════════════════════════════════

* Storage for the by-type figure: d_nd and d_def with 90% CI
foreach e of local expvars {
    foreach m in bnd blo_nd bhi_nd bdef blo_def bhi_def pdiff {
        matrix `m'_`e' = J(6, 1, 0)
    }
}

eststo clear

foreach e of local expvars {

    di as result _n "========================================================"
    di as result "PART B — EXPOSURE: `lbl_`e'' (by resolution type)"
    di as result "========================================================"
    di as result "h    d_nd      SE       d_def     SE       p(d_nd=d_def)"

    local elistB_`e'

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local row = `h' + 2

        capture xtscc dy_`h' onset_nd onset_def z_`e' Dnd_`e' Ddef_`e' ///
            `controls' i.year if sample == 1, fe lag(`lag')

        if _rc == 0 {
            * Equality test on the two amplification terms
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
di as result "  Same sign caveat as Part A. The two interpretable results are"
di as result "  claimsgov_assets (nexus, d_def<0 at h3-h4) and stdebt_share (rollover,"
di as result "  d_def<0 at h3-h4). With 90 tests in this file, ~4-5 significant cells are"
di as result "  expected by chance: treat any single panel as suggestive, not decisive."

* ── TABLE 6: exposure interactions by resolution type (Word/RTF) ─────────
local t6note "Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h. Both onset dummies and both exposure interactions enter jointly. 'Onset(type) x exposure' is the extra output effect per 1 SD of pre-crisis exposure in that crisis type. p(nd=def) is the p-value of the equality of the two interaction terms. Exposure measured at t-1. Country and year fixed effects; continuation years excluded. Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01. A negative coefficient is evidence of channel transmission ONLY for an exposure with a single economically sensible sign. credit, investment, FDI, bank claims on government and debt service are GDP ratios that proxy financial development as much as vulnerability, so their positive coefficients indicate the development reading dominates and are uninformative about transmission. The interpretable measures are short-term debt share (pure rollover risk, no development reading) and bank claims on government / bank assets (a portfolio share rather than a size measure). This file runs 9 exposures x 5 horizons x 2 panels = 90 tests; at the 5% level roughly 4-5 significant cells are expected by chance. No single panel should be read as decisive. "

local panel A
local writemode replace
local t6fail 0

foreach e of local expvars {
    if "`elistB_`e''" == "" {
        di as error "  ** Table 6: no estimates for exposure `e' — panel skipped"
        local t6fail 1
        continue
    }
    local t6extra
    if "`e'" == "fdi" local t6extra addnotes("`t6note'")

    capture esttab `elistB_`e'' using "$tabs/table6_exposure_by_resolution.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(Dnd_`e' Ddef_`e') order(Dnd_`e' Ddef_`e') ///
        coeflabel(Dnd_`e' "Onset(non-default) x exposure (per SD)" ///
                  Ddef_`e' "Onset(default) x exposure (per SD)") ///
        mtitles nonumber ///
        stats(pdiff N N_g, labels("p (nd = def)" "Observations" "Countries") fmt(3 0 0)) ///
        title("Table 6. Exposure amplification by resolution -- Panel `panel': `lbl_`e''") ///
        `t6extra'

    if _rc == 608 {
        di as error "  ** table6_exposure_by_resolution.rtf is OPEN IN WORD — close it and re-run."
        local t6fail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 6: esttab failed for exposure `e' (rc=" _rc ")"
        local t6fail 1
        continue
    }
    local writemode append
    if      "`panel'"=="A" local panel B
    else if "`panel'"=="B" local panel C
    else if "`panel'"=="C" local panel D
    else if "`panel'"=="D" local panel E
}
if `t6fail' == 0 di as result "Table 6 saved: $tabs/table6_exposure_by_resolution.rtf"
else di as error "Table 6 written with warnings (see messages above)."

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
    export delimited "$tabs/exposure_interactions.csv", replace
    di as result "Pooled exposure CSV saved: $tabs/exposure_interactions.csv"
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
    export delimited "$tabs/exposure_by_resolution.csv", replace
    di as result "By-type exposure CSV saved: $tabs/exposure_by_resolution.csv"
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
            ytitle("Onset x exposure (pp/SD)", size(small)) ///
            title("`lbl_`e''", size(small) color(navy)) ///
            legend(off) graphregion(color(white)) plotregion(color(white)) ///
            name(gA_`i', replace)
    restore
    local fignamesA `fignamesA' gA_`i'
    local ++i
}
graph combine `fignamesA', cols(3) ///
    title("Exposure Heterogeneity (Pooled): Output Loss by Pre-Crisis Exposure", ///
          size(medsmall) color(navy)) ///
    note("Interaction d_h (onset x standardized pre-crisis exposure), 90% CI. Below zero => output falls more where exposure was higher." ///
         "Driscoll-Kraay SE. Country & year FE. Same GDP outcome as the main LP.", size(vsmall)) ///
    graphregion(color(white)) xsize(11) ysize(7)
graph export "$figs/fig13f_exposure_interactions.pdf", replace
di as result "Figure saved: fig13f_exposure_interactions.pdf"
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
            ytitle("Onset x exposure (pp/SD)", size(small)) ///
            title("`lbl_`e''", size(small) color(navy)) ///
            legend(off) graphregion(color(white)) plotregion(color(white)) ///
            name(gB_`i', replace)
    restore
    local fignamesB `fignamesB' gB_`i'
    local ++i
}
graph combine `fignamesB', cols(3) ///
    title("Exposure Amplification by Resolution Type", size(medsmall) color(navy)) ///
    subtitle("Green triangles = non-default; red squares = default-linked", size(small)) ///
    note("Interaction d (onset-by-type x standardized pre-crisis exposure), 90% CI. More-negative in a type => stronger amplification there." ///
         "Driscoll-Kraay SE. Country & year FE. Same GDP outcome as the main LP.", size(vsmall)) ///
    graphregion(color(white)) xsize(11) ysize(7)
graph export "$figs/fig13g_exposure_by_resolution.pdf", replace
di as result "Figure saved: fig13g_exposure_by_resolution.pdf"
foreach nm of local fignamesB {
    capture graph drop `nm'
}

di as result _n "13b_exposure_heterogeneity.do complete."
