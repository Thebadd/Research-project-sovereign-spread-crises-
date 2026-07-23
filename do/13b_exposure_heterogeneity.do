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
  with a pre-determined (pre-crisis) measure of exposure to each channel:

      dy_{i,t+h} = a_i + l_t + b_h * onset_all
                   + d_h * (onset_all x Z_exposure_it)
                   + phi * Z_exposure_it + controls + e

  where the OUTCOME is the SAME cumulative GDP-per-capita change used in the
  main LP (dy_h), and Z_exposure is the standardized pre-crisis exposure.

  A negative and significant d_h means: a one-standard-deviation higher
  pre-crisis exposure to the channel makes the output loss DEEPER. That is
  direct, cross-sectional evidence that the channel transmits the shock to
  output -- much harder to dismiss as a common cause than a simple IRF.

  EXPOSURE MEASURES (one per channel mechanism)
  ---------------------------------------------
    credit       Private credit / GDP        -> credit / bank-dependence channel
    claims_govt  Bank claims on govt / GDP   -> sovereign-bank nexus channel
    inv          Investment / GDP            -> capital-accumulation channel
    fdi          FDI / GDP                   -> external-financing channel

  DESIGN CHOICES
  --------------
  * Exposure is PRE-DETERMINED: measured at t-1 (L.var), so it cannot be an
    outcome of the crisis it is supposed to amplify.
  * Exposure is STANDARDIZED (mean 0, SD 1) over the estimation sample, so d_h
    is "effect per 1 SD of exposure" and b_h is the onset effect at AVERAGE
    exposure (directly comparable to the Table 1 baseline).
  * Same controls, country + year FE, Driscoll-Kraay SE, continuation years
    excluded -- identical to the main LP, so only the interaction is new.

  OUTPUTS
  -------
    "$tabs/table5_exposure_interactions.rtf"  -- Word table (one panel/exposure)
    "$tabs/exposure_interactions.csv"         -- raw coefficients
    "$figs/fig13f_exposure_interactions.pdf"  -- d_h across horizons, 2x2 grid
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ── Exposure variables and readable labels ───────────────────────────────
local expvars    credit claims_govt inv fdi
local lbl_credit      "Private credit/GDP"
local lbl_claims_govt "Bank claims on govt/GDP"
local lbl_inv         "Investment/GDP"
local lbl_fdi         "FDI/GDP"

* ── Baseline controls (identical to the main output LP, 02_lp_all.do) ─────
* VIX and ust10y are absorbed by year FE; excluded as elsewhere.
local controls l1_gdpg l2_gdpg ca debt infl imf

* ══════════════════════════════════════════════════════════════════════════
* 1. BUILD PRE-DETERMINED, STANDARDIZED EXPOSURES AND INTERACTIONS
* ══════════════════════════════════════════════════════════════════════════

foreach e of local expvars {
    capture drop exp_`e' z_`e' Dz_`e'

    * Pre-crisis level: value at t-1 (predetermined w.r.t. onset at t)
    gen exp_`e' = L.`e'

    * Standardize over the estimation sample (mean 0, SD 1)
    quietly summarize exp_`e' if sample == 1
    gen z_`e' = (exp_`e' - r(mean)) / r(sd)
    label var z_`e' "Std. pre-crisis `lbl_`e''"

    * Interaction: onset x standardized exposure
    gen Dz_`e' = onset_all * z_`e'
    label var Dz_`e' "Onset x `lbl_`e''"

    quietly count if onset_all == 1 & sample == 1 & !missing(z_`e')
    di as result "  `e': " r(N) " onset obs with non-missing pre-crisis exposure"
}

* ══════════════════════════════════════════════════════════════════════════
* 2. ESTIMATE OUTPUT LP WITH EXPOSURE INTERACTION, BY EXPOSURE AND HORIZON
*    Coefficient of interest: Dz_`e' (= d_h). Negative => channel amplifies
*    the output loss where pre-crisis exposure was higher.
* ══════════════════════════════════════════════════════════════════════════

* Storage for the figure: d_h with 90% CI, one 5-row column per exposure
foreach e of local expvars {
    foreach m in dcoef dlo90 dhi90 {
        matrix `m'_`e' = J(5, 1, .)
    }
}

eststo clear

foreach e of local expvars {

    di as result _n "========================================================"
    di as result "EXPOSURE: `lbl_`e'' (z_`e')"
    di as result "  d_h = onset x exposure interaction (per 1 SD)"
    di as result "========================================================"
    di as result "h    b_onset   SE       d_interact  SE       p(d)"

    local elist_`e'

    forvalues h = 0/4 {
        local lag = max(1, `h'+1)
        local row = `h' + 1

        capture xtscc dy_`h' onset_all z_`e' Dz_`e' `controls' i.year ///
            if sample == 1, fe lag(`lag')

        if _rc == 0 {
            eststo x_`e'_`h', title("h=`h'")
            local elist_`e' `elist_`e'' x_`e'_`h'

            matrix dcoef_`e'[`row',1] = _b[Dz_`e']
            matrix dlo90_`e'[`row',1] = _b[Dz_`e'] - 1.645*_se[Dz_`e']
            matrix dhi90_`e'[`row',1] = _b[Dz_`e'] + 1.645*_se[Dz_`e']

            local pd = 2*(1 - normal(abs(_b[Dz_`e']/_se[Dz_`e'])))
            di "h=" `h' "   " %7.3f _b[onset_all] "  " %6.3f _se[onset_all] ///
               "   " %7.3f _b[Dz_`e'] "  " %6.3f _se[Dz_`e'] ///
               "   " %5.3f `pd'
        }
        else di as error "h=" `h' ": xtscc failed for exposure `e' (rc=" _rc ")"
    }
}

di as result _n "Interpretation:"
di as result "  d_h < 0 and significant => output falls MORE where pre-crisis"
di as result "  exposure to the channel was higher: evidence the channel transmits."

* ══════════════════════════════════════════════════════════════════════════
* 3. TABLE EXPORT — TABLE 5: Exposure interactions (Word/RTF)
*    One panel per exposure; columns = horizons. Rows: onset (effect at mean
*    exposure) and onset x exposure (the amplification term). Each esttab is
*    wrapped in capture so a locked file or missing estimate warns, not halts.
*    Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

local t5note "Dependent variable: cumulative change in log real GDP per capita (pp) from t-1 to t+h (same as Table 1). Each column adds one channel's standardized pre-crisis exposure and its interaction with crisis onset. 'Onset x exposure' is the effect per 1 SD of pre-crisis exposure; a negative value means the output loss is deeper where exposure was higher. Exposure measured at t-1. Country and year fixed effects; continuation years excluded. Driscoll-Kraay standard errors in parentheses. * p<0.10, ** p<0.05, *** p<0.01."

local panel A
local writemode replace
local t5fail 0

foreach e of local expvars {

    if "`elist_`e''" == "" {
        di as error "  ** Table 5: no estimates for exposure `e' — panel skipped"
        local t5fail 1
        continue
    }

    * Methodology note only on the last panel
    local t5extra
    if "`e'" == "fdi" local t5extra addnotes("`t5note'")

    capture esttab `elist_`e'' using "$tabs/table5_exposure_interactions.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_all Dz_`e') order(onset_all Dz_`e') ///
        coeflabel(onset_all "Spread-crisis onset (at mean exposure)" ///
                  Dz_`e' "Onset x pre-crisis exposure (per SD)") ///
        mtitles nonumber ///
        stats(N N_g, labels("Observations" "Countries") fmt(0 0)) ///
        title("Table 5. Exposure heterogeneity -- Panel `panel': `lbl_`e''") ///
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
    * advance panel letter A -> B -> C -> D
    if "`panel'" == "A"      local panel B
    else if "`panel'" == "B" local panel C
    else if "`panel'" == "C" local panel D
    else if "`panel'" == "D" local panel E
}

if `t5fail' == 0 di as result "Table 5 saved: $tabs/table5_exposure_interactions.rtf"
else di as error "Table 5 written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* 4. RAW CSV OF INTERACTION COEFFICIENTS (d_h) FOR ALL EXPOSURES
* ══════════════════════════════════════════════════════════════════════════

preserve
    clear
    local nobs = 5 * 4     // 5 horizons x 4 exposures
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
            replace horizon  = `h'                in `row'
            replace d_coef   = dcoef_`e'[`h'+1,1] in `row'
            replace d_lo90   = dlo90_`e'[`h'+1,1] in `row'
            replace d_hi90   = dhi90_`e'[`h'+1,1] in `row'
            local ++row
        }
    }
    order exposure horizon d_coef d_lo90 d_hi90
    export delimited "$tabs/exposure_interactions.csv", replace
    di as result "Exposure interactions CSV saved: $tabs/exposure_interactions.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 5. FIGURE — INTERACTION COEFFICIENT d_h ACROSS HORIZONS (2x2 GRID)
*    One panel per exposure. A path below zero => amplification.
* ══════════════════════════════════════════════════════════════════════════

local c_amp "157 36 73"

local fignames_h
local i = 1
foreach e of local expvars {

    preserve
        clear
        set obs 5
        gen horizon = _n - 1
        foreach m in dcoef dlo90 dhi90 {
            svmat `m'_`e', names(`m')
            rename `m'1 `m'
        }

        twoway ///
            (rarea dlo90 dhi90 horizon, color("`c_amp'%20") lwidth(none)) ///
            (connected dcoef horizon, ///
                lcolor("`c_amp'") lwidth(medthick) msymbol(circle) mcolor("`c_amp'")), ///
            yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
            xlabel(0(1)4, labsize(medsmall)) ///
            ylabel(, format(%4.1f) labsize(medsmall)) ///
            xtitle("Years after onset", size(small)) ///
            ytitle("Onset x exposure (pp per SD)", size(small)) ///
            title("`lbl_`e''", size(medsmall) color(navy)) ///
            legend(off) ///
            graphregion(color(white)) plotregion(color(white)) ///
            name(g_exp_`i', replace)
    restore

    local fignames_h `fignames_h' g_exp_`i'
    local ++i
}

graph combine `fignames_h', ///
    cols(2) rows(2) ///
    title("Exposure Heterogeneity: Output Loss by Pre-Crisis Channel Exposure", ///
          size(medium) color(navy)) ///
    note("Interaction coefficient d_h (onset x standardized pre-crisis exposure) with 90% CI." ///
         "Below zero => output falls more where pre-crisis exposure to the channel was higher." ///
         "Driscoll-Kraay SE. Country & year FE. Same GDP outcome as the main LP.", size(vsmall)) ///
    graphregion(color(white)) xsize(9) ysize(7)

graph export "$figs/fig13f_exposure_interactions.pdf", replace
di as result "Figure saved: fig13f_exposure_interactions.pdf"

foreach nm of local fignames_h {
    capture graph drop `nm'
}

di as result _n "13b_exposure_heterogeneity.do complete."
