/*===========================================================================
  02A_DESCRIPTIVE_FACTS.DO
  Stylised facts: the raw GDP and transmission-channel paths around a spread
  crisis, before any conditioning. Runs BEFORE the local projections and is
  the figure the Data section of the paper leads with.

  WHY THIS EXISTS
  ---------------
  Every estimate in this paper is conditional: country and year fixed effects,
  eight predetermined controls, an estimator chosen for its properties. A
  reader is entitled to see the object of interest before being asked to
  accept a specification. If the default gap is only visible after the
  controls go in, that is worth knowing; if it is visible in the raw data, the
  regressions are sharpening a fact rather than manufacturing one.

  This is the analog of the reference paper's descriptive tier (their Fig. 1/
  Section 2.2), which reports conditional MEANS of the same horizon-
  differenced outcome by restructuring type, with no regression at all. We
  follow that construction, extended from GDP alone to the transmission
  channels this project actually estimates with AIPW (13c_aipw_channels.do):
  bank credit, investment, claims on government, and the two sovereign-bank
  nexus shares. UNLIKE the reference paper's Figure 1, this file does NOT
  cover gross capital inflows or real lending interest rates -- this project
  has no source data for either, and building it would be new data
  construction, not a figure extension; that scope decision is stated here
  rather than left implicit.

  WHAT IS AND IS NOT DONE TO THE DATA
  -----------------------------------
  DONE: every outcome is country-demeaned. A cumulative-change outcome
  reflects systematic cross-country differences in trend, so a raw average
  across episodes would partly rank countries rather than describe crises.
  Subtracting each country's own mean over the estimation sample removes
  that and leaves the deviation from a country's normal path -- exactly what
  a country fixed effect does in the regressions, and the same demeaning the
  reference paper applies. Each channel outcome (ch_v_h = F h.v - L.v) is
  built here with the identical construction used in 13c_aipw_channels.do,
  since ch_v_h is not persisted in panel_lp.dta and this file runs before 13c.

  NOT DONE: no controls, no year effects, no weighting, no estimator. The
  plotted point at horizon h for group g is the simple average of the
  country-demeaned outcome over that group's onsets. No standard errors are
  shown, deliberately: this is a description, and attaching inference to it
  would invite reading it as a result.

  HORIZONS: same convention as the rest of the paper -- Year 1 is the crisis
  year, Year 0 is the hard-coded pre-crisis baseline, Year -1 is the single
  estimable pre-crisis point (on the same t-1 base as every other horizon).

  Output: $figs/fig0_descriptive_paths.pdf     (GDP, nd vs def -- the Data-section figure)
          $figs/fig0a_descriptive_all.pdf      (GDP, all episodes pooled)
          $figs/fig0_descriptive_<channel>.pdf (standalone figure, one per channel -- each
                                                 channel is its own independent graph, not
                                                 combined into a multi-panel figure)
          $tabs/descriptive_paths.csv          (all plotted numbers, GDP + channels)
          $tabs/descriptive_summary.csv        (pre-crisis characteristics by group)
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 0. CHANNEL OUTCOMES ch_v_h = F h.v - L.v (h=0..4), IDENTICAL CONSTRUCTION TO
*    13c_aipw_channels.do -- built fresh here since this file runs early in
*    the pipeline (before 13c) and ch_v_h is not persisted in panel_lp.dta.
*    Scope: the channels this project actually has data for and estimates
*    with AIPW (credit, inv, claims_govt, claimsgov_assets, claimpriv_assets)
*    -- NOT the reference paper's full Figure 1 (gross capital inflows, real
*    lending rates), for which this project has no source data. See this
*    file's header note on that scope decision.
* ══════════════════════════════════════════════════════════════════════════
foreach v in credit inv claims_govt claimsgov_assets claimpriv_assets {
    local src `v'
    if inlist("`v'","credit","inv") local src ln_r_`v'
    capture drop `v'_base
    gen double `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen double ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    * Year -1 placebo point, same t-1 base as every other horizon (mirrors
    * dy_m2's construction for GDP: change from t-2 to the t-1 base).
    capture drop ch_`v'_m2
    gen double ch_`v'_m2 = L2.`src' - `v'_base
}

* ══════════════════════════════════════════════════════════════════════════
* 1. COUNTRY-DEMEAN THE OUTCOMES
*    Demeaning is over sample==1 (onset + tranquil years, continuation and
*    carry-in rows excluded) so the reference path is a country's behaviour
*    in the same universe the regressions use, not one that includes the
*    continuation years the design deliberately drops.
* ══════════════════════════════════════════════════════════════════════════
foreach h in m2 0 1 2 3 4 {
    capture drop cmean_`h' dd_`h'
    quietly bysort cid: egen double cmean_`h' = mean(dy_`h') if sample==1
    quietly gen double dd_`h' = dy_`h' - cmean_`h' if sample==1
}
label var dd_0 "Country-demeaned cumulative GDP change, crisis year"

foreach v in credit inv claims_govt claimsgov_assets claimpriv_assets {
    foreach h in m2 0 1 2 3 4 {
        capture drop cmean_`v'_`h' dd_`v'_`h'
        quietly bysort cid: egen double cmean_`v'_`h' = mean(ch_`v'_`h') if sample==1
        quietly gen double dd_`v'_`h' = ch_`v'_`h' - cmean_`v'_`h' if sample==1
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 2. GROUP MEANS BY HORIZON  (rows: 1 = Year -1, 2 = Year 0, 3..7 = Years 1-5)
* ══════════════════════════════════════════════════════════════════════════
foreach g in all nd def {
    matrix desc_`g' = J(7, 1, .)
    matrix nobs_`g' = J(7, 1, .)
    matrix desc_`g'[2,1] = 0          // Year 0 = the baseline, zero by construction
    matrix nobs_`g'[2,1] = .
}

* Year -1: the pre-crisis placebo point, same t-1 base as every other horizon
foreach g in all nd def {
    quietly summarize dd_m2 if onset_`g'==1 & sample==1
    matrix desc_`g'[1,1] = r(mean)
    matrix nobs_`g'[1,1] = r(N)
}

* Years 1-5
forvalues h = 0/4 {
    local row = `h' + 3
    foreach g in all nd def {
        quietly summarize dd_`h' if onset_`g'==1 & sample==1
        matrix desc_`g'[`row',1] = r(mean)
        matrix nobs_`g'[`row',1] = r(N)
    }
}

di as result _n "════════════════════════════════════════════════════════════"
di as result "DESCRIPTIVE PATHS — country-demeaned mean cumulative GDP change"
di as result "No controls, no fixed effects beyond the demeaning, no estimator."
di as result "════════════════════════════════════════════════════════════"
di as result "Year      All (n)        Non-default (n)     Default-linked (n)"
forvalues r = 1/7 {
    local yr = `r' - 2
    if `r' == 2 {
        di "  " %2.0f `yr' "     0.000  (base)      0.000  (base)        0.000  (base)"
    }
    else {
        di "  " %2.0f `yr' "  " %8.3f desc_all[`r',1] " (" %3.0f nobs_all[`r',1] ")" ///
           "  " %8.3f desc_nd[`r',1]  " (" %3.0f nobs_nd[`r',1] ")" ///
           "    " %8.3f desc_def[`r',1] " (" %3.0f nobs_def[`r',1] ")"
    }
}
di as result _n "  Read the Year -1 row as the descriptive counterpart of the placebo"
di as result "  test: if the two groups were already diverging before the crisis it"
di as result "  shows up here, with no specification standing between the reader and"
di as result "  the data."

* ══════════════════════════════════════════════════════════════════════════
* 2b. CHANNEL GROUP MEANS BY HORIZON -- identical construction to Section 2,
*     applied to each active channel (credit, inv, claims_govt,
*     claimsgov_assets, claimpriv_assets). No regression, no controls.
* ══════════════════════════════════════════════════════════════════════════
foreach v in credit inv claims_govt claimsgov_assets claimpriv_assets {
    foreach g in all nd def {
        matrix descv_`v'_`g' = J(7, 1, .)
        matrix nobsv_`v'_`g' = J(7, 1, .)
        matrix descv_`v'_`g'[2,1] = 0
        matrix nobsv_`v'_`g'[2,1] = .
    }
    foreach g in all nd def {
        quietly summarize dd_`v'_m2 if onset_`g'==1 & sample==1
        matrix descv_`v'_`g'[1,1] = r(mean)
        matrix nobsv_`v'_`g'[1,1] = r(N)
    }
    forvalues h = 0/4 {
        local row = `h' + 3
        foreach g in all nd def {
            quietly summarize dd_`v'_`h' if onset_`g'==1 & sample==1
            matrix descv_`v'_`g'[`row',1] = r(mean)
            matrix nobsv_`v'_`g'[`row',1] = r(N)
        }
    }

    di as result _n "════════════════════════════════════════════════════════════"
    di as result "DESCRIPTIVE PATHS — country-demeaned mean cumulative `v' change"
    di as result "No controls, no fixed effects beyond the demeaning, no estimator."
    di as result "════════════════════════════════════════════════════════════"
    di as result "Year      All (n)        Non-default (n)     Default-linked (n)"
    forvalues r = 1/7 {
        local yr = `r' - 2
        if `r' == 2 {
            di "  " %2.0f `yr' "     0.000  (base)      0.000  (base)        0.000  (base)"
        }
        else {
            di "  " %2.0f `yr' "  " %8.3f descv_`v'_all[`r',1] " (" %3.0f nobsv_`v'_all[`r',1] ")" ///
               "  " %8.3f descv_`v'_nd[`r',1]  " (" %3.0f nobsv_`v'_nd[`r',1] ")" ///
               "    " %8.3f descv_`v'_def[`r',1] " (" %3.0f nobsv_`v'_def[`r',1] ")"
        }
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 3. PLOT DATASET + FIGURES
* ══════════════════════════════════════════════════════════════════════════
preserve
    clear
    set obs 7
    gen horizon = _n - 2
    foreach g in all nd def {
        svmat desc_`g', names(b_`g')
        rename b_`g'1 b_`g'
        svmat nobs_`g', names(n_`g')
        rename n_`g'1 n_`g'
    }
    foreach v in credit inv claims_govt claimsgov_assets claimpriv_assets {
        foreach g in all nd def {
            svmat descv_`v'_`g', names(b_`v'_`g')
            rename b_`v'_`g'1 b_`v'_`g'
            svmat nobsv_`v'_`g', names(n_`v'_`g')
            rename n_`v'_`g'1 n_`v'_`g'
        }
    }
    label var horizon "Year (Year 1 = crisis year)"
    export delimited "$tabs/descriptive_paths.csv", replace
    di as result "Descriptive paths saved: $tabs/descriptive_paths.csv"

    local c_nd  "0 84 166"
    local c_def "157 36 73"
    local c_all "23 55 94"

    * ── Figure 0: the Data-section figure — nd vs def, raw ────────────────
    twoway ///
        (connected b_nd horizon, ///
            lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) lwidth(medthick)) ///
        (connected b_def horizon, ///
            lcolor("`c_def'") mcolor("`c_def'") msymbol(square) ///
            lpattern(dash) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xline(0.5, lpattern(solid) lcolor(gs11) lwidth(thin)) ///
        xlabel(-1(1)5, labsize(medsmall)) ///
        ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year (Year 0 = pre-crisis baseline, Year 1 = crisis year)", size(small)) ///
        ytitle("Cumulative change in log real GDP (pp)", size(small)) ///
        title("Output around a spread crisis, by resolution", size(medium) color(navy)) ///
        subtitle("Country-demeaned means. No controls, no estimator.", size(small)) ///
        legend(order(1 "Non-default" 2 "Default-linked") size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    * Note text (kept as source comment, no longer rendered on the figure --
    * the legend, previously overlapping the plot at ring(0) pos(7), now
    * takes the bottom position this note used to occupy):
    * "Simple average of the country-demeaned cumulative change in log real GDP over each group's onsets.
    *  Demeaning is within country over the estimation sample and removes cross-country differences in trend
    *  growth; nothing else is done to the data. No confidence bands are shown because this is a description,
    *  not an estimate — the conditional versions with inference are Figure 2 and Table 2."
    graph export "$figs/fig0_descriptive_paths.pdf", replace
    capture graph export "$figs/fig0_descriptive_paths.png", replace width(1200)
    di as result "Figure saved: fig0_descriptive_paths.pdf"

    * ── Figure 0a: pooled, for the motivating paragraph ───────────────────
    twoway (connected b_all horizon, ///
            lcolor("`c_all'") mcolor("`c_all'") msymbol(circle) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xline(0.5, lpattern(solid) lcolor(gs11) lwidth(thin)) ///
        xlabel(-1(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year (Year 0 = pre-crisis baseline, Year 1 = crisis year)", size(small)) ///
        ytitle("Cumulative change in log real GDP (pp)", size(small)) ///
        title("Output around a spread crisis, all episodes", size(medium) color(navy)) ///
        subtitle("Country-demeaned means. No controls, no estimator.", size(small)) ///
        legend(off) ///
        note("All 61 identified onsets. Construction as in Figure 0.", size(vsmall)) ///
        graphregion(color(white)) plotregion(color(white))
    graph export "$figs/fig0a_descriptive_all.pdf", replace
    di as result "Figure saved: fig0a_descriptive_all.pdf"

    * ── Figures 0b-0f: channel descriptive paths, same construction as
    * Figure 0, one standalone figure per active channel (GDP already has
    * its own standalone Figure 0/0a above). Each channel gets its own
    * independent graph/PDF -- no combined small-multiple panel.
    local panellab_credit "Bank credit"
    local panellab_inv "Investment"
    local panellab_claims_govt "Claims on govt"
    local panellab_claimsgov_assets "Claims/assets (govt)"
    local panellab_claimpriv_assets "Claims/assets (private)"
    foreach v in credit inv claims_govt claimsgov_assets claimpriv_assets {
        twoway ///
            (connected b_`v'_nd horizon, lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) lwidth(medthick)) ///
            (connected b_`v'_def horizon, lcolor("`c_def'") mcolor("`c_def'") msymbol(square) lpattern(dash) lwidth(medthick)), ///
            yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) xline(0.5, lpattern(solid) lcolor(gs11) lwidth(thin)) ///
            xlabel(-1(1)5, labsize(small)) ylabel(, format(%4.1f) labsize(small)) ///
            xtitle("Years relative to crisis onset") ytitle("`v' (pp, country-demeaned)", size(small)) ///
            title("`panellab_`v''", size(medium)) ///
            legend(order(1 "Non-default" 2 "Default-linked") pos(6) size(small) rows(1)) ///
            name(gk_`v', replace) graphregion(color(white)) plotregion(color(white))
        graph export "$figs/fig0_descriptive_`v'.pdf", replace name(gk_`v')
        graph drop gk_`v'
    }
    di as result "Figures saved: fig0_descriptive_<channel>.pdf, one standalone graph per channel"
restore

* ══════════════════════════════════════════════════════════════════════════
* 4. PRE-CRISIS CHARACTERISTICS BY GROUP
*    The Data section needs a table showing what the two groups look like
*    going in. Everything is measured at t-1, i.e. the same predetermined
*    values the regressions condition on, so the table describes exactly the
*    variation the controls are asked to absorb. Tranquil years are included
*    as the reference column because they are the control group throughout.
* ══════════════════════════════════════════════════════════════════════════
local sumvars l1_gdpg l_debt l_banking_crisis l_govexp l_open ///
              l_credit_bank l_lninfl exchange2 l_spr_mean claimsgov_assets

di as result _n "════════════════════════════════════════════════════════════"
di as result "PRE-CRISIS CHARACTERISTICS (all at t-1) BY GROUP"
di as result "════════════════════════════════════════════════════════════"
di as result %-22s "Variable" "   Tranquil    Non-default   Default    p(nd=def)"

tempname S
tempfile sumf
postfile `S' str24 variable double m_tranq double m_nd double m_def double pdiff ///
    long n_tranq long n_nd long n_def using "`sumf'", replace

foreach v of local sumvars {
    capture confirm variable `v'
    if _rc continue

    quietly summarize `v' if sample==1 & onset_all==0
    local mt = r(mean)
    local nt = r(N)
    quietly summarize `v' if sample==1 & onset_nd==1
    local mn = r(mean)
    local nn = r(N)
    quietly summarize `v' if sample==1 & onset_def==1
    local md = r(mean)
    local nd = r(N)

    * two-sample t-test on the nd/def difference; the paper's identifying
    * comparison is between these two groups, so that is the difference worth
    * testing here rather than crisis-vs-tranquil.
    local pd = .
    capture ttest `v' if sample==1 & onset_all==1, by(nondefault)
    if _rc == 0 local pd = r(p)

    post `S' ("`v'") (`mt') (`mn') (`md') (`pd') (`nt') (`nn') (`nd')
    di as result %-22s "`v'" "  " %9.2f `mt' "  " %9.2f `mn' "  " %9.2f `md' "   " %6.3f `pd'
}
postclose `S'

preserve
    use "`sumf'", clear
    label var m_tranq "Mean, tranquil years"
    label var m_nd    "Mean, non-default onsets"
    label var m_def   "Mean, default-linked onsets"
    label var pdiff   "p-value, nd = def"
    export delimited "$tabs/descriptive_summary.csv", replace
    di as result _n "Summary table saved: $tabs/descriptive_summary.csv"
restore

di as result _n "  A significant p(nd=def) means the two groups entered their crises"
di as result "  differently on that dimension, which is the selection the propensity"
di as result "  models of Section 5 are built to address. Read this table alongside"
di as result "  the first-stage probit rather than as a balance test that has to pass."

di as result _n "02a_descriptive_facts.do complete."
