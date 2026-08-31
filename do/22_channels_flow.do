/*===========================================================================
  22_CHANNELS_FLOW.DO
  Transmission channels of a sovereign spread crisis — FLOW treatment

  WHY THIS FILE EXISTS
  --------------------
  11_channels.do re-runs the Act-1 output LP on six intermediate outcomes
  (credit, claims on govt/GDP, investment, government spending, primary
  balance, FDI) under ONSET coding: one treated row per episode; two further
  sovereign-bank NEXUS channels (claimsgov_assets, claimpriv_assets, bank
  claims scaled by bank assets rather than GDP) are added the same way by
  11b_nexus_channels.do/13c_aipw_channels.do. Both sets are ported here. This
  file
  asks the same question 20_lp_flow.do asks of GDP — what does each channel do
  WHILE a country is in a spread crisis, not just after one starts — using the
  same flow treatment (in_crisis, built in 18_transforms.do; 234 treated
  country-years, onset and continuation alike).

  Outcome construction is IDENTICAL to 11_channels.do and is not rebuilt: the
  log-real-level form for credit/inv/govexp, the ratio form for pb/fdi (both
  change sign, so no log is possible), and pre_<var> = L.src - L2.src as each
  channel's own pre-crisis change. That construction is about the outcome
  variable, not the treatment, so nothing about flow coding touches it.

  WHAT DOES CHANGE, MATCHING 20_lp_flow.do EXACTLY
  --------------------------------------------------
  * Treatment: in_crisis, on sample_flow (not sample).
  * Controls: EPISODE-DATED. $ctrl_flow (built in 18_transforms.do) plus, for
    each channel, an episode-dated version of its own pre_<var> built HERE —
    pre_<var> is the channel's own pre-crisis change and is exactly as much a
    mediator for a continuation row as l1_gdpg is (a row three years into an
    episode has a pre_credit that reflects two years of the crisis's own
    effect on credit), so it needs the same entry-year dating. Built in-file
    rather than added to 18's global $ctrl_flow, because pre_<var> itself is
    only defined once 11's channel-construction logic runs, and this file
    already reproduces that logic — no reason to make stage 18 aware of what
    a "channel" is.
  * FE: country AND year. This is single-stage, so it follows the same rule as
    20 (METHODOLOGY.md section 2), not the two-stage rule.
  * DK lag: max(2, h+3), matching 20 — the regressor is serially correlated
    within an episode under flow coding, which onset coding does not have.
  * Robustness row r_noyearfe, matching 20: the no-year-FE version, reported so
    a reader can separate the estimator from the FE choice if this file is ever
    compared against a flow AIPW channel file.

  HORIZON CONVENTION — as in 20: Year 0 = baseline (zero by construction),
  Year 1 = the crisis year in question. Year 1 here is NOT the same object as
  Year 1 in Table 3 (11_channels.do): the onset one is the first year of every
  episode, the flow one pools the first year of one episode with the fourth
  year of another. See 20_lp_flow.do's header for the full argument.

  NOT BUILT: a flow IPW/AIPW channel file (the analogue of 12's Spec B or of
  13c). Same reason the two-stage tier stops at 21_aipw_flow.do: a propensity
  model for a time-varying treatment is not well-posed here (METHODOLOGY.md
  section 1, and the headers of 20 and 21). This file's channel results should
  therefore be read the same way 20's GDP results are read — magnitude and
  state conditional on selection into crisis, not a second identified effect.

  Outputs
  -------
    "$tabs/table11_channels_flow.rtf"   one panel per channel, flow spec
    "$tabs/channels_flow.csv"           raw coefficients, all channels, all variants
    "$clean/irf_flow_ch_<var>.dta"      one IRF dataset per channel (8)
    "$figs/fig11_channels_flow.pdf/.png"  3x3 grid (8 panels), drawn here (self-contained)
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

foreach v in in_crisis in_crisis_nd in_crisis_def sample_flow ep_seq {
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

* ══════════════════════════════════════════════════════════════════════════
* 1. CHANNEL OUTCOMES — identical construction to 11_channels.do
* ══════════════════════════════════════════════════════════════════════════
* claimsgov_assets/claimpriv_assets (the "nexus" channels, 11b_nexus_channels.do)
* are asset-share ratios (bank claims / bank assets), not GDP ratios, so like
* pb/fdi they take no ln_r_ transform -- they fall through to the raw `src'
* branch below unchanged. Coverage is thin: IMF nexus data start 2001 and 6
* panel countries are entirely absent (China, Ecuador, El Salvador, India,
* Lebanon, Vietnam -- see 01c_merge_nexus.do), worse again once L2. (for
* pre_<var>) and episode-averaging (for epc_pre_<var>) are layered on; the
* coverage report right after this loop makes that visible per channel.
local channels credit claims_govt inv govexp pb fdi claimsgov_assets claimpriv_assets

foreach var of local channels {
    local src `var'
    if inlist("`var'","credit","inv","govexp") local src ln_r_`var'
    capture drop `var'_base
    gen `var'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`src' - `var'_base
        label var ch_`var'_`h' "Cum. change `var': F`h' vs t-1 (flow)"
    }
    capture drop pre_`var'
    gen pre_`var' = L.`src' - L2.`src'

    * Episode-dated version of the channel's own pre-crisis change, same
    * construction as $ctrl_flow's epc_* terms in 18_transforms.do. Unaffected
    * by that file's annual-criterion redefinition of `in_crisis' -- the
    * mechanism is identical, it just now entry-dates a somewhat smaller set of
    * rows (the ~13 old gap-bridged rows are no longer in_crisis==1 at all).
    capture drop epc_pre_`var' _ent_pre_`var'
    quietly bysort cid ep_seq: egen double _ent_pre_`var' = ///
        max(cond(onset_all==1, pre_`var', .))
    quietly gen double epc_pre_`var' = cond(in_crisis==1, _ent_pre_`var', pre_`var')
    label var epc_pre_`var' "pre_`var' at episode entry (tranquil rows keep own t-1)"
    quietly drop _ent_pre_`var'
    * Restore panel order INSIDE the loop, not just after it: the preceding
    * bysort silently re-sorts the physical dataset by cid ep_seq, and the
    * NEXT iteration's `gen \`var'_base = L.\`src'' would then be computed on
    * data no longer sorted cid year. Belt and suspenders with the sort after
    * the loop closes.
    sort cid year
}

* bysort above re-sorts the physical dataset by its by-list as a side effect
* (independent of what xtset declared), leaving it sorted by cid ep_seq rather
* than cid year. xtscc does its own sort check and errors r(5) "not sorted" if
* this is not restored before the next panel command.
sort cid year

di as result _n "=== FLOW COVERAGE AT h=0 (sample_flow==1 & in_crisis==1) ==="
foreach var of local channels {
    quietly count if in_crisis==1 & sample_flow==1 & !missing(ch_`var'_0)
    di as result "  `var': " r(N) " / 234 treated rows with non-missing data at h=0"
}

* ══════════════════════════════════════════════════════════════════════════
* 2. CONTROL SETS — $ctrl_flow (episode-dated common core) + each channel's
*    own episode-dated pre-trend; drop the core term equal to the channel's
*    own lagged level, mirroring 11_channels.do exactly.
* ══════════════════════════════════════════════════════════════════════════
* credit drops epc_l_credit_bank (same reasoning as 11: correlates 0.950 with
* the credit outcome, close to an LDV); govexp drops epc_l_govexp. Locals must
* be defined BEFORE use in `: list A - B' -- an undefined local silently
* evaluates as empty, which would leave the term in rather than drop it.
local epc_lc epc_l_credit_bank
local epc_lg epc_l_govexp
local lc0    l_credit_bank
local lg0    l_govexp

* CONTROL SET. flow_ctrl_variant: 0 = the ADOPTED flow-tier baseline,
* $ctrl_flow, built in 18_transforms.do from $ctrl_core_flowbase
* (l_banking_duration -> the l_banking_crisis DUMMY, l_ca -> tot_chg -> exchange2,
* l_hyperinfl -> l_lninfl -- see that file's "ADOPTED FLOW-TIER CORE CONTROL
* SET"). 1 = the adopted core PLUS the reference paper's own additional
* predictors, tot_chg (the term exchange2 replaced) and l_imf (18_transforms.do's "ALTERNATE FLOW
* CONTROL SET"). Default 0. This drives every control-set reference below,
* INCLUDING the identity check in section 3 -- because the ADOPTED set is
* no longer term-for-term identical to $ctrl_core, that check no longer
* reproduces 11_channels.do's published credit-channel coefficient under
* the default; it is a pure self-consistency check (flow must still
* collapse to onset coding under whichever control set is active), not a
* literal match to that published figure.
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
if `flow_ctrl_variant'==1 local ctrl_flow_base $ctrl_flow_flowplus
else                       local ctrl_flow_base $ctrl_flow
if `flow_ctrl_variant'==1 local ctrl_core_base $ctrl_core_flowplus
else                       local ctrl_core_base $ctrl_core_flowbase

* EXPLORATORY: set to 1 to drop year FE and match the reference paper's
* single-stage rule (country FE only), matching 20_lp_flow.do's toggle of the
* same name. Default 0 = current baseline. Drives section 4's estimation
* loop only -- NOT the identity check above, which deliberately runs both
* variants side by side regardless of this setting.
local drop_year_fe 0
local yearfe = cond(`drop_year_fe', "", "i.year")

local ctrl_credit      : list ctrl_flow_base - epc_lc
local ctrl_credit      `ctrl_credit' epc_pre_credit
local ctrl_claims_govt `ctrl_flow_base' epc_pre_claims_govt
local ctrl_inv         `ctrl_flow_base' epc_pre_inv
local ctrl_govexp      : list ctrl_flow_base - epc_lg
local ctrl_govexp      `ctrl_govexp' epc_pre_govexp
local ctrl_pb          `ctrl_flow_base' epc_pre_pb
local ctrl_fdi         `ctrl_flow_base' epc_pre_fdi
* Nexus channels: no core term dropped, same convention as 11b/13c (they are
* not literally identical to any $ctrl_core regressor, unlike credit/govexp).
local ctrl_claimsgov_assets `ctrl_flow_base' epc_pre_claimsgov_assets
local ctrl_claimpriv_assets `ctrl_flow_base' epc_pre_claimpriv_assets

local ctrl_row_credit      : list ctrl_core_base - lc0
local ctrl_row_credit      `ctrl_row_credit' pre_credit
local ctrl_row_claims_govt `ctrl_core_base' pre_claims_govt
local ctrl_row_inv         `ctrl_core_base' pre_inv
local ctrl_row_govexp      : list ctrl_core_base - lg0
local ctrl_row_govexp      `ctrl_row_govexp' pre_govexp
local ctrl_row_pb          `ctrl_core_base' pre_pb
local ctrl_row_fdi         `ctrl_core_base' pre_fdi

capture program drop _critvals
program define _critvals, rclass
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 {
        return scalar c90 = 1.645
        return scalar c95 = 1.960
    }
    else {
        return scalar c90 = invttail(`dfr', 0.05)
        return scalar c95 = invttail(`dfr', 0.025)
    }
end

capture program drop _pval
program define _pval, rclass
    args b se
    tempname dfr
    scalar `dfr' = e(df_r)
    if missing(`dfr') | `dfr' <= 0 return scalar p = 2*(1 - normal(abs(`b'/`se')))
    else                           return scalar p = 2*ttail(`dfr', abs(`b'/`se'))
end

capture program drop _nflowcount
program define _nflowcount, rclass
    syntax varname(numeric) , Outcome(varname) Controls(varlist) Samp(varname)
    tempvar esmp tagep tagcty
    quietly gen byte `esmp' = e(sample)
    quietly count if `esmp' == 1
    if r(N) == 0 {
        quietly replace `esmp' = (`samp' == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar nrow = r(N)
    quietly egen byte `tagep' = tag(cid ep_seq) if `varlist'==1 & `esmp'==1
    quietly count if `tagep'==1
    return scalar nep = r(N)
    quietly egen byte `tagcty' = tag(cid) if `varlist'==1 & `esmp'==1
    quietly count if `tagcty'==1
    return scalar ncty = r(N)
end

sort cid year   // belt-and-suspenders: guaranteed panel order before the first xtscc call

* ══════════════════════════════════════════════════════════════════════════
* 3. IDENTITY CHECK — credit channel, on sample==1. Same logic as 20's
*    Section 2: on sample==1 continuation rows are excluded, in_crisis IS
*    onset_all, and epc_* collapses to row-dated, so flow and onset must
*    return the SAME coefficient under whichever control set is active.
*    Under the ADOPTED default this is a pure self-consistency check, NOT a
*    literal match to 11_channels.do's published Table 3 coefficient (see
*    the `flow_ctrl_variant' comment above) -- the adopted set is no longer
*    term-for-term identical to $ctrl_core. Run once (credit); the same
*    construction governs the other channels.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "3. IDENTITY CHECK (credit channel) — flow on sample==1 must equal 11_channels"
di as result "════════════════════════════════════════════════════════════"

quietly xtscc ch_credit_0 in_crisis `ctrl_credit' if sample==1, fe lag(1)
local f0 = _b[in_crisis]
quietly xtscc ch_credit_0 onset_all `ctrl_row_credit' if sample==1, fe lag(1)
local o0 = _b[onset_all]
di as result "  no year FE   flow: " %9.6f `f0' "   onset: " %9.6f `o0'
if abs(`f0'-`o0') > 1e-6 di as error "  ** no-year-FE anchor FAILED (credit)."

quietly xtscc ch_credit_0 in_crisis `ctrl_credit' i.year if sample==1, fe lag(1)
local f0y = _b[in_crisis]
quietly xtscc ch_credit_0 onset_all `ctrl_row_credit' i.year if sample==1, fe lag(1)
local o0y = _b[onset_all]
di as result "  with year FE flow: " %9.6f `f0y' "   onset: " %9.6f `o0y' ///
             "   (self-consistency under the ACTIVE control set; not a literal" ///
             "    match to 11_channels.do's published Table 3 under the default)"
if abs(`f0y'-`o0y') > 1e-6 di as error "  ** with-year-FE anchor FAILED (credit)."
else                       di as result "  MATCH — treatment, sample flag and control dating are correct."

sort cid year   // belt-and-suspenders: guaranteed panel order before the estimation loop

* ══════════════════════════════════════════════════════════════════════════
* 4. ESTIMATION — one loop over all six channels, pooled flow LP
* ══════════════════════════════════════════════════════════════════════════
tempname F
tempfile flowf
postfile `F' str16 spec str24 channel int hdisp double b double se double p ///
    long nrow long nep long ncty long N using "`flowf'", replace

eststo clear
foreach ch of local channels {
    foreach m in b lo90 hi90 lo95 hi95 {
        matrix `m'_`ch' = J(7, 1, .)
        matrix `m'_`ch'[2,1] = 0
    }
}

di as result _n "════════════════════════════════════════════════════════════"
di as result "4. FLOW CHANNEL LP (Year 1 = the crisis year), country and year FE"
di as result "════════════════════════════════════════════════════════════"

foreach ch of local channels {
    di as result _n "--- CHANNEL: `ch' ---"
    forvalues h = 0/4 {
        local hd  = `h' + 1
        local row = `h' + 3
        local lag = max(2, `h' + 3)

        capture xtscc ch_`ch'_`h' in_crisis `ctrl_`ch'' `yearfe' if sample_flow==1, fe lag(`lag')
        if _rc {
            di as error "  h=`hd': failed for `ch' (rc=" _rc ")"
            continue
        }
        local bb = _b[in_crisis]
        local ss = _se[in_crisis]
        _pval `bb' `ss'
        local pp = r(p)
        _critvals
        local c90=r(c90)
        local c95=r(c95)
        _nflowcount in_crisis, outcome(ch_`ch'_`h') controls(`ctrl_`ch'') samp(sample_flow)
        local nr=r(nrow)
        local ne=r(nep)
        local nc=r(ncty)

        eststo t11f_`ch'_h`h'
        estadd scalar nep = `ne'
        estadd scalar ncty = `nc'
        local elist_`ch' `elist_`ch'' t11f_`ch'_h`h'

        matrix b_`ch'[`row',1]=`bb'
        matrix lo90_`ch'[`row',1]=`bb'-`c90'*`ss'
        matrix hi90_`ch'[`row',1]=`bb'+`c90'*`ss'
        matrix lo95_`ch'[`row',1]=`bb'-`c95'*`ss'
        matrix hi95_`ch'[`row',1]=`bb'+`c95'*`ss'

        di as result "  h=" %1.0f `hd' "  b=" %7.3f `bb' "  SE=" %6.3f `ss' ///
            "  p=" %5.3f `pp' "  rows=" %4.0f `nr' "  ep=" %3.0f `ne' "  cty=" %3.0f `nc'
        post `F' ("pooled") ("`ch'") (`hd') (`bb') (`ss') (`pp') (`nr') (`ne') (`nc') (e(N))

        * robustness: no year FE
        capture quietly xtscc ch_`ch'_`h' in_crisis `ctrl_`ch'' if sample_flow==1, fe lag(`lag')
        if _rc==0 post `F' ("r_noyearfe") ("`ch'") (`hd') (_b[in_crisis]) (_se[in_crisis]) ///
            (2*(1-normal(abs(_b[in_crisis]/_se[in_crisis])))) (.) (.) (.) (e(N))
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 5. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
local ptitle_credit      "Table 11. Transmission channels, flow specification -- Panel A: Private credit"
local ptitle_claims_govt "Panel B: Bank claims on govt/GDP"
local ptitle_inv         "Panel C: Investment"
local ptitle_govexp      "Panel D: Govt expenditure"
local ptitle_pb          "Panel E: Primary balance/GDP"
local ptitle_fdi         "Panel F: FDI/GDP"
local ptitle_claimsgov_assets "Panel G: Bank claims on govt / bank assets (nexus, doom-loop)"
local ptitle_claimpriv_assets "Panel H: Bank claims on private sector / bank assets (nexus, reallocation)"

local writemode replace
foreach ch of local channels {
    if "`elist_`ch''" == "" {
        di as error "  ** Table 11 panel `ch': no estimates stored — skipped"
        continue
    }
    local extra
    if "`ch'"=="fdi" local extra addnotes("Dependent variable: cumulative change in the channel variable from t-1 to t+h under FLOW treatment (in_crisis=1 in every year of an episode). Year 1 = the crisis year, NOT the same object as Year 1 in Table 3 (11_channels.do) -- see file header. Episode-dated controls plus the channel's own episode-dated pre-crisis change. Country and year FE. Driscoll-Kraay SE, lag max(2,h+3). Episodes/countries reported, not just observations. * p<0.10, ** p<0.05, *** p<0.01.")
    capture esttab `elist_`ch'' using "$tabs/table11_channels_flow.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(in_crisis) coeflabel(in_crisis "In a spread crisis") ///
        mtitles nonumber ///
        stats(nep ncty N, labels("Episodes" "Countries" "Observations") fmt(0 0 0)) ///
        title("`ptitle_`ch''") `extra'
    if _rc == 608 di as error "  ** table11_channels_flow.rtf OPEN IN WORD — close and re-run."
    else if _rc   di as error "  ** Table 11 panel `ch' failed (rc=" _rc ")"
    local writemode append
}
di as result "Table 11 saved: $tabs/table11_channels_flow.rtf"

postclose `F'
preserve
    use "`flowf'", clear
    label var spec "Specification variant"
    label var channel "Channel outcome"
    label var hdisp "Horizon (1 = crisis year)"
    label var nep "Distinct episodes"
    label var ncty "Distinct countries"
    label var nrow "Treated country-years"
    export delimited "$tabs/channels_flow.csv", replace
    di as result "Raw coefficients saved: $tabs/channels_flow.csv"
restore

* ── IRF datasets + figure ────────────────────────────────────────────────
foreach ch of local channels {
    preserve
        clear
        set obs 7
        gen horizon = _n - 2
        foreach m in b lo90 hi90 lo95 hi95 {
            svmat `m'_`ch', names(`m')
            rename `m'1 `m'
        }
        save "$clean/irf_flow_ch_`ch'.dta", replace
    restore
}
di as result "IRF datasets saved: irf_flow_ch_<channel>.dta (6)"

local c_flow "23 55 94"
local c_zero "150 150 150"
local title_credit "Private credit"
local title_claims_govt "Bank claims on govt"
local title_inv "Investment"
local title_govexp "Govt expenditure"
local title_pb "Primary balance"
local title_fdi "FDI"
local title_claimsgov_assets "Nexus: claims on govt"
local title_claimpriv_assets "Nexus: claims on private"

local i = 1
local gnames
foreach ch of local channels {
    preserve
        use "$clean/irf_flow_ch_`ch'.dta", clear
        keep if horizon >= 0
        twoway ///
            (rarea lo95 hi95 horizon, color("`c_flow'%15") lwidth(none)) ///
            (rarea lo90 hi90 horizon, color("`c_flow'%25") lwidth(none)) ///
            (connected b horizon, lcolor("`c_flow'") mcolor("`c_flow'") msymbol(circle) lwidth(medthick)), ///
            yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
            xlabel(0(1)5, labsize(small)) ylabel(, format(%4.1f) labsize(small)) ///
            xtitle("Year", size(small)) ytitle("", size(small)) ///
            title("`title_`ch''", size(small) color(navy)) legend(off) ///
            graphregion(color(white)) plotregion(color(white)) name(gf`i', replace)
    restore
    local gnames `gnames' gf`i'
    local ++i
}
graph combine `gnames', cols(3) ///
    title("Transmission Channels — Flow Specification", size(medsmall) color(navy)) ///
    subtitle("In a spread crisis vs tranquil years; 90/95% CIs; country and year FE", size(small)) ///
    graphregion(color(white)) xsize(11) ysize(7)
capture graph export "$figs/fig11_channels_flow.pdf", replace
if _rc di as error "  ** fig11_channels_flow.pdf export failed (rc=" _rc ")"
else {
    capture graph export "$figs/fig11_channels_flow.png", replace width(1200)
    di as result "Figure 11 saved: fig11_channels_flow.pdf/.png"
}
foreach nm of local gnames {
    capture graph drop `nm'
}

di as result _n "22_channels_flow.do complete."
