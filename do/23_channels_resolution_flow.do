/*===========================================================================
  23_CHANNELS_RESOLUTION_FLOW.DO
  Transmission channels by resolution type — FLOW treatment, Spec A only

  WHY THIS FILE EXISTS
  --------------------
  12_channels_resolution.do splits each of the six channels from
  11_channels.do into non-default vs default-linked, under ONSET coding. This
  file re-runs its Spec A (joint OLS, both types entered together, tranquil
  omitted — the reference-paper baseline) on the FLOW treatment built in
  18_transforms.do and used in 20_lp_flow.do and 22_channels_flow.do.

  NOT BUILT: a flow version of 12's Spec B (IPW, per-type-vs-tranquil, two
  first stages). Same reason the two-stage tier stops at 21_aipw_flow.do: for
  a time-varying treatment the strongest predictor of being in crisis at t is
  being in crisis at t-1, a post-treatment outcome of the same episode, so a
  propensity model here is not well-posed (METHODOLOGY.md section 1). This
  file's results are read the way 20's and 22's are — magnitude and state
  conditional on selection into crisis, not a second identified effect.

  SPECIFICATION
  -------------
  ch_<var>_h = a_i + lambda_t + b_nd*in_crisis_nd + b_def*in_crisis_def
               + X_flow (episode-dated common core + channel's own
                 episode-dated pre-trend) + e
  Tranquil country-years are the omitted category. Country and year FE
  (single-stage rule, matching 20 and 22). DK SE, lag max(2,h+3).

  THE DIFFERENCE IS TESTED WITH LINCOM, NOT THE CLOGG z THAT 12 USES.
  Under flow coding a country with BOTH a non-default and a default-linked
  episode contributes treated rows to BOTH arms (10 countries do: Argentina,
  Bulgaria, Dominican Republic, Ecuador, Ghana, Lebanon, Russia, Sri Lanka,
  Ukraine, Uruguay), so b_nd and b_def are not independent and the Clogg
  formula sqrt(se_nd^2+se_def^2) is the wrong SE for their difference. lincom
  uses the estimated covariance from the joint regression, which is correct.
  Same reasoning, same departure, as 20_lp_flow.do's Section 3.

  Outputs
  -------
    "$tabs/table12_channels_resolution_flow.rtf"  one panel per channel
    "$tabs/channels_resolution_flow.csv"           raw coefficients
    "$clean/irf_flow_ch_<var>_nd.dta / _def.dta"   12 IRF datasets
    "$figs/fig12_channels_resolution_flow.pdf/.png"  2x3 grid, nd/def overlay
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
* 1. CHANNEL OUTCOMES — identical construction to 11/22
* ══════════════════════════════════════════════════════════════════════════
local channels credit claims_govt inv govexp pb fdi

foreach var of local channels {
    local src `var'
    if inlist("`var'","credit","inv","govexp") local src ln_r_`var'
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
}

local epc_lc epc_l_credit_bank
local epc_lg epc_l_govexp
local ctrl_credit      : list ctrl_flow - epc_lc
local ctrl_credit      `ctrl_credit' epc_pre_credit
local ctrl_claims_govt $ctrl_flow epc_pre_claims_govt
local ctrl_inv         $ctrl_flow epc_pre_inv
local ctrl_govexp      : list ctrl_flow - epc_lg
local ctrl_govexp      `ctrl_govexp' epc_pre_govexp
local ctrl_pb          $ctrl_flow epc_pre_pb
local ctrl_fdi         $ctrl_flow epc_pre_fdi

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
end

* ══════════════════════════════════════════════════════════════════════════
* 2. ESTIMATION — joint nd/def per channel, tranquil omitted
* ══════════════════════════════════════════════════════════════════════════
tempname F
tempfile flowf
postfile `F' str16 spec str24 channel int hdisp double b double se double p ///
    long nep long N using "`flowf'", replace

eststo clear
foreach ch of local channels {
    foreach g in nd def {
        foreach m in b lo90 hi90 lo95 hi95 {
            matrix `m'_`ch'_`g' = J(7, 1, .)
            matrix `m'_`ch'_`g'[2,1] = 0
        }
    }
}

di as result _n "════════════════════════════════════════════════════════════"
di as result "FLOW CHANNELS BY RESOLUTION TYPE (Year 1 = the crisis year)"
di as result "Difference via lincom (covariance-correct) — see file header"
di as result "════════════════════════════════════════════════════════════"

foreach ch of local channels {
    di as result _n "--- CHANNEL: `ch' ---"
    forvalues h = 0/4 {
        local hd  = `h' + 1
        local row = `h' + 3
        local lag = max(2, `h' + 3)

        capture xtscc ch_`ch'_`h' in_crisis_nd in_crisis_def `ctrl_`ch'' i.year ///
            if sample_flow==1, fe lag(`lag')
        if _rc {
            di as error "  h=`hd': failed for `ch' (rc=" _rc ")"
            continue
        }
        local bnd=_b[in_crisis_nd]
        local snd=_se[in_crisis_nd]
        local bdef=_b[in_crisis_def]
        local sdef=_se[in_crisis_def]
        _critvals
        local c90=r(c90)
        local c95=r(c95)
        _pval `bnd' `snd'
        local pnd=r(p)
        _pval `bdef' `sdef'
        local pdef=r(p)

        _nflowcount in_crisis_nd, outcome(ch_`ch'_`h') controls(`ctrl_`ch'') samp(sample_flow)
        local nend=r(nep)
        _nflowcount in_crisis_def, outcome(ch_`ch'_`h') controls(`ctrl_`ch'') samp(sample_flow)
        local nedef=r(nep)

        quietly lincom in_crisis_def - in_crisis_nd
        local bdiff=r(estimate)
        local sdiff=r(se)
        local dfd=r(df)
        if !missing(`dfd') & `dfd'>0 local pdi=2*ttail(`dfd', abs(`bdiff'/`sdiff'))
        else                         local pdi=2*(1-normal(abs(`bdiff'/`sdiff')))

        eststo t12f_`ch'_h`h'
        estadd scalar nepnd = `nend'
        estadd scalar nepdef = `nedef'
        estadd scalar bdiff = `bdiff'
        estadd scalar pdiff = `pdi'
        local elist_`ch' `elist_`ch'' t12f_`ch'_h`h'

        foreach g in nd def {
            matrix b_`ch'_`g'[`row',1]    = `b`g''
            matrix lo90_`ch'_`g'[`row',1] = `b`g'' - `c90'*`s`g''
            matrix hi90_`ch'_`g'[`row',1] = `b`g'' + `c90'*`s`g''
            matrix lo95_`ch'_`g'[`row',1] = `b`g'' - `c95'*`s`g''
            matrix hi95_`ch'_`g'[`row',1] = `b`g'' + `c95'*`s`g''
        }

        di as result "  h=" %1.0f `hd' "  ND=" %7.3f `bnd' "(" %5.3f `pnd' "," %2.0f `nend' "ep)" ///
            "  DEF=" %7.3f `bdef' "(" %5.3f `pdef' "," %2.0f `nedef' "ep)" ///
            "  diff=" %7.3f `bdiff' " p=" %5.3f `pdi'

        post `F' ("nd")   ("`ch'") (`hd') (`bnd')  (`snd')  (`pnd')  (`nend')  (e(N))
        post `F' ("def")  ("`ch'") (`hd') (`bdef')  (`sdef') (`pdef') (`nedef') (e(N))
        post `F' ("diff") ("`ch'") (`hd') (`bdiff') (`sdiff') (`pdi')  (.)      (e(N))
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 3. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
local ptitle_credit      "Table 12. Channels by resolution, flow specification -- Panel A: Private credit"
local ptitle_claims_govt "Panel B: Bank claims on govt/GDP"
local ptitle_inv         "Panel C: Investment"
local ptitle_govexp      "Panel D: Govt expenditure"
local ptitle_pb          "Panel E: Primary balance/GDP"
local ptitle_fdi         "Panel F: FDI/GDP"

local writemode replace
foreach ch of local channels {
    if "`elist_`ch''" == "" {
        di as error "  ** Table 12 panel `ch': no estimates stored — skipped"
        continue
    }
    local extra
    if "`ch'"=="fdi" local extra addnotes("Both treatment dummies entered jointly; tranquil country-years omitted. Difference and p-value from lincom (covariance-correct — see file header, NOT the Clogg z 12_channels_resolution.do uses). Flow treatment: in_crisis_nd/def = 1 in every year of the respective episode type. Country and year FE, DK SE lag max(2,h+3). No flow IPW column — see file header.")
    capture esttab `elist_`ch'' using "$tabs/table12_channels_resolution_flow.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(in_crisis_nd in_crisis_def) ///
        coeflabel(in_crisis_nd "Non-default crisis" in_crisis_def "Default-linked crisis") ///
        mtitles nonumber ///
        stats(bdiff pdiff nepnd nepdef N, ///
              labels("Difference (def-nd)" "p (lincom)" "Episodes, nd" "Episodes, def" "Observations") ///
              fmt(3 3 0 0 0)) ///
        title("`ptitle_`ch''") `extra'
    if _rc == 608 di as error "  ** table12_channels_resolution_flow.rtf OPEN IN WORD."
    else if _rc   di as error "  ** Table 12 panel `ch' failed (rc=" _rc ")"
    local writemode append
}
di as result "Table 12 saved: $tabs/table12_channels_resolution_flow.rtf"

postclose `F'
preserve
    use "`flowf'", clear
    label var spec "nd / def / diff"
    label var channel "Channel outcome"
    label var hdisp "Horizon (1 = crisis year)"
    label var nep "Distinct episodes"
    export delimited "$tabs/channels_resolution_flow.csv", replace
    di as result "Raw coefficients saved: $tabs/channels_resolution_flow.csv"
restore

* ── IRF datasets + figure ────────────────────────────────────────────────
foreach ch of local channels {
    foreach g in nd def {
        preserve
            clear
            set obs 7
            gen horizon = _n - 2
            foreach m in b lo90 hi90 lo95 hi95 {
                svmat `m'_`ch'_`g', names(`m')
                rename `m'1 `m'
            }
            gen series = "`g'"
            save "$clean/irf_flow_ch_`ch'_`g'.dta", replace
        restore
    }
}
di as result "IRF datasets saved: irf_flow_ch_<channel>_nd/_def.dta (12)"

local c_nd "0 84 166"
local c_def "157 36 73"
local c_zero "150 150 150"
local title_credit "Private credit"
local title_claims_govt "Bank claims on govt"
local title_inv "Investment"
local title_govexp "Govt expenditure"
local title_pb "Primary balance"
local title_fdi "FDI"

local i = 1
local gnames
foreach ch of local channels {
    preserve
        use "$clean/irf_flow_ch_`ch'_nd.dta", clear
        append using "$clean/irf_flow_ch_`ch'_def.dta"
        keep if horizon >= 0
        twoway ///
            (connected b horizon if series=="nd",  lcolor("`c_nd'")  mcolor("`c_nd'")  msymbol(circle) lwidth(medthick)) ///
            (connected b horizon if series=="def", lcolor("`c_def'") mcolor("`c_def'") msymbol(square) lpattern(dash) lwidth(medthick)), ///
            yline(0, lpattern(dash) lcolor("`c_zero'") lwidth(thin)) ///
            xlabel(0(1)5, labsize(small)) ylabel(, format(%4.1f) labsize(small)) ///
            xtitle("Year", size(small)) ytitle("", size(small)) ///
            title("`title_`ch''", size(small) color(navy)) legend(off) ///
            graphregion(color(white)) plotregion(color(white)) name(gr`i', replace)
    restore
    local gnames `gnames' gr`i'
    local ++i
}
graph combine `gnames', cols(3) ///
    title("Channels by Resolution — Flow Specification", size(medsmall) color(navy)) ///
    subtitle("Navy circles = non-default; brick squares = default-linked", size(small)) ///
    graphregion(color(white)) xsize(11) ysize(7)
capture graph export "$figs/fig12_channels_resolution_flow.pdf", replace
if _rc di as error "  ** fig12_channels_resolution_flow.pdf export failed (rc=" _rc ")"
else {
    capture graph export "$figs/fig12_channels_resolution_flow.png", replace width(1200)
    di as result "Figure 12 saved: fig12_channels_resolution_flow.pdf/.png"
}
foreach nm of local gnames {
    capture graph drop `nm'
}

di as result _n "23_channels_resolution_flow.do complete."
