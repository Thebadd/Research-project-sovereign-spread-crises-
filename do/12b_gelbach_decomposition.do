/*===========================================================================
  12B_GELBACH_DECOMPOSITION.DO
  Gelbach (2016) decomposition: how much of the estimated GDP cost of a
  spread-crisis onset is "explained" by each transmission channel?

  For each channel and horizon h=0..4, compare two NESTED regressions on the
  SAME sample (same-sample comparison is required for a valid decomposition),
  each estimated with Table 2's EXACT headline spec -- xtreg, country FE
  only (no year FE), plain robust SE -- since that headline number is what
  is being decomposed (see 02_lp_all.do/03_lp_resolution.do's header for the
  full argument: this is the Stata-idiomatic equivalent of the reference
  paper's own reg ..., vce(robust) noconstant with explicit country dummies):
    baseline: xtreg dy_h onset_nd onset_def $ctrl_core, fe vce(robust)
    full:     xtreg dy_h onset_nd onset_def ch_<channel>_h $ctrl_core, fe vce(robust)
  explained_nd  = (b_nd_base  - b_nd_full)  / b_nd_base
  explained_def = (b_def_base - b_def_full) / b_def_base

  Interpretation: the share of the onset coefficient's magnitude that is
  statistically "absorbed" once the channel's own contemporaneous change is
  added as a control. This is NOT a causal mediation estimate (the channel
  variable is not exogenous conditional on controls -- it could itself be
  driven by the same crisis shock, or GDP could drive the channel rather
  than the reverse), but it is a real, quantitative upgrade over reporting
  the channel's own single-outcome regression in isolation, since it
  directly ties the channel to the GDP cost it is meant to help explain.

  Caveat printed per cell: on the thin default-linked sample (14-21 treated
  episodes depending on channel coverage), adding a channel control can
  inflate se_def_full well beyond se_def_base -- a large explained_def share
  built on an unstable full-spec SE is flagged, not silently reported as if
  it were precise.

  Output: $tabs/gelbach_decomposition.csv
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

local controls $ctrl_core

* ── Rebuild channel outcomes (same construction as 11_channels.do) ─────────
* OUTCOME SCALE. Strictly-positive GDP-ratio channels use the LOG REAL LEVEL
* (ln_r_*, built in 18_transforms) so the outcome is a cumulative percent
* change in the variable itself rather than in its ratio to a collapsing
* GDP. pb and fdi change sign, so no log is possible and they keep the
* ratio form.
foreach var in credit claims_govt inv govexp pb fdi {
    local src `var'
    if inlist("`var'","credit","inv","govexp") local src ln_r_`var'
    capture drop `var'_base
    gen `var'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`src' - `var'_base
    }
}

local channels    credit claims_govt inv govexp pb fdi
local chan_labels `" "Private credit/GDP" "Bank claims on govt/GDP" "Investment/GDP" "Govt expenditure/GDP" "Primary balance/GDP" "FDI/GDP" "'

di as result _n "════════════════════════════════════════════════════════════"
di as result "GELBACH DECOMPOSITION: share of onset_nd/onset_def coefficient"
di as result "'explained' by adding each channel as a control"
di as result "════════════════════════════════════════════════════════════"
di as result "  baseline: xtreg dy_h onset_nd onset_def \$ctrl_core, fe vce(robust)"
di as result "  full:     xtreg dy_h onset_nd onset_def ch_<channel>_h \$ctrl_core, fe vce(robust)"
di as result "  (Table 2's exact headline spec; both estimated on the identical sample, restricted to non-missing channel data)"

tempfile gelbach_out
tempname R
postfile `R' str20 channel horizon double b_nd_base double b_nd_full double explained_nd ///
    double b_def_base double b_def_full double explained_def double se_def_base double se_def_full ///
    long nep_nd long nep_def using "`gelbach_out'", replace

local i = 1
foreach ch of local channels {
    local lbl : word `i' of `chan_labels'
    di as result _n "--- CHANNEL: `ch' (`lbl') ---"
    di "h   b_nd_base  b_nd_full  expl_nd%   b_def_base  b_def_full  expl_def%   nep_nd  nep_def"

    forvalues h = 0/4 {
        local hd  = `h' + 1

        quietly count if onset_nd==1  & sample==1 & !missing(ch_`ch'_`h')
        local nepnd  = r(N)
        quietly count if onset_def==1 & sample==1 & !missing(ch_`ch'_`h')
        local nepdef = r(N)

        * baseline and full spec estimated on the SAME sample (both restricted
        * to non-missing channel data) -- required for a valid decomposition.
        * Matches Table 2's headline spec exactly (xtreg, country FE only,
        * no year FE, robust SE) -- this is the number being decomposed, so
        * the baseline here must be identical to it, not the separate
        * cluster-SE convention used by the IPW/AIPW stages.
        capture xtreg dy_`h' onset_nd onset_def `controls' ///
            if sample==1 & !missing(ch_`ch'_`h'), fe vce(robust)
        if _rc {
            di as error "h=`hd': baseline regression failed for `ch' (rc=" _rc ")"
            continue
        }
        local bnd_base    = _b[onset_nd]
        local bdef_base   = _b[onset_def]
        local se_def_base = _se[onset_def]

        capture xtreg dy_`h' onset_nd onset_def ch_`ch'_`h' `controls' ///
            if sample==1 & !missing(ch_`ch'_`h'), fe vce(robust)
        if _rc {
            di as error "h=`hd': channel-controlled regression failed for `ch' (rc=" _rc ")"
            continue
        }
        local bnd_full    = _b[onset_nd]
        local bdef_full   = _b[onset_def]
        local se_def_full = _se[onset_def]

        * Explained share; guard against dividing by a near-zero baseline
        * coefficient (share is uninformative/explosive when b_base ~ 0).
        local expl_nd  = .
        local expl_def = .
        if abs(`bnd_base')  > 0.05 local expl_nd  = 100*(`bnd_base'  - `bnd_full')  / `bnd_base'
        if abs(`bdef_base') > 0.05 local expl_def = 100*(`bdef_base' - `bdef_full') / `bdef_base'

        * Flag an unstable full-spec SE on the default-linked cell (thin sample):
        * a >50% SE inflation is a warning that expl_def should not be over-read.
        local flag ""
        if `se_def_full' > 1.5*`se_def_base' local flag "  ** se_def inflated >50% -- expl_def unstable"

        post `R' ("`ch'") (`hd') (`bnd_base') (`bnd_full') (`expl_nd') ///
            (`bdef_base') (`bdef_full') (`expl_def') (`se_def_base') (`se_def_full') (`nepnd') (`nepdef')

        local expl_nd_s  = cond(missing(`expl_nd'),  "    .", string(`expl_nd',  "%6.1f"))
        local expl_def_s = cond(missing(`expl_def'), "    .", string(`expl_def', "%6.1f"))

        di "h=" `hd' "  " %8.3f `bnd_base' "  " %8.3f `bnd_full' "  " "`expl_nd_s'" ///
           "     " %8.3f `bdef_base' "  " %8.3f `bdef_full' "  " "`expl_def_s'" ///
           "     " %4.0f `nepnd' "    " %4.0f `nepdef' "`flag'"
    }
    local ++i
}
postclose `R'

* ── Export ───────────────────────────────────────────────────────────────
preserve
    use "`gelbach_out'", clear
    order channel horizon b_nd_base b_nd_full explained_nd b_def_base b_def_full explained_def se_def_base se_def_full nep_nd nep_def
    export delimited "$tabs/gelbach_decomposition.csv", replace
    di as result _n "Gelbach decomposition saved: $tabs/gelbach_decomposition.csv"
restore

di as result _n "════════════════════════════════════════════════════════════"
di as result "12b_gelbach_decomposition.do complete."
di as result "Read explained_nd/explained_def as the % of the onset coefficient"
di as result "absorbed by that channel -- NOT a causal mediation estimate."
di as result "Cells flagged '** se_def inflated' should be read as unstable,"
di as result "not as strong evidence either way."
di as result "════════════════════════════════════════════════════════════"
