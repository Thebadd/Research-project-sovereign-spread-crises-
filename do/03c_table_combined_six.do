/*===========================================================================
  03C_TABLE_COMBINED_SIX.DO
  Single merged onset-tier OLS table, all SIX of the project's headline
  outcome variables (GDP, Bank credit, Investment, Bank claims on
  government/GDP, FDI, Real lending rate) as panels in ONE RTF -- Act 2,
  matching 03_lp_resolution.do's own CURRENT headline GDP spec
  (table2_output_resolution.rtf) and 12_channels_resolution.do's own
  channel spec (table4_channels_resolution.rtf) EXACTLY, panel for panel,
  with no change to either specification.

  WHY A SEPARATE FILE, NOT JUST ESTTAB-ING ACROSS 03 AND 12: confirmed
  this session that 12_channels_resolution.do's own `eststo clear' (its
  own line 121) wipes out 03_lp_resolution.do's t2_h* estimates when both
  run in the same Stata session (00_master.do runs 03 then 12) -- there is
  also no existing project precedent for `estimates save'/`estimates use'
  across do-files (grepped, none found). Every existing multi-panel table
  in this project re-estimates fresh within ONE do-file's own eststo/
  esttab block, so this file follows that same convention: it re-runs
  BOTH the GDP regression and the five channel regressions itself, in one
  file, so the eststo store is never cleared between them.

  SPEC, UPDATED to match 03_lp_resolution.do's CURRENT headline (previously
  this file used a JOINT regression + Wald F-test; that design was retired
  in 03_lp_resolution.do and this file had drifted out of sync with it --
  fixed here, not a new choice):
    TWO SEPARATE regressions per horizon, one per resolution type, EACH vs
    tranquil years with the RIVAL type dropped from that regression's own
    sample -- matching this project's AIPW design (08b_aipw.do's
    `_aipwpair') and 03_lp_resolution.do's own `_lpdiffboot'. Pooling the
    two arms in one joint regression lets the default-linked coefficient
    be estimated partly off non-default onset years (an implicit control,
    not an excluded observation) and forces a single residual variance
    across two populations this project treats as economically distinct --
    which is exactly why 03_lp_resolution.do dropped that design. See that
    file's own header for the full argument; unchanged here, only applied
    across all six outcomes instead of GDP alone.

    Because onset_nd and onset_def no longer share one regression, their
    difference has no directly estimable covariance from a single `test'
    command -- there is no joint VCE for a Wald F-test to use. The
    difference (default - non-default) and its SE/CI/p-value therefore
    come from the SAME paired row-bootstrap program 03_lp_resolution.do
    already defines, `_lpdiffboot' (copied here verbatim, not
    reimplemented, so both files share one tested mechanism): each of
    `nboot' replications resamples rows once (stratified so the tranquil
    pool and each arm's treated-row count are held fixed,
    `bsample, strata()') and refits BOTH arm regressions on that same
    resampled data, capturing the covariance the two arms share through
    their common tranquil control pool -- a naive sqrt(se_nd^2+se_def^2)
    would miss this. A Clogg et al. (1995) z is reported as a COMPANION
    statistic (independence-assuming, permissive); the bootstrap is the
    governing, adopted test throughout.

    xtreg <outcome> <dummy> <controls> if sample==1 & <rival>==0[&
    common_abcd==1], fe vce(robust) -- country FE only, no year FE, plain
    robust SE, tranquil years the omitted category, matching Asonuma et
    al.'s own Table I1 design, per arm. GDP + credit/inv/claims_govt run
    on the balanced A-D sample (common_abcd==1); FDI/real_lending keep
    their own best-available sample, unrestricted -- exactly as in both
    source files, unchanged by this fix.

  DISPLAY MECHANISM: esttab expects one eststo per column, each holding
  one regression's results -- it cannot natively show two coefficients
  that come from two separate regressions as if from one model.
  3_lp_resolution.do's own solution is reused verbatim: `ereturn post' a
  synthetic two-coefficient "model" (b_combo = (bnd, bdef), block-diagonal
  V_combo = diag(se_nd^2, se_def^2)) from each arm's own ORIGINAL
  (non-bootstrapped) point estimate and robust SE, purely so esttab's
  `keep(onset_nd onset_def)' can display both coefficients and their own
  SEs in one column, exactly as before. This posted V is NEVER used for
  inference -- the actual covariance the two arms share is what the
  paired bootstrap above already captures for the difference row.

  STATS ROW STANDARDIZED ACROSS ALL SIX PANELS to 03_lp_resolution.do's own
  Table 2 set: bdiff, sediff (bootstrap SE of the difference), lodiff/
  hidiff (95% bootstrap CI), pdiff (bootstrap p-value), cloggz/cloggp
  (Clogg z companion + its p-value), nepnd/nepdef (episodes per arm),
  nctynd/nctydef (countries per arm), nnd/ndef (observations per arm's own
  regression). The former bdiff/fdiff/pdiff/nepnd/nepdef/N/N_g set (Wald
  F-statistic, one shared N/N_g) no longer applies -- there is no joint
  regression to report a single N/N_g or Wald F from, and each arm now has
  its own N.

  Output: $tabs/table_combined_six_resolution.rtf.
  Run AFTER 18_transforms.do. Not wired into 00_master.do. Does not modify
  03_lp_resolution.do, 12_channels_resolution.do, or either of their own
  RTF outputs -- both keep producing their existing tables unchanged.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year
local controls $ctrl_core

* ── Channel outcomes ch_v_h = F h.v - L.v (h=0..4), same construction as
* 12_channels_resolution.do -- log real level for credit/inv, ratio for
* claims_govt/fdi/real_lending. Also builds pre_<v>, the pre-crisis-change
* control each channel's own ctrl_<v> local references below.
foreach v in credit inv claims_govt fdi real_lending {
    local src `v'
    if inlist("`v'","credit","inv") local src ln_r_`v'
    capture drop `v'_base
    gen double `v'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen double ch_`v'_`h' = F`h'.`src' - `v'_base
    }
    capture drop pre_`v'
    gen double pre_`v' = L.`src' - L2.`src'
}

* ── REPRODUCIBILITY: same seed/draw count as 03_lp_resolution.do's own
* `_lpdiffboot' calls, so this file's bootstrap draws are governed by the
* identical rationale (see that file's header) -- not re-derived here.
set seed 20260819
local nboot = 1000

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — _lpdiffboot, copied VERBATIM from 03_lp_resolution.do (see that
*   file's own header for the full rationale). Paired row-bootstrap of the
*   OLS def-nd DIFFERENCE, the OLS analog of 08b_aipw.do's `_aipwpair'.
*   _lpdiffboot , y() dnd() ifnd() ddef() ifdef() ctrlnd() ctrldef() reps()
* ══════════════════════════════════════════════════════════════════════════
capture program drop _lpdiffboot
program define _lpdiffboot, rclass
    syntax , Y(string) DND(string) IFND(string) DDEF(string) IFDEF(string) ///
             CTRLND(string) CTRLDEF(string) REPS(integer)

    capture xtreg `y' `dnd' `ctrlnd' if `ifnd', fe vce(robust)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bnd   = _b[`dnd']
    local send  = _se[`dnd']
    local nnd   = e(N)
    local ngnd  = e(N_g)

    capture xtreg `y' `ddef' `ctrldef' if `ifdef', fe vce(robust)
    if _rc {
        return scalar ok = 0
        exit
    }
    local bdef  = _b[`ddef']
    local sedef = _se[`ddef']
    local ndef  = e(N)
    local ngdef = e(N_g)

    local dh = `bdef' - `bnd'

    capture drop _pool
    quietly gen byte _pool = 0 if (`ifnd') | (`ifdef')
    quietly replace _pool = 1 if `dnd' == 1
    quietly replace _pool = 2 if `ddef' == 1

    tempname pf
    tempfile bf
    quietly postfile `pf' double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            quietly keep if !missing(_pool)
            quietly bsample, strata(_pool)
            capture regress `y' `dnd' `ctrlnd' i.cid if `ifnd', vce(robust)
            local t1 = cond(_rc==0, _b[`dnd'], .)
            capture regress `y' `ddef' `ctrldef' i.cid if `ifdef', vce(robust)
            local t2 = cond(_rc==0, _b[`ddef'], .)
            if !missing(`t1') & !missing(`t2') quietly post `pf' (`t2' - `t1')
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

    local pdiff = .
    if !missing(`se') & `se' > 0 local pdiff = 2*(1 - normal(abs(`dh'/`se')))

    local cloggz = .
    local cloggp = .
    if !missing(`send') & !missing(`sedef') & (`send'^2 + `sedef'^2) > 0 {
        local cloggz = `dh' / sqrt(`send'^2 + `sedef'^2)
        local cloggp = 2*(1 - normal(abs(`cloggz')))
    }

    return scalar ok    = 1
    return scalar bnd   = `bnd'
    return scalar send  = `send'
    return scalar nnd   = `nnd'
    return scalar ngnd  = `ngnd'
    return scalar bdef  = `bdef'
    return scalar sedef = `sedef'
    return scalar ndef  = `ndef'
    return scalar ngdef = `ngdef'
    return scalar dh    = `dh'
    return scalar se    = `se'
    return scalar lo    = `lo'
    return scalar hi    = `hi'
    return scalar nboot = `nd'
    return scalar p     = `pdiff'
    return scalar cloggz = `cloggz'
    return scalar cloggp = `cloggp'
end

* ── Helper: episode/country counts inside e(sample), matching
* 03_lp_resolution.do's own `_nepcount' exactly (copied verbatim), so all
* six panels report episodes AND countries per arm the same way.
capture program drop _nepcount
program define _nepcount, rclass
    syntax varname(numeric) , Outcome(varname) Controls(varlist)
    tempvar esmp
    quietly gen byte `esmp' = e(sample)
    quietly count if `esmp' == 1
    if r(N) == 0 {
        quietly replace `esmp' = (sample == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar n = r(N)
    tempvar tagcty
    quietly egen byte `tagcty' = tag(cid) if `varlist'==1 & `esmp'==1
    quietly count if `tagcty'==1
    return scalar ncty = r(N)
end

eststo clear

* ══════════════════════════════════════════════════════════════════════════
* SIX OUTCOMES, ONE LOOP: label | outcome-variable stem | own controls |
*   own eststo-name stem | balanced-sample flag. GDP kept as its own first
*   pass (own eststo names t2_h*, matching 03_lp_resolution.do's naming)
*   so the export block below can address it by name exactly as before;
*   channels loop with names t4_<ch>_h*, matching 12_channels_resolution.do.
* ══════════════════════════════════════════════════════════════════════════

local ctrl_credit       l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
local ctrl_inv          $ctrl_core pre_inv
local ctrl_claims_govt  $ctrl_core pre_claims_govt
local ctrl_fdi          $ctrl_core pre_fdi
local ctrl_real_lending $ctrl_core pre_real_lending

* ══════════════════════════════════════════════════════════════════════════
* PANEL A — GDP (mirrors 03_lp_resolution.do's own headline block)
* ══════════════════════════════════════════════════════════════════════════
di as result _n "=== PANEL A: GDP -- SEPARATE REGRESSIONS PER ARM + PAIRED BOOTSTRAP DIFFERENCE ==="
forvalues h = 0/4 {
    _lpdiffboot, y(dy_`h') ///
        dnd(onset_nd)  ifnd(sample==1 & onset_def==0 & common_abcd==1) ///
        ddef(onset_def) ifdef(sample==1 & onset_nd==0 & common_abcd==1) ///
        ctrlnd(`controls') ctrldef(`controls') reps(`nboot')

    if !r(ok) {
        di as error "  ** GDP h=" `h'+1 ": separate-regression estimate failed."
        continue
    }

    local bnd  = r(bnd)
    local bdef = r(bdef)
    local snd  = r(send)
    local sdef = r(sedef)
    local dh_    = r(dh)
    local se_    = r(se)
    local lo_    = r(lo)
    local hi_    = r(hi)
    local p_     = r(p)
    local cloggz_ = r(cloggz)
    local cloggp_ = r(cloggp)
    local nnd_   = r(nnd)
    local ndef_  = r(ndef)

    matrix b_combo = (`bnd', `bdef')
    matrix colnames b_combo = onset_nd onset_def
    matrix V_combo = diag((`snd'^2, `sdef'^2))
    matrix colnames V_combo = onset_nd onset_def
    matrix rownames V_combo = onset_nd onset_def
    ereturn post b_combo V_combo

    _nepcount onset_nd,  outcome(dy_`h') controls(`controls')
    local nepnd = r(n)
    local nctynd = r(ncty)
    _nepcount onset_def, outcome(dy_`h') controls(`controls')
    local nepdef = r(n)
    local nctydef = r(ncty)

    eststo t2_h`h', title("h=`=`h'+1'")
    estadd scalar bdiff   = `dh_'
    estadd scalar sediff  = `se_'
    estadd scalar lodiff  = `lo_'
    estadd scalar hidiff  = `hi_'
    estadd scalar pdiff   = `p_'
    estadd scalar cloggz  = `cloggz_'
    estadd scalar cloggp  = `cloggp_'
    estadd scalar nepnd   = `nepnd'
    estadd scalar nepdef  = `nepdef'
    estadd scalar nctynd  = `nctynd'
    estadd scalar nctydef = `nctydef'
    estadd scalar nnd     = `nnd_'
    estadd scalar ndef    = `ndef_'

    di "h=" `h'+1 ":  beta_nd=" %6.3f `bnd' "  beta_def=" %6.3f `bdef' ///
       "  diff(def-nd)=" %6.3f `dh_' "  boot SE=" %6.3f `se_' ///
       "  [" %6.3f `lo_' ", " %6.3f `hi_' "]  p=" %5.3f `p_' ///
       "  Clogg z=" %6.3f `cloggz_' " (p=" %5.3f `cloggp_' ")"
}

* ══════════════════════════════════════════════════════════════════════════
* PANELS B-F — channels (mirrors 12_channels_resolution.do's own loop,
*   same separate-regression + paired-bootstrap mechanism as GDP above)
* ══════════════════════════════════════════════════════════════════════════
local channels credit inv claims_govt fdi real_lending

foreach ch of local channels {
    local ctrl `ctrl_`ch''
    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

    di as result _n "=== PANEL: `ch' -- SEPARATE REGRESSIONS PER ARM + PAIRED BOOTSTRAP DIFFERENCE ==="
    forvalues h = 0/4 {
        _lpdiffboot, y(ch_`ch'_`h') ///
            dnd(onset_nd)  ifnd(sample==1 & onset_def==0`balflag') ///
            ddef(onset_def) ifdef(sample==1 & onset_nd==0`balflag') ///
            ctrlnd(`ctrl') ctrldef(`ctrl') reps(`nboot')

        if r(ok) {
            local bnd  = r(bnd)
            local bdef = r(bdef)
            local snd  = r(send)
            local sdef = r(sedef)
            local dh_    = r(dh)
            local se_    = r(se)
            local lo_    = r(lo)
            local hi_    = r(hi)
            local p_     = r(p)
            local cloggz_ = r(cloggz)
            local cloggp_ = r(cloggp)
            local nnd_   = r(nnd)
            local ndef_  = r(ndef)

            matrix b_combo = (`bnd', `bdef')
            matrix colnames b_combo = onset_nd onset_def
            matrix V_combo = diag((`snd'^2, `sdef'^2))
            matrix colnames V_combo = onset_nd onset_def
            matrix rownames V_combo = onset_nd onset_def
            ereturn post b_combo V_combo

            _nepcount onset_nd,  outcome(ch_`ch'_`h') controls(`ctrl')
            local nepnd = r(n)
            local nctynd = r(ncty)
            _nepcount onset_def, outcome(ch_`ch'_`h') controls(`ctrl')
            local nepdef = r(n)
            local nctydef = r(ncty)

            eststo t4_`ch'_`h', title("h=`=`h'+1'")
            estadd scalar bdiff   = `dh_'
            estadd scalar sediff  = `se_'
            estadd scalar lodiff  = `lo_'
            estadd scalar hidiff  = `hi_'
            estadd scalar pdiff   = `p_'
            estadd scalar cloggz  = `cloggz_'
            estadd scalar cloggp  = `cloggp_'
            estadd scalar nepnd   = `nepnd'
            estadd scalar nepdef  = `nepdef'
            estadd scalar nctynd  = `nctynd'
            estadd scalar nctydef = `nctydef'
            estadd scalar nnd     = `nnd_'
            estadd scalar ndef    = `ndef_'
            local elist_`ch' `elist_`ch'' t4_`ch'_`h'

            di "h=" `h'+1 ":  beta_nd=" %6.3f `bnd' "  beta_def=" %6.3f `bdef' ///
               "  diff(def-nd)=" %6.3f `dh_' "  boot SE=" %6.3f `se_' ///
               "  [" %6.3f `lo_' ", " %6.3f `hi_' "]  p=" %5.3f `p_' ///
               "  Clogg z=" %6.3f `cloggz_' " (p=" %5.3f `cloggp_' ")"
        }
        else di as error "  ** separate-regression estimate failed for `ch' h=" `h'+1
    }
}

* ══════════════════════════════════════════════════════════════════════════
* SINGLE MERGED ESTTAB — GDP panel first (replace), each channel appended
* ══════════════════════════════════════════════════════════════════════════
local statslabels labels("Difference (default - non-default)" "  Bootstrap SE (paired, def-nd)" ///
                          "  95% bootstrap CI, lower" "  95% bootstrap CI, upper" ///
                          "  p-value (paired row bootstrap)" ///
                          "  Clogg et al. (1995) z" "  p-value of Clogg z" ///
                          "Episodes (non-default)" "Episodes (default)" ///
                          "Countries (non-default arm)" "Countries (default arm)" ///
                          "Observations (non-default regression)" "Observations (default regression)")
local statsfmt fmt(3 3 3 3 3 3 3 0 0 0 0 0 0)

capture esttab t2_h0 t2_h1 t2_h2 t2_h3 t2_h4 using "$tabs/table_combined_six_resolution.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(onset_nd onset_def) order(onset_nd onset_def) ///
    coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
    mtitles nonumber ///
    stats(bdiff sediff lodiff hidiff pdiff cloggz cloggp nepnd nepdef nctynd nctydef nnd ndef, `statslabels' `statsfmt') ///
    title("Table X. Output cost by crisis resolution -- Panel A: GDP") ///
    addnotes("Dependent variable: cumulative change in the indicated outcome (pp) from t-1 to t+h. Non-default and default-linked onsets are each estimated in a SEPARATE regression vs tranquil years, with the rival resolution type dropped from that regression's own sample -- matching this project's AIPW design (08b_aipw.do) and 03_lp_resolution.do's own headline spec. Jorda (2005) local projections; country fixed effects only (no year FE); common-core controls plus (for channel panels) the channel's own pre-crisis change. Robust (heteroskedasticity-only) standard errors in parentheses, from each arm's own regression. Difference = beta(default) - beta(non-default); its standard error, 95% CI, and p-value come from a PAIRED ROW-BOOTSTRAP (1000 replications, seed 20260819) that refits both arm regressions on the same resampled rows each draw, preserving the covariance the two arms share through their common tranquil control pool. Clogg z is a companion statistic (permissive, assumes independence), not a replacement; the bootstrap governs. GDP, Bank credit, Investment, and Claims on government/GDP are estimated on the balanced A-D sample (common_abcd); FDI and Real lending rate use their own best-available sample. * p<0.10, ** p<0.05, *** p<0.01.")

local ptitle_credit       "Panel B: Bank credit"
local ptitle_inv          "Panel C: Investment"
local ptitle_claims_govt  "Panel D: Bank claims on government / GDP"
local ptitle_fdi          "Panel E: FDI"
local ptitle_real_lending "Panel F: Real lending rate"

foreach ch of local channels {
    if "`elist_`ch''" == "" {
        di as error "  ** table_combined_six_resolution.rtf: no estimates for `ch' -- panel skipped"
        continue
    }
    capture esttab `elist_`ch'' using "$tabs/table_combined_six_resolution.rtf", append ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_nd onset_def) order(onset_nd onset_def) ///
        coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
        mtitles nonumber ///
        stats(bdiff sediff lodiff hidiff pdiff cloggz cloggp nepnd nepdef nctynd nctydef nnd ndef, `statslabels' `statsfmt') ///
        title("`ptitle_`ch''")
}

if _rc == 608 di as error "  ** table_combined_six_resolution.rtf is OPEN IN WORD -- close it and re-run."
else if _rc di as error "  ** esttab failed on the last panel (rc=" _rc ")"
else di as result "Combined six-variable table saved: $tabs/table_combined_six_resolution.rtf"
