/*===========================================================================
  12_CHANNELS_RESOLUTION.DO
  Transmission Channels by Resolution Type (Non-Default vs. Default-Linked)

  Research question: do the transmission channels of spread crises differ
  between episodes resolved without default and those linked to default?

  STRATEGY:
  ---------
  For each of the 7 channels from 11_channels.do: TWO SEPARATE LPs per
  horizon, one per resolution type, EACH vs tranquil with the RIVAL type
  dropped from that regression's own sample:
    ch_var(h) = αi + β_nd(h)·onset_nd  + X·δ + ε      [if sample & onset_def==0]
    ch_var(h) = αi + β_def(h)·onset_def + X·δ + ε     [if sample & onset_nd==0]
  Country FE only, no year FE, plain robust SE -- the Stata-idiomatic
  equivalent of the reference paper's own reg ..., vce(robust) noconstant
  with explicit country dummies (see 02_lp_all.do's header).

  This REPLACES the previous headline, which entered onset_nd and onset_def
  in ONE joint regression on the full sample -- the same change made to
  03_lp_resolution.do's Table 2, for the same reason: matching this
  project's AIPW design (08b_aipw.do's `_aipwpair', each arm's own
  regression restricted to `if1'/`if2' excluding the rival type), and
  avoiding a joint regression in which one arm's own onset years contribute
  to identifying the other. Since the two arms no longer share one
  regression's VCE, the difference (def - nd) and its SE/CI/p-value now
  come from a PAIRED ROW-BOOTSTRAP (`_lpdiffboot' below, copied verbatim
  from 03_lp_resolution.do, mirroring how the AIPW files already share
  `_aipw'/`_aipwpair' program bodies across files): each of `nboot'
  replications resamples rows ONCE (stratified so the tranquil pool and
  each arm's treated-row count are held fixed) and refits BOTH arm
  regressions on that same resampled data, so the covariance the two arms
  inherit from their shared tranquil control pool is preserved in the
  bootstrap distribution of the difference -- unlike the retired Wald
  F-test's replacement candidate, treating the two arms as independent
  (sqrt(se_nd^2+se_def^2)), which would understate or misstate that
  covariance. Point estimates (b_nd, b_def) and their own robust SEs still
  come from the original, non-resampled fit of each arm's own regression;
  only the difference's uncertainty is bootstrapped.

  IPW REMOVED: this file used to carry a parallel IPW-weighted (Spec B)
  comparison, dropped project-wide once 08b_aipw.do's doubly-robust AIPW
  estimator superseded plain IPW as the estimator this project reports --
  see METHODOLOGY.md and 08b_aipw.do's header.

  BALANCED A-D SAMPLE (common_abcd, built in 18_transforms.do): the
  Investment, Bank credit, and Claims-on-govt columns are now ACTUALLY
  estimated on the common episode set shared with 03_lp_resolution.do's
  Table 2 and 08b_aipw.do's AIPW GDP result (Asonuma et al.'s own
  "balanced panels" convention -- see the per-channel loop below for the
  full argument). govexp, pb, fdi, and real_lending keep their own
  best-available sample, unrestricted.

  OUTPUTS:
  --------
  - fig12a_channels_ols.pdf   : 2×3 grid (nd=navy, def=brick)
  - channels_resolution.csv   : full coefficient table
===========================================================================*/

use "$clean/panel_lp.dta", clear
* safety: define the common core if this file is run standalone (master/18 also set it)
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"
sort cid year
xtset cid year

* ══════════════════════════════════════════════════════════════════════════
* 1. GENERATE CHANNEL OUTCOME VARIABLES
* ══════════════════════════════════════════════════════════════════════════

* OUTCOME SCALE. Strictly-positive GDP-ratio channels use the LOG REAL LEVEL
* (ln_r_*, built in 18_transforms) so the outcome is a cumulative percent change in
* the variable itself rather than in its ratio to a GDP that is collapsing — the
* reference paper's var2/var3 construction. pb and fdi change sign, so no log is
* possible and they keep the ratio; for a balance the ratio is the right object.
foreach var in credit claims_govt inv govexp pb fdi real_lending {
    local src `var'
    if inlist("`var'","credit","inv","govexp") local src ln_r_`var'
    capture drop `var'_base
    gen `var'_base = L.`src'
    forvalues h = 0/4 {
        capture drop ch_`var'_`h'
        gen ch_`var'_`h' = F`h'.`src' - `var'_base
    }
    * pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    * own-outcome pre-trend control, matching the AIPW spec in 13c.
    capture drop pre_`var'
    gen pre_`var' = L.`src' - L2.`src'
}

* ══════════════════════════════════════════════════════════════════════════
* 2. COVERAGE CHECK BY RESOLUTION TYPE
* ══════════════════════════════════════════════════════════════════════════

di as result _n "=== DATA COVERAGE AT ONSET BY RESOLUTION TYPE ==="
di as result "  Variable        nd (39)   def (22)"
foreach var in credit claims_govt inv govexp pb fdi real_lending {
    quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`var'_0)
    local n_nd = r(N)
    quietly count if onset_def == 1 & sample == 1 & !missing(ch_`var'_0)
    local n_def = r(N)
    di as result "  `var'" _col(20) `n_nd' " / 39" _col(32) `n_def' " / 22"
}

* ══════════════════════════════════════════════════════════════════════════
* 3. LP ESTIMATION BY CHANNEL
* ══════════════════════════════════════════════════════════════════════════

local channels   credit claims_govt inv govexp pb fdi real_lending

* Controls: common core ($ctrl_core) + each channel's own pre_<v> (same as 11_channels.do)
* Common-core controls (Asonuma-aligned $ctrl_core) + each channel's own pre_<v>;
* the core term measuring the channel's own lagged level is dropped from its own reg.
* For the credit channel that term is l_credit_bank: not literally the lag of the
* outcome series (credit = all financial corps, l_credit_bank = banks only), but the
* two correlate 0.950, so keeping it would be close to putting the lagged dependent
* variable on the RHS. It stays dropped.
local ctrl_credit      l1_gdpg l_debt l_banking_crisis l_govexp l_open l_lninfl exchange2 pre_credit
local ctrl_claims_govt $ctrl_core pre_claims_govt
local ctrl_inv         $ctrl_core pre_inv
local ctrl_govexp      l1_gdpg l_debt l_banking_crisis l_open l_credit_bank l_lninfl exchange2 pre_govexp
local ctrl_pb          $ctrl_core pre_pb
local ctrl_fdi         $ctrl_core pre_fdi
local ctrl_real_lending $ctrl_core pre_real_lending

* ── REPRODUCIBILITY: seed the bootstrap (matches 03_lp_resolution.do / 08b_aipw.do) ──
set seed 20260819
local nboot = 1000     // matches 08b_aipw.do's own G=1000

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — paired row-bootstrap of the def-nd DIFFERENCE, copied verbatim
*   from 03_lp_resolution.do's `_lpdiffboot' (see that file's header for the
*   full argument; mirrors how 08b_aipw.do/13c_aipw_channels.do already
*   share `_aipw'/`_aipwpair' program bodies across files). Point estimates
*   come from each arm's own original `xtreg ..., fe vce(robust)' fit; only
*   the difference's SE/CI/p come from the bootstrap.
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
end

* Initialize storage matrices
foreach ch of local channels {
    foreach grp in nd def {
        foreach m in b lo90 hi90 lo95 hi95 {
            matrix `m'_`grp'_`ch' = J(6, 1, 0)
        }
    }
    matrix pval_`ch' = J(6, 1, .)
}

* ── Loop over channels ───────────────────────────────────────────────────

eststo clear   // capture OLS resolution-split estimates for Table 4

foreach ch of local channels {

    local ctrl `ctrl_`ch''

    * BALANCED A-D SAMPLE (common_abcd, built in 18_transforms.do): Investment
    * and Bank credit (this file's own `inv'/`credit') plus Claims on
    * government (`claims_govt') are now ACTUALLY estimated on the same
    * common episode set as 03_lp_resolution.do's GDP and 08b_aipw.do's AIPW
    * GDP -- restricted to onsets where all four outcomes are non-missing at
    * every horizon h=0..4 (Asonuma et al.'s own "balanced panels"
    * convention). govexp, pb, fdi, and real_lending keep their own
    * best-available sample, unrestricted -- matching the reference paper's
    * own text that not every panel of a multi-panel result is balanced.
    local balflag
    if inlist("`ch'","credit","inv","claims_govt") local balflag " & common_abcd==1"

    di as result _n "========================================"
    di as result "CHANNEL: `ch'"
    di as result "========================================"
    di "h   b_nd     b_def    p(nd=def)"

    forvalues h = 0/4 {
        local row = `h' + 2

        * SEPARATE regressions per arm (rival type dropped from each), then
        * a paired row-bootstrap for the difference -- see this file's header
        * and `_lpdiffboot' above.
        _lpdiffboot, y(ch_`ch'_`h') ///
            dnd(onset_nd)  ifnd(sample==1 & onset_def==0`balflag') ///
            ddef(onset_def) ifdef(sample==1 & onset_nd==0`balflag') ///
            ctrlnd(`ctrl') ctrldef(`ctrl') reps(`nboot')

        if r(ok) {
            local bnd  = r(bnd)
            local bdef = r(bdef)
            local snd  = r(send)
            local sdef = r(sedef)

            matrix b_nd_`ch'[`row',1]    = `bnd'
            matrix lo90_nd_`ch'[`row',1] = `bnd'  - 1.645*`snd'
            matrix hi90_nd_`ch'[`row',1] = `bnd'  + 1.645*`snd'
            matrix lo95_nd_`ch'[`row',1] = `bnd'  - 1.960*`snd'
            matrix hi95_nd_`ch'[`row',1] = `bnd'  + 1.960*`snd'
            matrix b_def_`ch'[`row',1]   = `bdef'
            matrix lo90_def_`ch'[`row',1]= `bdef' - 1.645*`sdef'
            matrix hi90_def_`ch'[`row',1]= `bdef' + 1.645*`sdef'
            matrix lo95_def_`ch'[`row',1]= `bdef' - 1.960*`sdef'
            matrix hi95_def_`ch'[`row',1]= `bdef' + 1.960*`sdef'

            * Difference block: paired row-bootstrap SE/CI/p-value (NOT the
            * retired Wald F-test -- the two arms no longer share one
            * regression's VCE) + episode counts.
            local bdiff  = r(dh)
            local sediff = r(se)
            local lodiff = r(lo)
            local hidiff = r(hi)
            local pdiff_ = r(p)
            matrix pval_`ch'[`row',1] = `pdiff_'
            local pd_`ch'_`h' = `pdiff_'   // store before eststo (which can reset r())
            quietly count if onset_nd  == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepnd = r(N)
            quietly count if onset_def == 1 & sample == 1 & !missing(ch_`ch'_`h')
            local nepdef = r(N)
            local nnd_  = r(nnd)
            local ndef_ = r(ndef)

            * Synthesize a two-coefficient "model" via `ereturn post' from the
            * ORIGINAL (non-bootstrapped) per-arm point estimates and their
            * own robust SEs -- `_lpdiffboot' leaves e() holding the last
            * bootstrap-loop regression, not either arm's own fit. See
            * 03_lp_resolution.do's identical block for the full argument.
            matrix b_combo = (`bnd', `bdef')
            matrix colnames b_combo = onset_nd onset_def
            matrix V_combo = diag((`snd'^2, `sdef'^2))
            matrix colnames V_combo = onset_nd onset_def
            matrix rownames V_combo = onset_nd onset_def
            ereturn post b_combo V_combo

            eststo t4_`ch'_`h', title("h=`=`h'+1'")
            estadd scalar bdiff  = `bdiff'
            estadd scalar sediff = `sediff'
            estadd scalar lodiff = `lodiff'
            estadd scalar hidiff = `hidiff'
            estadd scalar pdiff  = `pdiff_'
            estadd scalar nepnd  = `nepnd'
            estadd scalar nepdef = `nepdef'
            estadd scalar nnd    = `nnd_'
            estadd scalar ndef   = `ndef_'
            local elist_`ch' `elist_`ch'' t4_`ch'_`h'
            local b_nd_o  = `bnd'
            local b_def_o = `bdef'
            local p_o     = `pd_`ch'_`h''
        }
        else {
            local b_nd_o  = .
            local b_def_o = .
            local p_o     = .
            di as error "regression failed for `ch' h=`=`h'+1'"
        }

        di "h=" `h'+1 "  " %7.3f `b_nd_o'  "  " %7.3f `b_def_o' "  " %5.3f `p_o'
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 3b. ROBUSTNESS: CREDIT CHANNEL EXCLUDING BULGARIA (LEAVE-ONE-OUT)
*
* Bulgaria's 1994 default-linked onset is followed by a real, well-documented
* banking-system collapse (WDI credit/GDP: 66% in 1993 -> 8.5% in 1997 -- a
* third of banks closed, hyperinflation, currency board adopted mid-1997),
* landing inside this LP's own h=3/h=4 outcome window. l_banking_crisis in
* the control set is LAGGED (t-1, predetermined) and cannot and should not
* net this out -- doing so would control away part of the very transmission
* channel (default -> banking distress -> credit collapse) this regression
* exists to measure, not a genuine confound. With only 13-14 default-linked
* episodes in the balanced A-D sample, though, one country's real crisis can
* set the def-arm average almost by itself -- this is a small-N
* generalizability question, not an omitted-variable one, checked directly
* here rather than argued informally. Diagnostic only: does NOT change the
* headline credit estimate above or the exported Table 4/IRF figures.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "========================================"
di as result "ROBUSTNESS: credit channel, Bulgaria excluded (leave-one-out)"
di as result "========================================"
di "h   b_nd(ex.Bulg)  b_def(ex.Bulg)  b_def(full sample)  p(nd=def, ex.Bulg)"
forvalues h = 0/4 {
    capture xtreg ch_credit_`h' onset_nd onset_def `ctrl_credit' ///
        if sample==1 & common_abcd==1 & country!="Bulgaria", fe vce(robust)
    if _rc == 0 {
        local bnd_lo1  = _b[onset_nd]
        local bdef_lo1 = _b[onset_def]
        test onset_nd = onset_def
        local p_lo1 = r(p)
        local bdef_full = b_def_credit[`h'+2,1]
        di "h=" `h'+1 "  " %9.3f `bnd_lo1' "  " %10.3f `bdef_lo1' ///
           "  " %14.3f `bdef_full' "  " %10.3f `p_lo1'
    }
    else di as error "  leave-one-out regression failed for credit h=" `h'+1 " (rc=" _rc ")"
}
di as result "  Read b_def(ex.Bulg) against b_def(full sample): if they stay close, the"
di as result "  headline credit coefficient is not just Bulgaria; if it collapses toward"
di as result "  zero, most of the def-arm signal was one country's crisis."

* ══════════════════════════════════════════════════════════════════════════
* 4. TABLE EXPORT — TABLE 4: Transmission channels by resolution type
*   Word/RTF, multi-panel: one panel per channel, columns = horizons h=0..4.
*   Each panel reports non-default and default-linked onset coefficients from the
*   JOINT regression (the reference paper's baseline), robust SE in parentheses,
*   plus the difference and the Wald F-test of their equality (test onset_nd =
*   onset_def) -- the covariance-correct answer for two coefficients estimated
*   in one regression, and the reference paper's own difference-test convention
*   (Table I1), not Clogg z or a bootstrap.
*   OLS spec (matches Tables 1-3). First panel replaces; the rest append.
*   Requires: ssc install estout
* ══════════════════════════════════════════════════════════════════════════

local t4note "Dependent variable: cumulative change in the channel variable (pp) from t-1 to t+h. Non-default and default-linked onsets are each estimated in a SEPARATE regression vs tranquil years, with the rival resolution type dropped from that regression's own sample (matching this project's AIPW design, 08b_aipw.do); the coefficient/SE shown for each arm come from that arm's own regression. Jorda (2005) local projections; country fixed effects only (no year FE); common-core controls plus the channel's own pre-crisis change; continuation years excluded. Robust (heteroskedasticity-only) standard errors in parentheses, from each arm's own regression. Difference = beta(default) - beta(non-default); its SE, 95% CI, and p-value come from a PAIRED ROW-BOOTSTRAP (1000 replications, seed 20260819): each replication resamples rows once, stratified so the tranquil pool and each arm's treated-row count are held fixed, and refits BOTH arm regressions on that same resampled data -- preserving the covariance the two arms inherit from their shared tranquil control pool, unlike treating them as independent. * p<0.10, ** p<0.05, *** p<0.01."

* Per-channel panel titles (Panel A carries the overall table caption)
local ptitle_credit      "Table 4. Channels by resolution (nd vs. def, joint) -- Panel A: Private credit/GDP"
local ptitle_claims_govt "Panel B: Bank claims on govt/GDP"
local ptitle_inv         "Panel C: Investment/GDP"
local ptitle_govexp      "Panel D: Govt expenditure/GDP"
local ptitle_pb          "Panel E: Primary balance/GDP"
local ptitle_fdi         "Panel F: FDI/GDP"
local ptitle_real_lending "Panel G: Real lending interest rate"

* Write one panel per channel to a single RTF. First successful panel uses
* "replace"; the rest "append". Each esttab is wrapped in capture so a locked
* file (open in Word) or a missing estimate warns and is skipped instead of
* halting the do-file.
local writemode replace
local t4fail 0

foreach ch in credit claims_govt inv govexp pb fdi real_lending {

    if "`elist_`ch''" == "" {
        di as error "  ** Table 4: no estimates for channel `ch' — panel skipped"
        local t4fail 1
        continue
    }

    local t4extra
    if "`ch'" == "real_lending" local t4extra addnotes("`t4note'")

    capture esttab `elist_`ch'' using "$tabs/table4_channels_resolution.rtf", `writemode' ///
        b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
        keep(onset_nd onset_def) order(onset_nd onset_def) ///
        coeflabel(onset_nd "Non-default onset" onset_def "Default-linked onset") ///
        mtitles nonumber ///
        stats(bdiff sediff lodiff hidiff pdiff nepnd nepdef nnd ndef, ///
              labels("Difference (default - non-default)" "  Bootstrap SE (paired, def-nd)" ///
                     "  95% bootstrap CI, lower" "  95% bootstrap CI, upper" ///
                     "  p-value (paired row bootstrap)" ///
                     "Episodes (non-default)" "Episodes (default)" ///
                     "Observations (non-default regression)" "Observations (default regression)") ///
              fmt(3 3 3 3 3 0 0 0 0)) ///
        title("`ptitle_`ch''") `t4extra'

    if _rc == 608 {
        di as error "  ** table4_channels_resolution.rtf is OPEN IN WORD — close it and re-run."
        local t4fail 1
        continue
    }
    else if _rc {
        di as error "  ** Table 4: esttab failed for panel `ch' (rc=" _rc ")"
        local t4fail 1
        continue
    }

    local writemode append
}

if `t4fail' == 0 di as result "Table 4 saved: $tabs/table4_channels_resolution.rtf"
else di as error "Table 4 written with warnings (see messages above)."

* ══════════════════════════════════════════════════════════════════════════
* 5. SAVE IRF DATASETS
* ══════════════════════════════════════════════════════════════════════════

foreach ch of local channels {
    foreach grp in nd def {
        preserve
            clear
            set obs 6
            gen horizon = _n - 1     // 0 (baseline), 1..5
            foreach m in b lo90 hi90 lo95 hi95 {
                svmat `m'_`grp'_`ch', names(`m')
                rename `m'1 `m'
            }
            gen channel = "`ch'"
            gen group   = "`grp'"
            save "$clean/irf_`grp'_`ch'.dta", replace
        restore
    }
}

* ══════════════════════════════════════════════════════════════════════════
* 6. EXPORT SUMMARY TABLE
* ══════════════════════════════════════════════════════════════════════════

preserve
    clear
    local nrows = 5 * 7   // 5 horizons × 7 channels
    set obs `nrows'
    gen channel  = ""
    gen horizon  = .
    foreach v in b_nd b_def p {
        gen `v' = .
    }

    local row = 1
    foreach ch of local channels {
        forvalues h = 0/4 {
            replace channel = "`ch'"              in `row'
            replace horizon = `h'                 in `row'
            replace b_nd    = b_nd_`ch'[`h'+2,1]  in `row'
            replace b_def   = b_def_`ch'[`h'+2,1] in `row'
            replace p       = pval_`ch'[`h'+2,1]  in `row'
            local ++row
        }
    }

    order channel horizon b_nd b_def p
    export delimited "$tabs/channels_resolution.csv", replace
    di as result "Table saved: $tabs/channels_resolution.csv"
restore

* ══════════════════════════════════════════════════════════════════════════
* 7. FIGURE — 2×3 MULTI-PANEL
* ══════════════════════════════════════════════════════════════════════════

* UNIFORM IRF STYLE (project-wide onset-tier convention): non-default =
* blue, default-linked = red, both solid, markers match line color, no
* legend (color already distinguishes the two lines).
local c_nd  "blue"
local c_def "red"

local titlelabels `" "Bank credit" "Bank claims on government" "Investment" "Government expenditure" "Primary balance" "FDI" "Real lending rate" "'

local i = 1
foreach ch of local channels {
    local tlab : word `i' of `titlelabels'

    * Y-axis title shown only on the leftmost panel of each row (cols(4)
    * rows(2): panels 1 and 5), matching the reference paper's own Figure 2
    * -- not repeated on every panel.
    local ytit ""
    if inlist(`i', 1, 5) local ytit "Cumulative percent change"

    use "$clean/irf_nd_`ch'.dta",  clear
    append using "$clean/irf_def_`ch'.dta"

    twoway ///
        (rarea lo90 hi90 horizon if group=="nd", ///
            color("`c_nd'%20") lwidth(none)) ///
        (rarea lo90 hi90 horizon if group=="def", ///
            color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if group=="nd", ///
            lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) ///
            lwidth(medthick) msize(small)) ///
        (connected b horizon if group=="def", ///
            lcolor("`c_def'") mcolor("`c_def'") msymbol(square) ///
            lwidth(medthick) msize(small)) ///
        , ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)5, labsize(medium)) ///
        ylabel(, format(%9.0f) labsize(medium) angle(horizontal)) ///
        xtitle("Year", size(medium)) ///
        ytitle("`ytit'", size(medsmall)) ///
        title("`tlab'", size(medlarge) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(ols_`i', replace)

    local ++i
}

* No note() and no overall combine title() -- paper-ready. Each panel's own
* {bf:...} title already names the channel, so a combine-level title would
* only repeat what the panel titles already say once several are merged
* into one figure; dropped rather than kept as redundant text.
graph combine ols_1 ols_2 ols_3 ols_4 ols_5 ols_6 ols_7, ///
    cols(4) rows(2) ///
    graphregion(color(white)) xsize(12) ysize(7)

graph export "$figs/fig12a_channels_ols.pdf", replace
di as result "Figure saved: fig12a_channels_ols.pdf"
forvalues i = 1/7 {
    capture graph drop ols_`i'
}

* ══════════════════════════════════════════════════════════════════════════
* 8. COMBINED 6-PANEL FIGURE: Panel A GDP, B Investment, C Bank credit,
*    D Claims on government, E FDI, F Real lending rate. GDP's own IRF is
*    NOT built in this file (it is Table 2's own headline result) -- read
*    from 03_lp_resolution.do's saved irf_nd.dta/irf_def.dta (horizon -1..5,
*    columns b/se/lo90/hi90/lo95/hi95/series; 03 runs before this file in
*    00_master.do, so those files already exist). Kept ALONGSIDE the 7-panel
*    figure above (Section 7), not a replacement for it.
* ══════════════════════════════════════════════════════════════════════════
local combo_vars   gdp inv credit claims_govt fdi real_lending
local combo_labels `" "Panel A: GDP" "Panel B: Investment" "Panel C: Bank credit" "Panel D: Claims on govt" "Panel E: FDI" "Panel F: Real lending rate" "'
local i = 1
foreach cv of local combo_vars {
    local clab : word `i' of `combo_labels'
    local ytit ""
    if inlist(`i', 1, 4) local ytit "Cumulative percent change"

    if "`cv'" == "gdp" {
        * rename series->group in BOTH files before appending -- append
        * matches columns by NAME, so renaming only the in-memory dataset
        * would leave irf_def.dta's own "series" column unmatched (a
        * separate, mostly-empty column) rather than merged into "group".
        use "$clean/irf_nd.dta", clear
        rename series group
        tempfile _gdpnd
        save `_gdpnd'
        use "$clean/irf_def.dta", clear
        rename series group
        append using `_gdpnd'
        keep if horizon >= 0
    }
    else {
        use "$clean/irf_nd_`cv'.dta",  clear
        append using "$clean/irf_def_`cv'.dta"
    }

    twoway ///
        (rarea lo90 hi90 horizon if group=="nd", ///
            color("`c_nd'%20") lwidth(none)) ///
        (rarea lo90 hi90 horizon if group=="def", ///
            color("`c_def'%20") lwidth(none)) ///
        (connected b horizon if group=="nd", ///
            lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) ///
            lwidth(medthick) msize(small)) ///
        (connected b horizon if group=="def", ///
            lcolor("`c_def'") mcolor("`c_def'") msymbol(square) ///
            lwidth(medthick) msize(small)) ///
        , ///
        yline(0, lcolor(gs10) lpattern(dash) lwidth(thin)) ///
        xlabel(0(1)5, labsize(large)) ///
        ylabel(, format(%9.0f) labsize(large) angle(horizontal)) ///
        xtitle("Year", size(large)) ///
        ytitle("`ytit'", size(large)) ///
        title("`clab'", size(medlarge) color(navy)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(comb_`i', replace)

    local ++i
}
graph combine comb_1 comb_2 comb_3 comb_4 comb_5 comb_6, ///
    cols(3) rows(2) graphregion(color(white)) xsize(10) ysize(7)
graph export "$figs/fig12_combined_ols.pdf", replace
di as result "Figure saved: fig12_combined_ols.pdf (Panel A-F: GDP, Investment, Bank credit, Claims on government, FDI, Real lending rate)"
forvalues i = 1/6 {
    capture graph drop comb_`i'
}

di as result _n "12_channels_resolution.do complete."
