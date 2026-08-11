/*===========================================================================
  19_SAMPLE_AUDIT.DO   —  DIAGNOSTIC ONLY, ESTIMATES NOTHING

  Why the headline regression uses 38 of 61 onsets, control by control.

  61 onsets -> 59 reach `sample' (2 lack ln_gdp_base) -> 38 survive listwise
  deletion on $ctrl_core. This file accounts for the missing ~21 so the decision
  about which controls to keep is made on evidence rather than on marginal
  coverage counts, which OVERSTATE the fixable part: the same country-year is
  usually missing several series at once, so dropping one control recovers far
  fewer onsets than its own coverage gap suggests.

  Sections:
    1. Coverage of each core control at onsets (all / nd / def)
    2. Leave-one-out: onsets each control UNIQUELY kills (the honest marginal cost)
    3. Candidate control sets: resulting onset counts, side by side
    4. l_credit vs l_credit_bank — is the swap free?
    5. Interpolation potential: how much missingness is interior gaps
    6. The dropped onsets, named, with the controls responsible

  Read-only: loads panel_lp.dta, writes one CSV, changes nothing on disk that
  the pipeline consumes. Not wired into 00_master.do — run it on its own.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"

local core   $ctrl_core
local outref dy_0        // reference horizon: h=0 has the largest sample

* Helper: onsets of a given type surviving `sample' + non-missing outcome + a control list
capture program drop _survive
program define _survive, rclass
    * ctrls optional so the no-controls ceiling spec can be priced too
    syntax [anything(name=ctrls id="controls")] , Outcome(varname) Treat(varname)
    tempvar t
    quietly gen byte `t' = (sample == 1)
    quietly markout `t' `outcome' `ctrls'
    quietly count if `treat' == 1 & `t' == 1
    return scalar n = r(N)
end

* ══════════════════════════════════════════════════════════════════════════
* 0. THE ARITHMETIC
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "19_SAMPLE_AUDIT — where the onsets go"
di as result "{hline 78}"

foreach t in onset_all onset_nd onset_def {
    quietly count if `t' == 1
    local n_tot = r(N)
    quietly count if `t' == 1 & sample == 1
    local n_smp = r(N)
    quietly count if `t' == 1 & sample == 1 & !missing(`outref')
    local n_out = r(N)
    _survive `core', outcome(`outref') treat(`t')
    local n_reg = r(n)
    di as result _n "  `t':"
    di as result "     in the crisis database            : `n_tot'"
    di as result "     reach sample==1                   : `n_smp'"
    di as result "     ... with a non-missing outcome    : `n_out'"
    di as result "     ... surviving the full core       : `n_reg'   <- what identifies the coefficient"
}

* ══════════════════════════════════════════════════════════════════════════
* 1. COVERAGE OF EACH CORE CONTROL AT ONSETS
*    Raw coverage. Read this as an upper bound on what dropping a control can
*    recover, never as the actual gain — see section 2 for that.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "1. RAW COVERAGE AT ONSETS (upper bound on the gain from dropping)"
di as result "{hline 78}"
di as result "  control" _col(26) "all/61" _col(38) "nd/40" _col(48) "def/21"

foreach v of local core {
    quietly count if onset_all==1 & !missing(`v')
    local ca = r(N)
    quietly count if onset_nd==1  & !missing(`v')
    local cn = r(N)
    quietly count if onset_def==1 & !missing(`v')
    local cd = r(N)
    di as result "  `v'" _col(26) %3.0f `ca' _col(38) %3.0f `cn' _col(48) %3.0f `cd'
}

* ══════════════════════════════════════════════════════════════════════════
* 2. LEAVE-ONE-OUT: WHAT EACH CONTROL UNIQUELY COSTS
*    An onset "uniquely killed" by v passes every OTHER core control and the
*    outcome, and fails only on v. That count IS the number of onsets you would
*    recover by dropping v — no more, no less. Controls with a big raw coverage
*    gap but a small unique cost are missing on the SAME country-years as
*    something else, so dropping them buys nothing.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "2. MARGINAL COST — onsets recovered by dropping this control ALONE"
di as result "{hline 78}"
di as result "  control" _col(26) "all" _col(34) "nd" _col(42) "def" _col(52) "(raw gap, all)"

tempname AUD
postfile `AUD' str24 control long n_all long n_nd long n_def ///
    using "$tabs/sample_audit_loo.dta", replace

foreach v of local core {
    local others : list core - v

    tempvar tk
    quietly gen byte `tk' = (sample == 1)
    quietly markout `tk' `outref' `others'

    * passes everything else, fails only on v
    quietly count if onset_all==1 & `tk'==1 & missing(`v')
    local ua = r(N)
    quietly count if onset_nd==1  & `tk'==1 & missing(`v')
    local un = r(N)
    quietly count if onset_def==1 & `tk'==1 & missing(`v')
    local ud = r(N)

    quietly count if onset_all==1 & missing(`v')
    local rawgap = r(N)

    di as result "  `v'" _col(26) %3.0f `ua' _col(34) %3.0f `un' _col(42) %3.0f `ud' ///
                 _col(52) "(" %2.0f `rawgap' ")"
    post `AUD' ("`v'") (`ua') (`un') (`ud')
    drop `tk'
}
postclose `AUD'

di as result _n "  A control whose marginal cost is far below its raw gap is missing on the"
di as result "  same country-years as another control — dropping it recovers nothing."

* ══════════════════════════════════════════════════════════════════════════
* 3. CANDIDATE CONTROL SETS
*    Each row is a complete alternative to $ctrl_core. The question is how many
*    onsets it buys and what it gives up — the trade-off is yours to make, this
*    only prices it.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "3. CANDIDATE CORES — onsets surviving, at h=0 and h=4"
di as result "{hline 78}"
di as result "  spec" _col(34) "h=0 all/nd/def" _col(56) "h=4 all/nd/def"

* Each entry: "label | varlist"
local spec1 "full core (current)          | $ctrl_core"
local spec2 "banking as a dummy           | l1_gdpg l_debt l_ca l_banking_crisis l_govexp l_open l_credit_bank l_hyperinfl"
local spec3 "drop the depth term          | l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_hyperinfl"
local spec4 "drop l_debt l_ca (Asonuma)   | l1_gdpg l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
local spec5 "Asonuma + l_credit           | l1_gdpg l_banking_duration l_govexp l_open l_credit l_hyperinfl"
local spec6 "lean (>=58/61 controls only) | l1_gdpg l_banking_duration l_hyperinfl"
local spec7 "no controls (ceiling)        | "

forvalues s = 1/7 {
    * split "label | varlist" on the pipe
    local full  `"`spec`s''"'
    local pos   = strpos(`"`full'"', "|")
    local lbl   = strtrim(substr(`"`full'"', 1, `pos'-1))
    local vars  = strtrim(substr(`"`full'"', `pos'+1, .))

    _survive `vars', outcome(dy_0) treat(onset_all)
    local a0 = r(n)
    _survive `vars', outcome(dy_0) treat(onset_nd)
    local n0 = r(n)
    _survive `vars', outcome(dy_0) treat(onset_def)
    local d0 = r(n)

    _survive `vars', outcome(dy_4) treat(onset_all)
    local a4 = r(n)
    _survive `vars', outcome(dy_4) treat(onset_nd)
    local n4 = r(n)
    _survive `vars', outcome(dy_4) treat(onset_def)
    local d4 = r(n)

    di as result "  `lbl'" _col(34) %2.0f `a0' " / " %2.0f `n0' " / " %2.0f `d0' ///
                 _col(56) %2.0f `a4' " / " %2.0f `n4' " / " %2.0f `d4'
}

di as result _n "  The def column is the one that matters: the default premium is the"
di as result "  paper's headline and it is identified off the smallest cell."

* ══════════════════════════════════════════════════════════════════════════
* 4. l_credit VS l_credit_bank
*    Both measure private-credit depth; `credit' is all financial corporations
*    (FS.AST.PRVT.GD.ZS), `credit_bank' is banks only (FD.AST.PRVT.GD.ZS). In EM
*    panels banks are nearly all of private credit, so if credit_bank is better
*    covered the swap is close to free.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "4. l_credit VS l_credit_bank"
di as result "{hline 78}"

foreach v in l_credit l_credit_bank {
    capture confirm variable `v'
    if _rc {
        di as error "  `v' not in panel_lp.dta — skipped"
        continue
    }
    quietly count if !missing(`v')
    local nrow = r(N)
    quietly count if onset_all==1 & !missing(`v')
    local ca = r(N)
    quietly count if onset_nd==1  & !missing(`v')
    local cn = r(N)
    quietly count if onset_def==1 & !missing(`v')
    local cd = r(N)
    di as result "  `v'" _col(24) "panel rows: " %5.0f `nrow' ///
                 "   onsets all/nd/def: " %2.0f `ca' " / " %2.0f `cn' " / " %2.0f `cd'
}

* Do they actually agree where both exist? A high correlation means the swap is
* a coverage gain, not a change of concept.
capture confirm variable l_credit_bank
if _rc == 0 {
    quietly count if !missing(l_credit) & !missing(l_credit_bank)
    di as result "  overlap (both non-missing): " r(N) " rows"
    quietly correlate l_credit l_credit_bank
    di as result "  correlation where both observed: " %5.3f r(rho)
    di as result "     (>0.95 => same object, so the swap is a coverage gain not a spec change)"
}

* ══════════════════════════════════════════════════════════════════════════
* 5. INTERPOLATION POTENTIAL
*    An INTERIOR gap (missing years bracketed by observed years within the same
*    country) can be linearly interpolated: these ratios are slow-moving stocks
*    and interpolation invents no variation, it just spans a reporting hole.
*    Missingness at the ENDS cannot — that would be extrapolation, and it is
*    exactly where crises cluster, so it must be left alone.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "5. INTERPOLATION POTENTIAL (interior gaps only)"
di as result "{hline 78}"

foreach v in credit credit_bank debt open govexp ca {
    capture confirm variable `v'
    if _rc continue

    tempvar fy ly
    quietly bysort cid: egen `fy' = min(cond(!missing(`v'), year, .))
    quietly bysort cid: egen `ly' = max(cond(!missing(`v'), year, .))

    quietly count if missing(`v') & sample==1 & inrange(year, `fy', `ly')
    local interior = r(N)
    quietly count if missing(`v') & sample==1 & !inrange(year, `fy', `ly')
    local edge = r(N)
    quietly count if onset_all==1 & missing(`v') & inrange(year, `fy', `ly')
    local ionset = r(N)

    di as result "  `v'" _col(20) "interior gaps: " %4.0f `interior' ///
                 "   at the ends: " %4.0f `edge' ///
                 "   onsets in interior gaps: " %2.0f `ionset'
    drop `fy' `ly'
}

di as result _n "  'onsets in interior gaps' is what interpolation would recover for that"
di as result "  series. If it is ~0 the missingness is structural (country never"
di as result "  reported, or stopped) and interpolation is not the answer."

* ══════════════════════════════════════════════════════════════════════════
* 6. THE DROPPED ONSETS, NAMED
*    Which episodes are actually being lost, and to what. Eyeball this: if the
*    dropped set is dominated by the recent, severe episodes (Lebanon, Sri Lanka,
*    Ghana, Zambia) then the surviving sample is not just smaller, it is
*    SELECTED, and that is a bigger problem than the loss of power.
* ══════════════════════════════════════════════════════════════════════════
di as result _n "{hline 78}"
di as result "6. ONSETS DROPPED BY THE FULL CORE (in sample, outcome present)"
di as result "{hline 78}"

tempvar tfull
quietly gen byte `tfull' = (sample == 1)
quietly markout `tfull' `outref' `core'

* count of core controls missing, per row, to show what is responsible
tempvar nmiss
quietly gen byte `nmiss' = 0
foreach v of local core {
    quietly replace `nmiss' = `nmiss' + missing(`v')
}

quietly count if onset_all==1 & sample==1 & !missing(`outref') & `tfull'==0
di as result "  dropped: " r(N) " onsets" _n

capture noisily list country year onset_nd onset_def `nmiss' `core' ///
    if onset_all==1 & sample==1 & !missing(`outref') & `tfull'==0, ///
    noobs sepby(country) abbrev(12)

* Export for offline inspection
preserve
    quietly keep if onset_all==1 & sample==1 & !missing(`outref') & `tfull'==0
    quietly keep country iso3 year onset_nd onset_def `core'
    quietly export delimited "$tabs/sample_audit_dropped_onsets.csv", replace
    di as result _n "  Dropped-onset list saved: $tabs/sample_audit_dropped_onsets.csv"
restore

di as result _n "{hline 78}"
di as result "19_sample_audit.do complete — nothing estimated, nothing changed."
di as result "{hline 78}"
