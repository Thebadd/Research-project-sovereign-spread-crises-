/*===========================================================================
  26_LP_DEBTCRISIS_FLOW.DO
  Output cost of a BROADENED sovereign debt crisis taxonomy — flow-coded,
  three treatment arms (non-default, preemptive default, post-default) —
  built on 17b_merge_at_full.do's AT-side columns.

  WHY THIS FILE EXISTS
  ---------------------
  Every existing flow-tier file in this project (20-25) uses the project's
  ORIGINAL treatment: a spread crisis (EMBIG-based), classified default-
  linked only when an AT restructuring coincides with THAT episode. This
  file answers a different, broader question: the cost of a sovereign DEBT
  crisis, where default-linked status can also come directly from the full
  AT database (17b) when a restructuring happened without producing a
  detected spread crisis. Non-default crises are still identified ONLY via
  the spread database, unchanged — there is no restructuring-record
  equivalent for "elevated borrowing cost, no default" (this asymmetry is
  deliberate: a default is a discrete, documented event AT can identify
  independently of spreads; a non-default crisis has no such independent
  marker, its only observable signature IS the spread itself).

  THE COMBINED TREATMENT IS BUILT HERE, SELF-CONTAINED — NOT IN
  18_TRANSFORMS.DO
  ---------------------------------------------------------------------------
  Same discipline as 25_aipw_nexus_split_flow.do building its own a_nexus
  rather than modifying 18_transforms.do: this broadened taxonomy is a
  design choice specific to this file's research question, not a change to
  the project's core treatment definition that every other file (20-25)
  should inherit. Those files are untouched by this one.

  BUILDING dc_in_crisis_* FROM in_crisis_* (spread) + at_default_year (AT)
  ---------------------------------------------------------------------------
  dc_in_crisis_nd    = in_crisis_nd, UNCHANGED (spread-only by design).
  dc_in_crisis       = in_crisis_nd | dc_default_year (below).
  dc_default_year    = in_crisis_def==1  (spread-triggered default-linked
                        year, the existing classification)
                        OR at_default_year==1 & in_crisis==0 (an AT-only
                        default-linked year: no spread crisis was detected,
                        but AT records a restructuring in progress — a
                        genuinely NEW treated year, entirely additive to the
                        existing 234 flow-treated rows)
  CONFLICT (flagged, not silently resolved): rows where the SPREAD
  classification says non-default (in_crisis_nd==1) but AT independently
  records a default/restructuring window covering that same row
  (at_default_year==1). This file's design explicitly trusts AT for default
  identification (the entire reason this file exists), so these rows are
  RECLASSIFIED to default-linked -- but every such row is printed in
  Section 1's diagnostic block, not silently reassigned.

  2-YEAR RESTRUCTURING-GAP RULE (motivated by Zambia: spread crisis onset
  2015, actual default Nov 2020 -- a 5-year gap during which the ORIGINAL
  in_crisis_def classification, inherited whole-episode from
  10_skeleton.do/18_transforms.do's onset-level nondefault flag, marks every
  year 2015-2024 default-linked even though no default had happened until
  2020). For any spread-based default episode that DOES have at least one
  row overlapping an AT window, at_episode_onset (AT's own restructuring-
  start date) gives an independent check: if it starts more than 2 years
  after the spread episode's own onset, the years BEFORE the AT onset are
  reclassified non-default here (the market-stress episode existed, but no
  default had actually begun), and only the years from the AT onset onward
  keep default-linked status and AT type. LIMITATION, stated plainly: this
  only fires where SOME AT overlap exists. Zambia, Ghana, Lebanon, Sri Lanka
  and Ukraine currently have ZERO AT overlap at all (AT's vintage cutoff is
  ~2020, before any of their actual defaults were recorded there) -- they
  are "Unclassified" for that reason, not a long-gap reason, so this rule
  does NOT split them; that requires their events to be added to
  17b_merge_at_full.do's at_supplement block first. Scope, confirmed with
  the user: this rule is applied ONLY here, not to the project's original
  in_crisis_def/nondefault classification used by every onset/flow file
  (02-13, 20-25).

  PREEMPTIVE vs POST-DEFAULT TYPE
  ---------------------------------------------------------------------------
  Comes directly from at_type where a row's own year falls inside an AT
  window. For a spread-triggered default-linked EPISODE where the onset (or
  some continuation) row does overlap an AT window but not every row does,
  the type is filled across the whole episode (same forward-fill idiom
  18_transforms.do already uses for nd_ep, `bysort cid ep_seq: egen ... max
  (...)'). Episodes with NO overlapping AT window anywhere (the existing
  Episode_Summary-based default classification was right that a
  restructuring happened, but the FULL AT database's own dating does not
  reach that row) are left "Unclassified (no AT type match)" and reported
  as a coverage diagnostic, not guessed.

  SPECIFICATION — MIRRORS 20_lp_flow.do's SECTION 3 EXACTLY
  ---------------------------------------------------------------------------
  dy_h on dc_in_crisis_nd, dc_in_crisis_preempt, dc_in_crisis_post (jointly,
  tranquil years omitted) + $ctrl_flow (episode-dated core, unchanged --
  this file does not redefine the control set) + country & year FE.
  Driscoll-Kraay SE at lag(max(2,h+3)), the established flow-tier rule
  (the regressor is itself flow-coded and serially correlated within an
  episode). Pairwise differences (post-vs-preempt, post-vs-nd, preempt-vs-nd)
  via `lincom' -- covariance-correct within the joint regression, the SAME
  reasoning 20/23 already use for their own multi-arm differences (this is
  plain OLS, no propensity model, so none of the AIPW files' bootstrap/
  Clogg-z machinery applies here).

  Same flow_ctrl_variant / drop_year_fe toggles as every other flow file,
  for consistency, though this file does not build its own control-set
  variant — it reads $ctrl_flow / $ctrl_flow_flowplus exactly as built in
  18_transforms.do.

  Outputs
  -------
    "$tabs/table_debtcrisis_flow.rtf"        three columns: nd / preemptive /
                                              post, esttab's standard layout
    "$tabs/table_debtcrisis_flow_layout.rtf" SAME regression, laid out to match
                                              the reference paper's Table 2:
                                              coefficient(stars)/SE/Observations-
                                              Countries-Episodes stacked under
                                              each arm, columns h=1..h=5. Hand-
                                              built (file write) since esttab's
                                              stats() rows can't interleave a
                                              different text row under each
                                              coefficient -- only one shared row
                                              per column at the table's bottom.
    "$tabs/debtcrisis_flow.csv"              raw coefficients + pairwise diffs +
                                              PER-ARM observations/episodes/
                                              countries at every horizon
                                              (nrow/nep/ncid columns)
    "$figs/fig_debtcrisis_flow.pdf"          three-line IRF, nd/preempt/post overlay

  PER-ARM N/EPISODES/COUNTRIES, AND A BUG FIXED ALONG THE WAY
  ---------------------------------------------------------------------------
  esttab's own N/N_g stats (in the RTF) are for the POOLED regression (all
  three arms + tranquil controls), not broken out by arm -- not useful for
  judging how much any one arm's coefficient rests on. _nflowcount computes
  the per-arm observation count, episode count and (now) country count
  directly from e(sample), and the estimation loop posts all three to
  debtcrisis_flow.csv and prints them to console after Section 3's table.
  Fixed while adding this: the preempt/post arms were being episode-tagged
  on `ep_seq` (the original spread-episode counter), under which every
  AT-only country's default rows share ep_seq==0 (onset_all is structurally
  0 for them, so ep_seq never increments) -- collapsing every one of that
  country's AT-only episodes into a single tagged row and silently
  undercounting. `_nflowcount` now takes an `epvar()` argument; the nd call
  passes ep_seq (correct for it, since non-default is always spread-based),
  the preempt/post calls pass dc_ep_seq (this file's own combined counter,
  which increments correctly for both spread-based and AT-only episodes).

  SELF-CONTAINED: reads only $clean/panel_lp.dta (which already carries
  17b's at_* columns, merged before 18_transforms.do ran and therefore
  carried through unchanged into panel_lp.dta).

  COUNTRY COVERAGE — THE ONE FILE WHERE IT DIFFERS FROM THE REST OF THE
  PROJECT
  ---------------------------------------------------------------------------
  10b_skeleton_atonly.do (between 10_skeleton.do and 11_weo.do) extends the
  panel with country-year rows for every AT-covered country outside this
  project's 52-country spread panel, flagged `atonly_country==1`. Per
  Pescatori & Sy (2007) and METHODOLOGY.md §5.3, AT's own default record is
  a primary authority that does not need spread confirmation to count, so
  these countries' default-linked years are legitimate additions here --
  but they can NEVER populate dc_in_crisis_nd (no spread data exists to
  establish a non-default finding for them). Every other file in the
  project (02-25) reads $clean/panel_lp.dta's `sample'/`sample_flow' flags,
  which EXCLUDE atonly_country==1 rows (18_transforms.do), so their
  published samples are completely unaffected by this expansion. THIS file
  alone builds its own broadened flag, `sample_flow_bc' (Section 0 below),
  that includes them -- deliberately, since giving these countries a real
  surrounding window of untreated years (not just their crisis years) is
  what lets them contribute genuine within-country identifying variation
  under the country+year FE design, rather than being fully absorbed by
  their own fixed effect.
===========================================================================*/

use "$clean/panel_lp.dta", clear
if "$ctrl_core"=="" global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"
sort cid year
xtset cid year

* BROADENED sample flag: unlike every other flow file, this one deliberately
* INCLUDES 10b_skeleton_atonly.do's AT-only countries (outside the 52-country
* spread panel) -- they can only ever appear in dc_in_crisis_preempt/_post
* (dc_in_crisis_nd is structurally 0 for them; see METHODOLOGY.md §5.3), but
* their surrounding untreated years are real, needed control observations for
* this file's own country+year FE regression. sample_flow itself (built in
* 18_transforms.do) EXCLUDES them, matching every other onset/flow file's
* existing sample -- this file reads panel_lp.dta's own atonly_country flag
* directly rather than sample_flow, so this is the ONLY file in the project
* where the broadened taxonomy's country coverage actually differs from the
* base panel's.
capture confirm variable atonly_country
if _rc gen byte atonly_country = 0
gen byte sample_flow_bc = !missing(ln_gdp_base) & carryin==0
label var sample_flow_bc "Broadened flow sample: sample_flow's own rule, WITHOUT the atonly_country exclusion"

foreach v in in_crisis in_crisis_nd in_crisis_def onset_all continuation ///
             ep_seq sample_flow at_default_year at_type {
    capture confirm variable `v', exact
    if _rc {
        di as error "  ** `v' not in panel_lp.dta — re-run 18_transforms.do (and 17b_merge_at_full.do" ///
                     " before it) first."
        exit 111
    }
}
if "$ctrl_flow" == "" {
    di as error "  ** \$ctrl_flow not set — re-run 18_transforms.do, then this file."
    exit 111
}

* ══════════════════════════════════════════════════════════════════════════
* HELPERS — identical to 20_lp_flow.do's
* ══════════════════════════════════════════════════════════════════════════
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
    syntax varname(numeric) , Outcome(varname) Controls(varlist) Samp(varname) [Epvar(varname)]
    * Epvar defaults to ep_seq (the original spread-episode counter), correct
    * for dc_in_crisis_nd since non-default rows are always spread-based. For
    * dc_in_crisis_preempt/_post the caller MUST pass epvar(dc_ep_seq) instead:
    * an AT-only country's default rows all carry ep_seq==0 (onset_all is
    * structurally 0 for them), so tag(cid ep_seq) would silently collapse
    * every one of that country's AT-only episodes into a single tagged row --
    * confirmed the hard way (this bug undercounted preempt/post episode
    * counts before dc_ep_seq, the file's own combined counter, was wired in
    * here as the fix).
    tempvar esmp tagep tagcid
    quietly gen byte `esmp' = e(sample)
    quietly count if `esmp' == 1
    if r(N) == 0 {
        quietly replace `esmp' = (`samp' == 1)
        quietly markout `esmp' `outcome' `controls'
    }
    quietly count if `varlist' == 1 & `esmp' == 1
    return scalar nrow = r(N)
    local epv = cond("`epvar'"=="", "ep_seq", "`epvar'")
    quietly egen byte `tagep' = tag(cid `epv') if `varlist'==1 & `esmp'==1
    quietly count if `tagep'==1
    return scalar nep = r(N)
    quietly egen byte `tagcid' = tag(cid) if `varlist'==1 & `esmp'==1
    quietly count if `tagcid'==1
    return scalar ncid = r(N)
end

* ══════════════════════════════════════════════════════════════════════════
* 1. BUILD THE BROADENED TREATMENT — dc_in_crisis_nd / _preempt / _post
* ══════════════════════════════════════════════════════════════════════════
di as result _n "════════════════════════════════════════════════════════════"
di as result "1. BROADENED DEBT-CRISIS TREATMENT — build and diagnostics"
di as result "════════════════════════════════════════════════════════════"

* Conflict rows: spread says non-default, AT independently records a default
* window over the same row. Trusted to AT (the reason this file exists),
* but printed, never silently resolved.
capture drop _dc_conflict
gen byte _dc_conflict = (in_crisis_nd==1 & at_default_year==1)
quietly count if _dc_conflict==1
if r(N) > 0 {
    di as error "  ** " r(N) " country-years classified non-default by the spread database"
    di as error "     but flagged as an AT default window -- RECLASSIFIED to default-linked."
    di as error "     Listed below, not silently overridden:"
    list cid country year at_episode_onset at_type if _dc_conflict==1, noobs sepby(cid)
}

* New, entirely additive AT-only default years: tranquil in the existing
* flow coding, but AT records a restructuring in progress.
capture drop _dc_atonly
gen byte _dc_atonly = (in_crisis==0 & at_default_year==1)
quietly count if _dc_atonly==1
di as result "  AT-only default years (tranquil under the existing spread coding," ///
    " now treated): `r(N)'"
if r(N) > 0 {
    quietly egen byte _dc_atonly_ep = tag(cid at_episode_onset) if _dc_atonly==1
    quietly count if _dc_atonly_ep==1
    di as result "  ... spanning `r(N)' distinct AT-only episodes."
    capture drop _dc_atonly_ep
    * Split by country origin: the 52-country spread panel vs. the
    * 10b_skeleton_atonly.do expansion. Only the latter can ever appear here
    * with dc_in_crisis_nd structurally 0 for the whole country (no spread
    * data exists for it) -- see METHODOLOGY.md §5.3.
    quietly count if _dc_atonly==1 & atonly_country==0
    local n_at52 = r(N)
    quietly count if _dc_atonly==1 & atonly_country==1
    local n_atnew = r(N)
    di as result "      of which `n_at52' from the original 52-country panel," ///
        " `n_atnew' from the AT-only country expansion (10b)."
    if `n_atnew' > 0 {
        quietly egen byte _dc_newctry_ep = tag(cid at_episode_onset) ///
            if _dc_atonly==1 & atonly_country==1
        quietly count if _dc_newctry_ep==1
        di as result "      expansion countries contributing: `r(N)' episodes."
        capture drop _dc_newctry_ep
    }
}

capture drop dc_default_year dc_in_crisis dc_in_crisis_nd
gen byte dc_default_year = (in_crisis_def==1) | (_dc_atonly==1) | (_dc_conflict==1)
gen byte dc_in_crisis_nd = (in_crisis_nd==1) & !_dc_conflict
gen byte dc_in_crisis    = dc_in_crisis_nd | dc_default_year
label var dc_default_year "Broadened default-linked treatment (spread OR AT-sourced)"
label var dc_in_crisis_nd "Broadened non-default treatment (spread-only, unchanged)"
label var dc_in_crisis    "Broadened debt-crisis treatment (nd or default-linked)"

* ── 2-YEAR RESTRUCTURING-GAP RULE ───────────────────────────────────────────
* If a spread-based default episode's earliest overlapping AT restructuring
* window starts more than 2 years after the spread episode's own onset, the
* years BEFORE that AT onset are reclassified non-default (the market-stress
* episode existed, but no default/restructuring had actually begun yet);
* years from the AT onset onward keep their default-linked status and AT
* type. Applies only to episodes with at least one AT-overlapping row --
* episodes with NO AT overlap anywhere (the "Unclassified" cases: Zambia,
* Ghana, Lebanon, Sri Lanka, Ukraine, as of this AT vintage) have no
* independent restructuring date to test against and are left as-is (see
* this file's header for why).
capture drop _dc_ep_onset_yr _dc_ep_at_yr _dc_gap _dc_predefault_reclass
bysort cid ep_seq: egen int _dc_ep_onset_yr = min(cond(onset_all==1, year, .)) ///
    if in_crisis_def==1
bysort cid ep_seq: egen int _dc_ep_at_yr = min(cond(at_default_year==1, at_episode_onset, .)) ///
    if in_crisis_def==1
gen int _dc_gap = _dc_ep_at_yr - _dc_ep_onset_yr if in_crisis_def==1
gen byte _dc_predefault_reclass = ///
    in_crisis_def==1 & !missing(_dc_gap) & _dc_gap > 2 & year < _dc_ep_at_yr

quietly egen byte _dc_split_ep_tag = tag(cid ep_seq) ///
    if in_crisis_def==1 & !missing(_dc_gap) & _dc_gap > 2
quietly count if _dc_split_ep_tag==1
if r(N) > 0 {
    di as result "  " r(N) " default-linked episode(s) split under the 2-year restructuring-gap rule:"
    list cid country _dc_ep_onset_yr _dc_ep_at_yr _dc_gap if _dc_split_ep_tag==1, noobs
    quietly count if _dc_predefault_reclass==1
    di as result "  ... reclassifying " r(N) " country-years from default-linked to non-default."
}
capture drop _dc_split_ep_tag

* Rebuild dc_default_year / dc_in_crisis_nd / dc_in_crisis with the
* reclassification applied -- everything downstream (the dc_ep_seq
* running-counter below, the dc_type mode-fill, the 3-arm regression) reads
* these generically, so no other logic needs to change: dc_default_year now
* correctly toggles off at the pre-restructuring years and back on at the AT
* onset year, and the EXISTING running-counter (increments whenever
* dc_default_year turns on after being off) automatically creates two
* separate episodes out of a split one.
drop dc_default_year dc_in_crisis_nd dc_in_crisis
gen byte dc_default_year = ///
    (in_crisis_def==1 & !_dc_predefault_reclass) | (_dc_atonly==1) | (_dc_conflict==1)
gen byte dc_in_crisis_nd = ///
    ((in_crisis_nd==1) | (in_crisis_def==1 & _dc_predefault_reclass==1)) & !_dc_conflict
gen byte dc_in_crisis = dc_in_crisis_nd | dc_default_year
label var dc_default_year "Broadened default-linked treatment (spread OR AT-sourced; 2yr-gap rule applied)"
label var dc_in_crisis_nd "Broadened non-default treatment (spread-only + pre-restructuring reclassified years)"
label var dc_in_crisis    "Broadened debt-crisis treatment (nd or default-linked)"
drop _dc_ep_onset_yr _dc_ep_at_yr _dc_gap _dc_predefault_reclass

* ── Episode grouping for the broadened default treatment ───────────────────
* Existing in_crisis_def rows keep their own ep_seq (already built in
* 18_transforms.do). AT-only default runs need their OWN episode grouping,
* since they are not part of any existing ep_seq. Build a combined counter
* the same running-sum idiom 18_transforms.do uses for ep_seq: increments at
* every row where dc_default_year turns on after being off in the prior year
* for that country (whether the trigger is spread-based or AT-only).
capture drop dc_default_onset dc_ep_seq
bysort cid (year): gen byte dc_default_onset = ///
    dc_default_year==1 & (dc_default_year[_n-1]!=1 | cid!=cid[_n-1])
bysort cid (year): gen int dc_ep_seq = sum(dc_default_onset)

* ── Type fill: preemptive/post-default, episode-level ───────────────────────
* Row-level at_type where a direct AT-window match exists; forward/backward-
* filled across the whole default episode from ANY row that does match,
* mirroring nd_ep's fill logic in 18_transforms.do. Episodes with no AT
* window match anywhere are left unclassified and reported, not guessed.
capture drop dc_type
gen str24 _at_type_nonmissing = at_type if dc_default_year==1 & !missing(at_type)
bysort cid dc_ep_seq: egen str24 dc_type = mode(_at_type_nonmissing) if dc_default_year==1 & dc_ep_seq>0
replace dc_type = "Unclassified (no AT type match)" if dc_default_year==1 & missing(dc_type)
drop _at_type_nonmissing

quietly egen byte _dc_ep_tag = tag(cid dc_ep_seq) if dc_default_year==1 & dc_ep_seq>0
di as result _n "  Broadened default-linked episodes, by type (episode count, not row count):"
foreach t in "Strictly preemptive" "Weakly preemptive" "Post-default" "Unclassified (no AT type match)" {
    quietly count if _dc_ep_tag==1 & dc_type=="`t'"
    di as result "    `t': " r(N)
}
capture drop _dc_ep_tag

capture drop dc_in_crisis_preempt dc_in_crisis_post
gen byte dc_in_crisis_preempt = dc_default_year==1 & ///
    inlist(dc_type, "Strictly preemptive", "Weakly preemptive")
gen byte dc_in_crisis_post = dc_default_year==1 & dc_type=="Post-default"
label var dc_in_crisis_preempt "Broadened treatment: preemptive default-linked (strictly or weakly)"
label var dc_in_crisis_post    "Broadened treatment: post-default"

quietly count if dc_default_year==1 & dc_type=="Unclassified (no AT type match)"
if r(N) > 0 di as error "  ** " r(N) " default-linked country-years have NO preemptive/post-default" ///
    " type (excluded from both dc_in_crisis_preempt and dc_in_crisis_post -- they" ///
    " still count in dc_default_year/dc_in_crisis, but not in the 3-arm regression below)."

* ══════════════════════════════════════════════════════════════════════════
* 2. CONTROLS — UNCHANGED, reads $ctrl_flow exactly as built in 18_transforms.do
* ══════════════════════════════════════════════════════════════════════════
local flow_ctrl_variant 0
if `flow_ctrl_variant'==1 & "$ctrl_core_flowplus"=="" {
    di as error "  ** flow_ctrl_variant==1 requested but \$ctrl_core_flowplus is empty (exchange2"
    di as error "     unavailable, exch missing) -- re-run 01_build_panel.do/12_wdi.do/18_transforms.do"
    di as error "     after confirming data/raw/officialexchangerate.xlsx is present, or use 0."
    exit 111
}
local controls = cond(`flow_ctrl_variant'==1, "$ctrl_flow_flowplus", "$ctrl_flow")

local drop_year_fe 0
local yearfe = cond(`drop_year_fe', "", "i.year")

* ══════════════════════════════════════════════════════════════════════════
* 3. ESTIMATION — three arms jointly, tranquil omitted
* ══════════════════════════════════════════════════════════════════════════
foreach g in nd preempt post {
    foreach m in b se lo90 hi90 {
        matrix `m'_dc_`g' = J(7, 1, .)
        matrix `m'_dc_`g'[2,1] = 0
    }
}

tempname F
tempfile dcf
postfile `F' str16 arm int hdisp double b double se double p long nrow long nep long ncid using "`dcf'", replace

eststo clear
di as result _n "════════════════════════════════════════════════════════════"
di as result "3. BROADENED DEBT-CRISIS LP — three arms, tranquil omitted"
di as result "════════════════════════════════════════════════════════════"
di as result "h    ND        SE      PREEMPT   SE      POST      SE      p(post=preempt) p(post=nd) p(preempt=nd)"

forvalues h = 0/4 {
    local hd  = `h' + 1
    local row = `h' + 3
    local lag = max(2, `h'+3)

    capture noisily xtscc dy_`h' dc_in_crisis_nd dc_in_crisis_preempt dc_in_crisis_post ///
        `controls' `yearfe' if sample_flow_bc==1, fe lag(`lag')
    if _rc {
        di as error "  ** h=`hd' failed (rc=" _rc ")"
        continue
    }

    local bnd = _b[dc_in_crisis_nd]
    local snd = _se[dc_in_crisis_nd]
    local bpr = _b[dc_in_crisis_preempt]
    local spr = _se[dc_in_crisis_preempt]
    local bpo = _b[dc_in_crisis_post]
    local spo = _se[dc_in_crisis_post]
    _critvals
    local c90 = r(c90)
    _pval `bnd' `snd'
    local pnd = r(p)
    _pval `bpr' `spr'
    local ppr = r(p)
    _pval `bpo' `spo'
    local ppo = r(p)

    quietly lincom dc_in_crisis_post - dc_in_crisis_preempt
    local pd_postpre = 2*ttail(r(df), abs(r(estimate)/r(se)))
    quietly lincom dc_in_crisis_post - dc_in_crisis_nd
    local pd_postnd = 2*ttail(r(df), abs(r(estimate)/r(se)))
    quietly lincom dc_in_crisis_preempt - dc_in_crisis_nd
    local pd_prend = 2*ttail(r(df), abs(r(estimate)/r(se)))

    eststo dc_h`h'
    estadd scalar pdpp = `pd_postpre'
    estadd scalar pdpn = `pd_postnd'
    estadd scalar pdrn = `pd_prend'

    _nflowcount dc_in_crisis_nd,      outcome(dy_`h') controls(`controls') samp(sample_flow_bc) epvar(ep_seq)
    local nend  = r(nep)
    local nendn = r(nrow)
    local nendc = r(ncid)
    _nflowcount dc_in_crisis_preempt, outcome(dy_`h') controls(`controls') samp(sample_flow_bc) epvar(dc_ep_seq)
    local nepr  = r(nep)
    local neprn = r(nrow)
    local neprc = r(ncid)
    _nflowcount dc_in_crisis_post,    outcome(dy_`h') controls(`controls') samp(sample_flow_bc) epvar(dc_ep_seq)
    local nepo  = r(nep)
    local nepon = r(nrow)
    local nepoc = r(ncid)

    matrix b_dc_nd[`row',1]      = `bnd'
    matrix se_dc_nd[`row',1]     = `snd'
    matrix lo90_dc_nd[`row',1]   = `bnd' - `c90'*`snd'
    matrix hi90_dc_nd[`row',1]   = `bnd' + `c90'*`snd'
    matrix b_dc_preempt[`row',1]    = `bpr'
    matrix se_dc_preempt[`row',1]   = `spr'
    matrix lo90_dc_preempt[`row',1] = `bpr' - `c90'*`spr'
    matrix hi90_dc_preempt[`row',1] = `bpr' + `c90'*`spr'
    matrix b_dc_post[`row',1]    = `bpo'
    matrix se_dc_post[`row',1]   = `spo'
    matrix lo90_dc_post[`row',1] = `bpo' - `c90'*`spo'
    matrix hi90_dc_post[`row',1] = `bpo' + `c90'*`spo'

    di as result "  " %1.0f `hd' "  " %7.3f `bnd' " " %6.3f `snd' "  " ///
        %7.3f `bpr' " " %6.3f `spr' "  " %7.3f `bpo' " " %6.3f `spo' "   " ///
        %5.3f `pd_postpre' "        " %5.3f `pd_postnd' "     " %5.3f `pd_prend'

    post `F' ("nd")      (`hd') (`bnd') (`snd') (`pnd') (`nendn') (`nend') (`nendc')
    post `F' ("preempt") (`hd') (`bpr') (`spr') (`ppr') (`neprn') (`nepr') (`neprc')
    post `F' ("post")    (`hd') (`bpo') (`spo') (`ppo') (`nepon') (`nepo') (`nepoc')
}
postclose `F'

di as result _n "  p(post=preempt)/p(post=nd)/p(preempt=nd) are covariance-correct lincom"
di as result "  tests within the same joint regression -- no bootstrap/Clogg-z machinery"
di as result "  needed here (plain OLS, no propensity model, matching 20/23's own design)."

di as result _n "  Per-arm N / episodes / countries, by horizon (from debtcrisis_flow.csv):"
di as result "  h    arm       N   episodes  countries"
preserve
    use "`dcf'", clear
    forvalues h = 1/5 {
        foreach a in nd preempt post {
            quietly summarize nrow if arm=="`a'" & hdisp==`h'
            local nn = r(mean)
            quietly summarize nep if arm=="`a'" & hdisp==`h'
            local ne = r(mean)
            quietly summarize ncid if arm=="`a'" & hdisp==`h'
            local nc = r(mean)
            di as result "  " %1.0f `h' "    " %-8s "`a'" %5.0f `nn' "     " %5.0f `ne' "      " %5.0f `nc'
        }
    }
restore

* ══════════════════════════════════════════════════════════════════════════
* 4. EXPORTS
* ══════════════════════════════════════════════════════════════════════════
capture esttab dc_h0 dc_h1 dc_h2 dc_h3 dc_h4 using "$tabs/table_debtcrisis_flow.rtf", replace ///
    b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(dc_in_crisis_nd dc_in_crisis_preempt dc_in_crisis_post) ///
    order(dc_in_crisis_nd dc_in_crisis_preempt dc_in_crisis_post) ///
    coeflabel(dc_in_crisis_nd "Non-default" dc_in_crisis_preempt "Preemptive default" ///
              dc_in_crisis_post "Post-default") ///
    mtitles("h=1" "h=2" "h=3" "h=4" "h=5") nonumber ///
    stats(pdpp pdpn pdrn N N_g, ///
          labels("p (post=preempt)" "p (post=nd)" "p (preempt=nd)" "Observations" "Countries") ///
          fmt(3 3 3 0 0)) ///
    title("Broadened debt-crisis LP (flow, 3 arms): output cost by type") ///
    addnotes("Dependent variable: cumulative change in log real GDP (pp) from t-1 to t+h. All three treatment dummies entered" ///
             "jointly, tranquil years omitted. Treatment = broadened debt-crisis taxonomy (17b_merge_at_full.do): non-default" ///
             "unchanged (spread-database only); default-linked = spread-triggered OR AT-sourced, sub-split preemptive/post-" ///
             "default from the full Asonuma-Trebesch (2016) database. Country and year FE; Driscoll-Kraay SE at" ///
             "lag(max(2,h+3)). Pairwise p-values are covariance-correct lincom tests within the joint regression." ///
             "* p<0.10, ** p<0.05, *** p<0.01.")
if _rc == 608 di as error "  ** table_debtcrisis_flow.rtf is OPEN IN WORD — close it and re-run."
else if _rc   di as error "  ** esttab failed (rc=" _rc ")"
else          di as result "Table saved: $tabs/table_debtcrisis_flow.rtf"

* ── 4b. LAYOUT MATCHING THE REFERENCE PAPER'S TABLE 2 ───────────────────────
* Coefficient (with stars) row, SE row, and an Observations/Countries/
* Episodes row directly beneath each arm's own SE row, columns = h=1..h=5 --
* esttab's stats() rows sit once at the BOTTOM of the whole table, indexed
* by MODEL (column) not by COEFFICIENT (row), so they cannot produce a
* different Obs/Countries/Episodes line per arm interleaved under that
* arm's own block. Written directly via file write. Tab-stopped RTF text
* (not a full \trowd grid table): simpler to generate correctly without a
* live Stata session to verify cell-width arithmetic against.
capture program drop _starstr
program define _starstr, rclass
    args pval
    local s ""
    if !missing(`pval') {
        if `pval' < .01       local s "***"
        else if `pval' < .05  local s "**"
        else if `pval' < .10  local s "*"
    }
    return local stars "`s'"
end

preserve
    use "`dcf'", clear
    forvalues h = 1/5 {
        foreach a in nd preempt post {
            quietly summarize b if arm=="`a'" & hdisp==`h'
            local b_`a'_`h' = r(mean)
            quietly summarize se if arm=="`a'" & hdisp==`h'
            local se_`a'_`h' = r(mean)
            quietly summarize p if arm=="`a'" & hdisp==`h'
            local p_`a'_`h' = r(mean)
            quietly summarize nrow if arm=="`a'" & hdisp==`h'
            local nrow_`a'_`h' = r(mean)
            quietly summarize nep if arm=="`a'" & hdisp==`h'
            local nep_`a'_`h' = r(mean)
            quietly summarize ncid if arm=="`a'" & hdisp==`h'
            local ncid_`a'_`h' = r(mean)
        }
    }
restore

capture file close dctab
file open dctab using "$tabs/table_debtcrisis_flow_layout.rtf", write replace
file write dctab "{\rtf1\ansi\deff0" _n
file write dctab "{\b Broadened debt-crisis LP (flow): output cost by type\par}" _n
file write dctab "{\i Observations/Countries/Episodes reported beneath each arm's standard error, matching the reference" _n
file write dctab " paper's Table 2 layout.\par}" _n
file write dctab "\par" _n
file write dctab "\tab h = 1\tab h = 2\tab h = 3\tab h = 4\tab h = 5\par" _n
file write dctab "\par" _n

foreach spec in "nd Non-default" "preempt Preemptive default" "post Post-default" {
    gettoken key lbl : spec
    file write dctab "{\b `lbl'}\par" _n

    local coefline ""
    local seline ""
    local cntline ""
    forvalues h = 1/5 {
        _starstr `p_`key'_`h''
        local st = r(stars)
        local bstr : display %5.2f `b_`key'_`h''
        local sestr : display %5.2f `se_`key'_`h''
        local nn : display %4.0f `nrow_`key'_`h''
        local cc : display %3.0f `ncid_`key'_`h''
        local ee : display %3.0f `nep_`key'_`h''
        local coefline "`coefline'\tab `bstr'`st'"
        local seline   "`seline'\tab (`sestr')"
        local cntline  "`cntline'\tab `nn'/`cc'/`ee'"
    }
    file write dctab "`coefline'\par" _n
    file write dctab "`seline'\par" _n
    file write dctab "Observations/Countries/Episodes`cntline'\par" _n
    file write dctab "\par" _n
}
file write dctab "{\i * p<0.10, ** p<0.05, *** p<0.01. Countries and episodes are counted per arm at each horizon" _n
file write dctab " -- a country/episode contributes to an arm's count only if it has at least one treated" _n
file write dctab " observation of that type surviving listwise deletion at that horizon, so these totals can" _n
file write dctab " shrink slightly across h=1 to h=5 as later-horizon outcomes go missing.\par}" _n
file write dctab "}" _n
file close dctab
di as result "Layout-matched table saved: $tabs/table_debtcrisis_flow_layout.rtf"

preserve
    use "`dcf'", clear
    label var arm  "nd / preempt / post"
    label var hdisp "Horizon (1 = crisis year)"
    label var nrow "Observations (country-years) in this arm at this horizon"
    label var nep  "Distinct episodes"
    label var ncid "Distinct countries"
    export delimited "$tabs/debtcrisis_flow.csv", replace
    di as result "Raw coefficients saved: $tabs/debtcrisis_flow.csv"
restore

* ── IRF datasets + figure ────────────────────────────────────────────────
foreach g in nd preempt post {
    preserve
        clear
        set obs 7
        gen horizon = _n - 2
        foreach m in b lo90 hi90 {
            svmat `m'_dc_`g', names(`m')
            rename `m'1 `m'
        }
        gen series = "`g'"
        save "$clean/irf_debtcrisis_flow_`g'.dta", replace
    restore
}

local c_nd "0 84 166"
local c_pr "230 159 0"
local c_po "157 36 73"
preserve
    use "$clean/irf_debtcrisis_flow_nd.dta", clear
    append using "$clean/irf_debtcrisis_flow_preempt.dta"
    append using "$clean/irf_debtcrisis_flow_post.dta"
    keep if horizon >= 0
    twoway ///
        (connected b horizon if series=="nd",      lcolor("`c_nd'") mcolor("`c_nd'") msymbol(circle) lwidth(medthick)) ///
        (connected b horizon if series=="preempt", lcolor("`c_pr'") mcolor("`c_pr'") msymbol(triangle) lpattern(dash) lwidth(medthick)) ///
        (connected b horizon if series=="post",    lcolor("`c_po'") mcolor("`c_po'") msymbol(square) lpattern(shortdash) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(gs8) lwidth(thin)) ///
        xlabel(0(1)5, labsize(medsmall)) ylabel(, format(%4.1f) labsize(medsmall)) ///
        xtitle("Year (Year 1 = crisis year)", size(medsmall)) ///
        ytitle("Cum. change in log real GDP (pp)", size(medsmall)) ///
        title("Output Cost of a Broadened Sovereign Debt Crisis", size(medium)) ///
        subtitle("Non-default vs preemptive default vs post-default. Spread crisis (existing) + AT-only" ///
                 " default years (17b_merge_at_full.do).", size(small)) ///
        legend(order(1 "Non-default" 2 "Preemptive default" 3 "Post-default") ring(0) pos(7) cols(1) size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    capture graph export "$figs/fig_debtcrisis_flow.pdf", replace
    if _rc di as error "  ** fig_debtcrisis_flow.pdf export failed (rc=" _rc ")"
    else {
        capture graph export "$figs/fig_debtcrisis_flow.png", replace width(1400)
        di as result "Figure saved: fig_debtcrisis_flow.pdf/.png"
    }
restore

di as result _n "26_lp_debtcrisis_flow.do complete."
di as result "  Compare the ND column to 20_lp_flow.do's own in_crisis_nd column: it should be"
di as result "  close (identical unless a spread/AT conflict row reclassified a handful of"
di as result "  observations, printed in Section 1 above) -- the sanity check that the"
di as result "  broadened design did not silently change what 'non-default' means."
