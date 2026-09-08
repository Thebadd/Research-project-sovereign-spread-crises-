/*===========================================================================
  18B_MINDUR_VARIANT.DO
  Robustness variant #1 of the annual-criterion episode-dating rule
  (in_crisis, 18_transforms.do:455): a crisis year only counts if it is part
  of a run of >=2 CONSECUTIVE crisis_any==1 years. A lone one-year spike is
  RECLASSIFIED AS TRANQUIL (folded into the control pool), not dropped --
  per the user's own confirmed choice, discussed and settled in this
  project's own planning conversation (AskUserQuestion: "Reclassify as
  tranquil" over "drop from the sample entirely").

  WHY THIS EXISTS: the project's baseline in_crisis (annual criterion, no
  persistence/gap bridging -- 18_transforms.do:455) treats every year that
  independently clears the spread criterion, including a single noisy
  one-year spike with no sustained crisis around it. This variant tests
  whether the project's headline results are being driven by such spikes,
  by removing them from the treated set entirely (not just widening their
  CI) and re-checking how many onsets/rows actually change status.

  STANDALONE: reads panel_lp.dta, builds NEW columns under distinct names
  (in_crisis_mindur, not in_crisis) -- does NOT modify 18_transforms.do,
  in_crisis, onset_all, continuation, or any existing baseline column.
  Does NOT re-estimate any LP/AIPW file in this pass -- diagnostic/
  column-construction only; wiring a robustness row into the estimation
  files is a follow-up once this column is validated against the console
  diagnostics below.

  Run AFTER 18_transforms.do. Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ── Run-length of consecutive crisis_any==1 years within each country ──────
* Relies on the same "bysort cid (year)" ordering convention already used
* for ep_seq (18_transforms.do:481-482) -- consistent with this project's
* existing assumption that the annual panel has no within-country year gaps
* once carryin scaffolding is accounted for.
capture drop _newrun _run_id _run_len
bysort cid (year): gen byte _newrun = (crisis_any==1) & (crisis_any[_n-1]!=1 | _n==1)
bysort cid (year): gen long _run_id = sum(_newrun) if crisis_any==1
bysort cid _run_id (year): gen long _run_len = _N if crisis_any==1

* ── MINIMUM-DURATION annual-criterion treatment ─────────────────────────────
* crisis_any==1 AND run length >=2. A lone one-year spike (_run_len==1) is
* reclassified tranquil (folded into the control pool), NOT dropped.
gen byte in_crisis_mindur = (crisis_any==1 & _run_len>=2 & carryin==0)
label var in_crisis_mindur "Annual-criterion treatment, single-year spikes (<2 consecutive yrs) reclassified tranquil (robustness)"

* ── Diagnostics ──────────────────────────────────────────────────────────
di as result _n "=== MINIMUM-DURATION ROBUSTNESS VARIANT (>=2 consecutive crisis_any years) ==="

quietly count if in_crisis==1
local n_base = r(N)
quietly count if in_crisis_mindur==1
local n_mindur = r(N)
quietly count if in_crisis==1 & in_crisis_mindur==0
local n_reclass = r(N)

di as result "Baseline in_crisis treated rows:            " `n_base'
di as result "Minimum-duration in_crisis_mindur treated:  " `n_mindur'
di as result "Reclassified treated -> tranquil (spikes):  " `n_reclass'

* By resolution type -- a spike-heavy arm losing disproportionately many
* episodes is worth knowing before any re-estimation on this sample.
foreach ndval in 1 0 {
    local lbl = cond(`ndval'==1, "Non-default", "Default-linked")
    quietly count if in_crisis==1 & in_crisis_mindur==0 & nd_ep==`ndval'
    di as result "  Reclassified, `lbl' arm: " r(N)
}

* Onset-level count: how many distinct episodes (ep_seq groups) are
* entirely single-year spikes and vanish under this variant.
quietly egen byte _spike_ep_tag = tag(cid ep_seq) if in_crisis==1 & in_crisis_mindur==0 & onset_all==1
quietly count if _spike_ep_tag==1
di as result "Onsets that are themselves single-year spikes (whole episode reclassified): " r(N)

* Full list of reclassified country-years, for manual auditing.
di as result _n "Reclassified country-years (treated under baseline, tranquil under min-duration):"
list country year crisis_any nd_ep if in_crisis==1 & in_crisis_mindur==0, noobs sepby(country)

capture drop _newrun _run_id _run_len _spike_ep_tag

di as result _n "18b_mindur_variant.do complete."
di as result "New column: in_crisis_mindur (not wired into any estimation file yet)."
