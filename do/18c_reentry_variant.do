/*===========================================================================
  18C_REENTRY_VARIANT.DO
  Robustness variant #2 of the annual-criterion episode-dating rule
  (in_crisis, 18_transforms.do:455): a crisis persists from its onset until
  the spread falls back to within 20% of a PRE-CRISIS REFERENCE LEVEL, not
  until crisis_any itself flips off. This is a genuinely different
  assumption from the baseline (which treats only years that individually
  clear the spread criterion) -- here, once a crisis starts, it is treated
  as ongoing through any elevated-but-sub-threshold years until the spread
  actually recovers to something close to what it was before the crisis.

  PRE-CRISIS REFERENCE LEVEL, per onset, per the user's own confirmed
  if/else rule (not a weighted blend):
    - LOCAL reference = mean spread over the TRANQUIL years strictly
      between the previous episode and this onset, if that gap has >=3
      tranquil years.
    - Else FALL BACK to the country's full-sample tranquil-year mean
      spread (all tranquil years anywhere in that country's history).
  Implementation note: ep_seq increments AT the onset row (running sum of
  onset_all, 18_transforms.do:481-482), so the gap immediately preceding
  onset k's own ep_seq value carries ep_seq = k-1 -- the local reference is
  therefore built from ep_seq-1's tranquil rows, not the onset's own
  ep_seq group (which would wrongly include this episode's own aftermath).

  RE-ENTRY THRESHOLD: spread <= 1.2 * pre-crisis reference ends the crisis
  (20% tolerance band, per the user's own confirmed choice) -- the current
  year's actual spread level (spr_mean), not a lag, since this is asking
  "has this year's spread recovered," not predicting next year's.

  STANDALONE: reads panel_lp.dta, builds NEW columns under distinct names
  (pre_ref, in_crisis_reentry, sample_reentry) -- does NOT modify
  18_transforms.do, in_crisis, onset_all, continuation, sample, or any
  other existing baseline column. Saves $clean/panel_lp_reentry.dta (does
  NOT overwrite panel_lp.dta). Does NOT re-estimate any LP/AIPW file in
  this pass -- diagnostic/column-construction only; wiring sample_reentry
  into an Act 2 re-estimation is a follow-up (see companion files).

  Run AFTER 18_transforms.do. Not wired into 00_master.do.
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ── Pre-crisis reference level per onset ────────────────────────────────
* Local reference: mean spr_mean over TRANQUIL years carrying ep_seq =
* (this onset's ep_seq - 1) -- exactly the gap after the previous episode
* and before this one. Needs >=3 such tranquil years; else falls back.
preserve
    keep if in_crisis==0 & carryin==0
    collapse (mean) _local_mean=spr_mean (count) _local_n=spr_mean, by(cid ep_seq)
    rename ep_seq _ep_seq_gap
    tempfile _localref
    save `_localref'
restore

capture drop _ep_seq_gap _local_mean _local_n
gen long _ep_seq_gap = ep_seq - 1
merge m:1 cid _ep_seq_gap using `_localref', keep(master match) nogen

capture drop _fullref _fullref_fill
bysort cid: egen double _fullref = mean(spr_mean) if in_crisis==0 & carryin==0
bysort cid: egen double _fullref_fill = mean(_fullref)

capture drop pre_ref
gen double pre_ref = cond(!missing(_local_n) & _local_n>=3, _local_mean, _fullref_fill)
label var pre_ref "Pre-crisis reference spread level (local gap mean if >=3 tranquil yrs, else country full-sample tranquil mean)"

di as result _n "=== RE-ENTRY-LEVEL ROBUSTNESS VARIANT (episode ends when spread <= 1.2x pre-crisis reference) ==="
quietly count if onset_all==1 & missing(pre_ref)
di as result "Onsets with NO usable pre-crisis reference (country has no tranquil history at all): " r(N)
quietly count if onset_all==1 & !missing(_local_n) & _local_n>=3
di as result "Onsets using the LOCAL gap reference (>=3 tranquil years immediately before): " r(N)
quietly count if onset_all==1 & !missing(pre_ref) & (missing(_local_n) | _local_n<3)
di as result "Onsets falling back to the country FULL-SAMPLE tranquil reference: " r(N)

* ── Re-entry-based episode persistence ──────────────────────────────────
* Starting at each onset, a row stays "in crisis" as long as spr_mean has
* remained continuously above 1.2*pre_ref since the onset (cumulative-AND /
* running-min: once a row drops into the band, every later row in that
* same ep_seq group stays out too, until the NEXT onset starts a new group).
capture drop _above _cum_above
gen byte _above = (spr_mean > 1.2*pre_ref) if ep_seq>=1 & !missing(pre_ref)
bysort cid ep_seq (year): replace _above = 1 if onset_all==1   // onset itself always starts the episode
bysort cid ep_seq (year): gen byte _cum_above = _above[1]
bysort cid ep_seq (year): replace _cum_above = min(_cum_above[_n-1], _above) if _n>1

capture drop in_crisis_reentry
gen byte in_crisis_reentry = (_cum_above==1) & carryin==0
replace in_crisis_reentry = 0 if missing(in_crisis_reentry)
label var in_crisis_reentry "Treatment persists from onset until spread <=1.2x pre-crisis reference (robustness)"

* ── Diagnostics: episode-length comparison vs baseline episode membership ──
capture drop _len_reentry _len_baseline
bysort cid ep_seq: egen long _len_reentry  = total(in_crisis_reentry) if ep_seq>=1
bysort cid ep_seq: egen long _len_baseline = total(in_crisis_episode) if ep_seq>=1

quietly count if in_crisis_reentry==1
di as result _n "Total treated rows, re-entry variant:    " r(N)
quietly count if in_crisis_episode==1
di as result "Total treated rows, baseline episode-membership (in_crisis_episode): " r(N)
quietly count if in_crisis==1
di as result "Total treated rows, baseline annual-criterion (in_crisis):           " r(N)

quietly count if onset_all==1 & _len_reentry > _len_baseline
di as result "Episodes LENGTHENED under re-entry rule (spread still above band after baseline's own end): " r(N)
quietly count if onset_all==1 & _len_reentry < _len_baseline
di as result "Episodes SHORTENED under re-entry rule (spread back in band before baseline's own end):      " r(N)
quietly count if onset_all==1 & _len_reentry == _len_baseline
di as result "Episodes UNCHANGED in length:                                                                " r(N)

di as result _n "Per-onset comparison (baseline episode length vs re-entry episode length):"
list country year pre_ref _len_baseline _len_reentry if onset_all==1, noobs sepby(country)

* ── sample_reentry: the onset-tier estimation SAMPLE under the re-entry
* episode definition ────────────────────────────────────────────────────
* CORRECTION (this session): re-entry was earlier claimed to have "no
* effect" on the onset-tier LP/AIPW design, on the reasoning that
* onset_nd/onset_def are 1 only at the single onset row and the outcome
* (F h.dy - L.dy) is measured forward from that row regardless of how
* long the crisis lasts -- both true, but incomplete. `sample' (18_
* transforms.do: (continuation==0) & !missing(ln_gdp_base) & carryin==0 &
* atonly_country==0) excludes bridged continuation years from the
* regression ENTIRELY, using the BASELINE bridging rule -- so every year
* the baseline currently treats as a TRANQUIL CONTROL year (because its
* own rule says the episode already ended) but the re-entry rule still
* counts as "in crisis" (spread still above 1.2x pre-crisis reference)
* contaminates the comparison group under the baseline sample; symmetrically
* a year the baseline still excludes as continuation==1 might have already
* recovered under re-entry and could validly join the tranquil pool. THIS
* is the real, non-trivial channel through which re-entry changes the
* onset-tier regression: not the treatment coding or the outcome, but the
* CONTROL-POOL COMPOSITION `sample' selects. sample_reentry keeps every
* onset row (still the treated row, in_crisis_reentry==1 there by
* construction) and every year NOT still "in crisis" under the re-entry
* rule; drops any post-onset year re-entry still counts as ongoing crisis,
* regardless of what the baseline's own continuation flag said about it.
capture drop sample_reentry
gen byte sample_reentry = ((in_crisis_reentry==0) | (onset_all==1)) & ///
    !missing(ln_gdp_base) & carryin==0 & atonly_country==0
label var sample_reentry "Onset-tier estimation sample under the re-entry episode definition (robustness)"

quietly count if sample==1 & sample_reentry==0
di as result _n "Rows in baseline sample but EXCLUDED under sample_reentry (still 'in crisis' per re-entry rule): " r(N)
quietly count if sample==0 & sample_reentry==1
di as result "Rows EXCLUDED from baseline sample but INCLUDED under sample_reentry (recovered earlier per re-entry rule): " r(N)
quietly count if sample==1 & sample_reentry==1
di as result "Rows in BOTH samples (unaffected by the swap): " r(N)

capture drop _ep_seq_gap _local_mean _local_n _fullref _fullref_fill _above _cum_above _len_reentry _len_baseline

sort cid year
save "$clean/panel_lp_reentry.dta", replace
di as result _n "Saved: $clean/panel_lp_reentry.dta (panel_lp.dta + pre_ref, in_crisis_reentry, sample_reentry)"

di as result _n "18c_reentry_variant.do complete."
di as result "New columns: pre_ref, in_crisis_reentry, sample_reentry."
