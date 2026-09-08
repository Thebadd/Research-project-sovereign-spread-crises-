/*===========================================================================
  13E_NEXUS_BARS.DO
  Descriptive companion to 13d_aipw_nexus_split.do: a country-level bar
  chart of PRE-CRISIS SOVEREIGN-BANK NEXUS EXPOSURE (a_nexus), one bar per
  crisis onset, split by resolution type (Non-default vs Default-linked),
  sorted descending, with a horizontal line at the median across ALL onsets
  pooled. This is the descriptive analog of Asonuma, Chamon, Erce & Sasahara
  (2024)'s own Panel C/D bar-chart figure (their split: weakly vs strictly
  preemptive, on bank credit/GDP) -- adapted here to THIS project's actual
  split variable and axis: the sovereign-bank nexus, split by non-default vs
  default-linked spread crisis, matching 13d's own headline design rather
  than the paper's preemptive-type distinction (this project does not use
  preemptive/strictly-preemptive as its main heterogeneity axis; the
  sovereign-bank "doom loop" nexus is this project's own contribution).

  AMPLIFIER: identical construction to 13d_aipw_nexus_split.do (lines
  166-173 there) -- copied verbatim, not re-derived: pre-crisis
  claims-on-government/bank-assets ratio (L.claimsgov_assets), country-mean
  filled where the t-1 value is missing (coverage is thin, 2001+). The
  median line is computed over ALL onsets pooled (sample==1 & onset_all==1),
  matching 13d's own cutoff population and the reference paper's own
  "median of all restructurings" line (computed before restricting to one
  panel's own subset).

  PURELY DESCRIPTIVE: this file does not estimate anything. It reads
  panel_lp.dta and reuses 13d's amplifier construction read-only; it does
  not touch 13d_aipw_nexus_split.do or change any AIPW result.

  Output: $figs/fig_nexus_bar_nd.pdf (non-default onsets), $figs/fig_nexus_
          bar_def.pdf (default-linked onsets).
  Run AFTER 18_transforms.do. Not wired into 00_master.do (standalone
  descriptive export, matching this project's convention for companion
  figure files).
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

* ── Amplifier: identical construction to 13d_aipw_nexus_split.do:166-173 ──
capture drop a_nexus a_nexus_cm
gen a_nexus = L.claimsgov_assets
bysort cid: egen a_nexus_cm = mean(claimsgov_assets)
replace a_nexus = a_nexus_cm if missing(a_nexus)

* Median over ALL onsets (matching 13d's own cutoff population, and the
* reference paper's own "median of all restructurings" line, computed
* before restricting to one panel's own subset).
quietly summarize a_nexus if sample==1 & onset_all==1, detail
local medall = r(p50)
di as result _n "=== NEXUS BAR CHART: pre-crisis sovereign-bank nexus by resolution type ==="
di as result "Median pre-crisis nexus, all onsets pooled = " %5.2f `medall'

* One row per onset.
keep if sample==1 & onset_all==1
drop if missing(a_nexus)
quietly count
local npooled = r(N)
gen str24 iso_year = country + " " + string(year)

foreach t in "nd" "def" {
    local tlab = cond("`t'"=="nd", "Non-default", "Default-linked")
    quietly count if onset_`t'==1
    local nthis = r(N)
    di as result "  `tlab': " `nthis' " onsets with non-missing nexus"
    preserve
        keep if onset_`t' == 1
        gsort -a_nexus
        local medstr : display %4.1f `medall'
        * text() x-position: roughly a third of the way in from the left,
        * where bars have already dropped near/below the median line, so
        * the label sits in open space to the RIGHT of that point (place(e))
        * and stays fully inside the plot -- not spanning left across the
        * tall bars/into the y-axis title, which is what clipped it before.
        quietly count
        local xtxt = max(1, round(`r(N)' * 0.30))
        local ytxt = `medall' + 5
        capture noisily graph bar a_nexus, ///
            over(iso_year, sort(a_nexus) descending label(angle(45) labsize(vsmall))) ///
            bar(1, color("121 168 208") lcolor(gs8)) ///
            yline(`medall', lcolor(navy)) ///
            text(`ytxt' `xtxt' "Median of all spread crises = `medstr'%", ///
                place(e) size(small) color(black) justification(left)) ///
            ytitle("Bank claims on government-to-total asset", size(small)) ///
            title("`tlab'", size(medium) color(navy)) ///
            graphregion(color(white)) bgcolor(white) ///
            ysize(3) xsize(5)
        if _rc == 0 {
            graph export "$figs/fig_nexus_bar_`t'.pdf", replace
            di as result "Figure saved: fig_nexus_bar_`t'.pdf"
        }
        else di as error "  ** fig_nexus_bar_`t' failed (rc=" _rc ")"
    restore
}

di as result _n "13e_nexus_bars.do complete."
di as result "fig_nexus_bar_nd.pdf / fig_nexus_bar_def.pdf: one bar per onset,"
di as result "pre-crisis sovereign-bank nexus, sorted descending, navy line ="
di as result "the median across all `npooled' onsets pooled (both types), matching"
di as result "13d_aipw_nexus_split.do's own median-split cutoff population."
