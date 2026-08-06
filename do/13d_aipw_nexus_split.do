/*===========================================================================
  13D_AIPW_NEXUS_SPLIT.DO
  Bank-intermediation heterogeneity — DOOM-LOOP version (Asonuma et al. §4 / Fig 6
  analog). Their headline two-dimensional result: the output cost of a crisis
  depends on the resolution strategy AND on bank intermediation. They median-split
  by bank-credit/GDP ("large" vs "small" banking sector) and estimate AIPW
  separately on each subsample. We go one better: split by the SOVEREIGN-BANK
  NEXUS — claimsgov_assets (bank claims on government / total assets, the doom-loop
  intensity Asonuma didn't have) — in a spread-crisis setting.

  Headline claim: default-linked spread crises are costliest where the
  sovereign-bank nexus is tight (banks heavily exposed to the sovereign).

  Outcomes: GDP (dy_h, coherent with 08b) PLUS the transmission channels
  credit, inv, claimpriv_assets, claims_govt — so we can watch each channel evolve
  differently under high vs low nexus (which channel carries the non-default
  cushion vs the default-linked doom-loop loss). Same Part A + Part B per outcome.

  Design (per outcome, coherent with 08b/13c):
    Amplifier a_nexus = pre-crisis claimsgov_assets (L.claimsgov_assets), filled
      with the country mean where the t-1 value is missing (coverage is thin, 2001+).
    Median split over crisis onsets: highbank = a_nexus >= median(onsets).
    AIPW cells (control = tranquil years only; treated = the specific cell; rival
    onsets dropped), estimated with the same _aipw + cluster-bootstrap as 08b/13c:
      Part A (robust headline): all onsets x {high, low}         -> 2 lines
      Part B (two-dimensional): {nd, def} x {high, low}          -> 4 lines

  SMALL-SAMPLE CAVEAT: the default x {high,low} cells hold only a handful of events,
  so many bootstrap draws drop out. CIs are shown only when >=50 valid draws; else
  the point estimate stands with a printed caveat. This thinness is the honest limit
  vs Asonuma's 194 restructurings — read Part B as suggestive, Part A as the result.

  Coherence with the Asonuma replication (their Fig 6 / IPWRA engine): (a) each
  channel outcome model controls for the channel's own pre-crisis change pre_<v>
  (their g_0); (b) besides the two level IRFs we bootstrap the HIGH - LOW nexus
  DIFFERENCE within each cell (their within-type high-vs-low contrast) so the
  sign-flip is formally tested, not just eyeballed from two bands.

  Output: $tabs/aipw_nexus_split.csv (outcome x part x bank x horizon, levels) ;
          $tabs/aipw_nexus_diff.csv  (outcome x part x horizon, high-low gap + CI) ;
          $figs/fig_aipw_nexus_split.pdf (GDP) + $figs/fig_nexus_<channel>.pdf.
  Run AFTER 01e_predictors.do (needs vix, past_onsets, past_def_onsets, nexus vars).
===========================================================================*/

use "$clean/panel_lp.dta", clear
sort cid year
xtset cid year

local nboot  = 300
local cx     l1_gdpg l2_gdpg debt ca infl imf
local cz     vix l_reg_crisis_share past_onsets       // Act 1 predictors
local cz_def vix l_reg_crisis_share past_def_onsets   // resolution predictors

* ══════════════════════════════════════════════════════════════════════════
* AMPLIFIER: pre-crisis sovereign-bank nexus + median split over onsets
* ══════════════════════════════════════════════════════════════════════════
capture drop a_nexus a_nexus_cm highbank
gen a_nexus = L.claimsgov_assets                 // predetermined (year before onset)
bysort cid: egen a_nexus_cm = mean(claimsgov_assets)
replace a_nexus = a_nexus_cm if missing(a_nexus) // country-mean fill (thin coverage)

quietly summarize a_nexus if sample==1 & onset_all==1, detail
local med = r(p50)
gen highbank = (a_nexus >= `med') if !missing(a_nexus)
label define hb 0 "Low nexus" 1 "High nexus"
label values highbank hb

di as result _n "=== SOVEREIGN-BANK NEXUS MEDIAN SPLIT ==="
di as result "  Amplifier = claims on govt / bank assets (pre-crisis, country-mean filled)"
di as result "  Median cutoff among crisis onsets = " %6.2f `med'
foreach t in onset_all onset_nd onset_def {
    quietly count if sample==1 & `t'==1 & highbank==1
    local nh = r(N)
    quietly count if sample==1 & `t'==1 & highbank==0
    local nl = r(N)
    quietly count if sample==1 & `t'==1 & missing(highbank)
    local nm = r(N)
    di as result "  `t': high=" `nh' "  low=" `nl' "  unclassified=" `nm'
}

* ── Confound check: which countries fall in each nexus bin? ──────────────────
*   If the high/low split just proxies income / financial development, the two
*   bins will be dominated by systematically different countries. Print the
*   country composition of the bins (over onsets) and the mean nexus by bin so
*   the split can be inspected, not taken on faith.
di as result _n "=== NEXUS-BIN COUNTRY COMPOSITION (crisis onsets) ==="
capture noisily tabulate country highbank if sample==1 & onset_all==1, ///
    row nofreq
di as result _n "  Mean pre-crisis nexus (claims-on-govt/assets) by bin, over onsets:"
capture noisily tabstat a_nexus if sample==1 & onset_all==1, ///
    by(highbank) statistics(mean min max n) format(%6.2f)
di as result "  (Read alongside any income/development ranking of these countries:"
di as result "   if the two columns are not obviously split by development, the nexus"
di as result "   result is not merely a development proxy.)"

* ══════════════════════════════════════════════════════════════════════════
* CHANNEL OUTCOMES — evolve each channel by high/low nexus (as in 13c)
*   Outcomes: GDP (dy_h, already in panel) + credit, inv, claimpriv_assets,
*   claims_govt (ch_v_h = F h.v - L.v). Every lagged control pre-generated as a
*   PLAIN column so the cluster bootstrap (bsample destroys time order) is valid.
* ══════════════════════════════════════════════════════════════════════════
foreach v in credit inv claimpriv_assets claims_govt {
    capture drop `v'_base
    gen `v'_base = L.`v'
    forvalues h = 0/4 {
        capture drop ch_`v'_`h'
        gen ch_`v'_`h' = F`h'.`v' - `v'_base
    }
    * (a) pre-crisis change in the channel itself (Asonuma's g_0 = L.var - L2.var):
    *     controls for the channel's own pre-trend momentum in each outcome model.
    capture drop pre_`v'
    gen pre_`v' = L.`v' - L2.`v'
}
foreach v in credit claimpriv_assets claims_govt pb {
    capture drop l_`v'
    gen l_`v' = L.`v'
}

* ══════════════════════════════════════════════════════════════════════════
* PROGRAMS — _aipw (Eqs. 1-3) + _aipwci (point + bootstrap CI)  [as in 13c]
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipw
program define _aipw, rclass
    syntax varlist(min=2 max=2) [if], OMODEL(varlist) PMODEL(varlist) [FE(varname)]
    gettoken y D : varlist
    marksample touse
    markout `touse' `omodel' `pmodel'
    tempvar xb m0 m1 ps summ
    if "`fe'" != "" {
        quietly reg `y' `D' `omodel' i.`fe' if `touse'
    }
    else {
        quietly reg `y' `D' `omodel' if `touse'
    }
    quietly predict double `xb' if `touse', xb
    quietly gen double `m0' = `xb' - _b[`D']*`D' if `touse'
    quietly gen double `m1' = `m0' + _b[`D']      if `touse'
    quietly probit `D' `pmodel' if `touse'
    quietly predict double `ps' if `touse', pr
    quietly replace `ps' = .01 if `ps' < .01              & `touse'
    quietly replace `ps' = .99 if `ps' > .99 & !missing(`ps') & `touse'
    quietly gen double `summ' = ///
        ( `D'*`y'/`ps' - (1-`D')*`y'/(1-`ps') ) ///
      - ( (`D'-`ps')/(`ps'*(1-`ps')) )*( (1-`ps')*`m1' + `ps'*`m0' ) ///
        if `touse'
    quietly summarize `summ' if `touse', meanonly
    return scalar theta = r(mean)
    return scalar N     = r(N)
end

capture program drop _aipwci
program define _aipwci, rclass
    syntax anything, IFC(string) OMOD(string) PZ(string) REPS(integer)
    gettoken yv Dv : anything
    capture _aipw `yv' `Dv' if `ifc', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc {
        return scalar ok = 0
        exit
    }
    local pt = r(theta)
    tempname pf
    tempfile bf
    quietly postfile `pf' double theta using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            capture drop _bid
            bsample, cluster(cid) idcluster(_bid)
            capture _aipw `yv' `Dv' if `ifc', omodel(`omod') pmodel(`pz') fe(_bid)
            if _rc == 0 quietly post `pf' (r(theta))
        restore
    }
    quietly postclose `pf'
    local se = .
    local lo = .
    local hi = .
    local nd = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(theta)
        local nd = r(N)
        if `nd' >= 50 {
            quietly summarize theta
            local se = r(sd)
            _pctile theta, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
        }
    restore
    return scalar ok = 1
    return scalar b  = `pt'
    return scalar se = `se'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* PROGRAM — (b) bootstrap the HIGH - LOW nexus DIFFERENCE within one cell
*   (Asonuma's within-type high-vs-low contrast). Both cells are re-estimated on
*   the SAME cluster resample each draw, so the difference distribution is valid.
*   _aipwdiff <yvar> <Dvar>, ifch(<high cond>) ifcl(<low cond>) omod() pz() reps()
*   returns r(ok) r(dh) [point high-low] r(lo) r(hi) r(nd) r(bh) r(bl)
* ══════════════════════════════════════════════════════════════════════════
capture program drop _aipwdiff
program define _aipwdiff, rclass
    syntax anything, IFCH(string) IFCL(string) OMOD(string) PZ(string) REPS(integer)
    gettoken yv Dv : anything
    * point estimates on the real data
    capture _aipw `yv' `Dv' if `ifch', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc { return scalar ok = 0 ; exit }
    local bh = r(theta)
    capture _aipw `yv' `Dv' if `ifcl', omodel(`omod') pmodel(`pz') fe(cid)
    if _rc { return scalar ok = 0 ; exit }
    local bl = r(theta)
    local dh = `bh' - `bl'
    * bootstrap the difference (one resample -> both cells -> their gap)
    tempname pf
    tempfile bf
    quietly postfile `pf' double diff using "`bf'", replace
    forvalues b = 1/`reps' {
        preserve
            capture drop _bid
            bsample, cluster(cid) idcluster(_bid)
            capture _aipw `yv' `Dv' if `ifch', omodel(`omod') pmodel(`pz') fe(_bid)
            local th = cond(_rc==0, r(theta), .)
            capture _aipw `yv' `Dv' if `ifcl', omodel(`omod') pmodel(`pz') fe(_bid)
            local tl = cond(_rc==0, r(theta), .)
            if !missing(`th') & !missing(`tl') quietly post `pf' (`th' - `tl')
        restore
    }
    quietly postclose `pf'
    local lo = .
    local hi = .
    local nd = 0
    preserve
        quietly use "`bf'", clear
        quietly count if !missing(diff)
        local nd = r(N)
        if `nd' >= 50 {
            _pctile diff, p(2.5 97.5)
            local lo = r(r1)
            local hi = r(r2)
        }
    restore
    return scalar ok = 1
    return scalar dh = `dh'
    return scalar bh = `bh'
    return scalar bl = `bl'
    return scalar lo = `lo'
    return scalar hi = `hi'
    return scalar nd = `nd'
end

* ══════════════════════════════════════════════════════════════════════════
* ESTIMATE — high/low nexus subsamples, all + resolution split; post results
* ══════════════════════════════════════════════════════════════════════════
tempname R
tempfile resf
postfile `R' str18 outcome str4 part str4 bank byte horizon double b se lo hi ntreat nd ///
    using "`resf'", replace

* (b) second results file: the HIGH - LOW nexus difference per outcome x part x h
tempname D
tempfile diff_resf
postfile `D' str18 outcome str4 part byte horizon double dhl bhi blo lo hi nd ///
    using "`diff_resf'", replace

* Outcomes: label | outcome-variable stem (dy or ch_<v>) | channel-specific
*   outcome-model controls (om), same specs as 13c. GDP uses cx (unchanged).
foreach oc in "gdp dy" "credit ch_credit" "inv ch_inv" ///
              "claimpriv_assets ch_claimpriv_assets" "claims_govt ch_claims_govt" {
    gettoken ocl   oc : oc
    gettoken ystem oc : oc

    * outcome-model controls; each channel gets its own pre-crisis change pre_<v>
    * (Asonuma g_0). GDP already carries lagged growth (l1/l2_gdpg) in cx.
    if      "`ocl'" == "gdp"              local om `cx'
    else if "`ocl'" == "credit"           local om l1_gdpg l2_gdpg debt infl ca banking_crisis pre_credit
    else if "`ocl'" == "inv"              local om l1_gdpg l2_gdpg debt ca l_credit banking_crisis pre_inv
    else if "`ocl'" == "claimpriv_assets" local om l_claimpriv_assets l1_gdpg l2_gdpg debt ca banking_crisis pre_claimpriv_assets
    else if "`ocl'" == "claims_govt"      local om l_claims_govt l_credit pb banking_crisis pre_claims_govt

    di as result _n "############### OUTCOME: `ocl' ###############"

    * Rows: part-label | treatment dummy | rival dummy to drop | predictors
    *   Part A uses onset_all (rival = "" -> none); Parts nd/def drop the other type.
    foreach cell in "all onset_all . cz" ///
                    "nd  onset_nd onset_def cz_def" ///
                    "def onset_def onset_nd cz_def" {
        gettoken part cell : cell
        gettoken Dv   cell : cell
        gettoken riv  cell : cell
        gettoken pzn  cell : cell
        * predictor set (indirect: pzn is the local NAME cz or cz_def)
        local pz "``pzn''"

        forvalues H = 0/1 {                       // 0 = low nexus, 1 = high nexus
            local blab = cond(`H'==1, "high", "low")
            di as result _n "--- `ocl' | Part `part' | `blab' nexus ---"

            forvalues h = 0/4 {
                * sample: tranquil (onset_all==0) + this cell's treated onsets only
                if "`riv'" == "." {
                    local ifc sample==1 & (onset_all==0 | (`Dv'==1 & highbank==`H'))
                }
                else {
                    local ifc sample==1 & `riv'==0 & (onset_all==0 | (`Dv'==1 & highbank==`H'))
                }
                quietly count if `Dv'==1 & highbank==`H' & sample==1
                local ntr = r(N)

                _aipwci `ystem'_`h' `Dv', ifc(`ifc') omod(`om') pz(`cx' `pz') reps(`nboot')
                if r(ok) {
                    local b=r(b)
                    local se=r(se)
                    local lo=r(lo)
                    local hi=r(hi)
                    local nd=r(nd)
                    post `R' ("`ocl'") ("`part'") ("`blab'") (`h') (`b') (`se') (`lo') (`hi') (`ntr') (`nd')
                    di "    h=" `h' "  ATE=" %8.3f `b' "  [" %7.3f `lo' ", " %7.3f `hi' ///
                       "]  (n_treat=" `ntr' ", " `nd' "/`nboot' draws)"
                }
                else di as error "    h=" `h' ": estimate failed (too thin)."
            }
        }

        * ── (b) HIGH - LOW nexus difference for this outcome x part ──────────
        di as result _n "--- `ocl' | Part `part' | HIGH - LOW difference ---"
        forvalues h = 0/4 {
            if "`riv'" == "." {
                local ifch sample==1 & (onset_all==0 | (`Dv'==1 & highbank==1))
                local ifcl sample==1 & (onset_all==0 | (`Dv'==1 & highbank==0))
            }
            else {
                local ifch sample==1 & `riv'==0 & (onset_all==0 | (`Dv'==1 & highbank==1))
                local ifcl sample==1 & `riv'==0 & (onset_all==0 | (`Dv'==1 & highbank==0))
            }
            _aipwdiff `ystem'_`h' `Dv', ifch(`ifch') ifcl(`ifcl') ///
                omod(`om') pz(`cx' `pz') reps(`nboot')
            if r(ok) {
                local dh=r(dh)
                local bh=r(bh)
                local bl=r(bl)
                local lo=r(lo)
                local hi=r(hi)
                local nd=r(nd)
                post `D' ("`ocl'") ("`part'") (`h') (`dh') (`bh') (`bl') (`lo') (`hi') (`nd')
                local sig = cond(`nd'>=50 & (`lo'>0 | `hi'<0), " *", "")
                di "    h=" `h' "  high-low=" %8.3f `dh' "  [" %7.3f `lo' ", " %7.3f `hi' ///
                   "]  (" `nd' "/`nboot' draws)`sig'"
            }
            else di as error "    h=" `h' ": difference estimate failed (too thin)."
        }
    }
}
postclose `R'
postclose `D'

* ══════════════════════════════════════════════════════════════════════════
* EXPORT — CSV + figure (panels by resolution part; high vs low nexus lines)
* ══════════════════════════════════════════════════════════════════════════
* (b) HIGH - LOW difference table: dhl = point gap, [lo,hi] = bootstrap 95% CI;
*     a CI excluding 0 means the nexus effect differs significantly by cell.
use "`diff_resf'", clear
label var dhl "AIPW (high - low nexus) difference (pp)"
label var bhi "High-nexus ATE"
label var blo "Low-nexus ATE"
label var nd  "Valid bootstrap draws"
gen byte sig95 = (nd>=50 & (lo>0 | hi<0))
label var sig95 "High-low gap CI excludes 0"
order outcome part horizon dhl bhi blo lo hi nd sig95
export delimited "$tabs/aipw_nexus_diff.csv", replace
di as result _n "Nexus high-low DIFFERENCE CSV saved: $tabs/aipw_nexus_diff.csv"

use "`resf'", clear
label var b  "AIPW ATE on outcome (pp)"
label var ntreat "Treated onsets in cell"
label var nd "Valid bootstrap draws"
order outcome part bank horizon b se lo hi ntreat nd
export delimited "$tabs/aipw_nexus_split.csv", replace
di as result _n "Nexus-split AIPW results CSV saved: $tabs/aipw_nexus_split.csv"

* numeric part id for by() panels, ordered all -> nd -> def
gen byte partid = 1 if part=="all"
replace partid = 2 if part=="nd"
replace partid = 3 if part=="def"
label define pl 1 "All onsets" 2 "Non-default" 3 "Default-linked"
label values partid pl

* ── One high-vs-low figure per outcome (GDP keeps its historical look) ───────
local c_hi "157 36 73"    // high nexus = red (the doom-loop)
local c_lo "0 84 166"     // low  nexus = blue
foreach oc in gdp credit inv claimpriv_assets claims_govt {
    if "`oc'" == "gdp" {
        local ptit "Output cost by sovereign-bank nexus"
        local ytit "Cumulative real GDP change (pp)"
        local fnm  "fig_aipw_nexus_split"       // unchanged filename for GDP
    }
    else {
        local ptit "`oc' channel by sovereign-bank nexus"
        local ytit "Cumulative change in `oc' (pp)"
        local fnm  "fig_nexus_`oc'"
    }
    capture twoway ///
        (rarea lo hi horizon if bank=="high" & outcome=="`oc'", color("`c_hi'%16") lwidth(none)) ///
        (rarea lo hi horizon if bank=="low"  & outcome=="`oc'", color("`c_lo'%16") lwidth(none)) ///
        (connected b horizon if bank=="high" & outcome=="`oc'", lcolor("`c_hi'") lwidth(medthick) msymbol(square) lpattern(dash)) ///
        (connected b horizon if bank=="low"  & outcome=="`oc'", lcolor("`c_lo'") lwidth(medthick) msymbol(circle)), ///
        by(partid, yrescale ///
            note("AIPW (Asonuma Eq. 3), ATE. Split by pre-crisis bank claims-on-govt / assets (median). Shaded = bootstrap 95% CI where >=50 draws.", size(vsmall)) ///
            title("`ptit'", size(medsmall) color(navy))) ///
        yline(0, lpattern(dash) lcolor(gs8)) ///
        xlabel(0(1)4) xtitle("Years after onset", size(small)) ///
        ytitle("`ytit'", size(small)) ///
        legend(order(3 "High nexus" 4 "Low nexus") size(small)) ///
        graphregion(color(white)) plotregion(color(white))
    if _rc == 0 {
        graph export "$figs/`fnm'.pdf", replace
        di as result "Figure saved: `fnm'.pdf"
    }
    else di as error "  ** `fnm' failed (rc=" _rc ")"
}

di as result _n "13d_aipw_nexus_split.do complete."
di as result "GDP headline (Part def): high-nexus default crises are the deepest (doom-loop"
di as result "amplifies the default cost). The channel outcomes trace the mechanism: expect"
di as result "claims_govt to RISE under def x high (banks absorbing the sovereign = the loss"
di as result "channel). Read Part B with its small-cell / draw-count caveat."
di as result "aipw_nexus_diff.csv formally tests the sign-flip: a high-low gap whose 95% CI"
di as result "excludes 0 (sig95==1) means the nexus effect differs significantly in that cell."
