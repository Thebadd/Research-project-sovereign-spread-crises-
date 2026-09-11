/*===========================================================================
  13E_NEXUS_MECHANISM_DIAGRAM.DO
  Conceptual 2x2 schematic of the sovereign-bank-nexus mechanism (Fig for the
  paper). NOT an estimation file — it hard-codes the qualitative pattern that
  13d produced so the reader can see, at a glance, why the GDP cost flips sign
  across the resolution x nexus grid:

      Columns = pre-crisis sovereign-bank nexus (Low | High)
      Rows    = crisis resolution (Non-default top | Default-linked bottom)

    nd x low  : costly (-5 to -9.5)  — credit/investment squeeze, external roll
    nd x high : CUSHION (~0)         — captive domestic buffer, banks hold calmly
    def x low : costly early, fades  — milder bank reallocation
    def x high: AMPLIFIER (-6 to -11)— doom loop: claims_govt UP, priv lending DOWN

  Numbers are the 13d AIPW ATE ranges on GDPpc; the arrow in the def x high cell
  depicts the within-bank reallocation (money moves firms -> government).

  Output: $figs/fig_nexus_mechanism_2x2.pdf.  Run any time after paths are set
  (no data dependency — draws on a blank canvas).
===========================================================================*/

* Blank 0-10 x 0-10 canvas (two invisible corner points fix the plot region).
clear
set obs 2
gen x = .
gen y = .
replace x = 0 in 1
replace y = 0 in 1
replace x = 10 in 2
replace y = 10 in 2

local c_cush "0 110 70"     // cushion  = green
local c_ampl "157 36 73"    // amplifier = red
local c_cost "60 60 60"     // ordinary cost = grey

capture twoway ///
    (scatter y x, msymbol(none)) ///
    (pcarrowi 7.5 8.6 7.5 6.4, lcolor("`c_ampl'") mcolor("`c_ampl'") lwidth(medthick)) ///
    , ///
    xline(5, lcolor(gs10) lwidth(medium)) ///
    yline(5, lcolor(gs10) lwidth(medium)) ///
    xscale(range(0 10) off) yscale(range(0 10) off) ///
    xlabel(none) ylabel(none) ///
    title("The sovereign-bank nexus: cushion vs. amplifier", size(medium) color(navy)) ///
    subtitle("AIPW output cost (GDPpc) and bank-portfolio response, by crisis resolution x pre-crisis nexus", ///
        size(vsmall)) ///
    note("Cells labelled with the 13d AIPW ATE range on GDPpc across h=0..4. Nexus = pre-crisis bank claims on govt / assets (median split)." ///
         "Def x High is the doom-loop cell: banks raise sovereign claims and cut private lending (arrow); Non-default x High is inert (captive buffer).", ///
        size(vsmall)) ///
    legend(off) graphregion(color(white)) plotregion(color(white) margin(zero)) ///
    text(10.3 2.5 "LOW nexus", size(small) color(navy) place(c)) ///
    text(10.3 7.5 "HIGH nexus", size(small) color(navy) place(c)) ///
    text(7.5 -0.6 "Non-default", size(small) color(navy) orient(vertical) place(c)) ///
    text(2.5 -0.6 "Default-linked", size(small) color(navy) orient(vertical) place(c)) ///
    /// ── Top-left: nd x low  (ordinary cost) ──
    text(4.55 2.5 "Costly", size(medsmall) color("`c_cost'") place(c)) ///
    text(3.9 2.5 "GDP -5 to -9.5" "credit & investment squeeze;" "rollover in external markets", ///
        size(vsmall) color("`c_cost'") place(c)) ///
    /// ── Top-right: nd x high (cushion) ──
    text(4.55 7.5 "CUSHION  (~0)", size(medsmall) color("`c_cush'") place(c)) ///
    text(3.9 7.5 "GDP ~ 0, insignificant" "banks calmly hold govt debt" "(captive funding buffer)", ///
        size(vsmall) color("`c_cush'") place(c)) ///
    /// ── Bottom-left: def x low (cost, fades) ──
    text(9.55 2.5 "Costly, then fades", size(medsmall) color("`c_cost'") place(c)) ///
    text(8.9 2.5 "GDP -6 -> -2 (loses sig.)" "milder bank reallocation", ///
        size(vsmall) color("`c_cost'") place(c)) ///
    /// ── Bottom-right: def x high (amplifier) ──
    text(9.55 7.5 "AMPLIFIER", size(medsmall) color("`c_ampl'") place(c)) ///
    text(9.0 7.5 "GDP -6 -> -11 (deepest)", size(vsmall) color("`c_ampl'") place(c)) ///
    text(6.9 6.4 "Gov claims UP" "(+9 to +12)", size(vsmall) color("`c_ampl'") place(c)) ///
    text(6.9 8.7 "Priv lending DOWN" "(-5 to -6)", size(vsmall) color("`c_ampl'") place(c)) ///
    text(8.15 7.5 "money reallocates firms -> govt", size(vsmall) color("`c_ampl'") place(c))

if _rc == 0 {
    graph export "$figs/fig_nexus_mechanism_2x2.pdf", replace
    di as result "Figure saved: fig_nexus_mechanism_2x2.pdf"
}
else di as error "  ** fig_nexus_mechanism_2x2 failed (rc=" _rc ")"

di as result _n "13e_nexus_mechanism_diagram.do complete."
