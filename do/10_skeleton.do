/*===========================================================================
  10_SKELETON.DO   —  FROM-SCRATCH REBUILD, STAGE 1
  Panel skeleton built ENTIRELY from David's own spread-crisis database
  (EM_Spread_Crisis_DB_FINAL.xlsx). This is the only "kept" data — the crisis
  definition (onsets, spreads, classification) is his own work. Every MACRO
  variable is rebuilt from official sources in stages 11+ and merged on iso3.

  Produces the country-year skeleton with a standard ISO3 key:
    iso3 country region year  |  spr_max spr_mean crisis_any onset_all
    continuation ep_id ep_status  |  classification nondefault onset_nd onset_def
    sample_base
  Saves: $clean/panel_build.dta   (macro merged on top of this in 11+)

  Source sheets (README of the workbook documents JP Morgan EMBIG Global):
    Panel_Annual    — full country-year panel, onset dating (D_it), spreads.
    Episode_Summary — 61 episodes, default vs non-default Classification.
===========================================================================*/

* ── 1. Import the annual crisis panel (row 1 = title, row 2 = header) ────────
import excel "$raw/EM_Spread_Crisis_DB_FINAL.xlsx", ///
    sheet("Panel_Annual") cellrange(A3) clear
* No firstrow: columns import as A..M. Rename by position (stable layout).
rename (A B C D E F G H I J K L M) ///
       (country region year spr_max spr_mean crit1 crit2 ///
        crisis_any ep_id ep_status onset_all continuation dq_note)

* numeric coercion (import may type some as string if blanks present)
foreach v in year spr_max spr_mean crit1 crit2 crisis_any onset_all continuation {
    capture confirm numeric variable `v'
    if _rc capture destring `v', replace force
}
drop dq_note
label var crit1 "Criterion 1 flag (spread > 1000bps)"
label var crit2 "Criterion 2 flag (99pct HRT)"
drop if missing(country) | missing(year)

* ── 2. ISO3 crosswalk (52 crisis countries → ISO3, the universal merge key) ──
gen str3 iso3 = ""
replace iso3="ARG" if country=="Argentina"
replace iso3="ARM" if country=="Armenia"
replace iso3="AZE" if country=="Azerbaijan"
replace iso3="BLZ" if country=="Belize"
replace iso3="BOL" if country=="Bolivia"
replace iso3="BRA" if country=="Brazil"
replace iso3="BGR" if country=="Bulgaria"
replace iso3="CHL" if country=="Chile"
replace iso3="CHN" if country=="China"
replace iso3="COL" if country=="Colombia"
replace iso3="CRI" if country=="Costa Rica"
replace iso3="CIV" if country=="Cote d'Ivoire"
replace iso3="DOM" if country=="Dominican Republic"
replace iso3="ECU" if country=="Ecuador"
replace iso3="EGY" if country=="Egypt"
replace iso3="SLV" if country=="El Salvador"
replace iso3="GHA" if country=="Ghana"
replace iso3="GTM" if country=="Guatemala"
replace iso3="HND" if country=="Honduras"
replace iso3="HUN" if country=="Hungary"
replace iso3="IND" if country=="India"
replace iso3="IDN" if country=="Indonesia"
replace iso3="JAM" if country=="Jamaica"
replace iso3="JOR" if country=="Jordan"
replace iso3="KAZ" if country=="Kazakhstan"
replace iso3="KEN" if country=="Kenya"
replace iso3="LBN" if country=="Lebanon"
replace iso3="MYS" if country=="Malaysia"
replace iso3="MEX" if country=="Mexico"
replace iso3="MAR" if country=="Morocco"
replace iso3="NAM" if country=="Namibia"
replace iso3="NGA" if country=="Nigeria"
replace iso3="PAK" if country=="Pakistan"
replace iso3="PAN" if country=="Panama"
replace iso3="PRY" if country=="Paraguay"
replace iso3="PER" if country=="Peru"
replace iso3="PHL" if country=="Philippines"
replace iso3="POL" if country=="Poland"
replace iso3="ROU" if country=="Romania"
replace iso3="RUS" if country=="Russia"
replace iso3="SEN" if country=="Senegal"
replace iso3="SRB" if country=="Serbia"
replace iso3="ZAF" if country=="South Africa"
replace iso3="LKA" if country=="Sri Lanka"
replace iso3="TTO" if country=="Trinidad and Tobago"
replace iso3="TUN" if country=="Tunisia"
replace iso3="TUR" if country=="Turkey"
replace iso3="UKR" if country=="Ukraine"
replace iso3="URY" if country=="Uruguay"
replace iso3="VEN" if country=="Venezuela"
replace iso3="VNM" if country=="Vietnam"
replace iso3="ZMB" if country=="Zambia"

quietly count if iso3==""
if r(N) > 0 {
    di as error "  ** `r(N)' rows with unmapped country — extend the crosswalk:"
    tabulate country if iso3=="", missing
}
label var iso3 "ISO3 country code (universal merge key)"

* ── 3. Merge default/non-default classification from Episode_Summary ─────────
preserve
    import excel "$raw/EM_Spread_Crisis_DB_FINAL.xlsx", ///
        sheet("Episode_Summary") cellrange(A3) clear
    * A=Country C=Onset Year H=Classification
    keep A C H
    rename (A C H) (country onset_year classification)
    capture destring onset_year, replace force
    drop if missing(country) | missing(onset_year)
    duplicates drop country onset_year, force
    rename onset_year year
    tempfile classif
    save `classif'
restore
merge m:1 country year using `classif', keep(master match) nogen

* ── 4. Build resolution dummies (mirrors the original 01_build logic) ────────
gen byte nondefault = .
replace nondefault = 1 if strpos(classification,"NON-DEFAULT") > 0
replace nondefault = 0 if !missing(classification) & strpos(classification,"NON-DEFAULT")==0
* Venezuela 2008: restructuring began 2017 (long lag) — treated as non-default onset.
replace nondefault = 1 if country=="Venezuela" & year==2008 & onset_all==1

gen byte onset_nd  = (onset_all==1 & nondefault==1)
replace  onset_nd  = 0 if missing(onset_nd)
gen byte onset_def = (onset_all==1 & nondefault==0)
replace  onset_def = 0 if missing(onset_def)

label var onset_all "Spread-crisis onset (all episodes)"
label var onset_nd  "Onset — non-default episodes"
label var onset_def "Onset — default-linked episodes"
label var nondefault "1 = non-default episode onset"

* ── 4b. LEGACY carry-over: IMF-program dummy ────────────────────────────────
*   Kept as-is from the prior dataset by decision — 98% coverage, gap-free at all
*   onsets, but its ORIGINAL SOURCE could not be recovered (undocumented upstream;
*   likely IMF MONA / History of Lending Arrangements, and it appears to include
*   precautionary FCLs e.g. Chile/Colombia/Mexico). This is the ONLY variable
*   carried from the old FINAL_DATASET_STATA_V2.csv; everything else is rebuilt
*   from official sources. Flagged accordingly in DATA_SOURCES.md.
preserve
    import delimited "$raw/FINAL_DATASET_STATA_V2.csv", varnames(1) clear
    keep country year imf_program
    rename imf_program imf
    capture destring year imf, replace force
    drop if missing(country) | missing(year)
    duplicates drop country year, force
    tempfile imfleg
    save `imfleg'
restore
merge m:1 country year using `imfleg', keep(master match) nogen
label var imf "IMF-program dummy (LEGACY carry-over; original source unrecovered)"

* ── 4c. CARRY-IN ROWS: lag scaffolding at each country's left edge ──────────
*   The panel is UNBALANCED on the left: it starts when a country enters the EMBIG
*   index, not when its macro data begins (Armenia ~2013, Kenya and Zambia ~2014).
*   Every downstream merge is keep(master match) with this skeleton as master, so
*   WEO's 1980-2031 series is thrown away for every pre-entry year. The lag
*   operators in 18 then fall off the start of each country's panel: a 2013 entrant
*   with a 2013 onset has no l1_gdpg (which needs GDP at t-1 AND t-2), which is why
*   l1_gdpg covers only 42/61 onsets while gdp_real covers 59/61. Likewise
*   l_banking_crisis covers 52/61 even though banking_crisis is 61/61 — those 9 are
*   lost purely because L. has no previous row, not because anything is missing.
*
*   Fix: add 3 rows immediately BEFORE each country's first panel year (3 covers the
*   deepest lag in use — L3.ln_gdp in dy_m2, L2.gdpg in l2_gdpg). Stages 11-16 merge
*   on iso3 x year / country x year, so they populate these rows automatically with
*   no per-source edits, and all transform logic stays in stage 18.
*
*   These rows exist ONLY to feed L. They are flagged carryin==1 and excluded from
*   sample_base here and from `sample' in 18, so they never enter an estimation.
*   They carry no EMBIG quote, so we assert NOTHING about whether those years were
*   tranquil — spr_max/spr_mean stay missing by design.
gen byte carryin = 0
preserve
    bysort iso3 (year): keep if _n == 1
    keep iso3 country region year
    rename year first_yr
    expand 3
    bysort iso3: gen int year = first_yr - _n     // first_yr-1, -2, -3
    drop first_yr
    gen byte carryin = 1
    tempfile preroll
    save `preroll'
restore
append using `preroll'

* Crisis variables are known-zero on carry-in rows (no episode can have started in a
* year the crisis DB does not cover). Spreads deliberately left missing.
foreach v in onset_all onset_nd onset_def crisis_any continuation crit1 crit2 {
    capture replace `v' = 0 if carryin == 1
}
label var carryin "1 = pre-entry scaffolding row (lag support only; never in sample)"

* ── 5. Preliminary estimation-sample flag (finalised in 18 once GDP merged) ──
gen byte sample_base = (continuation==0) & carryin==0
label var sample_base "Onset + tranquil years (excl. continuation & carry-in); GDP-availability added in 18"

* ── 6. Save skeleton ────────────────────────────────────────────────────────
order iso3 country region year spr_max spr_mean crit1 crit2 crisis_any ///
      onset_all onset_nd onset_def nondefault continuation ep_id ep_status ///
      classification carryin sample_base
sort iso3 year
* panel_build.dta = the growing panel; each stage (11+) merges macro on top of it.
save "$clean/panel_build.dta", replace

di as result _n "10_skeleton.do complete."
quietly levelsof iso3, local(cc)
di as result "  Countries: " `: word count `cc''  "   Rows: " _N
quietly count if carryin==1
di as result "    of which carry-in (lag scaffolding, never in sample): `r(N)'"
quietly summarize year if carryin==0
di as result "    EMBIG panel span: " r(min) "-" r(max) _c
quietly summarize year
di as result "   (with carry-in: " r(min) "-" r(max) ")"
quietly count if onset_all==1
di as result "  Onsets: all=" r(N) _c
quietly count if onset_nd==1
di as result "   nd=" r(N) _c
quietly count if onset_def==1
di as result "   def=" r(N)
