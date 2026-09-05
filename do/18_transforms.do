/*===========================================================================
  18_TRANSFORMS.DO   —  FROM-SCRATCH REBUILD, STAGE 8  (finalise panel_lp.dta)
  Build all econometric transforms IN CODE from the sourced levels, then save
  the analysis file $clean/panel_lp.dta consumed by 02..16 unchanged.

    gdpg              real GDP growth, % = 100*(ln gdp_real - ln L.gdp_real)
    l1_gdpg, l2_gdpg  first/second lag of gdpg (control momentum, = cx)
    ln_gdp_base       L.ln_gdp (headline LP baseline at t-1; ln_gdppc_base = per-capita robustness)
    dy_0..dy_4        cumulative % change in log TOTAL real GDP (ln gdp_real), F h vs t-1 (headline outcome)
    dy_pc_0..dy_pc_4  per-capita version (robustness)
    dy_m2             pre-trend placebo: h=-2 on the same t-1 base as dy_h (h=-1 is
                      trivially 0 under that base and is not estimated)
    l_hyperinfl (L.infl>50, core), l_lninfl (= ln(1+L.infl/100)) and l_infl (raw)
    built as continuous robustness alternatives, not in the core,
    l_govexp, l_open, l_credit (+l_credit_bank robustness), $ctrl_core  common-core ingredients
    tot_chg, exchange2, ex_dum1-5  robustness-tier controls (Asonuma additional)
    l_spr_mean/max    lagged EMBIG spread (balance table)
    sample            onset + tranquil years, excl. continuation, GDP base present
===========================================================================*/

use "$clean/panel_build.dta", clear
capture xtset cid year
sort cid year

* ── Real GDP growth + its lags (Asonuma gdpg2 analog; controls cx) ──────────
capture drop gdpg l1_gdpg l2_gdpg
gen double gdpg = 100*(ln(gdp_real) - ln(L.gdp_real)) ///
    if gdp_real > 0 & L.gdp_real > 0 & !missing(gdp_real, L.gdp_real)
label var gdpg "Real GDP growth, % (WEO NGDP_R, log-difference)"
gen double l1_gdpg = L.gdpg
gen double l2_gdpg = L2.gdpg
label var l1_gdpg "L1 real GDP growth"
label var l2_gdpg "L2 real GDP growth"

* ── HEADLINE LP outcome: cumulative % change in log TOTAL real GDP ──────────
*   Aligned with Asonuma et al. (their g_h = F h.ln(gdp_real) - L.ln(gdp_real)),
*   i.e. cumulative real-GDP growth relative to t-1 (NOT per capita).
capture drop ln_gdp ln_gdp_base dy_0 dy_1 dy_2 dy_3 dy_4 dy_m2
gen double ln_gdp = ln(gdp_real) if gdp_real > 0 & !missing(gdp_real)
gen double ln_gdp_base = L.ln_gdp
label var ln_gdp_base "log real GDP at t-1 (LP baseline)"
forvalues h = 0/4 {
    gen double dy_`h' = (F`h'.ln_gdp - ln_gdp_base) * 100
    label var dy_`h' "Cum. % change in log real GDP (total): F`h' vs t-1"
}
* ── Pre-trend placebo on the same (total real GDP) series, SAME BASE as dy_h ──
* dy_h = F h.ln_gdp - L.ln_gdp for h=0..4, always against the t-1 base. Plugging
* h=-1 into that identical formula gives L.ln_gdp - L.ln_gdp = 0 for every row —
* the base year trivially equals itself, so h=-1 is not an estimable placebo, it
* IS the anchor (matches the standard event-study convention, e.g. Ugarte-Ruiz
* 2025 WP 25-09 sec. 2.8/5.8, where h=-1 is normalized to zero, not estimated).
* The first REAL pre-crisis horizon under this same base is h=-2:
*   dy_m2 = F(-2).ln_gdp - L.ln_gdp = L2.ln_gdp - L.ln_gdp
* i.e. the change from t-2 to the t-1 base — never touches the onset year t.
* dy_m2 is algebraically -l1_gdpg (both are +/-100*(L.ln_gdp - L2.ln_gdp)), so the
* pre-trend regressions in 02/03 must drop l1_gdpg from the RHS — see the
* `controls_pre' local there. Controlling for the placebo outcome would guarantee
* a null test.
gen double dy_m2 = (L2.ln_gdp - L.ln_gdp) * 100
label var dy_m2 "Pre-trend h=-2 (same base as dy_h): GDP(t-2) - GDP(t-1)"

* ── ROBUSTNESS outcome: per-capita version (kept as dy_pc_*) ─────────────────
capture drop ln_gdppc_base dy_pc_0 dy_pc_1 dy_pc_2 dy_pc_3 dy_pc_4
gen double ln_gdppc_base = L.ln_gdppc
label var ln_gdppc_base "log real GDPpc at t-1 (per-capita LP baseline)"
forvalues h = 0/4 {
    gen double dy_pc_`h' = (F`h'.ln_gdppc - ln_gdppc_base) * 100
    label var dy_pc_`h' "Cum. % change in log real GDPpc: F`h' vs t-1 (robustness)"
}

* ── Lagged spreads (used in the balance table) ──────────────────────────────
capture drop l_spr_mean l_spr_max
gen double l_spr_mean = L.spr_mean
gen double l_spr_max  = L.spr_max
label var l_spr_mean "L1 EMBIG mean spread (bps)"
label var l_spr_max  "L1 EMBIG max spread (bps)"

* ── Asonuma-aligned COMMON-CORE control set (plain predetermined columns) ───
*   One control set used as the OUTCOME-model controls in every GDP + channel
*   regression, mirroring Asonuma's $convar: GDP momentum, gov spending, openness,
*   bank-credit depth, hyperinflation dummy, banking-crisis dummy + our debt/ca.
*   Plain lagged columns so the SAME set works in xtscc LP and bsample AIPW.
* Guard: these sources are built upstream — infl/govexp in stage 11 (WEO),
* open/credit_bank in stage 12 (WDI). If any is absent, panel_build.dta was not
* built by the full 10->18 chain (e.g. 18 was run alone on a stale file, or a
* WDI raw file was missing). Stop with a clear message instead of a cryptic halt.
foreach req in infl govexp open credit_bank {
    capture confirm variable `req'
    if _rc {
        di as error "  ** 18_transforms: required variable '`req'' is not in panel_build.dta."
        di as error "     Built upstream: infl/govexp -> stage 11 (11_weo); open/credit_bank -> stage 12 (12_wdi)."
        di as error "     Fix: run the FULL chain in order (do 00_master.do), with all data/raw files present."
        exit 111
    }
}
capture drop l_hyperinfl l_infl l_lninfl l_govexp l_open l_credit_bank l_credit l_debt l_ca l_banking_crisis l_banking_duration l_banking_duration_total
* INFLATION IN THE CORE = l_hyperinfl (David's decision), matching Asonuma et al.
* LAGGED like every other core control: built from L.infl, so it flags hyperinflation at
* t-1 and is predetermined relative to the onset year (the l_ prefix now says so).
* The `if !missing(L.infl)' guard matters: Stata treats missing as +infinity, so without
* it every missing-inflation row would be coded as a hyperinflation year.
* Raw inflation cannot be used directly here: across the 52 panel countries 1994-2026 the
* median is 9.0% while Venezuela 2018 is 65,374% (9 obs above 1000%), so a raw term would
* let a handful of country-years set the coefficient. The dummy sidesteps that entirely by
* flagging the hyperinflation regime instead of its magnitude — 63 country-years exceed the
* 50% threshold (Venezuela 16, Turkey 11, Bulgaria/Argentina/Lebanon 4 each, 17 countries
* in all), so the flag is well populated in the raw data.
*   NB it collapsed to a single estimation-sample observation in the pre-2026 runs, but
*   that was the OLD banking control truncating the panel at 2018, which killed the entire
*   recent cluster (Argentina 2019-24, Lebanon 2020-23, Turkey 2022-24, Venezuela 2019-26).
*   With the Laeven-Valencia rebuild extending the sample past 2018 those years return.
*   Worth re-checking the probit for "predicts failure perfectly" on the new sample.
gen byte   l_hyperinfl    = (L.infl > 50) if !missing(L.infl)
* Continuous alternatives, built but NOT in the core — available for a dummy-vs-continuous
* robustness spec, and a one-token swap in the 14 $ctrl_core definitions if wanted. Use
* l_lninfl rather than l_infl for that: the log tames the tail (Venezuela ~75x the median
* instead of ~7,000x). Safe because the minimum is -8.53% (Azerbaijan 1999), well above
* the -100% where ln(1+x) breaks.
gen double l_infl         = L.infl
gen double l_lninfl       = ln(1 + L.infl/100) if !missing(L.infl) & L.infl > -100
gen double l_govexp       = L.govexp
gen double l_open         = L.open
gen double l_credit_bank  = L.credit_bank
gen double l_credit       = L.credit
* Exchange rate — robustness-tier, not in $ctrl_core. reer_chg (12_wdi.do) is
* already a % change, so this is a plain lag, same pattern as every other
* control here.
capture confirm variable reer_chg
if !_rc {
    gen double l_reer_chg = L.reer_chg
    label var l_reer_chg "L1 real effective exchange rate, % change (predetermined; robustness alt., not in core)"
}
else di as error "  ** reer_chg not found — skipping l_reer_chg (add REER_INDEX.xlsx to build it)."
* debt, current account, and the banking-crisis flag lagged to t-1 so the ENTIRE
* common core is predetermined (t-1) — internally coherent and matching the
* reference paper's pre-crisis controls (avoids onset-year simultaneity/bad-control).
gen double l_debt         = L.debt
gen double l_ca           = L.ca
gen byte   l_banking_crisis = L.banking_crisis
* DURATION, not a flag — the reference paper's banking_duration_lv2018 analog.
* banking_duration counts how long the crisis has ALREADY lasted at that year
* (1,2,3,...), so L. of it is genuinely t-1 information: at t you know how long it
* had run as of t-1. See 16_banking.do step 5 for why years-so-far beats total
* episode length here. The dummy and the total-length version stay built as
* robustness alternatives, each one token from the core.
gen int    l_banking_duration       = L.banking_duration
gen int    l_banking_duration_total = L.banking_duration_total
* Global-push predictors lagged to t-1 (predetermined), matching the reference
* paper's federal_funds2 = L.federal_funds. Plain saved columns => bootstrap-safe
* in the AIPW (no L. operator at estimation time). Rates are pure time series, so
* L. under xtset cid year gives the year-(t-1) value for every country. vix/ust10y
* built too so a future push-variable swap is one line. Raw fedfunds/vix/ust10y
* columns are kept (coverage report below still uses them).
gen double l_fedfunds     = L.fedfunds
gen double l_vix          = L.vix
gen double l_ust10y       = L.ust10y
label var l_fedfunds "L1 US fed funds rate (predetermined; = Asonuma federal_funds2)"
label var l_vix      "L1 CBOE VIX (predetermined)"
label var l_ust10y   "L1 US 10y Treasury yield (predetermined)"
label var l_hyperinfl    "L1 hyperinflation dummy (L.infl > 50; predetermined) — common core"
label var l_infl         "L1 CPI inflation, % (raw; robustness alternative, not in core)"
label var l_lninfl       "L1 log gross inflation = ln(1+L.infl/100) (robustness alt., not in core)"
label var l_govexp       "L1 govt expenditure, % GDP"
label var l_open         "L1 trade openness, % GDP"
label var l_credit_bank  "L1 bank credit to private / GDP (financial depth, by-banks; COMMON CORE)"
label var l_credit       "L1 private credit / GDP (all fin. corps; robustness alt., NOT in core)"
label var l_debt         "L1 public debt, % GDP (predetermined)"
label var l_ca           "L1 current account, % GDP (predetermined)"
label var l_banking_crisis      "L1 systemic banking-crisis dummy (robustness alt., not in core)"
label var l_banking_duration    "L1 years the banking crisis had already lasted (predetermined) - COMMON CORE"
label var l_banking_duration_total "L1 total banking-crisis length (robustness alt., not in core)"

* $ctrl_core's real, asserted definition lives further down this file, after
* exchange2 (one of its terms) is built and confirmed to exist -- see
* "CORE CONTROL SET" below. Nothing between here and there depends on
* $ctrl_core being set yet.

* ── ROBUSTNESS-tier controls (Asonuma additional controls; NOT in the core) ──
*   terms-of-trade change and nominal-FX-change quantile dummies, built the
*   Asonuma way. Predetermined. Used only in robustness specs.
capture drop tot_chg exchange2 ex_dum1 ex_dum2 ex_dum3 ex_dum4 ex_dum5
* terms of trade (optional robustness source) — skip if not present
capture confirm variable tot
if !_rc {
    gen double tot_chg = 100*(ln(L.tot) - ln(L2.tot)) if L.tot>0 & L2.tot>0 & !missing(L.tot,L2.tot)
    label var tot_chg "L1 terms-of-trade log-change, % (WDI; robustness)"
}
else di as error "  ** tot not found — skipping tot_chg (robustness only; add termsoftrade.xlsx to build it)."
* official FX (optional robustness source) — build exchange2 + ex_dum bins, else skip
capture confirm variable exch
if !_rc {
    gen double exchange2 = ln(1+L.exch) - ln(1+L2.exch) if !missing(L.exch,L2.exch)
    label var exchange2 "L1 nominal exchange-rate log-change (robustness)"
    quietly summarize exchange2 if continuation==0 & !missing(ln_gdp_base), detail
    foreach k in 1 2 3 4 5 {
        gen byte ex_dum`k' = 0 if !missing(exchange2)
    }
    replace ex_dum1 = 1 if exchange2 < r(p5)
    replace ex_dum2 = 1 if exchange2 < r(p25)
    replace ex_dum2 = 0 if exchange2 < r(p5)
    replace ex_dum3 = 1 if exchange2 < r(p50)
    replace ex_dum3 = 0 if exchange2 < r(p25)
    replace ex_dum4 = 1 if exchange2 < r(p75)
    replace ex_dum4 = 0 if exchange2 < r(p50)
    replace ex_dum5 = 1 if exchange2 < r(p95)
    replace ex_dum5 = 0 if exchange2 < r(p75)
    label var ex_dum1 "FX-change bin < p5 (Asonuma ex_dum; robustness)"
}
else di as error "  ** exch not found — skipping exchange2/ex_dum (robustness only; add officialexchangerate.xlsx)."

* ── CHANNEL OUTCOMES ON THE PAPER'S SCALE: LOG REAL LEVELS ─────────────────
* The reference paper does NOT use ratios to GDP as channel outcomes. It uses log
* real LEVELS, built by multiplying the ratio back by real GDP:
*     var2 = ln(investment*gdp_real)*100      (real investment)
*     var3 = ln(gdp_real*credit_bank)*100     (real bank credit)
* exactly as its GDP outcome is ln(gdp_real)*100.
*
* This is not cosmetic. A change in X/GDP confounds what happened to X with what
* happened to GDP: in a crisis the denominator falls, so the ratio rises even when
* X is flat or falling, and the "channel" then reports the recession rather than
* the mechanism. Taking logs and differencing removes the denominator entirely —
* ln(X/GDP * GDP) = ln(X) up to a constant — so the outcome is the cumulative
* percent change in X itself, on the same scale as dy_h. That also makes the two
* directly comparable: "credit fell x% while GDP fell y%".
*
* Only strictly-positive GDP-ratio channels can take a log, so pb (primary
* balance), fdi (net inflows) and ca (current account) keep their ratio form —
* all three change sign, and for a balance the ratio is the object of interest
* anyway. The nexus variables (claimsgov_assets, claimpriv_assets) are shares of
* BANK ASSETS, not of GDP, so the denominator problem does not arise and they are
* left alone too.
*
* The *_base and pre_* terms in the channel files are built from these columns, so
* the levels flow through to the outcomes, the own pre-trend control and the AIPW.
capture drop ln_r_credit ln_r_credit_bank ln_r_inv ln_r_govexp ln_r_claims_govt
foreach v in credit credit_bank inv govexp claims_govt {
    capture confirm variable `v'
    if _rc {
        di as error "  ** 18_transforms: `v' not found — ln_r_`v' not built."
        continue
    }
    * ln(gdp_real * ratio)*100, the paper's form. The ratio is in percent, so this
    * is the real level up to a factor of 100 — a constant inside the log, which
    * differences away in F^h - L. Guarded so a zero or negative ratio yields
    * missing rather than a Stata log-of-nonpositive.
    gen double ln_r_`v' = ln(gdp_real * `v')*100 ///
        if `v' > 0 & gdp_real > 0 & !missing(`v', gdp_real)
    label var ln_r_`v' "Log real `v' level x100 (= paper's ln(gdp_real*x)*100)"
}

* ln_r_claims_govt is built above but NOT used as the claims_govt channel outcome.
* Bank claims on central government is close to a NET position (credit to govt minus
* govt deposits at banks), so the ratio sits near zero for many country-years and is
* negative for some — ln() requires strict positivity, so ln_r_claims_govt drops the
* low-claims country-years outright (342 missing in one run; onset coverage at h=0
* fell 55->47). That is selection on one tail of the outcome, not noise, so the
* channel file keeps claims_govt on its ratio-to-GDP form (ppt of GDP) instead.
* asinh gives a level-based alternative without that cost: it is defined at zero and
* for negative values, behaves like ln() in the tails, and reads approximately as a
* percent change, so it loses no observations to a positivity requirement.
capture drop as_r_claims_govt
capture confirm variable claims_govt
if !_rc {
    gen double as_r_claims_govt = asinh(gdp_real * claims_govt)
    label var as_r_claims_govt "asinh(gdp_real*claims_govt): level robustness for claims_govt (see note above)"
}

* ══════════════════════════════════════════════════════════════════════════
* FLOW TREATMENT — "being in a spread crisis" (consumed by 20_lp_flow.do)
*
* The headline design treats an episode as a point event: onset_all is 1 in the
* first year only and `sample' below drops every continuation year. That answers
* "what follows the START of an episode". It cannot answer "what is the output
* cost of BEING IN a spread crisis", because the 173 country-years during which
* countries were actually in their crises are discarded.
*
* This block builds the parallel objects for that second question. It ADDS
* columns; it changes nothing existing, so every estimate in 02-13 is untouched.
*
* TREATMENT = THE ANNUAL CRITERION, NOT EPISODE MEMBERSHIP.
* Earlier versions of this file (through commit e94309f's revert) treated a
* whole dated episode as the treated unit: onset_all | continuation, bridging
* 13 "gap" rows where the annual criterion (crisis_any = crit1|crit2, i.e. the
* EMBIG spread actually clearing 1000bp or the HRT21 QoQ rule) is not met that
* specific year. Checking the reference paper's own replication script and
* dataset (Figure 3's construction, dum1/dum2/dum3) showed their treatment
* flags are onset-year-only — an ONSET-CODED design, not a design that bridges
* intervening years the way this project's own flow tier does. That discovery
* does not by itself settle what the FLOW tier (a genuinely different,
* broader question this project asks beyond the onset tier — "what is the
* cost of BEING IN a crisis", not just "what follows its start") should bridge
* through, but on reflection the annual criterion is the more defensible
* choice for it: a year is treated if and only if it independently clears the
* same threshold that defines every other treated year, full stop. Bridging a
* row that itself does NOT clear the criterion asks the estimator to price a
* year as "in crisis" using no crisis evidence of its own — the episode-dating
* rule's gap clause is a useful device for NAMING a sequence of years as one
* episode (Decision 2 below keeps it for exactly that), but it is a much
* weaker basis for TREATMENT than for classification.
*
* in_crisis THEREFORE = crisis_any==1 & carryin==0 (annual criterion,
* previously available only as the unused `in_crisis_sp'), replacing episode
* membership as the default. The old episode-membership definition is kept,
* unchanged in formula, under `in_crisis_episode' — an explicit, named
* robustness/historical comparison, not deleted. See 20_lp_flow.do Section 4
* for where it is reported.
*
* CONCRETE EFFECT (see the diagnostic block below the assertions for the
* live-derived counts): the 13 gap rows below the annual criterion drop out
* of the treated set entirely and rejoin the ordinary tranquil/control pool,
* rather than being priced as treated with no crisis evidence of their own.
* They are not homogeneous:
*   9 rule-consistent one-year gaps (a crisis year on BOTH sides) — Brazil 2000,
*     Cote d'Ivoire 2009, Ecuador 1997, Nigeria 1997, Nigeria 2021, Pakistan
*     2010, Venezuela 1997, Venezuela 2000, Zambia 2017.
*   2 spanning a TWO-year gap — Brazil 1996, 1997 (crisis in 1995 and 1998).
*   2 TRAILING — Ukraine 2010, 2011 (no crisis follows; 2014 is a fresh onset).
* Documented in DATA_SOURCES.md section 1.
* ══════════════════════════════════════════════════════════════════════════
capture drop in_crisis ep_seq nd_ep in_crisis_nd in_crisis_def in_crisis_sp in_crisis_episode sample_flow gap_year gap_year_dropped

gen byte in_crisis = (crisis_any==1) & carryin==0
label var in_crisis "In a spread crisis (annual criterion: crit1|crit2 met that year)"

* Old definition, preserved as a named robustness/historical variant.
gen byte in_crisis_episode = (onset_all==1 | continuation==1) & carryin==0
label var in_crisis_episode "In a spread crisis (OLD definition: episode membership, onset or continuation — robustness/historical)"

* The rows that were treated under the old definition but are not under the
* new one: inside a dated episode, below the annual criterion.
gen byte gap_year_dropped = (in_crisis_episode==1 & in_crisis==0)
label var gap_year_dropped "Excluded from in_crisis under the annual-criterion redefinition (was treated under in_crisis_episode)"

* ── Episode sequence number: a CLASSIFICATION grouping (which dated crisis
* sequence a row's resolution type is inherited from), independent of the
* redefinition above — it is still built from onset_all/continuation, exactly
* as before, and still needed to forward-fill the resolution type onto every
* now-treated row, since the onset row is by construction always itself a
* crisis_any==1 row (the annual criterion is what defines an onset in
* 10_skeleton.do's own dating rule) and so is always one of the treated rows
* under the new definition too. ──
* nondefault is merged in 10_skeleton on country x ONSET YEAR, so it is present
* on the 61 onset rows and MISSING on all 173 continuation rows. Splitting the
* flow treatment by resolution without filling it first would silently drop
* every continuation row — the exact rows this block exists to add — and would
* still run without error, reporting plausible numbers. Hence the fill and the
* assertions below.
bysort cid (year): gen int ep_seq = sum(onset_all)
label var ep_seq "Running episode counter within country (0 = before first onset)"

bysort cid ep_seq: egen byte nd_ep = max(nondefault)
label var nd_ep "Resolution type of the episode, filled to all its years (1=non-default)"

* ── EPISODE-DATED CONTROLS ──────────────────────────────────────────────────
* Under flow coding every element of the flow tier's core control set,
* measured at a continuation row's OWN t-1, is an outcome of the crisis that
* row is already in: lagged growth most obviously, but debt, bank credit,
* government spending and banking-crisis status equally. A row three years
* into an episode would be conditioned on covariates the episode itself
* produced — mediators, on the right-hand side.
*
* The reference paper never faces this because every one of its treated rows is
* an onset, so its L./L2. controls are predetermined by construction. Dating
* the controls at the EPISODE's entry year reproduces exactly that timing:
* the core control set as it stood the year before the episode began, held
* fixed for all of its years; tranquil rows keep their own t-1, which is
* predetermined for them trivially. Onset rows are unchanged, since for them
* t-1 IS the entry year — so the onset and flow specifications share an
* identical control set at Year 1 and their Year-1 coefficients are directly
* comparable (see the row-dated $ctrl_core robustness column in the flow
* files for the plain, non-entry-dated version of the same set).
*
* Same construction as Callaway & Sant'Anna (2021), who measure covariates in
* each cohort's pre-treatment period for this reason.
*
* $ctrl_flow is the flow BASELINE (entry-dated); $ctrl_core (row-dated) is
* retained as the robustness column in 20_lp_flow.do, so the table shows
* what the entry-dating choice costs.
*
* ══════════════════════════════════════════════════════════════════════════
* CORE CONTROL SET — $ctrl_core / $ctrl_flow
*
* Used PROJECT-WIDE: every onset-tier file (02, 03, 06, 07, 08, 08b, 11,
* 11b, 12, 12b, 13, 13b, 13c, 13d, 13f, 19) reads $ctrl_core directly
* (row-dated, plain L. lags); every flow-tier file (20-26) reads $ctrl_flow
* (this set, entry-dated) as its outcome-model controls and $ctrl_core
* (row-dated) as its propensity/robustness set. One control-set definition,
* not two parallel ones.
*
* Three swaps relative to the set this project used before this adoption
* (l1_gdpg, l_debt, l_ca, l_banking_duration, l_govexp, l_open,
* l_credit_bank, l_hyperinfl):
*   l_banking_duration (years-so-far)   -> l_banking_crisis (0/1 dummy)
*   l_hyperinfl (L.infl>50 dummy)       -> l_lninfl (continuous log inflation)
*   l_ca (current account)              -> exchange2 (log FX change)
* The second swap was made after 21_aipw_flow.do's diagnostic history (its
* def-arm probit coefficient tables) repeatedly showed l_hyperinfl driving
* quasi-complete separation on its own -- a binary flag that is 1 almost
* exclusively on default-linked country-years functions close to a
* near-perfect predictor of the def arm by construction, not a genuine
* covariate. l_lninfl (= ln(1+L.infl/100), already built above as a
* robustness alternative before this adoption) keeps the same underlying
* inflation information as a continuous variable, which cannot by itself
* perfectly separate treated from control the way a threshold dummy can.
*
* THE THIRD SWAP, l_ca -> exchange2, IS ON LITERATURE GROUNDS, NOT PURELY
* EMPIRICAL PERFORMANCE -- state both sides plainly. The reference paper's
* own $convar carries exchange-rate depreciation (their ex_dum1-ex_dum5
* percentile bins) as a baseline control; the current account is not in
* their $convar at all. exchange2 (the continuous log FX change) is used
* here rather than their literal ex_dum1-ex_dum5 bins:
* 21b_first_stage_table_flow.do's flow_ctrl_variant testing found the bins
* separate on this project's much smaller panel exactly the way country FE
* and past_def_onsets did (ex_dum1, and in the def arm ex_dum2 as well, have
* zero outcome variation and Stata drops the dummy and every row in the bin
* automatically), so the continuous level is the workable proxy, not a
* literal reproduction of their construction. tot_chg (terms-of-trade
* log-change, this project's own earlier addition, never in the reference
* paper's $convar) is available as a robustness-tier alternative to
* exchange2 (see $ctrl_core_flowplus below) rather than folded into the
* core, since exchange2 is the term with the literature-fidelity claim.
*
* CONSEQUENCE FOR EVERY PUBLISHED NUMBER, STATED PLAINLY: this control set
* now governs Table 1, Table 2, Table 3, the AIPW tables (08b/13c/13d), and
* every channel table (11/12), not just the flow tier -- every one of those
* numbers changes relative to what this project reported before this
* adoption, because the regression's own control set changed. That is the
* intended effect of unifying the control set, not a side effect to correct.
* Whether l_hyperinfl's known separation risk in the FLOW AIPW's propensity
* model also affects the ONSET-tier AIPW files (08b/13c/13d) is a separate
* empirical question this adoption does not by itself answer -- it removes
* a KNOWN risk in one place, and the onset-tier propensity models should be
* checked directly (not assumed clean) after this change.
* ══════════════════════════════════════════════════════════════════════════
capture confirm variable exchange2
if _rc {
    di as error "  ** exchange2 not built (exch missing) -- \$ctrl_core needs it as a"
    di as error "     core term now, not an optional one. Add data/raw/officialexchangerate.xlsx"
    di as error "     and re-run 01_build_panel.do/12_wdi.do before this file."
    exit 111
}
global ctrl_core "l1_gdpg l_debt l_banking_crisis l_govexp l_open l_credit_bank l_lninfl exchange2"

global ctrl_flow ""
foreach X of global ctrl_core {
    capture drop epc_`X' _ent_`X'
    quietly bysort cid ep_seq: egen double _ent_`X' = max(cond(onset_all==1, `X', .))
    quietly gen double epc_`X' = cond(in_crisis==1, _ent_`X', `X')
    label var epc_`X' "`X' at episode entry (tranquil rows keep own t-1)"
    quietly drop _ent_`X'
    global ctrl_flow "$ctrl_flow epc_`X'"
}
di as result "  FLOW: episode-dated control set built -> $ctrl_flow"

* Onset rows must be untouched by the re-dating: for them t-1 IS the entry year.
local nbad = 0
foreach X of global ctrl_core {
    quietly count if onset_all==1 & carryin==0 & !missing(`X') & abs(epc_`X' - `X') > 1e-9
    local nbad = `nbad' + r(N)
}
if `nbad' != 0 di as error "  ** FLOW: episode-dated controls differ from row-dated on `nbad' ONSET rows — ep_seq is wrong"
else           di as result "  FLOW: episode-dated controls match row-dated on every onset row (correct)"

* ══════════════════════════════════════════════════════════════════════════
* ALTERNATE (BIGGER) CONTROL SET — $ctrl_core_flowplus / $ctrl_flow_flowplus
*
* Core PLUS two candidates not folded into the core itself:
*   - tot_chg: terms-of-trade log-change, the term exchange2 was preferred
*     over on literature-fidelity grounds (see "CORE CONTROL SET" above) --
*     kept available here so that choice can be tested rather than forced.
*   - l_imf: an IMF-supported-program dummy, lagged (predetermined). `imf' is
*     already in panel_build.dta (01_build_panel.do); this is the first place
*     it is lagged and entry-dated for use as a CONTROL rather than as the
*     treatment modifier 09_lp_imf.do studies it as.
* Freedom House civil-liberties/political-rights indices were also
* considered but are NOT built here -- no raw source file for them exists in
* this project yet; add them the same way GEO_CEPII.xlsx was added for the
* contagion predictor if that data becomes available.
*
* NOT used anywhere by default. Selected via each flow file's
* `flow_ctrl_variant' toggle (0=core, 1=this set), default 0. Same
* episode-dating mechanism as $ctrl_flow above, just parameterised on a
* different base list -- not a new pattern.
* ══════════════════════════════════════════════════════════════════════════
* The preceding bysort loop re-sorts the physical dataset by cid ep_seq as a
* side effect, same as every other bysort/egen loop in this file -- L.imf
* then needs the data back in cid year order to match what xtset declared,
* or it errors r(5) "not sorted" (confirmed the hard way; belt-and-suspenders
* here as elsewhere).
sort cid year
capture drop l_imf
gen byte l_imf = L.imf
label var l_imf "L1 IMF-supported program dummy (predetermined)"

capture confirm variable tot_chg
if _rc {
    di as error "  ** tot_chg not built (tot missing) -- \$ctrl_core_flowplus will exclude it."
    di as error "     add data/raw/termsoftrade.xlsx and re-run 01_build_panel.do first if wanted."
    global ctrl_core_flowplus "$ctrl_core l_imf"
}
else global ctrl_core_flowplus "$ctrl_core tot_chg l_imf"
global ctrl_flow_flowplus ""
foreach X of global ctrl_core_flowplus {
    * Dropped SEPARATELY: the 8 $ctrl_core terms already have epc_ built
    * from the $ctrl_flow loop above, so only _ent_`X' is guaranteed absent
    * for those (tot_chg/l_imf are new here and need both dropped).
    capture drop epc_`X'
    capture drop _ent_`X'
    quietly bysort cid ep_seq: egen double _ent_`X' = max(cond(onset_all==1, `X', .))
    quietly gen double epc_`X' = cond(in_crisis==1, _ent_`X', `X')
    label var epc_`X' "`X' at episode entry (tranquil rows keep own t-1) -- flowplus"
    quietly drop _ent_`X'
    global ctrl_flow_flowplus "$ctrl_flow_flowplus epc_`X'"
}
di as result "  FLOWPLUS (exploratory): episode-dated control set built -> $ctrl_flow_flowplus"

gen byte in_crisis_nd  = (in_crisis==1 & nd_ep==1)
gen byte in_crisis_def = (in_crisis==1 & nd_ep==0)
label var in_crisis_nd  "In a NON-DEFAULT spread crisis"
label var in_crisis_def "In a DEFAULT-LINKED spread crisis"

* ── Estimation sample ───────────────────────────────────────────────────────
* carryin==0 excludes the pre-EMBIG scaffolding rows added in 10_skeleton: they exist
* only so the L. operators above have a previous row to point at, and carry no spread
* data, so they must never enter an estimation.
* atonly_country==0 excludes 10b_skeleton_atonly.do's AT-only countries (outside
* this project's 52-country spread panel): they were never tested against the
* spread criterion, so including them here would silently add untested-for-
* spread control observations to EVERY existing onset/flow file's (02-25) sample
* -- exactly what METHODOLOGY.md §5.3 says does not happen. Only
* 26_lp_debtcrisis_flow.do builds its own broadened sample flag that includes
* them. capture confirm keeps this file runnable standalone if 10b never ran.
capture confirm variable atonly_country
if _rc gen byte atonly_country = 0
capture drop sample
gen byte sample = (continuation==0) & !missing(ln_gdp_base) & carryin==0 & atonly_country==0
label var sample "Estimation sample (onset + tranquil, excl. continuation, carry-in & AT-only countries, GDP base present)"

* Flow estimation sample: identical to `sample' EXCEPT that continuation years
* are kept. Required because every existing regression is `if sample==1', which
* by construction contains zero treated rows beyond each episode's onset year.
* sample is nested inside sample_flow (asserted below).
gen byte sample_flow = !missing(ln_gdp_base) & carryin==0 & atonly_country==0
label var sample_flow "Flow estimation sample (all episode years + tranquil, excl. carry-in)"

* ── Assertions on the flow objects ──────────────────────────────────────────
* These are cheap and they catch the two failure modes that would otherwise be
* invisible: a treatment that silently lost the continuation rows, and a
* resolution type that failed to propagate off the onset row.
* Under the annual-criterion redefinition, in_crisis's row count is no longer
* a fixed, derivable-in-advance number the way episode membership was (it
* depends on how many mid-episode years happen to clear crisis_any each
* year) — report it live rather than asserting a hard-coded expectation.
quietly count if in_crisis_episode==1
local n_old_ep = r(N)
if `n_old_ep' != 234 di as error "  ** FLOW: in_crisis_episode = `r(N)' rows, expected 234 — check onset_all/continuation/carryin"
quietly count if in_crisis==1
local n_new_ann = r(N)
di as result "  FLOW REDEFINITION: treated rows `n_old_ep' (old, episode membership) -> `n_new_ann' (new, annual criterion)"
quietly count if in_crisis_episode==1 & in_crisis==0
local n_dropped_cont = r(N)
quietly count if onset_all==1 & carryin==0 & in_crisis==0
local n_dropped_onset = r(N)
if `n_dropped_onset' != 0 di as error "  ** FLOW: `n_dropped_onset' onset rows do not clear crisis_any — onset dating is inconsistent with the annual criterion"
di as result "    dropped: `n_dropped_cont' continuation rows below the annual criterion (" %4.1f 100*`n_dropped_cont'/`n_old_ep' " pct of the old total)"
quietly count if gap_year_dropped==1 & !(in_crisis_episode==1 & in_crisis==0)
if r(N) != 0 di as error "  ** FLOW: gap_year_dropped disagrees with the direct in_crisis_episode/in_crisis count on `r(N)' rows"
else         di as result "  FLOW: gap_year_dropped matches the direct count (correct)"

quietly count if in_crisis==1 & missing(nd_ep)
if r(N) != 0 di as error "  ** FLOW: `r(N)' episode-years have no resolution type — the nd_ep fill FAILED"

* The fill must reproduce the source value on the onset row it came from. This
* is the check that actually bites: egen max() over the wrong group would still
* give a constant nd_ep, so testing constancy alone proves nothing.
quietly count if onset_all==1 & carryin==0 & nd_ep != nondefault
if r(N) != 0 di as error "  ** FLOW: nd_ep disagrees with nondefault on `r(N)' onset rows — ep_seq groups are wrong"

* Each episode must carry exactly one onset row; more means ep_seq merged two.
quietly preserve
    quietly keep if in_crisis==1
    quietly collapse (sum) nons=onset_all, by(cid ep_seq)
    quietly count if nons != 1
    local badep = r(N)
    quietly count
    local nep_flow = r(N)
quietly restore
if `badep' != 0 di as error "  ** FLOW: `badep' episodes do not have exactly one onset row"
if `nep_flow' != 61 di as error "  ** FLOW: `nep_flow' distinct episodes, expected 61"
else                di as result "  FLOW: 61 distinct episodes recovered from ep_seq (correct)"

quietly count if sample==1 & sample_flow==0
if r(N) != 0 di as error "  ** FLOW: `r(N)' rows in sample but not sample_flow — nesting violated"

quietly count if in_crisis_nd==1
local n_ndf = r(N)
quietly count if in_crisis_def==1
local n_deff = r(N)
di as result "  FLOW: crisis-years by resolution — non-default `n_ndf', default-linked `n_deff'"
di as result "        (expect 113 / 121 with the Venezuela-2008 override at 10_skeleton.do:119)"

* Gap-year breakdown — the rows inside a dated episode but below the annual
* criterion, now excluded from in_crisis (see the redefinition above).
quietly count if gap_year_dropped==1
local ngap = r(N)
di as result "  FLOW: gap years dropped (in old episode, criterion not met): `ngap' of `n_old_ep' old-definition treated rows"
if `ngap' > 0 {
    di as result "        listed below; 9 are one-year gaps, 2 span a two-year gap (Brazil)"
    di as result "        and 2 are trailing (Ukraine) — see the block above."
    quietly levelsof country if gap_year_dropped==1, local(gapc)
    foreach c of local gapc {
        quietly levelsof year if gap_year_dropped==1 & country=="`c'", local(gy)
        di as result "          `c': `gy'"
    }
}

* Treated share of each country's panel years — flow coding concentrates the
* treatment in a handful of chronic cases (Venezuela is in crisis in most of its
* panel), which matters for how much the country fixed effects can absorb.
di as result _n "  FLOW: countries with the highest treated share of panel years"
preserve
    quietly keep if carryin==0
    collapse (sum) ncris=in_crisis (count) nyr=year, by(country)
    gen double share = 100*ncris/nyr
    gsort -ncris
    quietly count
    local nshow = min(8, r(N))
    forvalues i = 1/`nshow' {
        di as result "        " %-16s country[`i'] %3.0f ncris[`i'] " of " %3.0f nyr[`i'] " years (" %4.1f share[`i'] "%)"
    }
restore

* ── Save the analysis file (drop-in replacement consumed by 02..16) ─────────
compress
sort cid year
xtset cid year
save "$clean/panel_lp.dta", replace

* ══════════════════════════════════════════════════════════════════════════
* VALIDATION SUMMARY
* ══════════════════════════════════════════════════════════════════════════
di as result _n "18_transforms.do complete — panel_lp.dta rebuilt from source."
quietly count if sample==1
di as result "  sample rows: `r(N)'"
quietly count if onset_all==1 & sample==1
di as result "  onsets in sample: `r(N)'"

* Carry-in scaffolding must never leak into an estimation sample.
quietly count if carryin==1
local ncar = r(N)
quietly count if sample==1 & carryin==1
if r(N) == 0 di as result "  carry-in rows: `ncar' present, 0 in sample (correct)"
else         di as error  "  ** `r(N)' CARRY-IN ROWS LEAKED INTO sample==1 — investigate"
di as result _n "Coverage at onsets (non-missing) for the key analysis variables:"
foreach v in dy_0 l1_gdpg debt ca infl l_hyperinfl l_lninfl imf credit fdi claims_govt inv govexp pb ///
             banking_crisis l_banking_crisis reer_chg l_reer_chg tot_chg exchange2 revenue_gdp open ///
             claimsgov_assets claimpriv_assets vix fedfunds ust10y {
    capture confirm variable `v'
    if !_rc {
        quietly count if onset_all==1 & !missing(`v')
        di as result "    `v': `r(N)' / 61"
    }
    else di as error "    `v': MISSING VARIABLE"
}

* ── Inflation control: confirm the log transform tamed the tail ─────────────
* l_hyperinfl (the retired dummy) is kept built as a robustness alternative,
* not part of $ctrl_core any more -- l_lninfl (continuous) is the core term
* project-wide now, see "CORE CONTROL SET" above. These lines report how many
* observations trip the 50% threshold and what the continuous alternative
* looks like, so the choice stays visible in the run log.
capture confirm variable l_lninfl
if !_rc {
    di as result _n "Inflation control (l_hyperinfl = L.infl > 50, robustness alt. only;"
    di as result "                    l_lninfl, continuous, is the CORE term project-wide):"
    quietly count if l_hyperinfl==1 & sample==1
    di as result "    l_hyperinfl==1 in estimation sample: `r(N)'"
    quietly summarize l_infl if sample==1, detail
    di as result "    raw  l_infl   median " %8.2f r(p50) "   max " %11.1f r(max)
    quietly summarize l_lninfl if sample==1, detail
    di as result "    log  l_lninfl median " %8.3f r(p50) "   max " %11.3f r(max) ///
                 "  min " %8.3f r(min) "  (core term, project-wide)"
}

* ── Core control set: distribution sanity check for every term ─────────────
* Verifies no term in $ctrl_core has a tail extreme enough to worry about
* (the reason l_hyperinfl -> l_lninfl and the raw ln(x/GDP) checks exist
* above) -- printed rather than assumed, same standard applied to inflation.
* min/max/skewness flagged if skewness exceeds 3 in absolute value, an
* informal threshold, not a hard rule.
di as result _n "Core control set (\$ctrl_core) — distribution check:"
foreach X of global ctrl_core {
    quietly summarize `X' if sample_flow==1, detail
    if r(N) > 0 {
        local skew = cond(r(sd) > 0, ///
            (r(mean) - r(p50)) / r(sd) * 3, .)   // Pearson's 2nd skewness coefficient, a cheap proxy
        local flag = cond(!missing(`skew') & abs(`skew') > 3, "  ** long tail, check", "")
        di as result "    " %-20s "`X'" "min " %9.2f r(min) "  median " %9.2f r(p50) ///
                     "  max " %9.2f r(max) "  skew~" %6.2f `skew' "`flag'"
    }
}
