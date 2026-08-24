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

global ctrl_core "l1_gdpg l_debt l_ca l_banking_duration l_govexp l_open l_credit_bank l_hyperinfl"

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
* TREATMENT = EPISODE MEMBERSHIP, NOT THE ANNUAL CRITERION FLAG.
* The workbook's episode rule (README sheet, after Detragiache & Spilimbergo
* 2001) is: onset = criterion met in t but not t-1; continuation = consecutive
* crisis years OR a one-year gap followed by a crisis year; end = two
* consecutive tranquil years. The gap clause exists to stop a single prolonged
* distress fragmenting into several "episodes".
*
* So `crisis_any' (the annual criterion flag) is NOT episode membership: it is 0
* on 13 rows the dating rule places INSIDE an episode. Using it would contradict
* the rule that generated the 61 episodes and would push mid-episode years into
* the tranquil CONTROL pool — worse than either including or dropping them.
* in_crisis therefore = onset_all | continuation (234 rows vs 221).
*
* The 13 gap rows are not homogeneous, and the breakdown is printed below rather
* than left implicit:
*   9 rule-consistent one-year gaps (a crisis year on BOTH sides) — Brazil 2000,
*     Cote d'Ivoire 2009, Ecuador 1997, Nigeria 1997, Nigeria 2021, Pakistan
*     2010, Venezuela 1997, Venezuela 2000, Zambia 2017.
*   2 spanning a TWO-year gap — Brazil 1996, 1997 (crisis in 1995 and 1998).
*   2 TRAILING — Ukraine 2010, 2011 (no crisis follows; 2014 is a fresh onset).
* The last four are exceptions to the "two consecutive tranquil years ends the
* episode" rule. They are retained deliberately (the dating is not revised), but
* note what strict application would imply: Brazil would split into two episodes,
* creating a 1998 onset that does not currently exist and taking the count to 62.
* Documented in DATA_SOURCES.md section 1.
* ══════════════════════════════════════════════════════════════════════════
capture drop in_crisis ep_seq nd_ep in_crisis_nd in_crisis_def in_crisis_sp sample_flow gap_year

gen byte in_crisis = (onset_all==1 | continuation==1) & carryin==0
label var in_crisis "In a spread crisis (episode membership: onset or continuation)"

* Annual-criterion variant, for the robustness column only.
gen byte in_crisis_sp = (crisis_any==1) & carryin==0
label var in_crisis_sp "In a spread crisis (annual criterion flag; robustness variant)"

* The gap rows: inside an episode but below the annual criterion.
gen byte gap_year = (in_crisis==1 & crisis_any==0)
label var gap_year "Mid-episode year below the annual crisis criterion"

* ── Episode sequence number, needed to forward-fill the resolution type ──────
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
* Under flow coding every element of $ctrl_core measured at a continuation
* row's OWN t-1 is an outcome of the crisis that row is already in: lagged
* growth most obviously, but debt, the current account, bank credit, government
* spending and banking-crisis duration equally. A row three years into an
* episode would be conditioned on covariates the episode itself produced —
* mediators, on the right-hand side.
*
* The reference paper never faces this because every one of its treated rows is
* an onset, so its L./L2. controls are predetermined by construction. Dating
* the controls at the EPISODE's entry year reproduces exactly that timing:
* $ctrl_core as it stood the year before the episode began, held fixed for all
* of its years; tranquil rows keep their own t-1, which is predetermined for
* them trivially. Onset rows are unchanged, since for them t-1 IS the entry
* year — so the onset and flow specifications share an identical control set at
* Year 1 and their Year-1 coefficients are directly comparable.
*
* Same construction as Callaway & Sant'Anna (2021), who measure covariates in
* each cohort's pre-treatment period for this reason.
*
* $ctrl_flow is the flow BASELINE; $ctrl_core is retained as the row-dated
* robustness column in 20_lp_flow.do, so the table shows what the choice costs.
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
* EXPLORATORY ALTERNATE FLOW CONTROL SET — $ctrl_core_flowalt / $ctrl_flow_flowalt
*
* NOT used anywhere by default. Drops l_debt/l_ca, swaps the banking-crisis
* DURATION (l_banking_duration) for the DUMMY (l_banking_crisis), and adds
* the exchange rate (l_reer_chg) and terms of trade (tot_chg) -- a one-off
* "what happens if" test of the flow tier's control set, requested to be kept
* separate from $ctrl_core/$ctrl_flow so every onset-tier file and every
* flow-file identity check that depends on matching $ctrl_core keeps working
* unchanged. Each flow file (20-24) reads this only if its own
* `use_flowalt_ctrl' toggle is set to 1; default is 0 (current baseline).
* Same episode-dating mechanism as $ctrl_flow above, just parameterised on a
* different base list -- not a new pattern.
global ctrl_core_flowalt "l1_gdpg l_govexp l_open l_credit_bank l_hyperinfl l_banking_crisis l_reer_chg tot_chg"
global ctrl_flow_flowalt ""
foreach X of global ctrl_core_flowalt {
    capture drop epc_`X' _ent_`X'
    quietly bysort cid ep_seq: egen double _ent_`X' = max(cond(onset_all==1, `X', .))
    quietly gen double epc_`X' = cond(in_crisis==1, _ent_`X', `X')
    label var epc_`X' "`X' at episode entry (tranquil rows keep own t-1) -- flowalt"
    quietly drop _ent_`X'
    global ctrl_flow_flowalt "$ctrl_flow_flowalt epc_`X'"
}
di as result "  FLOWALT (exploratory): episode-dated control set built -> $ctrl_flow_flowalt"

gen byte in_crisis_nd  = (in_crisis==1 & nd_ep==1)
gen byte in_crisis_def = (in_crisis==1 & nd_ep==0)
label var in_crisis_nd  "In a NON-DEFAULT spread crisis"
label var in_crisis_def "In a DEFAULT-LINKED spread crisis"

* ── Estimation sample ───────────────────────────────────────────────────────
* carryin==0 excludes the pre-EMBIG scaffolding rows added in 10_skeleton: they exist
* only so the L. operators above have a previous row to point at, and carry no spread
* data, so they must never enter an estimation.
capture drop sample
gen byte sample = (continuation==0) & !missing(ln_gdp_base) & carryin==0
label var sample "Estimation sample (onset + tranquil, excl. continuation & carry-in, GDP base present)"

* Flow estimation sample: identical to `sample' EXCEPT that continuation years
* are kept. Required because every existing regression is `if sample==1', which
* by construction contains zero treated rows beyond each episode's onset year.
* sample is nested inside sample_flow (asserted below).
gen byte sample_flow = !missing(ln_gdp_base) & carryin==0
label var sample_flow "Flow estimation sample (all episode years + tranquil, excl. carry-in)"

* ── Assertions on the flow objects ──────────────────────────────────────────
* These are cheap and they catch the two failure modes that would otherwise be
* invisible: a treatment that silently lost the continuation rows, and a
* resolution type that failed to propagate off the onset row.
quietly count if in_crisis==1
if r(N) != 234 di as error "  ** FLOW: in_crisis = `r(N)' rows, expected 234 — check onset_all/continuation/carryin"
else           di as result "  FLOW: in_crisis = 234 rows (61 onset + 173 continuation)"

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

* Gap-year breakdown — the 13 rows inside an episode but below the annual criterion.
quietly count if gap_year==1
local ngap = r(N)
di as result "  FLOW: gap years (in episode, criterion not met): `ngap' of 234 treated rows"
if `ngap' > 0 {
    di as result "        listed below; 9 are one-year gaps, 2 span a two-year gap (Brazil)"
    di as result "        and 2 are trailing (Ukraine) — see the block above."
    quietly levelsof country if gap_year==1, local(gapc)
    foreach c of local gapc {
        quietly levelsof year if gap_year==1 & country=="`c'", local(gy)
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
             banking_crisis l_banking_crisis reer_chg l_reer_chg tot_chg revenue_gdp open claimsgov_assets ///
             claimpriv_assets vix fedfunds ust10y {
    capture confirm variable `v'
    if !_rc {
        quietly count if onset_all==1 & !missing(`v')
        di as result "    `v': `r(N)' / 61"
    }
    else di as error "    `v': MISSING VARIABLE"
}

* ── Inflation control: confirm the log transform tamed the tail ─────────────
* The core carries l_hyperinfl. These lines report how many observations trip the
* 50% threshold and what the continuous alternatives look like, so the choice stays
* visible in the run log.
capture confirm variable l_lninfl
if !_rc {
    di as result _n "Inflation control (core uses l_hyperinfl = L.infl > 50):"
    quietly count if l_hyperinfl==1 & sample==1
    di as result "    l_hyperinfl==1 in estimation sample: `r(N)'"
    di as result "       (was 1 pre-rebuild, when the panel was truncated at 2018)"
    quietly summarize l_infl if sample==1, detail
    di as result "    raw  l_infl   median " %8.2f r(p50) "   max " %11.1f r(max)
    quietly summarize l_lninfl if sample==1, detail
    di as result "    log  l_lninfl median " %8.3f r(p50) "   max " %11.3f r(max)   "  (alt., not in core)"
}
