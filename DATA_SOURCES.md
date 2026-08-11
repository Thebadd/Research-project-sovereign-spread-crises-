# Data Sources

Every macro variable in this project is rebuilt from official sources in the
`10`–`18` build chain (`do/00_master.do`), keyed on **ISO3 × year**. The only
inputs kept from prior work are (a) the author's own **spread-crisis database**
(onsets, spreads, default/non-default classification) and (b) one **legacy**
column (`imf`) whose upstream source could not be recovered. Nothing is inherited
from the old opaque master CSV.

Each row below: variable → build stage → raw file → series code/sheet → transform.

## 1. Crisis definition — author's own spread-crisis DB
Source file: `data/raw/EM_Spread_Crisis_DB_FINAL.xlsx` (spreads from **JP Morgan
EMBIG Global**). Built in `10_skeleton.do`.

| Variable | Origin | Notes |
|---|---|---|
| `onset_all` | Panel_Annual, "Onset D_it" | first year of each spread-crisis episode |
| `onset_nd` / `onset_def` | derived from `classification` | non-default / default-linked onset |
| `nondefault` | Episode_Summary "Classification" | 1 = non-default episode |
| `spr_max`, `spr_mean` | Panel_Annual | EMBIG max / mean annual spread (bps) |
| `crit1`, `crit2` | Panel_Annual | 1000bps / 99pct-HRT criterion flags |
| `region`, `continuation`, `crisis_any` | Panel_Annual | — |
| `iso3` | crosswalk in `10_skeleton.do` | 52-country ISO3 key |
| `carryin` | added in `10_skeleton.do` | 1 = pre-entry scaffolding row (see below); **never in an estimation sample** |

### Carry-in rows — lag scaffolding at the panel's left edge
The panel is **unbalanced on the left**: it begins when a country enters the JP Morgan
EMBIG index, not when its macro data begins (Armenia ≈2013, Kenya and Zambia ≈2014). Every
macro merge in stages 11–16 is `keep(master match)` with this skeleton as master, so the
sources' earlier coverage — WEO runs 1980–2031, WDI from 1960 — was discarded for all
pre-entry years. The lag operators in `18_transforms.do` then fell off the start of each
country's panel: a 2013 entrant with a 2013 onset had no `l1_gdpg`, which needs GDP at
*t−1* **and** *t−2*. That is why `l1_gdpg` covered only 42 of 61 onsets while `gdp_real`
covered 59, and why `l_banking_crisis` covered 52 of 61 even though `banking_crisis` is
non-missing on every row.

`10_skeleton.do` therefore appends **3 rows immediately before each country's first panel
year** (3 covers the deepest lag in use: `L3.ln_gdp` in `dy_m2`, `L2.gdpg` in `l2_gdpg`).
Because stages 11–16 merge on `iso3 × year` or `country × year`, they populate these rows
automatically, and all transform logic stays in stage 18.

These rows are **scaffolding only**. They are flagged `carryin == 1`, forced to
`onset_all = onset_nd = onset_def = crisis_any = continuation = 0`, excluded from
`sample_base` (stage 10) and from `sample` (stage 18), and excluded from the regional
contagion aggregates in `17_predictors.do` so they cannot dilute the member denominator.
They carry **no EMBIG quote** — `spr_max` / `spr_mean` stay missing by design — so the
construction asserts nothing about whether those pre-entry years were tranquil. Stage 18
prints a leak check confirming zero carry-in rows reach `sample == 1`.

This is a data-plumbing fix: it changes no specification and adds no source, it only stops
observations being deleted for want of data already on disk.

## 2. IMF World Economic Outlook — `11_weo.do`
Source: `data/raw/WEOApr2026all.xlsx`, sheet "Countries", by `INDICATOR.ID`.

| Variable | WEO code | Definition |
|---|---|---|
| `gdp_real` / `ln_gdp` | NGDP_R | real GDP (const. prices); `ln_gdp = ln(NGDP_R)` — **headline LP outcome base** |
| `gddpc_real` / `ln_gdppc` | NGDPRPC | real GDP per capita (const. prices) — per-capita robustness outcome (`dy_pc_*`) |
| `gdp_nominal` | NGDP | nominal GDP (for ratios) |
| `pop` | LP | population (millions) |
| `infl` | PCPIPCH | CPI inflation, period avg, % |
| `ca` | BCA_NGDPD | current account, % GDP |
| `debt` | GGXWDG_NGDP | general govt gross debt, % GDP |
| `govexp` | GGX_NGDP | general govt expenditure, % GDP |
| `revenue_gdp` | GGR_NGDP | general govt revenue, % GDP |
| `inv` | NID_NGDP | total investment, % GDP |
| `pb` | GGXONLB_NGDP | primary balance, % GDP |
| `gdp_last_actual` | `LATEST_ACTUAL_ANNUAL_DATA` (col. N, NGDPRPC rows) | last **outturn** year for real GDP p.c.; later years are IMF projections |

### Outturns vs projections
The Apr-2026 WEO vintage carries values through 2031, but only years up to
`LATEST_ACTUAL_ANNUAL_DATA` are actual data — mostly **2024 or 2025**, earlier for a few
countries. Because the LP outcome `dy_h` is GDP at *t+h*, an onset in 2022 or later has its
h=3/h=4 outcome built partly from forecasts, and those are exactly the horizons where the
estimated cost and the default premium are largest.

`11_weo.do` therefore carries the boundary into the panel as `gdp_last_actual` (the 4-digit
year is parsed out with a regex, since some countries report fiscal years as `FY2024/25`),
and reports each run how many onsets have a projection-based outcome at each horizon.
`02_lp_all.do` and `03_lp_resolution.do` use it for an **outturns-only robustness cut**
(`(year + h) <= gdp_last_actual`), exported as Tables 1R and 2R. Those tables are a
falsification check on the long horizons, not a replacement for the headline: the sample
necessarily shrinks at h=3/h=4, so wider intervals there are expected and carry no
information on their own — only a change in sign or magnitude would.

## 3. World Bank WDI — `12_wdi.do`
Classic WDI wide files (Country Code = ISO3) + REER DataBank file.

| Variable | WDI code | Definition |
|---|---|---|
| `credit` | **FS.AST.PRVT.GD.ZS** | domestic credit to private sector, all financial corp., % GDP (channel variable) |
| `credit_bank` | FD.AST.PRVT.GD.ZS | domestic credit to private sector **by banks**, % GDP (robustness) |
| `fdi` | BX.KLT.DINV.WD.GD.ZS | FDI net inflows, % GDP |
| `claims_govt` | FS.AST.CGOV.GD.ZS | claims on central government, % GDP |
| `exp_gdp` / `imp_gdp` | NE.EXP.GNFS.ZS / NE.IMP.GNFS.ZS | exports / imports, % GDP |
| `open` | = `exp_gdp + imp_gdp` | trade openness, % GDP (Asonuma control) |
| `reer_chg` | REER_INDEX.xlsx (2010=100) | REER, % YoY change |

## 4. IMF Monetary & Financial Statistics — `13_ifs_nexus.do`
Source: `data/raw/Nexus_Sovereign_Bank.xlsx` — aggregate **Other Depository
Corporations** (deposit-taking banks ex central bank) balance sheet, LCU levels,
2001–2024. All three are shares of **total bank assets** (%).

| Variable | Definition |
|---|---|
| `claimsgov_assets` | total bank claims on **general government** / total assets (headline doom-loop) |
| `netclaimsgov_assets` | claims on govt **minus govt deposits** / total assets |
| `claimpriv_assets` | bank claims on the **private sector** / total assets |

Coverage: 2001+ only; absent for **China (mainland), Ecuador, El Salvador, India,
Lebanon, Vietnam** (not in the IMF ODC file).

## 5. World Bank International Debt Statistics — `14_ids.do`
Files: `International_Debt_Statistics.xlsx` (+ `Total_debt_service.xlsx`).

| Variable | IDS code | Definition |
|---|---|---|
| `stdebt_share` | DT.DOD.DSTC.ZS | short-term debt, % total external debt (rollover) |
| `reserves_extdebt` | FI.RES.TOTL.DT.ZS | reserves, % total external debt (buffer) |
| `intpay_gni` | DT.INT.DECT.GN.ZS | interest on external debt, % GNI |
| `debt_service` | DT.TDS.DECT.EX.ZS | total debt service, % of exports (liquidity) |

## 6. FRED — `15_rates.do`
Global-push predictors (merged on year; absorbed by year FE in the LP, used only
in the first-stage propensity model).

| Variable | FRED series | File |
|---|---|---|
| `ust10y` | GS10 | GS10.csv |
| `fedfunds` | FEDFUNDS | FEDFUNDS.xlsx (annual avg) |
| `vix` | VIXCLS | VIXCLS.xlsx (annual avg) |

## 7. Laeven–Valencia banking crises — `16_banking.do`
Source: `SYSTEMIC_BANKING_CRISES_DATABASE_2026.xlsx` — Laeven & Valencia, *Systemic
Banking Crises Database*, 2026 update. Sheet `Crisis Resolution and Outcomes`
(header row 1, data from row 2; `A`=Country, `C`=Start, `D`=End). **164 episodes,
120 countries, 1976–2023.**

| Variable | Definition |
|---|---|
| `banking_crisis` | systemic banking-crisis dummy: **1 in every year of a crisis, 0 otherwise** |

**Construction.** The source is an *episode list* (one row per crisis, with start and
end years), not a country-year panel. Each `[Start, End]` window is expanded to
country-years and set to 1; **every other panel country-year is set to 0**. This
zero-fill is the substantive point: absence from the L-V list means *no systemic
banking crisis*, not *unknown*.

**Why this replaced the earlier GFDD extract.** `banking_crisis` previously came from
`Banking_crisis_dummy_data_set.xlsx` (World Bank GFDD `GFDD.OI.19`), whose year columns
are entirely missing from 2018 on. The variable therefore existed only 1990–2017, and
`l_banking_crisis = L.banking_crisis` only 1991–2018. Since `l_banking_crisis` is in `$ctrl_core`,
that one control truncated the whole estimation sample at 2018 (every `xtscc` table
printed `2019–2026 (omitted)`), and listwise deletion cut the headline LP from 61 onsets
to roughly 20 — excluding Argentina 2018, Turkey 2018, Lebanon 2019, Zambia 2020,
Sri Lanka 2022, Ghana 2022 and COVID. The old file is retained in `data/raw/` but is no
longer read by any do-file.

**Conventions and caveats.**
- *Ongoing crises.* A blank `End` means the crisis had not ended at publication, so the
  country stays coded 1 through the last panel year (Vietnam 2022–, Sri Lanka 2023–).
- *Country coverage.* 43 of our 52 countries match L-V by exact name; **Türkiye** is the
  only naming mismatch (mapped to `Turkey`). The other 8 — Belize, Guatemala, Honduras,
  Namibia, Pakistan, Serbia, South Africa, Trinidad and Tobago — carry no L-V crisis at
  all and are correctly 0 throughout. Within the 1994–2026 panel window, **33** countries
  record at least one crisis and **123** country-years are flagged; the remaining 11
  matched countries (Chile, Côte d'Ivoire, Egypt, El Salvador, India, Jordan, Morocco,
  Panama, Peru, Senegal, Tunisia) were last hit before 1994, so their in-window zeros are
  also genuine.
- *Post-2023 years.* L-V observed through 2023; zeros for 2024–2026 assume no *new*
  systemic banking crisis began in those years (ongoing crises are still carried forward).
- *Merge key.* Country × year (the workbook has no ISO3 column), following the same
  `merge m:1 country year` pattern as `13_ifs_nexus.do`.
- The `Resolution Details` sheet is **not** used: its date cells are corrupted (Argentina's
  2001 crisis reads `2026-11-01`). `Crisis Resolution and Outcomes` has clean integer years.

## 8. Derived predictors & transforms — `17_predictors.do`, `18_transforms.do`
All computed in code from the sourced series.

| Variable | Formula |
|---|---|
| `l_reg_crisis_share` (Z2) | lagged leave-one-out regional onset share |
| `past_onsets` / `past_def_onsets` (Z3) | lagged cumulative own (default-linked) onsets |
| `gdpg` | `100·(ln gdp_real − ln L.gdp_real)` |
| `l1_gdpg`, `l2_gdpg` | first/second lag of `gdpg` |
| `ln_gdp_base` | `L.ln_gdp` (total real GDP baseline) |
| `dy_0…dy_4` | `(F h.ln_gdp − ln_gdp_base)·100` — **headline LP outcome (total real GDP growth, Asonuma-aligned)** |
| `dy_pc_0…dy_pc_4` | per-capita version (robustness) |
| `dy_m1`, `dy_m2` | pre-trend placebos (on total real GDP) |
| `l_spr_mean`, `l_spr_max` | lagged spreads |
| `l_hyperinfl` | `L.infl > 50` — lagged hyperinflation flag, **in `$ctrl_core`** |
| `l_infl` | `L.infl` — raw lagged CPI inflation, % (built, not in core) |
| `l_lninfl` | `ln(1 + L.infl/100)` — lagged log gross inflation (built, not in core) |
| `cid`, `sample` | numeric country id; estimation-sample flag |

### How inflation enters the core
The core carries the **hyperinflation dummy**, matching Asonuma et al. Raw inflation cannot
be used directly: across the 52 panel countries 1994–2026 the median is 9.0% while
Venezuela 2018 is 65,374%, with 9 observations above 1000%. A raw term sits ~7,000× beyond
the median, so the coefficient would be set by a handful of country-years rather than the
relationship being controlled for. The dummy avoids this by flagging the hyperinflation
*regime* instead of its magnitude.

The flag is well populated in the raw data: **63 country-years exceed 50%**, across 17
countries (Venezuela 16, Turkey 11, Bulgaria / Argentina / Lebanon 4 each).

*A caveat to re-check after the next run.* In the runs preceding the Laeven–Valencia
rebuild, the first-stage probit reported `l_hyperinfl != 0 predicts failure perfectly;
1 obs not used` — only **one** estimation-sample observation had it equal to 1, so it was
dropped from the propensity model while still drawing a "significant" LP coefficient
(+3.54, p=0.026 at h=0) that flipped sign by h=3. That was a symptom of the old banking
control truncating the panel at 2018, which removed the entire recent hyperinflation
cluster (Argentina 2019–24, Lebanon 2020–23, Turkey 2022–24, Venezuela 2019–26). With the
sample extended past 2018 those years return, and the dummy should carry real variation —
but the probit output is worth checking for that message again.

`l_infl` (raw) and `l_lninfl` (log) are built alongside as continuous alternatives for a
dummy-vs-continuous robustness spec; swapping either into the core is a one-token edit in
the 14 `$ctrl_core` definitions. Prefer `l_lninfl` over `l_infl` for that — the log cuts the
Venezuela ratio from ~7,000× to ~75× and is safe because the minimum is −8.53%
(Azerbaijan 1999), well above the −100% at which `ln(1+x)` breaks.

**All `$ctrl_core` definitions must stay byte-identical** across `00_master.do`,
`18_transforms.do` and the 12 standalone guards — a divergence between them is what caused
the earlier `variable l_debt not found  r(111)` crash.

## 9. Legacy / source unrecovered
| Variable | Status |
|---|---|
| `imf` | IMF-program dummy carried over from the prior dataset (`10_skeleton.do`). 98% coverage, gap-free at onsets, but the **original source could not be recovered** (likely IMF MONA / History of Lending Arrangements; appears to include precautionary FCLs). The only non-rebuilt macro variable. |

## Not (yet) built
- **Capital-flow inflows** (portfolio / other / FDI, IMF BOP BPM6) — Asonuma Panels
  D/E channel outcomes; optional extension, not in the current analysis.
- **Nominal exchange-rate FX-regime dummies** (`ex_dum`) — Asonuma control; would
  need `PA.NUS.FCRF`. Only for full Asonuma parity.
