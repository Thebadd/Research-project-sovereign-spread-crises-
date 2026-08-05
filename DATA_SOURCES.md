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

## 2. IMF World Economic Outlook — `11_weo.do`
Source: `data/raw/WEOApr2026all.xlsx`, sheet "Countries", by `INDICATOR.ID`.

| Variable | WEO code | Definition |
|---|---|---|
| `gdppc_real` / `ln_gdppc` | NGDPRPC | real GDP per capita (const. prices); `ln_gdppc = ln(NGDPRPC)` |
| `gdp_real` | NGDP_R | real GDP (const. prices) — used for `gdpg` |
| `gdp_nominal` | NGDP | nominal GDP (for ratios) |
| `pop` | LP | population (millions) |
| `infl` | PCPIPCH | CPI inflation, period avg, % |
| `ca` | BCA_NGDPD | current account, % GDP |
| `debt` | GGXWDG_NGDP | general govt gross debt, % GDP |
| `govexp` | GGX_NGDP | general govt expenditure, % GDP |
| `revenue_gdp` | GGR_NGDP | general govt revenue, % GDP |
| `inv` | NID_NGDP | total investment, % GDP |
| `pb` | GGXONLB_NGDP | primary balance, % GDP |

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
Source: `Banking_crisis_dummy_data_set.xlsx` (World Bank GFDD `GFDD.OI.19`).

| Variable | Definition |
|---|---|
| `banking_crisis` | systemic banking-crisis dummy (0/1). Kept as provided — no duration. |

## 8. Derived predictors & transforms — `17_predictors.do`, `18_transforms.do`
All computed in code from the sourced series.

| Variable | Formula |
|---|---|
| `l_reg_crisis_share` (Z2) | lagged leave-one-out regional onset share |
| `past_onsets` / `past_def_onsets` (Z3) | lagged cumulative own (default-linked) onsets |
| `gdpg` | `100·(ln gdp_real − ln L.gdp_real)` |
| `l1_gdpg`, `l2_gdpg` | first/second lag of `gdpg` |
| `ln_gdppc_base` | `L.ln_gdppc` |
| `dy_0…dy_4` | `(F h.ln_gdppc − ln_gdppc_base)·100` (LP outcome) |
| `dy_m1`, `dy_m2` | pre-trend placebos |
| `l_spr_mean`, `l_spr_max` | lagged spreads |
| `cid`, `sample` | numeric country id; estimation-sample flag |

## 9. Legacy / source unrecovered
| Variable | Status |
|---|---|
| `imf` | IMF-program dummy carried over from the prior dataset (`10_skeleton.do`). 98% coverage, gap-free at onsets, but the **original source could not be recovered** (likely IMF MONA / History of Lending Arrangements; appears to include precautionary FCLs). The only non-rebuilt macro variable. |

## Not (yet) built
- **Capital-flow inflows** (portfolio / other / FDI, IMF BOP BPM6) — Asonuma Panels
  D/E channel outcomes; optional extension, not in the current analysis.
- **Nominal exchange-rate FX-regime dummies** (`ex_dum`) — Asonuma control; would
  need `PA.NUS.FCRF`. Only for full Asonuma parity.
