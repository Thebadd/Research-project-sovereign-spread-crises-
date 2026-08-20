# Data (paper-ready draft)

> Continuous-prose draft of the Data section, written to sit immediately
> before the Results section drafted in `RESULTS_SECTION_DRAFT.md`.
> `DATA_SOURCES.md` is the full variable-by-variable provenance record —
> source file, series code and transform for every column — and
> `METHODOLOGY.md` carries the estimator and control-set detail. This
> section states what the data are, how a crisis is identified, and what
> the two groups look like before anything is estimated.

## The panel

The estimation sample is an annual, country-level panel of 52 emerging and
frontier economies running from 1994 to 2026. Coverage begins when a country
enters the J.P. Morgan EMBI Global index rather than at a common date, so the
panel is unbalanced on the left — Armenia enters around 2013, Kenya and
Zambia around 2014 — and the sample is defined by index membership rather
than by an income or region cutoff. That is the appropriate frame for this
question: the object of interest is the behaviour of a market-priced
sovereign spread, and a country with no traded external benchmark has no
spread to cross a threshold. It also fixes the external validity of
everything that follows. These are results about sovereigns that borrow in
international markets and are watched by them, not about low-income
sovereigns whose external liabilities are official and unpriced.

Macro variables are rebuilt from official sources for every country-year
rather than inherited: the IMF *World Economic Outlook* (April 2026 vintage)
supplies real GDP, the current account, general government debt, expenditure
and the primary balance, investment, inflation and population; World Bank
WDI supplies domestic credit to the private sector by banks, claims on
central government, trade openness and FDI; the sovereign-bank exposure
measures are built from the IMF's Monetary and Financial Statistics as shares
of aggregate deposit-taking-bank assets; and banking-crisis dates come from
the Laeven–Valencia systemic banking crises database. Spreads and the crisis
classification come from the author's own EMBIG-based spread-crisis database.
The full mapping is in `DATA_SOURCES.md`.

The exposure measures carry a coverage restriction that binds later: the IMF
depository-corporations file begins in 2001 and omits six of the 52 countries
entirely, so the sovereign-bank results in the Results section are estimated
on a strictly smaller panel than the output results, and the difference is
not random with respect to country size.

Two features of the panel construction have direct consequences for the
estimates and are therefore stated here rather than in a footnote. First,
because the panel starts at index entry rather than at the start of macro
coverage, the lag operators used to build the controls initially fell off the
left edge of each country's series: a country entering in 2013 with a 2013
onset had no lagged growth rate, even though the underlying WEO series runs
from 1980. Three scaffolding rows are appended before each country's first
panel year so the lags can be constructed from data already on disk. These
rows are flagged, carry no spread quote, are never treated, and are excluded
from every estimation sample; they add no source and change no specification.
Second, and less benignly, the WEO vintage runs to 2026 while the last
realised national-accounts observation is 2024 or 2025 for most countries.
Any onset from 2022 onward therefore has its fourth- and fifth-year outcome
built partly from IMF projections. Seven and nine episodes respectively are
affected, at precisely the horizons where the estimated cost is largest, and
the Results section reports what happens when the outcome is restricted to
realised data.

## Identifying a spread crisis

An episode is a spread crisis when the country's EMBIG spread crosses either
of two thresholds: an absolute one at 1,000 basis points, and a
country-specific one at the 99th percentile of that country's own historical
spread distribution. The absolute criterion alone would over-select
persistently high-spread sovereigns and under-select those whose spreads are
normally tight; the relative criterion alone would generate an onset in every
country's worst year regardless of whether the level was economically
meaningful. Requiring either, rather than both, keeps episodes that are
severe in level and episodes that are severe relative to the country's own
history, which is the union the definition is meant to capture.

This yields 61 onsets across the 52 countries. Each episode is classified by
how it was resolved: 40 are **non-default**, in which the spread crossed the
threshold and the sovereign continued to service its external debt, and 21
are **default-linked**, in which the episode coincided with or was followed
by a default or distressed restructuring. That split is the paper's central
partition and it is not a partition of severity — the classification is made
on the resolution of the episode, not on the height of the spread — which is
what allows the comparison in the Results section to be read as being about
resolution rather than about how bad the crisis was to begin with.

Treatment is defined narrowly: each episode contributes exactly one treated
observation, its onset year, and continuation years are dropped from the
estimation sample entirely rather than being counted as either treated or
untreated. The control group is therefore tranquil years — country-years in
which no episode is under way — and the multi-year dynamics of an episode are
traced through the horizon of the outcome rather than through additional
treated rows. The cost of this convention is that a five-year crisis
identifies the coefficient exactly as much as a one-year crisis, so episode
duration is never a source of variation and the paper cannot speak to whether
protracted crises are more costly than brief ones.

## What the raw data show

Before any specification is imposed, Figure 0 plots the cumulative change in
log real GDP around an onset, separately for the two resolution types. The
only operation performed on the data is country-demeaning: countries differ
systematically in trend growth, so an unadjusted average across episodes
would partly rank countries rather than describe crises, and subtracting each
country's own mean over the estimation sample removes exactly the variation a
country fixed effect removes in the regressions. Nothing else is done — no
controls, no year effects, no weighting, no estimator — and no confidence
bands are shown, because attaching inference to a description invites reading
it as a result.

| Year | All | Non-default | Default-linked |
|:--:|:--:|:--:|:--:|
| −1 | 0.29 | −0.92 | 2.49 |
| 0 | 0 (base) | 0 (base) | 0 (base) |
| 1 | −2.69 | −0.80 | −6.11 |
| 2 | −4.84 | −2.65 | −8.82 |
| 3 | −3.27 | −1.29 | −6.75 |
| 4 | −2.44 | −0.25 | −6.29 |
| 5 | −2.40 | −0.11 | −6.31 |

*Country-demeaned means, percentage points. n = 38 non-default and 21
default-linked onsets at Years 1–2, falling to 36 and 21 by Year 5. Year 0 is
the pre-crisis baseline and is zero by construction; Year −1 is the single
estimable pre-crisis point, measured on the same $t-1$ base.*

Two things are worth reading off this table, and one of them is a warning
rather than a preview. The divergence that the Results section estimates is
already visible without any conditioning: default-linked episodes lose six to
nine percentage points of output and do not recover within five years, while
non-default episodes lose under three points at their worst and are back to
their country's normal path by Year 4. The raw gap is in fact *larger* than
the conditional estimates reported later, so the controls attenuate the
finding rather than producing it — the regressions are sharpening a fact, not
manufacturing one. Note also that the descriptive table uses all 21
default-linked episodes, whereas listwise deletion on the control set leaves
15 of them identifying the headline coefficients; part of the difference in
magnitude is a difference in sample, not only in specification.

The warning is the Year −1 row. Default-linked episodes sit 2.5 points above
their country's normal path two years before onset and non-default episodes
sit 0.9 points below it. On the construction used here — output at $t-2$
minus output at the $t-1$ base — a positive value means output was *higher*
before the base year, that is, a pre-crisis slowdown. The two groups were
already on diverging growth paths going in. The Results section takes this up
directly as a placebo test and stress-tests the headline result against it.

## Pre-crisis characteristics

Table 1 reports the mean of each control, measured in the year before onset,
for tranquil country-years and for each resolution type, with a two-sample
test of the non-default/default-linked difference. Everything is measured at
$t-1$, so the table describes precisely the variation the controls are asked
to absorb.

| Variable ($t-1$) | Tranquil | Non-default | Default-linked | $p$ (nd = def) |
|:--|--:|--:|--:|--:|
| Real GDP growth (%) | 3.74 | 4.62 | 1.21 | **0.001** |
| Public debt (% GDP) | 49.90 | 50.94 | 60.79 | 0.308 |
| Current account (% GDP) | −1.79 | −4.44 | −4.12 | 0.861 |
| Banking-crisis duration (yrs) | 0.14 | 0.18 | 0.52 | 0.192 |
| Government expenditure (% GDP) | 27.91 | 26.20 | 28.16 | 0.446 |
| Trade openness (% GDP) | 72.96 | 62.07 | 65.51 | 0.707 |
| Bank credit to private sector (% GDP) | 47.06 | 37.76 | 37.20 | 0.945 |
| Hyperinflation flag | 0.01 | 0.03 | 0.10 | 0.292 |
| EMBIG spread, mean (bps) | 326.44 | 403.70 | 592.02 | **0.001** |
| Bank claims on government (% assets) | 13.56 | 15.53 | 12.56 | 0.554 |

The result of interest here is how little separates the two groups. Of the
ten pre-crisis characteristics, only two differ significantly: growth in the
year before onset, and the level of the spread. Public debt, the current
account, banking-sector distress, government spending, trade openness, credit
depth, hyperinflation and — notably — bank exposure to the sovereign are all
statistically indistinguishable between episodes that ended in default and
episodes that did not.

This should not be read as a balance test that the design has passed. Two
significant differences out of ten is not nothing, and both of them matter:
default-linked countries entered their crises with a third of the growth
momentum of non-default ones (1.21% against 4.62%) and with spreads roughly
190 basis points higher. The first is the pre-trend already visible in the
Year −1 row above, and the reason one lag of GDP growth is a non-negotiable
element of the control set rather than an optional addition. The second says
default-linked episodes were more severe as market events at the outset, so
the resolution comparison is not a clean comparison of two equally distressed
groups.

What the table does establish is the *shape* of the selection problem, which
is narrower than one might have expected. The two groups are not drawn from
different populations of fundamentals; they differ on the growth path they
were already on and on how severely the market had already priced them. That
is exactly the kind of selection a propensity-score design can address, and
it is why the Results section moves from the single-stage projections to
inverse-probability weighting and then to a doubly-robust estimator rather
than resting on the conditional projections alone. Read this table alongside
the first-stage probit reported there, not as a hurdle the design has to
clear before the estimates can be believed.

One further observation belongs here because it constrains a later section
rather than this one. Bank claims on government as a share of bank assets —
the sovereign-bank exposure measure the paper uses to test for amplification
— is if anything *lower* in default-linked episodes than in non-default ones
(12.6% against 15.5%), and the difference is nowhere near significant.
Whatever the sovereign-bank nexus does in these data, it does not operate by
default-linked countries having entered their crises with visibly more
exposed banking systems. It has to work through the interaction of exposure
with resolution, which is how the Results section tests it.
