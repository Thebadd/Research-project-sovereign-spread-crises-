# Introduction (paper-ready draft)

> Continuous-prose draft of the Introduction, written to the structure of
> Asonuma, Chamon, Erce & Sasahara (2024) — motivation and gap, design,
> results, mechanism, heterogeneity, then placement in the literature — but
> on this paper's own narrative rather than theirs. Every magnitude below is
> traceable to `EMPIRICAL_ANALYSIS.md`; the two paragraphs of country
> examples in the opening should be checked against `Episode_Summary` in the
> crisis database before submission, since the classification there is
> authoritative and this draft names episodes from memory of the sample.

Between 2020 and 2023 a series of emerging and frontier sovereigns saw their
external borrowing costs rise past levels that markets treat as pricing
default. Some of them defaulted: Argentina, Ecuador, Lebanon, Zambia, Sri
Lanka and Ghana all restructured or missed payments. Others crossed the same
thresholds and did not, servicing their debt through the episode and
eventually seeing spreads compress again. The literature on the cost of
sovereign debt crises has a great deal to say about the first group and
almost nothing to say about the comparison. It conditions on the credit event
— a default, a restructuring, a missed payment — and estimates what follows
it, typically against a control group of countries that were not in distress
at all. That design answers a question about defaults. It cannot answer the
question a finance minister facing a spread of a thousand basis points
actually has, which is not what happens after a default but what is gained by
avoiding one.

This paper takes the market event rather than the credit event as the unit of
analysis. We identify sovereign spread crises directly from secondary-market
pricing, classify each episode by how it was ultimately resolved, and ask
whether the output cost of a spread crisis is a consequence of the distress
itself or of the credit event that some episodes end in. The comparison group
for a default is therefore not a tranquil country but a country that reached
the same point of market distress and did not default. Holding the distress
fixed and varying only the resolution is what the design is for, and it is
what separates this paper from the literature it builds on.

Our sample covers 52 emerging and frontier economies over 1994–2026 and 61
spread-crisis episodes, of which 40 were resolved without default and 21
coincided with or were followed by default or distressed restructuring. An
episode begins when a country's J.P. Morgan EMBI Global spread crosses either
of two thresholds: a level threshold of a thousand basis points, following
Pescatori and Sy (2007), or a repricing threshold of 378 basis points on the
quarter-on-quarter change, which Horn, Reinhart and Trebesch (2021) calibrate
as the ninety-ninth percentile of pooled quarterly spread changes across 45
economies. Requiring either rather than both captures both sustained distress
and sudden repricing; a level criterion alone would select persistently
high-spread sovereigns and miss the violent repricing of countries whose
spreads are normally tight. Episodes are dated following Detragiache and
Spilimbergo (2001), which bridges a single non-crisis year so that one
prolonged distress is not fragmented into several artificial episodes. We
estimate impulse responses by local projections (Jordà, 2005), with tranquil
country-years as controls.

The pooled response is modest. Averaged across all 61 episodes, a spread
crisis is followed by a cumulative output shortfall of two to three
percentage points at the one- and two-year horizons, weakening to marginal
significance by the third year and indistinguishable from zero thereafter. On
its own that number would support a reassuring reading of market distress as
a transitory event. Splitting the sample by resolution shows that it does not
summarise the sample so much as conceal it. Episodes resolved without default
carry no reliable output cost at any horizon, with a coefficient that is
small, inconsistently signed and significant at one horizon out of five.
Default-linked episodes carry a cost of between four and seven percentage
points, significant at every horizon out to five years, and the formal test of
equality between the two rejects at four of five horizons. The output cost of
a sovereign debt crisis, on this evidence, is not a generic consequence of
crossing into distressed pricing. It is concentrated in, and largely specific
to, the episodes that end in a credit event.

Two features of the data qualify that finding and we state them at the outset
rather than in a limitations section. The first is a pre-trend. A placebo test
on pre-crisis growth returns a significant coefficient for default-linked
episodes and none for non-default ones, and its sign says those countries were
already decelerating going into their crises: growth in the year before onset
averages 1.2 per cent for default-linked episodes against 4.6 per cent for
non-default ones. A group already slowing will post weaker cumulative growth
afterwards whether or not the resolution matters. One lag of growth is
therefore in the control set throughout, and adding a second, non-overlapping
lag moves the default-linked coefficient by no more than a tenth of a
percentage point at any horizon, which addresses the version of the concern
that reaches further back than the baseline can see without disposing of the
concern itself. The second is selection. Because a country's ability to avoid
default depends on the conditions it faces, ordinary least squares is not
credible on its own, and we follow the convention of the local-projection
literature in reweighting on an estimated propensity score (Jordà & Taylor,
2016) and then in estimating the doubly-robust augmented inverse-probability
weighted analog, which requires only one of the propensity and outcome models
to be correctly specified. Under AIPW, with the difference between resolution
types bootstrapped directly on paired draws, the extra cost of default is
significant at the first three horizons and is not distinguishable from zero
at the fourth and fifth. A third check points the same way for an unrelated
reason: the panel runs to 2026 while realised national accounts end in 2024 or
2025 for most countries, so episodes from 2022 onward have their fourth- and
fifth-year outcomes built partly from projections, and restricting the outcome
to realised data nearly halves the default-linked coefficient at exactly those
horizons. Three independent checks converging on the same boundary is what
leads us to claim the resolution gap for Years 1 through 3 and to report, but
not build on, the longer horizons.

Turning to transmission, we re-estimate the same projection on six
intermediate outcomes and find private credit to be the clearest responder,
contracting by five to seven percentage points by the third and fourth years
with a delayed, compounding profile that matches the shape of the output
response itself. To move beyond parallel single-outcome regressions, which
establish co-movement and not mediation, we compute how much of the
default-linked output coefficient each channel statistically absorbs when
added to the headline specification on the identical sample, in the spirit of
Gelbach (2016). Private credit and government expenditure each account for
roughly forty per cent of the coefficient by the fourth and fifth years; the
primary balance and foreign direct investment account for essentially none of
it at any horizon. This is a ranking among candidate channels rather than a
decomposition into causal shares, and we are explicit that a coefficient which
shrinks when a channel is controlled for is equally consistent with reverse
causality, but it is a ranking the resolution-split regressions could not
deliver.

The paper's sharpest mechanism evidence concerns the sovereign-bank nexus, and
here we depart from the reference literature on measurement as well as on
question. Asonuma et al. proxy the role of banks by the size of the banking
system, private credit as a share of GDP, and split their sample at its
median. We find that measure to be uninterpretable in our panel: the
interaction between crisis onset and pre-crisis credit-to-GDP comes back
robustly positive, as do the equivalent interactions for investment and FDI,
and our own outcome equation shows high-credit countries growing more slowly
on trend and losing less in a crisis — the signature of an income and
development proxy rather than of exposure to a channel. We therefore measure
the nexus as banks' claims on the general government as a share of total bank
assets, a portfolio composition share that does not scale with financial
development and that is predetermined with respect to the crisis it is asked
to amplify. On that measure, one standard deviation of additional pre-crisis
exposure is associated with roughly four and a third additional percentage
points of lost output by the third and fourth years, and only in
default-linked episodes: the non-default interaction sits indistinguishably at
zero throughout and the equality test between the two rejects at both
horizons. The behavioural counterpart is visible directly. Where banks entered
the crisis heavily exposed and the crisis ended in default, their claims on
government rise by eight to seventeen percentage points of assets over three
consecutive years while their claims on the private sector fall, with no
comparable movement in any other cell. Notably, exposure does not differ
between the two resolution groups before their crises begin, so what the
evidence describes is an interaction between exposure and resolution rather
than defaulters simply having had more exposed banking systems to start with.

One result we report because it constrains the others is a null. Aguiar and
Gopinath (2006) imply that losing market access forces the current account
toward surplus, and it does: the external balance improves by one percentage
point of output in the crisis year and rises monotonically to two and a half
points by the fifth. But it does so identically in both groups, with no
horizon coming close to rejecting equality. Forced external adjustment appears
to be a generic consequence of losing market access rather than one of the
margins on which the default premium operates. Since output, credit, fiscal
policy and the sovereign-bank linkage all differ sharply by resolution, a
channel that does not is informative about where the default premium is not
located.

Our results speak to three literatures. They add to the work on the output
costs of sovereign default — Sturzenegger (2004), Tomz and Wright (2007),
Borensztein and Panizza (2009), Levy-Yeyati and Panizza (2011), Asonuma and
Trebesch (2016), Kuvshinov and Zimmermann (2019), Farah-Yacoub, Graf von
Luckner and Reinhart (2024) — by supplying the counterfactual that literature
mostly lacks. Where those papers estimate what follows a default relative to
countries not in distress, we estimate it relative to countries that reached
the same distressed pricing and avoided the credit event, and we find the cost
survives that far more demanding comparison. Our pre-crisis finding also
connects to Levy-Yeyati and Panizza's observation that output contraction
largely precedes the default that is supposed to have caused it, which in our
sample appears as a pre-crisis deceleration specific to the group that goes on
to default. Second, they relate to the theoretical and empirical work on
banks' sovereign exposure — Gennaioli, Martin and Rossi (2014), Bocola (2016),
Sosa-Padilla (2018), Broner, Erce, Martin and Ventura (2014), Erce and
Mallucci (2018) — by showing that direct balance-sheet exposure, and not the
size of the banking system, is the variable that conditions the output cost,
and by documenting the reallocation from private to sovereign claims that
accompanies it. Third, they contribute methodologically to the applied
local-projection literature (Jordà, 2005; Jordà and Taylor, 2016; Auerbach and
Gorodnichenko, 2016) by treating a binary market event with a genuine
zero-treatment control group, and by reporting, at each rung of the estimator
ladder, which conclusions survive the move from ordinary least squares to
reweighting to doubly-robust estimation and which do not.

Two limits bound everything that follows. Listwise deletion on the control set
leaves 31 non-default and 15 default-linked episodes identifying the headline
coefficients, so the precision of the resolution comparison rests on fifteen
treated observations, and the heterogeneity results rest on considerably
fewer; magnitudes should be read accordingly, and we report the within-sample
counts alongside every estimate rather than only the episode totals. And
throughout, the design identifies the association between how a crisis is
resolved and what follows it, conditional on observables and after correcting
for selection on those observables. It does not identify the effect of
choosing to default, and no estimate here should be read as one.

The remainder of the paper proceeds as follows. Section 2 describes the panel,
the identification of spread crises and their classification, and reports the
raw output paths and pre-crisis characteristics of the two groups. Section 3
sets out the estimator ladder and the reasoning behind each rung. Section 4
reports the pooled and by-resolution output responses and the robustness
checks that bound them. Section 5 turns to transmission channels and their
decomposition. Section 6 examines the sovereign-bank nexus. Section 7
concludes.
