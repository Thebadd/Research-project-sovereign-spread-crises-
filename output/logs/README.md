# Run logs

Two kinds of file live here.

**`YYYY-MM-DD_master.log`** — the full console output of a `00_master.do` run,
written automatically by the logging block at the top of the master. This is
the complete record: every regression, every coverage report, every warning.
Plain text so it greps and diffs.

**`YYYY-MM-DD_<file>.md`** — curated extracts, written by hand when a single
file was re-run outside the master and its numbers are cited in the write-up.
These keep the estimate tables and drop the code echo, which is already in git.

## Why these exist

Numbers in `EMPIRICAL_ANALYSIS.md` are quoted to three decimals. Until the
master logged its own output, the only record of where a given number came
from was whatever had been copied out of the Results window, and when an
estimate moved between runs there was no way to tell whether the code or the
session had changed. That ambiguity is what motivated seeding the bootstraps
(see the note in `13d_aipw_nexus_split.do`); logging is the other half of the
same fix.

## Current contents

- `2026-08-19_13_mechanisms.md` — Test 1 (credit supply vs demand), Test 3
  (current account, with the Clogg z added that day). Cited by §7b.
- `2026-08-19_13d_nexus_nboot1000.md` — the 1000-draw seeded nexus run, plus
  the comparison of the four contested cells across all three bootstrap runs.
  Cited by §8, and supersedes the earlier 300-draw runs.
