# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Topiary (pip package `topiary-asr`) is a Python framework for ancestral sequence
reconstruction (ASR). It wraps a pipeline of external phylogenetics tools (BLAST,
Muscle5, RAxML-NG, GeneRax, PastML, Open Tree of Life) behind a pandas-dataframe-based
API and a set of `topiary-*` CLI scripts, taking a user from a handful of seed
sequences to a full alignment, phylogenetic tree, and set of reconstructed ancestors
with an HTML report.

## Development setup

- Full install (creates external binaries, conda env, etc.) is `bash install.sh`.
  This is generally **not** needed for day-to-day development/testing — it does not
  require an NCBI API key or cluster access to work locally. It creates/updates a
  conda environment (default name `topiary`).
- For most development, the `topiary` conda environment should already be present
  at `~/miniconda3/envs/topiary/` with topiary itself installed via
  `pip install -e .`, plus `pytest`, `coverage`, `flake8`, `pytest-mock`, etc.
  Run Python/pytest/CLI commands using that environment's interpreter
  (`~/miniconda3/envs/topiary/bin/python`), not the system Python.
- `install.sh` additionally compiles custom-patched RAxML-NG and GeneRax
  (`dependencies/compile-raxml-ng.sh`, `dependencies/compile-generax.sh`) — only
  needed when those external tools themselves are being changed/rebuilt or on a
  fresh machine/cluster.

## Common commands

Run tests directly with pytest from the `topiary` conda env, e.g.:

```
pytest tests/topiary/quality
pytest tests/topiary/io/test_dataframes.py::test_specific_case
```

### Test tiers

Every test carries exactly one tier marker. Use these to pick how much you want
to run:

| tier | what it means | count | time |
| --- | --- | ---: | ---: |
| `unit` | hermetic: no network, no subprocess, no external binary | 245 | ~10s |
| `smoke` | real data and real file I/O, deterministic | 80 | ~25s |
| `integration` | needs an external binary or a live service | 51 | ~3min |

```
pytest tests -m unit           # fast feedback loop -- use this while iterating
pytest tests -m smoke
pytest tests                   # unit + smoke (integration is gated by flags)
```

How the tier is assigned:

- `integration` is **derived**, never written by hand — a test carrying any of
  the `--run-*`/`network` markers gets it automatically, so the two can never
  disagree.
- `smoke` is **explicit** (`@pytest.mark.smoke`), because "uses real data and
  takes a moment" is a judgement call.
- `unit` is the **default** for anything unmarked.

Because `unit` is the default, it is enforced rather than trusted: the
`block_subprocess` autouse fixture makes a `unit` test that shells out fail with
`SubprocessBlocked`. If you hit that, the test belongs in `smoke` (or behind a
`--run-*` flag if it needs a real binary) — don't remove the guard.

By default, pytest **skips** tests that hit external tools/servers. These are
opt-in via custom flags (see `tests/conftest.py`), and can be combined:

- `--run-raxml` — tests that shell out to RAxML-NG
- `--run-generax` — tests that shell out to GeneRax
- `--run-blast` — tests that run local/remote BLAST
- `--run-ncbi-server` — tests that hit the live NCBI server
- `--run-network` — tests that need live internet (Open Tree of Life, NCBI FTP,
  CDN fetches in report generation)

Tests that are **not** behind any opt-in flag have outbound connections blocked
by the `block_network` autouse fixture in `tests/conftest.py`. Tests carrying
`network`, `run_ncbi_server`, `run_blast`, `run_generax`, or `run_raxml` are
exempt (`_OPT_IN_MARKERS`) — the caller explicitly asked for them, so their
network use is declared rather than hidden.

The guard patches the parent process. On fork platforms (Linux) child processes
inherit it; on spawn platforms (macOS) they do not, so a test that reaches the
network only from a `multiprocessing` worker will pass locally on a Mac and fail
in Linux CI. Mock at the boundary rather than relying on the guard to catch it.

If a test starts failing with `NetworkAccessBlocked`, it has picked up a hidden
dependency on a remote service: either mock the call, or mark it
`@pytest.mark.network` (which makes it opt-in via `--run-network`). Do not "fix"
it by removing the guard.

Worth knowing: `run_generax` / `run_raxml` tests reach the Open Tree of Life API
even though they look purely local — setting up a GeneRax run annotates the
species tree from OTL. So `--run-generax` needs the internet, not just the
binary.

The working directory is likewise restored after every test by the `restore_cwd`
autouse fixture, so tests may `os.chdir` freely without a manual restore. Don't
add `current_dir = os.getcwd()` / `os.chdir(current_dir)` pairs back.

`--fd-report` tracks file descriptors opened and not closed by each test and
prints the worst offenders at the end of the session. It never fails a test; it
exists to locate the fd exhaustion seen when running the whole suite.

### CI

`.github/workflows/python-app.yml` is split by tier:

| job | what it needs | matrix | parallel? |
| --- | --- | --- | --- |
| `unit` | a pip install, nothing else | 2 OS x 3 python | yes |
| `smoke` | conda (muscle, blast) | 2 OS x 3 python | yes |
| `integration` | compiled RAxML-NG + GeneRax + live network | 1 config per push, full matrix nightly | no, serialized |

A manual `workflow_dispatch` run is treated exactly like the nightly —
full matrix, coverage, and the NCBI drift check — so a nightly-only failure
can be reproduced on demand instead of only at 07:00 UTC.

Only `integration` is serialized (`max-parallel: 1` plus a repo-wide
concurrency group) — it is the only job that reaches a live service (Open Tree
of Life), so it is the only one that can contend for one. That contention is why
the *whole* workflow used to be serialized; now it only costs the one job that
could cause it.

The `integration` job runs the entire suite rather than `-m integration`. That
is deliberate: it is the only job that exercises every test in a single process,
which is where cross-test contamination shows up (the file-descriptor
exhaustion in tests/BASELINE.md was invisible to any single test file).

`--run-ncbi-server` is **not run in CI at all**. Those tests query the live NCBI
Entrez server, which rate-limits unauthenticated traffic to ~3 requests/second
per IP; GitHub runners share IP space, so CI reliably gets HTTP 429. Serializing
does not help — the limit is per IP per unit time, not per concurrent job. The
logic is covered by `tests/topiary/ncbi/entrez/*_mock.py`. Run them against the
live server locally (`pytest tests --run-ncbi-server`); `NCBI_API_KEY` raises the
limit to ~10 req/s.

All Entrez call sites go through `topiary.ncbi.rate_limit()`, which waits only
for the remaining part of the interval since the last request. Add it before any
new Entrez call. (`ncbi/entrez/sequences.py` and `ncbi/blast/ncbi.py` are the
exception: they fan out over a multiprocessing pool and pass a lock down to the
workers, so they do their own throttling.)

**The one exception to "no NCBI in CI":** `tests/topiary/ncbi/entrez/test_ncbi_api_drift.py`
makes ~7 real calls to catch NCBI changing the *shape* of its responses, which
mocked tests cannot see. It runs nightly on a single matrix configuration only.
Its two rules keep it from reintroducing flakiness:

- It asserts on response **structure** — the specific keys topiary reads — never
  on values NCBI is entitled to change.
- A transport failure (429, 5xx, timeout, DNS) is a **skip**, not a failure.
  Only a response that parses differently than topiary expects fails.

If it fails, the message names the missing field and where topiary reads it.

`bash run_all_tests.sh` runs the full suite (flake8, test audit, coverage
with all of the above opt-in flags enabled, badge/report generation). This is slow
and touches the network/external binaries, so it should generally **not** be run by
Claude — but it's the reference for the exact test invocation syntax (including
which `--run-*` flags to pass) if that's needed.

Coverage settings live in `[tool.coverage.*]` in `pyproject.toml` —
`source_pkgs = ["topiary"]` and `branch = true`. Do not pass `--branch` or
`--source` on the command line; the config already handles it. Reports land in
`reports/coverage/` and `reports/htmlcov/`.

**Use `source_pkgs`, not `source = ["src/topiary"]`.** A directory only matches
when topiary is installed editable. CI installs non-editable, where the code
runs from site-packages, so a directory-based source traces nothing and reports
0% for everything. `tests/check_coverage_collected.py` guards against that
regression and runs before the floor check.

### Warnings

`filterwarnings = ["error", ...]` in `pyproject.toml`: **a warning that escapes a
test fails the test.** If a test expects a warning, assert on it —

```python
with pytest.warns(UserWarning,match="could not be assigned"):
    ott_list, species_list, results = species_to_ott(["Not a species"])
```

— which turns noise into a checked contract. If a *dependency* starts warning on
some platform, add a narrow ignore for that specific warning; don't delete the
`"error"` line.

Two ignores are in place, both for the same third-party pattern: ete4 parses
newick with `open(path).read()` and lets refcounting close the handle, which
raises `ResourceWarning` (and pytest's `PytestUnraisableExceptionWarning`
wrapper). It is noise, not a descriptor leak — fd counts stay flat. Use
`pytest --fd-report` to hunt real leaks; these filters don't hide those.

`run_all_tests.sh` enforces a **coverage floor of 92%** (`--fail-under=92`).
Ratchet it up as coverage improves; never lower it to make a change pass.

Test fixtures under `tests/data` must be **committed**. `.gitignore` has a
blanket `*.log` with a `!tests/data/**/*.log` negation for exactly this reason —
a test reading an uncommitted file passes locally and fails on a fresh checkout.
`run_all_tests.sh` fails if anything under `tests/data` is untracked. To check a
change the way CI sees it:

```bash
git checkout-index -a -f --prefix=/tmp/clean/ && cd /tmp/clean && pytest tests -m unit
```

Tests mirror the `src/topiary` package layout under `tests/topiary/` (e.g.
`src/topiary/quality/shrink.py` ↔ `tests/topiary/quality/test_shrink.py`).
`tests/audit_tests.py` reports tests that pass without verifying anything —
`pass`-bodied stubs, tests with no assertion, and tests silently shadowed by a
duplicate function name. All three are currently **zero**, and
`run_all_tests.sh` enforces that with `--max-stub 0 --max-noassert 0`, so a new
test that checks nothing fails the build. Run it as
`./tests/audit_tests.py tests`.

## Architecture

### Core data structure: the topiary dataframe

Nearly everything in topiary revolves around a single pandas dataframe format
(see `docs/source/data_structures.rst`). Required columns: `name`, `sequence`,
`species`. Topiary auto-adds `keep` (bool — sequences are never deleted, just
flagged `keep=False`), `uid` (unique 10-letter id, never edit by hand), and `ott`
(Open Tree of Life taxon id). An `alignment` column holds aligned sequences (all
same length) once alignment has occurred.

A **seed dataframe** (`name`, `species`, `aliases`, `sequence`) is the small,
human-curated starting point a user provides; `topiary.io.df_from_seed` expands it
into a full topiary dataframe via BLAST.

### API convention: functions take a df, return a modified copy

Most public topiary functions have the signature `df = topiary.something(df, ...)`,
returning a modified copy (new/changed columns, or `keep` flipped to `False`) rather
than mutating in place. This is what lets each pipeline step be run, written to CSV,
manually inspected/edited, and resumed independently. See `src/topiary/__init__.py`
for the full set of top-level re-exported functions and which submodule each comes
from (`quality`, `muscle`, `opentree`, `ncbi.blast`, `raxml`, `generax`, `io`).

### Package layout (`src/topiary/`)

- `ncbi/` — BLAST (local + NCBI remote) and NCBI Entrez access, reciprocal BLAST
  orthology calls.
- `opentree/` — Open Tree of Life species-tree lookups/resolution.
- `muscle/` — Muscle5 alignment wrapper.
- `quality/` — dataset size/redundancy reduction (`shrink_dataset`) and alignment
  polishing.
- `raxml/` — RAxML-NG wrappers: model selection, ML gene tree, ancestor inference,
  bootstraps.
- `generax/` — GeneRax wrapper for gene/species tree reconciliation.
- `pastml/` — gap-state ancestral inference (parsimony).
- `io/` — reading/writing topiary dataframes and alignments, paralog pattern
  compilation.
- `pipeline/` — the three top-level pipelines that chain the above modules together
  (see below).
- `reports/` — generates the self-contained HTML result reports (`reports/cards`
  for individual result cards, `reports/assets` for static assets bundled into
  output).
- `draw/` — tree/ancestor plotting (toytree/matplotlib based).
- `cli_scripts/` — thin argument-parsing wrappers around `pipeline/` functions;
  these are the `topiary-*` console entry points defined in `pyproject.toml`.
- `_private/` — internals not part of the public API, notably
  `_private/supervisor.py`'s `Supervisor` class, which manages the standardized
  `input/ working/ output/ run_parameters.json` run-directory layout used by every
  pipeline step that wraps external software. Also houses environment checks
  (`installed.py`, `environment.py`) and threading utilities. (topiary no longer
  uses MPI itself; GeneRax parallelism, if wanted, is invoked via a user-supplied
  launcher prefix — see below.)

### The three pipelines (`src/topiary/pipeline/`, exposed as `topiary-*` CLI scripts)

1. **seed_to_alignment** (`topiary-seed-to-alignment`): seed sequences → BLAST →
   reciprocal BLAST orthology calls → OTT/species-tree filtering → redundancy
   reduction → Muscle5 alignment → alignment polishing. Writes a numbered sequence
   of CSVs (`01_initial-dataframe.csv` … `06_alignment.fasta`) so each step's effect
   is inspectable/reversible. Meant to be followed by manual alignment inspection/
   editing (e.g. in AliView) before the next pipeline.
2. **alignment_to_ancestors** (`topiary-alignment-to-ancestors`): final alignment →
   RAxML-NG model selection → ML gene tree → ancestor inference (RAxML-NG +
   PastML for gaps) → GeneRax gene/species tree reconciliation (skipped for
   microbial-only datasets) → reconciled-tree ancestors → gene-tree bootstrap
   replicates. Typically run on a cluster; produces the `results/index.html` report.
3. **bootstrap_reconcile** (`topiary-bootstrap-reconcile`): takes the bootstrap
   replicates from step 2 and reruns GeneRax reconciliation on each one to compute
   branch supports on the reconciled tree. Computationally heavy, split out as its
   own step. It is a **re-entrant filesystem crawler** (`generax/_crawl.py`), not
   an MPI job: launch as many copies as you like (e.g. a SLURM job array, one task
   per node) pointing at the same run directory, and they coordinate through marker
   files — one builds the replicate directories, all of them claim-and-run
   replicates, one aggregates the supports. Re-running resumes (completed replicates
   are skipped; interrupted ones resume from GeneRax's own on-disk checkpoint). Each
   generax replicate runs single-core by default; to parallelize an individual
   reconciliation, pass a launcher prefix via `--generax-launch "mpirun -np N"`
   (GeneRax parallelism is MPI-only, and the user owns the launcher/allocation).

`seed_to_alignment` and `alignment_to_ancestors` support `--restart` to resume from
the last completed step in an existing output directory (state is tracked via
`Supervisor`/`run_parameters.json`); `bootstrap_reconcile` resumes automatically by
being re-run (no flag).

### Run directories

Any pipeline step that wraps external software (everything after seed-to-alignment)
writes output using the same directory shape: `input/`, `working/`, `output/`,
`run_parameters.json`, managed by `_private/supervisor.Supervisor`. This is what
makes `--restart` possible and keeps intermediate files inspectable.
