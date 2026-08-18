# Quantum Pilot

Two sibling, local-only quantum chemistry orchestration tools, plus a
small landing page that links to both:

- **[gaussbot](gaussbot/)** — Gaussian workflows: guess geometry →
  resilient PM6 pre-optimization → final-level reoptimization →
  reaction-mechanism study (TS search, IRC verification, energetics),
  or just a single geometry optimization.
- **[gamessbot](gamessbot/)** — GAMESS multireference workflows: RHF →
  CIS → CASSCF → XMCQDPT / TRANSITN (oscillator strengths), starting
  from an already-optimized geometry (including one gaussbot just
  produced).
- **[home/](home/)** — a static landing page ("Quantum Pilot") with
  one card per bot, linking to whichever local port each is actually
  running on.

Neither tool does any quantum chemistry itself — they build input
files, run Gaussian/GAMESS locally (or submit via PBS), and parse the
logs that come back. Gaussian/GAMESS do the actual work; this is glue
code, plus a GUI and CLI around it.

## Quick start

```bash
./install.sh
```

Interactive by default: pick gaussbot, gamessbot, or both; it looks
for a Gaussian environment already in your `~/.bashrc` (offers to use
it if found, otherwise asks); for gamessbot it also asks for the path
to `rungms` and a GAMESS scratch directory. Each bot gets its own venv
under `<bot>/venv/`, plus a `<bot>/.env` file holding whatever it
collected (nothing is written to your actual `~/.bashrc`).

Non-interactive (e.g. scripted/CI use) — set env vars and pass `--yes`:

```bash
BOT_INSTALL_CHOICE=3 \
BOT_G09ROOT=/home/you \
BOT_GAUSS_SCRDIR=/scratch \
BOT_RUNGMS_PATH=/path/to/gamess/rungms \
BOT_GAMESS_SCRATCH_DIR=/path/to/gamess_scratch \
./install.sh --yes
```

`BOT_INSTALL_CHOICE`: `1` = gaussbot only, `2` = gamessbot only, `3` =
both (default). See the comment block at the top of `install.sh` for
the full list of overrides.

The installer is safe to copy elsewhere and rerun — it resolves its
own location at runtime and never bakes an absolute path into anything
except the `.env` files it writes for you.

## Running a bot

Each bot has a portable launcher that sources its own `.env` and
starts its GUI:

```bash
gaussbot/run_gui.sh     # http://127.0.0.1:8765 by default
gamessbot/run_gui.sh    # http://127.0.0.1:8766 by default
```

Or, with the venv active, the CLI:

```bash
source gaussbot/venv/bin/activate && gaussbot --help
source gamessbot/venv/bin/activate && gamessbot --help
```

## The landing page

```bash
python3 home/serve.py   # http://127.0.0.1:8764
```

Opens a page with a card per bot. Each card's link is resolved live
(via `home/serve.py`'s `/config` endpoint reading `GAUSSBOT_GUI_PORT`/
`GAMESSBOT_GUI_PORT` out of each bot's own `.env`), so it stays correct
even if you change a bot's port after install — a card only works while
that bot's own server is actually running.

## `.env` files

Written once by `install.sh`, editable by hand any time afterward — no
reinstall needed for a path/port change to take effect, just restart
that bot's `run_gui.sh`.

- **`gaussbot/.env`**: `g09root`, `GAUSS_EXEDIR`, `GAUSS_SCRDIR`,
  optionally `GAUSSBOT_GUI_PORT`.
- **`gamessbot/.env`**: the same Gaussian variables (needed for its
  optional "guess geometry, optimize with Gaussian first" intake path),
  plus `GAMESSBOT_RUNGMS_PATH`, `GAMESSBOT_SCRATCH_DIR`, optionally
  `GAMESSBOT_GUI_PORT`.

## More detail

`gaussbot/README.md` has the fuller picture for that package —
pipeline stages, job types, the Python API, system dependencies
(Open Babel, `freqchk`, an optional `ANTHROPIC_API_KEY` for AI-assisted
structure refinement). gamessbot doesn't have its own README yet;
`gamessbot --help` and the in-GUI Tips panel cover the day-to-day
questions (active-space selection, convergence recovery, cost
tradeoffs).
