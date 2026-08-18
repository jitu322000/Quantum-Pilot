# gaussbot

A small orchestration layer around Gaussian for reaction mechanism
studies -- or just a single geometry optimization, when that's all
you need. For a mechanism study: guess geometry → resilient PM6
pre-optimization → user-chosen final-level reoptimization (plain
basis set or a mixed GenECP one) → TS guess + search (with a check
that the imaginary mode found is actually the reaction, not something
spurious) → optional IRC verification → barrier extraction. All
stages are built and run end to end locally, both against toy
molecules and against a real DFT/GenECP organometallic TS input.

This package does no quantum chemistry itself. It builds `.com`
files, will hand them to Gaussian via PBS, and will parse the logs
that come back. Gaussian does the actual work; this is glue code.

## Status

Built and tested against real `g09` runs (not just mocked):
- `gaussbot/geometry.py` — SMILES / file / PubChem → `Structure` (3D coords + charge + multiplicity); also `displace_structure()`, for nudging a geometry along a vibrational mode, and `rmsd()` (Kabsch-align + RMSD), the "how close are these two independently-obtained geometries" check used throughout the pipeline. `from_chemdraw()` (new): loads a `.cdx`/`.cdxml` ChemDraw file via `obabel -i{cdx|cdxml} <path> -omol -O <tmp>.mol` then `Chem.MolFromMolFile()` on the result -- deliberately *not* trusting the coordinates that come back from that conversion, since a real round trip (build a structure → export to `.cdxml` → reconvert to `.mol` → load) confirmed ChemDraw files carry flat 2D drawing-canvas coordinates (`z=0` on every atom, non-physical scale), not real geometry, even though bonds/bond-orders/formula all survive the round trip correctly. So only the atom/bond graph is kept, and it's re-embedded in 3D the same way `from_smiles()` already does -- that shared step (`Chem.AddHs()` → `AllChem.ETKDGv3()` embed → MMFF94 cleanup → lowest-energy conformer) was factored out of `from_smiles()` into `_embed_3d()` so both functions use exactly the same pipeline rather than two copies of it. `from_gaussian_log()` (new): loads a `Structure` straight from an already-finished Gaussian `.log` -- for trusting a geometry someone else already optimized (e.g. a reactant/product for the dedicated TS Search section below) without re-running PM6 on it. Requires `normal_termination` and `stationary_point_found` (raises clearly otherwise -- the whole point of this loader is that the caller doesn't have to independently verify the geometry) and a parseable `Charge = ... Multiplicity = ...` line (new in `parser.py`). `from_file()` dispatches `.log` here the same way it dispatches every other extension, so this is available anywhere a structure is uploaded, not just the new section.
- `gaussbot/input_builder.py` — `Structure`(s) + job type → Gaussian `.com` text. Every route includes `NoSymm` — a symmetry-constrained optimization can get trapped at a symmetric saddle point instead of a lower-symmetry minimum, and it keeps every geometry Gaussian reports in "Input orientation," i.e. the same Cartesian frame we sent in, instead of reorienting to a symmetry-adapted one. `%chk=` naming: every one of `pipeline.py`/`ts_search.py`/`irc.py`'s calls now passes `chk=<the .com path, minus ".com">` explicitly instead of relying on `build_input()`'s old default (`chk_name = chk or job_type`) -- that default meant every job of a given type, across *every study*, wrote to the same generic checkpoint filename (`pm6_opt_freq.chk`, etc.) in the process's own working directory, not even inside `jobs/<name>/`. Harmless while nothing read `.chk` files back; found while checking feasibility for G-corr (below), which is the first thing that actually needs the *right* checkpoint for the *right* molecule -- confirmed the collision was real (five stale, generically-named `.chk` files sitting loose in the project root from earlier sessions) before fixing it. Every attempt (including `_tryN` repairs) now gets its own uniquely-named `.chk` sitting right next to its `.com`/`.log`.
- `gaussbot/intake.py` — a reactant/product from an existing `.com` file, an `.xyz`/`.mol`/`.sdf`/`.pdb`/`.cdx`/`.cdxml` file, or (interactively) a SMILES string; flags an atom count/order mismatch before it becomes a QST2 error. `.cdx`/`.cdxml` dispatch straight through `geometry.from_file()` to the new `from_chemdraw()` -- `intake.load_structure()` itself needed no changes, since it already delegates anything that isn't `.com` to `from_file()`.
- `gaussbot/local_runner.py` — runs a `.com` file directly with `g09 < file.com > file.log` on this machine (no PBS), and checks the log for Gaussian's own "Normal termination" marker. `run_local()` takes an optional `cancel_event` (`threading.Event`) — runs Gaussian via `Popen` and polls it rather than blocking on `subprocess.run`, so a set event terminates the actual `g09` process (SIGTERM, then SIGKILL after a grace period) and raises `GaussianCancelledError` instead of waiting it out; this is what the GUI's "Stop" button ultimately reaches down into.
- `gaussbot/parser.py` — regex-parses a finished `.log`: SCF energy, vibrational frequencies (and which ones are imaginary), the final geometry, charge/multiplicity (new -- `Charge = ... Multiplicity = ...`, confirmed exact log text via `grep` on real logs under `jobs/`, e.g. `" Charge =  0 Multiplicity = 1"`; the last match in the file wins, since a restarted/reread job reprints the block), and — for repairing an imaginary-frequency result — the per-atom Cartesian displacement vector of that mode. Also `last_convergence_status()`: Gaussian's own four-criterion convergence table (Maximum/RMS Force, Maximum/RMS Displacement) from the most recent optimization step, converged or not — what `pipeline.py`'s repair loop now actually looks at before deciding how to retry, instead of reacting the same way regardless of how close an attempt actually got.
- `gaussbot/pipeline.py` — the resilience logic: `repair_and_optimize()` runs an Opt Freq at a given method/basis/resources and repairs it in place if it isn't a clean minimum. Non-convergence handling is now convergence-criteria-driven rather than one-size-fits-all: it reads `parser.last_convergence_status()` at the point an attempt gave up and, if at least 3 of the 4 criteria are already met (`_repair_strategy()`), continues from the last geometry with more cycles and the same small jitter as before (`"continue"`); if fewer are met, the small nudge clearly isn't enough on its own, so it uses a bigger one instead of repeating something that already didn't work (`"escalate"` — jitter scaled `CYCLE_BUMP_JITTER * 2**n`, capped at `MAX_ESCALATE_JITTER`). An imaginary frequency still gets nudged along that mode and reoptimized. Both repair budgets are now one user-facing `max_repair_attempts` (default 4, was two separate fixed numbers) — a geometry genuinely close to converging often just needs more attempts, not a different starting point. `preopt_with_escalation()` wraps `repair_and_optimize()` in a fallback ladder (PM6 on the given geometry → PM6 on a PubChem-refetched geometry → HF/STO-3G) for getting *some* trustworthy stationary point to hand off to the real study. Both functions take the same optional `cancel_event` `local_runner.run_local()` does, passed straight through.
- `gaussbot/level_select.py` — the GaussView-style picks for the real study: `prompt_for_study_type()` (single geometry optimization vs. a full reaction-mechanism study); `prompt_for_method_basis(elements)` (arrow-key list of common method/basis combos, type your own, *or* — if the elements involved make it relevant — "Mixed basis set (GenECP)", which walks into `prompt_for_genecp_assignment()`: checkbox which elements need an effective core potential, a name for each (elements sharing a name are grouped into one block, matching Gaussian's basis/pseudopotential format), and a conventional basis for everything else. Returns a `LevelChoice` (method + either `basis` or `basis_groups`/`ecp_groups`); `prompt_for_resources()` (`%nprocshared`/`%mem` default to a conservative 2/2, with the option to ask for more); `prompt_for_irc()` (whether to run the IRC at all, and how many points per direction, default 25); `prompt_for_ts_distortion_check()`/`prompt_for_irc_endpoint_reopt()` (whether to run `verification.py`'s two extra reactant/product cross-checks — both opt-in, off by default, since each adds two more Opt Freq jobs at the final level).
- `gaussbot/ts_search.py` — the TS stage: `run_ts_search()` tries a TS guess (given, or generated by `geometry.interpolate_structures()` — a crude linear-synchronous-transit interpolation between reactant and product), runs `Opt=(TS,CalcFC,NoEigenTest) Freq` (the `ts_highlevel` job type, GenECP-aware), and requires *exactly one* imaginary frequency whose displacement vector actually overlaps with the reactant→product atomic motion — not just any saddle point, the *right* one — retrying with a different interpolation fraction if it doesn't match. `geometry.align_to()` (Kabsch rigid-body fit) underpins both the interpolation and the mode-matching check, since two independently-optimized Gaussian jobs share no common frame otherwise. The successful `TSOutcome` also carries `imaginary_mode` (the same displacement vector, kept around for `verification.py`'s TS-mode-distortion check rather than reparsing the log a second time). `search_ts_staged()` (new, for the dedicated TS Search section below -- deliberately kept separate from `run_ts_search()` rather than refactored into it, so the existing reaction-mechanism study's TS-search behavior is completely unchanged): a more resilient two-stage search -- a cheap PM6 stage first (unless skipped), whose converged result becomes the fixed starting guess for the expensive final-level stage -- and, at each stage, an opt-in (`use_calcall_fallback`, off by default) `CalcAll` (new job route `ts_calcall`, full-Hessian-every-step, far more likely to locate a genuine saddle point than `CalcFC`'s Hessian-once, much more expensive) brute-force retry if the usual `CalcFC` route can't land a validated TS -- off by default since it's a real resource cost a user shouldn't pay without asking for it; with it off, exhausting the CalcFC attempts just fails that stage, logged as such rather than silently trying something expensive. Accepts `ts_guess`, or `reactant`+`product`, or both (at least one combination required) -- the reactant→product mode-overlap validation only runs when both `reactant`/`product` are given; with a TS-guess-only search, "exactly one imaginary frequency" (a genuine first-order saddle) is the only check available, since there's nothing to compare the mode against. Shares `run_ts_search()`'s small validated helpers (`_min_pairwise_distance`, `_cosine_similarity`, `_reaction_vector`) rather than duplicating them, via a new `_attempt_one_ts_job()`/`_search_ts_at_level()` pair of internal helpers.
- `gaussbot/irc.py` — `run_irc()` runs the IRC from the validated TS (`MaxPoints` per direction — from `prompt_for_irc()`, default 25 — same method/basis/resources) and checks it connects back to the optimized reactant and product by RMSD (`geometry.rmsd()`) in either direction pairing, since Gaussian's forward/reverse labeling is tied to the imaginary mode's arbitrary sign, not a physical direction. A path that hasn't fully relaxed within `MaxPoints` (increase it) looks the same in this check as a genuinely wrong TS (different problem) — either way it's flagged for a manual look rather than silently accepted. Opt-in (`cli.py`/`webapp.py` ask first) rather than always run. `reactant_trusted`/`product_trusted` (default `True`) support the endpoint-recovery path below: when one side has no real reference to check against (its own optimization failed), the usual dual-threshold gate can't apply to it, so `run_irc()` just reports both endpoints as found and leaves classifying *which* is which to `verification.py`'s elimination logic.
- `gaussbot/verification.py` — two independent, opt-in cross-checks that regenerate the reactant/product straight from the verified TS, rather than trusting the original guess-based optimization alone: `run_ts_distortion_check()` displaces the TS's converged geometry along its own imaginary mode (`ts_search.TSOutcome.imaginary_mode`) in the `+`/`-` directions (`x = x0 ± factor·mode`, default `factor=1.0` — a full, undamped step, since Gaussian's per-mode displacement table is already a ready-to-add Cartesian vector) and reoptimizes each side (through `pipeline.repair_and_optimize()`, so it gets the same non-convergence/imaginary-frequency repair logic the final-level reopt already has) to a nearby ground-state minimum; `run_irc_endpoint_reopt()` does the same to the IRC's forward/reverse endpoints, since those are just the last point of a fixed-step path, not necessarily a fully converged minimum. Either way, each reoptimized candidate is classified reactant-side vs. product-side by RMSD (whichever pairing minimizes total mismatch — same either-pairing logic `irc.py` already uses, since which physical direction is "+"/"forward" is arbitrary); if neither pairing is clearly better, the candidate is reported as ambiguous and left out rather than guessed. Both functions also take `reactant_trusted`/`product_trusted` (default `True`), threaded into `_assign_sides()`: when one side isn't trusted (its own optimization failed -- see endpoint recovery below), there's nothing to be ambiguous *between*, so classification switches to elimination against the side that *is* trusted instead of the paired-cost comparison. `select_best_candidate()` then picks whichever candidate (guess-based, TS-distortion, IRC-reopt) is actually lowest in energy (ZPE-corrected, same number `energetics.py` already uses) for reactant and product independently, and that's what the final barrier/reaction-energy report is built from -- a missing/failed `"guess"` candidate is simply skipped, which is what makes it double as the endpoint-recovery mechanism, not just a cross-check. Nothing about this is silent: `format_candidate_comparison()` logs every candidate's energy and which one won, in both the CLI output and the GUI's log panel. `recover_endpoints_from_ts()` (new, for the dedicated TS Search section below): given a found TS, regenerate the reactant/product straight from it. If both `reactant_ref`/`product_ref` are given, this is exactly `run_ts_distortion_check()` (their RMSD to each side is what classifies which is which). If neither is given -- a TS-guess-only search has no reference at all -- the two distorted-and-reoptimized sides come back genuinely *unclassified* (`UnclassifiedRecovery`, a new small dataclass: `side_a`/`side_b` instead of `reactant_candidate`/`product_candidate`) rather than guessed, since there's no way to tell which is which from the TS alone -- that's left to the user to determine chemically. Written to `ts_recover_a.com`/`ts_recover_b.com` in that case, vs. the existing `ts_distort_plus.com`/`ts_distort_minus.com` when there's a real reference to classify against.
- **Endpoint recovery** (in `cli.py`/`webapp.py`'s mechanism-study orchestration, opt-in, off by default via `recover_missing_endpoint`): if the reactant or product can't reach a clean minimum after every fallback -- at PM6 pre-opt *or* the final-level reopt -- the study no longer just aborts there. With recovery on, that side is marked missing (its best, non-converged attempt is kept as a rough TS-search reference -- `pipeline.py`'s `OptOutcome.structure` already holds this even on failure, it just wasn't used downstream before), the other side proceeds normally, and the TS search runs anyway. If the TS is found, `verification.py`'s TS-mode-distortion and/or IRC-endpoint-reopt cross-checks (if enabled -- recovery needs at least one to actually regenerate anything) independently reoptimize *both* sides from the TS, and `select_best_candidate()` picks up whichever one recovers the missing side, since a missing `"guess"` candidate is simply skipped rather than blocking the comparison. If neither recovery method was enabled or succeeded, the study fails with one clear message summarizing what *was* found (the TS, the other side's energy) instead of a bare error. Symmetric -- works for either side failing, not just the reactant, since the rest of the pipeline already treats them identically. Results carry `reactant_recovered`/`product_recovered` booleans so the GUI can say "recovered — original optimization failed" instead of the (different) "more stable than the original optimization" wording used when the guess-based optimization actually succeeded but a cross-check just did better.
- `gaussbot/thermochem.py` — optional, opt-in: G-corr (the thermal correction to Gibbs free energy) for a single structure at a user-chosen temperature, independently per molecule (no reference subtraction, unlike the Energy column). Gaussian's `freqchk` utility computes this from a `.chk` file, but it's interactive, not flag-driven -- there's no command-line way to hand it a temperature. `run_freqchk()` drives it with a verified prompt sequence (checkpoint file → "write Hyperchem files?" (no) → temperature → pressure (1 atm) → frequency scale factor (1, unscaled -- consistent with the rest of gaussbot never applying one) → "use principal isotope masses?" (yes -- skipping this prompt makes `freqchk` loop asking for every atom's isotope mass one at a time, discovered by actually triggering the runaway loop and killing it)) and regexes `Thermal correction to Gibbs Free Energy=` (Hartree) out of the output. `to_kelvin()` converts the GUI's K/°C/°F temperature field. Best-effort throughout: returns `None` rather than raising if `freqchk` is missing, times out, or the checkpoint has no frequency data -- callers skip the value rather than fail the study over it.
- `gaussbot/xyz_export.py` — optional, opt-in: `convert_log_to_xyz()` runs `obabel -ig09 <log> -oxyz -O <xyz>` (Open Babel) to hand back a plain `.xyz` of a finished log's final geometry, for dropping straight into a paper's SI. `convert_log_to_cdxml()` (new) does the same thing to `.cdxml` (`-ocdxml`) -- a ChemDraw-readable export of the final geometry, for dropping straight into a figure. Same best-effort shape both ways -- `None` on failure, never blocks the study.
- `gaussbot/energetics.py` — once the TS is verified: `prompt_for_energy_unit()` (kcal/mol, kJ/mol, eV, or Hartree), then `build_energy_report()` reports the barrier and reaction energy relative to the reactant (= 0), using the ZPE-corrected energy from each Opt Freq log (`LogResult.energy` in `parser.py` — "Sum of electronic and zero-point Energies," the standard 0 K reference for a barrier) rather than the bare electronic SCF energy.
- `gaussbot/cli.py` — the actual entry point. Asks what you want to do first (unless `-r`/`-p` already say): a single geometry optimization (`_run_single_study`), or a full reaction-mechanism study (`_run_mechanism_study` — pre-opt, final opt, TS search, optional TS-distortion cross-check, optional IRC, optional IRC-endpoint reopt, energetics). `gaussbot -r reactant.com -p product.xyz [-n job_name]` runs the mechanism study non-interactively; `gaussbot -r structure.com` (no `-p`) runs a single optimization; bare `gaussbot` asks interactively and, once done, offers to run another study before exiting. `--reactant-id`/`--product-id` give a name/SMILES to use as the PubChem fallback query if PM6 struggles; `--ts-guess` supplies a TS guess file instead of the auto-generated one. Also prompts for `max_repair_attempts` (`level_select.prompt_for_max_repair_attempts()`, default 4) alongside resources, for endpoint recovery (`prompt_for_endpoint_recovery()`, off by default) right before the PM6 pre-opt stage, whether to skip PM6 pre-optimization entirely and start directly at the final level (`prompt_for_skip_pm6()`, PM6 stays the default), and whether to export `.xyz`/`.cdxml` files / compute G-corr (`prompt_for_xyz_export()`/`prompt_for_cdxml_export()`/`prompt_for_gcorr()`, all off by default) at the very end. Filenames inside `jobs/<name>/` no longer repeat the job name as a prefix (`reactant_pm6opt.com`, not `my_reaction_reactant_pm6opt.com`) -- redundant, since the folder itself already says which job it is; `_run_preopt()`/`_run_final_opt()`'s `tag` argument is now a fixed role name (`"structure"`, `"reactant"`, `"product"`) instead of a job-name-prefixed one, and `verification.py`'s two cross-check functions dropped their `tag` parameter entirely in favor of fixed filenames (`ts_distort_plus.com`, `irc_forward_reopt.com`, etc.) for the same reason. The existing `%chk=` fix (checkpoint name derived from the `.com` path) means every `.chk` picks up the shorter name automatically, no separate change needed. When either extra cross-check produces a usable candidate, a comparison line (`verification.format_candidate_comparison()`) is printed for reactant and product before the final energy report, naming which candidate actually won; the final summary also prints the actual `.log` (and, if requested, `.xyz`/G-corr) file/value for reactant, TS, and product. No Stop/Restart here -- that's GUI-only (see below); the CLI is already interactively Ctrl+C-able, it just doesn't clean up the Gaussian child process on its own yet.
- `gaussbot/job_runner.py` — the PBS-side alternative to `local_runner.py`: `build_pbs_script()`/`write_pbs_script()` reproduce the site's own submission script line for line (Intel setvars, `g09.profile`, a per-job `GAUSS_SCRDIR` under `/scratch/$PBS_JOBID`, cleaned up after) — verified against the actual template text (the user supplied it directly: a `g09.sh` that writes exactly this script and `qsub`s it), not just eyeballed; `DEFAULT_G09ROOT`/`DEFAULT_INTEL_SETVARS`/`DEFAULT_PPN`/`DEFAULT_WALLTIME` match that one site's setup exactly, and will need editing for a different cluster (the GUI/CLI both warn about this explicitly wherever PBS mode is turned on). `submit_job()` wraps `qsub` and returns the job id; `job_status()`/`wait_for_job()` wrap `qstat -f` and poll until a job leaves the queue; `cancel_pbs_job()` wraps `qdel`. **Now wired in** via `run_pbs()`: the queue-based counterpart to `local_runner.run_local()` -- same contract (blocks until done, returns the log path, raises `local_runner.GaussianRunError`/`GaussianCancelledError` on failure/cancellation, `cancel_event` triggers a `qdel` instead of a `SIGTERM`) so it's a genuine drop-in, not a special case.

  **User-editable PBS script, not just a "go edit the .py constants" warning**: instead of only exposing the `DEFAULT_*` module constants, `default_pbs_template()` returns `build_pbs_script()`'s usual text with the job name/`.com` filename left as a placeholder (`__JOB__`, appearing in both `#PBS -N` and the `g09` line) -- the GUI's editable text box (see below) or the CLI's `prompt_for_pbs_template()` (opens it in `$EDITOR`, interactive sessions only -- scripted `-r`/`-p` runs never block on an editor, they just fall back to the module defaults) let the user rewrite *any* part of it (paths, `nodes`/`ppn`/`walltime`, extra module loads, anything), not just the four constants. `render_pbs_script(template, job_name)` substitutes every `__JOB__` occurrence with one job's actual name. Rather than threading the (possibly large) template string through every function that might run a Gaussian job, it's written **once**, per job/study, to `<out_dir>/pbs_template.sh` (`webapp.py`'s `_maybe_write_pbs_template()` / `cli.py`'s twin) -- `run_pbs()` checks for that file next to each `.com` it's about to submit and uses it (rendered per-job) instead of `build_pbs_script()`'s defaults if it's there. GUI: the "Enable PBS submission" button (see below) reveals a `<textarea>` pre-filled from `GET /api/pbs_template`, fetched once per session so an in-progress edit survives toggling the button off and back on.

  No real PBS system is reachable from here to submit a job to, so `run_pbs()` -- like the rest of this module -- is tested against a fake `qsub`/`qstat`/`qdel` (`tests/test_job_runner.py`: success, submission failure, missing-log-after-queue-exit, cancellation-triggers-qdel, template-substitution, and a full run using a `pbs_template.sh` written next to the `.com`) rather than validated live; the actual queue behavior of a real cluster is the one thing in this project that's never had a real end-to-end run. The template-editor path was however validated for real up to that boundary: a real GUI submission with PBS mode on and an edited template (`ppn`/`walltime` changed in the text box) correctly wrote the edited text to `jobs/<name>/pbs_template.sh` with both `__JOB__` placeholders intact, and failed with a clear `qsub: No such file or directory` (this machine genuinely has no PBS system) rather than silently doing the wrong thing.
- **`executor` switch (local vs. PBS)** — `gaussbot/executor.py`'s `run_com(com_path, executor="local"|"pbs", cancel_event=None)` is the one place that decides between `local_runner.run_local()` and `job_runner.run_pbs()`; every function that actually runs a Gaussian job (`pipeline.repair_and_optimize()`/`preopt_with_escalation()`, `ts_search.run_ts_search()`/`search_ts_staged()`, `irc.run_irc()`, `verification.run_ts_distortion_check()`/`run_irc_endpoint_reopt()`/`recover_endpoints_from_ts()`) now takes an `executor` parameter (default `"local"`, so nothing changes unless you opt in) and threads it straight through to `run_com()` instead of calling `run_local()` directly. `webapp.py`'s `JobRequest.executor` and `cli.py`'s `--executor {local,pbs}` flag (or, for a bare interactive `gaussbot`, a session-wide prompt asked once, not per-study) set it at the top and it flows all the way down. GUI: a small blue "Enable PBS submission" button next to "1. What do you want to do?" (turns green when on, `qsub`/`qstat`/`qdel` need to be on PATH) toggles it for every study submitted afterward, until clicked again -- and reveals the editable PBS-script text box described above, not just a static warning. Validated: the default `"local"` path re-run for real end to end after this refactor touched every call site (confirmed byte-identical behavior, same SCF energy as before); the `"pbs"` path itself can't be validated against a real queue from here, same limitation as `job_runner.py` above -- covered by its fake-qsub tests and code review instead.
- `gaussbot/webapp.py` + `gaussbot/web/` — a local browser GUI, `gaussbot-gui`. FastAPI backend that calls the exact same `pipeline.py`/`ts_search.py`/`irc.py`/`energetics.py` functions `cli.py` does — nothing about the actual chemistry is reimplemented, this only replaces the questionary/`input()` layer with a web form. Every submitted study runs on its own background thread with its own `JobState`, so several studies can run at once rather than one at a time; the frontend polls `GET /api/jobs/{id}` every 1.5s per study for a live-updating stage indicator and log panel, rendered as its own card in a "Studies" list. The form itself only ever starts one job per submit, but a "+ Start another study" button (present both while a study is running and after it finishes) reopens the form without disturbing the studies already in flight — click it, fill in a different structure/job name, submit again, and it runs alongside the others. Covers study type, structure input as either a file upload (`POST /api/analyze`) or a SMILES string (`POST /api/analyze_smiles`, `geometry.from_smiles` -- the GUI's counterpart to the CLI's SMILES intake option) -- both return the element list, used to drive the GenECP checkbox picker -- resources, a `max_repair_attempts` field, method/basis (including the same GenECP flow as the CLI), IRC opt-in, an energy-unit picker that's now a checkbox group rather than a single select (pick kcal/mol *and* eV, say, and the results panel reports the barrier and reaction energy in each unit chosen), two opt-in checkboxes for `verification.py`'s TS-mode-distortion and IRC-endpoint-reopt cross-checks (the second only enabled once "Run the IRC" is checked), a third opt-in checkbox for endpoint recovery (only useful alongside at least one of the other two, noted in its hint text), a checkbox to skip PM6 pre-optimization and start directly at the final level (PM6 stays the default), and opt-in checkboxes to export `.xyz` files (`xyz_export.convert_log_to_xyz()`) and/or `.cdxml` (ChemDraw) files (`xyz_export.convert_log_to_cdxml()`) and/or compute G-corr (`thermochem.run_freqchk()`, revealing a temperature field + K/°C/°F unit select) once the study finishes. Structure uploads (reactant/product/TS-guess) accept `.cdx`/`.cdxml` too, same as the CLI -- no extra backend work needed, `/api/analyze` already dispatches by extension straight to `geometry.from_file()`. Results render as two tables per study (`main`'s width bumped from 720px to 960px so there's room) instead of the old flat row list: an Energy table (Sr No / Compound / Energy per requested unit / optional G-corr per unit, one row each for Reactant, TS, Product) and a Files table (Sr No / Compound / Logs / optional XYZ / optional CDXML) sitting beside the reaction-coordinate plot in a wrapping flex row -- side by side on a wide viewport, stacked on a narrow one. The reaction-coordinate plot has a small print-icon button in its corner (`app.js`'s `wirePrintButtons()`) that grabs the live `<svg>`, opens a small new window with its own plain black-on-white stylesheet (not the app's dark theme), and calls `.print()` on it -- opened synchronously from the click so it isn't blocked as a popup. Filenames inside `jobs/<name>/` dropped the job-name prefix here too, same reasoning and same fixed role names (`"structure"`/`"reactant"`/`"product"`) as the CLI change above. A "Tips" button by the title opens a static panel of practical dos-and-don'ts (not math) -- reasonable starting geometries, GenECP element coverage, IRC point counts, charge/multiplicity sanity checks, TS-guess quality, resource sizing. Each running study card has a "Stop" button (`POST /api/jobs/{id}/cancel`, sets `JobState.cancel_event` -- reaches all the way down to `local_runner.run_local()` terminating the actual `g09` process) and, once cancelled, a "Restart" button (`POST /api/jobs/{id}/restart`) that replays the same request but skips any stage (`JobState.checkpoints`) that already finished cleanly before the cancellation -- pre-opt, final opt, and TS search are checkpointed; the opt-in extras (IRC, TS-distortion, IRC-endpoint reopt) are cheap enough by comparison to just re-run rather than also being resumable. Every study card also has a small trash-icon button (`DELETE /api/jobs/{id}`) that permanently removes its entry from the Studies list -- refused while the study is still running (stop it first), and never touches the actual `.com`/`.log` files under `jobs/<name>/`, only the in-memory `JobState` the GUI itself was tracking. Validated end to end in a real browser (Chromium via the dev-preview tooling, not just the API): uploaded real `.xyz` files through the actual file input, loaded structures from SMILES through the actual "or SMILES" field, drove the GenECP element checkboxes, toggled study type, ran three real mechanism/single studies concurrently through "+ Start another study" (confirmed independent live progress and independent results per card, including a kcal/mol+eV run whose reported barrier in each unit matched the expected conversion factor), ran full mechanism studies through the real UI from submit to a rendered results panel (HF/STO-3G NH3 inversion, barrier 9.53 kcal/mol, IRC "connected" badge), and -- for the Stop/Restart feature -- started a real mechanism study, confirmed via `ps` that a real `g09` process was running, clicked Stop, confirmed via `ps` again that the process was actually gone (not just that the UI stopped polling), clicked Restart, and confirmed the resumed run's log showed "Reusing already-completed ..." for the three stages that had genuinely finished before cancellation and for-real reran (with a matching, correct SCF energy) the one stage that hadn't, then completed the rest of the study normally.
- **GUI polish** — three small fixes/tweaks made alongside the TS Search work: (1) the G-corr temperature/unit fields showed on page load even before "Compute G-corr" was checked -- root cause was `.field-row { display: flex }` (author CSS) overriding the browser's own `[hidden]` rule (author stylesheets always win over the UA stylesheet regardless of selector specificity), fixed with an explicit `.field-row[hidden] { display: none; }`; (2) the reaction-plot's print button was printing its own printer *icon* instead of the plot -- `.reaction-plot` has two `<svg>` elements (the button's icon, then the actual chart) and `querySelector("svg")` grabbed whichever came first in the DOM (the icon); fixed by querying `.reaction-plot-svg` specifically; (3) "Reactant file"/"Product file" field labels are now bold (`#reactant-label`, `#product-label`) so they stand out from the surrounding format-hint text, which is deliberately lighter-weight; (4) the plot title ("Reaction coordinate (...)") and the TS point's energy label could overlap and render as illegible stacked text -- the TS is usually the highest point on the curve (`yMax`), so its marker sat right at the plot's top padding (`padTop`) with its value label offset just 10px above that, landing in the same few pixels as the title. Fixed by increasing `padTop` from 22 to 38, giving enough clearance between the title and the highest point's label regardless of which point that turns out to be. All four verified in the real browser (computed-style checks before/after, and intercepting `window.open` to confirm the print window's HTML now contains the actual chart markup with non-overlapping text, not the icon's path data).
- **Dedicated TS Search** — a third study type, alongside single-optimization and full mechanism study, for finding/verifying just the transition state without committing to a full study (no forced IRC/energetics). `level_select.prompt_for_study_type()` gained the third choice; `intake.prompt_for_ts_search_inputs()` (CLI) asks for a TS guess, a reactant+product pair, or both -- at least one of {TS guess} or {reactant AND product} required -- and warns (non-blockingly) if a reactant/product isn't already a `.log`, since this section, unlike the full mechanism study, does **not** PM6-preopt or final-level-reoptimize them; they're used exactly as given, straight into `ts_search.search_ts_staged()`. `cli.py`'s `_run_ts_search_study()` and `webapp.py`'s `_run_ts_search_job()` (`JobRequest.study_type = "ts_search"`, `reactant_token`/`product_token` now `Optional`, new `ts_recover_endpoints`/`ts_use_calcall` fields) mirror each other, both calling `search_ts_staged()` then, if requested, `verification.recover_endpoints_from_ts()`. GUI: a third radio button in section 1 reuses the *same* reactant/product/TS-guess upload fields and wiring the mechanism study already has (marked optional, with a small warning line -- `#reactant-log-warning`/`#product-log-warning` in `index.html`, toggled by `app.js`'s `updateLogWarning()` -- shown whenever the uploaded file isn't `.log`), a "skip PM6 stage" checkbox (reuses the existing `skip-pm6-preopt`), a CalcAll-fallback opt-in checkbox (`#ts-calcall-fallback`, unchecked by default), and a "recover reactant/product from the TS" checkbox, and section 6 ("Energy report") was split into an always-visible part (the `.xyz`/`.cdxml` export checkboxes -- previously only reachable for a mechanism study, since the whole section was hidden for "single" too; now visible there as well, a side-effect fix) and a `#energy-units-field`/`#gcorr-block` part hidden for TS Search (no barrier to report; G-corr isn't wired into `_run_ts_search_job()`). `renderStudyResults()` gained a `type === "ts_search"` branch: found/not-found, imaginary frequency, mode-overlap (or "n/a" when no reactant/product was given to check against), the TS `.log` path, and — if recovery was requested — either classified "Recovered reactant"/"Recovered product" rows or neutrally-labeled "Recovered side A"/"Recovered side B" rows, matching whichever of `DistortionOutcome`/`UnclassifiedRecovery` came back. Validated end to end for real: `search_ts_staged()` run directly against a real `g09` in all three input modes (TS-guess only, reactant+product only, all three together) on the vinyl-alcohol⇌acetaldehyde toy reaction -- confirmed the PM6 stage's converged TS correctly seeds the final-level guess, confirmed the overlap check is skipped (not just defaulted) when no reactant/product is given, confirmed `ts_pm6.com`/`ts_final.com` naming; `recover_endpoints_from_ts()` run directly for both the classified case (RMSD-classified reactant/product, correct energies) and the guess-only case (`UnclassifiedRecovery`, correctly labeled, correct honesty caveat in the log); a full real GUI submission (TS-guess-only mode, recovery enabled) confirmed the results panel renders every field correctly including the neutral Side A/B labeling.
- `gaussbot/llm_assist.py` — optional, off by default: reads a free-text description the user types next to a loaded reactant/product/TS-guess structure (a reaction name, a compound's common/IUPAC name, informal notes) and asks Claude to suggest a better starting point for it — a cleaner PubChem search query, a corrected SMILES if the given one looks inconsistent with the description, or charge/multiplicity if the description implies something non-default (a radical, an ion, a triplet). It never generates 3D coordinates itself — LLMs aren't reliable at that — geometry generation stays entirely with the existing deterministic tools (RDKit, PubChem); the LLM only ever picks a better lookup key or flags a likely mismatch. Wired into the GUI as `GET /api/refine/available` (a cheap, best-effort check for whether a credential is configured, used only to set an informational tooltip) and `POST /api/refine` (`webapp.py`): any PubChem hit or corrected SMILES the LLM suggests is re-embedded/fetched for real and formula-checked (`Structure.formula()`, Hill notation) against the original structure before it's ever offered back — a mismatched formula is silently dropped with an explanatory note instead of being offered. Nothing is applied automatically; the user sees the suggestion and clicks "Use this structure" (or doesn't). Needs `ANTHROPIC_API_KEY` (or another credential source the SDK resolves automatically) — without one, `is_available()` returns `False` and `/api/refine` returns a `503` with the real error message rather than pretending to work. Validated end to end except the live model call itself: `/api/refine`'s request/response wiring, the formula-match safety check (a real PubChem lookup for "vinyl alcohol" accepted, one for "benzene" against the same starting structure correctly rejected), and the graceful-failure path (no credential in this dev sandbox) were all exercised for real, with only the `anthropic.messages.parse()` call itself mocked or, in the browser, genuinely failing on the missing credential and surfacing that error in the UI as designed.

Known limitations, all found by actually running things rather than guessing:
- The auto-generated TS guess (Cartesian interpolation after Kabsch alignment) is degenerate for a linear-to-linear reaction — the rotation about the shared molecular axis is mathematically undetermined, so the alignment can pick a bad one and land two atoms on top of each other. `ts_search.py` catches this (checks the guess's minimum pairwise distance before ever handing it to Gaussian) and reports it clearly rather than letting g09 segfault on a "Small interatomic distances" error, but the real fix — a hand-supplied `--ts-guess` — is still on the user for that class of reaction. Found via the actual HCN⇌HNC toy reaction, which is exactly this pathological case (it's also why that reaction's own QST2 job fails, unrelated to this).
- `IRC=(CalcFC,Both)` — the syntax you'll see in most Gaussian references — is a syntax error on this installed G09 revision (RevC.01); `IRC=(CalcFC)` alone already runs both directions by default. `input_builder.py`'s route reflects this. Worth rechecking if this ever runs against a different Gaussian version.
- The site's original PBS script builds `. $g09root/g09/bsd/g09.profile` with a *double*-quoted `echo`, so `$g09root` there gets expanded immediately by the shell generating the script (the one running `qsub`), not the one that later executes it on the compute node — it only produces the right path if `g09root` happens to already be set in that outer shell's environment. `job_runner.build_pbs_script()` takes `g09root` as one explicit parameter used for both the export line and the source line, so the two are guaranteed to match regardless of environment. Everything else about the script is reproduced exactly as given.
- GenECP element assignment is interactive-only right now — no CLI flag equivalent for scripted/non-interactive mixed-basis runs (`-r`/`-p` skip *structure* intake, not the method/basis/resources/IRC prompts, which is also why those still show up even when `-r`/`-p` are given).
- Jittering a strictly *linear* molecule during the cycle-bump repair can trip Gaussian's own "Error in internal coordinate system" (a Berny-optimizer quirk for near-linear systems) — found on the same HCN test case that already breaks QST2 and the TS-guess alignment. Validated the repair fix on a non-linear molecule instead; a linear one hitting this would currently just exhaust its repair budget and report failure rather than recovering.
- The GUI's progress panel updates *stage by stage* (each PM6 pre-opt, each final opt, TS search, IRC), not attempt by attempt within a stage — `webapp.py` only reads each `outcome.log` once a stage function returns, deliberately, so this GUI pass didn't need to touch `pipeline.py`/`ts_search.py`/`irc.py`'s internals. A repair loop chewing through several retries will just look quiet until it finishes rather than streaming each attempt live. Threading a progress callback through those functions would fix this if it's worth the finer granularity later.
- The GUI's uploaded files land in `/tmp/gaussbot_uploads` (or `$GAUSSBOT_UPLOAD_DIR`) under random tokens and are never cleaned up automatically — fine for local single-user use, worth adding cleanup if this ever runs somewhere more shared.
- `llm_assist.py`'s actual Claude API call (`anthropic.messages.parse()`) has never run against a live credential — this dev environment has no `ANTHROPIC_API_KEY` configured. Everything around that call has been validated for real (the `/api/refine` request/response wiring, the formula-match safety check against real PubChem lookups, the graceful-failure path when there's no credential), but the model's actual suggestion quality — whether its PubChem queries and SMILES corrections are actually good — is unverified. Worth a real smoke test with a credential before relying on it for anything beyond a quick sanity check.
- `verification.py`'s reactant/product-side classification can come back genuinely ambiguous, and does in practice, not just in theory: run for real on the vinyl-alcohol⇌acetaldehyde case with a deliberately short IRC (15 points instead of the 25 default, chosen to keep the test fast), the IRC forward and reverse endpoints reoptimized to the *same* minimum (identical energy, identical RMSDs to both reactant and product) rather than one to each — because the IRC itself hadn't fully separated the two directions within that few points. `run_irc_endpoint_reopt()` correctly detected this (`CLASSIFY_AMBIGUOUS_MARGIN`) and reported neither candidate rather than picking one arbitrarily, but it's a real illustration that these cross-checks are only as good as what feeds them — a too-short IRC won't get a wrong answer out of this, but it can get *no* answer where a longer one would.
- The convergence-criteria-driven repair strategy (`pipeline._repair_strategy()`) is a heuristic, explicitly -- "continue if ≥3/4 criteria met, otherwise use a bigger kick" is a reasonable read of the real stuck case it was built against (`jobs/diels-2/..._try3.log`: Force converged, Displacement oscillating for 28 steps without trending down), but it hasn't been run across a large enough sample of genuinely difficult geometries to know its real-world success rate. The escalating-jitter magnitude (`CYCLE_BUMP_JITTER * 2**n`, capped at `MAX_ESCALATE_JITTER`) is likewise a reasonable-looking default, not a tuned one.
- Restart's stage-skipping is deliberately *stage-level*, not a true Gaussian checkpoint (`%chk`/`Opt=Restart`) warm restart -- a cancelled job resumes from the last fully completed stage (pre-opt, final opt, TS search), not from wherever mid-optimization-cycle it actually was killed. A cancellation that lands mid-way through, say, the final-level reactant optimization means that whole stage reruns from its pre-opt output on restart, not from its own partial progress. This was an explicit scope decision (a true warm restart is fragile across method/basis changes and a much bigger lift) rather than an oversight, but it's worth knowing before relying on it to save time on a very expensive stage killed near its end.
- Endpoint recovery's *failure-trigger* path (a side genuinely can't optimize, the study keeps going without it) was validated for real with a deliberately unoptimizable reactant geometry, and the elimination-classification math it depends on (`irc.py`/`verification.py`'s `*_trusted` flags) was validated in isolation against toy structures -- but a full live run all the way through "reactant fails → TS is found anyway using the fallback geometry → TS-mode distortion successfully recovers a proper reactant → the study finishes" was not observed end to end in this environment: forcing a real, moderate (not catastrophic) optimization failure that still leaves the geometry usable enough for a genuine TS search to succeed turned out to be hard to hit reliably by hand in the time available -- real molecules here either converged fine or failed so badly (clashing atoms) that no real TS existed to find between them and the product. The individual pieces are each proven; the full chain together is validated by code review, not a live run.
- The results panel's wrapping flex layout (tables beside the reaction plot on a wide viewport, stacked on a narrow one) was confirmed at desktop width via real `getBoundingClientRect()` measurements (same top offset, non-overlapping) -- the narrow/stacked side wasn't directly screenshotted, since the browser-resize tool wasn't reflecting into the page's own `window.innerWidth` in this environment at the time. The CSS itself (`display:flex; flex-wrap:wrap` with `flex-basis` values around 320-380px) is standard, well-understood behavior that reliably wraps once the container drops below their combined width, so this is a low-risk gap, but it's a real one -- worth an actual narrow-viewport screenshot next time the GUI is open somewhere that cooperates.
- ChemDraw export is `.cdxml` only, not the older binary `.cdx` -- confirmed via `obabel -L formats` that `.cdx` is read-only in this Open Babel build (`.cdx` upload works fine, since that's reading; there's just no way to write one back out). `.cdxml` is ChemDraw's own XML variant and opens in ChemDraw the same as a `.cdx` would.
- A ChemDraw upload's 3D geometry is only as good as RDKit's usual force-field embedding (ETKDGv3 + MMFF94), same caveat `from_smiles()` already has and already documented above (a force-field guess, not QM quality) -- since ChemDraw's own coordinates are discarded and can't be used as a better starting point even when the sketch was drawn with reasonable bond angles/geometry in mind. Nothing lost relative to typing the same molecule as SMILES instead; just worth knowing a careful ChemDraw sketch doesn't buy a better starting geometry than a SMILES string would.

## Structure

```
gaussbot/
  gaussbot/
    __init__.py
    geometry.py        # SMILES/file/PubChem -> Structure; displace_structure()
    input_builder.py    # Structure(s) + job_type -> .com text
    intake.py             # reactant/product file loading + interactive prompts
    local_runner.py        # run a .com file locally with g09
    job_runner.py           # PBS script generation + qsub/qstat -- not wired into cli.py yet
    parser.py                # .log -> SCF energy, frequencies, geometry, mode vectors
    pipeline.py                # repair_and_optimize / preopt_with_escalation
    level_select.py              # method/basis + resource picks (questionary)
    ts_search.py                   # TS guess generation, search, mode-match check
    irc.py                           # IRC + reactant/product connectivity check
    energetics.py                      # unit picker + relative-energy report
    cli.py                               # python -m gaussbot.cli -- the actual entry point
    webapp.py                             # FastAPI backend for the local GUI (gaussbot-gui)
    llm_assist.py                          # optional Claude-assisted structure refinement (advisory only)
    web/                                   # GUI frontend: index.html, app.js, style.css
  examples/
    build_pm6_and_qst2.py  # working demo of geometry.py + input_builder.py
  jobs/                     # per-job working directory (gitignore this)
  requirements.txt
```

## Install

```
cd gaussbot
python -m venv venv
source venv/bin/activate        # macOS/Linux; venv\Scripts\activate on Windows
pip install -r requirements.txt
pip install -e .
```

That last line registers `gaussbot` as a real terminal command (via
`pyproject.toml`'s entry point), so from then on:

```
gaussbot
```

runs the CLI from anywhere, as long as the venv is active. Next time
you open a terminal, you only need `source venv/bin/activate` again
-- no need to reinstall.

You'll also need Gaussian itself on PATH for `local_runner.py` to
find it -- source your site's Gaussian environment script first
(the same line your PBS scripts already use), e.g.:

```
source /path/to/g09/bsd/g09.profile
```

Two more system tools are needed, but only for the two GUI/CLI features
that are opt-in and off by default -- gaussbot runs fine without either
if you never check those boxes:
- [Open Babel](https://openbabel.org) (`obabel` on PATH) for the
  optional "export to .xyz"/"export to .cdxml" features
  (`xyz_export.py`) and for loading `.cdx`/`.cdxml` ChemDraw files as
  a reactant/product/TS-guess input (`geometry.from_chemdraw()`).
- `freqchk`, which ships with Gaussian itself (same `g09.profile`
  above puts it on PATH), for the optional G-corr feature
  (`thermochem.py`).

Uses RDKit for geometry embedding/parsing and the displacement helper,
`pubchempy` for the PubChem fallback lookup, `questionary` for the
CLI's method/basis/resource picker, `fastapi`/`uvicorn`/
`python-multipart` for the GUI (`gaussbot-gui`), and `anthropic` for
the GUI's optional AI-assisted structure refinement (`llm_assist.py`
-- only imported when that feature is actually used; the rest of
gaussbot has no dependency on it). `parser.py` is hand-rolled regex,
not `cclib` — the handful of fixed-format lines the pipeline actually
branches on haven't changed since G03, so a full log-format library
was more than this needed.

The AI-assisted refinement feature is entirely optional and off by
default in practice -- it just needs an `ANTHROPIC_API_KEY`
environment variable set (get one from
[console.anthropic.com](https://console.anthropic.com)) wherever
`gaussbot-gui` runs:

```
export ANTHROPIC_API_KEY=sk-ant-...
gaussbot-gui
```

Without it, the "Refine with AI" button still appears but returns a
clear error instead of silently doing nothing -- nothing else in
gaussbot is affected either way.

## Usage so far

```python
from gaussbot import from_smiles, from_file, build_input, write_com

reactant = from_smiles("C#N", label="HCN reactant")
product = from_smiles("[C-]#[NH+]", label="HNC product")

# QST2 requires the SAME atoms in the SAME order in both structures --
# it's one molecule's atoms moving, not two different molecules.
assert [a[0] for a in reactant.atoms] == [a[0] for a in product.atoms]

com_text = build_input("qst2", [reactant, product], method="PM6")
write_com(com_text, "jobs/ts_qst2.com")
```

A ChemDraw `.mol`/`.sdf` export works the same way:

```python
struct = from_file("my_molecule.mol", label="reactant")
```

A ChemDraw `.cdx`/`.cdxml` sketch also works, but through a different
path -- ChemDraw's own coordinates are 2D drawing-canvas positions, not
real geometry, so `from_file()` routes these two extensions to
`from_chemdraw()`, which discards the coordinates and re-embeds the
atom/bond graph in 3D the same way a SMILES string would be:

```python
struct = from_file("my_molecule.cdxml", label="reactant")
```

Job types currently supported by `build_input` (see
`input_builder.JOB_ROUTES` for the exact route strings):

| job_type         | structures needed | typical use                          |
|------------------|--------------------|---------------------------------------|
| `pm6_opt`        | 1                  | cheap pre-optimization                |
| `pm6_opt_freq`   | 1                  | pre-opt + frequency check             |
| `qst2`           | 2 (reactant, product) | PM6 TS guess from endpoints       |
| `qst3`           | 3 (reactant, product, guessed TS) | PM6 TS guess with a seed |
| `ts_highlevel`   | 1 (TS guess)       | DFT-level TS refinement               |
| `opt_highlevel`  | 1                  | DFT-level opt of reactant/product     |
| `irc`            | 1 (converged TS)   | confirm the TS connects R and P       |

## Run it

Reaction-mechanism study:
```
gaussbot -r reactant.com -p product.xyz -n my_reaction
```
Single geometry optimization (no `-p`):
```
gaussbot -r structure.com -n my_molecule
```
Bare `gaussbot` asks interactively which of the two you want, accepts
a SMILES string (the only path that does), and once a study finishes,
offers to run another before exiting.

For a mechanism study: loads the reactant/product (`.com`, `.xyz`,
`.mol`, `.sdf`, `.pdb`, `.cdx`, or `.cdxml` — picked by extension),
warns if the atom count/order doesn't match between them, then:

1. Runs the resilient PM6 pre-optimization on each (retrying through
   non-convergence and imaginary-frequency repairs, escalating to
   HF/STO-3G if PM6 can't get there), writing
   `jobs/my_reaction/reactant_pm6opt.com` and
   `product_pm6opt.com` (filenames inside `jobs/<name>/` don't repeat
   the job name -- the folder already says which job it is).
2. Prompts for the method/basis set to actually run the study at
   (arrow-key list, type your own, or a mixed GenECP basis — checkbox
   which elements need an effective core potential, e.g. a transition
   metal, name it, and pick a conventional basis for the rest) and how
   much of this machine to give Gaussian (`%nprocshared`/`%mem`
   default to a conservative 2/2, with the option to ask for more).
3. Reoptimizes both at that level through the same repair loop,
   writing `reactant_final.com` and `product_final.com`.
4. Prompts for a TS guess file (blank to auto-generate one by
   interpolating between the optimized reactant/product), runs the TS
   search, and checks the resulting imaginary frequency actually
   corresponds to the reactant↔product motion — writing
   `ts.com` once it finds one that matches.
5. Asks whether to also cross-check the reactant/product by
   distorting the converged TS along its own imaginary mode (`x = x0
   ± factor·mode`, default `factor=1.0`) and reoptimizing both sides
   (`ts_distort_plus.com`/`ts_distort_minus.com`) — off by default,
   since it adds two more Opt Freq jobs.
6. Asks whether to run the IRC, and if so how many points per
   direction (default 25) — writing `irc.com` and checking it
   connects back to the optimized reactant and product.
7. If the IRC ran: asks whether to also reoptimize its forward/reverse
   endpoints to proper minima (`irc_forward_reopt.com`/
   `irc_reverse_reopt.com`) — a third candidate, same idea as step 5,
   also off by default.
8. If either of steps 5/7 produced a usable candidate, compares its
   energy against the guess-based reactant/product (and against each
   other) and uses whichever is actually lowest for the final report —
   printing which one won and why.
9. Prompts for a unit (kcal/mol, kJ/mol, eV, Hartree) and prints the
   barrier and reaction energy relative to the reactant. (The GUI's
   equivalent step lets you pick more than one unit at once; the CLI
   still asks for a single one.)

A single optimization runs just step 1 (with the one structure instead
of a reactant/product pair) and step 2-3's final-level reopt, then
stops -- no TS/IRC/energetics, since those need a reaction. `gaussbot
--help` for the full flag list.

Note: `local_runner.py` only looks for `g09` on `PATH`. Your
interactive shell already has this via `.bashrc` sourcing
`g09.profile` (see Install below), but a non-interactive shell (cron,
a PBS batch script once that's wired back up) won't, and will fail
with a cryptic "No executable for file l1.exe" unless it sources
`g09.profile` itself first.

## Run the GUI

```
gaussbot-gui
```

Starts a local web server (`http://127.0.0.1:8765`) and opens it in
your browser. Same pipeline as the CLI, same resilience logic
underneath, just a web form instead of terminal prompts: pick study
type, upload your reactant/product (and optional TS guess) files
(`.com`, `.xyz`, `.mol`, `.sdf`, `.pdb`, `.cdx`, `.cdxml` -- shown
right on the label, and every one of them gets auto-converted to a
Gaussian input; a ChemDraw `.cdx`/`.cdxml` sketch is re-embedded in 3D
rather than trusting its 2D drawing coordinates -- see
`geometry.from_chemdraw()` above), set resources, pick the
method/basis (including the GenECP element picker), opt in to the
IRC, opt in to the two extra reactant/product cross-checks (TS-mode
distortion, and -- once IRC is on -- reoptimizing its endpoints; see
`verification.py` above), pick one or more energy units (check as
many as you want -- kcal/mol and eV together, say -- and the results
panel reports the barrier and reaction energy in each), opt in to
exporting the final structure(s) to `.xyz` and/or `.cdxml`
(ChemDraw), hit "Run study" -- then watch a live stage indicator and
log panel update as it works, ending in a results panel that names
which candidate (guess-based, TS-distortion, IRC-reopt) was actually
used if it wasn't the default one, the actual `.log` file each number
came from, and a small reaction-coordinate plot with a print button
in its corner for a quick black-on-white printout.
While a study is running, a "Stop" button on its card terminates the
actual Gaussian process; once stopped, a "Restart" button replays the
same study, skipping straight past whatever stages already finished
cleanly. Click "Tips" by the title for a quick panel of practical
dos-and-don'ts. Ctrl+C in the terminal it's running in to stop the
server.

Below each reactant/product structure field there's an optional
"describe this compound/reaction" text box and a "Refine with AI"
button (needs `ANTHROPIC_API_KEY` -- see Install above). Type
whatever you'd naturally call the thing -- a common name, informal
notes, "the enol tautomer of acetaldehyde" -- and Claude reads it
alongside the structure you already loaded and suggests a better
starting point: a PubChem query, a corrected SMILES, or
charge/multiplicity, if any of those seem warranted. Nothing is
applied automatically -- you review the suggestion and click "Use
this structure" yourself, or ignore it and keep what you had.

Once at least one study has been submitted, a "+ Start another
study" button stays available (next to the studies list, both while
studies are running and after they finish) -- click it to reopen the
form and submit a different study without waiting for the first one,
and both run side by side with their own live progress card.

## Run the demo

```
python examples/build_pm6_and_qst2.py
```

Builds PM6 pre-opt inputs for a toy HCN/HNC isomerization plus the
PM6 QST2 input, and prints the QST2 `.com` file. Not a real study —
just proof geometry.py and input_builder.py work together.

## Testing

```
pytest tests/
```

`tests/test_local_runner.py` and `tests/test_job_runner.py` are the
only real pytest coverage — both test against a fake stand-in
executable (`g09`, or `qsub`/`qstat`) rather than the real thing,
since neither a Gaussian install nor a PBS system was assumed to be
reachable in general. There's no pytest coverage for `geometry.py`,
`input_builder.py`, `intake.py`, `parser.py`, `pipeline.py`,
`level_select.py`, `ts_search.py`, `irc.py`, or `energetics.py` —
everything in those has instead been validated by actually running
real reactions through the real `g09` that *is* installed on this
particular machine:

- HCN reactant/product PM6 opt+freq, then HF/STO-3G reopt.
- A deliberately bad HCN guess, to force the non-convergence →
  MaxCycles-bump repair.
- A symmetric planar NH3 guess, a real PM6 saddle point, to force the
  imaginary-frequency → displacement repair.
- The NH3 umbrella inversion (two pyramidal minima related by
  inversion through the H3 plane), run through the *entire* pipeline
  end to end -- PM6 pre-opt, final opt, TS search, IRC, energetics:
  found the same planar D3h TS independently verified above
  (imaginary freq -962.5 cm⁻¹, matching the -961.6 cm⁻¹ found
  earlier), correctly flagged as matching the reactant→product motion
  (overlap -0.60); IRC connected cleanly to both minima; and the
  reported barrier (4.31 kcal/mol at PM6) landed in a believable
  ballpark next to the experimental NH3 inversion barrier (~5.8
  kcal/mol) for a semi-empirical method -- a real sanity check, not
  just "the code didn't crash."
- The *full* `cli.py` flow (all five stages back to back) on that
  same NH3 case, and separately on HCN⇌HNC (up through the TS stage)
  to confirm the degenerate-alignment guard fails cleanly instead of
  crashing g09 -- `level_select`'s and `energetics`'s prompts
  monkeypatched in both cases, since `questionary`'s arrow-key picker
  needs a real TTY.
- `job_runner.build_pbs_script()`'s output checked character for
  character against the site's actual PBS template.
- The cycle-bump repair fix (continue from the geometry actually
  reached, not the original guess): forced non-convergence on a
  distorted water guess, confirmed the `.com` geometry genuinely
  differs between `_try2`/`_try3` now, where it used to be identical.
- GenECP end to end for real: a mixed H/O basis (H at STO-3G, O with a
  LanL2DZ ECP) on water, generated `.com` matching a real reference
  organometallic GenECP input's format exactly, ran and converged.
  The interactive element-assignment prompt itself (checkbox which
  elements need an ECP, name each, group by shared name) tested with
  mocked `questionary` responses -- produces the right groupings.
- The restructured `cli.py` (study-type choice, single-opt path,
  opt-in IRC, repeat-or-exit loop): single-optimization path run for
  real; the reaction-mechanism path re-run for real post-restructure
  and reproduced the exact same NH3 numbers as before; the interactive
  study-type prompt and the repeat-loop's exit condition both
  exercised with mocked prompts.
- The GUI, in an actual Chromium browser (not just its FastAPI
  endpoints via `TestClient`, though those were checked first): real
  `.xyz` file uploads through the real file input (simulated via
  `DataTransfer`, the standard way to drive a file input without a
  native OS picker), the GenECP checkbox flow producing the exact
  `ecp_assignments`/`light_basis` payload the backend expects, the
  study-type toggle showing/hiding the right sections both directions,
  and a complete mechanism study run start to finish through the real
  UI -- submit, live-updating progress log, a rendered results panel
  (HF/STO-3G NH3 inversion, 9.53 kcal/mol barrier, IRC "connected").
- Three real studies (two vinyl-alcohol⇌acetaldehyde mechanism runs
  at B3LYP/6-31G(d), one water single-point optimization) submitted
  and run concurrently through "+ Start another study," each getting
  its own live-updating card with independent stage text and log --
  confirmed genuinely parallel (all three progressing through
  different stages at once, not queued), and each finished with its
  own correct results. One of the mechanism runs had both kcal/mol
  and eV checked as energy units; the rendered barrier in each
  (64.85 kcal/mol, 2.81 eV) matched the expected conversion factor.
- `/api/refine`'s full request/response wiring, both via `TestClient`
  (a mocked `llm_assist.refine_structure()` returning a canned
  suggestion, with the downstream PubChem lookup and formula check
  running for real -- confirmed a matching-formula suggestion
  ("vinyl alcohol") gets offered and a mismatched one ("benzene")
  gets correctly rejected) and in the real browser (no credential
  configured here, so the request genuinely fails and the real
  Anthropic SDK error -- "Could not resolve authentication method"
  -- surfaces gracefully in the UI instead of hanging or crashing).
- `verification.py`, run for real on the vinyl-alcohol⇌acetaldehyde
  case at B3LYP/6-31G(d), through both entry points:
  - The GUI, with TS-mode distortion, IRC (15 points, deliberately
    short to keep the test fast), and IRC-endpoint reopt all enabled
    together: both distortion directions converged and were correctly
    classified (`distort+1.00` -> product side, `distort-1.00` ->
    reactant side, by a clear RMSD margin); the IRC itself correctly
    flagged as not clearly connected at only 15 points; its forward
    and reverse endpoints both reoptimized to the *same* minimum
    (identical energy and RMSDs) -- exactly the scenario
    `CLASSIFY_AMBIGUOUS_MARGIN` exists for, and it correctly reported
    neither candidate rather than guessing one; the candidate
    comparison correctly kept the original guess-based reactant/
    product for both sides, since in this case they were already more
    stable than the TS-distortion candidates -- confirming
    `select_best_candidate()` doesn't just prefer the new candidates
    by default. Also unit-tested `renderStudyResults()` directly with
    a synthetic non-"guess" winner to confirm the source-note UI rows
    render correctly (not otherwise exercised live, since "guess" won
    in the real run).
  - The CLI (`cli._run_mechanism_study`, prompts monkeypatched,
    TS-distortion answered yes / IRC answered no): same reaction,
    same numbers as the GUI run (as expected -- identical underlying
    pipeline), comparison line printed clearly before the final energy
    report.
  - `distort_along_imaginary_mode()` and `select_best_candidate()`
    also isolated-tested against toy structures/fake `OptOutcome`s to
    confirm the `x = x0 ± factor·mode` displacement and the
    lowest-energy-wins selection logic directly, without needing a
    full Gaussian run for basic sanity checking.
- `parser.last_convergence_status()` tested against a real, already-
  on-disk log with a genuinely stuck optimization
  (`jobs/diels-2/diels-2_reactant_hf_sto3g_opt_try3.log` -- 28 internal
  steps, Force converged, Displacement never trending down): parses
  the real four-row table correctly (1/4 criteria met), and
  `pipeline._repair_strategy()` correctly returns `"escalate"` for it;
  a synthetic 3-of-4 case confirmed `"continue"` the other way.
- Stop/Restart, run for real in the GUI: started a mechanism study,
  confirmed via `ps aux` that a real `g09` process (plus its `l1.exe`/
  `l1110.exe` children) was actually running, clicked Stop, confirmed
  via `ps aux` again that every one of those processes was gone (not
  just that the browser stopped polling), confirmed the card showed
  "Cancelled" with the log intact up to that point. Clicked Restart:
  the new run's log showed "Reusing already-completed ..." for the
  three stages that had genuinely finished before cancellation
  (reactant pre-opt, product pre-opt, reactant final-level opt) and
  for-real reran the one that hadn't (product final-level opt,
  reconverging to the same correct SCF energy as an earlier from-
  scratch run of the identical reaction), then completed TS search,
  IRC, and energetics normally from there -- a full, real validation
  of the checkpoint-skip logic, not just that the button changes
  state.
- The study-entry delete/bin button: ran a real single-optimization
  study to completion in the GUI, confirmed the trash icon was
  disabled while running and enabled once done, clicked it, confirmed
  the card disappeared and the whole Studies panel hid itself (it was
  the only card) -- then confirmed on disk (`ls jobs/<name>/`) that
  every `.com`/`.log` file the study produced was still there,
  untouched; only the GUI's own entry was gone.
- Endpoint recovery: `_assign_sides(..., reactant_trusted=False)`
  isolated-tested against toy structures (one candidate placed close
  to a trusted product reference, one far) to confirm elimination
  classification correctly labels the closer one as the trusted side
  and the farther one as the recovered side. Run for real in the GUI
  with a deliberately unoptimizable reactant (a 3-atom `.xyz` with all
  atoms within 0.03 Å of each other, uploaded directly via
  `POST /api/analyze` to bypass RDKit's embedding, which would never
  produce something this broken): confirmed PM6 pre-opt failed with a
  real g09 error on every fallback rung, confirmed recovery correctly
  logged skipping the final-level reopt and continuing with the
  product, confirmed the product optimized normally and unaffected,
  and confirmed TS search proceeded (using the reactant's best failed
  attempt as its reference) rather than aborting -- it went on to fail
  there too, expected, since a genuinely nonsensical starting geometry
  doesn't have a real TS to find with the intact product. A second run
  with the recovery checkbox on but no actual failure (normal
  vinyl-alcohol⇌acetaldehyde case) reproduced byte-identical numbers to
  earlier confirmed-correct runs of that reaction, confirming the new
  flag adds no regression when it isn't needed. See Known limitations
  for what this did and didn't cover.
- The `%chk=` naming fix, `thermochem.run_freqchk()`, `xyz_export.convert_log_to_xyz()`,
  skip-PM6, and the two-table results UI, all validated together in one
  real GUI run (vinyl-alcohol⇌acetaldehyde, skip-PM6 + export-XYZ +
  G-corr-at-25°C all enabled): confirmed no PM6 `.com`/`.log` files were
  written at all (skip genuinely skipped, not just logged) and the
  final-level optimizations ran directly on the raw input structures;
  confirmed each of reactant/product/TS got its own uniquely-named
  `.chk` sitting right next to its `.com`/`.log` (the naming fix
  working in the full pipeline, not just the isolated test below);
  confirmed all three `.xyz` files were produced with correct atom
  counts/coordinates; confirmed 25°C was correctly converted to 298.15
  K (`gcorr_temperature_k` in the result) and G-corr values matched
  hand-running `freqchk` on the same checkpoint; confirmed the Energy
  and Files tables rendered with the right headers/rows/units
  (including a G-corr column per requested unit) and the reaction plot
  sat genuinely beside the tables (same top offset, non-overlapping,
  confirmed via `getBoundingClientRect()`) rather than just below them.
  `run_freqchk()` and `convert_log_to_xyz()` were also each isolated-
  tested directly against a real `.chk`/`.log` first, cross-checked
  against running `freqchk`/`obabel` by hand on the same files, before
  ever wiring them into the pipeline.
- ChemDraw I/O, the print button, and the shorter filenames: before
  writing `from_chemdraw()`, did an actual round trip by hand (a real
  structure → `.cdxml` via `obabel` → back to `.mol` via `obabel` →
  loaded in RDKit) to confirm what ChemDraw files actually carry --
  bonds/bond-orders/formula came through correctly (`C2H4O` preserved)
  but the coordinates were flat 2D drawing-canvas positions (`z=0` on
  every atom), not real geometry, which is what led to re-embedding in
  3D rather than trusting them. Then tested the real function the same
  way: built a structure, exported it to `.cdxml`, loaded it back
  through `from_chemdraw()`, confirmed the formula matched and the
  geometry was genuine 3D (nonzero z, Å-scale bond lengths, not the
  flat input). `pytest tests/` (18 tests) passed after all the changes.
  In the real GUI: uploaded a `.cdxml` reactant through the actual file
  input (simulated via a real `change` event with a `DataTransfer`-
  constructed `File`, since this environment's browser-automation tool
  can't drive a native OS file-picker dialog) and confirmed it parsed
  to the right atom count through the real `/api/analyze` call; ran a
  full mechanism study (HF/STO-3G vinyl-alcohol⇌acetaldehyde,
  reactant loaded from that `.cdxml`) through to completion and
  confirmed every file under `jobs/<name>/` used the new short naming
  (`reactant_pm6opt.com`, not `<job>_reactant_pm6opt.com`) with each
  `.chk` still sitting correctly next to its `.com`/`.log`; separately
  called `convert_log_to_xyz()`/`convert_log_to_cdxml()` directly
  against that run's real `reactant_final.log` and confirmed both
  produced correctly-named files; rendered the results panel with a
  mocked result carrying both `xyz_files` and `cdxml_files` and
  confirmed the Files table showed all three columns (Logs/XYZ/CDXML)
  with the right paths, and that the print button's click handler
  correctly grabbed the live `<svg>` and built a black-on-white,
  non-dark-theme print document from it. (The mechanism study itself
  had to be re-run once after the dev server it was running against
  was unexpectedly killed mid-run by an unrelated environment hiccup
  -- confirmed via the surviving `.log` files that everything up to
  that point, including the shortened naming, had been correct before
  switching to the more targeted verification above.)

## Next up

- Wiring `job_runner.py` into `cli.py`/`pipeline.py`/`ts_search.py`/
  `irc.py` behind a `--executor {local,pbs}` switch, once there's a
  real PBS system to validate submission against -- `job_runner.py`
  itself is built and unit-tested (fake `qsub`/`qstat`) but has never
  actually submitted a job. The one thing to carry over from
  `local_runner.py`'s early lesson: a PBS batch script isn't an
  interactive shell, so the generated script sources `g09.profile`
  itself explicitly rather than assuming it's already on `PATH` (the
  "l1.exe not found" segfault this session ran into at the very
  start, before that was fixed).
- Only Opt Freq is built as a job type. `prompt_for_study_type()` and
  the single-opt path exist specifically as the branch point for
  adding others later (single-point energy, NMR, relaxed/rigid
  scans) without disturbing the mechanism-study path -- not built yet,
  just deliberately left room for.
- The GUI's progress panel is stage-level, not attempt-level (see
  Known limitations) -- worth revisiting if repair-loop retries turn
  out to need finer-grained live feedback in practice.
- GenECP through the GUI, and non-interactive/scripted GenECP through
  the CLI, for the same reason noted above.
- A real smoke test of `llm_assist.py` against a live Claude API
  credential (see Known limitations) -- the wiring is validated, the
  model's actual suggestion quality isn't yet.
- The convergence-driven repair heuristic (`pipeline._repair_strategy()`)
  needs runs across more real difficult geometries to know its actual
  success rate -- it's grounded in one confirmed-real stuck case, not
  tuned against a broad sample yet.
- Stop/Restart's `cancel_event` plumbing is already threaded through
  every pipeline function (`pipeline.py`/`ts_search.py`/`irc.py`/
  `verification.py` all accept it) -- wiring an equivalent Ctrl+C-safe
  cleanup into `cli.py` would be straightforward if a killed CLI run
  leaving an orphaned `g09` process ever turns out to matter in
  practice.
- A true Gaussian checkpoint (`%chk`/`Opt=Restart`) warm restart, if
  stage-level restart (see Known limitations) turns out not to save
  enough time in practice on expensive stages cancelled near the end.
