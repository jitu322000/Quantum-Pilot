"""
cli.py

The actual entry point a person runs.

One flow: get a geometry (optimized Gaussian log / existing GAMESS
input / guess geometry optimized via gaussbot), run RHF, then
optionally CIS on top of it. CASSCF/XMCQDPT and UHF/ROHF are future
work (see gamess_input.py/rhf.py docstrings).

Interactive (`gamessbot` with no args) prompts for everything. Given
--gaussian-log / --gamess-inp / --smiles / --file it loads the geometry
from that instead of prompting, for scripting; the rest (basis, charge/
mult, SOSCF, CIS, GAMESS settings) is still prompted unless also given.

Run:
  gamessbot --gaussian-log water_opt.log -n water --gbasis "GBASIS=STO NGAUSS=3" --cis --nstate 5
  gamessbot                                              # fully interactive
"""

from __future__ import annotations

import argparse
import os
from typing import List, Optional

from .active_space import ActiveSpaceSuggestion, build_active_space, format_active_space_summary
from .casscf import (
    mo_source_from_casscf_outcome,
    run_casscf_staged,
    run_casscf_with_smaller_active_space_recovery,
    suggest_active_space_from_cis,
)
from .cis import run_cis
from .energetics import build_energy_table, format_latex_table
from .gamess_input import data_block_from_gamess_inp, data_block_from_gaussian_log
from .orbital_character import parse_orbital_character
from .intake import optimize_and_convert, prompt_for_geometry_source
from .level_select import (
    DEFAULT_RUNGMS_PATH,
    DEFAULT_SCRATCH_DIR,
    prompt_for_another_active_space_combo,
    prompt_for_casscf,
    prompt_for_casscf_mem_mwords,
    prompt_for_casscf_nstate,
    prompt_for_charge_mult,
    prompt_for_cis,
    prompt_for_executor,
    prompt_for_gamess_settings,
    prompt_for_gbasis,
    prompt_for_mo_source,
    prompt_for_pbs_template,
    prompt_for_soscf,
    prompt_for_transitn,
    prompt_for_xmcqdpt,
)
from .rhf import run_rhf_staged
from .transitn import run_transitn
from .xmcqdpt import run_xmcqdpt


def _confirm_active_space(suggestion: ActiveSpaceSuggestion, rhf_log_path: str) -> ActiveSpaceSuggestion:
    """Shows the suggested active space and lets the user accept it or
    replace it with their own choice of MOs -- per your request that the
    system propose an active space and ask before running with it."""
    import questionary

    print("\n" + format_active_space_summary(suggestion))
    active_mos = suggestion.occ_selected + suggestion.virt_selected
    if questionary.confirm(
        "Show orbital character (which atom/orbital each active MO is dominated by)?", default=False
    ).ask():
        character = parse_orbital_character(rhf_log_path, active_mos)
        for mo in active_mos:
            c = character.get(mo)
            print(f"  MO {mo}: {c.label} (|c|={c.coefficient:.4f})" if c else f"  MO {mo}: (not found in log)")
    use_suggested = questionary.confirm("Use this active space?", default=True).ask()
    if use_suggested or use_suggested is None:
        return suggestion

    occ_str = questionary.text(
        "Active occupied MOs (comma-separated):", default=",".join(map(str, suggestion.occ_selected))
    ).ask()
    virt_str = questionary.text(
        "Active virtual MOs (comma-separated):", default=",".join(map(str, suggestion.virt_selected))
    ).ask()
    occ = [int(x) for x in (occ_str or "").split(",") if x.strip()]
    virt = [int(x) for x in (virt_str or "").split(",") if x.strip()]
    return build_active_space(suggestion.n_occ, suggestion.norb, occ, virt)


def _maybe_write_pbs_template(out_dir: str, executor: str, pbs_template: Optional[str]) -> None:
    """If PBS mode is on and a template was picked (prompt_for_pbs_template()),
    write it once as <out_dir>/pbs_template.sh -- job_runner.run_pbs()
    picks it up automatically from there for every stage in this job."""
    if executor == "pbs" and pbs_template:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "pbs_template.sh"), "w") as f:
            f.write(pbs_template)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="gamessbot -- GAMESS RHF/CIS studies, starting from an optimized Gaussian "
        "log, an existing GAMESS input, or a guess geometry optimized with Gaussian first."
    )
    parser.add_argument("--gaussian-log", help="Already-optimized Gaussian .log file")
    parser.add_argument("--gamess-inp", help="Existing GAMESS input/punch file (has its own $DATA)")
    parser.add_argument("--smiles", help="SMILES for a guess geometry (optimized with Gaussian first)")
    parser.add_argument("--file", help="Geometry file (.xyz/.mol/.sdf/.pdb) for a guess geometry")
    parser.add_argument("-n", "--name", help="Job name (job folder under jobs/)")
    parser.add_argument("--charge", type=int, help="Molecular charge")
    parser.add_argument("--mult", type=int, help="Spin multiplicity")
    parser.add_argument("--gbasis", help='GBASIS line, e.g. "GBASIS=STO NGAUSS=3"')
    parser.add_argument("--diis", action="store_true", help="Use DIIS instead of SOSCF for RHF")
    parser.add_argument("--cis", action="store_true", help="Run CIS on top of RHF")
    parser.add_argument("--nstate", type=int, help="Number of CIS states (NSTATE), needs --cis")
    parser.add_argument("--casscf", action="store_true", help="Run CASSCF on top of RHF/CIS, needs --cis")
    parser.add_argument("--casscf-nstate", type=int, help="Number of CASSCF states (state-averaged), needs --casscf")
    parser.add_argument("--xmcqdpt", action="store_true", help="Run XMCQDPT on top of the converged CASSCF, needs --casscf")
    parser.add_argument("--rungms", default=None, help=f"Path to rungms (default: {DEFAULT_RUNGMS_PATH})")
    parser.add_argument("--scratch-dir", default=None, help=f"GAMESS scratch directory (default: {DEFAULT_SCRATCH_DIR})")
    parser.add_argument("--ncpus", type=int, default=None, help="Number of CPUs (default: 1)")
    parser.add_argument("--mem-mwords", type=int, default=None, help="Memory in MWORDS (default: 1)")
    parser.add_argument("--nprocs", type=int, default=2, help="Gaussian %%nprocshared for the guess-geometry optimization")
    parser.add_argument("--mem-gb", type=int, default=2, help="Gaussian %%mem (GB) for the guess-geometry optimization")
    parser.add_argument("--gaussian-method", default="B3LYP", help="Gaussian method for the final-level optimization of a guess geometry (after a PM6 pre-opt)")
    parser.add_argument("--gaussian-basis", default="6-31G(d)", help="Gaussian basis set for the final-level optimization of a guess geometry")
    parser.add_argument(
        "--executor", choices=["local", "pbs"], default=None,
        help="Where to actually run GAMESS jobs: 'local' (rungms directly) or 'pbs' (submit via "
        "qsub) -- prompted if not given for an interactive run, defaults to 'local' for a scripted one.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    print("gamessbot -- GAMESS RHF/CIS studies\n")

    scripted_geometry = args.gaussian_log or args.gamess_inp or args.smiles or args.file
    job_name = args.name or "job"
    out_dir = os.path.join("jobs", job_name)
    os.makedirs(out_dir, exist_ok=True)

    # --executor overrides; otherwise a bare interactive invocation asks
    # once per session (not per stage) -- a scripted invocation (one of
    # --gaussian-log/--gamess-inp/--smiles/--file given) just defaults
    # to "local" rather than prompting. Same split for the PBS script
    # template: only offered interactively (it can block on $EDITOR),
    # never for scripted/automated use -- scripted PBS submissions fall
    # back to job_runner.py's own DEFAULT_* constants instead.
    interactive = not scripted_geometry
    executor = args.executor
    if executor is None:
        executor = prompt_for_executor() if interactive else "local"
    pbs_template = prompt_for_pbs_template() if (interactive and executor == "pbs") else None

    if scripted_geometry:
        if args.gaussian_log:
            data_block = data_block_from_gaussian_log(args.gaussian_log, title=job_name)
        elif args.gamess_inp:
            data_block = data_block_from_gamess_inp(args.gamess_inp, title=job_name)
        else:
            from gaussbot.geometry import from_file, from_smiles

            structure = (
                from_smiles(args.smiles, label=job_name)
                if args.smiles
                else from_file(args.file, label=job_name)
            )
            print(f"Optimizing the guess geometry with Gaussian (PM6 pre-opt, then {args.gaussian_method}/{args.gaussian_basis}) before handing off to GAMESS:")
            data_block, pre, final = optimize_and_convert(
                structure, out_dir, job_name, method=args.gaussian_method, basis=args.gaussian_basis,
                nprocs=args.nprocs, mem_gb=args.mem_gb,
            )
            for line in pre.log:
                print(f"  {line}")
            for line in final.log:
                print(f"  {line}")
    else:
        job_name = args.name or input("\nName for this job (used for the job folder): ").strip() or "job"
        out_dir = os.path.join("jobs", job_name)
        os.makedirs(out_dir, exist_ok=True)
        data_block = prompt_for_geometry_source(out_dir, job_name, nprocs=args.nprocs, mem_gb=args.mem_gb)

    _maybe_write_pbs_template(out_dir, executor, pbs_template)

    charge = args.charge if args.charge is not None else None
    mult = args.mult if args.mult is not None else None
    if charge is None or mult is None:
        prompted_charge, prompted_mult = prompt_for_charge_mult()
        charge = charge if charge is not None else prompted_charge
        mult = mult if mult is not None else prompted_mult

    gbasis_line = args.gbasis or prompt_for_gbasis()
    use_soscf = not args.diis if args.diis else prompt_for_soscf()

    if args.cis:
        run_cis_flag, nstate = True, args.nstate or 5
    else:
        run_cis_flag, nstate = prompt_for_cis()

    if args.rungms and args.scratch_dir:
        rungms_path, scratch_dir = args.rungms, args.scratch_dir
        ncpus = args.ncpus or 1
        mem_mwords = args.mem_mwords or 1
    else:
        rungms_path, scratch_dir, ncpus, mem_mwords = prompt_for_gamess_settings()
        rungms_path = args.rungms or rungms_path
        scratch_dir = args.scratch_dir or scratch_dir
        ncpus = args.ncpus or ncpus
        mem_mwords = args.mem_mwords or mem_mwords

    print("\nRunning RHF:")
    rhf_outcome = run_rhf_staged(
        data_block, out_dir, charge=charge, mult=mult, gbasis_line=gbasis_line,
        rungms_path=rungms_path, scratch_dir=scratch_dir, ncpus=ncpus, mem_mwords=mem_mwords,
        use_soscf=use_soscf, executor=executor,
    )
    for line in rhf_outcome.trail:
        print(f"  {line}")

    if not rhf_outcome.success:
        print(f"\nRHF did not converge. See {rhf_outcome.log_path} for details.")
        return

    print(f"\nOK: RHF energy = {rhf_outcome.energy} Hartree (NORB={rhf_outcome.norb})")
    print(f"Log file: {rhf_outcome.log_path}")

    if not run_cis_flag:
        return

    print(f"\nRunning CIS ({nstate} states):")
    cis_outcome = run_cis(
        rhf_outcome, out_dir, nstate=nstate,
        rungms_path=rungms_path, scratch_dir=scratch_dir, ncpus=ncpus, mem_mwords=mem_mwords,
        executor=executor,
    )
    if not cis_outcome.success:
        print(f"\nCIS did not terminate normally. See {cis_outcome.log_path} for details.")
        return

    print(f"\nOK: {len(cis_outcome.states)} excited states found. Log file: {cis_outcome.log_path}")
    for state in cis_outcome.states:
        top = max(state.transitions, key=lambda t: abs(t.coefficient)) if state.transitions else None
        top_desc = f", dominant: {top.from_mo}->{top.to_mo} ({top.coefficient:.4f})" if top else ""
        print(f"  State {state.index}: E = {state.energy:.6f} Hartree, S = {state.spin}, sym = {state.space_sym}{top_desc}")

    if not (args.casscf if args.casscf else prompt_for_casscf()):
        return

    try:
        suggestion = suggest_active_space_from_cis(cis_outcome, rhf_outcome)
    except ValueError as e:
        print(f"\nCouldn't suggest an active space: {e}")
        return

    # Per your request: let the user compare several active-space/state
    # combinations against each other, not just run one -- each combo is
    # independently confirmed/adjusted, named cas-{e}{o}-sa{n} (and
    # xpt-.../trn-... for its own XMCQDPT/TRANSITN follow-ups, restarting
    # from THAT combo's own converged orbitals), with its own energy
    # table -- all combined into one LaTeX file at the end.
    combo_tables: List[str] = []
    combo_index = 0
    last_casscf_outcome = None
    while True:
        combo_index += 1
        print(f"\n--- Active-space/state combination #{combo_index} ---")
        active_space = _confirm_active_space(suggestion, rhf_outcome.log_path)

        casscf_nstate = (args.casscf_nstate or prompt_for_casscf_nstate(nstate)) if combo_index == 1 else prompt_for_casscf_nstate(nstate)
        casscf_mem_mwords = prompt_for_casscf_mem_mwords()

        # Only relevant from combination #2 onward -- the first combo has
        # no previous CASSCF result to continue from, so it always starts
        # fresh from the closed-shell (RHF) orbitals.
        if combo_index > 1 and last_casscf_outcome is not None:
            mo_source = mo_source_from_casscf_outcome(last_casscf_outcome) if prompt_for_mo_source() == "previous" else rhf_outcome
        else:
            mo_source = rhf_outcome

        n_electrons = 2 * active_space.ndoc
        n_orbitals = active_space.ndoc + active_space.nval
        cas_prefix = f"cas-{n_electrons}{n_orbitals}-sa{casscf_nstate}-c{combo_index}"

        print(f"\nRunning CASSCF ({casscf_nstate} states, NMCC={active_space.nmcc} "
              f"NDOC={active_space.ndoc} NVAL={active_space.nval}):")
        staged = run_casscf_staged(
            mo_source, active_space, out_dir, nstate=casscf_nstate,
            rungms_path=rungms_path, scratch_dir=scratch_dir, ncpus=ncpus, mem_mwords=casscf_mem_mwords,
            executor=executor, name_prefix=cas_prefix,
        )
        for line in staged.trail:
            print(f"  {line}")
        casscf_outcome = staged.outcome

        if not casscf_outcome.success and staged.exhausted:
            import questionary

            if questionary.confirm(
                "\nCASSCF still hasn't converged after more iterations. Try a smaller/easier "
                "active space first, then use its orbitals to restart this one?", default=True,
            ).ask():
                max_electrons = int((questionary.text(
                    "Smaller active space: max electrons:", default=str(min(4, 2 * active_space.ndoc)),
                ).ask() or "4").strip())
                max_orbitals = int((questionary.text(
                    "Smaller active space: max orbitals:", default=str(min(4, active_space.ndoc + active_space.nval)),
                ).ask() or "4").strip())
                print("\nRunning the smaller-active-space recovery:")
                casscf_outcome, active_space, recovery_trail = run_casscf_with_smaller_active_space_recovery(
                    mo_source, active_space, out_dir, nstate=casscf_nstate,
                    rungms_path=rungms_path, scratch_dir=scratch_dir,
                    smaller_max_electrons=max_electrons, smaller_max_orbitals=max_orbitals,
                    ncpus=ncpus, mem_mwords=casscf_mem_mwords, executor=executor, name_prefix=cas_prefix,
                )
                for line in recovery_trail:
                    print(f"  {line}")

        if not casscf_outcome.success:
            if casscf_outcome.result is not None and casscf_outcome.result.normal_termination:
                note = (" (ran to completion but the orbital optimization did not converge -- try a "
                         "different active space or a different starting geometry)")
            else:
                note = " (GAMESS exited abnormally -- check the log, e.g. for an unviable active space)"
            print(f"\nCASSCF did not succeed{note}. See {casscf_outcome.log_path} for details.")
            if not prompt_for_another_active_space_combo():
                break
            continue

        print(f"\nOK: CASSCF converged. Log file: {casscf_outcome.log_path}")
        for index, energy in casscf_outcome.result.state_energies:
            print(f"  State {index}: E = {energy:.6f} Hartree")
        last_casscf_outcome = casscf_outcome

        # active_space may have changed (regrown) after recovery -- the
        # tag/prefix for downstream stages must reflect what casscf_outcome
        # actually used, not the originally-requested combo.
        n_electrons = 2 * active_space.ndoc
        n_orbitals = active_space.ndoc + active_space.nval
        combo_tag = f"{n_electrons}{n_orbitals}-sa{casscf_nstate}-c{combo_index}"

        xmcqdpt_outcome = None
        if args.xmcqdpt if args.xmcqdpt else prompt_for_xmcqdpt():
            print(f"\nRunning XMCQDPT ({casscf_nstate} states):")
            xmcqdpt_outcome = run_xmcqdpt(
                casscf_outcome, out_dir, nstate=casscf_nstate,
                rungms_path=rungms_path, scratch_dir=scratch_dir, ncpus=ncpus, mem_mwords=casscf_mem_mwords,
                executor=executor, name=f"xpt-{combo_tag}",
            )
            if not xmcqdpt_outcome.success:
                print(f"\nXMCQDPT did not succeed. See {xmcqdpt_outcome.log_path} for details.")
                xmcqdpt_outcome = None
            else:
                print(f"\nOK: XMCQDPT complete. Log file: {xmcqdpt_outcome.log_path}")
                for index, energy in xmcqdpt_outcome.result.mcqdpt_state_energies:
                    print(f"  State {index}: E = {energy:.6f} Hartree")

        oscillator_strengths = None
        if prompt_for_transitn():
            print(f"\nRunning TRANSITN (oscillator strengths, {casscf_nstate} states):")
            transitn_outcome = run_transitn(
                casscf_outcome, out_dir, nstate=casscf_nstate,
                rungms_path=rungms_path, scratch_dir=scratch_dir, ncpus=ncpus, mem_mwords=casscf_mem_mwords,
                executor=executor, name=f"optical-{combo_tag}",
            )
            if not transitn_outcome.success:
                print(f"\nTRANSITN did not succeed. See {transitn_outcome.log_path} for details.")
            else:
                print(f"\nOK: TRANSITN complete. Log file: {transitn_outcome.log_path}")
                oscillator_strengths = transitn_outcome.oscillator_strengths
                for state_index, f_value in sorted(oscillator_strengths.items()):
                    print(f"  State {state_index} (S{state_index - 1}): f = {f_value:.4f}")

        rows = build_energy_table(
            casscf_outcome.result, xmcqdpt_outcome.result if xmcqdpt_outcome else None, oscillator_strengths,
        )
        print(f"\nVertical excitation energies (eV, relative to S0) -- CAS({n_electrons},{n_orbitals}), SA{casscf_nstate}:")
        header = f"  {'State':<6}{'CASSCF (eV)':>14}{'XMCQDPT (eV)':>15}"
        if oscillator_strengths:
            header += f"{'f':>10}"
        print(header)
        for row in rows:
            xmcqdpt_str = f"{row.xmcqdpt_ev:.2f}" if row.xmcqdpt_ev is not None else "--"
            row_str = f"  {row.label:<6}{row.casscf_ev:>14.2f}{xmcqdpt_str:>15}"
            if oscillator_strengths:
                f_str = f"{row.oscillator_strength:.4f}" if row.oscillator_strength is not None else "--"
                row_str += f"{f_str:>10}"
            print(row_str)

        combo_tables.append(format_latex_table(rows, nstate=casscf_nstate, active_space=active_space))

        if not prompt_for_another_active_space_combo():
            break

    if combo_tables:
        import questionary

        prompt = (
            f"\nSave all {len(combo_tables)} table(s) as a LaTeX file?"
            if len(combo_tables) > 1 else "\nSave this table as a LaTeX file?"
        )
        if questionary.confirm(prompt, default=False).ask():
            tex_filename = "energy_tables.tex" if len(combo_tables) > 1 else "energy_table.tex"
            tex_path = os.path.join(out_dir, tex_filename)
            with open(tex_path, "w") as f:
                f.write("\n\n".join(combo_tables))
            print(f"Wrote {tex_path}")


if __name__ == "__main__":
    main()
