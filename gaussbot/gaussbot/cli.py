"""
cli.py

The actual entry point a person runs.

First asks what this run is for (level_select.prompt_for_study_type,
skipped if --reactant/--product already say which): a single geometry
optimization, or a full reaction-mechanism study. Both share the same
resource + method/basis (including GenECP) prompts and the same
resilient PM6-pre-opt-then-final-opt machinery; the mechanism study
additionally searches for and verifies a TS, optionally runs the IRC,
and reports energetics.

Interactive (`gaussbot` with no args) loops, offering to run another
study, until you say no. Given --reactant/--product it runs once and
exits, for scripting.

Run:
  gaussbot -r reactant.com -p product.xyz [-n job_name]   # mechanism study
  gaussbot -r structure.com [-n job_name]                  # single opt
  gaussbot                                                   # interactive, asks which
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Optional

import questionary

from .intake import (
    load_structure,
    prompt_for_structure,
    prompt_for_reaction,
    prompt_for_ts_guess,
    prompt_for_ts_search_inputs,
    check_reaction_match,
)
from .level_select import (
    LevelChoice,
    prompt_for_method_basis,
    prompt_for_resources,
    prompt_for_irc,
    prompt_for_study_type,
    prompt_for_executor,
    prompt_for_pbs_template,
    prompt_for_ts_distortion_check,
    prompt_for_irc_endpoint_reopt,
    prompt_for_max_repair_attempts,
    prompt_for_endpoint_recovery,
    prompt_for_skip_pm6,
    prompt_for_ts_search_skip_pm6,
    prompt_for_ts_calcall_fallback,
    prompt_for_ts_recovery,
    prompt_for_xyz_export,
    prompt_for_cdxml_export,
    prompt_for_gcorr,
)
from .pipeline import preopt_with_escalation, repair_and_optimize, OptOutcome
from .ts_search import run_ts_search, search_ts_staged, TSOutcome
from .irc import run_irc
from .verification import (
    run_ts_distortion_check, run_irc_endpoint_reopt, select_best_candidate,
    format_candidate_comparison, recover_endpoints_from_ts, DistortionOutcome,
)
from .energetics import prompt_for_energy_unit, build_energy_report, format_energy_report, HARTREE_TO
from .geometry import Structure
from .thermochem import run_freqchk, to_kelvin
from .xyz_export import convert_log_to_xyz, convert_log_to_cdxml


def _maybe_write_pbs_template(out_dir: str, executor: str, pbs_template: Optional[str]) -> None:
    """If PBS mode is on and a template was picked (prompt_for_pbs_template()),
    write it once as <out_dir>/pbs_template.sh -- job_runner.run_pbs()
    picks it up automatically for every job this study submits."""
    if executor == "pbs" and pbs_template:
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "pbs_template.sh"), "w") as f:
            f.write(pbs_template)


def _run_preopt(
    structure: Structure,
    out_dir: str,
    tag: str,
    nprocs: int,
    mem_gb: int,
    max_repair_attempts: int = 4,
    pubchem_query: Optional[str] = None,
    exit_on_failure: bool = True,
    executor: str = "local",
) -> OptOutcome:
    """Run the resilient PM6 (-> PubChem -> HF/STO-3G) pre-optimization
    for one structure. Exits the program if every rung of the ladder
    fails and `exit_on_failure` is set (the default) -- there's
    normally no sound way to proceed from a geometry that isn't a
    genuine minimum. With `exit_on_failure=False` (endpoint recovery),
    the failed OptOutcome is returned instead, so the caller can
    continue with the other side and recover this one from the TS."""
    print(f"  optimizing {tag} ...")
    outcome = preopt_with_escalation(
        structure, out_dir, tag, pubchem_query=pubchem_query, nprocs=nprocs, mem_gb=mem_gb,
        max_repair_attempts=max_repair_attempts, executor=executor,
    )
    for line in outcome.log:
        print(f"    {line}")
    if not outcome.success:
        print(f"\nCouldn't get {tag} to a clean minimum after every fallback. See the trail above.")
        if exit_on_failure:
            sys.exit(1)
        print(f"Recovery is enabled -- continuing without {tag}; it may be recovered from the TS afterward.")
        return outcome
    print(f"    OK: {tag} is a minimum at {outcome.method}/{outcome.basis or 'no basis'}, "
          f"SCF energy = {outcome.result.scf_energy} a.u.")
    return outcome


def _run_final_opt(
    structure: Structure, com_path: str, tag: str, level: LevelChoice, nprocs: int, mem_gb: int,
    max_repair_attempts: int = 4,
    exit_on_failure: bool = True,
    executor: str = "local",
) -> OptOutcome:
    """Reoptimize an already PM6-minimized structure at the user-chosen
    final level, through the same repair loop. Exits the program if it
    can't reach a clean minimum and `exit_on_failure` is set (the
    default) -- unlike the PM6 pre-opt stage, there's no further
    escalation to fall back on here. With `exit_on_failure=False`, the
    failed OptOutcome is returned instead (see _run_preopt)."""
    level_desc = f"{level.method}/genecp" if level.basis_groups is not None else f"{level.method}/{level.basis}"
    print(f"  optimizing {tag} at the final level ...")
    outcome = repair_and_optimize(
        structure, com_path, method=level.method, basis=level.basis,
        basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, max_repair_attempts=max_repair_attempts, executor=executor,
    )
    for line in outcome.log:
        print(f"    {line}")
    if not outcome.success:
        print(f"\n{tag} didn't reach a clean minimum at {level_desc}. See the trail above.")
        if exit_on_failure:
            sys.exit(1)
        print(f"Recovery is enabled -- continuing without {tag}; it may be recovered from the TS afterward.")
        return outcome
    print(f"    OK: {tag} is a minimum at {level_desc}, SCF energy = {outcome.result.scf_energy} a.u.")
    return outcome


def _run_ts(
    reactant: Structure,
    product: Structure,
    com_path: str,
    level: LevelChoice,
    nprocs: int,
    mem_gb: int,
    ts_guess: Optional[Structure],
    executor: str = "local",
) -> TSOutcome:
    """Run the TS search and exit the program if no guess produced a TS
    whose imaginary mode actually matches the reactant<->product motion."""
    print("  searching for the TS ...")
    outcome = run_ts_search(
        reactant, product, com_path, method=level.method, basis=level.basis,
        basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, ts_guess=ts_guess, executor=executor,
    )
    for line in outcome.log:
        print(f"    {line}")
    if not outcome.success:
        print("\nCouldn't find a TS whose imaginary mode matches the reaction. See the trail above.")
        sys.exit(1)
    print(
        f"    OK: TS found, imaginary freq {outcome.result.imaginary_freqs[0]:.1f} cm^-1, "
        f"overlap with reactant->product motion = {outcome.match_overlap:.2f}"
    )
    return outcome


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="gaussbot -- resilient geometry optimization, and (given a reactant and "
        "product) a full reaction-mechanism study: TS search, IRC verification, energetics."
    )
    parser.add_argument("-r", "--reactant", help="Reactant/structure file: .com, .xyz, .mol, .sdf, .pdb, .cdx, .cdxml, or .log (already-optimized)")
    parser.add_argument("-p", "--product", help="Product file (omit for a single geometry optimization)")
    parser.add_argument("-n", "--name", help="Job name (job folder under jobs/); default derived from --reactant")
    parser.add_argument(
        "--reactant-id",
        help="Chemical name or SMILES for the reactant, used only as a PubChem fallback "
        "query if PM6 can't reach a minimum from the geometry you gave it",
    )
    parser.add_argument("--product-id", help="Same as --reactant-id, for the product")
    parser.add_argument(
        "--ts-guess",
        help="Transition-state guess file (.com/.xyz/.mol/.sdf/.pdb/.cdx/.cdxml/.log). If omitted, one is "
        "generated automatically by interpolating between the optimized reactant and product",
    )
    parser.add_argument(
        "--executor", choices=["local", "pbs"], default=None,
        help="Where to actually run Gaussian jobs: 'local' (g09 directly) or 'pbs' (submit via "
        "qsub -- see job_runner.py, and edit its DEFAULT_G09ROOT/DEFAULT_INTEL_SETVARS/DEFAULT_PPN/"
        "DEFAULT_WALLTIME constants for your own cluster's paths/settings first). If omitted, a bare "
        "interactive `gaussbot` asks once per session; a scripted -r/-p invocation defaults to 'local'.",
    )
    return parser.parse_args()


def _run_single_study(args: argparse.Namespace, executor: str = "local", pbs_template: Optional[str] = None) -> None:
    if args.reactant:
        structure = load_structure(args.reactant, label="structure")
        job_name = args.name or os.path.splitext(os.path.basename(args.reactant))[0]
    else:
        structure = prompt_for_structure("structure")
        job_name = args.name or input(
            "\nName for this job (used for the job folder): "
        ).strip() or "job"

    out_dir = os.path.join("jobs", job_name)
    _maybe_write_pbs_template(out_dir, executor, pbs_template)
    nprocs, mem_gb = prompt_for_resources()
    max_repair_attempts = prompt_for_max_repair_attempts()
    skip_pm6 = prompt_for_skip_pm6()
    export_xyz = prompt_for_xyz_export()
    export_cdxml = prompt_for_cdxml_export()
    gcorr = prompt_for_gcorr()

    starting_structure = structure
    if skip_pm6:
        print("\nSkipping PM6 pre-optimization -- starting directly at the final level.")
    else:
        print("\nPre-optimizing (PM6):")
        pre = _run_preopt(structure, out_dir, "structure", nprocs, mem_gb, max_repair_attempts, pubchem_query=args.reactant_id, executor=executor)
        starting_structure = pre.structure

    elements = sorted({a[0] for a in starting_structure.atoms})
    print("\nNow pick the level for the real calculation:")
    level = prompt_for_method_basis(elements)

    print("\nOptimizing at the final level:")
    final = _run_final_opt(
        starting_structure, os.path.join(out_dir, "structure_final.com"), "structure", level, nprocs, mem_gb,
        max_repair_attempts, executor=executor,
    )

    print(f"\nDone. Optimized geometry and frequencies are in jobs/{job_name}/.")
    print(f"Log file: {final.log_path}")

    if export_xyz:
        xyz_path = convert_log_to_xyz(final.log_path)
        print(f"XYZ file: {xyz_path}" if xyz_path else "Couldn't convert to .xyz (is Open Babel installed?).")

    if export_cdxml:
        cdxml_path = convert_log_to_cdxml(final.log_path)
        print(f"CDXML file: {cdxml_path}" if cdxml_path else "Couldn't convert to .cdxml (is Open Babel installed?).")

    if gcorr is not None:
        temp, unit = gcorr
        g_corr_hartree = run_freqchk(final.log_path.replace(".log", ".chk"), to_kelvin(temp, unit))
        if g_corr_hartree is not None:
            print(f"G-corr at {temp} {unit}: {g_corr_hartree:.6f} Hartree")
        else:
            print("Couldn't compute G-corr (is freqchk installed/working?).")


def _run_mechanism_study(args: argparse.Namespace, executor: str = "local", pbs_template: Optional[str] = None) -> None:
    if args.reactant and args.product:
        reactant = load_structure(args.reactant, label="reactant")
        product = load_structure(args.product, label="product")
        check_reaction_match(reactant, product)
        job_name = args.name or os.path.splitext(os.path.basename(args.reactant))[0]
        ts_guess = load_structure(args.ts_guess, label="ts_guess") if args.ts_guess else None
    else:
        reactant, product = prompt_for_reaction()
        job_name = args.name or input(
            "\nName for this reaction (used for the job folder): "
        ).strip() or "reaction"
        ts_guess = prompt_for_ts_guess()

    out_dir = os.path.join("jobs", job_name)
    _maybe_write_pbs_template(out_dir, executor, pbs_template)

    r_elems = [a[0] for a in reactant.atoms]
    p_elems = [a[0] for a in product.atoms]
    if r_elems != p_elems:
        print(
            "\nReactant/product atom count or order didn't match (see "
            "warning above) -- fix that before continuing, a TS search needs "
            "the same atoms in the same order on both ends. Stopping."
        )
        return

    nprocs, mem_gb = prompt_for_resources()
    max_repair_attempts = prompt_for_max_repair_attempts()
    recover = prompt_for_endpoint_recovery()
    skip_pm6 = prompt_for_skip_pm6()
    export_xyz = prompt_for_xyz_export()
    export_cdxml = prompt_for_cdxml_export()
    gcorr = prompt_for_gcorr()

    r_pre = p_pre = None
    if skip_pm6:
        print("\nSkipping PM6 pre-optimization for both reactant and product -- starting directly at the final level.")
    else:
        print("\nPre-optimizing reactant and product (PM6):")
        r_pre = _run_preopt(
            reactant, out_dir, "reactant", nprocs, mem_gb, max_repair_attempts,
            pubchem_query=args.reactant_id, exit_on_failure=not recover, executor=executor,
        )
        p_pre = _run_preopt(
            product, out_dir, "product", nprocs, mem_gb, max_repair_attempts,
            pubchem_query=args.product_id, exit_on_failure=not recover, executor=executor,
        )
        if not r_pre.success and not p_pre.success:
            print("\nNeither the reactant nor the product reached a clean PM6 minimum -- nothing usable to search a TS from. Stopping.")
            return

    r_preopt_ok = r_pre is None or r_pre.success
    p_preopt_ok = p_pre is None or p_pre.success

    elements = sorted(
        {a[0] for a in (r_pre.structure if r_pre else reactant).atoms}
        | {a[0] for a in (p_pre.structure if p_pre else product).atoms}
    )
    print("\nNow pick the level for the real study:")
    level = prompt_for_method_basis(elements)

    print("\nReoptimizing reactant and product at the final level:")
    r_final = None
    if r_preopt_ok:
        r_final = _run_final_opt(
            r_pre.structure if r_pre else reactant, os.path.join(out_dir, "reactant_final.com"), "reactant",
            level, nprocs, mem_gb, max_repair_attempts, exit_on_failure=not recover, executor=executor,
        )
    p_final = None
    if p_preopt_ok:
        p_final = _run_final_opt(
            p_pre.structure if p_pre else product, os.path.join(out_dir, "product_final.com"), "product",
            level, nprocs, mem_gb, max_repair_attempts, exit_on_failure=not recover, executor=executor,
        )

    reactant_missing = r_final is None or not r_final.success
    product_missing = p_final is None or not p_final.success
    if reactant_missing and product_missing:
        print("\nNeither the reactant nor the product reached a clean minimum at the final level -- nothing usable to search a TS from. Stopping.")
        return

    reactant_for_ts = (r_final.structure if r_final else None) or (r_pre.structure if r_pre else None) or reactant
    product_for_ts = (p_final.structure if p_final else None) or (p_pre.structure if p_pre else None) or product
    if reactant_missing or product_missing:
        print(f"\nUsing {'the reactant' if reactant_missing else 'the product'}'s best available (non-converged) "
              "geometry as the TS-search reference for that side -- not independently verified as a minimum.")
    else:
        print("\nBoth reactant and product are clean minima at the final level.")

    print("Searching for the TS:")
    ts_com = os.path.join(out_dir, "ts.com")
    ts_outcome = _run_ts(reactant_for_ts, product_for_ts, ts_com, level, nprocs, mem_gb, ts_guess, executor=executor)

    reactant_trusted, product_trusted = not reactant_missing, not product_missing

    distortion_outcome = None
    if prompt_for_ts_distortion_check():
        print("\nCross-checking via TS-mode distortion:")
        distortion_outcome = run_ts_distortion_check(
            ts_outcome, reactant_for_ts, product_for_ts, out_dir, level, nprocs, mem_gb,
            reactant_trusted=reactant_trusted, product_trusted=product_trusted, executor=executor,
        )
        for line in distortion_outcome.log:
            print(f"  {line}")

    irc_points = prompt_for_irc()
    irc_outcome = None
    if irc_points is not None:
        print(f"\nRunning the IRC ({irc_points} points/direction):")
        irc_com = os.path.join(out_dir, "irc.com")
        irc_outcome = run_irc(
            ts_outcome.structure, reactant_for_ts, product_for_ts, irc_com,
            method=level.method, basis=level.basis, basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
            nprocs=nprocs, mem_gb=mem_gb, maxpoints=irc_points,
            reactant_trusted=reactant_trusted, product_trusted=product_trusted, executor=executor,
        )
        for line in irc_outcome.log:
            print(f"  {line}")
        if not irc_outcome.success:
            print(
                "\nIRC didn't clearly connect the TS to both the reactant and product -- "
                f"see {irc_com.replace('.com', '.log')} and the RMSDs above before trusting this TS."
            )
    else:
        print("\nSkipping IRC verification.")

    irc_reopt_outcome = None
    if irc_outcome is not None and irc_outcome.forward_structure is not None and irc_outcome.reverse_structure is not None:
        if prompt_for_irc_endpoint_reopt():
            print("\nReoptimizing the IRC endpoints:")
            irc_reopt_outcome = run_irc_endpoint_reopt(
                irc_outcome, reactant_for_ts, product_for_ts, out_dir, level, nprocs, mem_gb,
                reactant_trusted=reactant_trusted, product_trusted=product_trusted, executor=executor,
            )
            for line in irc_reopt_outcome.log:
                print(f"  {line}")

    reactant_candidates = {
        "guess": r_final if (r_final is not None and r_final.success) else None,
        "ts_distortion": distortion_outcome.reactant_candidate.outcome if distortion_outcome and distortion_outcome.reactant_candidate else None,
        "irc_reopt": irc_reopt_outcome.reactant_candidate.outcome if irc_reopt_outcome and irc_reopt_outcome.reactant_candidate else None,
    }
    product_candidates = {
        "guess": p_final if (p_final is not None and p_final.success) else None,
        "ts_distortion": distortion_outcome.product_candidate.outcome if distortion_outcome and distortion_outcome.product_candidate else None,
        "irc_reopt": irc_reopt_outcome.product_candidate.outcome if irc_reopt_outcome and irc_reopt_outcome.product_candidate else None,
    }

    unrecovered = []
    r_winner = p_winner = None
    try:
        r_label, r_winner = select_best_candidate(reactant_candidates)
    except ValueError:
        unrecovered.append("reactant")
    try:
        p_label, p_winner = select_best_candidate(product_candidates)
    except ValueError:
        unrecovered.append("product")
    if unrecovered:
        print(
            f"\nCouldn't get a usable {' or '.join(unrecovered)} for the energy report even with recovery "
            f"enabled -- the TS was found (imaginary freq {ts_outcome.result.imaginary_freqs[0]:.1f} cm^-1) but "
            f"neither the original optimization nor TS-mode distortion/IRC-endpoint reopt produced a usable "
            f"{' or '.join(unrecovered)}. Try enabling TS-mode distortion and/or IRC-endpoint reopt, or supply "
            "a better starting geometry. Stopping."
        )
        return

    print()
    print(format_candidate_comparison("reactant", reactant_candidates, r_label))
    print(format_candidate_comparison("product", product_candidates, p_label))
    reactant_result, product_result = r_winner.result, p_winner.result
    reactant_log_path, product_log_path = r_winner.log_path, p_winner.log_path

    unit = prompt_for_energy_unit()
    report = build_energy_report(reactant_result, ts_outcome.result, product_result, unit)
    print()
    print(format_energy_report(report))
    print()
    print("Log files these numbers came from:")
    print(f"  Reactant: {reactant_log_path}")
    print(f"  TS:       {ts_outcome.log_path}")
    print(f"  Product:  {product_log_path}")
    if irc_outcome is not None:
        print(f"  IRC:      {irc_outcome.log_path}")

    labeled_logs = [("Reactant", reactant_log_path), ("TS", ts_outcome.log_path), ("Product", product_log_path)]

    if export_xyz:
        print()
        print("XYZ files:")
        for label, log_path in labeled_logs:
            xyz_path = convert_log_to_xyz(log_path)
            print(f"  {label}: {xyz_path if xyz_path else 'conversion failed (is Open Babel installed?)'}")

    if export_cdxml:
        print()
        print("CDXML files:")
        for label, log_path in labeled_logs:
            cdxml_path = convert_log_to_cdxml(log_path)
            print(f"  {label}: {cdxml_path if cdxml_path else 'conversion failed (is Open Babel installed?)'}")

    if gcorr is not None:
        temp, temp_unit = gcorr
        temp_k = to_kelvin(temp, temp_unit)
        print()
        print(f"G-corr at {temp} {temp_unit} ({unit}):")
        for label, log_path in labeled_logs:
            g_corr_hartree = run_freqchk(log_path.replace(".log", ".chk"), temp_k)
            if g_corr_hartree is None:
                print(f"  {label}: couldn't compute (is freqchk installed/working?)")
            else:
                print(f"  {label}: {g_corr_hartree * HARTREE_TO[unit]:.6f} {unit}")


def _run_ts_search_study(args: argparse.Namespace, executor: str = "local", pbs_template: Optional[str] = None) -> None:
    reactant, product, ts_guess = prompt_for_ts_search_inputs()

    job_name = input("\nName for this job (used for the job folder): ").strip() or "ts_search"
    out_dir = os.path.join("jobs", job_name)
    _maybe_write_pbs_template(out_dir, executor, pbs_template)

    nprocs, mem_gb = prompt_for_resources()
    max_repair_attempts = prompt_for_max_repair_attempts()
    skip_pm6 = prompt_for_ts_search_skip_pm6()
    use_calcall = prompt_for_ts_calcall_fallback()

    elements = sorted(
        {a[0] for a in (reactant.atoms if reactant else [])}
        | {a[0] for a in (product.atoms if product else [])}
        | {a[0] for a in (ts_guess.atoms if ts_guess else [])}
    )
    print("\nNow pick the level for the final TS search stage:")
    level = prompt_for_method_basis(elements)

    want_recovery = prompt_for_ts_recovery()
    export_xyz = prompt_for_xyz_export()
    export_cdxml = prompt_for_cdxml_export()

    print("\nSearching for the TS:")
    ts_outcome = search_ts_staged(
        out_dir, method=level.method, basis=level.basis,
        basis_groups=level.basis_groups, ecp_groups=level.ecp_groups,
        nprocs=nprocs, mem_gb=mem_gb, reactant=reactant, product=product, ts_guess=ts_guess,
        skip_pm6=skip_pm6, use_calcall_fallback=use_calcall, executor=executor,
    )
    for line in ts_outcome.log:
        print(f"  {line}")

    if not ts_outcome.success:
        print("\nCouldn't find a validated TS. See the trail above.")
        return

    overlap_desc = f", overlap with reactant->product motion = {ts_outcome.match_overlap:.2f}" if ts_outcome.match_overlap is not None else ""
    print(
        f"\nOK: TS found, imaginary freq {ts_outcome.result.imaginary_freqs[0]:.1f} cm^-1{overlap_desc}"
    )
    print(f"Log file: {ts_outcome.log_path}")

    if export_xyz:
        xyz_path = convert_log_to_xyz(ts_outcome.log_path)
        print(f"XYZ file: {xyz_path}" if xyz_path else "Couldn't convert to .xyz (is Open Babel installed?).")
    if export_cdxml:
        cdxml_path = convert_log_to_cdxml(ts_outcome.log_path)
        print(f"CDXML file: {cdxml_path}" if cdxml_path else "Couldn't convert to .cdxml (is Open Babel installed?).")

    if not want_recovery:
        return

    print("\nRecovering reactant/product from the TS:")
    recovery = recover_endpoints_from_ts(
        ts_outcome, out_dir, level, nprocs=nprocs, mem_gb=mem_gb,
        reactant_ref=reactant, product_ref=product, executor=executor,
    )
    for line in recovery.log:
        print(f"  {line}")

    if isinstance(recovery, DistortionOutcome):
        if recovery.reactant_candidate:
            print(f"Recovered reactant: {recovery.reactant_candidate.outcome.log_path}")
        if recovery.product_candidate:
            print(f"Recovered product: {recovery.product_candidate.outcome.log_path}")
    else:
        if recovery.side_a:
            print(f"Recovered side A (reactant or product -- not determined which): {recovery.side_a.log_path}")
        if recovery.side_b:
            print(f"Recovered side B (reactant or product -- not determined which): {recovery.side_b.log_path}")


def _run_one_study(args: argparse.Namespace, executor: str = "local", pbs_template: Optional[str] = None) -> None:
    if args.reactant and args.product:
        study_type = "mechanism"
    elif args.reactant:
        study_type = "single"
    else:
        study_type = prompt_for_study_type()

    if study_type == "single":
        _run_single_study(args, executor=executor, pbs_template=pbs_template)
    elif study_type == "ts_search":
        _run_ts_search_study(args, executor=executor, pbs_template=pbs_template)
    else:
        _run_mechanism_study(args, executor=executor, pbs_template=pbs_template)


def main() -> None:
    args = _parse_args()
    print("gaussbot -- geometry optimization, and reaction-mechanism studies (TS, IRC, energetics)\n")

    # --executor overrides; otherwise a bare interactive invocation asks
    # once per session (not per study) -- a non-interactive/scripted
    # invocation (-r/-p given) just defaults to "local" rather than
    # prompting, since there's nothing else to ask in that mode. Same
    # split for the PBS script template: only offered interactively
    # (it can block on $EDITOR), never for scripted/automated use --
    # scripted PBS submissions fall back to job_runner.py's own
    # DEFAULT_* constants instead.
    executor = args.executor
    interactive = not (args.reactant or args.product)
    if executor is None:
        executor = prompt_for_executor() if interactive else "local"
    pbs_template = prompt_for_pbs_template() if (interactive and executor == "pbs") else None

    if args.reactant or args.product:
        _run_one_study(args, executor=executor, pbs_template=pbs_template)
        return

    while True:
        _run_one_study(args, executor=executor, pbs_template=pbs_template)
        print()
        again = questionary.confirm("Run another study?", default=False).ask()
        if not again:
            break
        print()


if __name__ == "__main__":
    main()
