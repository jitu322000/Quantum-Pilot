"""
gamess_input.py

Turns a geometry source into a GAMESS $DATA block, and assembles full
GAMESS .inp files for the RHF, CIS, and CASSCF stages. Nothing here
talks to GAMESS itself -- see local_runner.py for that.
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from typing import List, Optional

from .active_space import ActiveSpaceSuggestion

_GROUP_RE_TEMPLATE = r"\${name}\b.*?\$END\b"

# Modern basis-function families GAMESS only defines as spherical
# harmonics -- confirmed by a real rungms abort: "MODERN BASIS SET
# FAMILIES SUCH AS CC, PCSEG, SPK, KTZ, MCP, IMCP, ZFK ARE INTENDED FOR
# USE ONLY AS SPHERICAL HARMONIC BASIS SETS. PLEASE SET ISPHER=1 IN THE
# $CONTRL GROUP." Older families (STO/N21/N31/N311/DZV/TZV) default to
# GAMESS's own cartesian ISPHER=-1 and are left alone.
_SPHERICAL_ONLY_GBASIS_PREFIXES = ("CC", "PCSEG", "SPK", "KTZ", "MCP", "IMCP", "ZFK")


def _needs_ispher(gbasis_line: str) -> bool:
    """Whether `gbasis_line` (the text between "$BASIS" and "$END") names
    a basis family that requires ISPHER=1. Must be applied identically to
    every $CONTRL built for the same job -- RHF's converged orbitals are
    read back into CIS/CASSCF/XMCQDPT via $VEC MOREAD, and a mismatched
    ISPHER setting changes the number of basis functions for any shell
    with d angular momentum or higher, corrupting that read."""
    m = re.search(r"GBASIS=(\S+)", gbasis_line, re.IGNORECASE)
    if not m:
        return False
    value = m.group(1).upper()
    return any(value.startswith(prefix) for prefix in _SPHERICAL_ONLY_GBASIS_PREFIXES)


def _normalize_group_lines(block: str) -> str:
    """GAMESS's card-image reader requires group markers ($DATA, $END,
    etc.) not to start in column 1 -- ensure every line beginning with
    '$' has a leading space, regardless of how the source file (e.g.
    obabel's gamin output, which omits it) formatted it. Confirmed by a
    real rungms failure ("ERROR, NO $DATA GROUP WAS FOUND") caused by a
    column-1 "$DATA" line coming straight out of the extraction regex."""
    return "\n".join(
        " " + line if line.startswith("$") else line
        for line in block.splitlines()
    )


def _extract_group(text: str, name: str) -> str:
    """Pull one `$NAME ... $END` group (inclusive) out of a GAMESS .inp
    or .dat file -- non-greedy, so it stops at the group's own $END
    rather than some later one in the same file."""
    m = re.search(_GROUP_RE_TEMPLATE.format(name=name), text, re.DOTALL | re.IGNORECASE)
    if m is None:
        raise ValueError(f"No ${name} group found")
    return _normalize_group_lines(m.group(0).strip())


def extract_data_block(text: str, title: Optional[str] = None) -> str:
    """
    Pull the $DATA ... $END block out of a GAMESS .inp or .dat file.
    If `title` is given, replaces the title card (the block's second
    line) -- useful since obabel's gamin output uses the source file's
    own path as the title, which isn't a great label to keep.
    """
    block = _extract_group(text, "DATA")
    if title is not None:
        lines = block.splitlines()
        if len(lines) < 2:
            raise ValueError(f"$DATA block is too short to have a title line: {block!r}")
        lines[1] = title
        block = "\n".join(lines)
    return block


def extract_vec_block(text: str) -> str:
    """Pull the $VEC ... $END block (the punched orbitals) out of a
    GAMESS .dat (PUNCH) file -- ignores the "--- CLOSED SHELL ORBITALS
    --- GENERATED AT ..." comment line GAMESS writes just above $VEC,
    since that line isn't part of the group itself (confirmed against a
    real punch file). See extract_vec_annotation() to grab that comment
    line separately, for carrying into a generated input as a plain
    comment above the orbital block."""
    return _extract_group(text, "VEC")


def extract_vec_annotation(text: str, group_name: str = "VEC") -> str:
    """
    Grabs the human-readable heading GAMESS prints just above a
    $VEC/$VEC1 punch block -- e.g. "--- CLOSED SHELL ORBITALS ---" for
    a plain RHF punch, or "--- OPTIMIZED MCSCF MO-S --- GENERATED AT
    ..." plus the filename and "E(MCSCF)=.../E(NUC)=..." lines that
    follow it for a CASSCF punch (confirmed against a real punch file
    and your own real TRANSITN reference input) -- purely informational
    text GAMESS's own card reader skips over between groups, carried
    into generated inputs so a reader can see at a glance which
    orbitals were copied in, per your request.

    Returns "" if no such heading is found immediately above the group
    (e.g. a hand-written $VEC with nothing before it) -- callers should
    treat that as "no annotation available", not an error.
    """
    m = re.search(_GROUP_RE_TEMPLATE.format(name=group_name), text, re.DOTALL | re.IGNORECASE)
    if m is None:
        return ""
    # the group marker is often indented (e.g. " $VEC"), so back up to
    # the start of ITS OWN line first -- text[:m.start()] alone would
    # otherwise end mid-line, on just that leading whitespace, and
    # never reach the real annotation lines above it
    line_start = text.rfind("\n", 0, m.start()) + 1
    prefix_lines = text[:line_start].splitlines()
    collected: List[str] = []
    for line in reversed(prefix_lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("$"):
            return ""  # walked off the top of the annotation without finding a "---" heading
        collected.append(line)
        if stripped.startswith("---"):
            break
    else:
        return ""
    return "\n".join(reversed(collected))


def extract_optimized_mcscf_vec_block(text: str) -> str:
    """
    Pull the $VEC ... $END block specifically from under the
    "OPTIMIZED MCSCF MO-S" heading in a CASSCF punch (.dat) file. A
    CASSCF punch has TWO $VEC blocks in sequence -- natural orbitals
    first, then optimized MCSCF MOs (confirmed against a real CASSCF
    .dat file) -- and XMCQDPT restarts from the latter, already in the
    correct active-space order with no further IORDER reordering needed
    (confirmed against your own real MCQDPT2 input's own comment to
    that effect).
    """
    idx = text.find("OPTIMIZED MCSCF MO-S")
    if idx == -1:
        raise ValueError('No "OPTIMIZED MCSCF MO-S" section found in this punch file')
    line_start = text.rfind("\n", 0, idx) + 1
    return extract_vec_block(text[line_start:])


def extract_optimized_mcscf_vec_annotation(text: str, group_name: str = "VEC") -> str:
    """Like extract_optimized_mcscf_vec_block(), but grabs the heading
    above that specific $VEC block (see extract_vec_annotation()) --
    used to carry the "--- OPTIMIZED MCSCF MO-S --- GENERATED AT ..."
    heading into CASSCF-restart inputs (XMCQDPT, TRANSITN)."""
    idx = text.find("OPTIMIZED MCSCF MO-S")
    if idx == -1:
        return ""
    line_start = text.rfind("\n", 0, idx) + 1
    return extract_vec_annotation(text[line_start:], group_name)


def relabel_vec_group(vec_block: str, new_name: str) -> str:
    """Renames a "$VEC...$END" block's opening line to $<new_name>
    (e.g. "$VEC1") -- GAMESS's own punch files always write plain
    "$VEC" for the OPTIMIZED MCSCF MO-S orbitals, but RUNTYP=TRANSITN
    specifically requires them presented as $VEC1 (confirmed against
    the manual's own $TRANST documentation: "Note that $GUESS is not
    read by this RUNTYP! Orbitals must be in $VEC1 and possibly $VEC2
    input groups.")."""
    lines = vec_block.splitlines()
    lines[0] = re.sub(r"\$VEC\b", f"${new_name}", lines[0], count=1, flags=re.IGNORECASE)
    return "\n".join(lines)


def data_block_from_gaussian_log(log_path: str, title: Optional[str] = None, timeout: int = 30) -> str:
    """
    Convert an already-optimized Gaussian .log to a GAMESS $DATA block
    via Open Babel (`obabel -ig09 <log> -ogamin -O <tmp>.inp`), then
    extract just the geometry ($DATA ... $END) -- obabel's own $CONTRL
    line is discarded, since the caller builds its own from scratch.

    Confirmed by hand: obabel's gamin output is atoms-only (COORD=CART,
    no inline basis), which pairs correctly with a separately-built
    $BASIS group rather than needing one inline.
    """
    with tempfile.NamedTemporaryFile(suffix=".inp", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        try:
            result = subprocess.run(
                ["obabel", "-ig09", log_path, "-ogamin", "-O", tmp_path],
                capture_output=True, text=True, timeout=timeout,
            )
        except FileNotFoundError as e:
            raise ValueError(
                "Open Babel ('obabel') isn't on PATH -- needed to convert a Gaussian log to GAMESS format."
            ) from e
        if result.returncode != 0:
            raise ValueError(f"Open Babel couldn't convert {log_path}: {result.stderr.strip()}")
        with open(tmp_path) as f:
            text = f.read()
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)

    return extract_data_block(text, title=title or os.path.splitext(os.path.basename(log_path))[0])


def data_block_from_gamess_inp(inp_path: str, title: Optional[str] = None) -> str:
    """Extract the $DATA ... $END block straight from an existing
    GAMESS input (or punch/.dat) file someone already has for this
    molecule -- no conversion needed, it's already in GAMESS format."""
    with open(inp_path) as f:
        text = f.read()
    return extract_data_block(text, title=title)


# ------------------------------------------------------------- input builders

def _annotated_vec_block(vec_block: str, annotation: str = "") -> str:
    """Prepends a human-readable heading (extract_vec_annotation()'s
    output, e.g. "--- OPTIMIZED MCSCF MO-S --- GENERATED AT ...") just
    above a $VEC/$VEC1 block in a generated input, purely as a comment
    -- GAMESS's card reader skips non-"$" lines between groups -- so a
    reader can see at a glance which orbitals were copied in, per your
    request. No-op if there's no annotation to add."""
    return f"{annotation}\n{vec_block}" if annotation else vec_block


def build_rhf_input(
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    mem_mwords: int = 1,
    use_soscf: bool = True,
    guess: str = "HUCKEL",
    vec_block: Optional[str] = None,
    norb: Optional[int] = None,
    vec_annotation: str = "",
) -> str:
    """
    Assemble a full GAMESS RHF input: $CONTRL/$SYSTEM/$BASIS/$DATA/$SCF/
    $GUESS(+ $VEC). `guess` is "HUCKEL" (a fresh start) or "MOREAD" (read
    `vec_block`, needs `norb`) -- see rhf.run_rhf_staged() for the
    SOSCF/DIIS/SOSCF(MOREAD) ladder this feeds. `gbasis_line` is the
    exact text that goes between "$BASIS" and "$END", e.g.
    "GBASIS=STO NGAUSS=3" or "GBASIS=CCT" -- see level_select.py for the
    curated GBASIS choices. `vec_annotation` (see extract_vec_annotation())
    is carried in as a plain comment just above the $VEC block, if given.
    """
    if guess == "MOREAD" and (vec_block is None or norb is None):
        raise ValueError("guess='MOREAD' needs both vec_block and norb")

    ispher = " ISPHER=1" if _needs_ispher(gbasis_line) else ""
    lines = [
        f" $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG={charge} MULT={mult}{ispher} $END",
        f" $SYSTEM MWORDS={mem_mwords} $END",
        f" $BASIS {gbasis_line} $END",
    ]
    if use_soscf:
        lines.append(" $SCF DIRSCF=.TRUE. SOSCF=.TRUE. DIIS=.FALSE. $END")
    else:
        lines.append(" $SCF DIRSCF=.TRUE. SOSCF=.FALSE. DIIS=.TRUE. $END")

    # $DATA second-to-last, $VEC (when reading orbitals) last -- keeps
    # the run's actual settings up top and the geometry/orbitals (the
    # bulky, less-often-read parts) out of the way at the bottom, for a
    # quick glance at what a calculation actually is.
    if guess == "MOREAD":
        lines.append(f" $GUESS GUESS=MOREAD NORB={norb} $END")
        lines.append(data_block)
        lines.append(_annotated_vec_block(vec_block, vec_annotation))
    else:
        lines.append(" $GUESS GUESS=HUCKEL $END")
        lines.append(data_block)

    return "\n".join(lines) + "\n"


def build_cis_input(
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    nstate: int,
    norb: int,
    vec_block: str,
    mem_mwords: int = 1,
    vec_annotation: str = "",
) -> str:
    """
    Assemble a full GAMESS CIS input, built on top of an already-
    converged RHF: reads `vec_block` (GUESS=MOREAD NORB=<norb>) instead
    of generating a fresh guess, and requests `nstate` singlet (MULT=1)
    excited states -- singlet only this round, since the RHF reference
    is closed-shell and UHF/ROHF-based triplet CIS is future work. See
    cis.run_cis(). `vec_annotation` (see extract_vec_annotation()) is
    carried in as a plain comment just above the $VEC block, if given.
    """
    ispher = " ISPHER=1" if _needs_ispher(gbasis_line) else ""
    lines = [
        f" $CONTRL SCFTYP=RHF RUNTYP=ENERGY CITYP=CIS ICHARG={charge} MULT={mult}{ispher} $END",
        f" $SYSTEM MWORDS={mem_mwords} $END",
        f" $BASIS {gbasis_line} $END",
        f" $GUESS GUESS=MOREAD NORB={norb} $END",
        f" $CIS NSTATE={nstate} MULT=1 $END",
        data_block,
        _annotated_vec_block(vec_block, vec_annotation),
    ]
    return "\n".join(lines) + "\n"


def _format_iorder(iorder: List[int]) -> str:
    """
    Renders a NORB-length permutation as individual GAMESS $GUESS
    IORDER statements, one per orbital that actually needs to move.
    active_space._build_iorder() only ever produces symmetric pairwise
    swaps (per your preference -- simpler and more conservative than a
    longer reordering cycle), so this naturally comes out as matched
    "IORDER(a)=b" / "IORDER(b)=a" pairs, e.g. "IORDER(10)=9" and
    "IORDER(9)=10" -- untouched (identity) positions are omitted
    entirely, since that's GAMESS's own IORDER default.
    """
    lines = [
        f"         IORDER({position})={value}"
        for position, value in enumerate(iorder, start=1)
        if value != position
    ]
    return "\n".join(lines)


def build_casscf_input(
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    vec_block: str,
    norb: int,
    active_space: ActiveSpaceSuggestion,
    nstate: int,
    wstate: Optional[List[float]] = None,
    mem_mwords: int = 1,
    maxit: int = 120,
    vec_annotation: str = "",
) -> str:
    """
    Assembles a full GAMESS FORS-CASSCF input on top of the RHF-
    converged orbitals: $GUESS GUESS=MOREAD + NORDER=1/IORDER (bringing
    `active_space`'s selected -- often scattered -- MOs into the single
    contiguous window $DRT's NMCC/NDOC/NVAL scheme requires), $MCSCF
    CISTEP=GUGA FULLNR=.T. FORS=.T., $DRT with the suggested/adjusted
    NMCC/NDOC/NVAL, and $GUGDIA/$GUGDM2 for `nstate` state-averaged
    roots. Format confirmed against theochem.cc's own CASSCF tutorial
    input and a real CASSCF input from this user's own prior work
    (including the IORDER reordering, for a case where the CIS-selected
    orbitals weren't already contiguous).

    `wstate` defaults to equal weights across all `nstate` states
    (state-averaged CASSCF); pass explicit weights to favor specific
    states instead.
    """
    weights = wstate if wstate is not None else [1] * nstate
    if len(weights) != nstate:
        raise ValueError(f"wstate must have exactly {nstate} entries, got {len(weights)}")
    if len(active_space.iorder) != norb:
        raise ValueError(
            f"active_space.iorder has {len(active_space.iorder)} entries, expected {norb} (norb)"
        )

    ispher = " ISPHER=1" if _needs_ispher(gbasis_line) else ""
    iorder_lines = _format_iorder(active_space.iorder)
    # NORDER=1 only when there's an actual reordering to apply -- an
    # empty IORDER block still parsed fine in testing, but the real
    # reference input this was checked against (a converged, real
    # production CASSCF+XMCQDPT run) never sets NORDER=1 unless MOs
    # genuinely needed reordering, so this matches that convention
    # instead of leaving a no-op NORDER=1 (and a stray blank line) in
    # every generated input.
    guess_line = (
        f" $GUESS  GUESS=MOREAD NORB={norb} NORDER=1\n{iorder_lines}\n $END"
        if iorder_lines
        else f" $GUESS  GUESS=MOREAD NORB={norb} $END"
    )
    lines = [
        f" $CONTRL RUNTYP=ENERGY SCFTYP=MCSCF ICHARG={charge} MULT={mult}{ispher} $END",
        f" $SYSTEM MWORDS={mem_mwords} $END",
        f" $BASIS {gbasis_line} $END",
        guess_line,
        f" $MCSCF  CISTEP=GUGA FULLNR=.T. MAXIT={maxit} FORS=.T. $END",
        f" $DRT    GROUP=C1 FORS=.T. NMCC={active_space.nmcc} NDOC={active_space.ndoc} "
        f"NVAL={active_space.nval} $END",
        f" $GUGDIA NSTATE={nstate} ITERMX={maxit} $END",
        " $GUGDM2 WSTATE(1)=" + ",".join(str(w) for w in weights) + " $END",
        data_block,
        _annotated_vec_block(vec_block, vec_annotation),
    ]
    return "\n".join(lines) + "\n"


def build_xmcqdpt_input(
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    vec_block: str,
    norb: int,
    active_space: ActiveSpaceSuggestion,
    nstate: int,
    wstate: Optional[List[float]] = None,
    mem_mwords: int = 1,
    maxit: int = 120,
    edshft: float = 0.04,
    xzero: bool = True,
    vec_annotation: str = "",
) -> str:
    """
    Assembles a full GAMESS XMCQDPT (MCQDPT2) input on top of a
    converged CASSCF: the same $DATA/$BASIS/$MCSCF/$DRT/$GUGDIA/$GUGDM2
    group shell as build_casscf_input() (same active space, re-running
    the CASSCF reference in this job too, per GAMESS's own design) plus
    $CONTRL MPLEVL=2, $MCSCF FINCI=MOS, and the two XMCQDPT-specific
    groups $MRMP MRPT=MCQDPT and $MCQDPT KSTATE(...)/XZERO/EDSHFT.

    `vec_block` should be the CASSCF run's own "OPTIMIZED MCSCF MO-S"
    block (gamess_input.extract_optimized_mcscf_vec_block()) -- already
    in the active space's final order, so unlike build_casscf_input()
    this does NOT set NORDER/IORDER (confirmed against theochem.cc's own
    XMCQDPT tutorial input and a real XMCQDPT input from this user's own
    prior work, neither of which reorder here).

    `edshft` (a small energy shift, Hartree) and `xzero` (the extended
    multistate formulation) default to the values both of those real
    references use -- EDSHFT=0.04 is also the manual's own suggested
    range for damping intruder-state singularities.
    """
    weights = wstate if wstate is not None else [1] * nstate
    if len(weights) != nstate:
        raise ValueError(f"wstate must have exactly {nstate} entries, got {len(weights)}")

    kstate = ",".join(["1"] * nstate)  # correct every state by default
    xzero_str = ".T." if xzero else ".F."

    ispher = " ISPHER=1" if _needs_ispher(gbasis_line) else ""
    lines = [
        f" $CONTRL SCFTYP=MCSCF RUNTYP=ENERGY ICHARG={charge} MULT={mult} MPLEVL=2{ispher} $END",
        f" $SYSTEM MWORDS={mem_mwords} $END",
        f" $BASIS {gbasis_line} $END",
        f" $GUESS  GUESS=MOREAD NORB={norb} $END",
        f" $MCSCF  CISTEP=GUGA FULLNR=.T. MAXIT={maxit} FORS=.T. FINCI=MOS $END",
        f" $DRT    GROUP=C1 FORS=.T. NMCC={active_space.nmcc} NDOC={active_space.ndoc} "
        f"NVAL={active_space.nval} $END",
        f" $GUGDIA NSTATE={nstate} ITERMX={maxit} $END",
        " $GUGDM2 WSTATE(1)=" + ",".join(str(w) for w in weights) + " $END",
        " $MRMP   MRPT=MCQDPT $END",
        f" $MCQDPT KSTATE(1)={kstate} XZERO={xzero_str} EDSHFT={edshft} $END",
        data_block,
        _annotated_vec_block(vec_block, vec_annotation),
    ]
    return "\n".join(lines) + "\n"


def build_transitn_input(
    data_block: str,
    charge: int,
    mult: int,
    gbasis_line: str,
    vec_block: str,
    norb: int,
    active_space: ActiveSpaceSuggestion,
    nstate: int,
    wstate: Optional[List[float]] = None,
    mem_mwords: int = 1,
    itermx: int = 120,
    vec_annotation: str = "",
) -> str:
    """
    Assembles a full GAMESS RUNTYP=TRANSITN input on top of a converged
    CASSCF: computes radiative transition moments and oscillator
    strengths between the CASSCF states (OPERAT=DM), via GAMESS's
    $TRANST group. Confirmed against the local GAMESS 2024 manual's own
    $TRANST documentation (docs-input.txt) and a real production
    TRANSITN input from this user's own prior work:

      - $DRT1 (not plain $DRT), with NFZC playing the same role as
        $DRT's NMCC.
      - $TRANST NFZC=<nmcc> NSTATE=<nstate> OPERAT=DM IROOTS(1)=<nstate>
        NOCC=<nmcc+ndoc+nval> -- NOCC is confirmed by the manual to be
        NFZC+NDOC+NVAL for a plain CAS-CI (not CI-SD/FOCI/SOCI);
        IROOTS(1) is how many of the NSTATE roots to analyze for
        transition moments -- set equal to nstate here, to get every
        requested state's oscillator strength relative to the ground
        state.
      - No MAXIT on $MCSCF here (confirmed absent in the real reference
        input too) -- this restarts from already-converged CASSCF
        orbitals purely to analyze transitions, not to re-optimize them.
      - Orbitals must be in $VEC1, not plain $VEC (the manual: "$GUESS
        is not read by this RUNTYP! Orbitals must be in $VEC1 ...") --
        `vec_block` should already be labeled $VEC1 (see
        gamess_input.relabel_vec_group()) before it's passed in here;
        this function does not relabel it itself, so it also works
        with an already-$VEC1-labeled block from elsewhere.

    `vec_block` should otherwise be the CASSCF run's own "OPTIMIZED
    MCSCF MO-S" block, same source as build_xmcqdpt_input().
    """
    weights = wstate if wstate is not None else [1] * nstate
    if len(weights) != nstate:
        raise ValueError(f"wstate must have exactly {nstate} entries, got {len(weights)}")

    nocc = active_space.nmcc + active_space.ndoc + active_space.nval
    ispher = " ISPHER=1" if _needs_ispher(gbasis_line) else ""
    lines = [
        f" $CONTRL SCFTYP=MCSCF CITYP=GUGA RUNTYP=TRANSITN ICHARG={charge} MULT={mult}{ispher} $END",
        f" $SYSTEM MWORDS={mem_mwords} $END",
        f" $BASIS {gbasis_line} $END",
        f" $GUESS  GUESS=MOREAD NORB={norb} $END",
        " $MCSCF  CISTEP=GUGA FULLNR=.T. FORS=.T. FINCI=MOS $END",
        f" $DRT1   GROUP=C1 FORS=.T. NFZC={active_space.nmcc} NDOC={active_space.ndoc} "
        f"NVAL={active_space.nval} STSYM=A $END",
        f" $TRANST NFZC={active_space.nmcc} NSTATE={nstate} OPERAT=DM IROOTS(1)={nstate} NOCC={nocc} $END",
        f" $GUGDIA NSTATE={nstate} ITERMX={itermx} $END",
        " $GUGDM2 WSTATE(1)=" + ",".join(str(w) for w in weights) + " $END",
        data_block,
        _annotated_vec_block(vec_block, vec_annotation),
    ]
    return "\n".join(lines) + "\n"
