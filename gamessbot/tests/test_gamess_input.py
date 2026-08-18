import os

import pytest

from gamessbot.active_space import build_active_space
from gamessbot.gamess_input import (
    build_casscf_input,
    build_cis_input,
    build_rhf_input,
    build_transitn_input,
    data_block_from_gaussian_log,
    extract_data_block,
    extract_optimized_mcscf_vec_annotation,
    extract_optimized_mcscf_vec_block,
    extract_vec_annotation,
    extract_vec_block,
    relabel_vec_group,
)

# Real files already on this machine from prior GAMESS/Gaussian work --
# used as ground truth rather than hand-written fixtures, same practice
# as the rest of this package's verification.
H2O_RHF_DAT = "/home/srr/jitendra/scratch/h2o-rhf.dat"
GAUSSBOT_LOG = "/home/srr/BOT/gaussbot/jobs/test12/product_final.log"

pytestmark = pytest.mark.skipif(
    not (os.path.exists(H2O_RHF_DAT) and os.path.exists(GAUSSBOT_LOG)),
    reason="real reference files not present on this machine",
)


def test_extract_vec_block_from_real_punch_file():
    with open(H2O_RHF_DAT) as f:
        text = f.read()
    vec_block = extract_vec_block(text)
    lines = vec_block.splitlines()
    assert lines[0].startswith(" $VEC")
    assert lines[-1].strip() == "$END"
    # every line in between is a group marker or a punched-coefficient row
    assert len(lines) > 2


def test_extract_data_block_leading_space_normalized():
    # obabel's gamin output omits the leading space GAMESS's card-image
    # reader requires on group markers -- confirmed by a real rungms
    # failure ("ERROR, NO $DATA GROUP WAS FOUND") before this was fixed.
    raw = "$CONTRL COORD=CART UNITS=ANGS $END\n$DATA\ntitle\nC1\nO 8.0 0.0 0.0 0.0\n $END\n"
    block = extract_data_block(raw, title="water")
    lines = block.splitlines()
    assert lines[0] == " $DATA"
    assert lines[1] == "water"
    assert lines[-1] == " $END"


def test_data_block_from_gaussian_log_real_obabel_conversion():
    data_block = data_block_from_gaussian_log(GAUSSBOT_LOG, title="acetaldehyde")
    lines = data_block.splitlines()
    assert lines[0] == " $DATA"
    assert lines[1] == "acetaldehyde"
    assert lines[2] == "C1"
    assert lines[-1] == " $END"
    assert len(lines) > 4  # more than just the header/footer -- real atom cards present


def test_build_rhf_input_huckel_guess():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    inp = build_rhf_input(data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3")
    assert " $CONTRL SCFTYP=RHF RUNTYP=ENERGY ICHARG=0 MULT=1 $END" in inp
    assert " $SCF DIRSCF=.TRUE. SOSCF=.TRUE. DIIS=.FALSE. $END" in inp
    assert " $GUESS GUESS=HUCKEL $END" in inp
    assert data_block in inp


def test_build_rhf_input_moread_guess_needs_vec_and_norb():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    with pytest.raises(ValueError):
        build_rhf_input(data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", guess="MOREAD")


def test_build_rhf_input_diis():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    inp = build_rhf_input(data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", use_soscf=False)
    assert " $SCF DIRSCF=.TRUE. SOSCF=.FALSE. DIIS=.TRUE. $END" in inp


def test_build_cis_input():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    inp = build_cis_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3",
        nstate=3, norb=7, vec_block=vec_block,
    )
    assert "CITYP=CIS" in inp
    assert " $GUESS GUESS=MOREAD NORB=7 $END" in inp
    assert vec_block in inp
    assert " $CIS NSTATE=3 MULT=1 $END" in inp


def test_build_rhf_input_sets_ispher_for_spherical_only_families():
    # Real rungms abort otherwise: "MODERN BASIS SET FAMILIES SUCH AS
    # CC, PCSEG, SPK, KTZ, MCP, IMCP, ZFK ARE INTENDED FOR USE ONLY AS
    # SPHERICAL HARMONIC BASIS SETS. PLEASE SET ISPHER=1 ..."
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    inp = build_rhf_input(data_block, charge=0, mult=1, gbasis_line="GBASIS=CCD")
    assert "ISPHER=1" in inp


def test_build_rhf_input_omits_ispher_for_cartesian_families():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    inp = build_rhf_input(data_block, charge=0, mult=1, gbasis_line="GBASIS=N31 NGAUSS=6 NDFUNC=1")
    assert "ISPHER" not in inp


def test_build_cis_input_sets_ispher_for_spherical_only_families():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    inp = build_cis_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=CCT",
        nstate=3, norb=7, vec_block=vec_block,
    )
    assert "ISPHER=1" in inp


def test_build_casscf_input_omits_norder_when_no_reordering_needed():
    # Real reference (a converged production CASSCF+XMCQDPT job) never
    # sets NORDER=1 unless MOs genuinely need reordering -- confirmed by
    # diffing our own generated input against it. occ/virt already
    # sitting in their target window (the common case) needs no swaps.
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    active_space = build_active_space(n_occ=5, norb=7, occ_selected=[5], virt_selected=[6])
    inp = build_casscf_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", vec_block=vec_block,
        norb=7, active_space=active_space, nstate=2,
    )
    assert "NORDER" not in inp
    assert " $GUESS  GUESS=MOREAD NORB=7 $END" in inp


def test_build_casscf_input_includes_norder_when_reordering_needed():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    # MOs 2 and 10 (far from the HOMO/LUMO boundary) force a real swap.
    active_space = build_active_space(n_occ=5, norb=10, occ_selected=[2], virt_selected=[10])
    inp = build_casscf_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", vec_block=vec_block,
        norb=10, active_space=active_space, nstate=2,
    )
    assert "NORDER=1" in inp
    assert "IORDER" in inp


def test_build_casscf_input_sets_itermx_matching_maxit():
    # Real reference input matches ITERMX to the outer MAXIT (both 120)
    # -- GAMESS's own $GUGDIA default (ITERMX=50) is fine for a tiny
    # active space's CI diagonalization but the reference's convention
    # is followed here for consistency/robustness on larger ones too.
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    active_space = build_active_space(n_occ=5, norb=7, occ_selected=[5], virt_selected=[6])
    inp = build_casscf_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", vec_block=vec_block,
        norb=7, active_space=active_space, nstate=2, maxit=200,
    )
    assert " $GUGDIA NSTATE=2 ITERMX=200 $END" in inp


# ------------------------------------------------- orbital-annotation carrying

ETHENE_RHF_DAT = "/home/srr/BOT/gamessbot/jobs/ethene-test-1/rhf_soscf.dat"
ETHENE_CASSCF_DAT = "/home/srr/BOT/gamessbot/jobs/Ethene-test-2/casscf_try1.dat"

pytestmark_annotation = pytest.mark.skipif(
    not (os.path.exists(ETHENE_RHF_DAT) and os.path.exists(ETHENE_CASSCF_DAT)),
    reason="real reference punch files not present on this machine",
)


@pytestmark_annotation
def test_extract_vec_annotation_real_rhf_punch():
    with open(ETHENE_RHF_DAT) as f:
        text = f.read()
    annotation = extract_vec_annotation(text)
    assert annotation.startswith("--- CLOSED SHELL ORBITALS ---")
    assert "E(RHF)=" in annotation


@pytestmark_annotation
def test_extract_optimized_mcscf_vec_annotation_real_casscf_punch():
    with open(ETHENE_CASSCF_DAT) as f:
        text = f.read()
    annotation = extract_optimized_mcscf_vec_annotation(text)
    assert annotation.startswith("--- OPTIMIZED MCSCF MO-S ---")
    assert "E(MCSCF)=" in annotation
    # the block extraction itself must still work unaffected
    vec_block = extract_optimized_mcscf_vec_block(text)
    assert vec_block.splitlines()[0].strip() == "$VEC"
    assert vec_block.splitlines()[-1].strip() == "$END"


def test_extract_vec_annotation_returns_empty_when_absent():
    text = " $DATA\nwater\nC1\n $END\n $VEC\n 1 1...\n $END\n"
    assert extract_vec_annotation(text) == ""


def test_relabel_vec_group_renames_opening_line_only():
    vec_block = " $VEC\n 1  1 1.0\n 1  2 2.0\n $END"
    relabeled = relabel_vec_group(vec_block, "VEC1")
    lines = relabeled.splitlines()
    assert lines[0].strip() == "$VEC1"
    assert lines[1] == " 1  1 1.0"
    assert lines[-1].strip() == "$END"


def test_annotated_vec_block_carried_into_rhf_moread_input():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    annotation = "--- CLOSED SHELL ORBITALS --- GENERATED AT some time\nwater\nE(RHF)= -1.0, E(NUC)= 1.0"
    inp = build_rhf_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", guess="MOREAD",
        vec_block=vec_block, norb=7, vec_annotation=annotation,
    )
    assert annotation in inp
    assert inp.index(annotation) < inp.index(" $VEC\n 1 1...")


def test_annotated_vec_block_carried_into_cis_input():
    data_block = " $DATA\nwater\nC1\nO 8.0 0.0 0.0 0.0\n $END"
    vec_block = " $VEC\n 1 1...\n $END"
    annotation = "--- CLOSED SHELL ORBITALS --- GENERATED AT some time"
    inp = build_cis_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3",
        nstate=3, norb=7, vec_block=vec_block, vec_annotation=annotation,
    )
    assert annotation in inp


# --------------------------------------------------------- build_transitn_input

def test_build_transitn_input_matches_real_reference_structure():
    # Confirmed against a real production TRANSITN input and the local
    # GAMESS 2024 manual's own $TRANST section (docs-input.txt):
    # NOCC = NFZC + NDOC + NVAL, IROOTS(1) = NSTATE, orbitals in $VEC1.
    data_block = " $DATA\ntitle\nC1\n $END"
    vec_block = " $VEC1\n 1  1 1.0\n $END"
    active_space = build_active_space(n_occ=101, norb=438, occ_selected=[100, 101], virt_selected=[102, 103])
    inp = build_transitn_input(
        data_block, charge=0, mult=1, gbasis_line="GBASIS=CCD", vec_block=vec_block,
        norb=438, active_space=active_space, nstate=8, mem_mwords=1000,
    )
    assert " $CONTRL SCFTYP=MCSCF CITYP=GUGA RUNTYP=TRANSITN ICHARG=0 MULT=1 ISPHER=1 $END" in inp
    assert " $DRT1   GROUP=C1 FORS=.T. NFZC=99 NDOC=2 NVAL=2 STSYM=A $END" in inp
    assert " $TRANST NFZC=99 NSTATE=8 OPERAT=DM IROOTS(1)=8 NOCC=103 $END" in inp
    assert " $GUGDIA NSTATE=8 ITERMX=120 $END" in inp
    assert "$VEC1" in inp
    assert "MAXIT" not in inp.split(" $DRT1")[0].split(" $MCSCF")[1]  # no MAXIT on $MCSCF, matching the real reference


def test_build_transitn_input_wstate_length_mismatch_raises():
    data_block = " $DATA\ntitle\nC1\n $END"
    vec_block = " $VEC1\n 1  1 1.0\n $END"
    active_space = build_active_space(n_occ=5, norb=10, occ_selected=[4, 5], virt_selected=[6, 7])
    with pytest.raises(ValueError):
        build_transitn_input(
            data_block, charge=0, mult=1, gbasis_line="GBASIS=STO NGAUSS=3", vec_block=vec_block,
            norb=10, active_space=active_space, nstate=2, wstate=[1],
        )
