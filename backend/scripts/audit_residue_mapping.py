#!/usr/bin/env python3

"""
RiceMetalAI
OsNRAMP5 residue-coordinate audit

Purpose:
    Reconcile the residue numbering used by:
        1. 8E6N reference structure
        2. global sequence alignment
        3. OsNRAMP5 target sequence
        4. AlphaFold OsNRAMP5 structure

This script does NOT perform docking.
It establishes the authoritative residue mapping that
future docking and AI modules will use.
"""

import json
from pathlib import Path

from Bio import SeqIO
from Bio.PDB import PDBParser


PROJECT_ROOT = Path(__file__).resolve().parents[2]

ALIGNMENT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_vs_8E6N_alignment.json"
)

MAPPING_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_metal_binding_mapping.json"
)

ALPHAFOLD_FILE = (
    PROJECT_ROOT
    / "data"
    / "structures"
    / "predicted"
    / "AF-Q8H4H5-F1-model_v6.pdb"
)

REFERENCE_SEQUENCE_FILE = (
    PROJECT_ROOT
    / "data"
    / "sequences"
    / "reference"
    / "8E6N_chainA.fasta"
)

TARGET_SEQUENCE_FILE = (
    PROJECT_ROOT
    / "data"
    / "sequences"
    / "proteins"
    / "OsNRAMP5_NP_001389970.1.fasta"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_residue_mapping_audit.json"
)


# Reference residues identified from the experimentally
# determined 8E6N metal environment analysis.
REFERENCE_SITES = [
    {
        "reference_position": 53,
        "reference_residue": "A",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "reference_position": 56,
        "reference_residue": "D",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "reference_position": 59,
        "reference_residue": "N",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "reference_position": 230,
        "reference_residue": "M",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "reference_position": 287,
        "reference_residue": "H",
        "metal": "MN",
        "metal_residue": 502,
    },
]


THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def load_json(path):
    """Load a JSON file."""

    if not path.exists():
        raise FileNotFoundError(f"Missing file:\n{path}")

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_sequence(path):
    """Load a FASTA sequence."""

    if not path.exists():
        raise FileNotFoundError(f"Missing FASTA file:\n{path}")

    record = SeqIO.read(path, "fasta")

    return str(record.seq).upper()


def load_alphafold_residues():
    """
    Read AlphaFold chain A and return:

        position -> {
            residue_1letter,
            residue_3letter,
            plddt
        }
    """

    if not ALPHAFOLD_FILE.exists():
        raise FileNotFoundError(
            f"Missing AlphaFold structure:\n{ALPHAFOLD_FILE}"
        )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(
        "OsNRAMP5",
        ALPHAFOLD_FILE
    )

    model = structure[0]
    chain = model["A"]

    residues = {}

    for residue in chain:

        hetflag, resseq, icode = residue.id

        # Ignore hetero residues.
        if hetflag != " ":
            continue

        residue_name = residue.resname.strip()

        if residue_name not in THREE_TO_ONE:
            continue

        one_letter = THREE_TO_ONE[residue_name]

        b_factors = [
            float(atom.get_bfactor())
            for atom in residue
        ]

        mean_plddt = (
            sum(b_factors) / len(b_factors)
            if b_factors
            else None
        )

        residues[int(resseq)] = {
            "residue_1letter": one_letter,
            "residue_3letter": residue_name,
            "mean_plddt": round(mean_plddt, 3)
            if mean_plddt is not None
            else None,
        }

    return residues


def build_alignment_map(alignment_data):
    """
    Build:

        reference_position
            ->
        target_position

    using the previously generated residue_mapping.
    """

    residue_mapping = alignment_data.get(
        "residue_mapping",
        []
    )

    mapping = {}

    for row in residue_mapping:

        target_position = row.get(
            "osnramp5_position"
        )

        reference_position = row.get(
            "reference_position"
        )

        if (
            target_position is None
            or reference_position is None
        ):
            continue

        mapping[int(reference_position)] = {
            "target_position": int(target_position),
            "target_residue": row.get(
                "osnramp5_residue"
            ),
            "reference_residue": row.get(
                "reference_residue"
            ),
        }

    return mapping


def build_mapping_file_map(mapping_data):
    """
    Load the previously generated metal-binding mapping
    as an independent comparison.
    """

    rows = mapping_data.get(
        "metal_binding_residue_mapping",
        []
    )

    result = {}

    for row in rows:

        reference_position = row.get(
            "reference_position"
        )

        if reference_position is None:
            continue

        result[int(reference_position)] = {
            "target_position": row.get(
                "osnramp5_position"
            ),
            "target_residue": row.get(
                "osnramp5_residue"
            ),
            "reference_residue": row.get(
                "reference_residue"
            ),
        }

    return result


def validate_reference_sequence(
    reference_sequence
):
    """Check reference residues against FASTA."""

    results = []

    for site in REFERENCE_SITES:

        position = site["reference_position"]
        expected = site["reference_residue"]

        if position > len(reference_sequence):

            results.append({
                "position": position,
                "expected": expected,
                "observed": None,
                "match": False,
                "status": "outside_sequence",
            })

            continue

        observed = reference_sequence[position - 1]

        results.append({
            "position": position,
            "expected": expected,
            "observed": observed,
            "match": observed == expected,
            "status": (
                "match"
                if observed == expected
                else "mismatch"
            ),
        })

    return results


def validate_target_sequence(
    target_sequence,
    position,
    expected
):
    """Validate target residue numbering."""

    if position is None:
        return {
            "position": None,
            "observed": None,
            "expected": expected,
            "match": False,
            "status": "unmapped",
        }

    if position < 1 or position > len(target_sequence):
        return {
            "position": position,
            "observed": None,
            "expected": expected,
            "match": False,
            "status": "outside_sequence",
        }

    observed = target_sequence[position - 1]

    return {
        "position": position,
        "observed": observed,
        "expected": expected,
        "match": observed == expected,
        "status": (
            "match"
            if observed == expected
            else "mismatch"
        ),
    }


def compare_site(
    site,
    alignment_map,
    mapping_file_map,
    target_sequence,
    alphafold_residues
):
    """
    Compare one reference residue through every
    coordinate system.
    """

    reference_position = site["reference_position"]
    reference_residue = site["reference_residue"]

    alignment_result = alignment_map.get(
        reference_position
    )

    mapping_result = mapping_file_map.get(
        reference_position
    )

    alignment_target_position = None

    if alignment_result:
        alignment_target_position = (
            alignment_result["target_position"]
        )

    mapping_target_position = None

    if mapping_result:
        mapping_target_position = (
            mapping_result["target_position"]
        )

    alignment_target_validation = validate_target_sequence(
        target_sequence,
        alignment_target_position,
        None,
    )

    mapping_target_validation = validate_target_sequence(
        target_sequence,
        mapping_target_position,
        None,
    )

    alphafold_position = mapping_target_position

    alphafold_result = None

    if alphafold_position is not None:

        alphafold_result = alphafold_residues.get(
            alphafold_position
        )

    # Compare the alignment-derived and mapping-file
    # target positions.
    coordinate_consistent = (
        alignment_target_position
        == mapping_target_position
    )

    # Compare AlphaFold residue with target sequence.
    alphafold_sequence_match = None

    if alphafold_result is not None:

        alphafold_sequence_match = (
            alphafold_result["residue_1letter"]
            == target_sequence[
                alphafold_position - 1
            ]
        )

    return {
        "metal": site["metal"],
        "metal_residue": site["metal_residue"],
        "reference": {
            "position": reference_position,
            "residue": reference_residue,
        },
        "alignment_mapping": {
            "target_position": alignment_target_position,
            "target_residue": (
                alignment_result["target_residue"]
                if alignment_result
                else None
            ),
            "target_sequence_check":
                alignment_target_validation,
        },
        "mapping_file": {
            "target_position": mapping_target_position,
            "target_residue": (
                mapping_result["target_residue"]
                if mapping_result
                else None
            ),
            "target_sequence_check":
                mapping_target_validation,
        },
        "alphafold": {
            "position": alphafold_position,
            "structure_residue": (
                alphafold_result["residue_1letter"]
                if alphafold_result
                else None
            ),
            "structure_residue_3letter": (
                alphafold_result["residue_3letter"]
                if alphafold_result
                else None
            ),
            "mean_plddt": (
                alphafold_result["mean_plddt"]
                if alphafold_result
                else None
            ),
            "sequence_match": alphafold_sequence_match,
        },
        "consistency": {
            "alignment_vs_mapping_file":
                coordinate_consistent,
            "mapping_to_alphafold":
                alphafold_sequence_match is True,
        },
    }


def summarize(results):
    """Create a concise audit summary."""

    total = len(results)

    alignment_mapping_consistent = sum(
        1
        for result in results
        if result["consistency"][
            "alignment_vs_mapping_file"
        ]
    )

    alphafold_matches = sum(
        1
        for result in results
        if result["consistency"][
            "mapping_to_alphafold"
        ]
    )

    return {
        "total_reference_sites": total,
        "alignment_mapping_consistent":
            alignment_mapping_consistent,
        "alphafold_sequence_matches":
            alphafold_matches,
        "all_alignment_mappings_consistent":
            alignment_mapping_consistent == total,
        "all_alphafold_positions_match_sequence":
            alphafold_matches == total,
    }


def main():

    print("=" * 70)
    print("RiceMetalAI - OsNRAMP5 Residue Mapping Audit")
    print("=" * 70)
    print()

    print("Loading alignment...")
    alignment_data = load_json(
        ALIGNMENT_FILE
    )

    print("Loading existing residue mapping...")
    mapping_data = load_json(
        MAPPING_FILE
    )

    print("Loading reference sequence...")
    reference_sequence = load_sequence(
        REFERENCE_SEQUENCE_FILE
    )

    print("Loading OsNRAMP5 sequence...")
    target_sequence = load_sequence(
        TARGET_SEQUENCE_FILE
    )

    print("Loading AlphaFold structure...")
    alphafold_residues = load_alphafold_residues()

    print()

    print(
        "Reference sequence length:",
        len(reference_sequence)
    )

    print(
        "OsNRAMP5 sequence length:",
        len(target_sequence)
    )

    print(
        "AlphaFold residues:",
        len(alphafold_residues)
    )

    print()

    print("Validating reference numbering...")

    reference_validation = (
        validate_reference_sequence(
            reference_sequence
        )
    )

    for result in reference_validation:

        status = "OK" if result["match"] else "ERROR"

        print(
            f"Reference {result['position']}: "
            f"{result['expected']} "
            f"-> "
            f"{result['observed']} "
            f"[{status}]"
        )

    print()

    alignment_map = build_alignment_map(
        alignment_data
    )

    mapping_file_map = (
        build_mapping_file_map(
            mapping_data
        )
    )

    print(
        "Building authoritative residue comparison..."
    )
    print()

    results = []

    for site in REFERENCE_SITES:

        result = compare_site(
            site=site,
            alignment_map=alignment_map,
            mapping_file_map=mapping_file_map,
            target_sequence=target_sequence,
            alphafold_residues=alphafold_residues,
        )

        results.append(result)

        reference_position = (
            result["reference"]["position"]
        )

        alignment_position = (
            result["alignment_mapping"][
                "target_position"
            ]
        )

        mapping_position = (
            result["mapping_file"][
                "target_position"
            ]
        )

        alphafold_position = (
            result["alphafold"]["position"]
        )

        print(
            f"Reference {reference_position}"
        )

        print(
            f"  Alignment target : "
            f"{alignment_position}"
        )

        print(
            f"  Mapping target   : "
            f"{mapping_position}"
        )

        print(
            f"  AlphaFold target : "
            f"{alphafold_position}"
        )

        print(
            "  Alignment/mapping consistent: "
            f"{result['consistency']['alignment_vs_mapping_file']}"
        )

        print(
            "  AlphaFold sequence match: "
            f"{result['consistency']['mapping_to_alphafold']}"
        )

        print()

    summary = summarize(results)

    print("=" * 70)
    print("AUDIT SUMMARY")
    print("=" * 70)

    print(
        "Reference sites:",
        summary["total_reference_sites"]
    )

    print(
        "Alignment/mapping consistent:",
        summary[
            "alignment_mapping_consistent"
        ],
    )

    print(
        "AlphaFold sequence matches:",
        summary[
            "alphafold_sequence_matches"
        ],
    )

    print()

    if (
        summary[
            "all_alignment_mappings_consistent"
        ]
        and
        summary[
            "all_alphafold_positions_match_sequence"
        ]
    ):
        status = "CONSISTENT"

    else:
        status = "REQUIRES_REVIEW"

    print("Overall mapping status:", status)
    print()

    result = {
        "project": "RiceMetalAI",
        "gene": "OsNRAMP5",
        "reference_structure": "8E6N",
        "alphafold_model":
            "AF-Q8H4H5-F1-model_v6",
        "purpose": (
            "Audit residue numbering and "
            "coordinate consistency before "
            "molecular docking."
        ),
        "reference_sites":
            REFERENCE_SITES,
        "reference_sequence_validation":
            reference_validation,
        "site_comparisons":
            results,
        "summary":
            summary,
        "final_status":
            status,
        "interpretation": (
            "This audit compares the previously "
            "generated alignment and residue mapping "
            "against the OsNRAMP5 sequence and "
            "AlphaFold structure. It is a coordinate "
            "quality-control step and does not "
            "establish metal binding."
        ),
    }

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with OUTPUT_FILE.open(
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            result,
            handle,
            indent=4
        )

    print(
        "Audit results saved to:"
    )

    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
