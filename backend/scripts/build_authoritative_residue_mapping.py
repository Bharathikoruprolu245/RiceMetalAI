
#!/usr/bin/env python3

"""
RiceMetalAI
Build authoritative residue mapping between:

8E6N PDB numbering
        ↓
8E6N chain-A sequence numbering
        ↓
OsNRAMP5 sequence
        ↓
OsNRAMP5 AlphaFold structure

This corrects for gaps between PDB residue numbering and
the continuous extracted reference sequence.

This script does NOT perform docking.
"""

import json
from pathlib import Path

from Bio import SeqIO
from Bio.PDB import PDBParser


PROJECT_ROOT = Path(__file__).resolve().parents[2]

PDB_FILE = (
    PROJECT_ROOT
    / "data"
    / "structures"
    / "pdb"
    / "8E6N.pdb"
)

REFERENCE_FASTA = (
    PROJECT_ROOT
    / "data"
    / "sequences"
    / "reference"
    / "8E6N_chainA.fasta"
)

TARGET_FASTA = (
    PROJECT_ROOT
    / "data"
    / "sequences"
    / "proteins"
    / "OsNRAMP5_NP_001389970.1.fasta"
)

ALIGNMENT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_vs_8E6N_alignment.json"
)

ALPHAFOLD_FILE = (
    PROJECT_ROOT
    / "data"
    / "structures"
    / "predicted"
    / "AF-Q8H4H5-F1-model_v6.pdb"
)

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_authoritative_residue_mapping.json"
)


# These are the residues identified in the 8E6N
# metal-proximity analysis.
REFERENCE_SITES = [
    {
        "pdb_residue_number": 53,
        "pdb_residue": "ALA",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "pdb_residue_number": 56,
        "pdb_residue": "ASP",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "pdb_residue_number": 59,
        "pdb_residue": "ASN",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "pdb_residue_number": 230,
        "pdb_residue": "MET",
        "metal": "MN",
        "metal_residue": 501,
    },
    {
        "pdb_residue_number": 287,
        "pdb_residue": "HIS",
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


def load_fasta(path):
    """Load one FASTA sequence."""

    if not path.exists():
        raise FileNotFoundError(
            f"Missing FASTA file:\n{path}"
        )

    record = SeqIO.read(path, "fasta")

    return str(record.seq).upper()


def load_alignment():
    """Load the previously generated global alignment."""

    if not ALIGNMENT_FILE.exists():
        raise FileNotFoundError(
            f"Missing alignment file:\n{ALIGNMENT_FILE}"
        )

    with ALIGNMENT_FILE.open(
        "r",
        encoding="utf-8"
    ) as handle:
        return json.load(handle)


def extract_pdb_chain_mapping():
    """
    Extract the actual PDB residue numbering and
    convert it to continuous sequence indexing.

    Returns:

        pdb_number -> {
            sequence_index,
            residue_3letter,
            residue_1letter
        }
    """

    if not PDB_FILE.exists():
        raise FileNotFoundError(
            f"Missing PDB structure:\n{PDB_FILE}"
        )

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "8E6N",
        PDB_FILE
    )

    model = structure[0]
    chain = model["A"]

    mapping = {}

    sequence_index = 0

    for residue in chain:

        hetflag, resseq, icode = residue.id

        # Ignore hetero residues.
        if hetflag != " ":
            continue

        residue_name = residue.resname.strip()

        if residue_name not in THREE_TO_ONE:
            continue

        sequence_index += 1

        mapping[int(resseq)] = {
            "pdb_residue_number": int(resseq),
            "insertion_code": str(icode),
            "sequence_index": sequence_index,
            "residue_3letter": residue_name,
            "residue_1letter":
                THREE_TO_ONE[residue_name],
        }

    return mapping


def build_reference_to_target_map(alignment_data):
    """
    The previous alignment was generated from:

        8E6N_chainA.fasta
        OsNRAMP5 FASTA

    Therefore the alignment's reference_position is
    a CONTINUOUS sequence position, not a PDB residue
    number.

    Returns:

        reference_sequence_position ->
        OsNRAMP5 position information
    """

    mapping = {}

    for row in alignment_data.get(
        "residue_mapping",
        []
    ):

        reference_position = row.get(
            "reference_position"
        )

        if reference_position is None:
            continue

        mapping[int(reference_position)] = {
            "target_position":
                row.get("osnramp5_position"),

            "target_residue":
                row.get("osnramp5_residue"),

            "reference_residue":
                row.get("reference_residue"),

            "match":
                row.get("match"),
        }

    return mapping


def extract_alphafold_residues():
    """
    Extract AlphaFold chain A residue numbering.

    AlphaFold model uses the OsNRAMP5 sequence numbering,
    so the resulting dictionary is:

        OsNRAMP5 position -> residue information
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

        if hetflag != " ":
            continue

        residue_name = residue.resname.strip()

        if residue_name not in THREE_TO_ONE:
            continue

        one_letter = THREE_TO_ONE[
            residue_name
        ]

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
            "residue_3letter": residue_name,
            "residue_1letter": one_letter,
            "mean_plddt": (
                round(mean_plddt, 3)
                if mean_plddt is not None
                else None
            ),
        }

    return residues


def build_authoritative_mapping():

    print("Loading 8E6N PDB numbering...")
    pdb_mapping = extract_pdb_chain_mapping()

    print("Loading 8E6N reference sequence...")
    reference_sequence = load_fasta(
        REFERENCE_FASTA
    )

    print("Loading OsNRAMP5 sequence...")
    target_sequence = load_fasta(
        TARGET_FASTA
    )

    print("Loading global alignment...")
    alignment_data = load_alignment()

    print("Loading AlphaFold structure...")
    alphafold_residues = (
        extract_alphafold_residues()
    )

    print()

    print(
        "8E6N PDB residues:",
        len(pdb_mapping)
    )

    print(
        "8E6N reference sequence:",
        len(reference_sequence)
    )

    print(
        "OsNRAMP5 sequence:",
        len(target_sequence)
    )

    print(
        "AlphaFold residues:",
        len(alphafold_residues)
    )

    print()

    alignment_map = (
        build_reference_to_target_map(
            alignment_data
        )
    )

    authoritative_sites = []

    print(
        "Building authoritative residue mapping..."
    )
    print()

    for site in REFERENCE_SITES:

        pdb_position = site[
            "pdb_residue_number"
        ]

        expected_pdb_residue = (
            site["pdb_residue"]
        )

        # ------------------------------------------------
        # STEP 1
        # PDB residue number → sequence index
        # ------------------------------------------------

        pdb_entry = pdb_mapping.get(
            pdb_position
        )

        if pdb_entry is None:

            print(
                f"PDB residue {pdb_position} "
                "NOT FOUND"
            )

            authoritative_sites.append({
                **site,
                "status": "pdb_residue_not_found",
            })

            continue

        sequence_index = pdb_entry[
            "sequence_index"
        ]

        observed_pdb_residue = (
            pdb_entry["residue_3letter"]
        )

        # ------------------------------------------------
        # STEP 2
        # Validate against reference FASTA
        # ------------------------------------------------

        fasta_residue = (
            reference_sequence[
                sequence_index - 1
            ]
        )

        expected_one_letter = (
            THREE_TO_ONE[
                expected_pdb_residue
            ]
        )

        reference_sequence_match = (
            fasta_residue
            == expected_one_letter
        )

        # ------------------------------------------------
        # STEP 3
        # Reference sequence index → OsNRAMP5
        # ------------------------------------------------

        alignment_entry = (
            alignment_map.get(
                sequence_index
            )
        )

        if alignment_entry is None:

            print(
                f"PDB {pdb_position} "
                f"→ sequence {sequence_index} "
                "→ TARGET UNMAPPED"
            )

            authoritative_sites.append({
                **site,
                "reference_sequence_position":
                    sequence_index,
                "reference_sequence_residue":
                    fasta_residue,
                "reference_sequence_match":
                    reference_sequence_match,
                "status":
                    "target_position_not_found",
            })

            continue

        target_position = (
            alignment_entry["target_position"]
        )

        target_residue = (
            alignment_entry["target_residue"]
        )

        # ------------------------------------------------
        # STEP 4
        # Validate target sequence
        # ------------------------------------------------

        target_sequence_residue = (
            target_sequence[
                target_position - 1
            ]
            if (
                target_position is not None
                and 1 <= target_position
                <= len(target_sequence)
            )
            else None
        )

        target_sequence_match = (
            target_sequence_residue
            == target_residue
        )

        # ------------------------------------------------
        # STEP 5
        # Validate AlphaFold position
        # ------------------------------------------------

        alphafold_entry = (
            alphafold_residues.get(
                target_position
            )
        )

        alphafold_sequence_match = False

        if alphafold_entry is not None:

            alphafold_sequence_match = (
                alphafold_entry[
                    "residue_1letter"
                ]
                == target_sequence_residue
            )

        # ------------------------------------------------
        # Final status
        # ------------------------------------------------

        fully_validated = all([
            observed_pdb_residue
            == expected_pdb_residue,

            reference_sequence_match,

            target_sequence_match,

            alphafold_entry is not None,

            alphafold_sequence_match,
        ])

        if fully_validated:
            status = "validated"
        else:
            status = "requires_review"

        result = {
            "metal": site["metal"],
            "metal_residue":
                site["metal_residue"],

            "reference": {
                "pdb_residue_number":
                    pdb_position,

                "pdb_residue":
                    observed_pdb_residue,

                "expected_pdb_residue":
                    expected_pdb_residue,

                "sequence_index":
                    sequence_index,

                "sequence_residue":
                    fasta_residue,

                "sequence_match":
                    reference_sequence_match,
            },

            "target": {
                "gene":
                    "OsNRAMP5",

                "position":
                    target_position,

                "residue":
                    target_residue,

                "sequence_residue":
                    target_sequence_residue,

                "sequence_match":
                    target_sequence_match,
            },

            "alphafold": {
                "model":
                    "AF-Q8H4H5-F1-model_v6",

                "position":
                    target_position,

                "residue":
                    (
                        alphafold_entry[
                            "residue_1letter"
                        ]
                        if alphafold_entry
                        else None
                    ),

                "residue_3letter":
                    (
                        alphafold_entry[
                            "residue_3letter"
                        ]
                        if alphafold_entry
                        else None
                    ),

                "mean_plddt":
                    (
                        alphafold_entry[
                            "mean_plddt"
                        ]
                        if alphafold_entry
                        else None
                    ),

                "sequence_match":
                    alphafold_sequence_match,
            },

            "validation": {
                "pdb_residue_valid":
                    (
                        observed_pdb_residue
                        == expected_pdb_residue
                    ),

                "reference_sequence_valid":
                    reference_sequence_match,

                "target_sequence_valid":
                    target_sequence_match,

                "alphafold_residue_valid":
                    alphafold_entry is not None,

                "alphafold_sequence_valid":
                    alphafold_sequence_match,

                "overall_status":
                    status,
            },
        }

        authoritative_sites.append(result)

        print(
            f"PDB {pdb_position} "
            f"({observed_pdb_residue}) "
            f"→ reference sequence "
            f"{sequence_index} "
            f"→ OsNRAMP5 "
            f"{target_position} "
            f"({target_residue}) "
            f"→ AlphaFold "
            f"{target_position} "
            f"[{status}]"
        )

    return {
        "project": "RiceMetalAI",

        "purpose": (
            "Authoritative residue-coordinate "
            "mapping for downstream structural "
            "analysis and molecular docking."
        ),

        "target": {
            "gene": "OsNRAMP5",
            "ncbi_accession":
                "NP_001389970.1",
            "uniprot_accession":
                "Q8H4H5",
            "sequence_length":
                len(target_sequence),
        },

        "reference": {
            "structure": "8E6N",
            "chain": "A",
            "reference_sequence_length":
                len(reference_sequence),
            "pdb_file":
                str(
                    PDB_FILE.relative_to(
                        PROJECT_ROOT
                    )
                ),
            "reference_fasta":
                str(
                    REFERENCE_FASTA.relative_to(
                        PROJECT_ROOT
                    )
                ),
        },

        "alignment": {
            "file":
                str(
                    ALIGNMENT_FILE.relative_to(
                        PROJECT_ROOT
                    )
                ),
            "method":
                alignment_data.get(
                    "alignment",
                    {}
                ).get(
                    "method"
                ),
            "mode":
                alignment_data.get(
                    "alignment",
                    {}
                ).get(
                    "mode"
                ),
        },

        "mapping_principle": {
            "step_1":
                "8E6N PDB residue number → 8E6N continuous sequence index",

            "step_2":
                "8E6N sequence index → OsNRAMP5 sequence position through global alignment",

            "step_3":
                "OsNRAMP5 sequence position → AlphaFold residue position",

            "important_note":
                "PDB residue numbers must not be treated as continuous sequence indices.",
        },

        "sites":
            authoritative_sites,

        "summary": {
            "total_sites":
                len(authoritative_sites),

            "validated_sites":
                sum(
                    1
                    for site in authoritative_sites
                    if site.get(
                        "validation",
                        {}
                    ).get(
                        "overall_status"
                    )
                    == "validated"
                ),

            "sites_requiring_review":
                sum(
                    1
                    for site in authoritative_sites
                    if site.get(
                        "validation",
                        {}
                    ).get(
                        "overall_status"
                    )
                    != "validated"
                ),
        },

        "interpretation": (
            "This mapping establishes the coordinate "
            "relationship between experimentally "
            "numbered 8E6N residues, the extracted "
            "8E6N sequence, the OsNRAMP5 sequence, "
            "and the OsNRAMP5 AlphaFold model. "
            "The mapping is a prerequisite for "
            "reliable downstream structural analysis "
            "and molecular docking."
        ),
    }


def main():

    print("=" * 70)
    print(
        "RiceMetalAI - Authoritative Residue Mapping"
    )
    print("=" * 70)
    print()

    result = build_authoritative_mapping()

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

    print()
    print("=" * 70)
    print("MAPPING SUMMARY")
    print("=" * 70)

    print(
        "Total sites:",
        result["summary"]["total_sites"]
    )

    print(
        "Validated:",
        result["summary"]["validated_sites"]
    )

    print(
        "Requires review:",
        result["summary"]["sites_requiring_review"]
    )

    print()

    print(
        "Authoritative mapping saved to:"
    )

    print(OUTPUT_FILE)

    print()

    if (
        result["summary"]["sites_requiring_review"]
        == 0
    ):
        print(
            "STATUS: ALL RESIDUE MAPPINGS VALIDATED"
        )
    else:
        print(
            "STATUS: RESIDUE MAPPING REQUIRES REVIEW"
        )


if __name__ == "__main__":
    main()
