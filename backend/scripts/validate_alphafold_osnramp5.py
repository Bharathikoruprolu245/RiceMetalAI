from pathlib import Path
from Bio.PDB import PDBParser
from Bio import SeqIO
import json
import statistics


BASE_DIR = Path(__file__).resolve().parents[2]

PDB_FILE = (
    BASE_DIR
    / "data"
    / "structures"
    / "predicted"
    / "AF-Q8H4H5-F1-model_v6.pdb"
)

PROTEIN_FILE = (
    BASE_DIR
    / "data"
    / "sequences"
    / "proteins"
    / "OsNRAMP5_NP_001389970.1.fasta"
)

OUTPUT_DIR = (
    BASE_DIR
    / "data"
    / "processed"
    / "structures"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "OsNRAMP5_AlphaFold_validation.json"
)

CANDIDATE_RESIDUES = {
    57: "A",
    60: "D",
    63: "N",
    232: "A",
    235: "M",
    337: "Q",
}


def load_structure():

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "OsNRAMP5",
        PDB_FILE
    )

    return structure


def load_target_sequence():

    record = SeqIO.read(
        PROTEIN_FILE,
        "fasta"
    )

    return str(record.seq)


def analyze_structure(structure):

    models = list(structure)

    total_atoms = 0
    chains = []
    residues = []

    confidence_values = []

    for model in models:

        for chain in model:

            chains.append(chain.id)

            for residue in chain:

                # Ignore hetero/water residues
                if residue.id[0] != " ":
                    continue

                residues.append({
                    "chain": chain.id,
                    "position": residue.id[1],
                    "residue": residue.resname
                })

                for atom in residue:

                    total_atoms += 1

                    # AlphaFold stores pLDDT in B-factor
                    confidence_values.append(
                        atom.get_bfactor()
                    )

    return {
        "model_count": len(models),
        "chains": sorted(set(chains)),
        "protein_residue_count": len(residues),
        "total_protein_atoms": total_atoms,
        "mean_plddt": (
            statistics.mean(confidence_values)
            if confidence_values
            else None
        ),
        "minimum_plddt": (
            min(confidence_values)
            if confidence_values
            else None
        ),
        "maximum_plddt": (
            max(confidence_values)
            if confidence_values
            else None
        ),
    }


def validate_sequence_length(
    structure_info,
    target_sequence
):

    structure_length = (
        structure_info["protein_residue_count"]
    )

    target_length = len(target_sequence)

    return {
        "target_sequence_length": target_length,
        "structure_residue_count": structure_length,
        "length_difference": (
            structure_length - target_length
        ),
        "length_match": (
            structure_length == target_length
        ),
    }


def analyze_candidate_residues(
    structure
):

    results = []

    model = next(structure.get_models())

    chain = model["A"]

    for position, expected_residue in (
        CANDIDATE_RESIDUES.items()
    ):

        if position not in chain:

            results.append({
                "position": position,
                "expected_residue": expected_residue,
                "status": "missing"
            })

            continue

        residue = chain[position]

        actual_three_letter = residue.resname

        three_to_one = {
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

        actual_one_letter = (
            three_to_one.get(
                actual_three_letter,
                "?"
            )
        )

        atom_confidence = [
            atom.get_bfactor()
            for atom in residue
        ]

        mean_confidence = (
            statistics.mean(atom_confidence)
            if atom_confidence
            else None
        )

        results.append({
            "position": position,
            "expected_residue": expected_residue,
            "structure_residue_3letter":
                actual_three_letter,
            "structure_residue_1letter":
                actual_one_letter,
            "residue_matches_expected": (
                actual_one_letter
                == expected_residue
            ),
            "mean_plddt": mean_confidence,
            "minimum_plddt": (
                min(atom_confidence)
                if atom_confidence
                else None
            ),
            "maximum_plddt": (
                max(atom_confidence)
                if atom_confidence
                else None
            ),
            "status": "present"
        })

    return results


def main():

    print("Loading OsNRAMP5 AlphaFold structure...")
    print(PDB_FILE)

    structure = load_structure()

    print("Structure loaded successfully.")
    print()

    target_sequence = load_target_sequence()

    structure_info = analyze_structure(
        structure
    )

    print(
        f"Models: "
        f"{structure_info['model_count']}"
    )

    print(
        f"Chains: "
        f"{', '.join(structure_info['chains'])}"
    )

    print(
        "Protein residues:",
        structure_info[
            "protein_residue_count"
        ]
    )

    print(
        "Protein atoms:",
        structure_info[
            "total_protein_atoms"
        ]
    )

    print(
        f"Mean pLDDT: "
        f"{structure_info['mean_plddt']:.2f}"
    )

    print(
        f"Minimum pLDDT: "
        f"{structure_info['minimum_plddt']:.2f}"
    )

    print(
        f"Maximum pLDDT: "
        f"{structure_info['maximum_plddt']:.2f}"
    )

    sequence_validation = (
        validate_sequence_length(
            structure_info,
            target_sequence
        )
    )

    print()
    print("Sequence validation:")
    print(
        "Target sequence length:",
        sequence_validation[
            "target_sequence_length"
        ]
    )

    print(
        "Structure residue count:",
        sequence_validation[
            "structure_residue_count"
        ]
    )

    print(
        "Length match:",
        sequence_validation[
            "length_match"
        ]
    )

    print()
    print(
        "Checking candidate residues..."
    )

    candidate_results = (
        analyze_candidate_residues(
            structure
        )
    )

    for candidate in candidate_results:

        position = candidate["position"]

        if candidate["status"] == "missing":

            print(
                f"{position}: MISSING"
            )

            continue

        residue = (
            candidate[
                "structure_residue_1letter"
            ]
        )

        plddt = candidate["mean_plddt"]

        matches = (
            candidate[
                "residue_matches_expected"
            ]
        )

        print(
            f"{residue}{position}"
            f" | expected "
            f"{CANDIDATE_RESIDUES[position]}"
            f" | pLDDT "
            f"{plddt:.2f}"
            f" | sequence match "
            f"{matches}"
        )

    result = {

        "target": {
            "gene": "OsNRAMP5",
            "uniprot": "Q8H4H5",
            "ncbi_accession":
                "NP_001389970.1",
            "sequence_file":
                str(
                    PROTEIN_FILE.relative_to(
                        BASE_DIR
                    )
                ),
        },

        "structure": {
            "model":
                "AF-Q8H4H5-F1-model_v6",
            "source":
                "AlphaFold Protein Structure Database",
            "file":
                str(
                    PDB_FILE.relative_to(
                        BASE_DIR
                    )
                ),
            **structure_info,
        },

        "sequence_validation":
            sequence_validation,

        "candidate_residues":
            candidate_results,

        "interpretation": (
            "AlphaFold predicted structure "
            "validated for chain, residue count, "
            "sequence length, and confidence. "
            "Candidate residues were checked "
            "for presence, residue identity, "
            "and local pLDDT confidence. "
            "Structural confidence does not "
            "by itself establish metal-binding "
            "function."
        ),
    }

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            result,
            handle,
            indent=4
        )

    print()
    print(
        "AlphaFold structure validation complete."
    )

    print()
    print("Results saved to:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
