from pathlib import Path
import json

from Bio import SeqIO
from Bio.Align import PairwiseAligner


# =========================================================
# File paths
# =========================================================

OSNRAMP5_FILE = Path(
    "data/sequences/proteins/OsNRAMP5_NP_001389970.1.fasta"
)

REFERENCE_FILE = Path(
    "data/sequences/reference/8E6N_chainA.fasta"
)

OUTPUT_DIR = Path(
    "data/processed/alignments"
)

OUTPUT_FILE = OUTPUT_DIR / "OsNRAMP5_vs_8E6N_alignment.json"


# =========================================================
# Load FASTA sequence
# =========================================================

def load_sequence(path: Path):

    if not path.exists():
        raise FileNotFoundError(
            f"Sequence file not found: {path}"
        )

    record = SeqIO.read(path, "fasta")

    return {
        "id": str(record.id),
        "description": str(record.description),
        "sequence": str(record.seq),
        "length": int(len(record.seq)),
    }


# =========================================================
# Create sequence aligner
# =========================================================

def create_aligner():

    aligner = PairwiseAligner()

    # Global alignment
    aligner.mode = "global"

    # Scoring
    aligner.match_score = 2
    aligner.mismatch_score = -1

    # Gap penalties
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    return aligner


# =========================================================
# Build residue mapping
# =========================================================

def build_residue_mapping(alignment):

    target = alignment.target
    query = alignment.query

    coordinates = alignment.coordinates

    mapping = []

    for i in range(len(coordinates[0]) - 1):

        target_start = int(coordinates[0][i])
        target_end = int(coordinates[0][i + 1])

        query_start = int(coordinates[1][i])
        query_end = int(coordinates[1][i + 1])

        target_step = target_end - target_start
        query_step = query_end - query_start

        # -------------------------------------------------
        # Both sequences contain residues
        # -------------------------------------------------

        if target_step > 0 and query_step > 0:

            shared_length = min(
                target_step,
                query_step
            )

            for j in range(shared_length):

                target_index = target_start + j
                query_index = query_start + j

                target_residue = str(
                    target[target_index]
                )

                query_residue = str(
                    query[query_index]
                )

                mapping.append(
                    {
                        "osnramp5_position": int(
                            target_index + 1
                        ),
                        "osnramp5_residue": target_residue,
                        "reference_position": int(
                            query_index + 1
                        ),
                        "reference_residue": query_residue,
                        "match": bool(
                            target_residue == query_residue
                        ),
                    }
                )

        # -------------------------------------------------
        # OsNRAMP5 has residues but reference has a gap
        # -------------------------------------------------

        elif target_step > 0 and query_step == 0:

            for j in range(target_step):

                target_index = target_start + j

                mapping.append(
                    {
                        "osnramp5_position": int(
                            target_index + 1
                        ),
                        "osnramp5_residue": str(
                            target[target_index]
                        ),
                        "reference_position": None,
                        "reference_residue": None,
                        "match": False,
                    }
                )

        # -------------------------------------------------
        # Reference has residues but OsNRAMP5 has a gap
        # -------------------------------------------------

        elif target_step == 0 and query_step > 0:

            for j in range(query_step):

                query_index = query_start + j

                mapping.append(
                    {
                        "osnramp5_position": None,
                        "osnramp5_residue": None,
                        "reference_position": int(
                            query_index + 1
                        ),
                        "reference_residue": str(
                            query[query_index]
                        ),
                        "match": False,
                    }
                )

    return mapping


# =========================================================
# Calculate sequence identity
# =========================================================

def calculate_identity(mapping):

    aligned_pairs = [
        item
        for item in mapping
        if (
            item["osnramp5_position"] is not None
            and item["reference_position"] is not None
        )
    ]

    identical = sum(
        1
        for item in aligned_pairs
        if item["match"]
    )

    aligned_count = len(aligned_pairs)

    if aligned_count > 0:
        identity = (
            identical / aligned_count
        ) * 100
    else:
        identity = 0.0

    return (
        int(aligned_count),
        int(identical),
        float(identity),
    )


# =========================================================
# Main analysis
# =========================================================

def main():

    print("Loading OsNRAMP5 sequence...")

    osnramp5 = load_sequence(
        OSNRAMP5_FILE
    )

    print(
        f"OsNRAMP5: {osnramp5['id']} "
        f"({osnramp5['length']} aa)"
    )

    print()

    print("Loading 8E6N reference sequence...")

    reference = load_sequence(
        REFERENCE_FILE
    )

    print(
        f"8E6N reference: {reference['id']} "
        f"({reference['length']} aa)"
    )

    print()

    print("Performing global sequence alignment...")

    aligner = create_aligner()

    alignments = aligner.align(
        osnramp5["sequence"],
        reference["sequence"],
    )

    if len(alignments) == 0:
        raise RuntimeError(
            "No sequence alignment was produced."
        )

    alignment = alignments[0]

    print("Alignment complete.")
    print()

    # -----------------------------------------------------
    # Alignment score
    # -----------------------------------------------------

    alignment_score = float(
        alignment.score
    )

    print(
        f"Alignment score: {alignment_score}"
    )

    # -----------------------------------------------------
    # Residue mapping
    # -----------------------------------------------------

    mapping = build_residue_mapping(
        alignment
    )

    (
        aligned_count,
        identical_count,
        identity_percentage,
    ) = calculate_identity(mapping)

    print(
        f"Aligned residue pairs: "
        f"{aligned_count}"
    )

    print(
        f"Identical residues: "
        f"{identical_count}"
    )

    print(
        f"Sequence identity: "
        f"{identity_percentage:.2f}%"
    )

    # -----------------------------------------------------
    # Create output directory
    # -----------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # -----------------------------------------------------
    # JSON result
    # -----------------------------------------------------

    result = {

        "target": {
            "gene": "OsNRAMP5",
            "accession": "NP_001389970.1",
            "sequence_length_aa": int(
                osnramp5["length"]
            ),
            "source_file": str(
                OSNRAMP5_FILE
            ),
        },

        "reference": {
            "structure": "8E6N",
            "chain": "A",
            "protein": "MntH",
            "organism": (
                "Deinococcus radiodurans"
            ),
            "sequence_length_aa": int(
                reference["length"]
            ),
            "source_file": str(
                REFERENCE_FILE
            ),
        },

        "alignment": {

            "method": (
                "Biopython PairwiseAligner"
            ),

            "mode": "global",

            "match_score": 2,

            "mismatch_score": -1,

            "open_gap_score": -5,

            "extend_gap_score": -0.5,

            "score": float(
                alignment_score
            ),

            "aligned_residue_pairs": int(
                aligned_count
            ),

            "identical_residues": int(
                identical_count
            ),

            "identity_percentage": round(
                float(identity_percentage),
                3
            ),
        },

        "residue_mapping": mapping,

        "interpretation": (
            "Global sequence alignment between "
            "OsNRAMP5 and the 8E6N MntH "
            "structural reference. The residue "
            "mapping can be used to transfer "
            "structural annotations from the "
            "reference structure to corresponding "
            "OsNRAMP5 positions."
        ),
    }

    # -----------------------------------------------------
    # Save JSON
    # -----------------------------------------------------

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            result,
            f,
            indent=4
        )

    print()

    print(
        "Alignment results saved to:"
    )

    print(OUTPUT_FILE)


# =========================================================
# Entry point
# =========================================================

if __name__ == "__main__":
    main()
