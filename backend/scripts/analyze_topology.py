from pathlib import Path
import json

from Bio import SeqIO


# ---------------------------------------------------------
# Input
# ---------------------------------------------------------

PROTEIN_FILE = Path(
    "data/sequences/proteins/"
    "OsNRAMP5_NP_001389970.1.fasta"
)

OUTPUT_DIR = Path(
    "data/processed/proteins"
)

OUTPUT_FILE = OUTPUT_DIR / (
    "OsNRAMP5_NP_001389970.1_topology.json"
)


# ---------------------------------------------------------
# Hydrophobic amino acids
# ---------------------------------------------------------

HYDROPHOBIC = set(
    "AVILMFWY"
)


# ---------------------------------------------------------
# Find hydrophobic windows
# ---------------------------------------------------------

def find_hydrophobic_windows(
    sequence,
    window_size=19,
    threshold=0.68
):

    candidates = []

    for start in range(
        0,
        len(sequence) - window_size + 1
    ):

        window = sequence[
            start:start + window_size
        ]

        hydrophobic_count = sum(
            aa in HYDROPHOBIC
            for aa in window
        )

        fraction = (
            hydrophobic_count
            / window_size
        )

        if fraction >= threshold:

            candidates.append(
                {
                    "start": start + 1,
                    "end": start + window_size,
                    "sequence": window,
                    "hydrophobic_fraction":
                        round(fraction, 3)
                }
            )

    return candidates


# ---------------------------------------------------------
# Merge overlapping windows
# ---------------------------------------------------------

def merge_candidates(candidates):

    if not candidates:
        return []

    merged = []

    current = {
        "start": candidates[0]["start"],
        "end": candidates[0]["end"],
        "max_hydrophobic_fraction":
            candidates[0]["hydrophobic_fraction"]
    }

    for candidate in candidates[1:]:

        if candidate["start"] <= current["end"] + 1:

            current["end"] = max(
                current["end"],
                candidate["end"]
            )

            current[
                "max_hydrophobic_fraction"
            ] = max(
                current[
                    "max_hydrophobic_fraction"
                ],
                candidate[
                    "hydrophobic_fraction"
                ]
            )

        else:

            merged.append(current)

            current = {
                "start": candidate["start"],
                "end": candidate["end"],
                "max_hydrophobic_fraction":
                    candidate[
                        "hydrophobic_fraction"
                    ]
            }

    merged.append(current)

    return merged


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    print(
        "Loading OsNRAMP5 protein..."
    )

    record = SeqIO.read(
        PROTEIN_FILE,
        "fasta"
    )

    sequence = str(record.seq)

    print(
        f"Protein: {record.id}"
    )

    print(
        f"Length: {len(sequence)} aa"
    )

    candidates = find_hydrophobic_windows(
        sequence
    )

    regions = merge_candidates(
        candidates
    )

    results = {

        "protein": {
            "accession":
                "NP_001389970.1",

            "gene":
                "OsNRAMP5",

            "sequence_length_aa":
                len(sequence)
        },

        "method": {
            "analysis":
                "hydrophobic-window screening",

            "window_size":
                19,

            "hydrophobic_fraction_threshold":
                0.68,

            "hydrophobic_residues":
                "AVILMFWY"
        },

        "candidate_hydrophobic_regions":
            regions,

        "interpretation":
            "Candidate hydrophobic regions "
            "identified computationally. "
            "These are not experimentally "
            "confirmed transmembrane helices."
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
            results,
            handle,
            indent=4
        )

    print()
    print(
        f"Candidate hydrophobic regions: "
        f"{len(regions)}"
    )

    for index, region in enumerate(
        regions,
        start=1
    ):

        print(
            f"Region {index}: "
            f"{region['start']}-"
            f"{region['end']} aa "
            f"(max hydrophobic fraction: "
            f"{region['max_hydrophobic_fraction']})"
        )

    print()
    print(
        "Topology screening complete."
    )

    print(
        "Results saved to:"
    )

    print(
        OUTPUT_FILE
    )


if __name__ == "__main__":
    main()
