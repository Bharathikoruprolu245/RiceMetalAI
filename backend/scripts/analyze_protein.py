from pathlib import Path
import json
from collections import Counter

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ---------------------------------------------------------
# Input / output files
# ---------------------------------------------------------

PROTEIN_FILE = Path(
    "data/sequences/proteins/"
    "OsNRAMP5_NP_001389970.1.fasta"
)

OUTPUT_DIR = Path(
    "data/processed/proteins"
)

OUTPUT_FILE = OUTPUT_DIR / (
    "OsNRAMP5_NP_001389970.1_properties.json"
)


# ---------------------------------------------------------
# Load protein sequence
# ---------------------------------------------------------

def load_protein():

    record = SeqIO.read(
        PROTEIN_FILE,
        "fasta"
    )

    sequence = str(record.seq)

    return record, sequence


# ---------------------------------------------------------
# Calculate protein properties
# ---------------------------------------------------------

def calculate_properties(sequence):

    analysis = ProteinAnalysis(sequence)

    amino_acid_counts = Counter(sequence)

    properties = {

        "sequence_length_aa": len(sequence),

        "molecular_weight_da":
            round(
                analysis.molecular_weight(),
                2
            ),

        "theoretical_pI":
            round(
                analysis.isoelectric_point(),
                3
            ),

        "aromaticity":
            round(
                analysis.aromaticity(),
                4
            ),

        "instability_index":
            round(
                analysis.instability_index(),
                3
            ),

        "gravy":
            round(
                analysis.gravy(),
                4
            ),

        "amino_acid_counts":
            dict(
                sorted(
                    amino_acid_counts.items()
                )
            ),

        "amino_acid_percentage":
            {
                aa: round(
                    percentage,
                    3
                )
                for aa, percentage
                in analysis.amino_acids_percent.items()
            }
    }

    return properties


# ---------------------------------------------------------
# Save results
# ---------------------------------------------------------

def save_results(results):

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


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    print(
        "Loading OsNRAMP5 protein..."
    )

    record, sequence = load_protein()

    print(
        f"Protein ID: {record.id}"
    )

    print(
        f"Protein length: {len(sequence)} aa"
    )

    properties = calculate_properties(
        sequence
    )

    results = {

        "protein": {
            "accession":
                "NP_001389970.1",

            "gene":
                "OsNRAMP5",

            "ncbi_gene_id":
                "4342859",

            "organism":
                "Oryza sativa Japonica Group",

            "source_file":
                str(PROTEIN_FILE)
        },

        "properties": properties
    }

    save_results(results)

    print()
    print(
        "Protein characterization complete."
    )

    print()
    print(
        f"Length: "
        f"{properties['sequence_length_aa']} aa"
    )

    print(
        f"Molecular weight: "
        f"{properties['molecular_weight_da']} Da"
    )

    print(
        f"Theoretical pI: "
        f"{properties['theoretical_pI']}"
    )

    print(
        f"Aromaticity: "
        f"{properties['aromaticity']}"
    )

    print(
        f"Instability index: "
        f"{properties['instability_index']}"
    )

    print(
        f"GRAVY: "
        f"{properties['gravy']}"
    )

    print()
    print(
        "Results saved to:"
    )

    print(
        OUTPUT_FILE
    )


if __name__ == "__main__":
    main()
