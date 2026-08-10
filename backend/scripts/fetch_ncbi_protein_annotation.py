from pathlib import Path

from Bio import Entrez, SeqIO


# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------

ENTREZ_EMAIL = "bharathikoruprolu25@gmail.com"

PROTEIN_ACCESSION = "NP_001389970.1"

OUTPUT_DIR = Path(
    "data/raw/ncbi"
)

OUTPUT_FILE = OUTPUT_DIR / (
    "OsNRAMP5_NP_001389970.1_protein_annotation.gb"
)


# ---------------------------------------------------------
# Fetch NCBI protein record
# ---------------------------------------------------------

def fetch_protein_record():

    Entrez.email = ENTREZ_EMAIL

    print(
        f"Fetching NCBI protein record "
        f"{PROTEIN_ACCESSION}..."
    )

    handle = Entrez.efetch(
        db="protein",
        id=PROTEIN_ACCESSION,
        rettype="gb",
        retmode="text"
    )

    record = SeqIO.read(
        handle,
        "genbank"
    )

    handle.close()

    return record


# ---------------------------------------------------------
# Extract annotated features
# ---------------------------------------------------------

def extract_features(record):

    annotations = []

    for feature in record.features:

        if feature.type not in {
            "Region",
            "Site",
            "Transmembrane",
            "Domain",
            "Signal peptide",
            "Mat_peptide"
        }:

            continue

        annotations.append(
            {
                "type": feature.type,

                "location":
                    str(feature.location),

                "qualifiers":
                    dict(feature.qualifiers)
            }
        )

    return annotations


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    record = fetch_protein_record()

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        SeqIO.write(
            record,
            handle,
            "genbank"
        )

    print()
    print("NCBI protein record retrieved.")

    print(
        f"Protein ID: {record.id}"
    )

    print(
        f"Description: {record.description}"
    )

    print(
        f"Sequence length: {len(record.seq)} aa"
    )

    annotations = extract_features(
        record
    )

    print()
    print(
        f"Annotated functional/topology features: "
        f"{len(annotations)}"
    )

    for feature in annotations:

        print()
        print(
            f"Type: {feature['type']}"
        )

        print(
            f"Location: "
            f"{feature['location']}"
        )

        print(
            f"Qualifiers: "
            f"{feature['qualifiers']}"
        )

    print()
    print(
        "Saved NCBI protein annotation:"
    )

    print(
        OUTPUT_FILE
    )


if __name__ == "__main__":
    main()
