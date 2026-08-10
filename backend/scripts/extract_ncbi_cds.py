from pathlib import Path

from Bio import Entrez, SeqIO


# ---------------------------------------------------------
# NCBI configuration
# ---------------------------------------------------------

Entrez.email = "bharathikoruprolu25@gmail.com"

MRNA_ACCESSION = "NM_001403041.1"

RAW_DIR = Path("data/raw/ncbi")
CDS_DIR = Path("data/sequences/cds")

GENBANK_FILE = (
    RAW_DIR / f"OsNRAMP5_{MRNA_ACCESSION}_genbank.gb"
)

CDS_FILE = (
    CDS_DIR / f"OsNRAMP5_{MRNA_ACCESSION}_CDS.fasta"
)


# ---------------------------------------------------------
# Fetch GenBank record from NCBI
# ---------------------------------------------------------

def fetch_genbank_record():
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    handle = Entrez.efetch(
        db="nuccore",
        id=MRNA_ACCESSION,
        rettype="gb",
        retmode="text",
    )

    record = SeqIO.read(handle, "genbank")

    handle.close()

    # Save the original NCBI GenBank record
    SeqIO.write(
        record,
        GENBANK_FILE,
        "genbank"
    )

    return record


# ---------------------------------------------------------
# Extract CDS from GenBank record
# ---------------------------------------------------------

def extract_cds(record):

    cds_features = [
        feature
        for feature in record.features
        if feature.type == "CDS"
    ]

    if not cds_features:
        raise RuntimeError(
            "No CDS feature found in NCBI record."
        )

    if len(cds_features) != 1:
        raise RuntimeError(
            f"Expected exactly 1 CDS feature, "
            f"but found {len(cds_features)}."
        )

    cds_feature = cds_features[0]

    # Extract CDS sequence
    cds_sequence = cds_feature.extract(record.seq)

    CDS_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # Save CDS FASTA
    with open(
        CDS_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        handle.write(
            f">{MRNA_ACCESSION}_CDS "
            f"OsNRAMP5 metal transporter Nramp5-like\n"
        )

        for i in range(
            0,
            len(cds_sequence),
            60
        ):
            handle.write(
                str(cds_sequence[i:i + 60])
                + "\n"
            )

    return cds_feature, cds_sequence


# ---------------------------------------------------------
# Verify CDS translation against NCBI annotation
# ---------------------------------------------------------

def verify_translation(
    cds_feature,
    cds_sequence
):

    # Translate CDS including the stop codon
    translated = str(
        cds_sequence.translate(
            to_stop=False
        )
    )

    # Get NCBI's annotated protein translation
    annotated_translation = (
        cds_feature.qualifiers
        .get("translation", [""])[0]
    )

    # Remove terminal stop symbol from our translation
    translated_without_stop = (
        translated.rstrip("*")
    )

    print("\nTranslation verification:")

    if (
        translated_without_stop
        == annotated_translation
    ):
        print(
            "PASS: CDS translation matches "
            "NCBI protein annotation."
        )
    else:
        print(
            "WARNING: CDS translation does not "
            "match NCBI protein annotation."
        )

    print(
        f"CDS nucleotide length: "
        f"{len(cds_sequence)} nt"
    )

    print(
        f"Translated length including stop: "
        f"{len(translated)}"
    )

    print(
        f"Protein amino-acid length: "
        f"{len(annotated_translation)} aa"
    )

    # Show whether the translation ends in a stop codon
    if translated.endswith("*"):
        print(
            "Terminal stop codon detected: YES"
        )
    else:
        print(
            "Terminal stop codon detected: NO"
        )

    return (
        translated_without_stop
        == annotated_translation
    )


# ---------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------

def main():

    print(
        "Fetching NCBI GenBank record..."
    )

    record = fetch_genbank_record()

    print(
        "Saved GenBank record:"
    )

    print(
        GENBANK_FILE
    )

    cds_feature, cds_sequence = extract_cds(
        record
    )

    # Get protein information from NCBI
    protein_id = (
        cds_feature.qualifiers
        .get("protein_id", ["unknown"])[0]
    )

    annotated_translation = (
        cds_feature.qualifiers
        .get("translation", [""])[0]
    )

    print()

    print(
        "CDS extracted successfully."
    )

    print(
        f"CDS length: "
        f"{len(cds_sequence)} nt"
    )

    print(
        f"CDS location: "
        f"{cds_feature.location}"
    )

    print(
        f"Protein ID: "
        f"{protein_id}"
    )

    print(
        f"Annotated protein length: "
        f"{len(annotated_translation)} aa"
    )

    print(
        f"Protein translation: "
        f"{annotated_translation[:30]}..."
    )

    # Verify CDS → protein relationship
    verification_passed = verify_translation(
        cds_feature,
        cds_sequence
    )

    print()

    if verification_passed:
        print(
            "BIOLOGICAL VALIDATION: PASS"
        )
        print(
            "CDS correctly translates to "
            "the NCBI annotated protein."
        )
    else:
        print(
            "BIOLOGICAL VALIDATION: FAILED"
        )
        print(
            "Please inspect the NCBI annotation "
            "before continuing."
        )

    print()

    print(
        "Saved CDS:"
    )

    print(
        CDS_FILE
    )


# ---------------------------------------------------------
# Entry point
# ---------------------------------------------------------

if __name__ == "__main__":
    main()
