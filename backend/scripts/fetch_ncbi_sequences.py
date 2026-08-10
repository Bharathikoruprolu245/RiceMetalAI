from pathlib import Path

from Bio import Entrez, SeqIO


Entrez.email = "bharathikoruprolu25@gmail.com"

MRNA_ACCESSION = "NM_001403041.1"
PROTEIN_ACCESSION = "NP_001389970.1"

OUTPUT_CDS = Path("data/sequences/cds/OsNRAMP5_NM_001403041.1.fasta")
OUTPUT_PROTEIN = Path("data/sequences/proteins/OsNRAMP5_NP_001389970.1.fasta")


def fetch_sequence(accession: str, output_file: Path):
    handle = Entrez.efetch(
        db="nuccore" if accession.startswith("NM_") else "protein",
        id=accession,
        rettype="fasta",
        retmode="text",
    )

    record = SeqIO.read(handle, "fasta")
    handle.close()

    output_file.parent.mkdir(parents=True, exist_ok=True)

    SeqIO.write(record, output_file, "fasta")

    return record


def main():
    print("Fetching OsNRAMP5 mRNA/CDS source record...")
    mrna_record = fetch_sequence(
        MRNA_ACCESSION,
        OUTPUT_CDS,
    )

    print(f"Saved: {OUTPUT_CDS}")
    print(f"Sequence length: {len(mrna_record.seq)} nt")

    print("\nFetching OsNRAMP5 protein...")
    protein_record = fetch_sequence(
        PROTEIN_ACCESSION,
        OUTPUT_PROTEIN,
    )

    print(f"Saved: {OUTPUT_PROTEIN}")
    print(f"Protein length: {len(protein_record.seq)} aa")


if __name__ == "__main__":
    main()
