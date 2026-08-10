from pathlib import Path

from Bio import Entrez


Entrez.email = "bharathikoruprolu25@gmail.com"

GENE_ID = "4342859"

OUTPUT_DIR = Path("data/raw/ncbi")
OUTPUT_FILE = OUTPUT_DIR / f"OsNRAMP5_gene_{GENE_ID}.txt"


def fetch_gene():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    handle = Entrez.efetch(
        db="gene",
        id=GENE_ID,
        rettype="gene_table",
        retmode="text",
    )

    data = handle.read()
    handle.close()

    OUTPUT_FILE.write_text(data, encoding="utf-8")

    return data


if __name__ == "__main__":
    print("Fetching NCBI Gene record...")

    data = fetch_gene()

    print(data)

    print("\nRaw NCBI record saved to:")
    print(OUTPUT_FILE)
