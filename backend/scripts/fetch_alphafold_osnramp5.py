from pathlib import Path
import json
import requests


# ============================================================
# Configuration
# ============================================================

UNIPROT_ACCESSION = "Q8H4H5"
GENE = "OsNRAMP5"

BASE_DIR = Path(__file__).resolve().parents[2]

OUTPUT_DIR = BASE_DIR / "data" / "structures" / "predicted"
METADATA_DIR = BASE_DIR / "data" / "structures" / "metadata"

PDB_FILE = OUTPUT_DIR / "AF-Q8H4H5-F1-model_v6.pdb"
METADATA_FILE = METADATA_DIR / "AF-Q8H4H5-F1-model_v6.json"

MODEL_URL = (
    "https://alphafold.ebi.ac.uk/files/"
    "AF-Q8H4H5-F1-model_v6.pdb"
)


# ============================================================
# Download AlphaFold model
# ============================================================

def download_model():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    METADATA_DIR.mkdir(parents=True, exist_ok=True)

    print("Downloading AlphaFold OsNRAMP5 structure...")
    print(f"UniProt: {UNIPROT_ACCESSION}")
    print(f"Model URL: {MODEL_URL}")
    print()

    response = requests.get(
        MODEL_URL,
        timeout=60
    )

    response.raise_for_status()

    if not response.content:
        raise RuntimeError(
            "Downloaded AlphaFold file is empty."
        )

    PDB_FILE.write_bytes(response.content)

    print("AlphaFold structure downloaded:")
    print(PDB_FILE)
    print(
        f"File size: {PDB_FILE.stat().st_size} bytes"
    )


# ============================================================
# Save metadata
# ============================================================

def save_metadata():

    metadata = {
        "gene": GENE,
        "uniprot_accession": UNIPROT_ACCESSION,
        "protein": "Metal transporter Nramp5",
        "organism": "Oryza sativa Japonica Group",
        "database": "AlphaFold Protein Structure Database",
        "model": "AF-Q8H4H5-F1-model_v6",
        "version": 6,
        "model_file": str(
            PDB_FILE.relative_to(BASE_DIR)
        ),
        "source_url": MODEL_URL,
        "structure_type": "AlphaFold predicted structure",
        "experimental": False,
        "note": (
            "Computationally predicted protein structure. "
            "It should not be treated as an experimentally "
            "determined structure."
        )
    }

    with open(
        METADATA_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            metadata,
            handle,
            indent=4
        )

    print()
    print("Metadata saved:")
    print(METADATA_FILE)


# ============================================================
# Main
# ============================================================

def main():

    download_model()
    save_metadata()

    print()
    print("AlphaFold OsNRAMP5 structure acquisition complete.")


if __name__ == "__main__":
    main()
