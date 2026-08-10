from pathlib import Path
from urllib.request import urlopen
import json


# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------

PDB_ID = "8E6N"

PDB_URL = (
    f"https://files.rcsb.org/download/{PDB_ID}.pdb"
)

OUTPUT_FILE = Path(
    f"data/structures/pdb/{PDB_ID}.pdb"
)

METADATA_FILE = Path(
    f"data/structures/metadata/{PDB_ID}.json"
)


# ---------------------------------------------------------
# Metadata
# ---------------------------------------------------------

METADATA = {
    "pdb_id": PDB_ID,
    "title": (
        "X-ray structure of the Deinococcus radiodurans "
        "Nramp/MntH divalent transition metal transporter "
        "G223W mutant in an outward-open, "
        "manganese-bound state"
    ),
    "organism": "Deinococcus radiodurans",
    "protein": "Divalent metal cation transporter MntH",
    "gene": "mntH",
    "uniprot": "Q9RTP8",
    "experimental_method": "X-ray diffraction",
    "resolution_angstrom": 2.40,
    "metal": "Manganese (II)",
    "state": "outward-open",
    "mutation": "G223W",
    "chain": "A",
    "modeled_residues": 398,
    "deposited_residues": 411,
    "source_database": "RCSB Protein Data Bank",
    "source_url": (
        "https://www.rcsb.org/structure/8E6N"
    )
}


# ---------------------------------------------------------
# Download
# ---------------------------------------------------------

def download_structure():

    print(f"Downloading PDB structure {PDB_ID}...")

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with urlopen(PDB_URL) as response:

        data = response.read()

    OUTPUT_FILE.write_bytes(data)

    print("Structure downloaded:")
    print(OUTPUT_FILE)

    print(
        f"File size: "
        f"{OUTPUT_FILE.stat().st_size} bytes"
    )


# ---------------------------------------------------------
# Save metadata
# ---------------------------------------------------------

def save_metadata():

    METADATA_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        METADATA_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            METADATA,
            handle,
            indent=4
        )

    print()
    print("Metadata saved:")
    print(METADATA_FILE)


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    download_structure()

    save_metadata()

    print()
    print("PDB structure acquisition complete.")


if __name__ == "__main__":
    main()
