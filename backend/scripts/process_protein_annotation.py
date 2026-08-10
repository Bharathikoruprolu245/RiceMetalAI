from pathlib import Path
import json

from Bio import SeqIO


# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------

ACCESSION = "NP_001389970.1"
GENE = "OsNRAMP5"

INPUT_FILE = Path(
    "data/raw/ncbi/OsNRAMP5_NP_001389970.1_protein_annotation.gb"
)

OUTPUT_FILE = Path(
    "data/processed/proteins/OsNRAMP5_NP_001389970.1_annotation.json"
)


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------

def get_first(qualifiers, key, default=None):
    values = qualifiers.get(key)

    if values:
        return values[0]

    return default


def location_to_residues(location):
    """
    Convert Biopython's zero-based, half-open coordinates
    into biological 1-based inclusive coordinates.
    """

    start = int(location.start) + 1
    end = int(location.end)

    return start, end


# ---------------------------------------------------------
# Main processing
# ---------------------------------------------------------

def main():

    print("Loading NCBI protein annotation...")

    record = SeqIO.read(INPUT_FILE, "genbank")

    print(f"Protein: {record.id}")
    print(f"Description: {record.description}")
    print(f"Length: {len(record.seq)} aa")

    domains = []
    transmembrane_regions = []
    disordered_regions = []

    # -----------------------------------------------------
    # Extract annotated features
    # -----------------------------------------------------

    for feature in record.features:

        start, end = location_to_residues(feature.location)

        qualifiers = feature.qualifiers

        # ---------------------------------------------
        # Domain / Region
        # ---------------------------------------------

        if feature.type == "Region":

            region_name = get_first(
                qualifiers,
                "region_name",
                ""
            )

            note = get_first(
                qualifiers,
                "note",
                ""
            )

            db_xref = qualifiers.get(
                "db_xref",
                []
            )

            # Detect SLC5-6 domain
            if region_name:

                domains.append(
                    {
                        "name": region_name,
                        "start": start,
                        "end": end,
                        "note": note,
                        "db_xref": db_xref
                    }
                )

            # Detect disordered region
            if "Disordered" in region_name:

                disordered_regions.append(
                    {
                        "start": start,
                        "end": end,
                        "evidence": note
                    }
                )

        # ---------------------------------------------
        # Transmembrane regions
        # ---------------------------------------------

        elif feature.type == "Site":

            site_type = get_first(
                qualifiers,
                "site_type",
                ""
            )

            note = get_first(
                qualifiers,
                "note",
                ""
            )

            if site_type == "transmembrane region":

                transmembrane_regions.append(
                    {
                        "start": start,
                        "end": end,
                        "site_type": site_type,
                        "evidence": note
                    }
                )

    # -----------------------------------------------------
    # Build output
    # -----------------------------------------------------

    annotation = {

        "protein": {
            "accession": ACCESSION,
            "gene": GENE,
            "description": record.description,
            "organism": "Oryza sativa Japonica Group",
            "sequence_length_aa": len(record.seq)
        },

        "domains": domains,

        "transmembrane_regions": transmembrane_regions,

        "disordered_regions": disordered_regions,

        "summary": {
            "domain_count": len(domains),
            "transmembrane_region_count": len(
                transmembrane_regions
            ),
            "disordered_region_count": len(
                disordered_regions
            )
        },

        "source": {
            "database": "NCBI RefSeq",
            "accession": ACCESSION,
            "source_file": str(INPUT_FILE)
        }
    }

    # -----------------------------------------------------
    # Save
    # -----------------------------------------------------

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            annotation,
            handle,
            indent=4
        )

    # -----------------------------------------------------
    # Report
    # -----------------------------------------------------

    print()
    print("Protein annotation processing complete.")
    print()

    print(
        f"Domains: "
        f"{len(domains)}"
    )

    print(
        f"Transmembrane regions: "
        f"{len(transmembrane_regions)}"
    )

    print(
        f"Disordered regions: "
        f"{len(disordered_regions)}"
    )

    print()
    print("Results saved to:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
