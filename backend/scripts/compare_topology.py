from pathlib import Path
import json


# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------

TOPOLOGY_FILE = Path(
    "data/processed/proteins/"
    "OsNRAMP5_NP_001389970.1_topology.json"
)

ANNOTATION_FILE = Path(
    "data/processed/proteins/"
    "OsNRAMP5_NP_001389970.1_annotation.json"
)

OUTPUT_FILE = Path(
    "data/processed/proteins/"
    "OsNRAMP5_NP_001389970.1_topology_comparison.json"
)


# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------

def load_json(path):
    """Load a JSON file."""

    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def calculate_overlap(start1, end1, start2, end2):
    """
    Calculate the number of overlapping residues
    between two 1-based inclusive regions.
    """

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start > overlap_end:
        return 0

    return overlap_end - overlap_start + 1


def compare_regions(screening_regions, ncbi_regions):
    """
    Compare computational hydrophobic regions against
    NCBI annotated transmembrane regions.
    """

    comparisons = []

    for screening_index, screening in enumerate(
        screening_regions,
        start=1
    ):

        screening_start = screening["start"]
        screening_end = screening["end"]

        best_match = None
        best_overlap = 0

        for ncbi_index, ncbi in enumerate(
            ncbi_regions,
            start=1
        ):

            ncbi_start = ncbi["start"]
            ncbi_end = ncbi["end"]

            overlap = calculate_overlap(
                screening_start,
                screening_end,
                ncbi_start,
                ncbi_end
            )

            if overlap > best_overlap:
                best_overlap = overlap
                best_match = {
                    "ncbi_region": ncbi_index,
                    "start": ncbi_start,
                    "end": ncbi_end,
                    "overlap_residues": overlap
                }

        # -------------------------------------------------
        # Calculate screening region length
        # -------------------------------------------------

        screening_length = (
            screening_end - screening_start + 1
        )

        # -------------------------------------------------
        # Calculate overlap percentage
        # -------------------------------------------------

        if screening_length > 0:
            overlap_percentage = round(
                (best_overlap / screening_length) * 100,
                2
            )
        else:
            overlap_percentage = 0.0

        # -------------------------------------------------
        # Determine agreement
        # -------------------------------------------------

        if best_overlap == 0:
            agreement = "no_overlap"

        elif overlap_percentage >= 50:
            agreement = "strong_overlap"

        else:
            agreement = "partial_overlap"

        comparisons.append(
            {
                "screening_region": screening_index,
                "screening_start": screening_start,
                "screening_end": screening_end,
                "screening_length": screening_length,
                "max_hydrophobic_fraction": screening[
                    "max_hydrophobic_fraction"
                ],
                "best_ncbi_match": best_match,
                "overlap_percentage": overlap_percentage,
                "agreement": agreement
            }
        )

    return comparisons


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():

    print("Loading topology screening results...")

    topology = load_json(TOPOLOGY_FILE)

    print("Loading NCBI protein annotation...")

    annotation = load_json(ANNOTATION_FILE)

    # -----------------------------------------------------
    # Extract regions
    # -----------------------------------------------------

    screening_regions = topology[
        "candidate_hydrophobic_regions"
    ]

    ncbi_regions = annotation[
        "transmembrane_regions"
    ]

    print()
    print(
        f"Computational hydrophobic regions: "
        f"{len(screening_regions)}"
    )

    print(
        f"NCBI annotated TM regions: "
        f"{len(ncbi_regions)}"
    )

    # -----------------------------------------------------
    # Compare
    # -----------------------------------------------------

    comparisons = compare_regions(
        screening_regions,
        ncbi_regions
    )

    # -----------------------------------------------------
    # Summary statistics
    # -----------------------------------------------------

    strong_matches = sum(
        1
        for item in comparisons
        if item["agreement"] == "strong_overlap"
    )

    partial_matches = sum(
        1
        for item in comparisons
        if item["agreement"] == "partial_overlap"
    )

    no_matches = sum(
        1
        for item in comparisons
        if item["agreement"] == "no_overlap"
    )

    total_screening_regions = len(
        screening_regions
    )

    if total_screening_regions > 0:

        agreement_percentage = round(
            (
                (
                    strong_matches
                    + partial_matches
                )
                / total_screening_regions
            )
            * 100,
            2
        )

    else:

        agreement_percentage = 0.0

    # -----------------------------------------------------
    # Output
    # -----------------------------------------------------

    result = {

        "protein": {
            "accession": "NP_001389970.1",
            "gene": "OsNRAMP5",
            "sequence_length_aa": 538
        },

        "methods": {

            "computational_screening": {
                "method": topology["method"]["analysis"],
                "window_size": topology["method"]["window_size"],
                "hydrophobic_fraction_threshold": topology[
                    "method"
                ][
                    "hydrophobic_fraction_threshold"
                ]
            },

            "reference_annotation": {
                "database": "NCBI RefSeq",
                "source": annotation["source"]
            }
        },

        "counts": {

            "computational_regions": len(
                screening_regions
            ),

            "ncbi_transmembrane_regions": len(
                ncbi_regions
            ),

            "strong_overlaps": strong_matches,

            "partial_overlaps": partial_matches,

            "no_overlaps": no_matches
        },

        "agreement": {
            "regions_with_any_overlap": (
                strong_matches
                + partial_matches
            ),

            "agreement_percentage": agreement_percentage
        },

        "region_comparisons": comparisons,

        "interpretation": (
            "Computational hydrophobic regions were "
            "compared against NCBI annotated "
            "transmembrane regions. Strong overlap "
            "indicates agreement between the independent "
            "hydrophobic screening method and the "
            "reference annotation. Partial or absent "
            "overlap indicates regions requiring further "
            "structural investigation."
        )
    }

    # -----------------------------------------------------
    # Save result
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
            result,
            handle,
            indent=4
        )

    # -----------------------------------------------------
    # Print summary
    # -----------------------------------------------------

    print()
    print("Topology comparison complete.")
    print()

    print(
        f"Strong overlaps: {strong_matches}"
    )

    print(
        f"Partial overlaps: {partial_matches}"
    )

    print(
        f"No overlaps: {no_matches}"
    )

    print(
        f"Overall agreement: "
        f"{agreement_percentage}%"
    )

    print()
    print("Results saved to:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
