"""
Map structurally relevant metal-binding residues from the 8E6N
NRAMP reference structure onto OsNRAMP5 using the previously
generated global sequence alignment.

Input:
    data/processed/alignments/OsNRAMP5_vs_8E6N_alignment.json

Output:
    data/processed/alignments/OsNRAMP5_metal_binding_mapping.json
"""

import json
from pathlib import Path


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

ALIGNMENT_FILE = Path(
    "data/processed/alignments/"
    "OsNRAMP5_vs_8E6N_alignment.json"
)

OUTPUT_FILE = Path(
    "data/processed/alignments/"
    "OsNRAMP5_metal_binding_mapping.json"
)


# ---------------------------------------------------------------------
# Reference metal-associated residues from 8E6N
# ---------------------------------------------------------------------
#
# These residues came from our previous distance-based analysis:
#
# MN 501:
#   ALA53
#   ASP56
#   ASN59
#   MET230
#
# MN 502:
#   HIS287
#
# Important:
# These are residues proximal to the detected manganese atoms.
# Proximity alone does NOT prove direct coordination.
# ---------------------------------------------------------------------

METAL_BINDING_RESIDUES = [
    {
        "metal_residue": 501,
        "metal": "MN",
        "reference_position": 53,
        "reference_residue": "ALA",
    },
    {
        "metal_residue": 501,
        "metal": "MN",
        "reference_position": 56,
        "reference_residue": "ASP",
    },
    {
        "metal_residue": 501,
        "metal": "MN",
        "reference_position": 59,
        "reference_residue": "ASN",
    },
    {
        "metal_residue": 501,
        "metal": "MN",
        "reference_position": 230,
        "reference_residue": "MET",
    },
    {
        "metal_residue": 502,
        "metal": "MN",
        "reference_position": 287,
        "reference_residue": "HIS",
    },
]


# ---------------------------------------------------------------------
# Load alignment
# ---------------------------------------------------------------------

def load_alignment():
    """Load the previously generated alignment JSON."""

    if not ALIGNMENT_FILE.exists():
        raise FileNotFoundError(
            f"Alignment file not found:\n{ALIGNMENT_FILE}\n\n"
            "Run align_nramp_reference.py first."
        )

    with ALIGNMENT_FILE.open("r") as handle:
        return json.load(handle)


# ---------------------------------------------------------------------
# Build reference-position mapping
# ---------------------------------------------------------------------

def build_reference_mapping(alignment):
    """
    Build a dictionary:

        reference_position -> alignment record

    Only reference residues that actually exist in the alignment
    are included.
    """

    mapping = {}

    for record in alignment["residue_mapping"]:

        reference_position = record.get("reference_position")

        if reference_position is None:
            continue

        mapping[int(reference_position)] = record

    return mapping


# ---------------------------------------------------------------------
# Determine residue relationship
# ---------------------------------------------------------------------

def classify_residue(reference_residue, target_residue):
    """
    Classify the mapped residue.

    Note:
    This is a simple sequence-level comparison and not a structural
    or biochemical equivalence prediction.
    """

    if target_residue is None:
        return "not_mapped"

    if reference_residue == target_residue:
        return "identical"

    return "different"


# ---------------------------------------------------------------------
# Map metal-binding residues
# ---------------------------------------------------------------------

def map_metal_binding_residues(alignment):
    """Map each reference metal-associated residue to OsNRAMP5."""

    reference_mapping = build_reference_mapping(alignment)

    results = []

    for site in METAL_BINDING_RESIDUES:

        reference_position = site["reference_position"]

        alignment_record = reference_mapping.get(reference_position)

        if alignment_record is None:

            result = {
                "metal_residue": site["metal_residue"],
                "metal": site["metal"],
                "reference_position": reference_position,
                "reference_residue": site["reference_residue"],
                "osnramp5_position": None,
                "osnramp5_residue": None,
                "mapping_status": "not_mapped",
                "residue_relationship": "not_mapped",
            }

        else:

            target_position = alignment_record.get(
                "osnramp5_position"
            )

            target_residue = alignment_record.get(
                "osnramp5_residue"
            )

            relationship = classify_residue(
                site["reference_residue"],
                target_residue,
            )

            result = {
                "metal_residue": site["metal_residue"],
                "metal": site["metal"],
                "reference_position": reference_position,
                "reference_residue": site["reference_residue"],
                "osnramp5_position": target_position,
                "osnramp5_residue": target_residue,
                "mapping_status": (
                    "mapped"
                    if target_position is not None
                    else "not_mapped"
                ),
                "residue_relationship": relationship,
            }

        results.append(result)

    return results


# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------

def create_summary(results):
    """Create simple mapping statistics."""

    mapped = [
        result
        for result in results
        if result["mapping_status"] == "mapped"
    ]

    identical = [
        result
        for result in mapped
        if result["residue_relationship"] == "identical"
    ]

    different = [
        result
        for result in mapped
        if result["residue_relationship"] == "different"
    ]

    not_mapped = [
        result
        for result in results
        if result["mapping_status"] == "not_mapped"
    ]

    return {
        "total_reference_sites": len(results),
        "mapped_sites": len(mapped),
        "identical_residues": len(identical),
        "different_residues": len(different),
        "unmapped_sites": len(not_mapped),
    }


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    print("Loading OsNRAMP5 / 8E6N alignment...")
    alignment = load_alignment()

    print(
        f"Target: {alignment['target']['gene']} "
        f"({alignment['target']['accession']})"
    )

    print(
        f"Reference: {alignment['reference']['structure']} "
        f"chain {alignment['reference']['chain']}"
    )

    print()

    print("Mapping reference metal-associated residues...")
    print()

    results = map_metal_binding_residues(alignment)

    for result in results:

        reference = (
            f"{result['reference_residue']}"
            f"{result['reference_position']}"
        )

        if result["mapping_status"] == "mapped":

            target = (
                f"{result['osnramp5_residue']}"
                f"{result['osnramp5_position']}"
            )

            print(
                f"MN {result['metal_residue']}: "
                f"{reference} -> {target} "
                f"({result['residue_relationship']})"
            )

        else:

            print(
                f"MN {result['metal_residue']}: "
                f"{reference} -> NOT MAPPED"
            )

    summary = create_summary(results)

    print()
    print("Mapping summary:")
    print(
        f"Reference sites: "
        f"{summary['total_reference_sites']}"
    )
    print(
        f"Mapped sites: "
        f"{summary['mapped_sites']}"
    )
    print(
        f"Identical residues: "
        f"{summary['identical_residues']}"
    )
    print(
        f"Different residues: "
        f"{summary['different_residues']}"
    )
    print(
        f"Unmapped sites: "
        f"{summary['unmapped_sites']}"
    )

    output = {
        "target": {
            "gene": alignment["target"]["gene"],
            "accession": alignment["target"]["accession"],
            "sequence_length_aa": alignment["target"][
                "sequence_length_aa"
            ],
        },
        "reference": {
            "structure": alignment["reference"]["structure"],
            "chain": alignment["reference"]["chain"],
            "protein": alignment["reference"]["protein"],
            "organism": alignment["reference"]["organism"],
            "sequence_length_aa": alignment["reference"][
                "sequence_length_aa"
            ],
        },
        "source_analysis": {
            "source": (
                "data/processed/structures/"
                "8E6N_metal_binding.json"
            ),
            "description": (
                "Residues identified as structurally proximal "
                "to detected manganese atoms in the 8E6N "
                "reference structure."
            ),
            "distance_cutoff_angstrom": 3.5,
        },
        "mapping_method": {
            "alignment_file": str(ALIGNMENT_FILE),
            "method": "Position mapping through global sequence alignment",
            "note": (
                "Mapped residues represent sequence correspondence "
                "to structurally proximal reference residues. "
                "They do not by themselves establish direct metal "
                "coordination or functional equivalence."
            ),
        },
        "summary": summary,
        "metal_binding_residue_mapping": results,
        "interpretation": (
            "Reference metal-associated residues were mapped "
            "onto OsNRAMP5 using the previously generated "
            "global sequence alignment. Identical residues "
            "represent sequence conservation at the mapped "
            "position, whereas different residues indicate "
            "sequence divergence. Structural and biochemical "
            "validation would be required before assigning "
            "metal-binding function to the corresponding "
            "OsNRAMP5 residues."
        ),
    }

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with OUTPUT_FILE.open("w") as handle:
        json.dump(
            output,
            handle,
            indent=4
        )

    print()
    print("Metal-binding residue mapping complete.")
    print()
    print("Results saved to:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
