import json
from pathlib import Path


# ============================================================
# PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]

MAPPING_FILE = (
    BASE_DIR
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_metal_binding_mapping.json"
)

OUTPUT_DIR = (
    BASE_DIR
    / "data"
    / "processed"
    / "alignments"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "OsNRAMP5_residue_conservation.json"
)


# ============================================================
# AMINO ACID CONVERSION
# ============================================================

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


# ============================================================
# AMINO ACID BIOCHEMICAL PROPERTIES
# ============================================================

AMINO_ACID_PROPERTIES = {

    "A": {
        "name": "Alanine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "R": {
        "name": "Arginine",
        "charge": "positive",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "N": {
        "name": "Asparagine",
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "D": {
        "name": "Aspartic acid",
        "charge": "negative",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "C": {
        "name": "Cysteine",
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobic": True,
    },

    "Q": {
        "name": "Glutamine",
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "E": {
        "name": "Glutamic acid",
        "charge": "negative",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "G": {
        "name": "Glycine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": False,
    },

    "H": {
        "name": "Histidine",
        "charge": "positive/neutral",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "I": {
        "name": "Isoleucine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "L": {
        "name": "Leucine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "K": {
        "name": "Lysine",
        "charge": "positive",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "M": {
        "name": "Methionine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "F": {
        "name": "Phenylalanine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "P": {
        "name": "Proline",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "S": {
        "name": "Serine",
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "T": {
        "name": "Threonine",
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobic": False,
    },

    "W": {
        "name": "Tryptophan",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },

    "Y": {
        "name": "Tyrosine",
        "charge": "neutral",
        "polarity": "polar",
        "hydrophobic": True,
    },

    "V": {
        "name": "Valine",
        "charge": "neutral",
        "polarity": "nonpolar",
        "hydrophobic": True,
    },
}


# ============================================================
# CHEMICAL GROUPS
# ============================================================

CONSERVATIVE_GROUPS = [

    # Hydrophobic
    {"A", "V", "I", "L", "M"},
    {"F", "W", "Y"},

    # Small/polar
    {"S", "T"},
    {"N", "Q"},

    # Acidic
    {"D", "E"},

    # Basic
    {"K", "R"},

]


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def same_conservative_group(residue1, residue2):
    """
    Determine whether two amino acids belong to
    the same conservative biochemical group.
    """

    if residue1 == residue2:
        return True

    for group in CONSERVATIVE_GROUPS:
        if residue1 in group and residue2 in group:
            return True

    return False


def classify_conservation(reference, target):
    """
    Classify the relationship between the reference residue
    and the OsNRAMP5 residue.
    """

    if reference is None or target is None:
        return "unmapped"

    if reference == target:
        return "identical"

    if reference not in AMINO_ACID_PROPERTIES:
        return "unknown"

    if target not in AMINO_ACID_PROPERTIES:
        return "unknown"

    reference_properties = AMINO_ACID_PROPERTIES[reference]
    target_properties = AMINO_ACID_PROPERTIES[target]

    # Conservative substitution
    if same_conservative_group(reference, target):
        return "conservative_substitution"

    # Same charge and polarity
    if (
        reference_properties["charge"]
        == target_properties["charge"]
        and
        reference_properties["polarity"]
        == target_properties["polarity"]
    ):
        return "chemically_similar"

    # Same polarity
    if reference_properties["polarity"] == target_properties["polarity"]:
        return "partially_conserved"

    return "chemically_different"


def compare_properties(reference, target):
    """
    Compare biochemical properties between two residues.
    """

    if reference is None or target is None:
        return {
            "charge_conserved": None,
            "polarity_conserved": None,
            "hydrophobicity_conserved": None,
        }

    reference_properties = AMINO_ACID_PROPERTIES[reference]
    target_properties = AMINO_ACID_PROPERTIES[target]

    return {
        "charge_conserved": (
            reference_properties["charge"]
            == target_properties["charge"]
        ),

        "polarity_conserved": (
            reference_properties["polarity"]
            == target_properties["polarity"]
        ),

        "hydrophobicity_conserved": (
            reference_properties["hydrophobic"]
            == target_properties["hydrophobic"]
        ),
    }


# ============================================================
# RESIDUE ANALYSIS
# ============================================================

def analyze_residue(mapping):

    reference_three = mapping.get("reference_residue")
    target = mapping.get("osnramp5_residue")

    # --------------------------------------------------------
    # Convert reference residue from 3-letter to 1-letter
    # --------------------------------------------------------

    reference = THREE_TO_ONE.get(reference_three)

    if reference is None:

        return {
            **mapping,

            "conservation": {
                "reference_residue_1letter": None,
                "classification": "unknown_reference_residue",
                "reference_properties": None,
                "osnramp5_properties": None,
                "charge_conserved": None,
                "polarity_conserved": None,
                "hydrophobicity_conserved": None,
            },
        }

    # --------------------------------------------------------
    # Handle unmapped target residue
    # --------------------------------------------------------

    if target is None:

        return {
            **mapping,

            "conservation": {
                "reference_residue_1letter": reference,
                "classification": "unmapped",
                "reference_properties": (
                    AMINO_ACID_PROPERTIES[reference]
                ),
                "osnramp5_properties": None,
                "charge_conserved": None,
                "polarity_conserved": None,
                "hydrophobicity_conserved": None,
            },
        }

    # --------------------------------------------------------
    # Validate target residue
    # --------------------------------------------------------

    if target not in AMINO_ACID_PROPERTIES:

        return {
            **mapping,

            "conservation": {
                "reference_residue_1letter": reference,
                "classification": "unknown_target_residue",
                "reference_properties": (
                    AMINO_ACID_PROPERTIES[reference]
                ),
                "osnramp5_properties": None,
                "charge_conserved": None,
                "polarity_conserved": None,
                "hydrophobicity_conserved": None,
            },
        }

    # --------------------------------------------------------
    # Get properties
    # --------------------------------------------------------

    reference_properties = AMINO_ACID_PROPERTIES[reference]
    target_properties = AMINO_ACID_PROPERTIES[target]

    # --------------------------------------------------------
    # Classification
    # --------------------------------------------------------

    classification = classify_conservation(
        reference,
        target,
    )

    # --------------------------------------------------------
    # Property comparison
    # --------------------------------------------------------

    property_comparison = compare_properties(
        reference,
        target,
    )

    # --------------------------------------------------------
    # Return result
    # --------------------------------------------------------

    return {
        **mapping,

        "conservation": {

            "reference_residue_1letter": reference,

            "classification": classification,

            "reference_properties": reference_properties,

            "osnramp5_properties": target_properties,

            "charge_conserved": (
                property_comparison["charge_conserved"]
            ),

            "polarity_conserved": (
                property_comparison["polarity_conserved"]
            ),

            "hydrophobicity_conserved": (
                property_comparison["hydrophobicity_conserved"]
            ),
        },
    }


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    print("Loading metal-binding residue mapping...")

    if not MAPPING_FILE.exists():

        raise FileNotFoundError(
            f"Mapping file not found:\n{MAPPING_FILE}"
        )

    with open(
        MAPPING_FILE,
        "r",
        encoding="utf-8",
    ) as f:

        mapping_data = json.load(f)

    target_info = mapping_data.get("target", {})
    reference_info = mapping_data.get("reference", {})

    print(
        f"Target: "
        f"{target_info.get('gene')} "
        f"({target_info.get('accession')})"
    )

    print(
        f"Reference: "
        f"{reference_info.get('structure')} "
        f"chain {reference_info.get('chain')}"
    )

    print()
    print("Analyzing residue conservation...")
    print()

    mappings = mapping_data.get(
        "metal_binding_residue_mapping",
        []
    )

    analyzed_residues = []

    # Counters
    identical = 0
    conservative = 0
    chemically_similar = 0
    partially_conserved = 0
    chemically_different = 0
    unmapped = 0
    unknown = 0

    # --------------------------------------------------------
    # Analyze each residue
    # --------------------------------------------------------

    for mapping in mappings:

        result = analyze_residue(mapping)

        analyzed_residues.append(result)

        reference_three = mapping.get(
            "reference_residue",
            "?"
        )

        reference_position = mapping.get(
            "reference_position",
            "?"
        )

        target_residue = mapping.get(
            "osnramp5_residue"
        )

        target_position = mapping.get(
            "osnramp5_position"
        )

        classification = result[
            "conservation"
        ]["classification"]

        print(
            f"{reference_three}{reference_position}"
            f" -> "
            f"{target_residue}{target_position}"
            f" : {classification}"
        )

        # Count classifications
        if classification == "identical":
            identical += 1

        elif classification == "conservative_substitution":
            conservative += 1

        elif classification == "chemically_similar":
            chemically_similar += 1

        elif classification == "partially_conserved":
            partially_conserved += 1

        elif classification == "chemically_different":
            chemically_different += 1

        elif classification == "unmapped":
            unmapped += 1

        else:
            unknown += 1

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    total = len(analyzed_residues)

    chemically_conserved = (
        identical
        + conservative
        + chemically_similar
    )

    if total > 0:

        conservation_percentage = round(
            (chemically_conserved / total) * 100,
            2,
        )

    else:

        conservation_percentage = 0.0

    print()
    print("Residue conservation analysis complete.")
    print()

    print("Summary:")
    print(f"Total mapped reference sites: {total}")
    print(f"Identical residues: {identical}")
    print(
        f"Conservative substitutions: "
        f"{conservative}"
    )
    print(
        f"Chemically similar: "
        f"{chemically_similar}"
    )
    print(
        f"Partially conserved: "
        f"{partially_conserved}"
    )
    print(
        f"Chemically different: "
        f"{chemically_different}"
    )
    print(f"Unmapped: {unmapped}")
    print(f"Unknown: {unknown}")

    print()
    print(
        "Chemically conserved sites: "
        f"{chemically_conserved}/{total}"
    )

    print(
        "Chemical conservation percentage: "
        f"{conservation_percentage}%"
    )

    # --------------------------------------------------------
    # Final JSON
    # --------------------------------------------------------

    result = {

        "target": target_info,

        "reference": reference_info,

        "source_analysis": mapping_data.get(
            "source_analysis",
            {}
        ),

        "mapping_method": mapping_data.get(
            "mapping_method",
            {}
        ),

        "analysis_method": {
            "name": (
                "Biochemical residue conservation analysis"
            ),

            "reference_residue_conversion": (
                "NCBI three-letter amino acid code "
                "converted to one-letter code"
            ),

            "classification": [
                "identical",
                "conservative_substitution",
                "chemically_similar",
                "partially_conserved",
                "chemically_different",
            ],
        },

        "summary": {

            "total_sites": total,

            "identical": identical,

            "conservative_substitutions": conservative,

            "chemically_similar": chemically_similar,

            "partially_conserved": partially_conserved,

            "chemically_different": chemically_different,

            "unmapped": unmapped,

            "unknown": unknown,

            "chemically_conserved_sites": (
                chemically_conserved
            ),

            "chemical_conservation_percentage": (
                conservation_percentage
            ),
        },

        "residue_conservation": analyzed_residues,

        "interpretation": (
            "Reference residues located near bound manganese "
            "in the 8E6N NRAMP/MntH structure were mapped "
            "onto OsNRAMP5 through sequence alignment and "
            "compared for biochemical property conservation. "
            "Identical or chemically similar substitutions "
            "indicate conservation of selected amino-acid "
            "properties, whereas chemically different "
            "substitutions indicate sequence divergence. "
            "This analysis does not by itself establish "
            "metal coordination or functional equivalence "
            "in OsNRAMP5."
        ),
    }

    # --------------------------------------------------------
    # Save output
    # --------------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8",
    ) as f:

        json.dump(
            result,
            f,
            indent=4,
        )

    print()
    print("Results saved to:")
    print(OUTPUT_FILE)


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    main()
