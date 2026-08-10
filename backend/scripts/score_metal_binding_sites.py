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

CONSERVATION_FILE = (
    BASE_DIR
    / "data"
    / "processed"
    / "alignments"
    / "OsNRAMP5_residue_conservation.json"
)

METAL_BINDING_FILE = (
    BASE_DIR
    / "data"
    / "processed"
    / "structures"
    / "8E6N_metal_binding.json"
)

OUTPUT_DIR = (
    BASE_DIR
    / "data"
    / "processed"
    / "alignments"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "OsNRAMP5_metal_binding_candidates.json"
)


# ============================================================
# SCORING WEIGHTS
# ============================================================

# Total score = 100

PROXIMITY_WEIGHT = 40
CONSERVATION_WEIGHT = 40
MAPPING_WEIGHT = 20


# ============================================================
# CONSERVATION SCORES
# ============================================================

CONSERVATION_SCORES = {

    "identical": 40,

    "conservative_substitution": 30,

    "chemically_similar": 25,

    "partially_conserved": 15,

    "chemically_different": 0,

    "unmapped": 0,

    "unknown": 0,

    "unknown_reference_residue": 0,

    "unknown_target_residue": 0,
}


# ============================================================
# LOAD JSON
# ============================================================

def load_json(path):

    if not path.exists():

        raise FileNotFoundError(
            f"Required file not found:\n{path}"
        )

    with open(
        path,
        "r",
        encoding="utf-8",
    ) as f:

        return json.load(f)


# ============================================================
# PROXIMITY SCORE
# ============================================================

def calculate_proximity_score(distance):

    """
    Convert reference metal distance into a
    normalized 0-40 score.

    Maximum evidence is assigned to residues
    at 2.0 Å or closer.

    Residues at or beyond 3.5 Å receive
    zero proximity score.

    This is a heuristic evidence score,
    NOT a physical binding probability.
    """

    cutoff = 3.5
    strong_distance = 2.0

    if distance is None:
        return 0.0

    distance = float(distance)

    if distance <= strong_distance:
        return float(PROXIMITY_WEIGHT)

    if distance >= cutoff:
        return 0.0

    fraction = (
        cutoff - distance
    ) / (
        cutoff - strong_distance
    )

    score = fraction * PROXIMITY_WEIGHT

    return round(
        max(0.0, min(PROXIMITY_WEIGHT, score)),
        2,
    )


# ============================================================
# CONSERVATION SCORE
# ============================================================

def calculate_conservation_score(
    classification
):

    raw_score = CONSERVATION_SCORES.get(
        classification,
        0,
    )

    # Normalize from 0-40
    return round(
        raw_score,
        2,
    )


# ============================================================
# BUILD REFERENCE DISTANCE LOOKUP
# ============================================================

def build_distance_lookup(
    metal_binding_data
):

    lookup = {}

    environments = metal_binding_data.get(
        "metal_binding_environments",
        [],
    )

    for environment in environments:

        metal = environment.get(
            "metal",
            {}
        )

        metal_residue = metal.get(
            "residue_number"
        )

        metal_name = metal.get(
            "residue_name"
        )

        for residue in environment.get(
            "nearby_residues",
            []
        ):

            reference_position = residue.get(
                "residue_number"
            )

            distance = residue.get(
                "distance_angstrom"
            )

            key = (
                metal_residue,
                reference_position,
            )

            lookup[key] = {
                "metal_residue": metal_residue,
                "metal": metal_name,
                "reference_position": reference_position,
                "reference_residue": residue.get(
                    "residue_name"
                ),
                "distance_angstrom": distance,
                "closest_atom": residue.get(
                    "closest_atom"
                ),
                "chain": residue.get(
                    "chain"
                ),
            }

    return lookup


# ============================================================
# BUILD CONSERVATION LOOKUP
# ============================================================

def build_conservation_lookup(
    conservation_data
):

    lookup = {}

    residues = conservation_data.get(
        "residue_conservation",
        []
    )

    for residue in residues:

        reference_position = residue.get(
            "reference_position"
        )

        metal_residue = residue.get(
            "metal_residue"
        )

        key = (
            metal_residue,
            reference_position,
        )

        conservation = residue.get(
            "conservation",
            {}
        )

        lookup[key] = {
            "classification": conservation.get(
                "classification",
                "unknown"
            ),

            "reference_residue_1letter": (
                conservation.get(
                    "reference_residue_1letter"
                )
            ),

            "reference_properties": (
                conservation.get(
                    "reference_properties"
                )
            ),

            "osnramp5_properties": (
                conservation.get(
                    "osnramp5_properties"
                )
            ),

            "charge_conserved": (
                conservation.get(
                    "charge_conserved"
                )
            ),

            "polarity_conserved": (
                conservation.get(
                    "polarity_conserved"
                )
            ),

            "hydrophobicity_conserved": (
                conservation.get(
                    "hydrophobicity_conserved"
                )
            ),
        }

    return lookup


# ============================================================
# SCORE ONE CANDIDATE
# ============================================================

def score_candidate(
    mapping,
    distance_lookup,
    conservation_lookup,
):

    metal_residue = mapping.get(
        "metal_residue"
    )

    reference_position = mapping.get(
        "reference_position"
    )

    target_position = mapping.get(
        "osnramp5_position"
    )

    target_residue = mapping.get(
        "osnramp5_residue"
    )

    reference_residue = mapping.get(
        "reference_residue"
    )

    mapping_status = mapping.get(
        "mapping_status",
        "unknown"
    )

    key = (
        metal_residue,
        reference_position,
    )

    # --------------------------------------------------------
    # Reference proximity
    # --------------------------------------------------------

    reference_data = distance_lookup.get(
        key,
        {}
    )

    distance = reference_data.get(
        "distance_angstrom"
    )

    proximity_score = (
        calculate_proximity_score(distance)
    )

    # --------------------------------------------------------
    # Conservation
    # --------------------------------------------------------

    conservation_data = (
        conservation_lookup.get(
            key,
            {}
        )
    )

    classification = (
        conservation_data.get(
            "classification",
            "unknown"
        )
    )

    conservation_score = (
        calculate_conservation_score(
            classification
        )
    )

    # --------------------------------------------------------
    # Mapping score
    # --------------------------------------------------------

    if mapping_status == "mapped":
        mapping_score = MAPPING_WEIGHT
    else:
        mapping_score = 0.0

    # --------------------------------------------------------
    # Total
    # --------------------------------------------------------

    total_score = (
        proximity_score
        + conservation_score
        + mapping_score
    )

    total_score = round(
        total_score,
        2,
    )

    # --------------------------------------------------------
    # Ranking interpretation
    # --------------------------------------------------------

    if total_score >= 80:

        confidence = "high"

    elif total_score >= 60:

        confidence = "moderate"

    elif total_score >= 40:

        confidence = "low"

    else:

        confidence = "very_low"

    return {

        "candidate": {

            "osnramp5_position": target_position,

            "osnramp5_residue": target_residue,

            "reference_position": reference_position,

            "reference_residue": reference_residue,

            "metal_residue": metal_residue,

            "metal": mapping.get(
                "metal"
            ),
        },

        "reference_evidence": {

            "distance_angstrom": distance,

            "closest_atom": (
                reference_data.get(
                    "closest_atom"
                )
            ),

            "chain": (
                reference_data.get(
                    "chain"
                )
            ),

            "proximity_score": proximity_score,
        },

        "sequence_evidence": {

            "mapping_status": mapping_status,

            "conservation_classification": (
                classification
            ),

            "conservation_score": (
                conservation_score
            ),

            "reference_residue_1letter": (
                conservation_data.get(
                    "reference_residue_1letter"
                )
            ),

            "target_residue": target_residue,

            "charge_conserved": (
                conservation_data.get(
                    "charge_conserved"
                )
            ),

            "polarity_conserved": (
                conservation_data.get(
                    "polarity_conserved"
                )
            ),

            "hydrophobicity_conserved": (
                conservation_data.get(
                    "hydrophobicity_conserved"
                )
            ),
        },

        "mapping_score": mapping_score,

        "total_candidate_score": total_score,

        "confidence": confidence,

        "interpretation": (
            "Heuristic candidate score combining "
            "reference metal proximity, sequence "
            "mapping, and biochemical conservation. "
            "This score is not a statistical binding "
            "probability and does not prove metal "
            "coordination in OsNRAMP5."
        ),
    }


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "Loading OsNRAMP5 metal-binding evidence..."
    )

    mapping_data = load_json(
        MAPPING_FILE
    )

    conservation_data = load_json(
        CONSERVATION_FILE
    )

    metal_binding_data = load_json(
        METAL_BINDING_FILE
    )

    target = mapping_data.get(
        "target",
        {}
    )

    reference = mapping_data.get(
        "reference",
        {}
    )

    print(
        f"Target: "
        f"{target.get('gene')} "
        f"({target.get('accession')})"
    )

    print(
        f"Reference: "
        f"{reference.get('structure')} "
        f"chain {reference.get('chain')}"
    )

    print()

    # --------------------------------------------------------
    # Build lookups
    # --------------------------------------------------------

    distance_lookup = (
        build_distance_lookup(
            metal_binding_data
        )
    )

    conservation_lookup = (
        build_conservation_lookup(
            conservation_data
        )
    )

    mappings = mapping_data.get(
        "metal_binding_residue_mapping",
        []
    )

    print(
        f"Candidate reference sites: "
        f"{len(mappings)}"
    )

    print()
    print(
        "Calculating candidate scores..."
    )
    print()

    candidates = []

    # --------------------------------------------------------
    # Score candidates
    # --------------------------------------------------------

    for mapping in mappings:

        result = score_candidate(
            mapping,
            distance_lookup,
            conservation_lookup,
        )

        candidates.append(
            result
        )

    # --------------------------------------------------------
    # Rank candidates
    # --------------------------------------------------------

    candidates.sort(
        key=lambda x: x[
            "total_candidate_score"
        ],
        reverse=True,
    )

    # --------------------------------------------------------
    # Assign rank
    # --------------------------------------------------------

    for index, candidate in enumerate(
        candidates,
        start=1
    ):

        candidate["rank"] = index

    # --------------------------------------------------------
    # Print results
    # --------------------------------------------------------

    print(
        "Ranked OsNRAMP5 candidate sites:"
    )

    print()

    for candidate in candidates:

        c = candidate[
            "candidate"
        ]

        score = candidate[
            "total_candidate_score"
        ]

        confidence = candidate[
            "confidence"
        ]

        distance = candidate[
            "reference_evidence"
        ][
            "distance_angstrom"
        ]

        conservation = candidate[
            "sequence_evidence"
        ][
            "conservation_classification"
        ]

        print(
            f"Rank {candidate['rank']}: "
            f"{c['osnramp5_residue']}"
            f"{c['osnramp5_position']}"
        )

        print(
            f"  Reference: "
            f"{c['reference_residue']}"
            f"{c['reference_position']}"
        )

        print(
            f"  Metal: "
            f"{c['metal']}"
            f" {c['metal_residue']}"
        )

        print(
            f"  Reference distance: "
            f"{distance} Å"
        )

        print(
            f"  Conservation: "
            f"{conservation}"
        )

        print(
            f"  Candidate score: "
            f"{score}/100"
        )

        print(
            f"  Evidence level: "
            f"{confidence}"
        )

        print()

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    high = sum(
        1
        for candidate in candidates
        if candidate["confidence"] == "high"
    )

    moderate = sum(
        1
        for candidate in candidates
        if candidate["confidence"] == "moderate"
    )

    low = sum(
        1
        for candidate in candidates
        if candidate["confidence"] == "low"
    )

    very_low = sum(
        1
        for candidate in candidates
        if candidate["confidence"] == "very_low"
    )

    # --------------------------------------------------------
    # Output JSON
    # --------------------------------------------------------

    result = {

        "target": target,

        "reference": reference,

        "method": {

            "name": (
                "Heuristic metal-binding candidate "
                "evidence scoring"
            ),

            "components": {

                "reference_metal_proximity": (
                    f"{PROXIMITY_WEIGHT}%"
                ),

                "sequence_conservation": (
                    f"{CONSERVATION_WEIGHT}%"
                ),

                "sequence_mapping": (
                    f"{MAPPING_WEIGHT}%"
                ),
            },

            "note": (
                "Scores are heuristic evidence scores "
                "and must not be interpreted as "
                "statistical binding probabilities."
            ),
        },

        "summary": {

            "total_candidates": len(
                candidates
            ),

            "high_evidence_candidates": high,

            "moderate_evidence_candidates": (
                moderate
            ),

            "low_evidence_candidates": low,

            "very_low_evidence_candidates": (
                very_low
            ),
        },

        "ranked_candidates": candidates,

        "interpretation": (
            "Candidate sites were ranked by combining "
            "distance to manganese in the reference "
            "NRAMP structure, successful sequence "
            "mapping, and biochemical conservation "
            "between reference and OsNRAMP5 residues. "
            "Higher scores indicate stronger support "
            "from the currently available computational "
            "evidence. These rankings do not establish "
            "direct metal coordination, binding affinity, "
            "or experimental probability of binding in "
            "OsNRAMP5."
        ),
    }

    # --------------------------------------------------------
    # Save
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

    print(
        "Metal-binding candidate scoring complete."
    )

    print()

    print(
        "Results saved to:"
    )

    print(
        OUTPUT_FILE
    )


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    main()
