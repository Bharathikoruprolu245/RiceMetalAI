#!/usr/bin/env python3

"""
Analyze the structural environment around OsNRAMP5
candidate metal-binding residues.

This is a structural screening step.

It does NOT prove metal coordination and does NOT perform docking.

Input:
    data/structures/predicted/AF-Q8H4H5-F1-model_v6.pdb

Candidate residues:
    A57
    D60
    N63
    A232
    M235
    Q337

Output:
    data/processed/structures/OsNRAMP5_metal_pocket.json
"""

import json
import math
from pathlib import Path

from Bio.PDB import PDBParser


# ============================================================
# PATHS
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

PDB_FILE = (
    PROJECT_ROOT
    / "data"
    / "structures"
    / "predicted"
    / "AF-Q8H4H5-F1-model_v6.pdb"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "structures"
)

OUTPUT_FILE = (
    OUTPUT_DIR
    / "OsNRAMP5_metal_pocket.json"
)


# ============================================================
# CANDIDATE RESIDUES
# ============================================================

CANDIDATES = {
    57: "A",
    60: "D",
    63: "N",
    232: "A",
    235: "M",
    337: "Q",
}


# ============================================================
# METAL-COORDINATING ATOM TYPES
# ============================================================

# Atoms that can potentially participate in metal coordination.
# This is a screening rule only.

DONOR_ATOMS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "ASN": ["OD1", "ND2"],
    "GLN": ["OE1", "NE2"],
    "HIS": ["ND1", "NE2"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "TYR": ["OH"],
    "CYS": ["SG"],
    "MET": ["SD"],
}


# ============================================================
# DISTANCE FUNCTION
# ============================================================

def distance(coord1, coord2):

    return math.sqrt(
        (coord1[0] - coord2[0]) ** 2
        + (coord1[1] - coord2[1]) ** 2
        + (coord1[2] - coord2[2]) ** 2
    )


# ============================================================
# LOAD STRUCTURE
# ============================================================

def load_structure():

    parser = PDBParser(QUIET=True)

    return parser.get_structure(
        "OsNRAMP5",
        str(PDB_FILE)
    )


# ============================================================
# GET CANDIDATE RESIDUES
# ============================================================

def get_candidates(structure):

    model = structure[0]
    chain = model["A"]

    candidates = []

    for position, expected in CANDIDATES.items():

        residue = chain[(" ", position, " ")]

        plddt_values = [
            float(atom.bfactor)
            for atom in residue
        ]

        mean_plddt = (
            sum(plddt_values) / len(plddt_values)
            if plddt_values
            else None
        )

        donor_atoms = []

        for atom_name in DONOR_ATOMS.get(
            residue.resname,
            []
        ):

            if atom_name in residue:

                donor_atoms.append(atom_name)

        candidates.append(
            {
                "position": position,
                "expected_residue": expected,
                "residue_3letter": residue.resname,
                "residue": residue,
                "mean_plddt": (
                    round(mean_plddt, 2)
                    if mean_plddt is not None
                    else None
                ),
                "potential_donor_atoms":
                    donor_atoms,
            }
        )

    return candidates


# ============================================================
# FIND NEARBY PROTEIN RESIDUES
# ============================================================

def find_nearby_residues(
    structure,
    candidate,
    cutoff=6.0
):

    model = structure[0]
    chain = model["A"]

    candidate_residue = candidate["residue"]

    nearby = {}

    for residue in chain:

        if residue.id[0] != " ":
            continue

        position = residue.id[1]

        if position == candidate["position"]:
            continue

        minimum_distance = float("inf")
        closest_atom_1 = None
        closest_atom_2 = None

        for atom1 in candidate_residue:

            for atom2 in residue:

                d = distance(
                    atom1.coord,
                    atom2.coord
                )

                if d < minimum_distance:

                    minimum_distance = d
                    closest_atom_1 = atom1.name
                    closest_atom_2 = atom2.name

        if minimum_distance <= cutoff:

            key = position

            nearby[key] = {
                "position": position,
                "residue_3letter": residue.resname,
                "distance_angstrom": round(
                    minimum_distance,
                    3
                ),
                "closest_candidate_atom":
                    closest_atom_1,
                "closest_neighbor_atom":
                    closest_atom_2,
                "potential_donor_atoms":
                    DONOR_ATOMS.get(
                        residue.resname,
                        []
                    ),
            }

    return sorted(
        nearby.values(),
        key=lambda x: x["distance_angstrom"]
    )


# ============================================================
# DONOR ATOM ANALYSIS
# ============================================================

def analyze_donor_environment(
    structure,
    candidate,
    cutoff=5.0
):

    model = structure[0]
    chain = model["A"]

    candidate_residue = candidate["residue"]

    donor_contacts = []

    for atom_name in candidate[
        "potential_donor_atoms"
    ]:

        if atom_name not in candidate_residue:
            continue

        candidate_atom = candidate_residue[
            atom_name
        ]

        for residue in chain:

            if residue.id[0] != " ":
                continue

            if residue.id[1] == candidate["position"]:
                continue

            donor_names = DONOR_ATOMS.get(
                residue.resname,
                []
            )

            for neighbor_atom_name in donor_names:

                if neighbor_atom_name not in residue:
                    continue

                neighbor_atom = residue[
                    neighbor_atom_name
                ]

                d = distance(
                    candidate_atom.coord,
                    neighbor_atom.coord
                )

                if d <= cutoff:

                    donor_contacts.append(
                        {
                            "candidate_atom":
                                atom_name,
                            "neighbor_position":
                                residue.id[1],
                            "neighbor_residue":
                                residue.resname,
                            "neighbor_atom":
                                neighbor_atom_name,
                            "distance_angstrom":
                                round(d, 3),
                        }
                    )

    return sorted(
        donor_contacts,
        key=lambda x: x[
            "distance_angstrom"
        ]
    )


# ============================================================
# POCKET CENTER
# ============================================================

def calculate_pocket_center(candidates):

    coordinates = []

    for candidate in candidates:

        residue = candidate["residue"]

        if "CA" in residue:
            coord = residue["CA"].coord

            coordinates.append([
                float(coord[0]),
                float(coord[1]),
                float(coord[2])
            ])

    if not coordinates:
        return None

    x = sum(
        coord[0]
        for coord in coordinates
    ) / len(coordinates)

    y = sum(
        coord[1]
        for coord in coordinates
    ) / len(coordinates)

    z = sum(
        coord[2]
        for coord in coordinates
    ) / len(coordinates)

    return [
        round(float(x), 3),
        round(float(y), 3),
        round(float(z), 3)
    ]


# ============================================================
# CANDIDATE DISTANCE FROM POCKET CENTER
# ============================================================

def distance_from_center(
    candidate,
    center
):

    if center is None:
        return None

    residue = candidate["residue"]

    if "CA" not in residue:
        return None

    d = distance(
        residue["CA"].coord,
        center
    )

    return round(d, 3)


# ============================================================
# CHEMICAL CONTRIBUTION
# ============================================================

def chemical_score(candidate):

    residue = candidate[
        "residue_3letter"
    ]

    donor_atoms = candidate[
        "potential_donor_atoms"
    ]

    score = 0

    # Strongest group for metal coordination
    if residue in {
        "ASP",
        "GLU",
        "HIS",
        "CYS"
    }:
        score += 30

    # Additional possible donor residues
    elif residue in {
        "ASN",
        "GLN",
        "SER",
        "THR",
        "TYR",
        "MET"
    }:
        score += 20

    else:
        score += 0

    if donor_atoms:
        score += 10

    return min(score, 40)


# ============================================================
# STRUCTURAL SCORE
# ============================================================

def structural_score(
    candidate,
    nearby_residues
):

    score = 0

    # Local density
    number_nearby = len(
        nearby_residues
    )

    if number_nearby >= 10:
        score += 30

    elif number_nearby >= 6:
        score += 25

    elif number_nearby >= 3:
        score += 20

    elif number_nearby >= 1:
        score += 10

    # AlphaFold confidence
    plddt = candidate[
        "mean_plddt"
    ]

    if plddt is not None:

        if plddt >= 90:
            score += 30

        elif plddt >= 80:
            score += 25

        elif plddt >= 70:
            score += 15

        elif plddt >= 50:
            score += 5

    return min(score, 60)


# ============================================================
# POCKET SCORE
# ============================================================

def calculate_pocket_score(
    candidate,
    nearby_residues
):

    chemistry = chemical_score(
        candidate
    )

    structure = structural_score(
        candidate,
        nearby_residues
    )

    total = chemistry + structure

    if total >= 75:
        evidence = "high"

    elif total >= 50:
        evidence = "moderate"

    elif total >= 25:
        evidence = "low"

    else:
        evidence = "very_low"

    return {
        "chemical_score": chemistry,
        "structural_score": structure,
        "total_score": total,
        "evidence_level": evidence,
    }


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "Loading OsNRAMP5 AlphaFold structure..."
    )

    print(PDB_FILE)
    print()

    if not PDB_FILE.exists():

        raise FileNotFoundError(
            f"Structure not found:\n{PDB_FILE}"
        )

    structure = load_structure()

    print(
        "Structure loaded successfully."
    )
    print()

    candidates = get_candidates(
        structure
    )

    print(
        f"Candidate residues: "
        f"{len(candidates)}"
    )

    print()

    # --------------------------------------------------------
    # Pocket center
    # --------------------------------------------------------

    pocket_center = calculate_pocket_center(
        candidates
    )

    print(
        "Candidate pocket center:"
    )

    print(pocket_center)
    print()

    # --------------------------------------------------------
    # Analyze candidates
    # --------------------------------------------------------

    candidate_results = []

    for candidate in candidates:

        position = candidate[
            "position"
        ]

        residue = candidate[
            "residue_3letter"
        ]

        print(
            f"Analyzing {residue}{position}..."
        )

        nearby = find_nearby_residues(
            structure,
            candidate,
            cutoff=6.0
        )

        donor_contacts = (
            analyze_donor_environment(
                structure,
                candidate,
                cutoff=5.0
            )
        )

        pocket_score = (
            calculate_pocket_score(
                candidate,
                nearby
            )
        )

        center_distance = (
            distance_from_center(
                candidate,
                pocket_center
            )
        )

        print(
            f"  pLDDT: "
            f"{candidate['mean_plddt']}"
        )

        print(
            f"  Nearby residues: "
            f"{len(nearby)}"
        )

        print(
            f"  Potential donor atoms: "
            f"{candidate['potential_donor_atoms']}"
        )

        print(
            f"  Donor contacts: "
            f"{len(donor_contacts)}"
        )

        print(
            f"  Pocket score: "
            f"{pocket_score['total_score']}/100"
        )

        print(
            f"  Evidence: "
            f"{pocket_score['evidence_level']}"
        )

        print()

        candidate_results.append(
            {
                "position": position,

                "residue": residue,

                "residue_1letter":
                    candidate[
                        "expected_residue"
                    ],

                "mean_plddt":
                    candidate[
                        "mean_plddt"
                    ],

                "potential_donor_atoms":
                    candidate[
                        "potential_donor_atoms"
                    ],

                "nearby_residues":
                    nearby,

                "potential_donor_contacts":
                    donor_contacts,

                "distance_from_candidate_center":
                    center_distance,

                "pocket_score":
                    pocket_score,
            }
        )

    # --------------------------------------------------------
    # Rank candidates
    # --------------------------------------------------------

    ranked_candidates = sorted(
        candidate_results,
        key=lambda x:
            x["pocket_score"][
                "total_score"
            ],
        reverse=True
    )

    for rank, candidate in enumerate(
        ranked_candidates,
        start=1
    ):

        candidate["rank"] = rank

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    high = sum(
        1
        for c in ranked_candidates
        if c["pocket_score"][
            "evidence_level"
        ] == "high"
    )

    moderate = sum(
        1
        for c in ranked_candidates
        if c["pocket_score"][
            "evidence_level"
        ] == "moderate"
    )

    low = sum(
        1
        for c in ranked_candidates
        if c["pocket_score"][
            "evidence_level"
        ] == "low"
    )

    very_low = sum(
        1
        for c in ranked_candidates
        if c["pocket_score"][
            "evidence_level"
        ] == "very_low"
    )

    # --------------------------------------------------------
    # Result
    # --------------------------------------------------------

    result = {

        "target": {
            "gene": "OsNRAMP5",
            "uniprot": "Q8H4H5",
            "ncbi_accession":
                "NP_001389970.1",
            "sequence_length_aa": 538
        },

        "structure": {
            "model":
                "AF-Q8H4H5-F1-model_v6",
            "source":
                "AlphaFold Protein Structure Database",
            "chain": "A",
            "file":
                "data/structures/predicted/"
                "AF-Q8H4H5-F1-model_v6.pdb"
        },

        "candidate_definition": {
            "residues": [
                "A57",
                "D60",
                "N63",
                "A232",
                "M235",
                "Q337"
            ],
            "source":
                "Mapped reference metal-associated "
                "residues and AlphaFold structural "
                "validation"
        },

        "analysis_method": {
            "name":
                "Structural metal-pocket screening",

            "nearby_residue_cutoff_angstrom":
                6.0,

            "donor_contact_cutoff_angstrom":
                5.0,

            "potential_donor_residues": [
                "ASP",
                "GLU",
                "ASN",
                "GLN",
                "HIS",
                "SER",
                "THR",
                "TYR",
                "CYS",
                "MET"
            ],

            "score_components": {
                "chemical_environment":
                    "potential metal-donor chemistry",

                "structural_environment":
                    "local residue density and AlphaFold pLDDT"
            }
        },

        "candidate_pocket_center": {
            "coordinate_system":
                "AlphaFold PDB coordinates",

            "center_xyz_angstrom":
                pocket_center
        },

        "summary": {
            "total_candidates":
                len(ranked_candidates),

            "high_evidence":
                high,

            "moderate_evidence":
                moderate,

            "low_evidence":
                low,

            "very_low_evidence":
                very_low
        },

        "ranked_candidates":
            ranked_candidates,

        "interpretation":
            (
                "Candidate residues were evaluated for "
                "their local structural environment, "
                "potential metal-donor atoms, nearby "
                "residue density, and AlphaFold confidence. "
                "Higher scores indicate stronger structural "
                "support for prioritization. This analysis "
                "does not establish direct metal coordination, "
                "metal-binding affinity, or experimental "
                "binding probability. Docking or additional "
                "structural validation is required."
            )
    }

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    OUTPUT_DIR.mkdir(
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

    print(
        "Metal pocket analysis complete."
    )

    print()

    print(
        "Ranked candidate sites:"
    )

    print()

    for candidate in ranked_candidates:

        print(
            f"Rank {candidate['rank']}: "
            f"{candidate['residue']}"
            f"{candidate['position']} "
            f"| Score "
            f"{candidate['pocket_score']['total_score']}"
            f"/100 "
            f"| "
            f"{candidate['pocket_score']['evidence_level']}"
        )

    print()

    print(
        "Results saved to:"
    )

    print(OUTPUT_FILE)


# ============================================================
# ENTRY POINT
# ============================================================

if __name__ == "__main__":
    main()
