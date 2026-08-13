#!/usr/bin/env python3

"""
Analyze the 3D geometry of candidate metal-binding residues
in the OsNRAMP5 AlphaFold structure.

Input:
    data/structures/predicted/AF-Q8H4H5-F1-model_v6.pdb

Candidates:
    A57, D60, N63, A232, M235, Q337

Outputs:
    data/processed/structures/OsNRAMP5_candidate_3d_geometry.json
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
    / "OsNRAMP5_candidate_3d_geometry.json"
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
# DISTANCE FUNCTION
# ============================================================

def calculate_distance(coord1, coord2):
    """
    Calculate Euclidean distance between two 3D coordinates.
    """

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

    structure = parser.get_structure(
        "OsNRAMP5",
        str(PDB_FILE)
    )

    return structure


# ============================================================
# EXTRACT CANDIDATES
# ============================================================

def get_candidate_residues(structure):

    model = structure[0]

    if "A" not in model:

        raise ValueError(
            "Chain A was not found in the AlphaFold structure."
        )

    chain = model["A"]

    candidate_residues = []

    for position, expected_residue in CANDIDATES.items():

        residue = chain[(" ", position, " ")]

        actual_residue = residue.resname

        candidate_residues.append(
            {
                "position": position,
                "expected_residue": expected_residue,
                "residue_3letter": actual_residue,
                "residue": residue,
            }
        )

    return candidate_residues


# ============================================================
# CA DISTANCE
# ============================================================

def calculate_ca_distances(candidates):

    results = []

    for i in range(len(candidates)):

        residue_a = candidates[i]["residue"]

        if "CA" not in residue_a:
            continue

        for j in range(i + 1, len(candidates)):

            residue_b = candidates[j]["residue"]

            if "CA" not in residue_b:
                continue

            distance = calculate_distance(
                residue_a["CA"].coord,
                residue_b["CA"].coord
            )

            results.append(
                {
                    "residue_1": candidates[i]["position"],
                    "residue_2": candidates[j]["position"],
                    "distance_angstrom": round(distance, 3)
                }
            )

    return results


# ============================================================
# CLOSEST ATOM DISTANCE
# ============================================================

def calculate_closest_atom_distances(candidates):

    results = []

    for i in range(len(candidates)):

        residue_a = candidates[i]["residue"]

        for j in range(i + 1, len(candidates)):

            residue_b = candidates[j]["residue"]

            minimum_distance = float("inf")

            atom_a_name = None
            atom_b_name = None

            for atom_a in residue_a:

                for atom_b in residue_b:

                    distance = calculate_distance(
                        atom_a.coord,
                        atom_b.coord
                    )

                    if distance < minimum_distance:

                        minimum_distance = distance
                        atom_a_name = atom_a.name
                        atom_b_name = atom_b.name

            results.append(
                {
                    "residue_1": candidates[i]["position"],
                    "residue_2": candidates[j]["position"],
                    "closest_atom_1": atom_a_name,
                    "closest_atom_2": atom_b_name,
                    "distance_angstrom": round(
                        minimum_distance,
                        3
                    )
                }
            )

    return results


# ============================================================
# IDENTIFY SPATIAL CLUSTERS
# ============================================================

def identify_clusters(candidates, ca_distances):

    """
    A simple structural clustering rule.

    Residues are considered spatially associated when their
    C-alpha atoms are <= 15 Å apart.

    This is a screening criterion, not a formal pocket detector.
    """

    threshold = 15.0

    graph = {
        candidate["position"]: set()
        for candidate in candidates
    }

    for pair in ca_distances:

        if pair["distance_angstrom"] <= threshold:

            residue_1 = pair["residue_1"]
            residue_2 = pair["residue_2"]

            graph[residue_1].add(residue_2)
            graph[residue_2].add(residue_1)

    clusters = []

    visited = set()

    for position in graph:

        if position in visited:
            continue

        cluster = set()

        stack = [position]

        while stack:

            current = stack.pop()

            if current in visited:
                continue

            visited.add(current)

            cluster.add(current)

            for neighbor in graph[current]:

                if neighbor not in visited:

                    stack.append(neighbor)

        clusters.append(
            sorted(cluster)
        )

    return clusters


# ============================================================
# LOCAL PLDDT
# ============================================================

def calculate_residue_plddt(residue):

    values = []

    for atom in residue:

        if hasattr(atom, "bfactor"):

            values.append(float(atom.bfactor))

    if not values:

        return None

    return round(
        sum(values) / len(values),
        2
    )


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    print("Loading OsNRAMP5 AlphaFold structure...")
    print(PDB_FILE)
    print()

    if not PDB_FILE.exists():

        raise FileNotFoundError(
            f"Structure file not found:\n{PDB_FILE}"
        )

    structure = load_structure()

    print("Structure loaded successfully.")
    print()

    model = structure[0]

    print(
        f"Model: {model.id}"
    )

    print(
        f"Chains: {len(list(model.get_chains()))}"
    )

    print()

    candidates = get_candidate_residues(
        structure
    )

    print("Candidate residues:")
    print()

    candidate_summary = []

    for candidate in candidates:

        residue = candidate["residue"]

        plddt = calculate_residue_plddt(
            residue
        )

        print(
            f"{candidate['expected_residue']}"
            f"{candidate['position']} "
            f"({candidate['residue_3letter']}) "
            f"| pLDDT {plddt}"
        )

        candidate_summary.append(
            {
                "position": candidate["position"],
                "expected_residue": candidate[
                    "expected_residue"
                ],
                "residue_3letter": candidate[
                    "residue_3letter"
                ],
                "plddt": plddt
            }
        )

    print()

    # --------------------------------------------------------
    # C-alpha distances
    # --------------------------------------------------------

    print("Calculating C-alpha distances...")

    ca_distances = calculate_ca_distances(
        candidates
    )

    print(
        f"Pairwise comparisons: "
        f"{len(ca_distances)}"
    )

    print()

    for pair in ca_distances:

        print(
            f"{pair['residue_1']}"
            f" -> "
            f"{pair['residue_2']}"
            f" : "
            f"{pair['distance_angstrom']}"
            f" Å"
        )

    print()

    # --------------------------------------------------------
    # Closest atom distances
    # --------------------------------------------------------

    print(
        "Calculating closest-atom distances..."
    )

    closest_atom_distances = (
        calculate_closest_atom_distances(
            candidates
        )
    )

    print()

    for pair in closest_atom_distances:

        print(
            f"{pair['residue_1']}"
            f" -> "
            f"{pair['residue_2']}"
            f" : "
            f"{pair['distance_angstrom']}"
            f" Å "
            f"("
            f"{pair['closest_atom_1']}"
            f" - "
            f"{pair['closest_atom_2']}"
            f")"
        )

    print()

    # --------------------------------------------------------
    # Spatial clustering
    # --------------------------------------------------------

    print("Identifying spatial clusters...")

    clusters = identify_clusters(
        candidates,
        ca_distances
    )

    print()

    for index, cluster in enumerate(
        clusters,
        start=1
    ):

        print(
            f"Cluster {index}: "
            + ", ".join(
                f"{CANDIDATES[position]}{position}"
                for position in cluster
            )
        )

    print()

    # --------------------------------------------------------
    # Nearest candidate relationships
    # --------------------------------------------------------

    nearest_neighbors = {}

    for position in CANDIDATES:

        related_pairs = [
            pair
            for pair in ca_distances
            if (
                pair["residue_1"] == position
                or pair["residue_2"] == position
            )
        ]

        if not related_pairs:

            continue

        nearest = min(
            related_pairs,
            key=lambda x: x[
                "distance_angstrom"
            ]
        )

        if nearest["residue_1"] == position:

            neighbor = nearest["residue_2"]

        else:

            neighbor = nearest["residue_1"]

        nearest_neighbors[str(position)] = {
            "nearest_candidate_position":
                neighbor,
            "distance_angstrom":
                nearest["distance_angstrom"]
        }

    # --------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------

    interpretation = (
        "Candidate residues were analyzed in the "
        "OsNRAMP5 AlphaFold structure using pairwise "
        "C-alpha and closest-atom distances. "
        "Spatially close residues may form a structural "
        "cluster and can be prioritized for further "
        "metal-binding pocket analysis. "
        "The 15 Å clustering threshold is a computational "
        "screening criterion and does not establish a "
        "metal-binding site or direct coordination."
    )

    # --------------------------------------------------------
    # Final result
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
            "chain": "A",
            "source":
                "AlphaFold Protein Structure Database",
            "file":
                "data/structures/predicted/"
                "AF-Q8H4H5-F1-model_v6.pdb"
        },

        "candidate_residues":
            candidate_summary,

        "distance_analysis": {

            "ca_distance_method":
                "C-alpha Euclidean distance",

            "closest_atom_method":
                "minimum heavy/protein atom "
                "distance between candidate residues",

            "pairwise_ca_distances":
                ca_distances,

            "pairwise_closest_atom_distances":
                closest_atom_distances
        },

        "spatial_clustering": {

            "threshold_angstrom": 15.0,

            "clusters": [
                {
                    "cluster_id": index + 1,
                    "residues": [
                        {
                            "position": position,
                            "residue":
                                CANDIDATES[position]
                        }
                        for position in cluster
                    ]
                }
                for index, cluster
                in enumerate(clusters)
            ]
        },

        "nearest_candidate_relationships":
            nearest_neighbors,

        "interpretation":
            interpretation
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
        "3D candidate geometry analysis complete."
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
