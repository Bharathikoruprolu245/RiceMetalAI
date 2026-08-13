import json
import math
from pathlib import Path

from Bio.PDB import PDBParser


# ============================================================
# INPUT / OUTPUT
# ============================================================

PDB_FILE = Path(
    "data/structures/predicted/AF-Q8H4H5-F1-model_v6.pdb"
)

OUTPUT_FILE = Path(
    "data/processed/structures/OsNRAMP5_focused_metal_pocket.json"
)


# ============================================================
# TARGET CANDIDATES
# ============================================================

# Primary donor-capable candidates identified from
# the previous pocket analysis.
PRIMARY_CANDIDATES = {
    60: {
        "residue": "ASP",
        "expected": "D",
        "donor_atoms": ["OD1", "OD2"]
    },
    63: {
        "residue": "ASN",
        "expected": "N",
        "donor_atoms": ["OD1", "ND2"]
    },
    235: {
        "residue": "MET",
        "expected": "M",
        "donor_atoms": ["SD"]
    },
    337: {
        "residue": "GLN",
        "expected": "Q",
        "donor_atoms": ["OE1", "NE2"]
    }
}


# Supporting residues identified in the broad pocket
SUPPORTING_RESIDUES = {
    57: "ALA",
    232: "ALA"
}


# Distance used for donor-atom neighborhood screening
DONOR_NEIGHBOR_CUTOFF = 6.0


# ============================================================
# HELPERS
# ============================================================

def atom_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms."""

    c1 = atom1.coord
    c2 = atom2.coord

    dx = float(c1[0]) - float(c2[0])
    dy = float(c1[1]) - float(c2[1])
    dz = float(c1[2]) - float(c2[2])

    return math.sqrt(
        dx * dx +
        dy * dy +
        dz * dz
    )


def get_residue(chain, position):

    for residue in chain:

        if residue.id[0] != " ":
            continue

        if residue.id[1] == position:
            return residue

    return None


def get_atom(residue, atom_name):

    if atom_name in residue:
        return residue[atom_name]

    return None


def get_plddt(residue):

    values = []

    for atom in residue.get_atoms():
        values.append(
            float(atom.get_bfactor())
        )

    if not values:
        return None

    return round(
        sum(values) / len(values),
        2
    )


def residue_center(residue):

    coords = []

    for atom in residue.get_atoms():

        coord = atom.coord

        coords.append([
            float(coord[0]),
            float(coord[1]),
            float(coord[2])
        ])

    if not coords:
        return None

    x = sum(c[0] for c in coords) / len(coords)
    y = sum(c[1] for c in coords) / len(coords)
    z = sum(c[2] for c in coords) / len(coords)

    return [
        round(float(x), 3),
        round(float(y), 3),
        round(float(z), 3)
    ]


def centroid(points):

    if not points:
        return None

    x = sum(p[0] for p in points) / len(points)
    y = sum(p[1] for p in points) / len(points)
    z = sum(p[2] for p in points) / len(points)

    return [
        round(float(x), 3),
        round(float(y), 3),
        round(float(z), 3)
    ]


def distance_from_point(point, atom):

    coord = atom.coord

    dx = float(point[0]) - float(coord[0])
    dy = float(point[1]) - float(coord[1])
    dz = float(point[2]) - float(coord[2])

    return math.sqrt(
        dx * dx +
        dy * dy +
        dz * dz
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print(
        "Loading OsNRAMP5 AlphaFold structure..."
    )

    print(PDB_FILE)

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "OsNRAMP5",
        str(PDB_FILE)
    )

    model = structure[0]
    chain = model["A"]

    print()
    print(
        "Structure loaded successfully."
    )

    print()
    print(
        "Primary donor candidates:"
    )

    for position, info in PRIMARY_CANDIDATES.items():

        print(
            f"  {info['residue']}{position} "
            f"({', '.join(info['donor_atoms'])})"
        )

    # --------------------------------------------------------
    # Collect donor atoms
    # --------------------------------------------------------

    donor_atoms = []

    candidate_results = []

    for position, info in PRIMARY_CANDIDATES.items():

        residue = get_residue(
            chain,
            position
        )

        if residue is None:

            print(
                f"WARNING: residue {position} "
                f"not found."
            )

            continue

        actual_residue = (
            residue.get_resname()
            .strip()
        )

        plddt = get_plddt(
            residue
        )

        center = residue_center(
            residue
        )

        residue_donors = []

        for atom_name in info["donor_atoms"]:

            atom = get_atom(
                residue,
                atom_name
            )

            if atom is None:
                continue

            coord = [
                float(atom.coord[0]),
                float(atom.coord[1]),
                float(atom.coord[2])
            ]

            residue_donors.append({
                "atom": atom_name,
                "coordinates": [
                    round(coord[0], 3),
                    round(coord[1], 3),
                    round(coord[2], 3)
                ]
            })

            donor_atoms.append({
                "position": position,
                "residue": actual_residue,
                "atom": atom_name,
                "atom_object": atom,
                "coordinates": coord
            })

        candidate_results.append({
            "position": position,
            "expected_residue": info["expected"],
            "structure_residue": actual_residue,
            "plddt": plddt,
            "residue_center": center,
            "potential_donor_atoms":
                residue_donors
        })

    # --------------------------------------------------------
    # Donor atom centroid
    # --------------------------------------------------------

    donor_points = [
        item["coordinates"]
        for item in donor_atoms
    ]

    donor_centroid = centroid(
        donor_points
    )

    print()
    print(
        "Potential donor atoms detected:",
        len(donor_atoms)
    )

    print()
    print(
        "Donor-atom centroid:"
    )

    print(
        donor_centroid
    )

    # --------------------------------------------------------
    # Distance of each donor from centroid
    # --------------------------------------------------------

    donor_geometry = []

    for item in donor_atoms:

        d = distance_from_point(
            donor_centroid,
            item["atom_object"]
        )

        donor_geometry.append({
            "residue_position":
                int(item["position"]),

            "residue":
                item["residue"],

            "atom":
                item["atom"],

            "distance_from_donor_centroid_angstrom":
                round(float(d), 3)
        })

    donor_geometry.sort(
        key=lambda x:
        x["distance_from_donor_centroid_angstrom"]
    )

    # --------------------------------------------------------
    # Pairwise donor distances
    # --------------------------------------------------------

    donor_pairwise = []

    for i in range(
        len(donor_atoms)
    ):

        for j in range(
            i + 1,
            len(donor_atoms)
        ):

            item1 = donor_atoms[i]
            item2 = donor_atoms[j]

            d = atom_distance(
                item1["atom_object"],
                item2["atom_object"]
            )

            donor_pairwise.append({

                "residue_1":
                    int(item1["position"]),

                "atom_1":
                    item1["atom"],

                "residue_2":
                    int(item2["position"]),

                "atom_2":
                    item2["atom"],

                "distance_angstrom":
                    round(float(d), 3)
            })

    donor_pairwise.sort(
        key=lambda x:
        x["distance_angstrom"]
    )

    # --------------------------------------------------------
    # Supporting residues
    # --------------------------------------------------------

    supporting_results = []

    for position, expected in SUPPORTING_RESIDUES.items():

        residue = get_residue(
            chain,
            position
        )

        if residue is None:
            continue

        supporting_results.append({

            "position":
                int(position),

            "residue":
                residue.get_resname().strip(),

            "expected_residue":
                expected,

            "plddt":
                get_plddt(residue),

            "residue_center":
                residue_center(residue)
        })

    # --------------------------------------------------------
    # Donor neighborhood
    # --------------------------------------------------------

    donor_neighborhood = []

    all_residues = [
        residue
        for residue in chain
        if residue.id[0] == " "
    ]

    for donor in donor_atoms:

        nearby = []

        for residue in all_residues:

            if residue.id[1] == donor["position"]:
                continue

            if "CA" not in residue:
                continue

            ca = residue["CA"]

            d = atom_distance(
                donor["atom_object"],
                ca
            )

            if d <= DONOR_NEIGHBOR_CUTOFF:

                nearby.append({

                    "position":
                        int(residue.id[1]),

                    "residue":
                        residue.get_resname().strip(),

                    "distance_angstrom":
                        round(float(d), 3)
                })

        nearby.sort(
            key=lambda x:
            x["distance_angstrom"]
        )

        donor_neighborhood.append({

            "donor_residue":
                int(donor["position"]),

            "donor_atom":
                donor["atom"],

            "nearby_residues":
                nearby
        })

    # --------------------------------------------------------
    # Basic structural interpretation
    # --------------------------------------------------------

    donor_residue_positions = sorted(
        set(
            item["position"]
            for item in donor_atoms
        )
    )

    if len(donor_residue_positions) >= 3:

        interpretation = (
            "The focused OsNRAMP5 pocket contains multiple "
            "candidate donor-bearing residues with measurable "
            "three-dimensional proximity. D60, N63, M235, "
            "and Q337 should therefore be prioritized for "
            "explicit metal interaction modeling. The donor "
            "centroid is only a geometric reference and is "
            "not an experimentally determined metal position."
        )

    else:

        interpretation = (
            "The focused pocket contains limited donor-residue "
            "support and requires additional structural analysis."
        )

    # --------------------------------------------------------
    # Save JSON
    # --------------------------------------------------------

    result = {

        "target": {

            "gene": "OsNRAMP5",

            "uniprot": "Q8H4H5",

            "ncbi_accession":
                "NP_001389970.1",

            "structure":
                "AF-Q8H4H5-F1-model_v6",

            "chain": "A",

            "sequence_length_aa": 538
        },

        "analysis": {

            "name":
                "Focused metal-binding pocket geometry",

            "primary_candidates":
                list(
                    PRIMARY_CANDIDATES.keys()
                ),

            "supporting_residues":
                list(
                    SUPPORTING_RESIDUES.keys()
                ),

            "donor_neighbor_cutoff_angstrom":
                DONOR_NEIGHBOR_CUTOFF
        },

        "primary_candidate_residues":
            candidate_results,

        "potential_donor_atoms":
            [
                {
                    "position":
                        int(item["position"]),

                    "residue":
                        item["residue"],

                    "atom":
                        item["atom"],

                    "coordinates":
                        [
                            round(
                                float(item["coordinates"][0]),
                                3
                            ),
                            round(
                                float(item["coordinates"][1]),
                                3
                            ),
                            round(
                                float(item["coordinates"][2]),
                                3
                            )
                        ]
                }

                for item in donor_atoms
            ],

        "donor_atom_centroid":
            donor_centroid,

        "donor_geometry":
            donor_geometry,

        "pairwise_donor_distances":
            donor_pairwise,

        "supporting_residues":
            supporting_results,

        "donor_neighborhood":
            donor_neighborhood,

        "interpretation":
            interpretation
    }

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        OUTPUT_FILE,
        "w"
    ) as handle:

        json.dump(
            result,
            handle,
            indent=4
        )

    # --------------------------------------------------------
    # Console summary
    # --------------------------------------------------------

    print()
    print(
        "Focused metal-binding pocket analysis complete."
    )

    print()
    print(
        "Donor atoms:"
    )

    for item in donor_atoms:

        print(
            f"  {item['residue']}"
            f"{item['position']} - "
            f"{item['atom']}"
        )

    print()
    print(
        "Closest donor-atom relationships:"
    )

    for item in donor_pairwise[:10]:

        print(
            f"  {item['residue_1']}"
            f"{item['atom_1']} -> "
            f"{item['residue_2']}"
            f"{item['atom_2']} : "
            f"{item['distance_angstrom']} Å"
        )

    print()
    print(
        "Results saved to:"
    )

    print(
        OUTPUT_FILE.resolve()
    )


if __name__ == "__main__":
    main()
