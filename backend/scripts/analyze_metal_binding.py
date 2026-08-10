from pathlib import Path
import json
import math

from Bio.PDB import PDBParser


# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------

PDB_FILE = Path("data/structures/pdb/8E6N.pdb")

OUTPUT_DIR = Path("data/processed/structures")
OUTPUT_FILE = OUTPUT_DIR / "8E6N_metal_binding.json"

METAL_NAMES = {
    "MN",
    "MG",
    "FE",
    "ZN",
    "CA",
    "CO",
    "NI",
    "CU",
}

# Distance threshold for identifying nearby residues.
# 3.5 Å is a reasonable first-pass structural screening cutoff.
DISTANCE_CUTOFF = 3.5


# ---------------------------------------------------------
# Helper functions
# ---------------------------------------------------------

def distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms."""
    return math.sqrt(
        sum(
            (a - b) ** 2
            for a, b in zip(atom1.coord, atom2.coord)
        )
    )


def find_metal_atoms(structure):
    """Find metal hetero atoms in the structure."""
    metals = []

    for model in structure:
        for chain in model:
            for residue in chain:

                residue_name = residue.get_resname().strip().upper()

                if residue_name in METAL_NAMES:
                    for atom in residue:
                        metals.append({
                            "model": model.id,
                            "chain": chain.id,
                            "residue_number": residue.id[1],
                            "residue_name": residue_name,
                            "atom_name": atom.get_name(),
                            "coord": [
                                float(atom.coord[0]),
                                float(atom.coord[1]),
                                float(atom.coord[2]),
                            ],
                        })

    return metals


def find_nearby_residues(structure, metal_atom):
    """
    Find protein residues with at least one atom
    within DISTANCE_CUTOFF Å of the metal atom.
    """

    nearby = []

    metal_coord = metal_atom["coord"]

    for model in structure:
        for chain in model:

            for residue in chain:

                # Ignore hetero residues.
                if residue.id[0] != " ":
                    continue

                min_distance = None
                closest_atom = None
                closest_residue_atom = None

                for atom in residue:

                    d = math.sqrt(
                        sum(
                            (a - b) ** 2
                            for a, b in zip(
                                atom.coord,
                                metal_coord
                            )
                        )
                    )

                    if min_distance is None or d < min_distance:
                        min_distance = d
                        closest_atom = atom.get_name()
                        closest_residue_atom = atom

                if min_distance is not None and min_distance <= DISTANCE_CUTOFF:

                    nearby.append({
                        "chain": chain.id,
                        "residue_number": residue.id[1],
                        "residue_name": residue.get_resname().strip(),
                        "closest_atom": closest_atom,
                        "distance_angstrom": round(
                            float(min_distance),
                            3
                        ),
                    })

    nearby.sort(
        key=lambda x: x["distance_angstrom"]
    )

    return nearby


# ---------------------------------------------------------
# Main analysis
# ---------------------------------------------------------

def main():

    print("Loading PDB structure...")
    print(PDB_FILE)

    if not PDB_FILE.exists():
        raise FileNotFoundError(
            f"PDB file not found: {PDB_FILE}"
        )

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "NRAMP_8E6N",
        PDB_FILE
    )

    print("Structure loaded successfully.")

    metals = find_metal_atoms(structure)

    print()
    print(f"Metal atoms detected: {len(metals)}")

    if not metals:
        print("No metal atoms detected.")
        return

    results = []

    for index, metal in enumerate(metals, start=1):

        print()
        print(
            f"Analyzing metal {index}: "
            f"{metal['residue_name']} "
            f"chain {metal['chain']} "
            f"residue {metal['residue_number']}"
        )

        nearby = find_nearby_residues(
            structure,
            metal
        )

        print(
            f"Nearby protein residues "
            f"within {DISTANCE_CUTOFF} Å: "
            f"{len(nearby)}"
        )

        for residue in nearby:

            print(
                f"  {residue['residue_name']}"
                f"{residue['residue_number']} "
                f"chain {residue['chain']} "
                f"{residue['distance_angstrom']} Å "
                f"via {residue['closest_atom']}"
            )

        results.append({
            "metal": metal,
            "distance_cutoff_angstrom": DISTANCE_CUTOFF,
            "nearby_residues": nearby,
        })

    # -----------------------------------------------------
    # Save results
    # -----------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    output = {
        "structure": {
            "pdb_id": "8E6N",
            "file": str(PDB_FILE),
            "description": (
                "NRAMP/MntH divalent metal transporter "
                "reference structure"
            ),
        },
        "analysis": {
            "method": (
                "Distance-based screening of protein "
                "residues surrounding detected metal atoms"
            ),
            "distance_cutoff_angstrom": DISTANCE_CUTOFF,
        },
        "metals_detected": len(metals),
        "metal_binding_environments": results,
        "interpretation": (
            "Residues listed here are structurally proximal "
            "to detected metal atoms within the selected "
            "distance cutoff. Proximity alone does not prove "
            "direct metal coordination."
        ),
    }

    with open(
        OUTPUT_FILE,
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            output,
            handle,
            indent=4
        )

    print()
    print("Metal-binding analysis complete.")
    print()
    print("Results saved to:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
