from pathlib import Path
from Bio.PDB import PDBParser


PDB_FILE = Path(
    "data/structures/pdb/8E6N.pdb"
)


def main():

    print("Loading PDB structure...")

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "8E6N",
        PDB_FILE
    )

    models = list(structure)

    print(
        f"Models: {len(models)}"
    )

    for model in models:

        chains = list(model)

        print(
            f"Model {model.id}: "
            f"{len(chains)} chain(s)"
        )

        for chain in chains:

            residues = list(chain)

            print(
                f"  Chain {chain.id}: "
                f"{len(residues)} residues"
            )

    # -----------------------------------------------------
    # Count atoms and hetero atoms
    # -----------------------------------------------------

    atom_count = 0
    hetero_residues = []

    for model in structure:

        for chain in model:

            for residue in chain:

                atom_count += len(
                    list(residue)
                )

                hetflag = residue.id[0]

                if hetflag != " ":

                    hetero_residues.append(
                        (
                            chain.id,
                            residue.id,
                            residue.resname
                        )
                    )

    print()
    print(
        f"Total atoms: {atom_count}"
    )

    print(
        f"Hetero residues: "
        f"{len(hetero_residues)}"
    )

    # -----------------------------------------------------
    # Search for manganese
    # -----------------------------------------------------

    manganese = []

    for item in hetero_residues:

        chain_id, residue_id, resname = item

        if resname.upper() in {
            "MN",
            "MNG"
        }:

            manganese.append(item)

    print()

    if manganese:

        print(
            "Manganese detected:"
        )

        for item in manganese:

            print(
                f"  Chain {item[0]}, "
                f"Residue {item[1]}, "
                f"Name {item[2]}"
            )

    else:

        print(
            "WARNING: Manganese was "
            "not detected."
        )

    print()
    print(
        "PDB validation complete."
    )


if __name__ == "__main__":
    main()
