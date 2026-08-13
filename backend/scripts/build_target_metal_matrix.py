#!/usr/bin/env python3

"""
RiceMetalAI
Build the scalable Target × Metal analysis matrix.
"""

import json
import sys
from pathlib import Path

sys.path.insert(
    0,
    str(Path(__file__).resolve().parent)
)

from load_target_config import load_targets
from load_metal_config import load_metals


PROJECT_ROOT = Path(__file__).resolve().parents[2]

OUTPUT_FILE = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "target_metal_matrix.json"
)


def build_matrix():

    targets = load_targets()
    metals = load_metals()

    matrix = []

    for target in targets:

        for metal_symbol in target["metals"]:

            if metal_symbol not in metals:
                raise ValueError(
                    f"Metal '{metal_symbol}' "
                    f"configured for "
                    f"{target['gene']} but not "
                    f"defined in metals.json."
                )

            metal = metals[metal_symbol]

            matrix.append(
                {
                    "gene": target["gene"],
                    "protein": target["protein"],
                    "ncbi_accession":
                        target["ncbi_accession"],
                    "uniprot_accession":
                        target["uniprot_accession"],
                    "sequence_length":
                        target["sequence_length"],
                    "target_status":
                        target["status"],

                    "metal": {
                        "symbol":
                            metal["symbol"],
                        "name":
                            metal["name"],
                        "atomic_number":
                            metal["atomic_number"],
                        "primary_ion":
                            metal["primary_ion"],
                        "oxidation_states":
                            metal[
                                "common_oxidation_states"
                            ],
                        "category":
                            metal["category"],
                        "essential_element":
                            metal[
                                "essential_element"
                            ]
                    },

                    "analysis_status":
                        "configured"
                }
            )

    return matrix


def main():

    print("=" * 70)
    print(
        "RiceMetalAI - Target × Metal Matrix"
    )
    print("=" * 70)
    print()

    matrix = build_matrix()

    print(
        f"Target-metal combinations: "
        f"{len(matrix)}"
    )
    print()

    for item in matrix:

        print(
            f"{item['gene']:12s} × "
            f"{item['metal']['symbol']:2s} "
            f"({item['metal']['name']})"
        )

    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result = {
        "project": "RiceMetalAI",

        "description": (
            "Scalable matrix connecting "
            "configured rice target proteins "
            "with their configured metal species."
        ),

        "total_combinations":
            len(matrix),

        "combinations":
            matrix
    }

    with OUTPUT_FILE.open(
        "w",
        encoding="utf-8"
    ) as handle:

        json.dump(
            result,
            handle,
            indent=4
        )

    print()
    print(
        "Matrix saved to:"
    )
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
