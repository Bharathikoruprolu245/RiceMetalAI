#!/usr/bin/env python3

"""
RiceMetalAI - Target Configuration Loader

Loads the master rice metal target configuration and provides
a reusable representation for downstream pipeline components.
"""

import csv
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = PROJECT_ROOT / "data" / "config" / "rice_metal_targets.csv"


REQUIRED_COLUMNS = [
    "gene",
    "protein",
    "ncbi_accession",
    "uniprot_accession",
    "sequence_length",
    "metals",
    "status",
]


def load_targets():
    """Load all target genes from the master CSV configuration."""

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(
            f"Target configuration not found:\n{CONFIG_FILE}"
        )

    targets = []

    with CONFIG_FILE.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise ValueError("Target configuration has no header.")

        missing = [
            column
            for column in REQUIRED_COLUMNS
            if column not in reader.fieldnames
        ]

        if missing:
            raise ValueError(
                f"Missing required columns: {', '.join(missing)}"
            )

        for row in reader:
            if not row.get("gene"):
                continue

            metals = [
                metal.strip()
                for metal in row["metals"].split(";")
                if metal.strip()
            ]

            target = {
                "gene": row["gene"].strip(),
                "protein": row["protein"].strip(),
                "ncbi_accession": row["ncbi_accession"].strip(),
                "uniprot_accession": row["uniprot_accession"].strip(),
                "sequence_length": int(row["sequence_length"]),
                "metals": metals,
                "status": row["status"].strip(),
            }

            targets.append(target)

    return targets


def get_target(gene):
    """Return one target by gene name."""

    targets = load_targets()

    for target in targets:
        if target["gene"].lower() == gene.lower():
            return target

    raise ValueError(
        f"Target gene '{gene}' was not found in configuration."
    )


def get_metal_targets(metal):
    """Return all genes associated with a particular metal."""

    targets = load_targets()

    metal = metal.strip()

    return [
        target
        for target in targets
        if metal in target["metals"]
    ]


def print_target(target):
    """Print a target in a readable format."""

    print("Target")
    print("-" * 50)

    print(f"Gene:              {target['gene']}")
    print(f"Protein:            {target['protein']}")
    print(f"NCBI accession:     {target['ncbi_accession']}")
    print(f"UniProt accession:  {target['uniprot_accession']}")
    print(f"Sequence length:    {target['sequence_length']} aa")
    print(f"Metals:             {', '.join(target['metals'])}")
    print(f"Status:             {target['status']}")


def main():

    print("RiceMetalAI Target Configuration")
    print("=" * 60)

    print(f"Configuration file:")
    print(CONFIG_FILE)
    print()

    targets = load_targets()

    print(f"Targets loaded: {len(targets)}")
    print()

    for target in targets:
        print_target(target)
        print()

    print("Configuration loading successful.")


if __name__ == "__main__":
    main()
