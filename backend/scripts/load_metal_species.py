#!/usr/bin/env python3

"""
RiceMetalAI
Metal chemical-species configuration loader.
"""

import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = (
    PROJECT_ROOT
    / "data"
    / "config"
    / "metal_species.json"
)


def load_species():
    """Load all configured chemical species."""

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(
            f"Metal species configuration not found:\n"
            f"{CONFIG_FILE}"
        )

    with CONFIG_FILE.open(
        "r",
        encoding="utf-8"
    ) as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(
            "Metal species configuration "
            "must be a JSON object."
        )

    return data


def get_species(metal):
    """Return all configured species for a metal."""

    species = load_species()

    if metal not in species:
        available = ", ".join(species.keys())

        raise KeyError(
            f"Unknown metal '{metal}'. "
            f"Available metals: {available}"
        )

    return species[metal]


def get_default_species(metal):
    """Return the species marked as default for docking."""

    candidates = get_species(metal)

    defaults = [
        item
        for item in candidates
        if item.get("default_for_docking") is True
    ]

    if len(defaults) == 0:
        raise ValueError(
            f"No default docking species configured "
            f"for {metal}."
        )

    if len(defaults) > 1:
        raise ValueError(
            f"Multiple default docking species "
            f"configured for {metal}."
        )

    return defaults[0]


def get_species_by_id(metal, species_id):
    """Return one species by ID."""

    candidates = get_species(metal)

    for item in candidates:
        if item["species_id"] == species_id:
            return item

    available = ", ".join(
        item["species_id"]
        for item in candidates
    )

    raise KeyError(
        f"Species '{species_id}' not found for "
        f"{metal}. Available: {available}"
    )


def print_species(metal, species):
    """Print one species."""

    print(
        f"{metal} → "
        f"{species['species_id']}"
    )

    print(
        f"  Name:              "
        f"{species['name']}"
    )

    print(
        f"  Formula:           "
        f"{species['formula']}"
    )

    print(
        f"  Oxidation state:   "
        f"{species['oxidation_state']}"
    )

    print(
        f"  Charge:            "
        f"{species['charge']}"
    )

    print(
        f"  Type:              "
        f"{species['species_type']}"
    )

    print(
        f"  Default docking:   "
        f"{species['default_for_docking']}"
    )


def main():

    print("=" * 70)
    print(
        "RiceMetalAI - Metal Species Configuration"
    )
    print("=" * 70)
    print()

    species = load_species()

    total = sum(
        len(items)
        for items in species.values()
    )

    print(
        f"Metals configured: {len(species)}"
    )

    print(
        f"Species configured: {total}"
    )

    print()

    for metal, candidates in species.items():

        for item in candidates:
            print_species(
                metal,
                item
            )
            print()

    print(
        "Metal species configuration "
        "loading successful."
    )


if __name__ == "__main__":
    main()
