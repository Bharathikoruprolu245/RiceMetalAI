#!/usr/bin/env python3

"""
RiceMetalAI
Generic metal configuration loader.

Loads the centralized metal definitions from:

data/config/metals.json
"""

import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = (
    PROJECT_ROOT
    / "data"
    / "config"
    / "metals.json"
)


def load_metals():
    """Load all configured metals."""

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(
            f"Metal configuration not found:\n{CONFIG_FILE}"
        )

    with CONFIG_FILE.open(
        "r",
        encoding="utf-8"
    ) as handle:

        metals = json.load(handle)

    if not isinstance(metals, dict):
        raise ValueError(
            "Metal configuration must contain "
            "a JSON object."
        )

    return metals


def get_metal(symbol):
    """Return configuration for one metal."""

    metals = load_metals()

    if symbol not in metals:
        available = ", ".join(metals.keys())

        raise KeyError(
            f"Unknown metal '{symbol}'. "
            f"Available metals: {available}"
        )

    return metals[symbol]


def get_configured_metal_symbols():
    """Return configured metal symbols."""

    return list(load_metals().keys())


def print_metal(metal):
    """Print one metal configuration."""

    print("Metal")
    print("-" * 50)

    print(f"Name:                 {metal['name']}")
    print(f"Symbol:               {metal['symbol']}")
    print(f"Atomic number:        {metal['atomic_number']}")
    print(
        "Oxidation states:     "
        f"{', '.join(map(str, metal['common_oxidation_states']))}"
    )
    print(f"Primary ion:          {metal['primary_ion']}")
    print(f"Docking charge:       {metal['docking_charge']}")
    print(f"Category:             {metal['category']}")
    print(
        f"Essential element:    "
        f"{metal['essential_element']}"
    )


def main():

    print("RiceMetalAI Metal Configuration")
    print("=" * 60)

    print("Configuration file:")
    print(CONFIG_FILE)
    print()

    metals = load_metals()

    print(
        f"Metals loaded: {len(metals)}"
    )
    print()

    for symbol, metal in metals.items():

        print_metal(metal)
        print()

    print(
        "Metal configuration loading successful."
    )


if __name__ == "__main__":
    main()
