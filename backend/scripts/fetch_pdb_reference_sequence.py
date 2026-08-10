from pathlib import Path

from Bio.PDB import PDBParser, PPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

PDB_FILE = Path("data/structures/pdb/8E6N.pdb")

OUTPUT_DIR = Path("data/sequences/reference")
OUTPUT_FILE = OUTPUT_DIR / "8E6N_chainA.fasta"


def main():

    print("Loading 8E6N structure...")

    if not PDB_FILE.exists():
        raise FileNotFoundError(
            f"PDB file not found: {PDB_FILE}"
        )

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "8E6N",
        PDB_FILE
    )

    model = structure[0]
    chain = model["A"]

    print(f"Chain: {chain.id}")

    ppb = PPBuilder()

    peptides = ppb.build_peptides(chain)

    if not peptides:
        raise RuntimeError(
            "Could not extract protein sequence from chain A."
        )

    sequences = []

    for peptide in peptides:
        sequences.append(str(peptide.get_sequence()))

    sequence = Seq("".join(sequences))

    print(f"Extracted sequence length: {len(sequence)} aa")

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    record = SeqRecord(
        sequence,
        id="8E6N_chainA",
        name="8E6N_chainA",
        description=(
            "Deinococcus radiodurans MntH "
            "NRAMP reference structure chain A"
        ),
    )

    SeqIO.write(
        record,
        OUTPUT_FILE,
        "fasta"
    )

    print()
    print("Reference sequence saved:")
    print(OUTPUT_FILE)


if __name__ == "__main__":
    main()
