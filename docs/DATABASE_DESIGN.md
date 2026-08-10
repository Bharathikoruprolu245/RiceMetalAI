# RiceMetalAI Database Design

Version: 1.0

Status: Draft

Author: K.Bharathi


---

# Overview

RiceMetalAI is an AI-assisted structural bioinformatics platform for
analyzing heavy metal transport proteins in rice and designing CRISPR
strategies to reduce heavy metal accumulation.

The database follows the biological workflow:

Heavy Metal
        ↓
Gene
        ↓
Protein
        ↓
Protein Structure
        ↓
Docking
        ↓
Binding Residues
        ↓
Mutation
        ↓
CRISPR Target
        ↓
AI Prediction

---

# Entity Relationship Overview

HeavyMetal

↓

Gene

↓

Protein

↓

ProteinStructure

↓

DockingResult

↓

BindingResidue

↓

Mutation

↓

CRISPRTarget

↓

AIPrediction

Publications connect to Genes and Proteins.

---

# Table 1 : heavy_metals

Purpose

Stores all heavy metals analyzed by RiceMetalAI.

Columns

id

UUID

Primary Key

name

VARCHAR(100)

Example:

Cadmium

symbol

VARCHAR(10)

Example:

Cd

ionic_form

VARCHAR(50)

Example:

Cd2+

category

VARCHAR(50)

Example:

Toxic

description

TEXT

created_at

TIMESTAMP

updated_at

TIMESTAMP

---

# Table 2 : genes

Purpose

Stores rice heavy metal transporter genes.

Columns

id

UUID

Primary Key

symbol

VARCHAR(50)

Example

OsNRAMP5

full_name

TEXT

gene_family

VARCHAR(50)

Example

NRAMP

species

VARCHAR(100)

Default

Oryza sativa

ncbi_gene_id

VARCHAR(50)

chromosome

VARCHAR(20)

strand

VARCHAR(5)

gene_length

INTEGER

function

TEXT

priority

INTEGER

created_at

TIMESTAMP

updated_at

TIMESTAMP

---

# Table 3 : proteins

Purpose

Protein information for every gene.

Columns

id

UUID

gene_id

FK -> genes

accession

VARCHAR(100)

sequence

TEXT

protein_length

INTEGER

molecular_weight

FLOAT

isoelectric_point

FLOAT

localization

VARCHAR(100)

function

TEXT

created_at

TIMESTAMP

updated_at

TIMESTAMP

---

# Table 4 : protein_structures

Purpose

AlphaFold / PDB structures.

Columns

id

UUID

protein_id

FK -> proteins

structure_source

VARCHAR

Examples

AlphaFold

PDB

structure_id

VARCHAR

confidence_score

FLOAT

resolution

FLOAT

download_url

TEXT

created_at

TIMESTAMP

---

# Table 5 : docking_results

Purpose

Stores docking experiments.

Columns

id

UUID

protein_id

FK

metal_id

FK

software

VARCHAR

binding_energy

FLOAT

binding_site

TEXT

docking_date

TIMESTAMP

notes

TEXT

---

# Table 6 : binding_residues

Purpose

Stores amino acids interacting with metals.

Columns

id

UUID

docking_result_id

FK

residue

VARCHAR

Example

ASP

position

INTEGER

interaction_type

VARCHAR

Example

Hydrogen Bond

Metal Coordination

distance

FLOAT

---

# Table 7 : mutations

Purpose

Stores designed amino acid mutations.

Columns

id

UUID

binding_residue_id

FK

wild_type

VARCHAR

mutant

VARCHAR

mutation

VARCHAR

Example

D145A

predicted_effect

TEXT

confidence

FLOAT

---

# Table 8 : crispr_targets

Purpose

Stores CRISPR guide RNAs.

Columns

id

UUID

mutation_id

FK

guide_sequence

TEXT

pam

VARCHAR

editing_type

VARCHAR

Examples

Knockout

Base Editing

Prime Editing

efficiency

FLOAT

off_target_score

FLOAT

---

# Table 9 : ai_predictions

Purpose

Stores AI model predictions.

Columns

id

UUID

protein_id

FK

mutation_id

FK

model_name

VARCHAR

prediction

TEXT

confidence

FLOAT

created_at

TIMESTAMP

---

# Table 10 : publications

Purpose

Scientific references.

Columns

id

UUID

title

TEXT

authors

TEXT

journal

TEXT

year

INTEGER

doi

VARCHAR

url

TEXT
