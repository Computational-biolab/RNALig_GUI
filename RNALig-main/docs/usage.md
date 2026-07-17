# RNALig Documentation

## Overview

RNALig is a machine learning (ML)-based structure-driven scoring function designed to predict RNA–ligand binding affinity (ΔG, kcal·mol⁻¹). The framework integrates structural, physicochemical, and interaction-based features derived from RNA–ligand complexes.

The RNALig pipeline consists of two sequential stages:

1. Feature extraction from RNA–ligand structures  
2. Binding affinity prediction using a trained ML model  

---

## Input Requirements

RNALig accepts structural input files in the following formats:

- PDB (`.pdb`)
- mmCIF (`.cif`)

Each input file must contain:

- RNA macromolecule coordinates  
- Bound ligand coordinates  
- Standard atomic representation without critical missing atoms  

---

## Workflow

### Step 1: Feature Extraction

Structural features are computed using the `Features_RNALig.py` script. This module extracts RNA, ligand, and complex-level descriptors, including geometric, physicochemical, and interaction features.

Example:

```bash
python Features_RNALig.py \
  --pdb input.pdb \
  --outcsv features.csv \
  --outdir results/
Output:
- features.csv → tabulated feature representation of the complex
- results/ → HTML-based visualization reports for RNA and ligand

### Step 2: Bindingg Affinity Prediction
The extracted features are used as input to the trained ML model.

Inputs:

Trained model (.pkl)
Feature file (features.csv)

The model predicts binding affinity (ΔG) for each complex.

Output:

Predictions_from_model.csv
| Column                                | Description                   |
| ------------------------------------- | ----------------------------- |
| Complex_ID                            | Identifier of input structure |
| Predicted_binding_affinity (kcal/mol) | Predicted binding free energy |
