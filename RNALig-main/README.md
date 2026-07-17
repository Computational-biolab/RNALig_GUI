# RNALig: ML-driven Structure Based Scoring Function for RNA–Ligand Binding Affinity

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](#requirements)
[![Colab CleanPDB](https://img.shields.io/badge/Colab-Clean%20PDB-black?logo=googlecolab)](https://colab.research.google.com/drive/1bKYbTiqtdPYGR4hRPqJp8IeC7GptRW_s)
[![Colab Predictor](https://img.shields.io/badge/Colab-Binding%20Affinity%20Predictor-black?logo=googlecolab)](https://colab.research.google.com/drive/1_BvsKnuCg0B3oYg2CS1ucm9YviEf2h__#scrollTo=w3m8xFJWyOKm)


---

## Overview

**RNALig** predicts RNA–small molecule binding free energy (ΔG, kcal·mol⁻¹) using structure-aware features derived from RNA, ligand, and complex geometries.
It performs **automated feature extraction**, **binding affinity prediction**, and **validation** with interpretable metrics.

**Supported input formats:** **PDB (`.pdb`)** and **mmCIF (`.cif`)**. Multi-model files are auto-split; hetero-ligands are detected automatically.

This repository contains:

* **Feature extraction (Linux CLI)** – batch processing of PDB/mmCIF files with automatic Results folder creation.
* **Binding-affinity predictor (Google Colab)** – ΔG prediction and visualization. Link: https://colab.research.google.com/drive/1_BvsKnuCg0B3oYg2CS1ucm9YviEf2h__#scrollTo=w3m8xFJWyOKm
* **Clean-PDB utility (Google Colab)** – structure cleanup and ligand chain handling. Link : https://colab.research.google.com/drive/1bKYbTiqtdPYGR4hRPqJp8IeC7GptRW_s

> Every feature is physically validated, unit-checked, and summarized via HTML visualization reports for transparency and reproducibility.

---

## Repository Layout

```
RNALig/
├── Codes                           # Unified script for features extraction for RNA, Ligand, and Complex features (PDB/mmCIF)
├── models/                         # Trained model
├── notebooks/                      # Colab notebooks (Predictor, Clean-PDB)                          
├── Data.zip                        # Results of training and testing for all the tools
├── LICENSE
├── README.md
└── environment.yml
```

---

## Quick Start (Feature Extraction)

### Step 1: Create Environment

```bash
conda create -n rnalig python=3.10 -y
conda activate rnalig
conda install -c conda-forge rdkit openbabel mdanalysis biopython freesasa -y
conda install -c bioconda viennarna -y
pip install numpy pandas scipy py3Dmol scikit-learn joblib openpyxl
```

Alternatively, use the provided environment file:

```bash
conda env create -f environment.yml
conda activate rnalig
```

---

### Step 2: Prepare Input Folder

Place your **PDB** or **mmCIF** files in a directory:

```bash
mkdir input_structures
cp path/to/1f27.pdb input_structures/  or
cp path/to/1f27.cif input_structures/
```

---

### Step 3: Run Feature Extraction

#### Batch Mode (PDB + mmCIF)

Extract features for all structures in a directory:

```bash
python Features_RNALig.py \
  --indir ./input_structures \
  --outcsv final_features.csv \
  --outdir ./viz \
  --viz_rna --viz_ligand
```

#### Single Structure Mode

Run on one RNA–ligand complex (PDB or mmCIF):

```bash
# PDB example
python Features_RNALig.py \
  --pdb ./input_structures/1f27.pdb \
  --outcsv final_features.csv \
  --outdir ./viz \
  --viz_rna --viz_ligand

# mmCIF example
python Features_RNALig.py \
  --pdb ./input_structures/1f27.cif \
  --outcsv final_features.csv \
  --outdir ./viz \
  --viz_rna --viz_ligand
```

---

### Step 4: Features Output Folder Structure

After execution, a **viz/** directory is automatically created:

```
viz/
├── final_features.csv              # Combined features for all complexes
├── RNA/                            # RNA-specific HTML reports
│   ├── 1f27_RNA_features.html
│   └── ...
├── Ligand/                         # Ligand-specific HTML reports
│   ├── 1f27_Ligand_features.html
│   └── ...
└── logs/ (optional)                # Log files if errors occur
```

Each **HTML file** visually summarizes feature extraction (pocket detection, SASA, shape metrics, etc.).

---

## Advanced Options

Use these flags for control over extraction and visualization:

* `--viz_rna`, `--viz_ligand` – generate per-structure HTML reports.
* `--pocket_cutoff 5.0` – distance (Å) for defining RNA pocket residues.
* `--pocket_sasa 0.05` – pocket SASA threshold.
* `--rna_label_topk 5` – top residue labels shown near the binding site.
* `--min_heavy 4` – minimum heavy atoms required for ligands.
* `--no-require_carbon` – allow ions or metal cofactors as ligands.
* `--keep_ions` – preserve ions in system.
* `--cutoff 5.0` – atomic contact cutoff for RNA–ligand interactions.

Example:

```bash
python Features_RNALig.py \
  --indir Training \
  --outdir Training/Results \
  --outcsv All_Features.csv \
  --results_rna --results_ligand \
  --pocket_cutoff 5.0 --pocket_sasa 0.05 \
  --rna_label_topk 5 --min_heavy 4 \
  --no-require_carbon --keep_ions --cutoff 5.0
```

---

## Features Summary

All RNA, Ligand, and Complex features are computed by **`Features_RNALig.py`**.

## Step 5: Binding Affinity Prediction
Upload :
- Trained model (.pkl)
- Features file (final_features.csv)
Output:
Prediction_from_model.csv
### Output

| Complex_ID | Predicted_binding_affinity (kcal/mol) |
|------------|--------------------------------------|
| 1f27       | -7.85                                |

### Pipeline Summary
RNA–Ligand Structure (PDB/CIF)
↓
Feature Extraction (Features_RNALig.py)
↓
Feature Table (CSV)
↓
RNALig ML Model (.pkl)
↓
Predicted Binding Affinity (ΔG)


### RNA Features

* **Shape metrics:** radius of gyration (Rg), κ² anisotropy, asphericity, acylindricity.
* **Sequence composition:** base counts (A, C, G, U), GC%.
* **Folding stability:** MFE (kcal·mol⁻¹) via ViennaRNA.
* **Solvent exposure:** SASA (Å²) using FreeSASA.
* **Pocket metrics:** residue count, pocket depth, mean and max distances.

### Ligand Features

* **Geometry/shape:** asphericity, eccentricity, inertia eigenvalues.
* **Physicochemical:** LogP, TPSA, HBD/HBA, charge stats, rotatable bonds, MMFF94 energy.
* **Surface:** SASA (polar/nonpolar), vdW volume.

### Complex Features

* **Contacts:** contact count, H-bonds, hydrophobic and vdW interactions.
* **Buried surface area (BSA)** and pocket COM distances.

---

## Validation & Reporting

* All outputs are unit-checked and logged in **HTML reports**.
* Each computed feature is summarized for inspection.
* Users can verify per-structure reports in `results/RNA/` and `results/Ligand/`.

---

## Troubleshooting

* **Missing freesasa:** ensure it’s installed via conda (`conda install -c conda-forge freesasa`).
* **CIF parsing errors:** use valid mmCIF files or re-export via the PDB website.
* **RDKit/OpenBabel conflicts:** avoid installing `rdkit-pypi`; always use conda-forge RDKit.
* **User-site conflicts:** if pip warns `Defaulting to user installation`, set `export PYTHONNOUSERSITE=1`.
## Polymeric Amino Acid Handling
* RNALig is designed for RNA–small molecule complexes. When a PDB structure contains polymeric amino acid residues (e.g., ARG) instead of a discrete small-molecule ligand, RNALig retains RNA-level features only and skips ligand and complex descriptors to preserve chemical validity.
* Such cases are explicitly annotated in the output:
&Ligand_Type = polymeric_amino_acid
&Feature_Extraction_Status = RNA_only_features
---
