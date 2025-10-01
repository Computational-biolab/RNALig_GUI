# final_GUI.py
import streamlit as st
import pandas as pd
import numpy as np
import tempfile
import os
import requests
import joblib
import subprocess
import shutil
import zipfile
import io
import json
import re
from Bio.PDB import PDBParser
import MDAnalysis as mda
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit.Chem import Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
import pubchempy as pcp
import py3Dmol

# --------------------
# Configuration + Caches
# --------------------
@st.cache_resource
def load_model_and_preprocessors():
    """
    Load ML model and preprocessing objects. Update paths if your model files live somewhere else.
    """
    try:
        model = joblib.load("best_rf_model.pkl")
        imputer = joblib.load("imputer.pkl")
        scaler = joblib.load("scaler.pkl")
        target_scaler = joblib.load("target_scaler.pkl")
    except Exception as e:
        # bubble up so caller can display a message
        raise e
    return model, imputer, scaler, target_scaler

# ordering of features expected by model
features = [
    "nucleotide_composition_A", "nucleotide_composition_U",
    "nucleotide_composition_G", "nucleotide_composition_C",
    "gc_content", "minimum_free_energy", "Molecular Weight (g/mol)",
    "Hydrogen Bond Donors", "Hydrogen Bond Acceptors", "LogP", "TPSA (Å²)",
    "Atom Count", "Chiral Atom Count", "num_electrostatic_contacts",
    "avg_electrostatic_distance", "num_vdw_contacts", "avg_vdw_distance"
]

# Feature groups for display
rna_feats_list = [
    "nucleotide_composition_A", "nucleotide_composition_U",
    "nucleotide_composition_G", "nucleotide_composition_C",
    "gc_content", "minimum_free_energy"
]
ligand_feats_list = [
    "Molecular Weight (g/mol)", "Hydrogen Bond Donors", "Hydrogen Bond Acceptors",
    "LogP", "TPSA (Å²)", "Atom Count", "Chiral Atom Count"
]
complex_feats_list = [
    "num_electrostatic_contacts", "avg_electrostatic_distance",
    "num_vdw_contacts", "avg_vdw_distance"
]

# small curated lists
common_ions = {'NA', 'K', 'MG', 'ZN', 'CA', 'CL', 'MN', 'FE', 'CU', 'CO', 'BR', 'IOD', 'NI', 'HG', 'AG', 'CD', 'AU', 'PB', 'RB', 'HOH'}
known_smiles = {
    "BTN": "O=C1CCC(N2C(C1=O)C(C(=O)NC2=O)C(=O)O)C",
    "ACT": "CC(=O)O",
    "P13": "C1CCCCC1",
    "AP7": "C1=CC=CC=C1",
    "ISI": "CCO",
    "MGR": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "P14": "CC(=O)NCCC(=O)O",
    "HPA": "C1=CC=C(C=C1)O"
}

parser = PDBParser(QUIET=True)

# atomic masses (approx, in g/mol)
ELEMENT_MASSES = {
    "H": 1.0079, "C": 12.0107, "N": 14.0067, "O": 15.999, "P": 30.9738, "S": 32.065,
    "CL": 35.453, "BR": 79.904, "F": 18.998, "B": 10.811, "I": 126.90, "ZN": 65.38,
    "MG": 24.305, "CA": 40.078, "FE": 55.845, "NA": 22.99, "K": 39.0983
}

# --------------------
# Utilities: fetch PDB, RNA features, ligand detection, ligand features, complex features
# --------------------
def fetch_pdb_from_rcsb(pdb_id):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp.write(response.content)
            return tmp.name
    return None

def extract_rna_features(pdb_file, pdb_id):
    structure = parser.get_structure(pdb_id, pdb_file)
    chain_seqs = {}
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                resn = residue.get_resname().strip()
                if resn in {"A", "U", "G", "C"}:
                    seq += resn
            if seq:
                prev = chain_seqs.get(chain.id, "")
                if len(seq) > len(prev):
                    chain_seqs[chain.id] = seq
    features_out = {}
    if not chain_seqs:
        return features_out
    longest_chain_id = max(chain_seqs.keys(), key=lambda k: len(chain_seqs[k]))
    rna_sequence = chain_seqs[longest_chain_id]
    L = len(rna_sequence)
    if L > 0:
        features_out['nucleotide_composition_A'] = rna_sequence.count('A') / L
        features_out['nucleotide_composition_U'] = rna_sequence.count('U') / L
        features_out['nucleotide_composition_G'] = rna_sequence.count('G') / L
        features_out['nucleotide_composition_C'] = rna_sequence.count('C') / L
        features_out['gc_content'] = (rna_sequence.count('G') + rna_sequence.count('C')) / L
    else:
        features_out['nucleotide_composition_A'] = features_out['nucleotide_composition_U'] = 0.0
        features_out['nucleotide_composition_G'] = features_out['nucleotide_composition_C'] = 0.0
        features_out['gc_content'] = 0.0
        features_out['minimum_free_energy'] = None
        return features_out

    mfe_val = None
    rnafold_path = shutil.which("RNAfold")
    if rnafold_path:
        try:
            proc = subprocess.run(
                [rnafold_path, "--noPS"],
                input=rna_sequence,
                capture_output=True,
                text=True,
                check=False,
                timeout=20
            )
            stdout = proc.stdout.strip().splitlines()
            if len(stdout) >= 2:
                struct_line = stdout[1].strip()
                if '(' in struct_line and ')' in struct_line:
                    try:
                        last_open = struct_line.rfind('(')
                        last_close = struct_line.rfind(')')
                        mfe_str = struct_line[last_open+1:last_close].strip()
                        mfe_num = mfe_str.split()[0]
                        mfe_val = float(mfe_num)
                    except Exception:
                        mfe_val = None
            else:
                combined = " ".join(stdout)
                m = re.search(r'\((-?\d+\.?\d*)\)', combined)
                if m:
                    try:
                        mfe_val = float(m.group(1))
                    except Exception:
                        mfe_val = None
        except Exception:
            mfe_val = None
    else:
        mfe_val = None

    features_out['minimum_free_energy'] = mfe_val
    return features_out

def detect_ligands_with_chain(pdb_file):
    structure = parser.get_structure("Structure", pdb_file)
    seen = set()
    ligands = []
    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                try:
                    resnum = residue.id[1]
                except Exception:
                    resnum = None
                if resname and resname not in {'A', 'U', 'G', 'C'} and resname.upper() not in common_ions:
                    key = (resname, chain.id, resnum)
                    if key not in seen:
                        seen.add(key)
                        ligands.append({"resname": resname, "chain": chain.id, "resnum": resnum})
    return ligands

def _element_from_atom(atom):
    try:
        el = atom.element.strip().upper()
        if el:
            return el
    except Exception:
        pass
    name = atom.get_name().strip().upper()
    if len(name) >= 2 and name[:2] in ELEMENT_MASSES:
        return name[:2]
    if len(name) >= 1 and name[0] in ELEMENT_MASSES:
        return name[0]
    return None

def _fetch_rcsb_component_smiles(resname):
    try:
        res = requests.get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname}")
        if res.status_code == 200:
            data = res.json()
            chem_comp = data.get('chem_comp', {})
            identifiers = chem_comp.get('chem_comp_identifier', [])
            for ident in identifiers:
                if ident.get('type') and ident['type'].upper() in ('SMILES', 'ISOMERIC_SMILES'):
                    return ident.get('identifier')
            sdf_url = f"https://files.rcsb.org/ligands/view/{resname}_ideal.sdf"
            sdf_res = requests.get(sdf_url)
            if sdf_res.status_code == 200:
                mol = Chem.MolFromMolBlock(sdf_res.text, sanitize=False)
                if mol:
                    Chem.SanitizeMol(mol)
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
                    return smiles
    except Exception:
        pass
    return None

def extract_ligand_features(ligand_name, pdb_file=None, selected_ligand_desc=None):
    feats = {}
    name_upper = ligand_name.strip().upper()
    smiles = known_smiles.get(name_upper)
    if not smiles:
        try:
            compounds = pcp.get_compounds(ligand_name, 'name')
            if compounds:
                smiles = compounds[0].isomeric_smiles
        except Exception:
            smiles = None
    if not smiles:
        smiles = _fetch_rcsb_component_smiles(name_upper)
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                feats['Molecular Weight (g/mol)'] = CalcExactMolWt(mol)
                feats['Hydrogen Bond Donors'] = Lipinski.NumHDonors(mol)
                feats['Hydrogen Bond Acceptors'] = Lipinski.NumHAcceptors(mol)
                feats['LogP'] = Crippen.MolLogP(mol)
                feats['TPSA (Å²)'] = rdMolDescriptors.CalcTPSA(mol)
                feats['Atom Count'] = mol.GetNumAtoms()
                feats['Chiral Atom Count'] = len([a for a in mol.GetAtoms() if a.HasProp("_ChiralityPossible")])
                return feats
        except Exception:
            pass

    if pdb_file and selected_ligand_desc:
        try:
            structure = parser.get_structure("STRUCT", pdb_file)
            found_atoms = []
            found = False
            for model in structure:
                if selected_ligand_desc["chain"] in [ch.id for ch in model]:
                    chain = model[selected_ligand_desc["chain"]]
                    for residue in chain:
                        resname = residue.get_resname().strip()
                        try:
                            resnum = residue.id[1]
                        except Exception:
                            resnum = None
                        if (resname == selected_ligand_desc["resname"] and resnum == selected_ligand_desc["resnum"]):
                            found = True
                            for atom in residue:
                                el = _element_from_atom(atom)
                                if el:
                                    found_atoms.append(el)
                            break
                if found:
                    break
            if not found:
                try:
                    u = mda.Universe(pdb_file)
                    ligand_atoms_md = u.select_atoms(f"resname {selected_ligand_desc['resname']} and chain {selected_ligand_desc['chain']}")
                    if len(ligand_atoms_md) > 0:
                        for a in ligand_atoms_md.atoms:
                            el = getattr(a, "element", None)
                            if el:
                                found_atoms.append(el.strip().upper())
                            else:
                                n = a.name.strip().upper()
                                if len(n) >= 2 and n[:2] in ELEMENT_MASSES:
                                    found_atoms.append(n[:2])
                                else:
                                    found_atoms.append(n[0])
                except Exception:
                    pass
            if len(found_atoms) > 0:
                total_mass = 0.0
                n_count = 0
                o_count = 0
                atom_count = 0
                for el in found_atoms:
                    atom_count += 1
                    mass = ELEMENT_MASSES.get(el, None)
                    if mass is None:
                        if len(el) >= 1:
                            mass = ELEMENT_MASSES.get(el[0], 0.0)
                        else:
                            mass = 0.0
                    total_mass += mass
                    if el == 'N':
                        n_count += 1
                    if el == 'O':
                        o_count += 1
                feats['Molecular Weight (g/mol)'] = float(total_mass) if total_mass > 0 else None
                feats['Hydrogen Bond Donors'] = int(n_count)
                feats['Hydrogen Bond Acceptors'] = int(o_count + n_count)
                feats['Atom Count'] = int(atom_count)
                feats['Chiral Atom Count'] = 0
                feats['LogP'] = None
                feats['TPSA (Å²)'] = None
                return feats
        except Exception:
            pass
    return feats

def extract_complex_features(pdb_file, selected_ligand_desc=None):
    feats = {}
    try:
        u = mda.Universe(pdb_file)
    except Exception:
        return None
    rna_atoms = u.select_atoms("resname A C G U")
    if len(rna_atoms) == 0:
        return None
    ligand_coords = []
    ligand_atom_names = []
    if selected_ligand_desc:
        structure = parser.get_structure("STRUCT", pdb_file)
        found = False
        for model in structure:
            if selected_ligand_desc["chain"] in [ch.id for ch in model]:
                chain = model[selected_ligand_desc["chain"]]
                for residue in chain:
                    resname = residue.get_resname().strip()
                    try:
                        resnum = residue.id[1]
                    except Exception:
                        resnum = None
                    if (resname == selected_ligand_desc["resname"] and resnum == selected_ligand_desc["resnum"]):
                        found = True
                        for atom in residue:
                            try:
                                ligand_coords.append(atom.get_vector().get_array())
                                ligand_atom_names.append(atom.get_name())
                            except Exception:
                                continue
                        break
            if found:
                break
        if not found:
            ligand_atoms_md = u.select_atoms(f"resname {selected_ligand_desc['resname']}")
            if len(ligand_atoms_md) == 0:
                ligand_atoms_md = u.select_atoms("not resname A C G U")
            if len(ligand_atoms_md) == 0:
                return None
            ligand_coords = ligand_atoms_md.positions.tolist()
            ligand_atom_names = [a.name for a in ligand_atoms_md.atoms]
    else:
        ligand_atoms_md = u.select_atoms("not resname A C G U")
        if len(ligand_atoms_md) == 0:
            return None
        ligand_coords = ligand_atoms_md.positions.tolist()
        ligand_atom_names = [a.name for a in ligand_atoms_md.atoms]

    ligand_coords = np.array(ligand_coords)
    if ligand_coords.size == 0:
        return None

    rna_positive = rna_atoms.select_atoms("name N*")
    ligand_oxygen_idx = [i for i, n in enumerate(ligand_atom_names) if n.upper().startswith("O")]
    if len(rna_positive) == 0 or len(ligand_oxygen_idx) == 0:
        electrostatic_distances = np.array([])
        contacts = np.array([])
    else:
        ligand_oxygen_coords = ligand_coords[ligand_oxygen_idx, :]
        electrostatic_distances = cdist(rna_positive.positions, ligand_oxygen_coords)
        contacts = electrostatic_distances < 5.0

    feats['num_electrostatic_contacts'] = int(np.sum(contacts)) if electrostatic_distances.size else 0
    feats['avg_electrostatic_distance'] = float(np.mean(electrostatic_distances[contacts])) if electrostatic_distances.size and np.any(contacts) else None

    distances = cdist(rna_atoms.positions, ligand_coords)
    vdw_contacts = distances < 4.0
    feats['num_vdw_contacts'] = int(np.sum(vdw_contacts))
    feats['avg_vdw_distance'] = float(np.mean(distances[vdw_contacts])) if np.any(vdw_contacts) else None

    return feats

# --------------------
# Visualization helpers (py3Dmol)
# --------------------
def show_structure_3d_widget(pdb_data, title="Structure", height=450, width=700):
    """Return html for py3Dmol view from pdb string content"""
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, "pdb")
    view.setStyle({}, {"hidden": True})
    view.addStyle({"resn": ["A", "U", "G", "C"]}, {"cartoon": {"color": "spectrum"}})
    view.addStyle({"not": {"resn": ["A", "U", "G", "C", "HOH", "NA", "K", "MG", "ZN", "CL", "CA"]}}, {"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()
    view.setBackgroundColor("white")
    view.spin(True)  # auto-rotate demo viewer
    return view._make_html()

# --------------------
# Streamlit app (two pages) with persistent header
# --------------------
st.set_page_config(layout="wide", page_title="RNALig — RNA–Ligand Affinity Predictor")

# Header display (logo + title) - shown on all pages
def display_header(logo_filename="logo.png", title="RNALig — The RNA–Ligand Feature Extractor & Binding Affinity Predictor"):
    # try multiple locations for logo: local folder, /mnt/data
    possible_paths = [
        os.path.join(os.getcwd(), logo_filename),
        os.path.join("/mnt/data", logo_filename),
        logo_filename
    ]
    logo_path = None
    for p in possible_paths:
        if os.path.exists(p):
            logo_path = p
            break

    cols = st.columns([1, 3, 1])
    with cols[1]:
        if logo_path:
            try:
                st.image(logo_path, width=320)  # larger logo
            except Exception:
                st.image(logo_path, use_column_width=True)
        st.markdown(f"## {title}")

display_header()  # show logo + title on every page

# Sidebar navigation
st.sidebar.title("RNALig Navigation")
page = st.sidebar.radio("Go to", ["🏠 Home", "🔮 Predict"])

# ----------------- HOMEPAGE -----------------
if page == "🏠 Home":
    st.markdown(
        """
        Welcome to **RNALig** 🎉

        This tool extracts structural and physicochemical features from RNA–ligand complexes (PDB files)
        and predicts binding affinity using a trained ML model.

        **What you can do:**
        - Upload single PDBs or a ZIP of PDBs.
        - Or enter a PDB ID to fetch from RCSB.
        - Choose ligand if multiple are present.
        - Visualize structure and download predictions.

        ---
        """
    )

    st.info("Tip: If a PDB contains many small solvent/ion entries, they are ignored. "
            "The ligand selector lists non-ion ligand residues as `RES — chain X — res N`.")

    # --- RNALig paragraph ---
    st.markdown(
        """
        **RNALig** is an AI-driven platform designed to accelerate RNA–ligand interaction studies by seamlessly 
        combining feature extraction, structural analysis, and binding affinity prediction. It leverages advanced 
        computational pipelines to analyze RNA secondary structure, ligand physicochemical properties, and RNA–ligand 
        complex interactions from experimental PDB files or RCSB entries. With an intuitive interface, RNALig empowers 
        researchers to upload single or multiple structures, automatically detect ligands, and visualize complexes in 
        interactive 3D. By integrating machine learning–based affinity prediction with domain-specific feature 
        engineering, RNALig provides a reliable and user-friendly tool to support RNA-targeted drug discovery, 
        aptamer design, and therapeutic development.
        """
    )

    # --- Movable 3D PDB structure viewer (homepage demo) ---
    demo_names = ["1f27.pdb", os.path.join("/mnt/data", "1f27.pdb")]
    demo_path = None
    for n in demo_names:
        if os.path.exists(n):
            demo_path = n
            break

    if demo_path:
        try:
            with open(demo_path, "r") as fh:
                pdb_text = fh.read()
            # subtitle removed as requested
            html = show_structure_3d_widget(pdb_text, title="Demo structure", height=520, width=900)
            st.components.v1.html(html, height=520)
        except Exception as e:
            st.warning("Failed to render demo PDB viewer.")
            st.exception(e)
    else:
        st.warning("Demo PDB not found (`1f27.pdb`). Place the file next to final_GUI.py or in /mnt/data.")

    # --- Contact info ---
    st.markdown("---")
    st.markdown("📧 **Contact:** computationalBioLab@gmail.com")

# ----------------- PREDICT PAGE -----------------
elif page == "🔮 Predict":
    st.title("Predict Binding Affinity with RNALig")
    st.write("Upload single PDB(s), a ZIP of many PDBs, or fetch a PDB from RCSB. If a PDB contains multiple ligands, choose the ligand using chain + residue number.")

    # Try to load model & preprocessors (lazy load)
    model = None
    imputer = None
    scaler = None
    target_scaler = None
    model_loading_error = None
    try:
        model, imputer, scaler, target_scaler = load_model_and_preprocessors()
    except Exception as e:
        model_loading_error = e
        st.warning("Model or preprocessing files not found or failed to load. Prediction will be disabled until model files are available.")
        st.write("Expected files: best_rf_model.pkl, imputer.pkl, scaler.pkl, target_scaler.pkl")
        # continue so user can still use upload/visualization

    pdb_id_input = st.text_input(" Enter PDB ID (optional to fetch from RCSB):")
    uploaded_files = st.file_uploader(" Upload one or more PDB files", type=["pdb"], accept_multiple_files=True)
    uploaded_zip = st.file_uploader(" Or upload a ZIP containing multiple PDB files", type=["zip"], accept_multiple_files=False)

    files_to_process = []

    # Fetch by PDB ID if provided
    if pdb_id_input:
        fetched_file = fetch_pdb_from_rcsb(pdb_id_input)
        if fetched_file:
            files_to_process.append((f"{pdb_id_input}.pdb", fetched_file))
        else:
            st.error("Failed to fetch PDB from RCSB. Check PDB ID or network.")

    # Process individual uploaded PDB files
    for uploaded_file in uploaded_files:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
            tmp.write(uploaded_file.read())
            files_to_process.append((uploaded_file.name, tmp.name))

    # Process uploaded ZIP file (extract .pdb files)
    if uploaded_zip:
        try:
            zip_bytes = uploaded_zip.read()
            z = zipfile.ZipFile(io.BytesIO(zip_bytes))
            pdb_names = [n for n in z.namelist() if n.lower().endswith(".pdb")]
            if len(pdb_names) == 0:
                st.warning("No PDB files found in the uploaded ZIP.")
            else:
                st.info(f"Found {len(pdb_names)} PDB file(s) in ZIP — extracting...")
                for name in pdb_names:
                    with z.open(name) as file_in_zip:
                        content = file_in_zip.read()
                        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
                            tmp.write(content)
                            base = os.path.basename(name)
                            files_to_process.append((base, tmp.name))
        except Exception as e:
            st.error("Failed to process ZIP file. Make sure it's a valid ZIP containing .pdb files.")
            st.exception(e)

    all_results = []

    if files_to_process:
        st.markdown(f"**Processing {len(files_to_process)} file(s).**")
        for file_index, (file_name, pdb_path) in enumerate(files_to_process, start=1):
            st.header(f" Processing {file_name}  ({file_index}/{len(files_to_process)})")
            pdb_id = os.path.splitext(os.path.basename(file_name))[0]

            # Detect ligands with chain+resid
            ligand_descs = detect_ligands_with_chain(pdb_path)
            rna_feats = extract_rna_features(pdb_path, pdb_id)

            # allow batch mode: predict for all ligands automatically
            batch_mode = st.checkbox("Predict for ALL detected ligands (batch mode)", value=False, key=f"batch_{pdb_id}_{file_index}")

            # Selection UI (chain-aware)
            selected_ligand_desc = None
            non_ion_ligands = [d for d in ligand_descs if d["resname"].upper() not in common_ions]
            if len(non_ion_ligands) == 0:
                st.info("No non-ion ligands detected in this PDB.")
            elif len(non_ion_ligands) == 1:
                selected_ligand_desc = non_ion_ligands[0]
                st.info(f"Automatically selected ligand: {selected_ligand_desc['resname']} — chain {selected_ligand_desc['chain']} — res {selected_ligand_desc['resnum']}")
            else:
                labels = [f"{d['resname']} — chain {d['chain']} — res {d['resnum']}" for d in non_ion_ligands]
                choice = st.selectbox(
                    "Multiple ligands detected — choose ligand (resname — chain — resnum):",
                    options=labels,
                    index=0,
                    key=f"ligand_select_{pdb_id}_{file_index}"
                )
                idx = labels.index(choice)
                selected_ligand_desc = non_ion_ligands[idx]

            # If batch mode: iterate all non_ion_ligands and make a small results table
            if batch_mode and len(non_ion_ligands) > 0:
                batch_rows = []
                for desc in non_ion_ligands:
                    ligand_feats = extract_ligand_features(desc["resname"], pdb_file=pdb_path, selected_ligand_desc=desc)
                    complex_feats = extract_complex_features(pdb_path, desc)
                    combined = {**(rna_feats or {}), **(ligand_feats or {}), **(complex_feats or {})}
                    if combined and model is not None:
                        df_tmp = pd.DataFrame([combined])
                        for feat in features:
                            if feat not in df_tmp.columns:
                                df_tmp[feat] = np.nan
                        X = df_tmp[features]
                        try:
                            X_imp = imputer.transform(X)
                            X_scaled = scaler.transform(X_imp)
                            y_pred_scaled = model.predict(X_scaled)
                            y_pred = target_scaler.inverse_transform(y_pred_scaled.reshape(-1, 1)).flatten()[0]
                        except Exception:
                            y_pred = None
                    else:
                        y_pred = None
                    batch_rows.append({
                        "PDB": pdb_id,
                        "Ligand": f"{desc['resname']}|chain{desc['chain']}|res{desc['resnum']}",
                        "Predicted Affinity": y_pred
                    })
                if batch_rows:
                    st.subheader("Batch predictions for all ligands")
                    st.table(pd.DataFrame(batch_rows))
            else:
                ligand_feats = {}
                if selected_ligand_desc:
                    ligand_feats = extract_ligand_features(selected_ligand_desc["resname"], pdb_file=pdb_path, selected_ligand_desc=selected_ligand_desc)
                else:
                    if len(ligand_descs) == 1:
                        selected = ligand_descs[0]
                        ligand_feats = extract_ligand_features(selected["resname"], pdb_file=pdb_path, selected_ligand_desc=selected)
                        selected_ligand_desc = selected

                complex_feats = extract_complex_features(pdb_path, selected_ligand_desc)
                combined = {**(rna_feats or {}), **(ligand_feats or {}), **(complex_feats or {})}

                if combined:
                    df = pd.DataFrame([combined])
                    for feat in features:
                        if feat not in df.columns:
                            df[feat] = np.nan

                    X = df[features]
                    if model is None:
                        st.warning("Model not loaded — prediction unavailable. You can still view features and visualizations.")
                        y_pred = None
                    else:
                        try:
                            X_imp = imputer.transform(X)
                            X_scaled = scaler.transform(X_imp)
                            y_pred_scaled = model.predict(X_scaled)
                            y_pred = target_scaler.inverse_transform(y_pred_scaled.reshape(-1, 1)).flatten()[0]
                        except Exception as e:
                            st.error("Prediction failed. Check your model and preprocessing pipeline.")
                            st.exception(e)
                            y_pred = None

                    if y_pred is not None:
                        df["Predicted Binding Affinity"] = y_pred
                    df["PDB ID"] = pdb_id
                    df["Selected Ligand"] = f"{selected_ligand_desc['resname']}|chain{selected_ligand_desc['chain']}|res{selected_ligand_desc['resnum']}" if selected_ligand_desc else ""

                    if y_pred is not None:
                        st.metric(label=" Predicted Binding Affinity (kcal/mol)", value=f"{y_pred:.2f}")
                    else:
                        st.info("No prediction available for this system (model missing or prediction failed).")

                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.subheader(" RNA Features")
                        display_rna = df[rna_feats_list].T
                        display_rna.columns = ["value"]
                        st.table(display_rna)
                    with col2:
                        st.subheader(" Ligand Features")
                        display_lig = df[ligand_feats_list].T
                        display_lig.columns = ["value"]
                        st.table(display_lig)
                    with col3:
                        st.subheader(" Complex Interaction Features")
                        display_comp = df[complex_feats_list].T
                        display_comp.columns = ["value"]
                        st.table(display_comp)

                    all_results.append(df)

                    # Visualizations
                    try:
                        with open(pdb_path, "r") as fh:
                            pdb_text = fh.read()
                        html = show_structure_3d_widget(pdb_text, title="RNA + Ligand Complex", height=450, width=800)
                        st.components.v1.html(html, height=450)
                        # optionally show separate RNA / ligand views (commented if too slow)
                        # show_rna_only(pdb_path)  # you can implement similar helper if needed
                        # show_ligand_only(pdb_path, selected_ligand_desc)
                    except Exception as e:
                        st.warning("Visualization failed for this file.")
                        st.exception(e)
                else:
                    st.warning(f" Could not extract features from {pdb_id}.")

        if all_results:
            final_df = pd.concat(all_results, ignore_index=True)
            csv = final_df.to_csv(index=False).encode("utf-8")
            st.download_button(" Download All Results as CSV", data=csv, file_name="rna_ligand_predictions.csv", mime="text/csv")
    else:
        st.info("Upload a PDB file, upload a ZIP of PDBs, or enter a PDB ID to begin.")

    st.markdown("---")
    st.markdown("Tip: If a PDB contains many small solvent/ion entries, they are ignored. The ligand selector lists non-ion ligand residues as `RES — chain X — res N`.")
