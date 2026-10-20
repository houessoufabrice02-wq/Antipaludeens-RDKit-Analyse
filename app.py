import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, FilterCatalog
import py3Dmol
from stmol import showmol

# --- CONFIG ---
st.set_page_config(page_title="Anima Palu", layout="wide")
st.title("💊 Anima Palu - Mon modèle interactif")
st.write("Par Fabrice Houessou - Évaluation antipaludéens et toxicité")
st.divider()

# --- 4 MOLÉCULES ANTIPALUDÉENNES ---
molecules = {
    "Chloroquine": "CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl",
    "Artémisinine": "C[C@@H]1CC[C@H]2[C@@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C",
    "Quinine": "COC1=CC2=C(C=CN=C2C=C1)[C@H]([C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)O",
    "Méfloquine": "FC(F)(F)c1cc(cc2c1nccc2[C@@H](O)[C@H]3N4CCCCC4CCC3)C(F)(F)F"
}

# --- FILTRES TOXICITÉ ---
@st.cache_resource
def get_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)

catalog = get_catalog()

def calculer_toxicite(mol):
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    pains = catalog.GetMatches(mol)
    
    if len(pains) > 0 or logp > 5:
        risque, color = "Élevé ⚠️", "red"
    elif logp > 3.5 or tpsa < 40:
        risque, color = "Modéré 🟡", "orange"
    else:
        risque, color = "Faible ✅", "green"
    
    lipinski = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
    
    return {"MW": round(mw,1), "LogP": round(logp,2), "TPSA": round(tpsa,1), "Risque": risque, "Color": color, "Lipinski": lipinski, "PAINS": pains}

# --- INTERFACE ---
choix = st.selectbox("Choisis une molécule :", list(molecules.keys()))
smiles = molecules[choix]

mol = Chem.MolFromSmiles(smiles)
mol_3d = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol_3d)

tox = calculer_toxicite(mol)

col1, col2 = st.columns([1,1.2])

with col1:
    st.subheader(f"🧬 {choix} - Vue 3D")
    block = Chem.MolToMolBlock(mol_3d)
    view = py3Dmol.view(width=450, height=450)
    view.addModel(block, 'mol')
    view.setStyle({'stick': {}})
    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()
    showmol(view, height=450, width=450)
    st.code(smiles, language="text")

with col2:
    st.subheader("📊 Propriétés Physico-chimiques")
    c1, c2, c3 = st.columns(3)
    c1.metric("Poids Mol.", f"{tox['MW']} Da")
    c2.metric("LogP", tox['LogP'])
    c3.metric("TPSA", f"{tox['TPSA']} Å²")
    
    st.divider()
    st.subheader("🧪 Profil Toxicité Prédit")
    
    if tox["Color"] == "green":
        st.success(f"**Risque Global : {tox['Risque']}**")
    elif tox["Color"] == "orange":
        st.warning(f"**Risque Global : {tox['Risque']}**")
    else:
        st.error(f"**Risque Global : {tox['Risque']}**")

    if tox["Lipinski"]:
        st.success("Règle de Lipinski : PASS ✅ - Bonne biodisponibilité orale")
    else:
        st.warning("Règle de Lipinski : FAIL ❌ - Mauvaise absorption possible")
    
    if len(tox["PAINS"]) == 0:
        st.success("Filtre PAINS : PASS ✅ - Aucune structure toxique/réactive détectée")
        st.success("Filtre BRENK : PASS ✅ - Structure propre")
    else:
        for m in tox["PAINS"]:
            st.error(f"ALERTE PAINS/BRENK : {m.GetDescription()} - Potentiellement toxique")
    
    st.info(f"""
    **Interprétation pour {choix} :**
    - LogP = {tox['LogP']} : {'Trop lipophile -> risque foie' if tox['LogP'] > 5 else 'Bon équilibre'}
    - TPSA = {tox['TPSA']} : {'Bonne perméabilité' if tox['TPSA'] < 90 else 'Faible passage membranaire'}
    """)
