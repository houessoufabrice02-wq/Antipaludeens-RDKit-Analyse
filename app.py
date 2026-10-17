import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import py3Dmol

st.set_page_config(page_title="Anima Palu", layout="wide")
st.title("💊 Anima Palu - Mon modèle interactif")
st.caption("Par Fabrice Houessou")

mols = {
    "Chloroquine": "CCN(CC)CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl",
    "Amodiaquine": "CCN(CC)CC1=C(C=CC(=C1)NC2=C3C=CC(=CC3=NC=C2)Cl)O",
    "Quinine": "COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O",
    "Artemether": "COC1OC2C(C)C3CCC4C(C)(C)CCC(C(C)(O)OC)C43C2OO1"
}

choice = st.selectbox("Choisis une molécule :", list(mols.keys()))
smiles = mols[choice]
mol = Chem.MolFromSmiles(smiles)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

col1, col2 = st.columns(2)
with col1:
    mblock = Chem.MolToMolBlock(mol)
    view = py3Dmol.view(width=450, height=450)
    view.addModel(mblock, 'mol')
    view.setStyle({'stick':{}})
    view.setBackgroundColor('0xeeeeee')
    view.zoomTo()
    st.components.v1.html(view._make_html(), height=500)

with col2:
    st.metric("Poids Moléculaire", f"{Descriptors.MolWt(mol):.1f}")
    st.metric("LogP", f"{Descriptors.MolLogP(mol):.2f}")
    st.info(smiles)