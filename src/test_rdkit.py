from rdkit import Chem
from rdkit.Chem import Draw

mol = Chem.MolFromSmiles('CCO')
print("RDKit marche ! Molecule:", mol)

Draw.MolToFile(mol, "../data/molecule.png")
print("Image sauvegardee: molecule.png")