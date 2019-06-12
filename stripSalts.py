from rdkit.Chem.SaltRemover import SaltRemover
from rdkit import Chem

remover = SaltRemover()

out=open("data/fragment_library.smi","w")
with open('data/diamond_fragment_library.smi', 'r') as reader:
    for smiles in reader.readlines():
        smiles=smiles.strip()

        rdkit_mol=Chem.MolFromSmiles(smiles)
        res = remover.StripMol(rdkit_mol)
        sm=Chem.MolToSmiles(res)
        out.write(sm+"\n")
out.close()
