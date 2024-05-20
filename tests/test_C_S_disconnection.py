from chemicalbd.bond_disconnector import C_S_disconnection
from rdkit import Chem

def C_S_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert C_S_disconnection(mol) == 0

