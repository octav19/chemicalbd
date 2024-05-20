from chemicalbd.bond_disconnector import C_S_disconnection
from rdkit import Chem

def test_C_S_absent():
    mol = Chem.MolFromSmiles('CSCC')
    assert C_S_disconnection(mol) == 0

