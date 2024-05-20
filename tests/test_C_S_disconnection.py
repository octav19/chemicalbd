from chemicalbd.bond_disconnector import C_S_disconnection
from rdkit import Chem
import pytest

def test_C_S_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert C_S_disconnection(mol) == 0

def test_C_S_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        C_S_disconnection(mol)
    

