from chemicalbd.bond_disconnector import C_S_disconnection
from rdkit import Chem
from IPython.display import display
import pytest

def test_C_S_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert C_S_disconnection(mol) == [0]

def test_C_S_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        C_S_disconnection(mol)

def test_C_S_present():
    mol = Chem.MolFromSmiles('CSCC')
    reactant_list = C_S_disconnection(mol)[1]
    assert Chem.MolToSmiles(reactant_list[0]) == Chem.MolToSmiles(Chem.MolFromSmiles('CI'))
    assert Chem.MolToSmiles(reactant_list[1]) == Chem.MolToSmiles(Chem.MolFromSmiles('CCS'))
    assert Chem.MolToSmiles(reactant_list[2]) == Chem.MolToSmiles(Chem.MolFromSmiles('C(=O)([O-])[O-].[K+].[K+]'))
    assert Chem.MolToSmiles(reactant_list[3]) == Chem.MolToSmiles(Chem.MolFromSmiles('CCI'))
    assert Chem.MolToSmiles(reactant_list[4]) == Chem.MolToSmiles(Chem.MolFromSmiles('CS'))
    assert Chem.MolToSmiles(reactant_list[5]) == Chem.MolToSmiles(Chem.MolFromSmiles('C(=O)([O-])[O-].[K+].[K+]'))