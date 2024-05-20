from chemicalbd.bond_disconnector import C_S_disconnection
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
from unittest.mock import patch
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
    with patch('IPython.display.display') as mock_display:
        mock_display.side_effect = lambda *args, **kwargs: None
        reactant_list = C_S_disconnection(mol)[1]
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list]
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CS')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCS')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C(=O)([O-])[O-].[K+].[K+]')) in smiles_list
    