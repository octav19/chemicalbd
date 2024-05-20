from chemicalbd.bond_disconnector import ester_disconnection
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

'''
The ester_disconnection function is
tested.
'''
#Test in the case the ester is absent
def test_ester_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert ester_disconnection(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_ester_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        ester_disconnection(mol)

#Test in the case an ester is present
def test_ester_present():
    mol = Chem.MolFromSmiles('CC(=O)OCC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = ester_disconnection(mol)[1] #List of reactants that could form the disconnected esters
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('[H+]')) in smiles_list
    