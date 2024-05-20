from chemicalbd.bond_disconnector import alpha_dicarbonyl_alkylation
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

'''
The alpha_dicarbonyl_alkylation function is
tested.
'''
#Test in the case the pattern is absent
def test_alpha_dicarbonyl_alkylation_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert alpha_dicarbonyl_alkylation(mol) == []

#Test in the case the function receives an input different from a mol object
def test_alpha_dicarbonyl_alkylation_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        alpha_dicarbonyl_alkylation(mol)

#Test in the case the pattern is present
def test_alpha_dicarbonyl_alkylation_present():
    mol = Chem.MolFromSmiles('CC(=O)C(C)(CC)C(=O)C')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = alpha_dicarbonyl_alkylation(mol) #List of reactants that could form the pattern that was disconnected 
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)C(C)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC[O-].[Na+]')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC(C(C)=O)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CI')) in smiles_list