from chemicalbd.bond_disconnector import swern_oxidation
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

'''
The swern_oxidation function is
tested.
'''
#Test in the case an aldehyde is absent
def test_aldehyde_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert swern_oxidation(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_aldehyde_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        swern_oxidation(mol)

#Test in the case an aldehyde is present
def test_aldehyde_present():
    mol = Chem.MolFromSmiles('CCC=O')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = swern_oxidation(mol)[1] #List of reactants that could form the aldehyde that was disconnected
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CS(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCN(CC)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('O=C(Cl)C(=O)Cl')) in smiles_list