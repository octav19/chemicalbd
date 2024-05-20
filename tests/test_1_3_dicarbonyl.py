from chemicalbd.bond_disconnector import dicarbonyl_1_3
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

'''
The dicarbonyl_1_3 function is
tested.
'''
#Test in the case the pattern is absent
def test_dicarbonyl_1_3_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert dicarbonyl_1_3(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_dicarbonyl_1_3_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        dicarbonyl_1_3(mol)

#Test in the case the pattern is present
def test_dicarbonyl_1_3_present():
    mol = Chem.MolFromSmiles('CC(=O)CC(=O)CCC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = dicarbonyl_1_3(mol)[1] #List of reactants that could form the 1,3 dicarbonyl pattern that was disconnected
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C=C(C)N1CCCCC1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C1CCNCC1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(=O)Cl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(=O)CC(C)=[N+]1CCCCC1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C=C(CCC)N1CCCCC1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)Cl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(CC(C)=O)=[N+]1CCCCC1')) in smiles_list