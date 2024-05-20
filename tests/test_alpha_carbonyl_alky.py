from chemicalbd.bond_disconnector import alpha_carbonyl_alkylation
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

'''
The alpha_carbonyl_alkylation function is
tested.
'''
#Test in the case the pattern is absent
def test_alpha_carbonyl_alkylation_absent():
    mol = Chem.MolFromSmiles('CCC')
    assert alpha_carbonyl_alkylation(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_alpha_carbonyl_alkylation_invalid ():
    mol = 'Ups'
    with pytest.raises(TypeError):
        alpha_carbonyl_alkylation(mol)

#Test in the case the pattern is present
def test_alpha_carbonyl_alkylation_present():
    mol = Chem.MolFromSmiles('CC(=O)C(CCC)C(=O)CCC(=O)C(CC)(C)C')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        #A list of lists of reactants that could form the disconnected pattern is returned; each list within the list contains the reactants for one type
        # that could participate in one of the three types of carbonyl alkylation
        reactant_list_of_list = alpha_carbonyl_alkylation(mol)[1]
        #The bidimensional list is transformed in an unidimensional list of reactants
    reactant_list = []
    for list_of_reactants in reactant_list_of_list:
        for reactant in list_of_reactants:
            reactant_list.append(reactant)
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC(C)(C)C(=O)CCC(=O)CC(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC[O-].[Na+]')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCOC(=O)CC(=O)C(C)(C)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C(C)=O)C(=O)CI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C(C)=O)C(=O)CC(C(=O)OCC)C(=O)C(C)(C)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C(C)=O)C(=O)CCC(=O)C(C)C')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCN(CC)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C(C)=O)C(=O)CCC(=O)C(C)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CI')) in smiles_list