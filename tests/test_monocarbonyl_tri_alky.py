"""
Tests for the `alpha_monocarbonyl_alkylation_trisubstituted` function from `chemicalbd.bond_disconnector` module.

The `alpha_monocarbonyl_alkylation_trisubstituted` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the pattern is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case the pattern is present.

"""
from chemicalbd.bond_disconnector import alpha_monocarbonyl_alkylation_trisubstituted
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the pattern is absent
def test_alpha_monocarbonyl_alkylation_trisubstituted_absent():
    """
    Test the behavior of `alpha_monocarbonyl_alkylation_trisubstituted` function when the pattern is absent.

    Verifies that the function returns an empty list when the specified pattern is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert alpha_monocarbonyl_alkylation_trisubstituted(mol) == []

#Test in the case the function receives an input different from a mol object
def test_alpha_monocarbonyl_alkylation_trisubstituted_invalid ():
    """
    Test the behavior of `alpha_monocarbonyl_alkylation_trisubstituted` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        alpha_monocarbonyl_alkylation_trisubstituted(mol)

#Test in the case the pattern is present
def test_alpha_monocarbonyl_alkylation_trisubstituted_present():
    """
    Test the behavior of `alpha_monocarbonyl_alkylation_trisubstituted` function when the pattern is present.

    Verifies that the function correctly identifies the specified pattern and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CC(=O)C(C)(CCC)CC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = alpha_monocarbonyl_alkylation_trisubstituted(mol) #List of reactants that could form the pattern that was disconnected 
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(CC)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCN(CC)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC(C)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCI')) in smiles_list