"""
Tests for the `alpha_dicarbonyl_alkylation` function from `chemicalbd.bond_disconnector` module.

The `alpha_dicarbonyl_alkylation` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the pattern is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case the pattern is present.

"""

from chemicalbd.bond_disconnector import alpha_dicarbonyl_alkylation
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the pattern is absent
def test_alpha_dicarbonyl_alkylation_absent():
    """
    Test the behavior of `alpha_dicarbonyl_alkylation` function when the pattern is absent.

    Verifies that the function returns an empty list when the alpha-dicarbonyl alkylation pattern is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert alpha_dicarbonyl_alkylation(mol) == []

#Test in the case the function receives an input different from a mol object
def test_alpha_dicarbonyl_alkylation_invalid ():
    """
    Test the behavior of `alpha_dicarbonyl_alkylation` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        alpha_dicarbonyl_alkylation(mol)

#Test in the case the pattern is present
def test_alpha_dicarbonyl_alkylation_present():
    """
    Test the behavior of `alpha_dicarbonyl_alkylation` function when the pattern is present.

    Verifies that the function correctly identifies the alpha-dicarbonyl alkylation pattern and returns the expected list of reactants.

    Returns:
    None
    """
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