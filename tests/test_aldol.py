"""
Tests for the `aldol` function from `chemicalbd.bond_disconnector` module.

The `aldol` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case an alpha-beta unsaturated carbonyl is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case an alpha-beta unsaturated carbonyl is present.

"""

from chemicalbd.bond_disconnector import aldol
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case an alpha-beta unsaturated carbonyl is absent
def test_aldol_absent():
    """
    Test the behavior of `aldol` function when an alpha-beta unsaturated carbonyl is absent.

    Verifies that the function returns [0] when the pattern is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert aldol(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_aldol_invalid ():
    """
    Test the behavior of `aldol` function when an invalid input is provided.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        aldol(mol)

#Test in the case an alpha-beta unsaturated carbonyl is present
def test_aldol_present():
    """
    Test the behavior of `aldol` function when an alpha-beta unsaturated carbonyl is present.

    Verifies that the function correctly identifies the pattern and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CC(=O)C=CC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = aldol(mol)[1] #List of reactants that could form the alpha-beta unsaturated carbonyl that was disconnected
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C=C(C)O[Si](C)(C)C')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C[Si](C)(C)Cl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCN(CC)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('Cl[Ti](Cl)(Cl)Cl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)CC(C)O')) in smiles_list