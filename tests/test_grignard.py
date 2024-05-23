"""
Tests for the `grignard` function from `chemicalbd.bond_disconnector` module.

The `grignard` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the alcohol is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case the alcohol is present.

"""

from chemicalbd.bond_disconnector import grignard
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the bond is absent
def test_C_O_absent():
    """
    Test the behavior of `grignard` function when the alcohol is absent.

    Verifies that the function returns [0] when the ether is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert grignard(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_C_O_invalid ():
    """
    Test the behavior of `grignard` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        grignard(mol)

#Test in the case the bond is present
def test_C_O_present():
    """
    Test the behavior of `grignard` function when the alcohol is present.

    Verifies that the function correctly identifies the alcohol and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CC(O)CC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = grignard(mol)[1] #List of reactants that could form the disconnected alcohol
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC[Mg]Br')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C1CCOC1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C[Mg]Br')) in smiles_list