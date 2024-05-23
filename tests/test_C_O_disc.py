"""
Tests for the `C_O_disconnection` function from `chemicalbd.bond_disconnector` module.

The `C_O_disconnection` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the ether is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case the ether is present.

"""

from chemicalbd.bond_disconnector import C_O_disconnection
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the bond is absent
def test_C_O_absent():
    """
    Test the behavior of `C_O_disconnection` function when the ether is absent.

    Verifies that the function returns [0] when the ether is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert C_O_disconnection(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_C_O_invalid ():
    """
    Test the behavior of `C_O_disconnection` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        C_O_disconnection(mol)

#Test in the case the bond is present
def test_C_O_present():
    """
    Test the behavior of `C_O_disconnection` function when the ether is present.

    Verifies that the function correctly identifies the ether and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('COCC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = C_O_disconnection(mol)[1] #List of reactants that could form the disconnected ether
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('[Na+].[OH-]')) in smiles_list