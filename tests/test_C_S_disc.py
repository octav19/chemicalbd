"""
Tests for the `C_S_disconnection` function from `chemicalbd.bond_disconnector` module.

The `C_S_disconnection` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the bond is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case the bond is present.

"""

from chemicalbd.bond_disconnector import C_S_disconnection
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the bond is absent
def test_C_S_absent():
    """
    Test the behavior of `C_S_disconnection` function when the C-S bond is absent.

    Verifies that the function returns [0] when the C-S bond is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert C_S_disconnection(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_C_S_invalid ():
    """
    Test the behavior of `C_S_disconnection` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        C_S_disconnection(mol)

#Test in the case the bond is present
def test_C_S_present():
    """
    Test the behavior of `C_S_disconnection` function when the C-S bond is present.

    Verifies that the function correctly identifies the C-S bond and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CSCC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = C_S_disconnection(mol)[1] #List of reactants that could form the disconnected C(sp3)-S bonds
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CS')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCS')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C(=O)([O-])[O-].[K+].[K+]')) in smiles_list
    