"""
Tests for the `acyl_chloride_disconnection` function from `chemicalbd.bond_disconnector` module.

The `acyl_chloride_disconnection` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the acyl chloride is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case an acyl chloride is present.

"""

from chemicalbd.bond_disconnector import acyl_chloride_disconnection
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the ester is absent
def test_acyl_chloride_absent():
    """
    Test the behavior of `acyl_chloride_disconnection` function when the acyl chloride is absent.

    Verifies that the function returns [0] when the ester bond is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert acyl_chloride_disconnection(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_acyl_chloride_invalid ():
    """
    Test the behavior of `acyl_chloride_disconnection` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        acyl_chloride_disconnection(mol)

#Test in the case an ester is present
def test_acyl_chloride_present():
    """
    Test the behavior of `acyl_chloride_disconnection` function when the acyl_chloride is present.

    Verifies that the function correctly identifies the amide and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CC(=O)Cl')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = acyl_chloride_disconnection(mol)[1] #List of reactants that could form the disconnected acyl chloride
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('O=S(Cl)Cl')) in smiles_list