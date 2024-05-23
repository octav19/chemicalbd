"""
Tests for the `amino_alcohol_1_2` function from `chemicalbd.bond_disconnector` module.

The `amino_alcohol_1_2` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the 1,2 amino-alcohol is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case the 1,2 amino-alcohol is present.

"""

from chemicalbd.bond_disconnector import amino_alcohol_1_2
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the pattern is absent
def test_amino_alcohol_1_2_absent():
    """
    Test the behavior of `amino_alcohol_1_2` function when the 1,2 amino-alcohol is absent.

    Verifies that the function returns [0] when the pattern is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert amino_alcohol_1_2(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_amino_alcohol_1_2_invalid ():
    """
    Test the behavior of `amino_alcohol_1_2` function when an invalid input is provided.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        amino_alcohol_1_2(mol)

#Test in the case the pattern is present
def test_amino_alcohol_1_2_present():
    """
    Test the behavior of `amino_alcohol_1_2` function when the 1,2 amino-alcohol is present.

    Verifies that the function correctly identifies the pattern and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CC(O)CNCCC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = amino_alcohol_1_2(mol)[1] #List of reactants that could form the 1,2 amino-alcohol that was disconnected 
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC1CO1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCN')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('[NH2-].[Na+]')) in smiles_list