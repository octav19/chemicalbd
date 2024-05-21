"""
Tests for the `swern_oxidation` function from `chemicalbd.bond_disconnector` module.

The `swern_oxidation` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case an aldehyde is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case an aldehyde is present.

"""
from chemicalbd.bond_disconnector import swern_oxidation
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case an aldehyde is absent
def test_aldehyde_absent():
    """
    Test the behavior of `swern_oxidation` function when the aldehyde is absent.

    Verifies that the function returns a list containing 0 when the aldehyde functional group is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert swern_oxidation(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_aldehyde_invalid ():
    """
    Test the behavior of `swern_oxidation` function with invalid input.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        swern_oxidation(mol)

#Test in the case an aldehyde is present
def test_aldehyde_present():
    """
    Test the behavior of `swern_oxidation` function when the aldehyde is present.

    Verifies that the function correctly identifies the aldehyde functional group and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC=O')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = swern_oxidation(mol)[1] #List of reactants that could form the aldehyde that was disconnected
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CS(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCN(CC)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('O=C(Cl)C(=O)Cl')) in smiles_list