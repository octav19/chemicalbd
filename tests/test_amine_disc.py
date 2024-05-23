"""
Tests for the `amine_disconnection` function from `chemicalbd.bond_disconnector` module.

The `amine_disconnection` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case an amine is absent.
Test in the case the function receives an input different from a Mol object.
Test in the case an amine is present.

"""

from chemicalbd.bond_disconnector import amine_disconnection
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case an alpha-beta unsaturated carbonyl is absent
def test_amine_absent():
    """
    Test the behavior of `amine_disconnection` function when an amine is absent.

    Verifies that the function returns [0] when the pattern is not found in the molecule.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCC')
    assert amine_disconnection(mol) == [0]

#Test in the case the function receives an input different from a mol object
def test_amine_invalid ():
    """
    Test the behavior of `amine_disconnection` function when an invalid input is provided.

    Verifies that the function raises a TypeError when a non-Mol object is passed as input.

    Returns:
    None
    """
    mol = 'Ups'
    with pytest.raises(TypeError):
        amine_disconnection(mol)

#Test in the case an alpha-beta unsaturated carbonyl is present
def test_amine_present():
    """
    Test the behavior of `amine_disconnection` function when an amine is present.

    Verifies that the function correctly identifies the pattern and returns the expected list of reactants.

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCN(C)CCC')
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = amine_disconnection(mol)[1] #List of reactants that could form the amine that was disconnected
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCNC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(=O)Cl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('c1ccncc1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCN(C)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('[AlH4-].[Li+]')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCNCC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('O=CCl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCN(C=O)CC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCNC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC(=O)Cl')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCC(=O)N(C)CC')) in smiles_list