"""
Tests for the `flatten_list` function from `chemicalbd.bond_disconnector` module.

The `disconnections` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the function receives an input that is not a list.
Test in the case the function receives an input that is a nested list, but not of only Mol objects.
Test in the case the function receives an valid input

"""
from chemicalbd.bond_disconnector import flatten_list
from rdkit import Chem
import pytest

#Test in the case the input is not a list
def test_not_list ():
    """
    Test behaviour of the 'flatten_list' function with a non-list input

    Verifies that the function raises a TypeError when a non-list input is passed

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCCC')
    with pytest.raises(TypeError):
        flatten_list(mol)

#Test in the case the input is a nested list, but not of Mol objects
def test_not_mol():
    """
    Test behaviour of the 'flatten_list' function with a nested list whose elements are
    not only lists or Mol objects as input

    Verifies that the function raises a TypeError when a nested list whose elements are
    not only lists or Mol objects input is passed

    Returns:
    None
    """
    mol = Chem.MolFromSmiles('CCCC')
    mol_string = 'CC'
    nested_list = [[], [mol, mol_string]]
    with pytest.raises(TypeError):
        flatten_list(nested_list)

#Test in the case of a valid input
def test_valid ():
    """
    Test behaviour of the 'flatten_list' function with a valid input

    Verifies that the function returns a flattened list when a valid input
    is passed

    Returns:
    None
    """
    mol_1 = Chem.MolFromSmiles('CCC')
    mol_2 = Chem.MolFromSmiles('CC')
    mol_3 = Chem.MolFromSmiles('C')
    nested_list = [[], [mol_1], [[mol_2, mol_3]]]
    flattened_list = flatten_list(nested_list)
    smiles_list = [Chem.MolToSmiles(mol) for mol in flattened_list]
    assert smiles_list == ['CCC', 'CC', 'C']