"""
Tests for the `unique_list_reactants` function from `chemicalbd.bond_disconnector` module.

The `disconnections` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the function receives an input that is not a tuple.
Test in the case the function receives an input that is a nested tuple, but not of Mol objects.
Test in the case the function receives an valid input

"""
from chemicalbd.bond_disconnector import unique_list_reactants
from rdkit import Chem
import pytest

#Test in the case the functions receives an input that is not a tuple
def test_not_tuple ():
    """
    Test behaviour of the 'unique_list_reactants' function with a non-tuple input

    Verifies that the function raises a TypeError when a non-tuple input is passed

    Returns:
    None
    """
    mol_1 = Chem.MolFromSmiles('CC')
    mol_2 = Chem.MolFromSmiles('C')
    tuple_1 = (mol_1,)
    tuple_2 = (mol_2,)
    list_of_tuples = [tuple_1, tuple_2]
    with pytest.raises(TypeError):
        unique_list_reactants(list_of_tuples)

#Test in the case the functions receives an input that is a nested tuple, but contains elements other than tuples
def test_not_all_tuple ():
    """
    Test behaviour of the 'unique_list_reactants' function with a nested tuple with not all the elements tuples input

    Verifies that the function raises a TypeError when a nested tuple with not all the elements tuples input is passed

    Returns:
    None
    """
    mol_1 = Chem.MolFromSmiles('CC')
    mol_2 = Chem.MolFromSmiles('C')
    tuple_1 = (mol_1,)
    list_1 = [mol_2]
    nested_tuple = [tuple_1, list_1]
    with pytest.raises(TypeError):
        unique_list_reactants(nested_tuple)
        
#Test in the case the functions receives an input that is a tuple of tuples of non-Mol objects
def test_not_mol ():
    """
    Test behaviour of the 'unique_list_reactants' function with a nested tuple of non-Mol objects as the input

    Verifies that the function raises a TypeError when a nested tuple of non-Mol objects is passed as the input

    Returns:
    None
    """
    tuple_1 = ('CC',)
    tuple_2 = ('C',)
    nested_tuple = (tuple_1, tuple_2)
    with pytest.raises(TypeError):
        unique_list_reactants(nested_tuple)

#Test in the case the function receives a valid input
def test_valid ():
    """
    Test behaviour of the 'unique_list_reactants' function with a valid input

    Verifies that the function returns a unique nested list of reactant molecules when a valid input is provided

    Returns:
    None
    """
    mol_1 = Chem.MolFromSmiles('CC')
    mol_2 = Chem.MolFromSmiles('C')
    mol_3 = Chem.MolFromSmiles('CC')
    mol_4 = Chem.MolFromSmiles('C')
    tuple_1 = (mol_1, mol_2)
    tuple_2 = (mol_3, mol_4)
    nested_tuple = (tuple_1, tuple_2)
    reactant_list = unique_list_reactants(nested_tuple) #A list containing a list of mol objects is created
    reactant_list_smiles = [Chem.MolToSmiles(reactant) for reactant in reactant_list[0]] #The SMILES of the mol objects are added in a new list
    assert reactant_list_smiles == ['CC', 'C'] #It is verified that a list of unique reactants was created