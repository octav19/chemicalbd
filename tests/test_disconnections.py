"""
Tests for the `disconnections` function from `chemicalbd.bond_disconnector` module.

The `disconnections` function is tested under different scenarios to ensure its correctness.

Test Cases:
Test in the case the function receives an input different from a str.
Test in the case the function receives an invalid SMILES string.
Test in the case no known is present.
Test in the case known disconnections are present

"""
from chemicalbd.bond_disconnector import disconnections
from rdkit import Chem
from IPython.display import display
from unittest.mock import patch
import pytest

#Test in the case the function receives an input different from a str
def test_disconnections_invalid_not_string ():
    """
    Test the behavior of `disconnections` function with a non-str input.

    Verifies that the function raises a TypeError when a non-str object is passed as input.

    Returns:
    None
    """
    mol_smiles = 100
    with pytest.raises(TypeError):
        disconnections(mol_smiles)

#Test in the case the function receives an input that is not a valid SMILES string
def test_disconnections_invalid_not_valid_SMILES ():
    """
    Test the behavior of `disconnections` function with an invalid SMILES str as input.

    Verifies that the function raises a ValueError when an invalid SMILES str is passed as input.

    Returns:
    None
    """
    mol_smiles = 'Hey'
    with pytest.raises(ValueError):
        disconnections(mol_smiles)

#Test in the case the molecule contains no known disconnections
def test_disconnections_absent ():
    """
    Test the behavior of `disconnections` function when all known disconnections are absent.

    Verifies that the function returns an empty list when no known disconnections are found in the molecule.

    Returns:
    None
    """
    mol_smiles = 'CCC'
    assert disconnections(mol_smiles) == []

#Test in the case the molecule contains known disconnections
def test_disconnections_present():
    """
    Test the behavior of `disconnections` function when known disconnections are present.

    Verifies that the function correctly identifies the known disconnections and returns the expected list of reactants.

    Returns:
    None
    """
    mol_smiles = 'CC(O)CC(CCC)C(=O)C'
    with patch('IPython.display.display') as mock_display: #display is a function for IPython; as the tests are done in Python, the function is transformed
        mock_display.side_effect = lambda *args, **kwargs: None #in a function that does nothing, so that errors do not appear
        reactant_list = disconnections(mol_smiles) #List of reactants that could form the patterns that were disconnected
    smiles_list = [Chem.MolToSmiles(reactant) for reactant in reactant_list] #The list of reactants is transformed in a list corresponding to their smiles
    #It is checked if each smiles that should appear in the list of smiles indeed appears
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C[Mg]Br)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C1CCOC1')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(CC=O)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('C[Mg]Br')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCOC(=O)C(CC(C)O)C(C)=O')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC[O-].[Na+]')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCO')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(CC(C)O)(C(C)=O)C(=O)OCC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(C(C)=O)C(=O)OCC')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CC(O)CI')) in smiles_list
    assert Chem.MolToSmiles(Chem.MolFromSmiles('CCCC(CC(C)O)(C(C)=O)C(=O)OCC')) in smiles_list