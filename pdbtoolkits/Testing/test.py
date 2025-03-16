"""
Tests for the alignment module.
"""

import os
import pytest
from rdkit import Chem
from pdbtoolkits.alignment import align_molecules, save_aligned_molecule

# Path to test data directory
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')

# Make sure the test data directory exists
os.makedirs(TEST_DATA_DIR, exist_ok=True)

# Sample PDB content for testing
SAMPLE_PDB_1 = """HEADER    PROTEIN
ATOM      1  N   ALA A   1      -0.966   0.248   0.974  1.00  0.00           N  
ATOM      2  CA  ALA A   1       0.257   0.243  -0.001  1.00  0.00           C  
ATOM      3  C   ALA A   1       1.391  -0.231   0.850  1.00  0.00           C  
ATOM      4  O   ALA A   1       1.391  -0.064   2.050  1.00  0.00           O  
ATOM      5  CB  ALA A   1       0.108   0.531  -1.487  1.00  0.00           C  
END
"""

SAMPLE_PDB_2 = """HEADER    PROTEIN
ATOM      1  N   ALA A   1      -0.666   0.548   0.674  1.00  0.00           N  
ATOM      2  CA  ALA A   1       0.557   0.443  -0.101  1.00  0.00           C  
ATOM      3  C   ALA A   1       1.691  -0.131   0.750  1.00  0.00           C  
ATOM      4  O   ALA A   1       1.691   0.036   1.950  1.00  0.00           O  
ATOM      5  CB  ALA A   1       0.408   0.731  -1.587  1.00  0.00           C  
END
"""

@pytest.fixture(scope="module")
def sample_pdbs():
    """Create sample PDB files for testing."""
    ref_pdb_path = os.path.join(TEST_DATA_DIR, "reference.pdb")
    target_pdb_path = os.path.join(TEST_DATA_DIR, "target.pdb")
    
    with open(ref_pdb_path, "w") as f:
        f.write(SAMPLE_PDB_1)
    
    with open(target_pdb_path, "w") as f:
        f.write(SAMPLE_PDB_2)
    
    return ref_pdb_path, target_pdb_path

def test_align_molecules(sample_pdbs):
    """Test basic alignment functionality."""
    ref_pdb_path, target_pdb_path = sample_pdbs
    
    # Test alignment with all atoms
    aligned_mol, rmsd, atom_mapping = align_molecules(ref_pdb_path, target_pdb_path)
    
    # Check that we got a valid molecule back
    assert aligned_mol is not None
    
    # Check that RMSD is reasonable (should be close to zero for very similar structures)
    assert rmsd < 1.0
    
    # Check that we got atom mappings
    assert len(atom_mapping) > 0
    
    # Check that the mappings contain the expected fields
    first_mapping = atom_mapping[0]
    required_fields = ['ref_atom_name', 'ref_element', 'target_atom_name', 'target_element']
    for field in required_fields:
        assert field in first_mapping
    
def test_align_with_specific_atoms(sample_pdbs):
    """Test alignment using specific atoms."""
    ref_pdb_path, target_pdb_path = sample_pdbs
    
    # Align using only nitrogen atoms
    aligned_mol, rmsd, atom_mapping = align_molecules(
        ref_pdb_path, 
        target_pdb_path,
        alignment_atoms={'N': 1}
    )
    
    # Check that we got a valid result
    assert aligned_mol is not None
    assert rmsd is not None
    
def test_save_aligned_molecule(sample_pdbs):
    """Test saving an aligned molecule."""
    _, target_pdb_path = sample_pdbs
    
    # Load a molecule
    mol = Chem.MolFromPDBFile(target_pdb_path)
    
    # Save it to a new file
    output_path = os.path.join(TEST_DATA_DIR, "saved_molecule.pdb")
    save_aligned_molecule(mol, output_path)
    
    # Check that the file exists and is not empty
    assert os.path.exists(output_path)
    assert os.path.getsize(output_path) > 0
    
    # Check that we can load it back
    loaded_mol = Chem.MolFromPDBFile(output_path)
    assert loaded_mol is not None