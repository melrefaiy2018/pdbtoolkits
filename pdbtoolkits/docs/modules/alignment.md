# Alignment Module

The alignment module provides tools for superimposing molecular structures and analyzing the results.

## Overview

Molecular alignment (or superposition) is a fundamental operation in structural biology and computational chemistry. This module enables accurate alignment of PDB structures with detailed atom mapping.

## Functions

### `align_molecules`

```python
from pdbtoolkits.alignment.superimpose import align_molecules

aligned_mol, rmsd, atom_mapping = align_molecules(ref_pdb_path, target_pdb_path, alignment_atoms=None)
```

Aligns a target molecule to a reference molecule and returns the aligned molecule, RMSD value, and atom mapping.

**Parameters:**
- `ref_pdb_path` (str): Path to the reference PDB file
- `target_pdb_path` (str): Path to the target PDB file to be aligned
- `alignment_atoms` (dict, optional): Dictionary specifying atoms to use for alignment
  - Keys are element symbols
  - Values are the number of atoms of that element to use
  - Example: `{'Mg': 1, 'N': 4}` uses 1 magnesium atom and 4 nitrogen atoms

**Returns:**
- `aligned_mol`: RDKit molecule object of the aligned target
- `rmsd`: Root Mean Square Deviation value (float)
- `atom_mapping`: List of dictionaries with atom mapping information

**Example:**
```python
# Align using all heavy atoms
aligned_mol, rmsd, mapping = align_molecules("reference.pdb", "target.pdb")
print(f"RMSD: {rmsd:.4f} Å")

# Align using specific atoms
aligned_mol, rmsd, mapping = align_molecules(
    "reference.pdb", 
    "target.pdb", 
    alignment_atoms={'C': 4, 'N': 2}
)
```

### `save_aligned_molecule`

```python
from pdbtoolkits.alignment.superimpose import save_aligned_molecule

save_aligned_molecule(aligned_mol, output_path)
```

Saves an aligned molecule to a PDB file.

**Parameters:**
- `aligned_mol`: RDKit molecule object to save
- `output_path` (str): Path where the PDB file should be saved

**Returns:**
- Path to the saved file

**Example:**
```python
# Align a molecule and save the result
aligned_mol, rmsd, mapping = align_molecules("reference.pdb", "target.pdb")
save_aligned_molecule(aligned_mol, "aligned_target.pdb")
```

## Implementation Details

The alignment module uses RDKit's molecular alignment capabilities with additional functionality for:

1. **Custom Atom Selection**: Specify which atoms to use for alignment with the `alignment_atoms` parameter
2. **Comprehensive Atom Mapping**: Generate detailed mapping between reference and target structures
3. **Detailed Output**: Provide RMSD values and mapping information for analysis

## Example Use Cases

### Basic Structure Alignment

```python
from pdbtoolkits.alignment.superimpose import align_molecules, save_aligned_molecule

# Align two PDB structures
aligned_mol, rmsd, mapping = align_molecules("protein_native.pdb", "protein_model.pdb")
print(f"Alignment RMSD: {rmsd:.4f} Å")

# Save the aligned structure
save_aligned_molecule(aligned_mol, "protein_model_aligned.pdb")
```

### Chlorophyll Alignment Using Metal Centers

```python
from pdbtoolkits.alignment.superimpose import align_molecules

# Align using magnesium centers and nitrogen atoms in chlorophyll
aligned_mol, rmsd, mapping = align_molecules(
    "chlorophyll_a.pdb",
    "chlorophyll_b.pdb",
    alignment_atoms={'Mg': 1, 'N': 4}
)
print(f"Chlorophyll alignment RMSD: {rmsd:.4f} Å")
```

### Analyzing the Atom Mapping

```python
# Align structures and examine the mapping
aligned_mol, rmsd, mapping = align_molecules("ref.pdb", "target.pdb")

# Print summary of atom mapping
print(f"{'Reference':<30} {'Target':<30}")
print("-" * 60)
for atom_map in mapping[:10]:  # First 10 atoms
    ref = f"{atom_map['ref_atom_name']} ({atom_map['ref_element']}) {atom_map['ref_residue']} {atom_map['ref_residue_number']}"
    target = f"{atom_map['target_atom_name']} ({atom_map['target_element']}) {atom_map['target_residue']} {atom_map['target_residue_number']}"
    print(f"{ref:<30} {target:<30}")
```

## Tips for Successful Alignment

1. **Choose appropriate reference and target structures**
   - Structures should have similar overall topology
   - Significant conformational differences may require special handling

2. **Specify alignment atoms for complex structures**
   - For metal-containing structures, using the metal center can improve alignment
   - For proteins, consider using backbone atoms
   - For ligands, use rigid core atoms

3. **Interpret RMSD values in context**
   - Lower RMSD values indicate better alignment
   - Typical values depend on the structure type:
     - < 1.0 Å: Excellent alignment
     - 1.0-2.0 Å: Good alignment
     - > 2.0 Å: Poor alignment or significantly different structures

## Troubleshooting

### Common Issues

1. **"No sub-structure match found between the probe and query mol"**
   - The structures are too different for automatic alignment
   - Try specifying alignment atoms explicitly
   - Check that the PDB files are properly formatted

2. **High RMSD values**
   - Try different sets of alignment atoms
   - Check if your structures have the expected conformation
   - For flexible molecules, consider aligning only rigid substructures

3. **Warning about hydrogen atoms**
   - These warnings are typically benign - RDKit is warning about hydrogen handling
   - The alignment should still proceed correctly