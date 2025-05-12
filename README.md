# pdb-toolKits

A comprehensive collection of Python tools for working with PDB (Protein Data Bank) files. This toolkit provides utilities for structural alignment, analysis, visualization of molecular structures and other useful capapilities that would be released soon.

## Features

- **Molecular Alignment**: Align molecules using specified atoms or all heavy atoms
- **Atom Mapping**: Generate detailed atom mapping between aligned structures
- **Visualization**: Generate configuration for NGL Viewer to visualize molecular structures

## Installation

```bash
# Clone the repository
git clone https://github.com/melrefaiy2018/pdb-toolKits.git
cd pdb-toolKits

# Install dependencies
conda create -n pdbtoolkits
conda activate pdbtoolkits

# Install in development mode
pip install -e .
```

## Dependencies

- RDKit: For molecular manipulation and alignment
- Logging: For error handling and debugging

## Usage Examples

### Molecular Alignment

```python
from pdbtoolkits.alignment import align_molecules, save_aligned_molecule

# Align target to reference using all heavy atoms
aligned_mol, rmsd, atom_mapping = align_molecules(
    "reference.pdb", 
    "target.pdb"
)
print(f"RMSD: {rmsd:.4f} Å")

# Save the aligned molecule
save_aligned_molecule(aligned_mol, "aligned_target.pdb")

# Align using specific atoms (e.g., Mg and 4 nitrogen atoms)
specific_alignment = align_molecules(
    "reference.pdb",
    "target.pdb",
    alignment_atoms={'Mg': 1, 'N': 4}
)
```

### Visualization with NGL

```python
from pdbtoolkits.visualization import prepare_ngl_visualization

# Generate NGL configuration for visualization
ngl_config = prepare_ngl_visualization("reference.pdb", "aligned_target.pdb")

# Use the config with NGL Viewer in a web application
```

## Project Structure

```
pdb-toolKits/
├── pdbtoolkits/            # Main package
│   ├── alignment/          # Alignment tools
│   ├── visualization/      # Visualization tools
├── tests/                  # Test directory
└── docs/                   # Documentation
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
