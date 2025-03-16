#!/usr/bin/env python
"""
Example script demonstrating how to use pdb-toolKits to align molecules
and generate visualization configuration.

Usage:
  python align_and_visualize.py reference.pdb target.pdb --output aligned.pdb
  python align_and_Visualize.py CLA_602_reference.pdb CLA_602_target.pdb --output aligned_CLA_602_target.pdb
"""

import argparse
import json
import logging
import os
import sys

from pdbtoolkits.alignment.superimpose import align_molecules, save_aligned_molecule
from pdbtoolkits.visualization.ngl import prepare_ngl_visualization

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description='Align PDB molecules and generate visualization')
    parser.add_argument('reference', help='Reference PDB file')
    parser.add_argument('target', help='Target PDB file to align')
    parser.add_argument('--output', '-o', default='aligned_output.pdb', help='Output path for aligned PDB')
    parser.add_argument('--vis_config', '-v', default='ngl_config.json', help='Output path for NGL visualization config')
    parser.add_argument('--atoms', '-a', nargs='+', help='Specific atoms to use for alignment (format: Element:Count, e.g., "Mg:1 N:4")')
    args = parser.parse_args()

    # Parse atom specifications if provided
    alignment_atoms = None
    if args.atoms:
        alignment_atoms = {}
        for atom_spec in args.atoms:
            element, count = atom_spec.split(':')
            alignment_atoms[element] = int(count)
        logger.info(f"Using specific atoms for alignment: {alignment_atoms}")

    try:
        # Perform alignment
        logger.info(f"Aligning {args.target} to {args.reference}")
        aligned_mol, rmsd, atom_mapping = align_molecules(
            args.reference,
            args.target,
            alignment_atoms
        )
        logger.info(f"Alignment complete with RMSD: {rmsd:.4f} Å")

        # Save aligned molecule
        save_aligned_molecule(aligned_mol, args.output)
        logger.info(f"Saved aligned molecule to {args.output}")

        # Generate visualization config
        vis_config = prepare_ngl_visualization(args.reference, args.output)
        with open(args.vis_config, 'w') as f:
            json.dump(vis_config, f, indent=2)
        logger.info(f"Saved visualization config to {args.vis_config}")

        # Print summary of atom mapping
        print("\nAtom Mapping Summary:")
        print(f"{'Reference':<30} {'Target':<30}")
        print("-" * 60)
        for mapping in atom_mapping[:20]:  # Show first 10 mappings
            ref_str = f"{mapping['ref_atom_name']} ({mapping['ref_element']}) {mapping['ref_residue']} {mapping['ref_residue_number']}"
            target_str = f"{mapping['target_atom_name']} ({mapping['target_element']}) {mapping['target_residue']} {mapping['target_residue_number']}"
            print(f"{ref_str:<30} {target_str:<30}")
        
        if len(atom_mapping) > 10:
            print(f"... and {len(atom_mapping) - 10} more atoms")

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        return 1

    return 0

if __name__ == "__main__":
    sys.exit(main())