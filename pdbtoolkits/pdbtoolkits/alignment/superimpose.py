"""
Molecular alignment tools for PDB structures.

This module provides functions for aligning PDB structures using RDKit.
"""

import logging
from rdkit import Chem
from rdkit.Chem import rdMolAlign

# Configure logger
logger = logging.getLogger(__name__)

def align_molecules(ref_pdb_path, target_pdb_path, alignment_atoms=None):
    """
    Align target molecule to reference molecule using specified atoms,
    then generate mapping information for all atoms.
    
    Parameters:
    - ref_pdb_path: Path to reference PDB file
    - target_pdb_path: Path to target PDB file to be aligned
    - alignment_atoms: Dictionary with atom symbols as keys and minimum counts as values
                       e.g., {'Mg': 1, 'N': 4} for aligning using Mg and 4 N atoms
                       
    Returns:
    - Tuple containing (aligned molecule, RMSD value, atom mapping information)
    """
    try:
        ref_mol = Chem.MolFromPDBFile(ref_pdb_path)
        target_mol = Chem.MolFromPDBFile(target_pdb_path)
        
        if not ref_mol or not target_mol:
            raise ValueError("Failed to load PDB files")

        # Prepare the alignment atom map if specific atoms are provided
        alignment_atom_map = []
        if alignment_atoms:
            for element, count in alignment_atoms.items():
                ref_indices = []
                target_indices = []
                
                # Get all atoms of this element type
                for atom in ref_mol.GetAtoms():
                    if atom.GetSymbol() == element:
                        ref_indices.append(atom.GetIdx())
                
                for atom in target_mol.GetAtoms():
                    if atom.GetSymbol() == element:
                        target_indices.append(atom.GetIdx())
                
                if len(ref_indices) < count or len(target_indices) < count:
                    raise ValueError(f"Not enough {element} atoms (need {count}, ref has {len(ref_indices)}, target has {len(target_indices)})")
                
                # Map the specified number of atoms for alignment
                for i in range(min(count, len(ref_indices), len(target_indices))):
                    alignment_atom_map.append((target_indices[i], ref_indices[i]))
        
        # Perform alignment using the specified atoms or all heavy atoms if none specified
        if alignment_atoms:
            rmsd = rdMolAlign.AlignMol(target_mol, ref_mol, atomMap=alignment_atom_map)
        else:
            rmsd = rdMolAlign.AlignMol(target_mol, ref_mol)
        
        # Generate comprehensive atom mapping (all atoms)
        complete_atom_mapping = []
        
        # Process reference molecule atoms
        ref_atoms = list(ref_mol.GetAtoms())
        for ref_atom in ref_atoms:
            ref_info = ref_atom.GetPDBResidueInfo()
            if ref_info:
                ref_mapping = {
                    'ref_atom_name': ref_info.GetName().strip(),
                    'ref_element': ref_atom.GetSymbol(),
                    'ref_residue': ref_info.GetResidueName(),
                    'ref_residue_number': ref_info.GetResidueNumber(),
                    'ref_chain': ref_info.GetChainId(),
                    'target_atom_name': '',
                    'target_element': '',
                    'target_residue': '',
                    'target_residue_number': '',
                    'target_chain': ''
                }
                
                # Find corresponding atom in target molecule by element and name similarity
                best_match = None
                for target_atom in target_mol.GetAtoms():
                    target_info = target_atom.GetPDBResidueInfo()
                    if target_info and target_atom.GetSymbol() == ref_atom.GetSymbol():
                        # If we have an exact name match, use that
                        if target_info.GetName().strip() == ref_info.GetName().strip():
                            best_match = target_atom
                            break
                        # Otherwise keep track of possible matches
                        elif best_match is None:
                            best_match = target_atom
                
                # If we found a match, add its details
                if best_match:
                    target_info = best_match.GetPDBResidueInfo()
                    ref_mapping.update({
                        'target_atom_name': target_info.GetName().strip(),
                        'target_element': best_match.GetSymbol(),
                        'target_residue': target_info.GetResidueName(),
                        'target_residue_number': target_info.GetResidueNumber(),
                        'target_chain': target_info.GetChainId()
                    })
                
                complete_atom_mapping.append(ref_mapping)
        
        return target_mol, rmsd, complete_atom_mapping
        
    except Exception as e:
        logger.error(f"Alignment failed: {str(e)}")
        raise

def save_aligned_molecule(mol, output_path):
    """
    Save the aligned molecule to a PDB file.
    
    Parameters:
    - mol: RDKit molecule object
    - output_path: Path to save the aligned PDB file
    
    Returns:
    - Path to the saved file
    """
    try:
        Chem.MolToPDBFile(mol, output_path)
        return output_path
    except Exception as e:
        logger.error(f"Failed to save aligned molecule: {str(e)}")
        raise