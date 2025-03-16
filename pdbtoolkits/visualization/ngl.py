"""
NGL-based visualization tools for PDB structures.

This module provides functions for generating visualization configurations
for the NGL viewer to display PDB structures.
"""

import logging

# Configure logger
logger = logging.getLogger(__name__)

def prepare_ngl_visualization(ref_pdb_path, aligned_pdb_path):
    """
    Generate NGL visualization settings for reference and aligned molecules.
    
    Parameters:
    - ref_pdb_path: Path to reference PDB file
    - aligned_pdb_path: Path to aligned PDB file
    
    Returns:
    - Dictionary with NGL viewer settings
    """
    try:
        # Create NGL viewer configuration
        ngl_config = {
            'files': [
                {
                    'path': ref_pdb_path,
                    'display_name': 'Reference',
                    'representation': [
                        {
                            'type': 'licorice',
                            'params': {
                                'colorScheme': 'element',
                                'colorValue': '#2E86C1',
                                'radius': 0.3
                            }
                        }
                    ]
                },
                {
                    'path': aligned_pdb_path,
                    'display_name': 'Aligned',
                    'representation': [
                        {
                            'type': 'licorice',
                            'params': {
                                'colorScheme': 'element',
                                'colorValue': '#E74C3C',
                                'radius': 0.3
                            }
                        }
                    ]
                }
            ],
            'stage_params': {
                'backgroundColor': 'white'
            }
        }
        
        return ngl_config
        
    except Exception as e:
        logger.error(f"NGL visualization preparation failed: {str(e)}")
        raise

def create_comparison_view(molecules, labels=None, colors=None, representations=None):
    """
    Create a more customizable comparison view for multiple molecules.
    
    Parameters:
    - molecules: List of paths to PDB files
    - labels: List of display names for each molecule
    - colors: List of colors for each molecule
    - representations: List of representation types ('cartoon', 'licorice', 'ball+stick', etc.)
    
    Returns:
    - Dictionary with NGL viewer settings
    """
    if labels is None:
        labels = [f"Molecule {i+1}" for i in range(len(molecules))]
    
    if colors is None:
        # Default color palette
        default_colors = ['#2E86C1', '#E74C3C', '#2ECC71', '#F1C40F', '#9B59B6']
        colors = [default_colors[i % len(default_colors)] for i in range(len(molecules))]
    
    if representations is None:
        representations = ['licorice'] * len(molecules)
    
    files_config = []
    for i, molecule_path in enumerate(molecules):
        files_config.append({
            'path': molecule_path,
            'display_name': labels[i],
            'representation': [
                {
                    'type': representations[i],
                    'params': {
                        'colorScheme': 'element',
                        'colorValue': colors[i],
                        'radius': 0.3
                    }
                }
            ]
        })
    
    ngl_config = {
        'files': files_config,
        'stage_params': {
            'backgroundColor': 'white'
        }
    }
    
    return ngl_config