# Visualization Module

The visualization module provides tools for generating configuration settings for NGL Viewer to visualize PDB structures.

## Overview

Visualization is an essential component of structural biology and molecular analysis. The visualization module helps generate configuration for the NGL Viewer, a powerful WebGL-based molecular visualization tool.

## Functions

### `prepare_ngl_visualization`

```python
from pdbtoolkits.visualization.ngl import prepare_ngl_visualization

ngl_config = prepare_ngl_visualization(ref_pdb_path, aligned_pdb_path)
```

Generates NGL visualization settings for reference and aligned molecules.

**Parameters:**
- `ref_pdb_path` (str): Path to the reference PDB file
- `aligned_pdb_path` (str): Path to the aligned PDB file

**Returns:**
- Dictionary with NGL viewer settings that can be used in web applications

**Example:**
```python
# Generate config for a reference and aligned structure
config = prepare_ngl_visualization("reference.pdb", "aligned_target.pdb")

# Save configuration to a JSON file
import json
with open("ngl_config.json", "w") as f:
    json.dump(config, f, indent=2)
```

### `create_comparison_view`

```python
from pdbtoolkits.visualization.ngl import create_comparison_view

ngl_config = create_comparison_view(molecules, labels=None, colors=None, representations=None)
```

Creates a more customizable comparison view for multiple molecules.

**Parameters:**
- `molecules` (list): List of paths to PDB files
- `labels` (list, optional): List of display names for each molecule
- `colors` (list, optional): List of colors for each molecule
- `representations` (list, optional): List of representation types

**Returns:**
- Dictionary with NGL viewer settings for a multi-molecule comparison

**Example:**
```python
# Create a comparison view of multiple conformations
molecules = ["conf1.pdb", "conf2.pdb", "conf3.pdb"]
labels = ["State A", "State B", "State C"]
colors = ["#2E86C1", "#E74C3C", "#2ECC71"]
representations = ["cartoon", "licorice", "ball+stick"]

config = create_comparison_view(molecules, labels, colors, representations)
```

## Configuration Structure

The generated NGL configuration has the following structure:

```json
{
  "files": [
    {
      "path": "reference.pdb",
      "display_name": "Reference",
      "representation": [
        {
          "type": "licorice",
          "params": {
            "colorScheme": "element",
            "colorValue": "#2E86C1",
            "radius": 0.3
          }
        }
      ]
    },
    {
      "path": "aligned.pdb",
      "display_name": "Aligned",
      "representation": [
        {
          "type": "licorice",
          "params": {
            "colorScheme": "element",
            "colorValue": "#E74C3C",
            "radius": 0.3
          }
        }
      ]
    }
  ],
  "stage_params": {
    "backgroundColor": "white"
  }
}
```

## Available Representation Types

The visualization module supports various representation types for molecular visualization:

- `"licorice"`: Stick representation showing bonds
- `"cartoon"`: Ribbon representation for proteins
- `"ball+stick"`: Atoms as balls and bonds as sticks
- `"spacefill"`: Van der Waals spheres
- `"surface"`: Molecular surface representation
- `"line"`: Simple line representation

## Example Web Integration

The generated configuration can be used with NGL Viewer in a web application. Here's a simple example of HTML/JavaScript code:

```html
<!DOCTYPE html>
<html>
<head>
  <title>PDB Viewer</title>
  <script src="https://unpkg.com/ngl@2.0.0-dev.37/dist/ngl.js"></script>
  <script src="ngl_config.json" type="application/json" id="ngl-config"></script>
  <style>
    #viewport { width: 800px; height: 600px; }
  </style>
</head>
<body>
  <div id="viewport"></div>
  <script>
    // Load configuration
    const config = JSON.parse(document.getElementById('ngl-config').textContent);
    
    // Create NGL Stage
    const stage = new NGL.Stage('viewport', { backgroundColor: config.stage_params.backgroundColor });
    
    // Load and display molecules
    config.files.forEach(fileConfig => {
      stage.loadFile(fileConfig.path).then(component => {
        // Add representations
        fileConfig.representation.forEach(rep => {
          component.addRepresentation(rep.type, rep.params);
        });
        
        // Set component name
        component.name = fileConfig.display_name;
        
        // Center view on all structures
        stage.autoView();
      });
    });
  </script>
</body>
</html>
```

## Future Enhancements

Planned enhancements for the visualization module include:

1. **Interactive interface** for manipulating visualization settings
2. **Additional representation types** for specialized visualization needs
3. **Animation capabilities** for showing conformational changes
4. **Selective visualization** of specific parts of structures
5. **Measurement tools** for distances, angles, and dihedrals