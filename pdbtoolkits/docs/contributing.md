# Contributing to pdb-toolKits

Thank you for your interest in contributing to pdb-toolKits! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

Please be respectful and considerate of others when contributing to this project. We aim to foster an inclusive and welcoming community.

## Getting Started

### Prerequisites

- Python 3.6 or higher
- Git
- RDKit and other dependencies listed in requirements.txt

### Development Setup

1. Fork the repository on GitHub
2. Clone your fork locally
   ```bash
   git clone https://github.com/yourusername/pdb-toolKits.git
   cd pdb-toolKits
   ```
3. Create a branch for your feature or bugfix
   ```bash
   git checkout -b feature/my-new-feature
   ```
4. Install development dependencies
   ```bash
   pip install -e ".[dev]"
   ```

## Development Guidelines

### Code Style

We follow the PEP 8 style guide for Python code. Please ensure your code adheres to these standards.

- Use 4 spaces for indentation
- Maximum line length of 88 characters (follows Black formatter)
- Use descriptive variable and function names
- Include docstrings for all functions, classes, and modules

### Documentation

All new functions, classes, and modules should be properly documented:

- Include docstrings in Google or NumPy format
- Update relevant documentation files in the `docs` directory
- Add example usage where appropriate

### Testing

All new features or bug fixes should include tests:

1. Write tests that cover your new functionality
2. Make sure all tests pass before submitting a pull request
3. Run the tests with:
   ```bash
   pytest
   ```

### Commit Messages

Write clear, concise commit messages that explain what changes were made and why:

- Use the present tense ("Add feature" not "Added feature")
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
- Reference issue numbers if applicable

## Adding New Functionality

### Adding a New Module

To add a new module to pdb-toolKits:

1. Create a new subdirectory under `pdbtoolkits/`
2. Add an `__init__.py` file to expose the module's functionality
3. Implement your functionality in appropriate Python files
4. Update documentation to describe the new module
5. Add example scripts to demonstrate usage
6. Add tests for the new functionality

### Example Module Structure

```
pdbtoolkits/
└── new_module/
    ├── __init__.py           # Expose functionality
    ├── core.py               # Core implementation
    └── utilities.py          # Helper functions
```

## Pull Request Process

1. Update the documentation with details of your changes
2. Update the README.md if appropriate
3. Make sure all tests pass
4. Submit a pull request to the main repository
5. The pull request will be reviewed by maintainers, who may request changes

## Feature Requests and Bug Reports

If you have ideas for features or have found a bug:

1. Check if the feature/bug is already in the issues
2. If not, create a new issue with a clear description
3. For bugs, include steps to reproduce, expected behavior, and actual behavior
4. For feature requests, describe what you want to accomplish and why it's valuable

## Repository Structure

Understanding the project structure helps when contributing:

```
pdb-toolKits/
├── pdbtoolkits/            # Main package
│   ├── alignment/          # Alignment module
│   └── visualization/      # Visualization module
├── docs/                   # Documentation
├── Examples/               # Example scripts
├── Testing/                # Test directory
├── setup.py                # Package setup
├── requirements.txt        # Dependencies
└── README.md               # Main README
```

## License

By contributing to pdb-toolKits, you agree that your contributions will be licensed under the project's MIT License.