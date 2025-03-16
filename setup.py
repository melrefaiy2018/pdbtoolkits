# setup.py
from setuptools import setup, find_packages

setup(
    name="pdbtoolkits",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "rdkit>=2020.09.1",
        "pytest>=6.0.0",
    ],
    author="Mohamed Elrefaiy",
    author_email="moelrefaiy@gmail.com",
    description="A toolkit for working with PDB files",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/pdb-toolKits",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
)