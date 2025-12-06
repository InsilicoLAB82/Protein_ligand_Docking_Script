# **Protein-Ligand Docking Script by InsilicoLAB**
<img width="3000" height="1800" alt="Vina_Results" src="https://github.com/user-attachments/assets/0ff0d030-e85d-4e2c-84e4-4e3a889c0407" />

# **Overview:**

This Python script automates the process of protein-ligand docking using AutoDock Vina and calculates the RMSD (Root Mean Square Deviation) of the docking poses compared to a reference ligand using LS-align. It generates docking results, including binding affinities and RMSD values, and produces a visualization of the docking results.

This script is useful for computational drug design and molecular docking studies, where it helps in predicting the interaction of a ligand with a receptor.

**Prerequisites:**

Before running this script, you need the following software and files:

AutoDock Vina executable (vina_1.2.7_linux_x86_64)

LS-align executable

Open Babel (for converting molecular file formats)

PDB and PDBQT files for the receptor and ligand(s)

Grid parameters file (txt format) for docking setup

Reference ligand file in PDBQT format for RMSD calculation

**Requirements:**

Python (>=3.x)

Subprocess module: For running external commands such as AutoDock Vina and LS-align.

Matplotlib and NumPy: For plotting the docking results.

Open Babel: For converting receptor and ligand files into the PDBQT format.

**Usage:**

Upload Required Files:

The following files need to be uploaded for the docking process:

Receptor File (in PDB format)

Ligand File (in PDB format)

Grid Parameters File (txt format containing grid center and size)

Reference Ligand File (in PDBQT format)

**Provide Docking Parameters:**

When prompted, input the desired docking parameters such as:

Exhaustiveness: Determines the search quality (default: 32).

Scoring Function: Choose between vina, vinardo, or ad4 (default: vina).

Number of Poses: How many poses to generate (default: 20).

Energy Range: The energy window for the docking poses (default: 3.0 kcal/mol).

Grid Spacing: The spacing of the grid in Angstroms (default: 0.375 Å).

RMSD Threshold: A value to filter poses based on RMSD (default: 999, for all poses).

Convert Files to PDBQT Format:

The receptor, ligand, and reference ligand files are converted to PDBQT format using Open Babel.

Run AutoDock Vina:

AutoDock Vina is executed with the provided parameters, and the docking process is carried out.

The script prints real-time output from the Vina docking process.

Calculate RMSD Using LS-align:

For each docking pose, the RMSD is calculated by comparing the ligand pose to the reference ligand using LS-align.

The RMSD values and binding affinities are stored for further analysis.

Results Visualization:

A scatter plot is generated, showing the RMSD (Å) versus the binding affinity (kcal/mol).

The best docking pose is marked on the plot with a red star, and the plot is saved as Vina_Results.png.

# **Final Output:**

The best docking pose (with the lowest RMSD and highest binding affinity) is saved in the final PDBQT format as result_docking_final.pdbqt.

Detailed results are printed, including the best pose's affinity and RMSD.
# +++++++++++++++++++++++++++++++++++++
# **Protein-Ligand Docking Pipeline for Google Colab**
# +++++++++++++++++++++++++++++++++++++
<img width="5253" height="3680" alt="docking_results" src="https://github.com/user-attachments/assets/9e3a7e2d-50d1-49d6-9ed1-b11e20f52ac8" />
# Overview
This repository contains a comprehensive Jupyter Notebook for automated protein-ligand docking using AutoDock Vina in Google Colab. Designed by InsilicoLAB, this pipeline streamlines molecular docking workflows with minimal setup requirements, making computational drug discovery accessible to researchers, students, and enthusiasts.
# Features
Automated File Processing: Handles multiple molecular formats (PDB, PDBQT, MOL2, SDF) with automatic conversion
Flexible Grid Configuration: Customizable binding site definition via parameter files or manual input
Multiple Scoring Functions: Supports Vina, Vinardo, and AD4 scoring functions
RMSD Analysis: Calculates structural similarity against reference ligands
Interactive Visualization: Generates comprehensive plots and summary statistics
Colab-Optimized: Runs entirely in Google Colab with automatic dependency installation
User-Friendly Interface: Step-by-step workflow with interactive parameter configuration
# Workflow
# Step 1: Environment Setup
The notebook automatically installs all required dependencies:
AutoDock Vina (v1.2.7)
OpenBabel (for file conversion)
Python libraries (Pandas, Matplotlib, NumPy, Seaborn, RDKit)
LS-align (for RMSD calculations)

# Step 2: Input Preparation
Upload your files through the interactive interface:
Receptor: Protein structure in PDB format
Ligand: Small molecule in PDB/MOL2/SDF format
Grid Parameters: Text file defining the binding site (optional - can be created manually)
Reference Ligand: For RMSD calculations (optional)

# Step 3: Docking Configuration
Customize docking parameters:
Grid Center & Size: Define the search space
Exhaustiveness: Search thoroughness (8-64)
Scoring Function: Choose between vina, vinardo, or ad4
Number of Poses: Output poses to generate
Energy Range: Maximum energy difference between poses

# Step 4: Run Docking
Execute the pipeline with one click:
File format conversion to PDBQT
AutoDock Vina execution
Pose extraction and scoring
RMSD calculation (if reference provided)

# Step 5: Results Visualization
Automatically generated outputs include:
Affinity vs RMSD Scatter Plot
Affinity Distribution Histogram
Pose Ranking Bar Chart
Interactive Summary Dashboard
CSV Export with detailed pose information

# Input File Formats
Supported Formats:
Receptor: .pdb, .pdbqt
Ligand: .pdb, .pdbqt, .mol2, .sdf
Grid Parameters: .txt with key-value pairs
