# Gibbs Free Energy Calculator for RNA Structures

This tool calculates the Gibbs free energy for RNA structures based on provided PDB files and precomputed score files.
This code is designed to work specifically with PDB (Protein Data Bank) files that contain RNA (Ribonucleic Acid) structures.

## Installation

```bash
git clone https://github.com/usama-ak/RNA-Folding-Problem.git
```

## Usage

Navigate to the directory:

    ```bash
    cd RNA-Folding-Project
    ```

Run the main script:

    ```bash
    python gibbs_estimation.py --pdb-file path/to/your/pdb/file.pdb --scores-dir path/to/score/files/
    ```

    The `--pdb-file` argument specifies the path to the PDB file containing the RNA structure for which you want to calculate the Gibbs free energy.  
    You can also specify the `--scores-dir` argument to provide a custom directory containing score files. If not specified, it defaults to `data/scores/`.

Alternatively, if you don't specify a PDB file, the script will process example files located in the data/examples/ directory:

    ```bash
    python gibbs_estimation.py
    ```

## File Structure

- `src/`: Directory containing utility modules
    - `utils/`: Directory containing utility functions for distance calculation, file I/O, and score calculation
- `data/`: Directory containing data files
    - `examples/`: Directory containing example files
    - `plots/`: Directory containing plots
    - `scores/`: Directory containing score files
    - `train-data`: Directory containing PDB files for training the objective function

## Scripts 

- `src/plot_scoring_profiles.r`: R script for plotting scoring profiles
- `src/obj_funct_training.py`: Python script for training the objective function, generating scores files
- `gibbs_estimation.py`: Main python script for calculating Gibbs free energy of an RNA structure 


## Authors

- [Usama AKHTAR](https://github.com/usama-ak)
- [Ali YOUNCHA](https://github.com/MrAli1582)

