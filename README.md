# Butina Molecular Clustering Overview 

A Python tool for clustering molecules from a CSV file using RDKit's implementation of the **Butina clustering algorithm**.

The script reads molecular structures supplied as SMILES, generates molecular fingerprints, calculates pairwise similarities, converts similarities to distances, and groups compounds into clusters. Optional outputs include a cluster-size bar chart, chemical structures of molecules in the largest clusters, and an Excel spreadsheet containing cluster members.



## Features

* Reads compound identifiers and SMILES from a CSV file.
* Supports two fingerprinting methods:

  * **MACCS**
  * **Morgan** (radius 2, 2048 bits as configured in the script).
* Supports two similarity coefficients:

  * **Tanimoto** similarity.
  * **Dice** similarity.
* Uses Butina clustering with a user-specified distance cutoff.
* Can optionally:

  * Plot the number of molecules in each cluster.
  * Display chemical structures from the largest clusters.
  * Export cluster members to an `.xlsx` file.



## Technologies

* Python 3.11
* RDKit
* Pandas
* Matplotlib
* Openpyxl



## Input 

A CSV file containing two columns:

* **ID**: Compound identifier.
* **SMILES**: Molecular SMILES representation.



## Output

Depending on the options selected, the script can produce:

* A bar chart showing the number of molecules in each cluster.
* Chemical structures of compounds in the largest clusters, number of clusters to examine can be specified.
* An Excel file showing the members of each cluster.



## Installation

Create the Conda environment:
conda env create -f butina\_clustering.yml



Activate it:
conda activate butina\_clustering



## Run the script

Follow the prompts to provide:

1. The path to your input CSV file.
2. The fingerprint method: `MACCS` or `Morgan`.
3. The similarity coefficient: `Tanimoto` or `Dice`.
4. The distance cutoff, between `0` and `1`.
5. Whether to generate a cluster-size bar chart (`y`/`n`).
6. Whether to display structures from the largest clusters (`y`/`n`).
7. Whether to export the cluster members to Excel (`y`/`n`).

### 

### Note: Choosing the distance cutoff

The script converts similarity to distance using:


distance = 1 - distance cutoff



Therefore, a **smaller distance cutoff is more restrictive** and groups compounds only when they are more similar. A larger cutoff permits less similar compounds to be grouped together. The appropriate value depends on the fingerprint and similarity coefficient selected.



## Citation

RDKit: Open-source cheminformatics. https://www.rdkit.org

Butina clustering: Butina, *J. Chem. Inf. Comput. Sci.*, **1999**, 747

