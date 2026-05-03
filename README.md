# Structural proximity and enrichment analysis pipeline

This repository contains two connected pipelines:

1. **Python structure-distance pipeline**

   * Reads protein structure files in **PDB** or **mmCIF** format.
   * Extracts chain sequences and CA coordinates.
   * Computes pairwise distances between protein-chain centers of mass.
   * Matches chain sequences to protein names using a FASTA reference.
   * Produces per-pair distance tables and summary statistics.

2. **R integration and correlation analysis pipeline**

   * Reads the structural distance summaries.
   * Reorders protein pairs to match the proximity-labeling naming convention.
   * Merges structural distances with proximity-labeling enrichment values.
   * Computes correlation statistics and generates the final scatter plot.

The pipeline was used to support the analyses in:

**Structural Insights into Bromodomain-Containing Complexes from *Trypanosoma cruzi* Revealed by Proximity Labeling and Stoichiometric Space Exploration** (https://www.biorxiv.org/content/10.64898/2026.03.22.713544v2)

---

## Repository contents

Recommended input folders and files:

* `cifs_af3/`

  * Input structure files used in the paper.
  * Contains the AF3-derived PDB/CIF models used for the distance analysis.

* `results_af3/`

  * Output distance tables generated from the structure files.
  * Includes `distances_averages_*.tsv` and related summary files.

* `Fasta files`

  * FASTA files containing the protein sequences used to map structure chains to protein names.

* `network.txt`

  * Proximity-labeling enrichment table used in the R analysis.

---

## Python pipeline: structure distance analysis

### What it does

The Python script scans one structure file or a directory of structure files and performs the following steps for each model:

1. Detects whether each file is PDB or mmCIF.
2. Extracts chain sequences from CA atoms.
3. Uses the CA coordinates of each chain to compute a chain center of mass.
4. Calculates all unique pairwise Euclidean distances between chain centers within each structure.
5. Matches each chain sequence to a protein name using the FASTA file.
6. Writes:

   * a **per-pair per-structure distance table**, and
   * an **aggregated summary table** with mean, standard deviation, median, and IQR.

### Inputs

* `structure_path`: a single structure file or a directory containing structure files.
* `fasta_file`: FASTA file with protein sequences and names.

### Outputs

Two TSV files are written to the selected output directory:

* `distances_<timestamp>.tsv`
* `distances_averages_<timestamp>.tsv`

### Example execution

```bash
python pdb_distance_analyzer.py cifs_af3 combined.fasta --output-dir results_af3 --recursive
```

### Command-line arguments

```text
usage: python pdb_distance_analyzer.py <structure_path> <fasta_file> [options]

Arguments:
  structure_path      Path to a PDB/CIF file or directory containing structure files
  fasta_file          Path to FASTA file with protein sequences and names

Options:
  --output-dir        Directory for output files (default: current directory)
  --recursive         Search for structure files recursively in subdirectories
```

### Notes

* The script uses **CA atoms only** as a proxy for chain position.
* The center of mass is computed as the arithmetic mean of all CA coordinates in a chain.
* If a chain sequence is not found in the FASTA file, a placeholder label such as `unknown_seq_1` is assigned.

---

## R pipeline: integration with proximity-labeling data

### What it does

The R script takes the structural distance summaries and combines them with proximity-labeling enrichment values from `network.txt`.

It performs the following steps:

1. Loads the enrichment table.
2. Computes a combined Euclidean enrichment difference from the GFP and NLS bait contrasts.
3. Loads the structural distance summary table from the Python pipeline.
4. Reorders `SOURCE` and `TARGET` so they match the naming scheme used in the enrichment dataset.
5. Joins the structural and enrichment data by protein pair.
6. Computes Spearman and Pearson correlations between structural distance and enrichment.
7. Generates and saves the correlation plot.

### Inputs

* `../bdfs_enrichment/interactors_and_proximal/network.txt`

  * Proximity-labeling enrichment table.
* `./results_af3/distances_averages_<timestamp>.tsv`

  * Summary distance table produced by the Python pipeline.

### Output

* `correlation_distance_vs_enrichment.png`

---

## Reproducing the paper analysis

To reproduce the data analysis reported in the paper, use the following resources provided in the repository:

* `cifs_af3/`: structure files used as input for the distance analysis.
* `combined.fasta`: FASTA sequences used to map structural chains to protein names.
* `results_af3/`: output distance tables generated from the structures.
* `network.txt`: proximity-labeling enrichment table used for the correlation analysis.

### Recommended reproduction workflow

1. Place the AF3-derived structure files in `cifs_af3/`.
2. Select the FASTA file containing the protein sequences for the analyzed system.
3. Run the Python script to generate chain-distance summaries.
4. Use the R script to merge the structural results with the enrichment table.
5. Generate the correlation plot and inspect the correlation statistics.

### Expected data flow

```text
cifs_af3/ + fasta_files/  -->  Python distance pipeline  -->  results_af3/distances_averages_*.tsv
results_af3/distances_averages_*.tsv + network.txt  -->  R pipeline  -->  correlation plot and statistics
```

---

## Dependencies

### Python

* Python 3
* `numpy`

### R

* R 4.x or compatible
* `tidyverse`

Install the R package with:

```r
install.packages("tidyverse")
```

Install the Python dependency with:

```bash
pip install numpy
```

---

## File format expectations

### FASTA

The FASTA file should contain protein sequences with headers that can be matched to the sequence labels used in the structures.

### Structure files

Supported formats:

* `.pdb`
* `.ent`
* `.cif`
* `.mmcif`

### Enrichment table

The enrichment file must contain at least the following columns used by the R script:

* `SOURCE`
* `TARGET`
* `N: Student's T-test Difference Bait_GFP`
* `N: Student's T-test Difference Bait_NLS`

---

## Output columns

### Python per-structure distance table

* `SOURCE`
* `TARGET`
* `seq1`
* `seq2`
* `chain1`
* `chain2`
* `distance`
* `structure_file`

### Python summary table

* `SOURCE`
* `TARGET`
* `seq1`
* `seq2`
* `mean_distance`
* `SD_distance`
* `median_distance`
* `IQR_distance`

### R merged table

* enrichment differences (`diff_GFP`, `diff_NLS`, `euclidean_diff`)
* structural distance summaries
* bait assignment (`BAIT`)

---

## Citation

If this repository is used in a publication or derivative analysis, please cite the associated paper:

**Structural Insights into Bromodomain-Containing Complexes from *Trypanosoma cruzi* Revealed by Proximity Labeling and Stoichiometric Space Exploration**
