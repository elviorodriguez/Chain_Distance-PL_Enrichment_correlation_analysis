This repository provides a lightweight Python-based pipeline to quantify spatial relationships between protein chains from structural models and integrate them with downstream experimental data. The tool was developed to support proximity labeling and interactome analyses by enabling systematic comparison between structural proximity and biochemical enrichment.

Overview

The pipeline parses protein structure files in PDB or mmCIF format, reconstructs chain sequences, and computes distances between chains based on their centers of mass (approximated from Cα atoms). Protein identities are assigned by matching reconstructed sequences to a user-provided FASTA reference. The output consists of both individual pairwise distances and aggregated statistics across multiple structures.

An accompanying R workflow is included to integrate these structural measurements with enrichment datasets and evaluate correlations between spatial proximity and experimental signals.

Key Features
Supports both PDB and mmCIF structure formats
Extracts chain sequences and coordinates directly from structure files
Computes center-of-mass distances between all chain pairs
Maps chains to protein identities using a FASTA reference
Handles unmatched sequences with automatic labeling
Aggregates results across multiple structures, reporting:
Mean distance
Standard deviation
Median
Interquartile range
Provides R scripts for:
Data integration with enrichment datasets
Pair reformatting and matching
Correlation analysis and visualization
Typical Use Case
Predict or obtain protein complex structures (e.g., from AlphaFold-based pipelines).
Run the analyzer on a directory of structure files.
Obtain a table of pairwise distances and summary statistics.
Integrate with experimental data (e.g., proximity labeling, proteomics).
Assess whether structural proximity correlates with enrichment or interaction strength.
Usage
python pdb_distance_analyzer.py <structure_path> <fasta_file> --output-dir results --recursive
Outputs
distances_<timestamp>.tsv: All pairwise chain distances
distances_averages_<timestamp>.tsv: Aggregated statistics per protein pair
Dependencies
Python ≥ 3.8
NumPy
R (for downstream analysis) with:
tidyverse
ggplot2
Notes
Distances are computed using Cα atoms as a proxy for chain centers, which provides a robust and efficient approximation for large-scale analyses.
Sequence matching requires exact identity; mismatches will be labeled as unknown sequences.
The pipeline is designed to scale to large structural datasets generated from combinatorial modeling approaches.
