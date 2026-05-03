#!/usr/bin/env python3
"""
PDB/CIF Chain Distance Analyzer

This script computes distances between protein chains' centers of mass in PDB or CIF files.
It matches sequences to protein names using a FASTA file and generates detailed
distance tables for proximity labeling analysis.

Usage:
    python pdb_distance_analyzer.py <structure_path> <fasta_file> [options]

Arguments:
    structure_path: Path to a PDB/CIF file or directory containing structure files
    fasta_file: Path to FASTA file with protein sequences and names

Options:
    --output-dir: Directory for output files (default: current directory)
    --recursive: Search for structure files recursively in subdirectories
"""

import os
import sys
import argparse
import re
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import numpy as np


class StructureParser:
    """Base parser for extracting chain sequences and coordinates from structure files."""
    
    def __init__(self, structure_file):
        self.structure_file = structure_file
        self.chains = {}
        
    def parse(self):
        """Extract chain sequences and CA atom coordinates. Override in subclasses."""
        raise NotImplementedError
    
    def compute_center_of_mass(self, chain_id):
        """Compute the center of mass for a chain using CA coordinates."""
        if chain_id not in self.chains:
            return None
        coords = self.chains[chain_id]['coordinates']
        return np.mean(coords, axis=0)


class PDBParser(StructureParser):
    """Parser for extracting chain sequences and coordinates from PDB files."""
    
    def parse(self):
        """Extract chain sequences and CA atom coordinates from PDB file."""
        chain_coords = defaultdict(list)
        chain_seqs = defaultdict(list)
        chain_residues = defaultdict(set)
        
        # Three-letter to one-letter amino acid code mapping
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        
        try:
            with open(self.structure_file, 'r') as f:
                for line in f:
                    # Parse ATOM records for CA atoms (center of mass approximation)
                    if line.startswith('ATOM'):
                        atom_name = line[12:16].strip()
                        if atom_name == 'CA':  # Alpha carbon
                            chain_id = line[21]
                            residue_name = line[17:20].strip()
                            residue_num = int(line[22:26].strip())
                            
                            # Extract coordinates
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            
                            chain_coords[chain_id].append([x, y, z])
                            
                            # Build sequence (only add each residue once)
                            if residue_num not in chain_residues[chain_id]:
                                chain_residues[chain_id].add(residue_num)
                                if residue_name in aa_map:
                                    chain_seqs[chain_id].append(aa_map[residue_name])
        
        except Exception as e:
            print(f"Error parsing {self.structure_file}: {e}", file=sys.stderr)
            return {}
        
        # Convert to final format
        for chain_id in chain_coords:
            self.chains[chain_id] = {
                'sequence': ''.join(chain_seqs[chain_id]),
                'coordinates': np.array(chain_coords[chain_id])
            }
        
        return self.chains


class CIFParser(StructureParser):
    """Parser for extracting chain sequences and coordinates from mmCIF files."""
    
    def parse(self):
        """Extract chain sequences and CA atom coordinates from CIF file."""
        chain_coords = defaultdict(list)
        chain_seqs = defaultdict(list)
        chain_residues = defaultdict(set)
        
        # Three-letter to one-letter amino acid code mapping
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        
        try:
            with open(self.structure_file, 'r') as f:
                in_atom_site = False
                col_indices = {}
                
                for line in f:
                    line = line.strip()
                    
                    # Check if we're entering the atom_site loop
                    if line.startswith('_atom_site.'):
                        in_atom_site = True
                        # Parse column name
                        col_name = line.split('.')[1].split()[0]
                        col_idx = len(col_indices)
                        col_indices[col_name] = col_idx
                        continue
                    
                    # Check if we're leaving the atom_site section
                    if in_atom_site and (line.startswith('#') or line.startswith('_')):
                        in_atom_site = False
                        continue
                    
                    # Parse atom data
                    if in_atom_site and line and not line.startswith('_'):
                        # Split by whitespace, handling quoted strings
                        parts = []
                        in_quote = False
                        current = []
                        
                        for char in line:
                            if char in ('"', "'"):
                                in_quote = not in_quote
                            elif char.isspace() and not in_quote:
                                if current:
                                    parts.append(''.join(current))
                                    current = []
                            else:
                                current.append(char)
                        if current:
                            parts.append(''.join(current))
                        
                        if len(parts) < len(col_indices):
                            continue
                        
                        # Extract relevant fields
                        try:
                            group_PDB = parts[col_indices.get('group_PDB', 0)]
                            atom_name = parts[col_indices.get('label_atom_id', 2)]
                            residue_name = parts[col_indices.get('label_comp_id', 5)]
                            chain_id = parts[col_indices.get('label_asym_id', 6)]
                            residue_num = int(parts[col_indices.get('label_seq_id', 8)])
                            x = float(parts[col_indices.get('Cartn_x', 10)])
                            y = float(parts[col_indices.get('Cartn_y', 11)])
                            z = float(parts[col_indices.get('Cartn_z', 12)])
                        except (KeyError, ValueError, IndexError):
                            continue
                        
                        # Only process ATOM records with CA atoms
                        if group_PDB == 'ATOM' and atom_name == 'CA':
                            chain_coords[chain_id].append([x, y, z])
                            
                            # Build sequence (only add each residue once)
                            if residue_num not in chain_residues[chain_id]:
                                chain_residues[chain_id].add(residue_num)
                                if residue_name in aa_map:
                                    chain_seqs[chain_id].append(aa_map[residue_name])
        
        except Exception as e:
            print(f"Error parsing {self.structure_file}: {e}", file=sys.stderr)
            return {}
        
        # Convert to final format
        for chain_id in chain_coords:
            self.chains[chain_id] = {
                'sequence': ''.join(chain_seqs[chain_id]),
                'coordinates': np.array(chain_coords[chain_id])
            }
        
        return self.chains


def get_parser(structure_file):
    """Return appropriate parser based on file extension."""
    suffix = Path(structure_file).suffix.lower()
    if suffix in ['.cif', '.mmcif']:
        return CIFParser(structure_file)
    elif suffix in ['.pdb', '.ent']:
        return PDBParser(structure_file)
    else:
        print(f"Warning: Unknown file type {suffix}, attempting PDB format", file=sys.stderr)
        return PDBParser(structure_file)


class FASTAParser:
    """Parser for FASTA files to map sequences to protein names."""
    
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.seq_to_name = {}
        
    def parse(self):
        """Parse FASTA file and create sequence-to-name mapping."""
        try:
            with open(self.fasta_file, 'r') as f:
                current_name = None
                current_seq = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence
                        if current_name and current_seq:
                            seq = ''.join(current_seq)
                            self.seq_to_name[seq] = current_name
                        
                        # Extract protein name (first word after '>')
                        # Stop at space, comma, semicolon, or pipe
                        header = line[1:]  # Remove '>'
                        match = re.match(r'^([^\s,;|]+)', header)
                        current_name = match.group(1) if match else header
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Don't forget the last sequence
                if current_name and current_seq:
                    seq = ''.join(current_seq)
                    self.seq_to_name[seq] = current_name
        
        except Exception as e:
            print(f"Error parsing FASTA file {self.fasta_file}: {e}", file=sys.stderr)
            return {}
        
        return self.seq_to_name


class DistanceAnalyzer:
    """Analyzer for computing and aggregating chain distances."""
    
    def __init__(self, fasta_parser):
        self.seq_to_name = fasta_parser.seq_to_name
        self.unknown_counter = 0
        self.unknown_seq_map = {}
        self.all_distances = []
        
    def get_protein_name(self, sequence):
        """Get protein name for a sequence, assign unknown_seq_i if not found."""
        if sequence in self.seq_to_name:
            return self.seq_to_name[sequence]
        
        # Check if we've already assigned an unknown name to this sequence
        if sequence in self.unknown_seq_map:
            return self.unknown_seq_map[sequence]
        
        # Assign new unknown identifier
        self.unknown_counter += 1
        unknown_name = f"unknown_seq_{self.unknown_counter}"
        self.unknown_seq_map[sequence] = unknown_name
        return unknown_name
    
    def analyze_structure(self, structure_file):
        """Analyze a single structure file and compute all pairwise chain distances."""
        parser = get_parser(structure_file)
        chains = parser.parse()
        
        if not chains:
            print(f"Warning: No chains found in {structure_file}", file=sys.stderr)
            return
        
        chain_ids = sorted(chains.keys())
        
        # Compute pairwise distances (only once per pair)
        for i, chain1 in enumerate(chain_ids):
            for chain2 in chain_ids[i+1:]:
                seq1 = chains[chain1]['sequence']
                seq2 = chains[chain2]['sequence']
                
                # Skip if sequences are empty
                if not seq1 or not seq2:
                    continue
                
                # Compute centers of mass
                com1 = parser.compute_center_of_mass(chain1)
                com2 = parser.compute_center_of_mass(chain2)
                
                # Compute Euclidean distance
                distance = np.linalg.norm(com1 - com2)
                
                # Get protein names
                protein1 = self.get_protein_name(seq1)
                protein2 = self.get_protein_name(seq2)
                
                # Store result
                self.all_distances.append({
                    'SOURCE': protein1,
                    'TARGET': protein2,
                    'seq1': seq1,
                    'seq2': seq2,
                    'chain1': chain1,
                    'chain2': chain2,
                    'distance': distance,
                    'structure_file': str(structure_file)
                })
    
    def compute_averages(self):
        """Compute average distances for each protein pair across all structure files."""
        # Group by protein pair
        pair_distances = defaultdict(list)
        pair_seqs = {}
        
        for record in self.all_distances:
            # Ensure consistent ordering (alphabetical by protein name)
            proteins = sorted([record['SOURCE'], record['TARGET']])
            pair_key = tuple(proteins)
            
            pair_distances[pair_key].append(record['distance'])
            
            # Store sequences (maintain consistent ordering)
            if pair_key not in pair_seqs:
                if proteins[0] == record['SOURCE']:
                    pair_seqs[pair_key] = (record['seq1'], record['seq2'])
                else:
                    pair_seqs[pair_key] = (record['seq2'], record['seq1'])
        
        # Compute statistics
        averages = []
        for pair_key, distances in pair_distances.items():
            distances_array = np.array(distances)
            seq1, seq2 = pair_seqs[pair_key]
            
            averages.append({
                'SOURCE': pair_key[0],
                'TARGET': pair_key[1],
                'seq1': seq1,
                'seq2': seq2,
                'mean_distance': np.mean(distances_array),
                'SD_distance': np.std(distances_array, ddof=1) if len(distances_array) > 1 else 0.0,
                'median_distance': np.median(distances_array),
                'IQR_distance': np.percentile(distances_array, 75) - np.percentile(distances_array, 25)
            })
        
        return averages


def find_structure_files(path, recursive=True):
    """Find all PDB/CIF files in a directory or return single file."""
    path = Path(path)
    
    if path.is_file():
        if path.suffix.lower() in ['.pdb', '.ent', '.cif', '.mmcif']:
            return [path]
        else:
            print(f"Error: {path} is not a PDB or CIF file", file=sys.stderr)
            return []
    
    if path.is_dir():
        if recursive:
            structure_files = (list(path.rglob('*.pdb')) + list(path.rglob('*.ent')) + 
                             list(path.rglob('*.cif')) + list(path.rglob('*.mmcif')))
        else:
            structure_files = (list(path.glob('*.pdb')) + list(path.glob('*.ent')) +
                             list(path.glob('*.cif')) + list(path.glob('*.mmcif')))
        return structure_files
    
    print(f"Error: {path} does not exist", file=sys.stderr)
    return []


def write_tsv(data, filename, columns):
    """Write data to TSV file."""
    with open(filename, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')
        
        # Write data rows
        for record in data:
            row = [str(record[col]) for col in columns]
            f.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Compute distances between protein chains in PDB/CIF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('structure_path', help='Path to PDB/CIF file or directory')
    parser.add_argument('fasta_file', help='Path to FASTA file with protein sequences')
    parser.add_argument('--output-dir', default='.', help='Output directory (default: current directory)')
    parser.add_argument('--recursive', action='store_true', 
                       help='Search for structure files recursively in subdirectories')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.structure_path):
        print(f"Error: Structure path '{args.structure_path}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.fasta_file):
        print(f"Error: FASTA file '{args.fasta_file}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse FASTA file
    print("Parsing FASTA file...")
    fasta_parser = FASTAParser(args.fasta_file)
    fasta_parser.parse()
    print(f"Found {len(fasta_parser.seq_to_name)} protein sequences in FASTA file")
    
    # Find structure files
    print("\nSearching for structure files...")
    structure_files = find_structure_files(args.structure_path, recursive=args.recursive)
    print(f"Found {len(structure_files)} structure files")
    
    if not structure_files:
        print("No structure files found. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    # Analyze all structure files
    print("\nAnalyzing structure files...")
    analyzer = DistanceAnalyzer(fasta_parser)
    
    for i, structure_file in enumerate(structure_files, 1):
        print(f"Processing {i}/{len(structure_files)}: {structure_file.name}")
        analyzer.analyze_structure(structure_file)
    
    print(f"\nTotal pairwise distances computed: {len(analyzer.all_distances)}")
    
    # Generate output filenames with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    distances_file = output_dir / f"distances_{timestamp}.tsv"
    averages_file = output_dir / f"distances_averages_{timestamp}.tsv"
    
    # Write individual distances
    print(f"\nWriting individual distances to {distances_file}")
    distance_columns = ['SOURCE', 'TARGET', 'seq1', 'seq2', 'chain1', 'chain2', 'distance', 'structure_file']
    write_tsv(analyzer.all_distances, distances_file, distance_columns)
    
    # Compute and write averages
    print(f"Computing average distances...")
    averages = analyzer.compute_averages()
    print(f"Writing average distances to {averages_file}")
    average_columns = ['SOURCE', 'TARGET', 'seq1', 'seq2', 'mean_distance', 
                      'SD_distance', 'median_distance', 'IQR_distance']
    write_tsv(averages, averages_file, average_columns)
    
    print(f"\nAnalysis complete!")
    print(f"Total protein pairs analyzed: {len(averages)}")
    print(f"Unknown sequences assigned: {analyzer.unknown_counter}")


if __name__ == '__main__':
    main()