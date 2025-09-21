#!/usr/bin/env python3

from sys import argv
from GTF import parse as parse_gtf_line
from os.path import basename, exists
from collections import defaultdict
import os

def split_gtf_by_gene_id(input_gtf, output_prefix=None, genes_per_file=100):
    """
    Split GTF file by gene_id, handling unsorted gene_ids.
    Groups genes into files with specified number of genes per file.
    """
    if output_prefix is None:
        output_prefix = basename(input_gtf).split('.')[0]

    # First pass: collect all lines for each gene_id
    gene_lines = defaultdict(list)
    header_lines = []

    print(f"Reading GTF file: {input_gtf}")

    with open(input_gtf, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                header_lines.append(line)
                continue

            try:
                fields = parse_gtf_line(line)
                gene_id = fields.get('gene_id')

                if gene_id is None:
                    print(f"Warning: Could not parse gene_id from line {line_num}")
                    continue

                gene_lines[gene_id].append(line)

            except Exception as e:
                print(f"Warning: Error parsing line {line_num}: {e}")
                continue

    print(f"Found {len(gene_lines)} unique gene_ids")

    # Create mappings
    gene_to_file = {}
    file_to_genes = defaultdict(list)

    # Group genes into files
    file_index = 0
    current_gene_count = 0

    for gene_id in sorted(gene_lines.keys()):  # Sort for consistent output
        if current_gene_count >= genes_per_file:
            file_index += 1
            current_gene_count = 0

        filename = f"{output_prefix}-{file_index}.gtf"
        gene_to_file[gene_id] = filename
        file_to_genes[filename].append(gene_id)
        current_gene_count += 1

    print(f"Will create {file_index + 1} output files")

    # Write files
    for filename, genes_in_file in file_to_genes.items():
        with open(filename, 'w') as f:
            # Write header lines first
            f.writelines(header_lines)

            # Write all lines for genes in this file
            total_lines = 0
            for gene_id in genes_in_file:
                lines = gene_lines[gene_id]
                f.writelines(lines)
                total_lines += len(lines)

            print(f"Created {filename} with {len(genes_in_file)} genes and {total_lines} lines")

    # Print mapping summary
    print(f"\nGene to file mapping:")
    for gene_id, filename in sorted(gene_to_file.items()):
        print(f"  {gene_id} -> {filename}")

def main():
    if len(argv) < 2:
        print("Usage: python split_gtf_by_gene_id.py <input.gtf> [output_prefix] [genes_per_file]")
        print("Splits GTF file by gene_id, creating files with specified number of genes.")
        print("Handles unsorted gene_ids properly.")
        print("Default: 100 genes per file")
        return

    input_gtf = argv[1]
    output_prefix = argv[2] if len(argv) > 2 else None
    genes_per_file = int(argv[3]) if len(argv) > 3 else 100

    if not exists(input_gtf):
        print(f"Error: Input file {input_gtf} does not exist")
        return

    split_gtf_by_gene_id(input_gtf, output_prefix, genes_per_file)

if __name__ == '__main__':
    main()