import re
import os
from sys import argv


def extract_numbers_from_file(file_path):
    # Open the file and read it from the end
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Look for the specific patterns in the last few lines of the file
    for line in reversed(lines[-8:]):  # Only check the last 6 lines, assuming info is at the end
        if "gene count" in line:
            gene_count = int(re.search(r'\d+', line).group())
        elif "single isoform gene" in line:
            single_isoform_gene = int(re.search(r'\d+', line).group())
        elif "big gene with 100+ isoforms" in line:
            big_gene_count = int(re.search(r'isoforms \d+', line).group()[9:])
        elif "infeasiible count" in line:
            infeasible_count = int(re.search(r'\d+', line).group())
        elif "feasible counts" in line:
            feasible_counts = [int(x) for x in re.findall(r'\d+', line)]

    # Return the extracted values
    return gene_count, single_isoform_gene, big_gene_count, infeasible_count, feasible_counts

if __name__ == "__main__":
    directory = argv[1]
    stats_files = [f for f in os.listdir(directory) if f.endswith('.stats')]

    total_gene_count = 0
    total_single_isoform_gene = 0
    total_big_gene_count = 0
    total_infeasible_count = 0
    total_feasible_counts = []

    
    for file in stats_files:
        file_path = os.path.join(directory, file)
        gene_count, single_isoform_gene, big_gene_count, infeasible_count, feasible_counts = extract_numbers_from_file(file_path)

        total_gene_count += gene_count
        total_single_isoform_gene += single_isoform_gene
        total_big_gene_count += big_gene_count
        total_infeasible_count += infeasible_count
        if not total_feasible_counts:
            total_feasible_counts = feasible_counts  # Initialize if empty
        else:
            total_feasible_counts = [x + y for x, y in zip(total_feasible_counts, feasible_counts)]

        print(f"Gene Count: {gene_count}")
        print(f"Single Isoform Gene: {single_isoform_gene}")
        print(f"Big Gene Count: {big_gene_count}")
        print(f"Infeasible Count: {infeasible_count}")
        print(f"Feasible Counts: {feasible_counts}")
    
    print(f"Total Gene Count: {total_gene_count}")
    print(f"Total Single Isoform Gene: {total_single_isoform_gene}")
    print(f"Total Big Gene Count: {total_big_gene_count}")
    print(f"Total Infeasible Count: {total_infeasible_count}")
    print(f"Total Feasible Counts: {total_feasible_counts}")
    print(f"Total examined gene count: {total_gene_count - total_single_isoform_gene - total_big_gene_count}")