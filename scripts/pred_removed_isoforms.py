#!/usr/bin/env python3
"""
Script to artificially remove isoforms from a reference GTF for benchmarking
transcript assembly accuracy.

Removes N isoforms from each multi-isoform transcript group if:
1. The removal does not remove the shortest transcript
2. The removal does not reduce total number of exons

The removed isoforms should theoretically be recoverable by assembly algorithms.
"""

import argparse
import random
import os
import sys
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove isoforms from GTF for benchmarking transcript recovery"
    )
    parser.add_argument(
        "--refgtf", "-r",
        required=True,
        help="Path to reference GTF file"
    )
    parser.add_argument(
        "--tennis", "-t",
        default="../programs/tennis.py",
        help="Path to tennis.py (default: ../programs/tennis.py)"
    )
    parser.add_argument(
        "--num_removal", "-n",
        type=int,
        default=1,
        help="Number of isoforms to remove per gene (default: 1)"
    )
    parser.add_argument(
        "--output_dir", "-o",
        default="removal_retrieval",
        help="Output directory (default: removal_retrieval)"
    )
    parser.add_argument(
        "--output_prefix", "-p",
        default="Ref_removal",
        help="Output file prefix (default: Ref_removal)"
    )
    parser.add_argument(
        "--seed", "-s",
        type=int,
        default=2025,
        help="Random seed for reproducibility"
    )
    return parser.parse_args()


def setup_tennis_import(tennis_path):
    """Resolve tennis.py path and add its directory to sys.path"""
    tennisSrc = tennis_path
    while os.path.islink(tennisSrc):
        newTarget = os.readlink(tennisSrc)
        if not os.path.isabs(newTarget):
            newTarget = os.path.join(os.path.dirname(os.path.abspath(tennisSrc)), newTarget)
        tennisSrc = newTarget
    tennisSrcDir = os.path.dirname(tennisSrc)
    sys.path.insert(0, tennisSrcDir)

    from tennis import Transcriptome
    return Transcriptome


def remove_isoforms(matrix, matrix_gids, num_removal):
    """
    Remove isoforms from transcript groups.

    Args:
        matrix: List of isoform groups (each group is a list of exon chains)
        matrix_gids: List of gene/group IDs corresponding to matrix
        num_removal: Number of isoforms to remove per group

    Returns:
        tuple: (removed_transcripts, remaining_transcripts, considered_gene_count)
    """
    considered_gene_count = 0
    removed_transcripts = []
    remaining_transcripts = []

    # Need at least num_removal + 2 isoforms to remove num_removal
    # (keep at least 2 remaining, and can't remove shortest)
    min_isoforms_required = num_removal + 2

    # Track duplicate count for reporting
    total_duplicates_removed = 0

    for idx in range(len(matrix)):
        group_id = matrix_gids[idx]

        # Deduplicate isoforms within the group
        seen_chains = set()
        group_isoforms = []
        for iso in matrix[idx]:
            chain_tuple = tuple(iso)
            if chain_tuple not in seen_chains:
                seen_chains.add(chain_tuple)
                group_isoforms.append(iso)
            else:
                total_duplicates_removed += 1

        if len(group_isoforms) < min_isoforms_required:
            continue

        # Count exons and get shortest transcript length
        shortest_isoform_length = min(len(iso) for iso in group_isoforms)

        # Track current exon counts (will be updated as we remove isoforms)
        exon_counter = Counter()
        for isoform in group_isoforms:
            if len(isoform) > 2:  # Must be multi-exon transcripts
                exon_counter.update(isoform)

        # Try to remove num_removal isoforms
        removed_from_group = []
        random.shuffle(group_isoforms)

        for isoform in list(group_isoforms):  # Iterate over a copy
            if len(removed_from_group) >= num_removal:
                break

            # Condition 1: Not the shortest transcript
            if len(isoform) == shortest_isoform_length:
                continue

            # Condition 2: Doesn't reduce total exon count
            is_exon_num_reduced = False
            for exon in isoform:
                if exon_counter[exon] <= 1:
                    is_exon_num_reduced = True
                    break
            if is_exon_num_reduced:
                continue

            # This isoform can be removed
            removed_from_group.append(isoform)
            group_isoforms.remove(isoform)

            # Update exon counter to reflect removal
            for exon in isoform:
                exon_counter[exon] -= 1

        # Only count this gene if we successfully removed the desired number
        if len(removed_from_group) == num_removal:
            for iso in removed_from_group:
                removed_transcripts.append((group_id, iso))
                print(f"Removed isoform {iso} exons from gene {group_id}")
            remaining_transcripts.append((group_id, group_isoforms))
            considered_gene_count += 1 

    if total_duplicates_removed > 0:
        print(f"Warning: Removed {total_duplicates_removed} duplicate isoforms from input")

    return removed_transcripts, remaining_transcripts, considered_gene_count


def save_gtf_files(tsm, removed_transcripts, remaining_transcripts, output_dir, output_prefix):
    """Save removed and remaining transcripts to GTF files"""
    # Assert no overlap between removed and remaining transcripts
    removed_chains = {tuple(chain) for gid, chain in removed_transcripts}
    remaining_chains = {tuple(chain) for gid, chains_list in remaining_transcripts for chain in chains_list}
    assert removed_chains.isdisjoint(remaining_chains), "Overlap found between removed and remaining transcripts!"

    os.makedirs(output_dir, exist_ok=True)
    save_basename = os.path.join(output_dir, output_prefix)

    # Save remaining transcripts
    remaining_file = save_basename + ".remaining.gtf"
    if os.path.exists(remaining_file):
        print(f"{remaining_file} exists. Removing it")
        os.remove(remaining_file)

    for gid, chains in remaining_transcripts:
        info = dict()
        tsm.save_novel_gene_in_gtf(chains, info, gid, "pexon_chain", addlN=9999, file=remaining_file)

    # Save removed transcripts
    removed_file = save_basename + ".artificial_removed.gtf"
    if os.path.exists(removed_file):
        print(f"{removed_file} exists. Removing it")
        os.remove(removed_file)

    for gid, chain in removed_transcripts:
        gchains = [chain]  # Put transcript in a gene
        info = dict()
        tsm.save_novel_gene_in_gtf(gchains, info, gid, "pexon_chain", addlN=1, file=removed_file)

    return remaining_file, removed_file


def main():
    args = parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    # Setup tennis import
    Transcriptome = setup_tennis_import(args.tennis)

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load transcriptome
    save_basename = os.path.join(args.output_dir, args.output_prefix)
    tsm = Transcriptome(
        args.refgtf,
        f'{save_basename}_tmp.stats',
        f'{save_basename}_tmp.pred.gtf'
    )

    # Get matrix representation
    matrix = tsm.get_chain_matrix('pexon_chain', 'tsstes_level')

    # Count multi-isoform genes
    multi_iso_gene = sum(1 for i in tsm.exonMatrixCollection if len(i) >= 2)
    print(f"Total multi-isoform genes: {multi_iso_gene}")

    # Remove isoforms
    removed_transcripts, remaining_transcripts, considered_gene_count = remove_isoforms(
        matrix, tsm.matrix_gids, args.num_removal
    )

    print(f"Removed {args.num_removal} isoform(s) from {considered_gene_count} transcript groups")
    print(f"Total removed transcripts: {len(removed_transcripts)}")

    # Show examples
    print("\nExamples:")
    shown = set()
    for gid, iso in removed_transcripts[:10]:
        if gid not in shown:
            print(f"  Gene {gid}: removed isoform with {len(iso)} exons")
            shown.add(gid)

    # Save GTF files
    remaining_file, removed_file = save_gtf_files(
        tsm, removed_transcripts, remaining_transcripts,
        args.output_dir, args.output_prefix
    )

    print(f"\nOutput files:")
    print(f"  Remaining transcripts: {remaining_file}")
    print(f"  Removed transcripts:   {removed_file}")
    print(f"\nUse '{remaining_file}' as input and '{removed_file}' as ground truth for benchmarking.")


if __name__ == "__main__":
    main()
