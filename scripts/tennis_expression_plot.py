#!/usr/bin/env python3
"""
Plot validated TENNIS isoform expression vs average transcript group expression.
"""

import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import GTF


def parse_gtf(gtf_file):
    """Parse GTF and compute grouping key (first_exon_end, last_intron_start) per transcript."""
    transcript_exons = defaultdict(list)

    for fields in GTF.lines(gtf_file):
        if fields['feature'] != 'exon':
            continue

        chrom = fields['seqname']
        start = int(fields['start'])
        end = int(fields['end'])
        strand = fields['strand']
        transcript_id = fields['transcript_id']

        if transcript_id:
            transcript_exons[transcript_id].append((chrom, start, end, strand))

    # Compute grouping key per transcript
    transcript_to_group = {}
    for tid, exons in transcript_exons.items():
        if len(exons) < 2:
            # Single exon, use exon boundaries as group key
            continue
        else:
            # Sort exons by position
            exons_sorted = sorted(exons, key=lambda x: x[1])
            chrom = exons_sorted[0][0]
            strand = exons_sorted[0][3]
            first_exon_end = exons_sorted[0][2]
            last_intron_start = exons_sorted[-2][2]  # end of second-to-last exon
            transcript_to_group[tid] = (chrom, first_exon_end, last_intron_start, strand)
    # print(transcript_to_group)
    return transcript_to_group


def parse_tmap(tmap_file):
    """Parse tmap file and return validated TENNIS-to-annotation mapping."""
    tennis_to_annot = {}

    df = pd.read_csv(tmap_file, sep='\t')
    # Filter for exact matches
    matched = df[df['class_code'] == '=']

    for _, row in matched.iterrows():
        ref_id = row['ref_id']      # annotation ID (TCONS_*)
        qry_id = row['qry_id']      # TENNIS ID
        tennis_to_annot[qry_id] = ref_id

    return tennis_to_annot


def load_quant_files(quant_pattern):
    """Load quantification files matching pattern, return dict of sample -> DataFrame."""
    files = glob.glob(quant_pattern)
    if not files:
        raise ValueError(f"No files found matching pattern: {quant_pattern}")

    sample_data = {}
    for f in files:
        df = pd.read_csv(f, sep='\t')
        sample_data[f] = df

    return sample_data


def compute_cpm(sample_data) -> dict:
    """Compute CPM for each sample and return combined DataFrame with average CPM."""
    all_cpm = []

    for sample_name, df in sample_data.items():
        total_reads = df['num_reads'].sum()
        cpm = df[['tname']].copy()
        cpm['cpm'] = (df['num_reads'] / total_reads) * 1e6
        cpm = cpm.set_index('tname')['cpm']
        all_cpm.append(cpm)

    # Combine and average
    combined = pd.concat(all_cpm, axis=1)
    avg_cpm = combined.mean(axis=1)

    return avg_cpm.to_dict()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quant', required=True, help='Glob pattern for quantification files (e.g., "sample.*.quant")')
    parser.add_argument('--tmap', required=True, help='tmap file for TENNIS-to-annotation mapping')
    parser.add_argument('--gtf', required=True, help='GTF file with annotation transcripts')
    parser.add_argument('--output', default='tennis_vs_group_expression.pdf', help='Output plot filename')
    parser.add_argument('--min_cpm', type=float, default=0.01, help='Minimum CPM threshold for plotting (default: 0.01)')

    args = parser.parse_args()

    # Parse inputs
    print("Parsing GTF...")
    transcript_to_group = parse_gtf(args.gtf)

    print("Parsing tmap...")
    tennis_to_annot = parse_tmap(args.tmap)

    print(f"Loading quantification files matching: {args.quant}")
    sample_data = load_quant_files(args.quant)
    print(f"  Found {len(sample_data)} sample files")

    print("Computing CPM...")
    avg_cpm = compute_cpm(sample_data)

    # Build group -> list of transcript CPMs
    group_to_transcripts = defaultdict(list)

    # Add annotation transcripts to groups
    for tid, group_key in transcript_to_group.items():
        if tid in avg_cpm:
            group_to_transcripts[group_key].append((tid, avg_cpm[tid], 'annotation'))

    # Map TENNIS transcripts to groups via their matched annotation
    # Note: quant file uses annotation IDs, so we use annot_id to look up CPM
    tennis_in_groups = {}  # tennis_id -> (group_key, annot_id)
    for tennis_id, annot_id in tennis_to_annot.items():
        assert annot_id in transcript_to_group
        assert annot_id in avg_cpm
        if annot_id in transcript_to_group and annot_id in avg_cpm:
            group_key = transcript_to_group[annot_id]
            # Mark this annotation transcript as a validated TENNIS isoform
            tennis_in_groups[tennis_id] = (group_key, annot_id)

    # Compute plot data
    x_vals = []  # Average group CPM
    y_vals = []  # TENNIS isoform CPM

    for tennis_id, (group_key, annot_id) in tennis_in_groups.items():

        tennis_cpm = avg_cpm[annot_id]  # Use annot_id to look up CPM

        # Average CPM across all isoforms in the group
        group_cpms = [cpm for _, cpm, _ in group_to_transcripts[group_key]]
        avg_group_cpm = np.mean(group_cpms)

        if tennis_cpm >= args.min_cpm and avg_group_cpm >= args.min_cpm:
            x_vals.append(avg_group_cpm)
            y_vals.append(tennis_cpm)

    print(f"Plotting {len(x_vals)} validated TENNIS isoforms...")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.scatter(x_vals, y_vals, alpha=0.5, s=20, c='steelblue', edgecolors='none')

    # Add diagonal line (y = x, meaning TENNIS equals group average)
    min_val = min(min(x_vals), min(y_vals)) if x_vals else 0.01
    max_val = max(max(x_vals), max(y_vals)) if x_vals else 100
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='y = x (equal to avg)')

    # Add lines for different ratios
    for ratio, color in [(2, 'gray'), (0.5, 'lightgray')]:
        ax.plot([min_val, max_val], [min_val * ratio, max_val * ratio],
                '--', color=color, alpha=0.5, label=f'y = {ratio}x avg')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Average transcript group CPM (log scale)')
    ax.set_ylabel('TENNIS isoform CPM (log scale)')
    ax.set_title('Validated TENNIS isoform expression vs average group expression')
    ax.legend(loc='lower right')

    # Equal aspect ratio for log-log
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {args.output}")

    # Print summary statistics
    if x_vals:
        ratios = np.array(y_vals) / np.array(x_vals)
        print(f"\nSummary statistics:")
        print(f"  Number of validated TENNIS isoforms plotted: {len(x_vals)}")
        print(f"  Median TENNIS/total ratio: {np.median(ratios):.3f}")
        print(f"  Mean TENNIS/total ratio: {np.mean(ratios):.3f}")
        print(f"  TENNIS isoforms >50% of group: {np.sum(ratios > 0.5)} ({100*np.sum(ratios > 0.5)/len(ratios):.1f}%)")
        print(f"  TENNIS isoforms >10% of group: {np.sum(ratios > 0.1)} ({100*np.sum(ratios > 0.1)/len(ratios):.1f}%)")


if __name__ == '__main__':
    main()
