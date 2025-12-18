#!/usr/bin/env python3
"""
Plot validated TENNIS isoform expression vs transcript group expression.
Generates two plots: one comparing to group total, one comparing to group average.
"""

import argparse
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import GTF

# Set global font size
plt.rcParams.update({'font.size': 16})


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
    parser.add_argument('--min_cpm_isoform', type=float, default=0.01, help='Minimum CPM threshold for scatter plot points (default: 0.01)')
    parser.add_argument('--min_cpm_group', type=float, default=0.01, help='Minimum CPM threshold for histogram groups (default: 0.01)')

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
    x_avg_vals = []  # Average group CPM (for validated TENNIS isoforms)
    x_total_vals = []  # Total group CPM (for validated TENNIS isoforms)
    y_vals = []  # TENNIS isoform CPM
    all_group_avg_cpms = []  # Average CPM for ALL groups (for histogram)
    all_group_total_cpms = []  # Total CPM for ALL groups (for histogram)

    # Compute average and total CPM for all groups (for histogram)
    for group_key, transcripts in group_to_transcripts.items():
        group_cpms = [cpm for _, cpm, _ in transcripts]
        avg_group_cpm = np.mean(group_cpms)
        total_group_cpm = np.sum(group_cpms)
        if avg_group_cpm >= args.min_cpm_group:
            all_group_avg_cpms.append(avg_group_cpm)
        if total_group_cpm >= args.min_cpm_group:
            all_group_total_cpms.append(total_group_cpm)

    # Compute scatter plot data (for validated TENNIS isoforms)
    for tennis_id, (group_key, annot_id) in tennis_in_groups.items():

        tennis_cpm = avg_cpm[annot_id]  # Use annot_id to look up CPM

        # Average and total CPM across all isoforms in the group
        group_cpms = [cpm for _, cpm, _ in group_to_transcripts[group_key]]
        avg_group_cpm = np.mean(group_cpms)
        total_group_cpm = np.sum(group_cpms)

        if tennis_cpm >= args.min_cpm_isoform and avg_group_cpm >= args.min_cpm_isoform:
            x_avg_vals.append(avg_group_cpm)
            x_total_vals.append(total_group_cpm)
            y_vals.append(tennis_cpm)

    print(f"Plotting {len(y_vals)} validated TENNIS isoforms...")

    # Determine output filenames
    base_output = args.output.rsplit('.', 1)[0] if '.' in args.output else args.output
    ext = args.output.rsplit('.', 1)[1] if '.' in args.output else 'pdf'

    # ===== Plot 1: TENNIS vs Group Total =====
    fig1, (ax_main1, ax_hist1) = plt.subplots(
        2, 1, figsize=(8, 10),
        gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05},
        sharex=True
    )

    ax_main1.scatter(x_total_vals, y_vals, alpha=0.5, s=20, c='steelblue', edgecolors='none')

    min_val1 = min(min(x_total_vals), min(y_vals)) if y_vals else 0.01
    max_val1 = max(max(x_total_vals), max(y_vals)) if y_vals else 100

    for pct, color in [(0.95, 'black'), (0.05, 'lightgray')]:
        ax_main1.plot([min_val1, max_val1], [min_val1 * pct, max_val1 * pct],
                '--', color=color, alpha=0.5, label=f'y = {int(pct*100)}% of x')

    ax_main1.set_xscale('log')
    ax_main1.set_yscale('log')
    ax_main1.set_ylabel('TENNIS isoform CPM (log scale)')
    # ax_main1.set_title('Validated TENNIS isoform expression vs total group expression')
    ax_main1.legend(loc='lower right')

    log_bins1 = np.logspace(np.log10(args.min_cpm_group), np.log10(max(all_group_total_cpms)), 50)
    ax_hist1.hist(all_group_total_cpms, bins=log_bins1, color='gray', alpha=0.7, edgecolor='none')
    ax_hist1.set_xscale('log')
    ax_hist1.set_xlabel('Total transcript group CPM (log scale)')
    ax_hist1.set_ylabel('Count of groups')

    plt.tight_layout()
    output_total = f"{base_output}_total.{ext}"
    plt.savefig(output_total, dpi=300, bbox_inches='tight')
    plt.close(fig1)
    print(f"Saved plot to {output_total}")

    # ===== Plot 2: TENNIS vs Group Average =====
    fig2, (ax_main2, ax_hist2) = plt.subplots(
        2, 1, figsize=(8, 10),
        gridspec_kw={'height_ratios': [4, 1], 'hspace': 0.05},
        sharex=True
    )

    ax_main2.scatter(x_avg_vals, y_vals, alpha=0.5, s=20, c='steelblue', edgecolors='none')

    min_val2 = min(min(x_avg_vals), min(y_vals)) if y_vals else 0.01
    max_val2 = max(max(x_avg_vals), max(y_vals)) if y_vals else 100
    ax_main2.plot([min_val2, max_val2], [min_val2, max_val2], 'k--', alpha=0.5, label='y = x (equal to avg)')

    for ratio, color in [(10, 'gray'), (0.1, 'gray')]:
        ax_main2.plot([min_val2, max_val2], [min_val2 * ratio, max_val2 * ratio],
                '--', color=color, alpha=0.5, label=f'y = {ratio}x avg')

    ax_main2.set_xscale('log')
    ax_main2.set_yscale('log')
    ax_main2.set_ylabel('TENNIS isoform CPM (log scale)')
    # ax_main2.set_title('Validated TENNIS isoform expression vs average group expression')
    ax_main2.legend(loc='lower right')

    log_bins2 = np.logspace(np.log10(args.min_cpm_group), np.log10(max(all_group_avg_cpms)), 50)
    ax_hist2.hist(all_group_avg_cpms, bins=log_bins2, color='gray', alpha=0.7, edgecolor='none')
    ax_hist2.set_xscale('log')
    ax_hist2.set_xlabel('Average transcript group CPM (log scale)')
    ax_hist2.set_ylabel('Count of groups')

    plt.tight_layout()
    output_avg = f"{base_output}_average.{ext}"
    plt.savefig(output_avg, dpi=300, bbox_inches='tight')
    plt.close(fig2)
    print(f"Saved plot to {output_avg}")

    # Print summary statistics
    if y_vals:
        # Total stats
        ratios_total = np.array(y_vals) / np.array(x_total_vals)
        print(f"\nSummary statistics (vs Total):")
        print(f"  Number of validated TENNIS isoforms plotted: {len(y_vals)}")
        print(f"  Median TENNIS/total ratio: {np.median(ratios_total):.3f}")
        print(f"  Mean TENNIS/total ratio: {np.mean(ratios_total):.3f}")
        print(f"  TENNIS isoforms >95% of total: {np.sum(ratios_total > 0.95)} ({100*np.sum(ratios_total > 0.95)/len(ratios_total):.1f}%)")
        print(f"  TENNIS isoforms >5% of total: {np.sum(ratios_total > 0.05)} ({100*np.sum(ratios_total > 0.05)/len(ratios_total):.1f}%)")

        # Average stats
        ratios_avg = np.array(y_vals) / np.array(x_avg_vals)
        print(f"\nSummary statistics (vs Average):")
        print(f"  Median TENNIS/avg ratio: {np.median(ratios_avg):.3f}")
        print(f"  Mean TENNIS/avg ratio: {np.mean(ratios_avg):.3f}")
        print(f"  TENNIS isoforms > 10x avg: {np.sum(ratios_avg > 10)} ({100*np.sum(ratios_avg > 10)/len(ratios_avg):.1f}%)")
        print(f"  TENNIS isoforms > avg: {np.sum(ratios_avg > 1)} ({100*np.sum(ratios_avg > 1)/len(ratios_avg):.1f}%)")
        print(f"  TENNIS isoforms > 0.1x avg: {np.sum(ratios_avg > 0.1)} ({100*np.sum(ratios_avg > 0.1)/len(ratios_avg):.1f}%)")


if __name__ == '__main__':
    main()
