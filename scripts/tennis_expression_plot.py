#!/usr/bin/env python3
"""
Plot validated TENNIS isoform expression vs total transcript group expression.
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


def compute_cpm(sample_data):
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

