#!/usr/bin/env python3
"""
Filter GTF file based on matching first intron start and last intron end positions.

Given two GTF files:
1. Reference GTF (first file) - extract intron boundaries for each transcript
2. Query GTF (second file) - filter transcripts that match the reference boundaries

For each transcript, we:
- Calculate first intron start position (end of first exon)
- Calculate last intron end position (start of last exon)
- Keep transcripts in query GTF that match any transcript in reference GTF
"""

import sys
from collections import defaultdict


def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string into a dictionary."""
    attributes = {}
    for item in attr_string.strip().split(';'):
        item = item.strip()
        if item:
            parts = item.split(' ', 1)
            if len(parts) == 2:
                key = parts[0]
                value = parts[1].strip('"')
                attributes[key] = value
    return attributes


def get_intron_boundaries(gtf_file):
    """
    Extract first intron start and last intron end for each transcript.

    Returns:
        dict: {transcript_id: (chrom, strand, first_intron_start, last_intron_end)}
    """
    # First, collect all exons for each transcript
    transcripts = defaultdict(list)

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'exon':
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_gtf_attributes(fields[8])

            transcript_id = attributes.get('transcript_id')
            if not transcript_id:
                continue

            transcripts[transcript_id].append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'strand': strand
            })

    # Calculate intron boundaries for each transcript
    intron_boundaries = {}

    for transcript_id, exons in transcripts.items():
        if len(exons) < 2:
            # Need at least 2 exons to have introns
            continue

        chrom = exons[0]['chrom']
        strand = exons[0]['strand']

        # Sort exons by start position
        exons_sorted = sorted(exons, key=lambda x: x['start'])

        # First intron start = end of first exon
        first_intron_start = exons_sorted[0]['end']

        # Last intron end = start of last exon
        last_intron_end = exons_sorted[-1]['start']

        intron_boundaries[transcript_id] = (chrom, strand, first_intron_start, last_intron_end)

    return intron_boundaries


def create_boundary_index(intron_boundaries):
    """
    Create an index for quick lookup of intron boundaries.

    Returns:
        dict: {(chrom, strand, first_intron_start, last_intron_end): [transcript_ids]}
    """
    boundary_index = defaultdict(list)

    for transcript_id, (chrom, strand, first_start, last_end) in intron_boundaries.items():
        key = (chrom, strand, first_start, last_end)
        boundary_index[key].append(transcript_id)

    return boundary_index


def filter_gtf_by_boundaries(query_gtf, boundary_index):
    """
    Filter query GTF file to keep only transcripts with matching intron boundaries.

    Args:
        query_gtf: Path to query GTF file
        boundary_index: Index of valid intron boundaries

    Returns:
        list: Filtered GTF lines
    """
    # First, get intron boundaries for query transcripts
    query_boundaries = get_intron_boundaries(query_gtf)

    # Determine which query transcripts to keep
    transcripts_to_keep = set()

    for transcript_id, (chrom, strand, first_start, last_end) in query_boundaries.items():
        key = (chrom, strand, first_start, last_end)
        if key in boundary_index:
            transcripts_to_keep.add(transcript_id)

    # Filter the GTF file
    filtered_lines = []

    with open(query_gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                filtered_lines.append(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            attributes = parse_gtf_attributes(fields[8])
            transcript_id = attributes.get('transcript_id')

            if transcript_id in transcripts_to_keep:
                filtered_lines.append(line)

    return filtered_lines


def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_gtf_by_intron_boundaries.py <reference.gtf> <query.gtf> <output.gtf>")
        print("\nDescription:")
        print("  Filters query GTF to keep transcripts with matching intron boundaries from reference GTF.")
        print("  Matching is based on: chromosome, strand, first intron start, and last intron end.")
        print("\nArguments:")
        print("  reference.gtf  - Reference GTF file to extract intron boundaries from")
        print("  query.gtf      - Query GTF file to filter")
        print("  output.gtf     - Output filtered GTF file")
        sys.exit(1)

    reference_gtf = sys.argv[1]
    query_gtf = sys.argv[2]
    output_gtf = sys.argv[3]

    print(f"Reading reference GTF: {reference_gtf}")
    ref_boundaries = get_intron_boundaries(reference_gtf)
    print(f"Found {len(ref_boundaries)} transcripts with introns in reference")

    print(f"\nCreating boundary index...")
    boundary_index = create_boundary_index(ref_boundaries)
    print(f"Found {len(boundary_index)} unique intron boundary combinations")

    print(f"\nFiltering query GTF: {query_gtf}")
    filtered_lines = filter_gtf_by_boundaries(query_gtf, boundary_index)

    print(f"\nWriting filtered GTF to: {output_gtf}")
    with open(output_gtf, 'w') as f:
        for line in filtered_lines:
            f.write(line)

    # Count filtered transcripts
    filtered_transcripts = set()
    for line in filtered_lines:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                attributes = parse_gtf_attributes(fields[8])
                transcript_id = attributes.get('transcript_id')
                if transcript_id:
                    filtered_transcripts.add(transcript_id)

    print(f"\nResults:")
    print(f"  Kept {len(filtered_transcripts)} transcripts in output")
    print(f"  Wrote {len(filtered_lines)} lines to output file")


if __name__ == '__main__':
    main()
