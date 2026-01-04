import pandas as pd
import matplotlib.pyplot as plt
import argparse
import re

def parse_gtf_attributes(attr_str):
    """Parse GTF attribute string to extract key-value pairs."""
    attrs = {}
    for item in attr_str.strip().split(';'):
        item = item.strip()
        if not item:
            continue
        match = re.match(r'(\S+)\s+"?([^"]+)"?', item)
        if match:
            attrs[match.group(1)] = match.group(2)
    return attrs

def load_gtf_class_codes(gtf_file):
    """Load GTF and extract transcript_id -> class_code mapping."""
    class_codes = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get('transcript_id')
            class_code = attrs.get('class_code')
            if transcript_id and class_code:
                class_codes[transcript_id] = class_code
    return class_codes

def select_row(group):
    """Select one row per query_id: prefer best_gtf, else highest inframe_length_pi."""
    best_gtf = group[group['notes'] == 'best_gtf']
    if len(best_gtf) > 0:
        return best_gtf.iloc[0]
    else:
        return group.loc[group['inframe_length_pi'].idxmax()]

def main():
    parser = argparse.ArgumentParser(description='Plot histogram of inframe_length_pi')
    parser.add_argument('--gtf', type=str, help='Optional GTF file to split by class_code')
    parser.add_argument('--gtf', type=str, help='Optional GTF file to split by class_code')
    args = parser.parse_args()

    # Read the TSV file
    df = pd.read_csv('orf.tsv', sep='\t')

    # For each query_id, select one entry
    selected = df.groupby('query_id').apply(select_row).reset_index(drop=True)

    if args.gtf:
        # Load class_code from GTF
        class_codes = load_gtf_class_codes(args.gtf)

        # Map class_code to selected transcripts
        selected['class_code'] = selected['query_id'].map(class_codes)

        # Split into matched (=) and not matched (!= '=')
        matched = selected[selected['class_code'] == '=']
        not_matched = selected[selected['class_code'] != '=']

        # Create two subplots
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        axes[0].hist(matched['inframe_length_pi'], bins=50, edgecolor='black', alpha=0.7, color='green')
        axes[0].set_xlabel('inframe_length_pi')
        axes[0].set_ylabel('Count')
        axes[0].set_title(f'Matched (class_code = "=") [n={len(matched)}]')

        axes[1].hist(not_matched['inframe_length_pi'], bins=50, edgecolor='black', alpha=0.7, color='orange')
        axes[1].set_xlabel('inframe_length_pi')
        axes[1].set_ylabel('Count')
        axes[1].set_title(f'Not Matched (class_code != "=") [n={len(not_matched)}]')

        plt.tight_layout()
        plt.savefig('inframe_length_pi_histogram_by_class.png', dpi=150)
        plt.show()

        print(f"Total unique query_ids: {len(selected)}")
        print(f"Matched: {len(matched)}, Not Matched: {len(not_matched)}")
        print(f"\nMatched summary:\n{matched['inframe_length_pi'].describe()}")
        print(f"\nNot Matched summary:\n{not_matched['inframe_length_pi'].describe()}")
    else:
        # Single histogram (original behavior)
        plt.figure(figsize=(10, 6))
        plt.hist(selected['inframe_length_pi'], bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('inframe_length_pi')
        plt.ylabel('Count')
        plt.title('Distribution of inframe_length_pi (one per query_id)')
        plt.tight_layout()
        plt.savefig('inframe_length_pi_histogram.png', dpi=150)
        plt.show()

        print(f"Total unique query_ids: {len(selected)}")
        print(f"Summary statistics:\n{selected['inframe_length_pi'].describe()}")

if __name__ == '__main__':
    main()
