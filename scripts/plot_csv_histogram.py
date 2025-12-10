#!/usr/bin/env python3
"""
Plot histogram of a specified column from stats.csv, colored by feasibility.
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description='Plot histogram of a column, colored by real_additional_nodes feasibility'
    )
    parser.add_argument(
        'csv_file',
        help='Path to the CSV file'
    )
    parser.add_argument(
        '--column',
        '-c',
        required=True,
        help='Column name to plot histogram for'
    )
    parser.add_argument(
        '--output',
        '-o',
        help='Output file path (optional, will show plot if not specified)'
    )
    parser.add_argument(
        '--bins',
        '-b',
        type=int,
        default=-1,
        help='Number of bins for histogram (default: 30)'
    )
    parser.add_argument(
        '--nostack',
        action='store_true',
        help='Do not stack feasible/infeasible bars (show side-by-side instead)'
    )

    args = parser.parse_args()

    # Read CSV file
    try:
        df = pd.read_csv(args.csv_file)
    except FileNotFoundError:
        print(f"Error: File '{args.csv_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # Check if required columns exist
    if args.column not in df.columns:
        print(f"Error: Column '{args.column}' not found in CSV", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)

    if 'real_additional_nodes' not in df.columns:
        print("Error: Column 'real_additional_nodes' not found in CSV", file=sys.stderr)
        sys.exit(1)

    # Create label for legend: 'infeasible' for negative values, actual number for non-negative
    df['node_label'] = df['real_additional_nodes'].apply(
        lambda x: 'infeasible' if x < 0 else str(int(x))
    )

    # Create histogram
    fig, ax = plt.subplots(figsize=(10, 6))

    # Get unique values sorted (infeasible last)
    unique_values = sorted(df['real_additional_nodes'].unique())
    non_negative = [v for v in unique_values if v >= 0]
    negative = [v for v in unique_values if v < 0]

    # Prepare data and labels for each category
    data_list = []
    labels = []
    colors = []

    # Use a colormap for non-negative values
    cmap = plt.cm.Greens
    num_non_neg = len(non_negative)

    for i, val in enumerate(non_negative):
        data_list.append(df[df['real_additional_nodes'] == val][args.column])
        labels.append(str(int(val)))
        # Scale colors from light to dark green
        colors.append(cmap(0.3 + 0.6 * i / max(num_non_neg - 1, 1)))

    # Add infeasible (all negative values grouped together)
    if negative:
        data_list.append(df[df['real_additional_nodes'] < 0][args.column])
        labels.append('infeasible')
        colors.append('#e74c3c')

    # Plot with different colors
    if args.bins > 0:
        bins = args.bins
    else:
        bins = list(range(0, 11)) + [np.inf]
    ax.hist(
        data_list,
        bins=bins,
        label=labels,
        color=colors,
        alpha=0.7,
        edgecolor='black',
        linewidth=0.5,
        stacked=not args.nostack
    )

    # ax.set_xlabel(args.column, fontsize=12)
    if args.column != "computed_upper_bound":
        ax.set_xlabel(args.column, fontsize=12)
    else:
        ax.set_xlabel('MST-Computed Upper Bound', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    # ax.set_title(f'Histogram of {args.column} by Feasibility', fontsize=14)
    ax.legend(title='TENNIS-SAT additional nodes', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save or show
    if args.output:
        plt.savefig(args.output, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    else:
        plt.show()


if __name__ == '__main__':
    main()
