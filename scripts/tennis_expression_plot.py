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

