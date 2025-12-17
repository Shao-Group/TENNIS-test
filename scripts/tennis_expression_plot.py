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
