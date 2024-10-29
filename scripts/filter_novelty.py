import sys
import GTF
import pandas as pd

def read_gtf(gtf):
    table = []
    with open(gtf, 'r') as f:
        for line in f.readlines():
            if line.startswith('#'): continue
            fields = GTF.parse(line)
            table.append(fields)
    return table


def extract_introns(gtf_table):
    """
    Extracts introns from a GTF table.
    
    Parameters:
    gtf_table (list): List of dictionaries containing GTF data.
    
    Returns:
    pd.DataFrame: DataFrame containing intron data.
    """
    introns = []
    transcripts = {}

    # Group exons by transcript_id
    for entry in gtf_table:
        if entry['feature'] == 'exon':
            transcript_id = entry['transcript_id']
            if transcript_id not in transcripts:
                transcripts[transcript_id] = []
            transcripts[transcript_id].append(entry)

    # Sort exons and extract introns
    for transcript_id, exons in transcripts.items():
        exons = sorted(exons, key=lambda x: int(x['start']))
        for i in range(1, len(exons)):
            intron = {
                'seqname': exons[i]['seqname'],
                'start': int(exons[i-1]['end']) + 1,
                'end': int(exons[i]['start']) - 1,
                'strand': exons[i]['strand'],
                'transcript_id': transcript_id                
            }
            introns.append(intron)
    return introns
    
# Return: transcripts_with_novel_introns, transcripts_without_novel_introns
def split_gtf_wrt_novel_introns(gtf, ref_gtf):
    introns = extract_introns(gtf)
    ref_introns = extract_introns(ref_gtf)

    ref_introns_set = set([(i['seqname'], i['start'], i['end'], i['strand']) for i in ref_introns])
    transcripts_id_with_novel_introns = set()
    for i in introns:
        if (i['seqname'], i['start'], i['end'], i['strand']) not in ref_introns_set:
            transcripts_id_with_novel_introns.add(i['transcript_id'])
    
    transcripts_with_novel_introns = []
    transcripts_without_novel_introns = []
    for line in gtf:
        if line['transcript_id'] in transcripts_id_with_novel_introns:
            transcripts_with_novel_introns.append(line)
        else:
            transcripts_without_novel_introns.append(line)
    return transcripts_with_novel_introns, transcripts_without_novel_introns 

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_novelty.py <reference_gtf> <input_gtf>")
        sys.exit(1)

    reference_gtf_path = sys.argv[1]
    input_gtf_path = sys.argv[2]

    reference_gtf_table = read_gtf(reference_gtf_path)
    input_gtf_table     = read_gtf(input_gtf_path)

    tx_with_novel_introns, tx_without_novel_introns = split_gtf_wrt_novel_introns(input_gtf_table, reference_gtf_table)
    
    GTF.write_file(tx_with_novel_introns, "tx_novel_introns.gtf")
    GTF.write_file(tx_without_novel_introns, "tx_existing_introns.gtf")
    

