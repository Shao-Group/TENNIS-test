from sys import argv
from GTF import parse as parse_gtf_line
from os.path import basename, exists
from collections import defaultdict


def quick_scan(gtf):
    # gid_span = {gid:(first appearance, last appearance)}
    # we split file at lines that is not inside any gid_span
    gid_span = dict() 

    f = open(gtf, 'r')
    line_counter = 0
    splits = []
    for line in f.readlines():
        if line.startswith('#'): continue
        line_counter += 1
        fields = parse_gtf_line(line)
        gid = fields['gene_id']
        if gid not in gid_span:
            gid_span[gid] = (line_counter, line_counter)
        else:
            s, t = gid_span[gid]
            gid_span[gid] = (min(s, line_counter), max(t, line_counter))
    
    sorted_span = sorted(list(gid_span.values()))
    sorted_start = sorted([x[0] for x in list(gid_span.values())])
    sorted_end   = sorted([x[1] for x in list(gid_span.values())], reverse=True)
    
    GENE_NUM_PER_FILE = 500
    for x in range(1, len(sorted_span), GENE_NUM_PER_FILE):
        s,t = sorted_span[x]
        assert s == sorted_start[x]
        if t == sorted_end[x]:
            splits.append(t)
        elif t < sorted_end[x]:
            print('previous gene spans into this one')
        else:
            print('this gene spans into later ones')
    f.close()
    return splits


# assuming gtf is sorted by gene ids
def split(gtf):
    global file_index
    line_batch = []
    gene_ids_permanent = set()
    gene_ids_batch = set()
    GENE_NUM_PER_FILE = 500

    splits = quick_scan(gtf)

    f = open(gtf, 'r')
    line_counter = 0
    split_counter = 0
    for line in f.readlines():
        if line.startswith('#'): continue
        line_counter += 1
        if line_counter > splits[split_counter]:
            # flush out previous genes
            assert (line_counter - 1 ==splits[split_counter])
            split_counter += 1
            write_file(line_batch)
            line_batch = []
            gene_ids_batch = set()
            file_index += 1
        fields = parse_gtf_line(line)
        gid = fields['gene_id']
        # new gid
        if gid not in gene_ids_batch:  
            # record new genes
            if gid in gene_ids_permanent:
                print(gid)
            assert  gid not in gene_ids_permanent# otherwise gtf is not sorted, a gene might be in different files
            gene_ids_permanent.add(gid)
            gene_ids_batch.add(gid)
        line_batch.append(line)
    write_file(line_batch)
    f.close()

def write_file(lines):
    global file_index
    global output_prefix
    if (file_index % 10 == 0 ): 
        print(file_index)
    with open(output_prefix + "-" + str(file_index) + ".gtf",'w') as f:
        f.writelines(lines)


if __name__ == '__main__':
    gtf = argv[1]
    print(gtf)
    # some dangerous GLOBAL variables
    file_index = 0
    output_prefix = basename(gtf) 
    output_prefix = (argv[2]) if len(argv) >= 3 else basename(gtf)
    print(output_prefix)
    split(gtf)
