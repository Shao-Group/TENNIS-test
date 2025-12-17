import sys
import sys
import os
if len(sys.argv) >= 5:
    p = sys.argv[4]
    sys.path.insert(0, p)
from GTF import get_attr


# get class code from .tmap
# return {tid: code}
def get_class_code(tmap):
    tid2class_code = {}
    f = open(tmap, 'r')
    line_count = 0
    for line in f.readlines():
        if line_count == 0:
            # print("skip header")
            line_count += 1
            continue
        fields = line.strip().split('\t')
        assert len(fields) >= 5
        tid2class_code[fields[4]] = fields[2]
        line_count += 1
    f.close()
    return tid2class_code

def precision_recall_curv(gtf, tmap, attr_name, output_file=None, filter_k=None, filter_v=None):
    if attr_name == 'pctIn_PSI_score':
        tid2attr = get_attr(gtf, 'PctIn', 'PSI_score', filter_k=filter_k, filter_v=filter_v)
    else:
        tid2attr = get_attr(gtf, attr_name, filter_k=filter_k, filter_v=filter_v)
    tid2class_code = get_class_code(tmap)

    sorted_tid2attr = sorted(tid2attr.items(), key=lambda x: x[1], reverse=True)
    cur_total = 0
    cur_true  = 0
    cur_pct   = None
    info = []
    for k,v in sorted_tid2attr:
        cur_total += 1
        if cur_pct is None:
            cur_pct = v
        assert v <= cur_pct
        cur_pct = v
        if tid2class_code[k] == '=':
            cur_true += 1
        cur_info = (cur_pct, cur_true, cur_true/cur_total, cur_total)
        # print(cur_info)
        info.append(cur_info)

    if output_file is None:
        output_file = 'precision_recall_info.txt'
    with open(output_file, 'w') as f:
        f.write("\t".join([attr_name, '#True', '%True', 'Total']) + "\n")
        for line in info:
            l = '\t'.join([str(x) for x in line])
            f.write(f"{l}\n")

if __name__ == "__main__":
    gtf = sys.argv[1]
    tmap = sys.argv[2]
    attr = sys.argv[3]
    if len(sys.argv) >= 6:
        output_file = sys.argv[4]
        filter_k, filter_v = sys.argv[5].split('=')
        precision_recall_curv(gtf, tmap, attr, output_file, filter_k, filter_v)
    elif len(sys.argv) >= 5:
        output_file = sys.argv[4]
        precision_recall_curv(gtf, tmap, attr, output_file)
    else:
        precision_recall_curv(gtf, tmap, attr)



