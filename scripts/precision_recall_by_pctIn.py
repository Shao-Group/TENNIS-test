import sys
import sys
import os
if len(sys.argv) >= 5:
    p = sys.argv[4]
    sys.path.insert(0, p)
from GTF import get_PctIn


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

def precision_recall_curv(gtf, tmap, output_file=None):
    tid2PctIn = get_PctIn(gtf)
    tid2class_code = get_class_code(tmap)
    
    sorted_tid2PctIn = sorted(tid2PctIn.items(), key=lambda x: x[1], reverse=True)
    cur_total = 0
    cur_true  = 0
    cur_pct   = 1
    info = []
    for k,v in sorted_tid2PctIn:
        cur_total += 1
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
        f.write("\t".join(['PctIn', '#True', '%True', 'Total']) + "\n")
        for line in info:
            l = '\t'.join([str(x) for x in line])
            f.write(f"{l}\n")

if __name__ == "__main__":
    gtf = sys.argv[1]
    tmap = sys.argv[2]
    if len(sys.argv) >= 4:
        output_file = sys.argv[3]
        precision_recall_curv(gtf, tmap, output_file)
    else:
        precision_recall_curv(gtf, tmap)



