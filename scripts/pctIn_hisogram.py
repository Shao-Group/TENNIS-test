import matplotlib.pyplot as plt
import numpy as np
from sys import argv

precision_recall_file = argv[1]
outputName = "hist.pdf" if len(argv) <= 2 else argv[2]

with open(precision_recall_file, 'r') as file:
    lines = file.readlines()
    pctIn = []
    pctIn_correct = []
    recall = []
    precision = []
    for line in lines[1:]:  # Skip the header
        values = line.strip().split('\t')
        pctIn.append(float(values[0]))
        if len(recall) == 0:
            if int(values[1]) > 0:
                print(int(values[1]), line)
                assert int(values[1]) == 1
                pctIn_correct.append((float(values[0])))
        elif int(values[1]) - recall[-1] > 0:           # increase recall, so current is correct
            pctIn_correct.append((float(values[0])))
        recall.append(int(values[1]))
        precision.append(float(values[2]))


plt.figure(figsize=(8, 6))
plt.hist(pctIn, bins=50, range=(0,1), color='skyblue', label='total')
plt.hist(pctIn_correct, bins=50, range=(0,1), color='peachpuff', label='supported')
plt.xlabel('PctIn', fontsize=16)
plt.ylabel('# predicted transcripts', fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(False)
plt.legend(fontsize=14, loc='upper right')
plt.savefig(outputName, dpi=600, bbox_inches='tight')
plt.show()


exit()