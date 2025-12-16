import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

precision_recall_file = argv[1] # from precision_recall_by_pctIn.py
rec1, prec1 = float(argv[2]), float(argv[3])/100 # Rand1 mean
stdrec1, stdprec1 = float(argv[4]), float(argv[5])/100 # Rand1 std
rec2, prec2 = float(argv[6]), float(argv[7])/100 # RandX mean
stdrec2, stdprec2 = float(argv[8]), float(argv[9])/100 # RandX std
rec3, prec3 = float(argv[10]), float(argv[11])/100 # PSI1
rec4, prec4 = float(argv[12]), float(argv[13])/100 # PSIX
outputName = "plot.pdf" if len(argv) <= 14 else argv[14]

# Read the data from the file
with open(precision_recall_file, 'r') as file:
    lines = file.readlines()
    pctIn = []
    recall = []
    precision = []
    for line in lines[10:]:  # Skip the header and a few points that are not informative
        values = line.strip().split('\t')
        pctIn.append(float(values[0].strip('()').split(',')[0]) if values[0].startswith('(') else float(values[0]))
        recall.append(int(values[1]))
        precision.append(float(values[2]))

plt.errorbar(rec1, prec1, xerr=stdrec1, yerr=stdprec1, fmt='^', color='red', markersize=6, capsize=3)
plt.errorbar(rec2, prec2, xerr=stdrec2, yerr=stdprec2, fmt='s', color='green', markersize=5, capsize=3)
plt.scatter([], [], color='red', marker='^', s=36, label='Rand1')
plt.scatter([], [], color='green', marker='s', s=30, label='RandX')
plt.scatter(rec3, prec3, color='blue', marker='d', s=32, label='PSI1')
plt.scatter(rec4, prec4, color='orange', marker='o', s=32, label='PSIX')

colors_list = ['#FF00FF',       
               '#FFE700',      
               '#00A693',       
               '#004165',      
               '#002040']       

custom_cmap = colors.LinearSegmentedColormap.from_list("custom", colors_list)
norm = colors.Normalize(vmin=min(pctIn), vmax=max(pctIn))

for i, (x, y, val) in enumerate(zip(recall, precision, pctIn)):
    if val < 0.5:  # circle
        plt.plot(x, y, 'o', color='black', fillstyle='none', markersize=10)
        break

for i, (x, y, val) in enumerate(zip(recall, precision, pctIn)):
    if val < 0.3333:  # circle
        plt.plot(x, y, 'o', color='black', fillstyle='none', markersize=10)
        break

plt.scatter(recall, precision, c=pctIn, cmap=custom_cmap, norm=norm, marker='.', s=1)
scatter = plt.scatter(recall, precision, c=pctIn, 
                     cmap=custom_cmap,
                     norm=norm,
                     marker='.', 
                     s=1)

cbar = plt.colorbar(scatter, label='PctIn', aspect=50)
cbar.ax.tick_params(labelsize=14) 
cbar.set_label('PctIn', size=16)  


plt.grid(False)
if  'rmrtr' in outputName:
    plt.ylim(0, 1)
else:    
    plt.ylim(0, 0.6)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('# mathching transcripts', fontsize=16)
plt.ylabel(r'% matched', fontsize=16)
# plt.title('Precision-Recall Curve', fontsize=18)
plt.legend(fontsize=12, loc='upper right')  

plt.savefig(outputName, dpi=600, bbox_inches='tight')

# plt.show()



