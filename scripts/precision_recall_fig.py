import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

precision_recall_file = argv[1]
rec1, prec1 = int(argv[2]), float(argv[3])/100 # Rand1
rec2, prec2 = int(argv[4]), float(argv[5])/100 # RandX
outputName = "plot.pdf" if len(argv) <= 6 else argv[6]

# Read the data from the file
with open(precision_recall_file, 'r') as file:
    lines = file.readlines()
    pctIn = []
    recall = []
    precision = []
    for line in lines[10:]:  # Skip the header and a few points that are not informative
        values = line.strip().split('\t')
        pctIn.append(float(values[0]))
        recall.append(int(values[1]))
        precision.append(float(values[2]))

plt.scatter(rec1, prec1, color='red', marker='^', s=32, label='Rand1') 
plt.scatter(rec2, prec2, color='green', marker='s', s=24, label='RandX')

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
plt.xlabel('# supported isoforms', fontsize=16)
plt.ylabel(r'% supported', fontsize=16)
# plt.title('Precision-Recall Curve', fontsize=18)
plt.legend(fontsize=12, loc='upper right')  

plt.savefig(outputName, dpi=600, bbox_inches='tight')

# plt.show()



