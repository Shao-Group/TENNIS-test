import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from sys import argv

# produce dm6 long read plot
# python ../analyses/precision_recall_fig.py precision_recall_dm6_longread.txt 149 23.2 171 18.3 dm6.lr.pdf

# removel retrival
# 1 exact
# python ../analyses/precision_recall_fig.py exact.retrival.ROC.txt 97 19.5 107 14.7 exact.rmrtr.pdf
# 2 both
# python ../analyses/precision_recall_fig.py both.retrival.ROC.txt 168 33.8 205 28.2 both.rmrtr.pdf

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
    for line in lines[10:]:  # Skip the header
        values = line.strip().split('\t')
        pctIn.append(float(values[0]))
        recall.append(int(values[1]))
        precision.append(float(values[2]))

# Calculate the average precision and recall
# NOT PLOTTED
avg_recall = 0
avg_total = sum(pctIn)
for i in range(len(recall)):
    is_true = (recall[i] >= recall[i-1] + 1) if i >= 1 else (recall[i] > 0)
    if is_true:
        avg_recall += pctIn[i]
avg_precision = avg_recall/avg_total
# print(avg_recall, avg_precision)

plt.scatter(rec1, prec1, color='red', marker='^', s=32, label='Rand1') 
plt.scatter(rec2, prec2, color='green', marker='s', s=24, label='RandX')
# plt.scatter(avg_recall, avg_precision, color='blue', marker='*', s=100, label='Avg') 


print(f"pctIn range: {min(pctIn)} to {max(pctIn)}")

# OR Approach 2: Normalize
colors_list = ['#FF00FF',      # magenta
               '#FFE700',      # yellow
               '#00A693',      # teal
               '#004165',      # blue
               '#002040']      # dark blue

custom_cmap = colors.LinearSegmentedColormap.from_list("custom", colors_list)
norm = colors.Normalize(vmin=min(pctIn), vmax=max(pctIn))

for i, (x, y, val) in enumerate(zip(recall, precision, pctIn)):
    if val < 0.5:  # Adjust this range as needed
        circle_color = custom_cmap(norm(val))  # This gets the exact same color as the point
        plt.plot(x, y, 'o', color='black', fillstyle='none', markersize=10)
        break

for i, (x, y, val) in enumerate(zip(recall, precision, pctIn)):
    if val < 0.3333:  # Adjust this range as needed
        circle_color = custom_cmap(norm(val))  # This gets the exact same color as the point
        plt.plot(x, y, 'o', color='black', fillstyle='none', markersize=10)
        break

plt.scatter(recall, precision, c=pctIn, cmap=custom_cmap, norm=norm, marker='.', s=1)
scatter = plt.scatter(recall, precision, c=pctIn, 
                     cmap=custom_cmap,
                     norm=norm,
                     marker='.', 
                     s=1)

# Rest of your plotting code
cbar = plt.colorbar(scatter, label='PctIn', aspect=50)
cbar.ax.tick_params(labelsize=14)  # Change tick font size
cbar.set_label('PctIn', size=16)   # Change label font size


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
plt.legend(fontsize=12, loc='upper right')  # or 'upper left', 'center', etc. # plt.legend(fontsize=12, bbox_to_anchor=(1.15, 1))

plt.savefig(outputName, dpi=600, bbox_inches='tight')


# plt.show()



