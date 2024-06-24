import matplotlib.pyplot as plt
import numpy as np

print('a')
labels = ['G1', 'G2', 'G3', 'G4', 'G5']
men_means = [20, 34, 30, 35, 27]
women_means = [25, 32, 34, 20, 25]
aa = [13,4,44,55,66]
bb = [45,4,44,55,66]
groups = [men_means, women_means, aa, bb]
groupsLabe = ['aa','bb','cc', 'dd']
x = np.arange(len(labels))  # the label locations
print(x)
width = 1.5/(len(labels)+len(groups))  # the width of the bars
span = 0.3/(len(labels)+len(groups)-1)

fig, ax = plt.subplots(figsize=(5,5))
axx =len(groups)*(width)+(len(groups)-1)*span/2
a= x  + width/2 -(len(groups)*(width)+(len(groups)-1)*span/2)
xPos = [axx/2-span/2 + a+(width+span)*i for i in range(0,len(groups)+0)]
print(xPos)
for i in range(len(groups)):
    rects1 = ax.bar(xPos[i], groups[i], width, label=groupsLabe[i])
    #rects2 = ax.bar(xPos, groups[i], width, label='Women')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Scores')
ax.set_title('Scores by group and gender')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

ax.bar_label(rects1, padding=3)
#ax.bar_label(rects2, padding=3)

fig.tight_layout()


plt.savefig("a.png",format='png')