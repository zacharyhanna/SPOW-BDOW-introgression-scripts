import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
# putting the data into an array section
dataTable = sys.argv[1]

ancestry = []
coverage = []
category = []

def grab_data(file_in):
    with open(file_in, 'r') as infile:
        for line in infile:
            line = line.strip()
            splitline = line.split(",")
            category.append(splitline[0].replace("_"," "))
            ancestry.append(float(splitline[1]))
            coverage.append(float(splitline[2]))

grab_data(dataTable)

x = {'ances': ancestry,
'cover': coverage,
'categ': category}
print len(ancestry)
print len(coverage)
df = pd.DataFrame(x)
print df.ances
print df.cover



markers={
"Eastern Barred Owl": "o",
"Western Barred Owl": "o",
"Siskiyou Barred Owl": "o",
"Spotted Owl Pre": "o",
"Spotted Owl": "o",
"Sparred Owl": "o"}

colors ={
"Eastern Barred Owl": "#d55e00",
"Western Barred Owl": "#cc79a7",
"Siskiyou Barred Owl": "#e69f00",
"Spotted Owl Pre": "#009e73",
"Spotted Owl": "#0072b2",
"Sparred Owl": "#f0e442"}

for kind in markers:
    d = df[df.categ==kind]
    #plt.scatter(d.ances, d.cover)
    plt.scatter(d.ances, d.cover, c=colors[kind], marker=markers[kind], lw=0, s=60)
plt.xlim(0,1)
plt.ylim(0,7)
plt.grid(linestyle='-')
plt.ylabel("Coverage (X)")
plt.xlabel("Percentage Spotted Owl Ancestry")

handles =[]
labels = []
#line1 = matplotlib.lines.Line2D(color="blue",label='blue data')
point1 = plt.scatter([],[],color="#cc79a7",marker='o',label='Barred Owl (western)',lw=0, s=60)
point2 = plt.scatter([],[],color="#e69f00",marker='o',label='Barred Owl (Siskiyou)',lw=0, s=60)
point3 = plt.scatter([],[],color="#d55e00",marker='o',label='Barred Owl (eastern)',lw=0, s=60)
point4 = plt.scatter([],[],color="#009e73",marker='o',label='Spotted Owl (pre-contact)',lw=0, s=60)
point5 = plt.scatter([],[],color="#0072b2",marker='o',label='Spotted Owl (post-contact)',lw=0, s=60)
point6 = plt.scatter([],[],color="#f0e442",marker='o',label='Hybrid',lw=0, s=60)
#plt.legend(handles, labels)
plt.legend(handles=[point1,point2,point3,point4,point5,point6],scatterpoints=1)
plt.savefig('foo.eps', bbox_inches='tight')
plt.show()
