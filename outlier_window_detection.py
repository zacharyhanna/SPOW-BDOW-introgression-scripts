import sys
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom
import csv
import math
dataTable = sys.argv[1]

sample_dict = {
4:("ASG007","Barred Owl (Siskiyou)"),
5:("ASG017","Barred Owl (western)"),
6:("ASG037","Barred Owl (Siskiyou)"),
7:("CAS87569","Spotted Owl (pre-contact)"),
8:("CMC40824","Barred Owl (eastern)"),
9:("CMCB40819","Barred Owl (eastern)"),
10:("CMCB41533","Barred Owl (eastern)"),
11:("CMCB41566","Barred Owl (eastern)"),
12:("CU51478","Barred Owl (eastern)"),
13:("HaigHybrid4","Spotted Owl (post-contact)"),
14:("Hoopa20005","Barred Owl (western)"),
15:("Hoopa20011","Barred Owl (western)"),
16:("Hoopa20014","Barred Owl (western)"),
17:("Hoopa20017","Barred Owl (western)"),
18:("Hoopa20018","Barred Owl (western)"),
19:("Hoopa20019","Barred Owl (western)"),
20:("JMR920","Barred Owl (Siskiyou)"),
21:("JPD386","Barred Owl (Siskiyou)"),
22:("LCW405","Barred Owl (western)"),
23:("LCW443","Barred Owl (western)"),
24:("LCW491","Barred Owl (western)"),
25:("MEF404","Barred Owl (Siskiyou)"),
26:("MEF432","Barred Owl (Siskiyou)"),
27:("MEF435","Barred Owl (Siskiyou)"),
28:("MEF457","Barred Owl (eastern)"),
29:("MK1012","Barred Owl (Siskiyou)"),
30:("MK1020","Barred Owl (Siskiyou)"),
31:("MK968","Barred Owl (Siskiyou)"),
32:("MK987","Barred Owl (Siskiyou)"),
33:("MK994","Barred Owl (western)"),
34:("MK998","Barred Owl (western)"),
35:("NSO138799040","Spotted Owl (post-contact)"),
36:("NSO168709365","Spotted Owl (post-contact)"),
37:("Sequoia","Spotted Owl (pre-contact)"),
38:("UWBM53433","Spotted Owl (post-contact)"),
39:("UWBM62061","Spotted Owl (pre-contact)"),
40:("UWBM65055","Barred Owl (western)"),
41:("UWBM67015","Barred Owl (western)"),
42:("UWBM74078","Barred Owl (western)"),
43:("UWBM76815","Barred Owl (western)"),
44:("UWBM79007","Barred Owl (western)"),
45:("UWBM79049","Barred Owl (western)"),
46:("UWBM79141","Barred Owl (western)"),
47:("UWBM91379","Spotted Owl (post-contact)"),
48:("UWBM91380","Spotted Owl (post-contact)"),
49:("UWBM91382","Barred Owl (western)"),
50:("UWBM91392","Spotted Owl (post-contact)"),
51:("UWBM91393","Spotted Owl (post-contact)"),
52:("UWBM91408","Spotted Owl (post-contact)"),
53:("ZRH455","Spotted Owl (post-contact)"),
54:("ZRH602","Barred Owl (western)"),
55:("ZRH604","Barred Owl (western)"),
56:("ZRH607","Hybrid"),
57:("ZRH962","Hybrid")
}

color_dict ={
"Barred Owl (eastern)": "#d55e00",
"Barred Owl (western)": "#cc79a7",
"Barred Owl (Siskiyou)": "#e69f00",
"Spotted Owl (pre-contact)": "#009e73",
"Spotted Owl (post-contact)": "#0072b2",
"Hybrid": "#f0e442"}

def lookoverlap(dict_in):
    dict_scafs = {}
    lendict = 0
    for samp in dict_in: #loop over all samples (keys in the first-level dictionary)
        #print lendict
        for scaf in dict_in[samp]: #loop over all scaffolds (keys in this sub-dictionary) that have windows for this sample
            for win_coords in dict_in[samp][scaf]: #loop down the list of window coordinates that is the value for this scaf key
                if scaf in dict_scafs:
                    if (win_coords[0],win_coords[1]) in dict_scafs[scaf]:
                        dict_scafs[scaf][(win_coords[0],win_coords[1])].append(samp)
                    else:
                        dict_scafs[scaf][(win_coords[0],win_coords[1])]=[samp]
                else:
                    dict_scafs[scaf]={(win_coords[0],win_coords[1]):[samp]}
    return dict_scafs

def test_pval(prob,samp,calc_vers):
    if sample_dict[samp][1] == "Barred Owl (Siskiyou)" or sample_dict[samp][1] == "Barred Owl (western)" or sample_dict[samp][1] == "Barred Owl (eastern)":
        if calc_vers == 1: # hard cutoff calculation
            if prob >= 0.4:
                return True
        elif calc_vers == 2:
            compare = prob
    elif sample_dict[samp][1] == "Spotted Owl (post-contact)" or sample_dict[samp][1] == "Spotted Owl (pre-contact)":
        if calc_vers == 1:
            if prob <= 0.6:
                return True
        elif calc_vers == 2:
            compare = 1-prob
    if calc_vers ==2:
        if compare < 0.05:
            return False
        else:
            return True

def calc_pval(avg_anc,reads,samp,calc_vers):
    if calc_vers == 1:
        prob = avg_anc
        return test_pval(prob,samp,1)
    elif calc_vers == 2:
        prob_success = 0.538 # taken from F1 hybrid average ancestry
        num_suc = avg_anc * reads
        # binom.cdf(successes, attempts, chance_of_success_per_attempt)
        prob = binom.cdf(int(round(num_suc)),reads,prob_success) # rounding number of successes if is not an integer
        #print prob
        #pval = prob_binom_cdf(int(round(num_suc)),prob_success,reads)
        return test_pval(prob,samp,2)

def grab_data(file_in):
    wind_dict = {}
    numwind_dict = {}
    numoutliers_dict = {}
    with open(file_in, 'r') as infile:
        for line in infile:
            line = line.strip()
            splitline = line.split("\t")
            if splitline[0][0] == "#": # get rid of lines that begin with a pound sign (scaffold summary lines)
                continue
            else:
                scaf = splitline[0]
                st_idx = int(splitline[1])
                end_idx = int(splitline[2])
                if end_idx-st_idx < 49999: # make sure window is at least 50kb
                    continue
                else:
                    for samp in range(4,len(splitline)):
                        if samp != 10 and samp != 13 and samp != 37 and samp != 56 and samp != 57: # we don't need the hybrid samples, the Haig Hybrid, or the reference samples
                            splitsamp = splitline[samp].split(":")
                            if len(splitsamp) == 3: #make sure have ancestry, snp, and read values
                                if splitsamp[0] == "-1" or splitsamp[0] == "-": # make sure ancestry value is not "-1" or "-"
                                    continue
                                else:
                                    avg_anc = float(splitsamp[0])
                                    snps = int(splitsamp[1])
                                    reads = int(splitsamp[2])
                                    if snps > 9:
                                        if samp in numwind_dict:
                                            numwind_dict[samp]+=1
                                        else:
                                            numwind_dict[samp]=1
                                        if calc_pval(avg_anc,reads,samp,1) == True: # check for at least 10 snps and cumulative binomial probability
                                            #print "Yes", samp, sample_dict[samp], scaf, snps, reads, avg_anc
                                            if samp in numoutliers_dict:
                                                numoutliers_dict[samp]+=1
                                            else:
                                                numoutliers_dict[samp]=1
                                            if samp in wind_dict:
                                                if scaf in wind_dict[samp]:
                                                    if wind_dict[samp][scaf][-1][0] < st_idx <= wind_dict[samp][scaf][-1][1]+1:
                                                        wind_dict[samp][scaf][-1][1]=end_idx
                                                    else:
                                                        wind_dict[samp][scaf].append([st_idx,end_idx])
                                                else:
                                                    wind_dict[samp][scaf]=[]
                                                    wind_dict[samp][scaf].append([st_idx,end_idx])
                                            else:
                                                wind_dict[samp]={}
                                                wind_dict[samp][scaf]=[]
                                                wind_dict[samp][scaf].append([st_idx,end_idx])
    return wind_dict, numwind_dict, numoutliers_dict

def make_hist_dict(wind_dict):
    hist_dict = {}
    for sample in wind_dict:
        hist_dict[sample]=[]
        for scaffold in wind_dict[sample]:
            for element in wind_dict[sample][scaffold]:
                len_wind = int(element[1])-int(element[0])+1
                if len_wind >= 50000:
                    hist_dict[sample].append(len_wind)
    return hist_dict

def hist_overlap_dict(in_dict):
    dict_num_samps = {}
    total_windows = 0.0
    for scaf in in_dict:
        for wind in in_dict[scaf]:
            total_windows += 1
            num_samps = len(in_dict[scaf][wind])
            if num_samps in dict_num_samps:
                dict_num_samps[num_samps]+=1
            else:
                dict_num_samps[num_samps]=1
    return dict_num_samps, total_windows

def hist_scaf_overlap(in_dict):
    dict_num_wind = {}
    total_scafs = 0.0
    for scaf in in_dict:
        total_scafs += 1
        num_winds = len(in_dict[scaf])+1
        if num_winds in dict_num_wind:
            dict_num_wind[num_winds] +=1
        else:
            dict_num_wind[num_winds]=1
        dict_num_wind[scaf]=num_winds
    return dict_num_wind, total_scafs


def write_csv(in_dict, filename):
    outfile = filename+".csv"
    with open(outfile,'w') as f:
        w = csv.writer(f)
        w.writerows(in_dict.items())

def write_out_dict(in_dict, filename):
    outfile = filename+".pydict"
    with open(outfile, "w") as file:
        file.write(str(in_dict))

filled_wind_dict = grab_data(dataTable)
filled_hist_dict = make_hist_dict(filled_wind_dict[0])
#write_csv(filled_wind_dict[0],"filled_wind_dict_0")
#write_csv(filled_wind_dict[1],"filled_wind_dict_1")
#write_csv(filled_wind_dict[2],"filled_wind_dict_2")
#write_out_dict(filled_wind_dict,"filled_wind_dicts")

# graph analyzed windows vs outlier windows
dict_analyzed_winds = filled_wind_dict[1]
dict_outlier_winds = filled_wind_dict[2]
def graph_w_v_w(analyzed,outlier):
    analyzed_ls = []
    outlier_ls = []
    samp_ls = []
    for samp in analyzed:
        samp_ls.append(samp)
        analyzed_ls.append(analyzed[samp])
        if samp in outlier:
            outlier_ls.append(outlier[samp])
        else:
            outlier_ls.append(0)
    f4 = plt.figure()
    fig4, ax4 = plt.subplots()

    rects4 = ax4.scatter(analyzed_ls,outlier_ls)
    #for i, txt in enumerate(samp_ls):
    #    ax4.annotate(
    #        txt,
    #        xy=(analyzed_ls[i],outlier_ls[i]),
    #        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    #ax1.set_xticks(xtick_range)
    ax4.set_xlabel("Number of Analyzed Windows")
    ax4.set_ylabel("Number of Outlier Windows")
    ax4.set_title("Number of Outlier vs Analyzed Windows in Each Sample")
    ax4.set_ylim([0,100])
    #for rect1 in rects1:
        #height1 = rect1.get_height()
        #scaled_back_ht = int(height*filled_hist_overlap[1])
        #ax1.text(rect1.get_x() + rect1.get_width()/2., height1, '%s' % scaled_back_ht, ha='center',va='bottom')
    fig4.savefig('outliers_vs_analyzed_windows.png', bbox_inches='tight',pad_inches=0)#,pad_inches=0)

graph_an_v_out = graph_w_v_w(dict_analyzed_winds,dict_outlier_winds)


filled_overlap = lookoverlap(filled_wind_dict[0])
#write_csv(filled_overlap,"filled_overlap")
#write_out_dict(filled_overlap,"filled_overlap")

filled_hist_overlap = hist_overlap_dict(filled_overlap)
#write_csv(filled_hist_overlap[0],"filled_hist_overlap_0")
#write_out_dict(filled_hist_overlap,"filled_hist_overlap")

filled_scaf_overlap = hist_scaf_overlap(filled_overlap)
#write_csv(filled_scaf_overlap[0],"filled_scaf_overlap_0")
#write_out_dict(filled_scaf_overlap,"filled_scaf_overlap")

f2x = []
f2y = []
for f2key in filled_hist_overlap[0]:
    f2x.append(f2key)
    f2y.append(filled_hist_overlap[0][f2key]/filled_hist_overlap[1]) #make the height a proportion of total overlapping windows, rathern than a raw number
f2x_sorted = f2x.sort()
f2 = plt.figure()
fig2, ax2 = plt.subplots()
rects = ax2.bar(f2x,f2y,align='center',width=1)
ax2.set_ylim([0,1])
if f2x and f2y:
    xtick_range = range(min(f2x), max(f2x)+1)
    ax2.set_xticks(xtick_range)
ax2.set_xlabel("Number of Samples")
ax2.set_ylabel("Proportion of Outlier Windows")
ax2.set_title("Number of Samples in an Outlier Window")
for rect in rects:
    height = rect.get_height()
    scaled_back_ht = int(height*filled_hist_overlap[1])
    ax2.text(rect.get_x() + rect.get_width()/2., height, '%s' % scaled_back_ht, ha='center',va='bottom')

fig2.savefig('number_samples_in_outliers.png', bbox_inches='tight',pad_inches=0)#,pad_inches=0)

statsfile = "outlier_window_stats.txt"
outstring1 = "number of scaffolds with outlier windows = " + str(len(filled_overlap)+1) + "\n"
with open(statsfile, "w") as file:
    file.write(outstring1)

ls_win_array = []
sampkey_array = []
hist_values = []

len_filled_hist = len(filled_hist_dict)
zpoints = [item * 10 for item in range(len_filled_hist)]

pop_order ={
"Barred Owl (eastern)": 2,
"Barred Owl (western)": 4,
"Barred Owl (Siskiyou)": 3,
"Spotted Owl (post-contact)": 1,
"Spotted Owl (pre-contact)": 5,
"Hybrid": 6}
sample_order_num = {}
for sample1 in sample_dict: #ordering the samples
    sample_order_num[sample1]=pop_order[sample_dict[sample1][1]]

# sort the samples
sorted_sample_ls = sorted(sample_dict, key=sample_order_num.__getitem__)
sorted_sample_ls
for sampkey in sorted_sample_ls:
    if sampkey in filled_hist_dict:
        yedges, xedges = np.histogram(np.asarray(filled_hist_dict[sampkey]), bins=[50000,100000,150000,200000])
        new_yedges = np.append(yedges,[0]) # add an extra zero onto yedges
        for eley in new_yedges:
            if eley == 0: # take out zero values from histogram array
                pass
        hist_values.append([sampkey,xedges,new_yedges,sample_dict[int(sampkey)][1],color_dict[sample_dict[int(sampkey)][1]],sample_order_num[sampkey]])

for nodivval in hist_values:
    nodivval[2] = np.divide(nodivval[2],float(filled_wind_dict[1][nodivval[0]]))
new_hist_vals = sorted(hist_values, key=lambda x: (x[5],x[2][0])) #sorts hist_values list by x[5] (sample_order_num[sampkey]]) first and then x[2][0] (the first element of new_yedges, which is the height of the 50kb bin)
new2_hist_vals = new_hist_vals
start_z_ax = 0
for newele in new2_hist_vals:
    newele.append(start_z_ax)
    list_trash = []
    for new2, zet in enumerate(newele[2]):
        if zet == 0 or zet == float(0):
            list_trash.append(new2)
    new1elem = np.delete(newele[1],list_trash)
    newele[1] = new1elem
    new2elem = np.delete(newele[2],list_trash)
    newele[2] = new2elem
    start_z_ax += 10

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
for elem in new2_hist_vals:
    ax.bar(elem[1],elem[2],elem[6],zdir='y',label=elem[3],color=elem[4],edgecolor=elem[4],width=50000,alpha=0.7)
ax.set_zlabel("Proportion of Total Windows",labelpad=10)
ax.set_ylabel("Samples",labelpad=-5)
ax.set_xlabel("Outlier Window Lengths (nt)")#,labelpad=10
x_axis_labs = ['50,000','100,000','150,000']
ax.set_xlim([50000, 200000])
ax.set_xticks([50000,100000,150000])
ax.set_xticklabels(x_axis_labs)
ax.tick_params(axis='x',which='major',pad=-2)
ax.set_yticklabels([])

point1 = plt.Rectangle((0,0),1,1,color="#cc79a7",edgecolor="#cc79a7",label='Barred Owl (western)')
point2 = plt.Rectangle((0,0),1,1,color="#e69f00",edgecolor="#e69f00",label='Barred Owl (Siskiyou)')
point3 = plt.Rectangle((0,0),1,1,color="#d55e00",edgecolor="#d55e00",label='Barred Owl (eastern)')
point5 = plt.Rectangle((0,0),1,1,color="#0072b2",edgecolor="#0072b2",label='Spotted Owl (post-contact)')

ax.legend([point1,point2,point3,point5],['Barred Owl (western)','Barred Owl (Siskiyou)','Barred Owl (eastern)','Spotted Owl (post-contact)'],loc=2,fontsize='small')
fig.tight_layout()
fig.savefig('outlier_windows_histograms_by_sample.png', bbox_inches='tight',pad_inches=0)#,pad_inches=0)
