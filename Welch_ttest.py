import sys
from scipy import stats
import numpy as np
# input file is a comma-separated value (csv) file with the population in the
# first field and the ancestry value in the second field
dataFile = sys.argv[1]

pop_dict = {}
def write_out_val(in_val, filename, title):
    global pop_dict
    outfile = filename+".txt"
    header_val = "calculated t-statistic, two-tailed p-value\n"
    with open(outfile, "a") as file:
        file.write(title)
        file.write(header_val)
        file.write(str(in_val))
        file.write("\n")

def grab_data(file_in):
    with open(file_in, 'r') as infile:
        for line in infile:
            line = line.strip()
            splitline = line.split(",")
            if splitline[1] in pop_dict:
                pop_dict[str(splitline[1])].append(float(splitline[0]))
            else:
                pop_dict[str(splitline[1])]=[float(splitline[0])]
                #pop_dict[str(splitline[1])].append(float(splitline[0]))
            #print splitline[0]
            #print splitline[1]


grab_data(dataFile)

vals_ary3 = np.asarray(pop_dict["Western_Barred_Owl"])
vals_ary4 = np.asarray(pop_dict["Siskiyou_Barred_Owl"])

ttest_val2 = stats.ttest_ind(vals_ary4,vals_ary3, equal_var = False) # perform Welch's t-test
title2 = "Siskiyou vs Western Barred Owls\n"
write_out_val(ttest_val2, "Welch_ttest_out", title2)

We_Sis_combo = pop_dict["Western_Barred_Owl"] + pop_dict["Siskiyou_Barred_Owl"]
vals_ary5 = np.asarray(We_Sis_combo)
vals_ary6 = np.asarray(pop_dict["Eastern_Barred_Owl"])

ttest_val3 = stats.ttest_ind(vals_ary5,vals_ary6, equal_var = False) # perform Welch's t-test
title3 = "AllWestern vs Eastern Barred Owls\n"
write_out_val(ttest_val3, "Welch_ttest_out", title3)

vals_ary7 = np.asarray(pop_dict["Spotted_Owl"])
vals_ary8 = np.asarray(pop_dict["Spotted_Owl_Pre"])

ttest_val4 = stats.ttest_ind(vals_ary8,vals_ary7, equal_var = False) # perform Welch's t-test
title4 = "Pre vs. Post Spotted Owls\n"
write_out_val(ttest_val4, "Welch_ttest_out", title4)

SPOWcombo = pop_dict["Spotted_Owl"] + pop_dict["Spotted_Owl_Pre"]
vals_ary9 = np.asarray(SPOWcombo)

All_BADOs = pop_dict["Western_Barred_Owl"] + pop_dict["Siskiyou_Barred_Owl"] + pop_dict["Eastern_Barred_Owl"]
vals_ary10 = np.asarray(All_BADOs)

ttest_val7 = stats.ttest_ind(vals_ary10,vals_ary9, equal_var = False) # perform Welch's t-test

title7 = "AllBarred vs AllSpotted Owls\n"
write_out_val(ttest_val7, "Welch_ttest_out", title7)
