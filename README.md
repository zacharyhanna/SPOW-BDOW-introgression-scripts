# SPOW-BDOW-introgression-scripts
We here provide the scripts that we developed for analyzing introgression in whole genome sequences obtained from XXX. 

## Contents
* [Introduction](#introduction)  
* [Sliding Window Analyses Pipeline](#sliding-window-analyses-pipeline)  
  * [1. Filter raw vcf file](#1-filter-raw-vcf-file)
  * [2. Site coverage calculation](#2-site-coverage-calculation)
  * [3. Exclude sites with excessive coverage](#3-exclude-sites-with-excessive-coverage)
  * [4. Calculation of coverage depth per sample](#4-calculation-of-coverage-depth-per-sample)
  * [5 Allele depth calculation](#5-allele-depth-calculation)
    * [5.1 Allele depth only calculation](#51-allele-depth-only-calculation)
    * [5.2 Extended allele depth calculation](#52-extended-allele-depth-calculation)
  * [6. Sliding window calculation](#6-sliding-window-calculation)
    * [6.1 sliding_window.sh](#61-sliding\_windowsh)
    * [6.2 ext_fmt_sliding_window_reads.sh](#62-ext\_fmt\_sliding\_window\_readssh)
  * [7. Compute means and standard deviations on allele depth file](#7-compute-means-and-standard-deviations-on-allele-depth-file)
    * [7.1 Get sample names from vcf](#71-get-sample-names-from-vcf)
    * [7.2 Merge sample names from vcf with means and standard deviations](#72-merge-sample-names-from-vcf-with-means-and-standard-deviations)
  * [8. Check for outlier windows](#8-check-for-outlier-windows)
* [Welch's _t_-test](#welchs-t-test)
* [Nucleotide diversity and FST calculation](#nucleotide-diversity-and-fst-calculation)
* [Citing the Repository](#citing-the-repository)  
* [References](#references)  

## Introduction
These scripts take a variant call format (vcf) file as input. We created our vcf file using UnifiedGenotyper from the Genome Analysis Toolkit (GATK) version 3.4-46 (DePristo et al. 2011; McKenna et al. 2010; Van der Auwera et al. 2013). Due to the specifics of our sample set and our analyses, you will not be able to directly use most of these scripts without making some modifications to the code. VCF files produced by other variant callers may require modifications to these scripts for correct parsing of the variant files. We are providing our code here for the purposes of documentation and with the hope that some of our methods may prove useful to others in their own genome-scale analyses.  
#### Specifics to our analyses that have affected the code of all scripts
- VCF file output from UnifiedGenotyper: raw_variants.vcf  
- Fields 16 and 43 in the vcf corresponded with the reference Barred and Spotted Owl samples, respectively.  

## Sliding Window Analyses Pipeline

### 1. Filter raw vcf file
Usage example:  
```
$ ./vcf_qual_filter.sh raw_variants.vcf > filtered_variants.vcf  
```  
vcf_qual_filter.sh shell script requirements:  
MAWK, if available - we used MAWK version 1.2 (Brennan 1994).  
If MAWK is not available, uses GNU Awk (GAWK) - we used GAWK version 4.0.1 (Free Software Foundation 2012)
  
Specifics to our analyses that have affected the code:  
1. We examined only biallelic variant sites fixed for alternative alleles between the reference Spotted and Barred Owl samples.  
   * Field 43 in the vcf corresponded with the reference Spotted Owl sample.  
   * Field 16 in the vcf corresponded with the reference Barred Owl sample.  
2. Our reference genome had some contaminant or mitochondrial scaffolds, which we removed from our analyses.  
   * You could just take out the line that throws out the following scaffolds: C7961234, C7963448, C7970814, C8091874, scaffold3674.  
3. We only examined sites with a Phred-scaled probability >50 that a polymorphism exists at that site.  
   * This was the "QUAL" field, which was field 6 of the vcf 
4. We required that the Phred-scaled genotype quality, which states the confidence in the genotype of a particular sample, must be >=30 for both the Spotted Owl and Barred Owl reference samples.
5. We required that the Barred Owl reference sample had zero reads that supported the Spotted Owl allele at a given variant site.
   * We examined the unfiltered allele depth (field AD) to make this determination.
6. We required that the Spotted Owl reference sample had zero reads that supported the Barred Owl allele at a given variant site and >=10 reads in support of the Spotted Owl allele.
   * We examined the unfiltered allele depth (field AD) to make this determination.

### 2. Site coverage calculation
Calculate the mean total coverage at a site and the standard deviation (σ).  
We are using the filtered set of sites for calculation of the average coverage.  
  
Usage example:  
```
$ cat filtered_variants.vcf | ./dp_cov_script.sh 
```  
Example output:  
```
meanDP = 129.19,stdevDP = 34.2783,number of sites = 5821431  
```
### 3. Exclude sites with excessive coverage
We excluded sites with coverage in excess of the mean + 5σ (we only kept sites with coverage <301 nt in our data set), as suggested by the GATK documentation (https://software.broadinstitute.org/gatk/guide/article?id=3225).  
  
Usage example:  
```
$ cat filtered_variants.vcf | ./vcf_dp_filter.sh >filtered_variants2.vcf  
```
### 4. Calculation of coverage depth per sample
This step is not actually required in the pipeline, but we calculated the mean and standard deviation (σ) of the coverage depth for each sample across all of the sites in our final filtered set of SNPs.  
  
Usage example:  
```
$ cat filtered_variants2.vcf | ./DP_means_std_dev.sh | head -1 >filtered_variants2_dp_means_stdev.txt  
```
Example output:
```
$ cat filtered_variants2_dp_means_stdev.txt
0.745503,0.913129 0.34499,0.599295 0.996806,1.05703 1.07951,1.10409 0.627326,0.826752 1.0912,1.11578 15.5489,5.81112 1.52076,1.35264 0.960982,1.04229 0.35384,0.604826 0.44408,0.680869 0.297244,0.553016 0.515357,0.735895 0.421296,0.662121 0.354241,0.608053 0.277709,0.536011 0.947581,1.03082 0.0355873,0.189918 0.751939,0.899872 0.482803,0.711412 0.722469,0.87494 1.74519,1.40872 0.944969,1.01582 2.16733,1.5994 6.40654,3.13418 2.07245,1.60605 0.975839,1.04753 0.11785,0.348755 0.713162,0.880505 0.59747,0.812452 0.853734,0.985254 1.38034,1.22845 1.98464,1.49242 60.8151,15.9385 0.168379,0.418069 0.495495,0.731765 0.377467,0.625025 0.72346,0.897134 0.510975,0.735475 0.75375,0.903665 0.0810105,0.286488 0.59008,0.797603 0.759363,0.90252 0.448514,0.697064 0.0166811,0.13067 0.262981,0.551261 0.0439719,0.211284 0.261251,0.522411 1.10659,1.13334 0.850381,0.964539 3.49167,2.20528 4.54619,2.81694 1.62449,1.39522 2.20901,1.61496
```
Each space-separated field in the output is "mean coverage,standard deviation".  
### 5. Allele depth calculation
There are two versions of this script whose output function with different downstream scripts.  
Both "AD_pct.sh" and "AD_pct_ex.sh" require:  
GNU Awk - we used GNU Awk version 4.0.1 (Free Software Foundation, 2012)  
cut (GNU coreutils) - we used cut (GNU coreutils) version 8.21 (Ihnat et al. 2013)
#### 5.1 Allele depth only calculation
Usage example:  
```
$ cat filtered_variants2.vcf | ./AD_pct.sh >ad_pct.txt  
```
Example output:
```
$ head -n 2 ad_pct.txt
scaffold18 506 -1 -1 0 1 -1 0 0 -1 0 -1 0 -1 -1 -1 -1 -1 0 -1 0 0 -1 0 0 0 0 0 -1 -1 -1 0 0 1 1 1 1 1 -1 0 -1 -1 -1 0 -1 -1 -1 -1 -1 -1 -1 1 -1 0 1 1
scaffold18 520 0 -1 -1 1 -1 0 0 -1 0 -1 0 -1 -1 -1 -1 -1 0 -1 0 0 -1 0 0 0 0 0 -1 -1 -1 0 0 1 -1 1 -1 1 -1 0 -1 -1 -1 0 -1 -1 -1 -1 -1 1 -1 1 -1 0 1 0.67
```
Each space-separated field after the scaffold and position fields is "Spotted Owl ancestry percentage".  
Ancestry percentage is between 0 and 1 and is the percentage of the sequences that support the Spotted Owl allele at this site fixed between our Spotted and Barred Owl reference sequences.  
"-1:0" is output for no data for a sample at a site 
#### 5.2 Extended allele depth calculation
This script returns the number of reads for a sample at a site in addition to the percentage ancestry.  
Usage example:  
```
$ cat filtered_variants2.vcf | ./AD_pct_ex.sh >ad_pct_ex.txt  
```  
Example output:  
```
$ head -n 2 ad_pct_ex.txt
scaffold18 506 -1:0 -1:0 0:1 1:2 -1:0 0:1 0:17 -1:0 0:1 -1:0 0:1 -1:0 -1:0 -1:0 -1:0 -1:0 0:2 -1:0 0:2 0:1 -1:0 0:2 0:2 0:4 0:6 0:2 -1:0 -1:0 -1:0 0:1 0:2 1:2 1:2 1:52 1:1 1:2 -1:0 0:1 -1:0 -1:0 -1:0 0:1 -1:0 -1:0 -1:0 -1:0 -1:0 -1:0 -1:0 1:3 -1:0 0:6 1:1 1:2 
scaffold18 520 0:1 -1:0 -1:0 1:2 -1:0 0:1 0:19 -1:0 0:3 -1:0 0:1 -1:0 -1:0 -1:0 -1:0 -1:0 0:2 -1:0 0:2 0:1 -1:0 0:2 0:2 0:4 0:7 0:1 -1:0 -1:0 -1:0 0:1 0:1 1:2 -1:0 1:52 -1:0 1:2 -1:0 0:1 -1:0 -1:0 -1:0 0:1 -1:0 -1:0 -1:0 -1:0 -1:0 1:1 -1:0 1:3 -1:0 0:6 1:1 0.67:3 
```
Each space-separated field after the scaffold and position fields is "Spotted Owl ancestry percentage:coverage for sample at site (number of sequences)".  
Ancestry percentage is between 0 and 1 and is the percentage of the sequences that support the Spotted Owl allele at this site fixed between our Spotted and Barred Owl reference sequences.  
"-1:0" is output for no data for a sample at a site  

### 6. Sliding window calculation  

#### 6.1 sliding_window.sh
This script will run a sliding window analysis on output produced from ad_pct.txt, which does not contain any read coverage information.  
  
Usage example:  
```
$ cat ad_pct.txt | ./sliding_window.sh 50000 >wnd_50k_noovlp.txt  
```
The above example calculates 50,000 bp windows with no overlap.  
  
"sliding_window.sh" shell script requirements:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  

##### Example of running full sliding window pipeline
As above, this example calculates 50,000 bp windows with no overlap.
```  
$ ./vcf_qual_filter.sh raw_variants.vcf | ./AD_pct.sh >ad_pct.txt  
$ cat ad_pct.txt | ./sliding_window.sh 50000 >wnd_50k_noovlp.txt  
```
The first pipe starts with the original vcf file and outputs to the ad_pct.txt file, which is then used in the second pipe for the sliding window analysis.  
These scripts enable you to pipe all the way through, but, since multiple runs of the sliding window routine are likely, we saved the ad_pct.txt file and worked from it for multiple sliding window analyses.  

Example output:  
```
$ head -n 3 wnd_50k_noovlp.txt
#scaffold12     552     292215  1415    0.074745 0.0392377 0.0733255 0.992282 0.0574987 0.0531025 0 0.0655141 0.0604773 0.997758 0.0667287 0.0647343 0.0300879 0.0442817 0.0614458 0.0856269 0.0607843 0.0909091 0.0390326 0.0925926 0.0659242 0.0550549 0.0649718 0.0442348 0.0648399 0.0712263 0.0426153 0.0657895 0.0610762 0.0487805 0.0492584 0.997751 0.999543 1 1 0.996109 0.0423661 0.0238786 0.0276817 0.0643652 0.075188 0.0613905 0.0581956 0.997938 1 0.0596817 1 0.997455 1 1 0.0544884 0.0670834 0.536692 0.497809 
scaffold12      1       50000   225     0.0737705 0.142857 0.0707965 1 0.0777778 0.0328639 0 0.0736196 0.0518519 1 0.166667 0.0980392 0.127119 0.0406977 0.0545455 0.0869565 0.0535714 0.2 0.0804598 0.095 0.0882353 0.0500936 0.0651927 0.0857143 0.0872328 0.0678392 0.0950704 0.222222 0.122222 0.111111 0.0731061 1 1 1 1 1 0.0555556 0.0507246 0.119565 0.0571429 0.0967742 0.0555556 0.0840336 1 - 0.0461538 1 1 1 1 0.0670194 0.0696182 0.491082 0.528932 
scaffold12      50001   100000  199     0.0258621 0.0447761 0.0396175 1 0.040293 0.0307692 0 0.0166667 0.0121655 1 0.0338983 0.0428571 0.0175439 0.0294118 0.0220588 0.0434783 0.010582 0 0.0373832 0.0178571 0.0245098 0.0418638 0.0420168 0.0368901 0.0356963 0.045098 0.00833333 0 0.0195195 0.05 0.0645161 0.995679 1 1 1 1 0 0.025 0.0238095 0.00724638 0.08 0.039604 0.0304348 1 1 0.0181818 1 1 1 1 0.00970018 0.0297659 0.512691 0.544737
```
Output is tab-separated.
The line beginning with "#" is a summary line for the scaffold.  
  
Field 1 | Field 2 | Field 3 | Field 4 | Fields 5-  
--- | --- | --- | --- | ---  
scaffold / contig | position of first SNP in scaffold | position of last SNP in the scaffold | number of SNPs in window (across all samples) | Average ancestry for sample across all windows in scaffold  
  
The successive lines provide:  

Field 1 | Field 2 | Field 3 | Field 4 | Fields 5-  
--- | --- | --- | --- | ---  
scaffold / contig | start position of window | end position of window | number of SNPs in window (across all samples) | Average ancestry for sample in window  

##### Example of a running a different sliding window
This example calculates 50,000 base windows sliding 5,000 bases at a time.  
```  
$ ./sliding_window.sh ad_pct.txt 50000 5000 >wnd_50k_5k_slide.txt  
```

#### 6.2 ext_fmt_sliding_window_reads.sh
This is modified version of the sliding_window.sh script that will work with the output from AD_pct_ex.sh. The sliding window analysis needs to be conducted with this script to enable outlier detection with the outlier_window_detection.py script. Usage is identical to that of sliding_window.sh.  
  
Usage example:  
```
$ cat ad_pct_ex.txt | ./ext_fmt_sliding_window_reads.sh >wnd_50k_noovlp_ext.txt  
```
The above example calculates 50,000 bp windows with no overlap.  
  
"ext_fmt_sliding_window_reads.sh" shell script requirements:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014) 

Example output:  
```
$ head -n 3 wnd_50k_noovlp_ext.txt
#scaffold12	552	292215	1415	0.075:817:4	0.039:446:0	0.073:851:7	0.99:907:16	0.058:629:3	0.053:924:13	0:1401:210	0.066:1099:19	0.06:824:11	1:446:5	0.067:537:5	0.065:345:2	0.03:493:11	0.044:478:7	0.061:415:7	0.086:327:1	0.061:850:19	0.091:33:1	0.039:696:12	0.093:576:9	0.066:651:6	0.055:1078:22	0.065:885:21	0.044:1234:38	0.065:1378:101	0.071:1236:24	0.043:896:20	0.066:152:0	0.061:796:13	0.049:697:10	0.049:854:8	1:978:22	1:1166:21	1:1401:693	1:330:7	1:514:11	0.042:417:7	0.024:862:19	0.028:578:6	0.064:659:17	0.075:133:1	0.061:676:11	0.058:726:5	1:485:6	1:11:0	0.06:377:4	1:83:2	1:393:3	1:920:11	1:780:9	0.054:1332:45	0.067:1355:85	0.54:1084:22	0.5:1216:32	
scaffold12	1	50000	225	0.074:122:170	0.14:56:68	0.071:113:176	1:138:220	0.078:90:107	0.033:142:277	0:225:3672	0.074:163:312	0.052:135:213	1:77:96	0.17:71:92	0.098:34:41	0.13:59:83	0.041:86:124	0.055:55:62	0.087:46:51	0.054:112:178	0.2:5:5	0.08:87:143	0.095:100:135	0.088:85:112	0.05:178:379	0.065:147:244	0.086:210:567	0.087:217:1541	0.068:199:557	0.095:142:234	0.22:18:21	0.12:135:202	0.11:99:139	0.073:132:210	1:143:271	1:189:540	1:225:13079	1:50:60	1:72:87	0.056:54:69	0.051:138:223	0.12:92:121	0.057:105:141	0.097:31:39	0.056:117:157	0.084:119:170	1:78:98	0:0:0	0.046:65:81	1:9:9	1:66:72	1:152:278	1:132:196	0.067:216:952	0.07:217:944	0.49:177:373	0.53:198:518	
scaffold12	50001	100000	199	0.026:116:166	0.045:67:87	0.04:122:228	1:126:208	0.04:91:128	0.031:130:208	0:199:3446	0.017:150:332	0.012:137:226	1:61:78	0.034:59:70	0.043:35:40	0.018:57:65	0.029:68:98	0.022:68:77	0.043:46:51	0.011:126:218	0:4:4	0.037:107:159	0.018:84:119	0.025:102:156	0.042:151:334	0.042:119:201	0.037:162:453	0.036:196:1327	0.045:170:493	0.0083:120:195	0:26:26	0.02:111:189	0.05:80:123	0.065:124:222	1:162:299	1:175:446	1:199:12419	1:52:63	1:79:96	0:57:71	0.025:120:190	0.024:84:107	0.0073:92:119	0.08:25:28	0.04:101:167	0.03:115:150	1:58:69	1:6:6	0.018:55:75	1:10:13	1:34:45	1:120:216	1:109:166	0.0097:189:700	0.03:193:845	0.51:153:340	0.54:171:465	
```
Output is tab-separated.
The line beginning with "#" is a summary line for the scaffold.

Field 1 | Field 2 | Field 3 | Field 4 | Fields 5-  
--- | --- | --- | --- | ---  
scaffold / contig | position of first SNP in scaffold | position of last SNP in the scaffold | number of SNPs in window (across all samples) | Average ancestry for sample across all windows in the scaffold : total number of SNPs for this sample across all windows in the scaffold : this read subfield is meaningless in the summary line   
  
The successive lines provide:  

Field 1 | Field 2 | Field 3 | Field 4 | Fields 5-  
--- | --- | --- | --- | ---  
scaffold / contig | start position of window | end position of window | number of SNPs in window (across all samples) | Average ancestry for sample in window : total number of SNPs for this sample in the window : number of reads across all variant sites in the scaffold  

### 7. Compute means and standard deviations on allele depth file
AD_pct.txt : allele depth file from above examples  
We are keeping the output for further use in the file "means_stdevs_ad.txt".  
  
Usage example:  
```
$ ./compute_ad_mean_stdev.sh ad_pct.txt >means_stdevs_ad.txt    
```  
Example output:  
```
$ cat means_stdevs_ad.txt
0.0685271,0.247087 0.0706799,0.25342 0.0689756,0.246139 0.995341,0.0658135 0.0655275,0.242632 0.0680178,0.243384 0,0 0.0662537,0.238591 0.068544,0.245222 0.995418,0.0664862 0.069415,0.250616 0.0706324,0.253674 0.0703329,0.25157 0.0702486,0.252144 0.0691084,0.250728 0.0697873,0.252425 0.0689068,0.246394 0.0715925,0.257479 0.0706598,0.250252 0.0705867,0.25238 0.070515,0.25017 0.0701403,0.242641 0.0697367,0.247418 0.0701523,0.240936 0.0695239,0.227327 0.0692743,0.240828 0.0683806,0.245452 0.0702735,0.254588 0.0700564,0.249522 0.0682737,0.247658 0.0727917,0.252138 0.994586,0.0691271 0.994712,0.0669258 1,0 0.995202,0.0685055 0.990622,0.0949053 0.0705875,0.253107 0.0689451,0.247426 0.0695106,0.250165 0.0690981,0.248127 0.0700916,0.254692 0.0703358,0.250722 0.0703353,0.249793 0.995633,0.0649109 0.994928,0.0709095 0.0671105,0.248352 0.995189,0.069063 0.995029,0.0695 0.99575,0.062973 0.995606,0.064032 0.0695541,0.235453 0.069667,0.23387 0.359312,0.427612 0.53809,0.376608
```
Each space-separated field in the output is "mean Spotted Owl ancestry percentage,standard deviation".  
  
The shell script requires:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  

#### 7.1 Get sample names from vcf
Usage example:  
```
$ grep -v "^#" raw_variants.vcf -B1 | head -1  
```  
The command requires:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  
head (GNU coreutils) - we used head (GNU coreutils) version 8.21 (Ihnat et al. 2013)  

#### 7.2 Merge sample names from vcf with means and standard deviations
Usage example:  
```
$ awk 'NR==1{b=1;for(i=10;i<=NF;i++)nm[b++]=$i}
     NR>1{for(i=1;i<=NF;i++){split($i,ms,",");
          mminus = (ms[1]-ms[2] > 0) ? ms[1]-ms[2] : 0;
          mplus  = (ms[1]+ms[2] > 1) ? 1 : ms[1]+ms[2]
          print i, nm[i], ms[1],ms[2],mminus,mplus}}' <(grep -v "^#" raw_variants.vcf -B1 |head -1) means_stdevs_ad.txt | sort -k3,3n | awk '{printf"%2s\t%12s\t",$1,$2;for(i=3;i<=NF;i++)printf "%-9s\t",$i;print ""}' | cat <(echo -e "Sample#\t Sample_Name\tmean\t\tstdev\t\tmean-stdev\tmean+stdev") -
```  
The command requires:  
cat (GNU coreutils) - we used cat (GNU coreutils) version 8.21 (Granlund & Stallman 2013)  
echo (GNU coreutils) - we used echo (GNU coreutils) version 8.21 (Fox & Ramey 2013)  
GNU Awk - we used GNU Awk version 4.0.1 (Free Software Foundation, 2012)  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  
head (GNU coreutils) - we used head (GNU coreutils) version 8.21 (Ihnat et al. 2013)  
sort (GNU coreutils) - we used sort (GNU coreutils) version 8.21 (Haertel & Eggert 2013)  

### 8. Check for outlier windows
Usage example:  
```
$ python outlier_window_detection.py wnd_50k_noovlp_ext.txt  
```  
The script produces four outputs (three graphical and one text file):  
  
**outlier_windows_histograms_by_sample.png** : graphs a histogram of outlier window lengths for each sample
![alt text](https://github.com/zrhanna/SPOW-BDOW-introgression-scripts/blob/master/images/outlier_windows_histograms_by_sample.png)  
  
**number_samples_in_outliers.png** : graphs a histogram of the number of outlier windows that are either unique or shared by more than one sample
![alt text](https://github.com/zrhanna/SPOW-BDOW-introgression-scripts/blob/master/images/number_samples_in_outliers.png)  
  
**outliers_vs_analyzed_windows.png** : for each sample, this graph plots the number of outlier windows versus the number of windows analyzed  
![alt text](https://github.com/zrhanna/SPOW-BDOW-introgression-scripts/blob/master/images/outliers_vs_analyzed_windows.png)  
  
**outlier_window_stats.txt** : this output is a single line giving the number of scaffolds that have outlier windows across all samples  

We ran this script utilizing the following required software (other versions of these will probably also work):  
Python version 2.7.12 (Python Software Foundation, 2016)  
matplotlib version 1.5.1 (Hunter, 2007; Matplotlib Development Team, 2016)  
NumPy version 1.10.4 (van der Walt et al., 2011; NumPy Developers, 2016)  
SciPy version 0.17.0 (Jones et al., 2001; van der Walt et al., 2011; SciPy developers 2016)  

## Welch's _t_-test
We ran Welch's _t_-tests to look for significant differences in the spotted owl ancestry of populations.  
  
Input file format is a comma-separated value (CSV) file with the ancestry value in the first field and the population name in the second field:  
```
0.068544,Eastern_Barred_Owl
0.0695239,Eastern_Barred_Owl
0.0683806,Siskiyou_Barred_Owl
```
Usage example:  
  
```
$ python Welch_ttest.py ancestry_value.csv
```
The script produces an output file named "Welch_ttest_out.txt" with output like this:
```
$ cat Welch_ttest_out.txt
Siskiyou vs Western Barred Owls
calculated t-statistic, two-tailed p-value
Ttest_indResult(statistic=-0.77077671527931058, pvalue=0.44863619081421191)
AllWestern vs Eastern Barred Owls
calculated t-statistic, two-tailed p-value
Ttest_indResult(statistic=2.9678999739315834, pvalue=0.035617969491983037)
Pre vs. Post Spotted Owls
calculated t-statistic, two-tailed p-value
Ttest_indResult(statistic=-0.93088180034592116, pvalue=0.52193033083919926)
AllBarred vs AllSpotted Owls
calculated t-statistic, two-tailed p-value
Ttest_indResult(statistic=-1913.4944317390673, pvalue=3.0081979197771044e-43)
```
There are three lines for each _t_-test. The first gives the names of the populations compared. The second gives a description of the values on the third line of the output. The third line gives the _t_-test t-value and p-value.  
  
We ran this script utilizing the following required software (other versions of these will probably also work):  
Python version 2.7.12 (Python Software Foundation, 2016)  
NumPy version 1.10.4 (van der Walt et al., 2011; NumPy Developers, 2016)  
SciPy version 0.17.0 (Jones et al., 2001; van der Walt et al., 2011; SciPy developers 2016)    

## Nucleotide diversity and FST calculation
countFstPi is a C script (provided executable compiled on an Ubuntu system) that takes an input file (of # of reads supporting ref vs. alt allele for each individual, without chr or pos information) and outputs pi and Fst. Calculations choose a random read from each individual at each site. A seed file called "seedms" that provides the initial seeds of the random numbers must be present for the program to work.  
  
Usage example:
```
$ cat <input_file> | ./countFstPi k pop1 pop2
```
k = number of individuals in input file (each line should have 2k fields)
pop1 and pop2 are 2-line files giving the number of individuals in the population on the first line and the sample numbers corresponding to the members of the population  
Example of pop file:  
```
$ cat pop1
6
1 2 3 4 5 6
```
line 1: 6 individuals in pop1  
line 2: individuals 1-6 are in pop1 (note this is not designating the input file columns for pop1 as there are 2 columns per individual)  
  
Output:
piW1 = pi within pop1  
piW2 = pi within pop2  
piB = pi between pop1 and pop2  
Fst = fixation index between pop1 and pop2
  
## Citing the repository

#### Authorship
Code authors: James B. Henderson, jhenderson@calacademy.org; Zachary R. Hanna, zachanna@berkeley.edu; Jeffrey D. Wall  
README.md author: Zachary R. Hanna  

#### Version 1.0.0

Please cite this repository as follows:  

Hanna ZR, Henderson JB, Wall JD. 2017. SPOW-BDOW-introgression-scripts. Version 1.0.0. Zenodo. DOI:  

## References
DePristo MA., Banks E., Poplin R., Garimella KV., Maguire JR., Hartl C., Philippakis AA., del Angel G., Rivas MA., Hanna M., McKenna A., Fennell TJ., Kernytsky AM., Sivachenko AY., Cibulskis K., Gabriel SB., Altshuler D., Daly MJ. 2011. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics 43:491–498. DOI: 10.1038/ng.806.  
  
Fox B., Ramey C. 2013. echo (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Free Software Foundation. 2012. GNU Awk. Version 4.0.1. Available at <https://www.gnu.org/software/gawk/>.  
  
Free Software Foundation 2014. GNU Grep. Version 2.16. Available at <https://www.gnu.org/software/grep/>.  
  
Granlund T., Stallman RM. 2013. cat (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Haertel M., Eggert P. 2013. sort (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Hunter JD. 2007. Matplotlib: A 2D graphics environment. Computing In Science & Engineering. 9:90–95. doi: 10.1109/MCSE.2007.55.  
  
Ihnat DM., MacKenzie D., Meyering J. 2013. cut (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Jones E, Oliphant T, Peterson P, SciPy developers. 2001. SciPy: Open Source Scientific Tools for Python. [Accessed 2017 Sep 29]. Available from: http://www.scipy.org.  
  
MacKenzie D. 2013. fold (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
MacKenzie D., Meyering J. 2013. head (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
McKenna A., Hanna M., Banks E., Sivachenko A., Cibulskis K., Kernytsky A., Garimella K., Altshuler D., Gabriel S., Daly M., DePristo MA. 2010. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research 20:1297–1303. DOI: 10.1101/gr.107524.110.  
  
Matplotlib Development Team. 2016. matplotlib. Version 1.5.1. [Accessed 2017 Sep 29]. Available from: https://github.com/matplotlib/matplotlib. doi: 10.5281/zenodo.44579.  
  
NumPy Developers. 2016. NumPy. Version 1.10.4. [Accessed 2017 Sep 29]. Available from: https://github.com/numpy/numpy.  
  
Python Software Foundation 2016. Python. Version 2.7.12. [Accessed 2016 Oct 1]. Available from: https://www.python.org.  
  
SciPy developers. 2016. SciPy. Version 0.17.0. [Accessed 2017 Sep 29]. Available from: https://github.com/scipy/scipy.  
  
van der Auwera GA., Carneiro MO., Hartl C., Poplin R., del Angel G., Levy-Moonshine A., Jordan T., Shakir K., Roazen D., Thibault J., Banks E., Garimella KV., Altshuler D., Gabriel S., DePristo MA. 2013. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics 11:11.10.1-11.10.33. DOI: 10.1002/0471250953.bi1110s43.  
  
van der Walt S, Colbert SC, Varoquaux G. 2011. The NumPy Array: A Structure for Efficient Numerical Computation. Computing in Science & Engineering. 13:22–30. doi: 10.1109/MCSE.2011.37.  
