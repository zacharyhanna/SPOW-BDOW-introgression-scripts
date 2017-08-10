# SPOW-BADO-introgression-scripts
We here provide the scripts that we developed for analyzing introgression in whole genome sequences obtained from Barred and Spotted Owls. 

## Contents
* [Introduction](#introduction)  
* [Sliding Window Analyses Pipeline](#sliding-window-analyses-pipeline)  
  * [1. Filter raw vcf file](#1-filter-raw-vcf-file)
  * [2. Site coverage calculation](#2-site-coverage-calculation)
  * [3. Exclude sites with excessive coverage](#3-exclude-sites-with-excessive-coverage)
* [Citing the Repository](#citing-the-repository)  
* [References](#references)  

## Introduction
These scripts take a variant call format (vcf) file as input. We created our vcf file using UnifiedGenotyper from the Genome Analysis Toolkit (GATK) version 3.4-46 (DePristo et al. 2011; McKenna et al. 2010; Van der Auwera et al. 2013). Due to the specifics of our sample set and our analyses, you will not be able to directly use most of these scripts without making some modifications to the code. VCF files produced by other variant callers may require modifications to these scripts for correct parsing of the variant files. We are providing our code here for the purposes of documentation and with the hope that some of our methods may prove useful to others in their own genome-scale analyses.  
#### Specifics to our analyses that have affected the code
- VCF file output from UnifiedGenotyper: raw_variants.vcf  
- We set hard filters filtered Our reference individual sample.
- Fields 16 and 43 in the vcf corresponded with the reference Barred and Spotted Owl samples, respectively.

## Sliding Window Analyses Pipeline

### 1. Filter raw vcf file
Usage example:  
$ ./vcf_qual_filter.sh raw_variants.vcf > filtered_variants.vcf  
  
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
$ cat filtered_variants.vcf | ./dp_cov_script.sh 
  
Example output:  
```
meanDP = 129.19,stdevDP = 34.2783,number of sites = 5821431  
```
### 3) Exclude sites with excessive coverage
We excluded sites with coverage in excess of the mean + 5σ (>300 nt in our data set), as suggested by the GATK documentation (https://software.broadinstitute.org/gatk/guide/article?id=3225).  
  
Usage example:  
$ cat filtered_variants.vcf | ./vcf_dp_filter.sh >filtered_variants2.vcf  

### 4) Calculation of coverage depth per sample
This step is not actually required in the pipeline, but we calculated the mean and standard deviation (σ) of the coverage depth for each sample across all of the sites in our final filtered set of SNPs.  
  
Usage example:  
$ cat filtered_variants2.vcf | ./DP_means_std_dev.sh | head -1 >filtered_variants2_dp_means_stdev.txt  

### 4) Allele depth calculation
Usage example:  
$ cat filtered_variants2.vcf | ./AD_pct.sh >ad_pct.txt  

"AD_pct.sh" shell script requirements:  
GNU Awk - we used GNU Awk version 4.0.1 (Free Software Foundation, 2012)  
cut (GNU coreutils) - we used cut (GNU coreutils) version 8.21 (Ihnat et al. 2013)

### 5) Sliding window calculation  
Usage example:  
$ cat ad_pct.txt | ./sliding_window.sh 40000 >wnd_40k_noovlp.txt  
The above example calculates 40,000 bp windows with no overlap.  
  
"sliding_window.sh" shell script requirements:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  

#### Example of running full sliding window pipeline
As above, this example calculates 40,000 bp windows with no overlap.
  
$ ./vcf_qual_filter.sh raw_variants.vcf | ./AD_pct.sh >ad_pct.txt  
$ cat ad_pct.txt | ./sliding_window.sh 40000 >wnd_40k_noovlp.txt  

The first pipe starts with the original vcf file and outputs to the ad_pct.txt file, which is then used in the second pipe for the sliding window analysis.  
These scripts enable you to pipe all the way through, but, since multiple runs of the sliding window routine are likely, we saved the ad_pct.txt file and worked from it for multiple sliding window analyses.  

##### Example of a running a different sliding window
This example calculates 40,000 base windows sliding 5,000 bases at a time.  
  
$ ./sliding_window.sh ad_pct.txt 40000 5000 >wnd_40k_5k_slide.txt  

### Compute means and standard deviations on allele depth file (AD_pct.txt)
We are keeping the output for further use in the file "means_stdevs_ad.txt".  
  
Usage example:  
$ compute_ad_mean_stdev.sh ad_pct.txt >means_stdevs_ad.txt    
  
The shell script requires:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  

### Check for outliers in the window
Usage example:  
$ outliers.sh wnd_40k_noovlp.txt  

### Get column (sample) names from vcf
Usage example:  
$ grep -v "^#" raw_variants.vcf -B1 | head -1  
  
The command requires:  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  
head (GNU coreutils) - we used head (GNU coreutils) version 8.21 (Ihnat et al. 2013)  

### Merge column (sample) names from vcf with means and standard deviations
Usage example:  
$ awk 'NR==1{b=1;for(i=10;i<=NF;i++)nm[b++]=$i}
     NR>1{for(i=1;i<=NF;i++){split($i,ms,",");
          mminus = (ms[1]-ms[2] > 0) ? ms[1]-ms[2] : 0;
          mplus  = (ms[1]+ms[2] > 1) ? 1 : ms[1]+ms[2]
          print i, nm[i], ms[1],ms[2],mminus,mplus}}' <(grep -v "^#" raw_variants.vcf -B1 |head -1) means_stdevs_ad.txt | sort -k3,3n | awk '{printf"%2s\t%12s\t",$1,$2;for(i=3;i<=NF;i++)printf "%-9s\t",$i;print ""}' | cat <(echo -e "Sample#\t Sample_Name\tmean\t\tstdev\t\tmean-stdev\tmean+stdev") -
  
The command requires:  
cat (GNU coreutils) - we used cat (GNU coreutils) version 8.21 (Granlund & Stallman 2013)  
echo (GNU coreutils) - we used echo (GNU coreutils) version 8.21 (Fox & Ramey 2013)  
GNU Awk - we used GNU Awk version 4.0.1 (Free Software Foundation, 2012)  
GNU Grep - we used GNU Grep version 2.16 (Free Software Foundation, 2014)  
head (GNU coreutils) - we used head (GNU coreutils) version 8.21 (Ihnat et al. 2013)  
sort (GNU coreutils) - we used sort (GNU coreutils) version 8.21 (Haertel & Eggert 2013)    

## Citing the repository

#### Authorship
Code authors: James B. Henderson, jhenderson@calacademy.org; Zachary R. Hanna, zachanna@berkeley.edu; Jeffrey D. Wall  
README.md author: Zachary R. Hanna  

#### Version 1.0.0

Please cite this repository as follows:  

Hanna ZR, Henderson JB, Wall JD. 2017. SPOW-BADO-introgression-scripts. Version 1.0.0. Zenodo. DOI:  

## References
DePristo MA., Banks E., Poplin R., Garimella KV., Maguire JR., Hartl C., Philippakis AA., del Angel G., Rivas MA., Hanna M., McKenna A., Fennell TJ., Kernytsky AM., Sivachenko AY., Cibulskis K., Gabriel SB., Altshuler D., Daly MJ. 2011. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics 43:491–498. DOI: 10.1038/ng.806.  
  
Fox B., Ramey C. 2013. echo (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Free Software Foundation. 2012. GNU Awk. Version 4.0.1. Available at <https://www.gnu.org/software/gawk/>.  
  
Free Software Foundation 2014. GNU Grep. Version 2.16. Available at <https://www.gnu.org/software/grep/>.  
  
Granlund T., Stallman RM. 2013. cat (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Haertel M., Eggert P. 2013. sort (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
Ihnat DM., MacKenzie D., Meyering J. 2013. cut (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
MacKenzie D. 2013. fold (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
MacKenzie D., Meyering J. 2013. head (GNU coreutils). Version 8.21. Available at <http://www.gnu.org/software/coreutils/coreutils.html>.  
  
McKenna A., Hanna M., Banks E., Sivachenko A., Cibulskis K., Kernytsky A., Garimella K., Altshuler D., Gabriel S., Daly M., DePristo MA. 2010. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research 20:1297–1303. DOI: 10.1101/gr.107524.110.  
  
Van der Auwera GA., Carneiro MO., Hartl C., Poplin R., del Angel G., Levy-Moonshine A., Jordan T., Shakir K., Roazen D., Thibault J., Banks E., Garimella KV., Altshuler D., Gabriel S., DePristo MA. 2013. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics 11:11.10.1-11.10.33. DOI: 10.1002/0471250953.bi1110s43.  
