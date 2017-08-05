# SPOW-BADO-introgression-scripts
We here provide the scripts that we developed for analyzing introgression in whole genome sequences obtained from Barred and Spotted Owls. 

## Introduction
These scripts take a variant call format (vcf) file as input. We created our vcf file using UnifiedGenotyper from the Genome Analysis Toolkit (GATK) version 3.4-46 (DePristo et al. 2011; McKenna et al. 2010; Van der Auwera et al. 2013). Due to the specifics of our sample set and our analyses, you will not be able to directly use most of these scripts without making some modifications to the code. VCF files produced by other variant callers may require modifications to these scripts for correct parsing of the variant files. We are providing our code here for the purposes of documentation and with the hope that some of our methods may prove useful to others in their own genome-scale analyses.  
#### Specifics to our analyses that have affected the code
VCF file output from UnifiedGenotyper: raw_variants.vcf
We set hard filters filtered Our reference individual sample

## Sliding Window Analyses Pipeline
Pipe from original vcf file to the ad pct file and use this for the sliding window.  
You could pipe all the way through, but, since multiple runs of the sliding window routine are likely, we saved the ad_pct.txt file and worked from it.  
  
#### Example of pipeline usage (40,000 bp windows with no overlap):  
$ ./vcf_qual_filter.sh raw_variants.vcf | ./AD_pct.sh >ad_pct.txt  
$ cat ad_pct.txt | ./sliding_window.sh 40000 >wnd_40k_noovlp.txt  

#### can run additional sliding windows, e.g., 40,000 base window sliding 5,000 bases at a time:  
$ ./sliding_window.sh ad_pct.txt 40000 5000 >wnd_40k_5k_slide.txt  

### to compute mean and std dev on AD pct file and keep this info for further use do this
$ compute_ad_mean_stdev.sh ad_pct.txt >means_stdevs_ad.txt  

### to check for outliers in the window do this:
$ outliers.sh wnd_40k_noovlp.txt  

### to get column names
$ grep -v "^#" raw_variants.vcf -B1 | head -1  
### to merge column names with means and stddev
$ awk 'NR==1{b=1;for(i=10;i<=NF;i++)nm[b++]=$i}
     NR>1{for(i=1;i<=NF;i++){split($i,ms,",");
          mminus = (ms[1]-ms[2] > 0) ? ms[1]-ms[2] : 0;
          mplus  = (ms[1]+ms[2] > 1) ? 1 : ms[1]+ms[2]
          print i, nm[i], ms[1],ms[2],mminus,mplus}}' <(grep -v "^#" 2016Feb18_2maskedaligned_UnifGen_raw_variants.vcf -B1 |head -1) means_stdevs_ad.txt \
 | sort -k3,3n | awk '{printf"%2s\t%12s\t",$1,$2;for(i=3;i<=NF;i++)printf "%-9s\t",$i;print ""}' \
 | cat <(echo -e "Sample#\t Sample_Name\tmean\t\tstdev\t\tmean-stdev\tmean+stdev") -

### Citing the repository

#### Authorship
Code authors: Zachary R. Hanna, zachanna@berkeley.edu; James B. Henderson, jhenderson@calacademy.org; Jeffrey D. Wall  
README.md author: Zachary R. Hanna  

#### Version 1.0.0

Please cite this repository as follows:  

Hanna ZR, Henderson JB, Wall JD. (2017). SPOW-BADO-introgression-scripts. Version 1.0.0. Zenodo. DOI:  

### References
DePristo MA., Banks E., Poplin R., Garimella KV., Maguire JR., Hartl C., Philippakis AA., del Angel G., Rivas MA., Hanna M., McKenna A., Fennell TJ., Kernytsky AM., Sivachenko AY., Cibulskis K., Gabriel SB., Altshuler D., Daly MJ. 2011. A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics 43:491–498. DOI: 10.1038/ng.806.  
  
McKenna A., Hanna M., Banks E., Sivachenko A., Cibulskis K., Kernytsky A., Garimella K., Altshuler D., Gabriel S., Daly M., DePristo MA. 2010. The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research 20:1297–1303. DOI: 10.1101/gr.107524.110.  
  
Van der Auwera GA., Carneiro MO., Hartl C., Poplin R., del Angel G., Levy-Moonshine A., Jordan T., Shakir K., Roazen D., Thibault J., Banks E., Garimella KV., Altshuler D., Gabriel S., DePristo MA. 2013. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. Current Protocols in Bioinformatics 11:11.10.1-11.10.33. DOI: 10.1002/0471250953.bi1110s43.  
