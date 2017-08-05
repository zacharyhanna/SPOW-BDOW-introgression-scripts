# SPOW-BADO-introgression-scripts

## Sliding Window Analyses Pipeline
Pipe from original vcf file to the ad pct file and use this for the sliding window.  
You could pipe all the way through, but, since multiple runs of the sliding window routine are likely, we saved the ad_pct.txt file and worked from it.  
  
Example of pipeline usage:  
$ ./vcf_qual_filter.sh 2016Feb18_2maskedaligned_UnifGen_raw_variants.vcf | ./AD_pct.sh >ad_pct.txt  
$ cat ad_pct.txt | ./sliding_window.sh 40000 >wnd_40k_noovlp.txt  

### can run additional sliding windows, e.g., 40,000 base window sliding 5,000 bases at a time:  
$ ./sliding_window.sh ad_pct.txt 40000 5000 >wnd_40k_5k_slide.txt  

### to compute mean and std dev on AD pct file and keep this info for further use do this
$ compute_ad_mean_stdev.sh ad_pct.txt >means_stdevs_ad.txt  

### to check for outliers in the window do this:
$ outliers.sh wnd_40k_noovlp.txt  

### to get column names
$ grep -v "^#" 2016Feb18_2maskedaligned_UnifGen_raw_variants.vcf -B1 | head -1  
### to merge column names with means and stddev
$ awk 'NR==1{b=1;for(i=10;i<=NF;i++)nm[b++]=$i}
     NR>1{for(i=1;i<=NF;i++){split($i,ms,",");
          mminus = (ms[1]-ms[2] > 0) ? ms[1]-ms[2] : 0;
          mplus  = (ms[1]+ms[2] > 1) ? 1 : ms[1]+ms[2]
          print i, nm[i], ms[1],ms[2],mminus,mplus}}' <(grep -v "^#" 2016Feb18_2maskedaligned_UnifGen_raw_variants.vcf -B1 |head -1) means_stdevs_ad.txt \
 | sort -k3,3n | awk '{printf"%2s\t%12s\t",$1,$2;for(i=3;i<=NF;i++)printf "%-9s\t",$i;print ""}' \
 | cat <(echo -e "Sample#\t Sample_Name\tmean\t\tstdev\t\tmean-stdev\tmean+stdev") -
