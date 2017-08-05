# SPOW-BADO-introgression-scripts

## Sliding Window Analyses Pipeline
Pipe from original vcf file to the ad pct file and use this for the sliding window.  
You could pipe all the way through, but, since multiple runs of the sliding window routine are likely, we saved the ad_pct.txt file and worked from it.  
  
Example of pipeline usage:  
$ ./vcf_qual_filter.sh 2016Feb18_2maskedaligned_UnifGen_raw_variants.vcf | ./AD_pct.sh >ad_pct.txt  
cat ad_pct.txt | ./sliding_window.sh 40000 >wnd_40k_noovlp.txt  
