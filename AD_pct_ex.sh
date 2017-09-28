#!/bin/bash
# 2017Se21 change to output sum of reads at a snp after allele depth percentage
# new output format is adp:sum
# change to see which field AD is in to accommodate various VCF formats
# presumption is that AD is in same place everywhere as it is in record 1
#
# also original version had an error for this data "./.:.:1" e.g from scaffold12 snp 756
# where using [3] and [4] had it reporting 1 instead of -1
#
# first awk will inject a single number as first line of input to second awk and for other lines
# will perform the equivalent of cut -f1,2,10-
#
# reason we change technique is so we use stdin to pipe in the file

#ad_pos=$(head -1 $1|awk '{if(split($9,a,":")>1){for(i=1;i<=length(a);i++){if(a[i]=="AD")ad_ix=i}}else if($9=="AD")ad_ix=1;}END{print ad_ix}')
#cut -f1,2,10- $1 | \
vcf_input=$1
awk 'BEGIN{OFS="\t"; OFMT="%.2g";CONVFMT=OFMT}
     NR==1{if(split($9,a,":")>1){for(i=1;i<=length(a);i++){if(a[i]=="AD")ad_ix=i}}else if($9=="AD")ad_ix=1; print ad_ix}
     {printf "%s\t%s", $1,$2; for(f=10;f<=NF;f++)printf "\t%s", $f; printf "\n" }' $vcf_input | \
  awk -v AD_pos=$ad_pos \
       'BEGIN{OFMT="%.2g";CONVFMT=OFMT}
        NR==1{if(NF==1 && $1 > 0){AD_pos=$1};
              if(AD_pos){m="AD in subfield "AD_pos" of field 9"}else{m="AD subfield not found in field 9"}; print m > "/dev/stderr";
              if(!AD_pos)exit; if(NF==1)next;
       }
       { printf "%s %s ", $1, $2;
        for(f=3; f<=NF; f++) {
           adp = -1;
           if(split($f,a,":") >= AD_pos) {
              if (split(a[AD_pos],ad,",")==2) {
                sum=ad[1]+ad[2]; if(sum!=0) adp=ad[1]/sum
              }
           } else sum = 0
           printf "%s:%s ", adp,sum
         }
         printf "\n"
       }'

exit 0
# original script below (advantage of only 1 split but depends on AD being in 2nd field with 2 parm 1st field)
cut -f1,2,10- $1 | \
  awk '{printf "%s %s ", $1, $2;
        for(f=3; f<=NF; f++) {
           adp = -1;
           if(split($f,a,/[,:\/]/) > 1) {
              sm=a[3]+a[4]; if(sm!=0) adp=a[3]/sm
           }
           printf "%s ", adp
         }
         printf "\n"
        }'
