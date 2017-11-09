#!/bin/bash
# modified vcf_qual_filter.sh for pi calculations - only examining biallelic sites with probability >50
# 14Mar2017 incorporated vcf_bl_filter.sh so this should do all the quality filtering in this script
# add removal of any comments so can be used on original vcf and change to use 41 for file
# and allow a missing $1 to let stdin by used instead of hardwiring vcf_snps.vcf
# Use mawk instead of awk if it is available, approx 10% faster
in_vcf=$1
awkpgm="mawk"
[ -z "$(which mawk)" ] && awkpgm="awk"

$awkpgm '/^#/ { next }
     $6 > 50 && $4 ~ "^[ACGT]$" && $5 ~ "^[ACGT]$"' $in_vcf \
| $awkpgm '$1=="C7961234"||$1=="C7963448"||$1=="C7970814"||$1=="C8091874"||$1=="scaffold3674"{next}
    $43 ~ "1/1"{next}
    {for(i=1; i<=(split($8,zfd,";"));i++) if(split(zfd[i],numDP,"DP=")==2 && numDP[2]<301) {print}}'
