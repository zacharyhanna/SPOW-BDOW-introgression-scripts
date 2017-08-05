#!/bin/bash
# 14Mar2017 incorporated vcf_bl_filter.sh so this should do all the quality filtering in this script
# add removal of any comments so can be used on original vcf and change to use 41 for file
# and allow a missing $1 to let stdin by used instead of hardwiring vcf_snps.vcf
# Use mawk instead of awk if it is available, approx 10% faster
in_vcf=$1
awkpgm="mawk"
[ -z "$(which mawk)" ] && awkpgm="awk"

$awkpgm '/^#/ { next }
     $16 ~ "1/1" && $43 ~ "0/0" && $6 > 50 && $4 ~ "^[ACGT]$" && $5 ~ "^[ACGT]$"' $in_vcf \
| $awkpgm 'BEGIN{gqual=".*:.*:.*:[3-9][0-9]*";boAD="./.:0,";soAD="./.:[1-9][0-9]*,0"}
     $1=="C7961234"||$1=="C7963448"||$1=="C7970814"||$1=="C8091874"||$1=="scaffold3674"{next}
     {if(match($16,gqual) && match($43,gqual) && match($16,boAD) && match($43, soAD)){print}}'
