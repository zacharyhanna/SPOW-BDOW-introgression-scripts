#!/bin/bash
geno_file=$1
window_file=$2
cat $geno_file $window_file |\
awk '
BEGIN{OFS="\t"}
{if(NR==1) {
        last_geno = NF+4;
        for (i=1; i<=NF; i++) {
          geno[i+4]=$i;
          if (geno[i+4]==0) {cutoff[i+4]=0.4}
          else if (geno[i+4]==1) {cutoff[i+4]=0.6}
        }
      }
    else if ($1 !~ /^#/) {
        printf "%s", $1
        for (d=2; d<=4; d++) {
          printf "\t%s", $d
        };
        for (c=5; c<=NF; c++) {
          if ($c != "-" && geno[c] == 0) {
            if ($c > 0.4) {
              printf "\t%s", $c
            }
            else {printf "\t-1"}
          }
          else if ($c != "-" && geno[c] == 1) {
            if ($c < 0.6) {
              printf "\t%s", $c
            }
            else {printf "\t-1"}
          }
        else {printf "\t-1"}
        }
        {printf "\n"}
      }
}
'
