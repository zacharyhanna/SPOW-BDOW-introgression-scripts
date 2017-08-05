#!/bin/bash
awk '   { maxi=NF-2; for(f=3; f<=NF; f++){i=f-2;if($f > -1){sum[i] += $f; sqrsum[i] += ($f)^2; snps[i]++}} }
     END{ for (i=1; i<=maxi; i++){
             div=snps[i]; mean = -1; stdev = -1;
             if (div > 0) {
                mean = sum[i]/div;
                stdev = sqrt((sqrsum[i]-sum[i]^2/div)/div);
             }
             printf "%s,%s ", mean, stdev
          }
          printf "\n" }' $1
