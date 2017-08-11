#!/bin/bash
in_vcf=$1
awk '{for(i=1; i<=(split($8,zfd,";"));i++) if(split(zfd[i],numDP,"DP=")==2) {sumDP+=numDP[2]; sqrsumDP += (numDP[2])^2; sites++}}
    END{div=sites;
        if(div > 0) {
            mean = sumDP/div;
            stdev = sqrt((sqrsumDP-sumDP^2/div)/div)
        }
        printf "meanDP = %s,stdevDP = %s,number of sites = %s", mean, stdev, div
        printf "\n"
      }' $in_vcf
