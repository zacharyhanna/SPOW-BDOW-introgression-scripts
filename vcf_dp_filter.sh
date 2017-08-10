#!/bin/bash
in_vcf=$1
awk '/^#/ { next }
    {for(i=1; i<=(split($8,zfd,";"));i++) if(split(zfd[i],numDP,"DP=")==2 && numDP[2]<301) {print}}' $in_vcf
