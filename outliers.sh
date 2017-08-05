#!/bin/bash
msfile="means_stdevs_ad.txt"
#wndfile="wnd_50k_1000.out"
wndfile=$1

if [ ! -f $msfile ]; then
  >&2 printf "\n\"$msfile\" does not exist. Run\n   ./compute_ad_mean_stdev.sh ad_pct.txt >$msfile\nto create it and then rerun this script.\n\n"
  exit 1
fi

awk 'function min(a,b){ return (a<=b)?a:b } function max(a,b){ return (a>=b)?a:b }
     NR==1{ num=NF;for(i=1;i<=num;i++){split($i,a,",");means[i]=a[1];stdev[i]=a[2];low[i]=max(0,a[1]-a[2]);hi[i]=min(1,a[1]+a[2])} }
     /^#/ { next }
     FNR!=NR{outl=0;line="";
             for(f=5;f<=NF;f++) {
                i=f-4; val= ($f == "-") ? -1 : $f;
                if ( val >= 0 && (val < low[i] || val > hi[i]) ) {line=line " " val; outl++}
                else {line=line " -"}
             }
             if (outl > 0) {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" outl "\t" line}
     }
     # END {for(i=1; i<=num;i++){print i, means[i], stdev[i], low[i], hi[i]}} # this prints them with index and one per line so can do sort -k2,2n
    ' $msfile $wndfile
