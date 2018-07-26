#!/bin/bash
# This is a version of the sliding_window.sh script that writes number of SNPs for each individual in
#    a column after the window val for the individual is written. This number has an s after it to make it easier to distinguish
# Also so this uses tabs between each field to make it easier to Paste into a spreadsheet
# And format %.2g used for some rounding of val output

# the input file should have ScaffName SnpPos adpct_ind1 adpct_ind2...
# e.g.: scaffold12 35061 0 -1 -1 -1 -1 0 0 0 -1 -1 -1 -1 -1 -1 -1 -1 1 -1 -1 0 0 0 0 0 0.625 0 0 -1 -1 -1 0 -1 1 1 -1 -1 -1 -1 0 0 -1 0 0 -1 -1 1 -1 1 1 -1 0 0 -1 0.25
# where -1 is placeholder for no value on input, on output if no sums for an ind then output is "-"
# this version v4 calculates the sums for each sliding window as does v3 but for the first record of
# a scaffold shows results for the entire snp set of the scaffold and shows first and last snp loc in fld 2 and 3, resp.
# the first line od a scaffold that summarizes its values for all snps begins with a # char so you
# can easily grep "^#" or grep -v "^#" to get the summary lines or exclude the summary lines

# accomodate usage where args given, a filename, one or two numeric args with or without the file
# arguments can be in any order. first 2 numeric args are used and last non-numeric arg used for file

# changed to file argument required or if missing uses stdin
# also changed num_slides_wnd to 1 so just 1 numeric arg defines a non-overlapping window

# changed to let first argument be window size and second argument be slide size, but checks to
# make sure that the window size is a multiple of the slide size (we require this since
# it lets use use a faster algorithm and this property seems to be met in all papers I read)

file=""  # was ad_pct.txt
wnd_size=50000
slide_size=0 # 0 means use window size

num_ints=0
numre="^[0-9]+$"
for arg in "$@" ; do
   if [[ $arg =~ $numre ]] ; then
      num_ints=$((num_ints+1))
      [[ ( $num_ints == 1 ) ]] && wnd_size=$arg
      [[ ( $num_ints == 2 ) ]] && slide_size=$arg
   else
      file=$arg
   fi
done

(($slide_size == 0)) && slide_size=$wnd_size
mod=$(($wnd_size % $slide_size))
if (( mod > 0 )); then
   >&2 echo "Slide of $slide_size must be a multiple of window size $wnd_size" 
   exit 1
fi
num_slides_in_wnd=$(($wnd_size / $slide_size))

awk -v SLIDE_SIZE=$slide_size -v SLIDES_IN_WND=$num_slides_in_wnd 'BEGIN {WndSize = SLIDE_SIZE * SLIDES_IN_WND;
            m1 = SLIDE_SIZE-1; wnd_inc = WndSize-1; wnd_slides_inc = SLIDES_IN_WND-1;
            print "Window of " WndSize " bp sliding " SLIDE_SIZE " bp at a time."  > "/dev/stderr"
     }
     function calc_all_scaffold_snps() {
        delete adpct_sum; delete snps_used;
        for (sld_ix = 1; sld_ix <= max_sld_ix; sld_ix++) {
            if (slide_snps[sld_ix] > 0) { # there are snps in this slide portion of the window
               snps_in_wnd += slide_snps[sld_ix]

               for (ind=1; ind <= NUM_IND; ind++) {  # we have from 1 to 54 individuals to calculate values
                  if (aSlideSnps[sld_ix][ind] > 0) { # this indiviudal had one or more snps in this slide slice
                     snps_used[ind] += aSlideSnps[sld_ix][ind]
                     adpct_sum[ind] += aSlideSums[sld_ix][ind] # this is the sum for individual ind in the current slide slice of the window
                  }
               }
           }
        }

        printf "#%s\t%s\t%s\t%s\t", cur_scaff, first_snp, last_snp, snps_in_wnd
        for (ind=1; ind <= NUM_IND; ind++) {
           if (snps_used[ind]==0) val = "-"; else {val = adpct_sum[ind] / snps_used[ind]}
           printf "%.2g:%d:%d\t", val, snps_used[ind], read_sum[ind];
        }
        printf "\n"
     }

     function calc_one_window_by_sld_inf(sldsum_ary, sldsnpcount_ary, sldreadcount_ary, wnd_start_slide, wnd_end_slide) {
        w_strt = ((wnd_start_slide - 1) * SLIDE_SIZE) + 1
        w_end =  (wnd_end_slide < max_sld_ix) ? w_strt + WndSize - 1 : last_snp #print last_snp pos if window ends there

        delete adpct_sum; delete snps_used; delete read_sum;
        snps_in_wnd = 0;
        for (sld_ix = wnd_start_slide; sld_ix <= wnd_end_slide; sld_ix++) {
            if (slide_snps[sld_ix] > 0) { # there are snps in this slide portion of the window
               snps_in_wnd += slide_snps[sld_ix]

               for (ind=1; ind <= NUM_IND; ind++) {  # we have from 1 to 54 individuals to calculate values
                  if (aSlideSnps[sld_ix][ind] > 0) { # this indiviudal had one or more snps in this slide slice
                     snps_used[ind] += aSlideSnps[sld_ix][ind]
                     adpct_sum[ind] += aSlideSums[sld_ix][ind] # this is the sum for individual ind in the current slide slice of the window
                     read_sum[ind] += aSlideReads[sld_ix][ind] # this is the sum for individual reads in the current slide slice of the window
                  }
               }
           }
        }
        if (snps_in_wnd == 0)
           return #nothing in this window to report

        printf "%s\t%s\t%s\t%s\t", cur_scaff, w_strt, w_end, snps_in_wnd
        for (ind=1; ind <= NUM_IND; ind++) {
           if (snps_used[ind]==0) val = "-"; else {val = adpct_sum[ind] / snps_used[ind]}
           printf "%.2g:%d:%d\t", val, snps_used[ind], read_sum[ind]; #output form is ancestry_value:snps_for_window:reads_for_window
        }
        printf "\n"
     }

     function calc_scaff_windows(sldsum_ary, sldsnpcount_ary, sldreadcount_ary, last_slide_ix) {
        calc_all_scaffold_snps()

        start_wnd = 1; end_wnd = start_wnd + wnd_inc

        for (;;start_wnd += SLIDE_SIZE ) {
           wnd_start_slide = int((start_wnd + m1)/SLIDE_SIZE)
           wnd_end_slide   = wnd_start_slide + wnd_slides_inc
           # if (wnd_start_slide > last_slide_ix) break; # done with this scaffold. using this test allows shrinking final windows
           if (wnd_end_slide > last_slide_ix && wnd_start_slide>1) break; # done with this scaffold. using this test stops after last snp is at end of window

           calc_one_window_by_sld_inf(sldsum_ary, sldsnpcount_ary, sldreadcount_ary, wnd_start_slide, wnd_end_slide)
        }
     }

     function new_scaff() {cur_scaff=$1; scaff_count++; first_snp=$2; last_snp=0; max_sld_ix=0; delete aSlideSums; delete aSlideSnps; delete aSlideReads; delete slide_snps}

     BEGIN{OFMT="%.2g"}
     NR==1 {NUM_IND = NF -2; fst_ind = 3; new_scaff()}
     cur_scaff != $1 { calc_scaff_windows(aSlideSums, aSlideSnps, aSlideReads, max_sld_ix); new_scaff() }
     # scaff_count > 2 { exit }
     { sld_ix = int(($2 + m1)/SLIDE_SIZE); max_sld_ix=sld_ix; last_snp=$2; slide_snps[sld_ix]++;
       fnum=1;for(f=fst_ind;f<=NF;f++){
          ad=$f+0; if(ad >= 0){aSlideSums[sld_ix][fnum] += ad; aSlideSnps[sld_ix][fnum] += 1}
          if(split($f,ards,":")==2){aSlideReads[sld_ix][fnum] += ards[2]}
          fnum++
       }
     }
     END{calc_scaff_windows(aSlideSums, aSlideSnps, aSlideReads, max_sld_ix)}
' $file
