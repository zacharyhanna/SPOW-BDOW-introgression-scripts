#!/bin/bash
vcf_input=$1
awk 'BEGIN{OFS="\t"}
     NR==1{if(split($9,a,":")>1){for(i=1;i<=length(a);i++){if(a[i]=="DP")dp_ix=i}}else if($9=="DP")dp_ix=1; print dp_ix}
     {printf "%s\t%s", $1,$2; for(f=10;f<=NF;f++)printf "\t%s", $f; printf "\n" }' $vcf_input | \
awk 'NR==1{if(NF==1 && $1 > 0) {dp_pos=$1};
            if(dp_pos){m="DP in subfield "dp_pos" of field 9"}else{m="DP subfield not found in field 9"};
              print m > "/dev/stderr";
            if(!dp_pos)exit; if(NF==1)next;
          }
      {maxi=NF-2; snps++;
        for(f=3; f<=NF; f++) {
          if(split($f,a,":") >= dp_pos) {
            sum[f] += a[dp_pos]; sqrsum[f] += (a[dp_pos])^2;
          }
        }
      }
      END{for (i=3; i<=maxi; i++){
            div=snps; mean = -1; stdev = -1;
            if(div > 0){
              mean = sum[i]/div;
              stdev = sqrt((sqrsum[i]-sum[i]^2/div)/div);
            }
            printf "%s,%s ", mean, stdev
          }
            printf "\n"
      }'
