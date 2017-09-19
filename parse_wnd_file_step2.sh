#!/bin/bash
in_file=$1
awk '
BEGIN{OFS="\t"}
{
  numfld=NF
  for (i=5; i<=NF; i++) {
    if ($i != -1) {
      {cnt[i] += 1}
    }
  }
}
END{{if (cnt[5] >= 1) {
      printf "%s", cnt[5]}
    else {printf "0"}
    }
    {for (d=6; d<=numfld; d++) {
        if (cnt[d] >= 1) {
            printf "\t%s", cnt[d]
        }
      else {printf "\t0"}
      }
    }
    {printf "\n"}
}' $in_file
