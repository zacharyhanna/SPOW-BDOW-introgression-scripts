#!/bin/bash
window_file=$1
awk '
BEGIN{OFS="\t"}
{for (a=1; a<=1; a++){
  printf "%s", $a
  };
for (b=2; b<=4; b++){
  printf "\t%s", $b
  };
for (c=5; c<=NF; c+=2){
  compare = $(c+1)
  sub(/s/,"",compare)
  if (compare+0 > 9) {
    printf "\t%s", $c
  }
  else {printf "\t-1"}
};
{printf "\n"}
}' $window_file
