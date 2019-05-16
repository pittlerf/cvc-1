#!/bin/bash

nts=8

conf_start=0
conf_step=8
conf_end=1240

conf_missing=""
conf_incomplete=""

for cid in $(seq $conf_start $conf_step $conf_end); do
  c4d=$(printf %04d $cid)
  echo -n $c4d
  files="$(ls $c4d/*.h5)"
  if [ $? -ne 0 ]; then
    conf_missing+=( $c4d )
    echo " missing"
  else
    nh5=$( echo $files | wc | awk '{print $2}' )
    echo "  $nh5"
    if [ $nh5 -ne $nts ]; then
      conf_incomplete+=( $c4d )
    fi
  fi
done

if [ -f missing.txt ]; then
  rm missing.txt
fi

if [ -f incomplete.txt ]; then
  rm incomplete.txt
fi

for cid in ${conf_missing[@]}; do
  echo $cid >> missing.txt
done

for cid in ${conf_incomplete[@]}; do
  echo $cid >> incomplete.txt
done


