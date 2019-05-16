#!/bin/bash
cstart=0
cstep=8
cend=1240

ts_rng_seed=13713032
Nt=64
nts=8
dts=$(( $Nt / $nts ))

arch=pascal

echo "cid src_ts" > src_timeslices.dat

for cid in $(seq ${cstart} ${cstep} ${cend}); do
  echo $cid

  c4id=$(printf %04d ${cid})

  if [ ! -d ${c4id} ]; then
    mkdir ${c4id}
  fi

  yaml_input=${c4id}/definitions.yaml
  tmlqcd_infile=${c4id}/invert.input
  cvc_infile=${c4id}/cpff.input
  jscript=${c4id}/job.cmd

  cp templates/definitions.yaml.template ${yaml_input}

  cp templates/invert.input.${arch}.template ${tmlqcd_infile}
  sed -i "s/_NSTORE_/${cid}/g" ${tmlqcd_infile}

  cp templates/job_${arch}.cmd.template ${jscript}
  sed -i "s/_NSTORE_/${c4id}/g" ${jscript}

  cp templates/cpff.input.template ${cvc_infile}

  # use rng.py to genereate a starting time slice and use maximal
  # separation in time to generate $nts time slices from this
  if [ ! -f rng.py ]; then
    echo rng.py could not be found!
    exit 1
  fi
  ts1=$(./rng.py --seed ${ts_rng_seed} --cid ${cid} --Nt ${Nt} )
  src_ts=( $ts1 )
  for i in $(seq 1 $(( $nts - 1 )) ); do
    new_ts=$(( (${ts1} + ${i}*${dts}) % ${Nt} ))  
    src_ts+=( $new_ts )
  done
  for i in $(seq 0 $(( $nts -1 )) ); do
    echo $cid ${src_ts[${i}]} >> src_timeslices.dat
    echo "source_coords = ${src_ts[${i}]},0,0,0" >> ${cvc_infile}
  done
done

