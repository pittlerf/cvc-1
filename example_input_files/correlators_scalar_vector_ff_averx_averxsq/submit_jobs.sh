#!/bin/bash
cstart=0
cstep=8
cend=1240

for cid in $(seq ${cstart} ${cstep} ${cend}); do
  c4id=$(printf %04d $cid)
  pushd ${c4id}
  stdout=$(sbatch job.cmd)
  echo $stdout
  jobid=$(echo $stdout | awk '{print $4}')
  #scontrol update job $jobid priority=65000
  popd
done
