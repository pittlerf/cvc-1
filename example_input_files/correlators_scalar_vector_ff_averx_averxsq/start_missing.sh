#!/bin/bash
while IFS= read -r c4d
do
  pushd $c4d
  sbatch job.cmd
  popd
done < missing.txt
