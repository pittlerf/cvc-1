#!/bin/bash
while IFS= read -r c4d
do
  pushd $c4d
  rm *.h5
  sbatch job.cmd
  popd
done < incomplete.txt
