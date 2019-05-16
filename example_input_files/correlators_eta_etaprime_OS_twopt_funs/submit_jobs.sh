cstart=0
cstep=8
cend=1240

for cid in $(seq ${cstart} ${cstep} ${cend}); do
  c4id=$(printf %04d $cid)
  pushd ${c4id}
  sbatch job.cmd
  popd
done
