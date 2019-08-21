#!/bin/bash

CORRBIN=""
if [ -z "$1" ]; then
  echo "The first argument of run_integration_test.sh must be the path to the 'correlators' executable!"
else
  CORRBIN="$1/correlators"
fi

LIMEDIR=""
if [ -z "$2" ]; then
  echo "The second argument of run_integration_test.sh must be the path of the c-lime installation directory!" 
else
  LIMEDIR="$2"
fi

if [ -z "$CORRBIN" -o -z "$LIMEDIR" ]; then
  exit 1
fi

# remove any existing h5 files, just in case (the code does not overwrite!)
rm *.h5

mpiproc=$(nproc)
# since we have parallelisation set up only in T, we can't have more than 8 procs
# in order to not cause to much load, let's limit ourselves to 4
if [ $mpiproc -gt 4 ]; then
  mpiproc=4
fi
mpirun -np ${mpiproc} $CORRBIN

corrfiles=( $(ls *.h5) )
if [ ${#corrfiles[@]} -eq 0 ]; then
  echo "Correlator files were not produced!"
  exit 1
fi

# the absolute comparison is somewhat meaningless as correlators may drop down to very small
# values
#corr_diff_abs_thresh="1e-12"
#echo "Analysing differences in correlator data (absolute deviation threshold: $corr_diff_abs_thresh)"

# The relative comparison is more meaningful, and a relative deviation no larger than the seventh
# decimal place should be more than stringent enough for most purposes!
corr_diff_rel_thresh="1e-7"
echo "Analysing differences in correlator data (relative deviation threshold: $corr_diff_rel_thresh)"

corr_differences=( )
for corrfile in ${corrfiles[@]}; do
  echo Comparing: $corrfile "<->" reference/$corrfile
  # the absolute comparison is not really very meaningful
  #h5diff -d $corr_diff_abs_thresh $corrfile reference/$corrfile
  h5diff -p $corr_diff_rel_thresh $corrfile reference/$corrfile
  if [ $? -ne 0 ]; then
    corr_differences+=( $corrfile )
  fi
done
if [ ${#corr_differences[@]} -ne 0 ]; then
  echo "Correlator differences found in:"
  for i in $(seq 0 $(( ${#corr_differences[@]} - 1 )) ); do
    echo ${corr_differences[$i]}
  done
else
  echo "No differenes found in correlator data!"
fi
echo ------------------------------------------------
echo

echo "Extracting all propagator data"
./extract_txtprop.sh
echo ------------------------------------------------
echo

txtprop_files=( $(ls *.txtprop) )
if [ ${#txtprop_files[@]} -eq 0 ]; then
  echo "txtprop files were not produced!"
  exit 1
fi

# for the propagators, we require the relative deviation to be no larger than
# in the 8th decimal place
prop_diff_rel_thresh="1e-8"
echo "Analysing differences in propagator data (relative deviation threshold: $prop_diff_rel_thresh)"
prop_differences=( )
for prop in ${txtprop_files[@]}; do
  echo Comparing: $prop "<->" reference/$prop
  numdiff -r $prop_diff_rel_thresh -X1:1 -X2:1 $prop reference/$prop
  echo ------------------------------------------------
  if [ $? -ne 0 ]; then
    prop_differences+=( $prop )
  fi
done
if [ ${#prop_differences[@]} -ne 0 ]; then
  echo "Propagator differences found in:"
  for i in $(seq 0 $(( ${#prop_differences[@]} - 1 )) ); do
    echo ${prop_differences[$i]}
  done
else
  echo "No differences found in propagator data!"
fi 
echo ------------------------------------------------
echo

if [ ${#differences[@]} -ne 0 -o ${#prop_differences[@]} -ne 0 ]; then
  echo "Errors found!"
  exit 1
else
  echo "No errors!"
  echo ------------------------------------------------
  echo
  exit 0
fi

