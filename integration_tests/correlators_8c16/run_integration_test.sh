#!/bin/bash

if [ -z "$CORRBIN" ]; then
  echo "The CORRBIN environment variable must be defined and point to the 'correlators' executable!"
fi

if [ -z "$LIMEDIR" ]; then
  echo "The LIMEDIR environment variable must be defined and point to the c-lime installation directory!" 
fi

if [ -z "$CORRBIN" -o -z "$LIMEDIR" ]; then
  exit 1
fi


# remove any existing h5 files, just in case (the code does not overwrite!)
rm *.h5

mpirun -np 1 $CORRBIN

echo "Analysing differences in correlator data"

corrfiles=( $(ls *.h5) )
if [ ${#corrfiles[@]} -eq 0 ]; then
  echo "Correlator files were not produced!"
  exit 1
fi
echo ${corrfiles[@]}

corr_differences=( )
for corrfile in ${corrfiles[@]}; do
  h5diff -p 1e-6 $corrfile reference/$corrfile
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

echo "Analysing differences in propagator data"
prop_differences=( )
for prop in ${txtprop_files[@]}; do
  numdiff -r1e-10 -X1:1 -X2:1 $prop reference/$prop
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

