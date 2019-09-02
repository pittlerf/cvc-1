LIMEDIR=""
if [ -z "$1" ]; then
  echo "The first argument of extract_txtprops.sh must be the path of the c-lime installation directory!" 
  exit 1
else
  LIMEDIR="$1"
  if [ ! -d "${LIMEDIR}" ]; then
    "${LIMEDIR} could not be found!"
    exit 1
  fi
fi

echo LIMEDIR="${LIMEDIR}"

propfiles=( $(ls *.lime) )
if [ ${#propfiles[@]} -eq 0 ]; then
  echo "Propagator files were not produced!"
  exit 1
fi

for i in ${propfiles[@]}; do
  filebase=$(basename $i .lime)
  echo $i
  ${LIMEDIR}/bin/lime_extract_record $i 2 1 $filebase.bin
  od --format f8 --endian big $filebase.bin > $filebase.txtprop
  rm $filebase.bin
done

