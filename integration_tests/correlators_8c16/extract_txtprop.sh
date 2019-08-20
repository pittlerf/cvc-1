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

