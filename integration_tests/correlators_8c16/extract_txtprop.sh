for i in *.lime; do
  filebase=$(basename $i .lime)
  echo $i
  ~/build/lime/broadwell/examples/lime_extract_record $i 2 1 $filebase.bin
  od --format f8 --endian big $filebase.bin > $filebase.txtprop
  rm $filebase.bin
done

