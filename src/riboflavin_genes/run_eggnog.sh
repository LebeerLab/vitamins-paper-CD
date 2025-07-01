dout=results
mkdir -p $dout

while read assembly; do

  outf="$(basename $assembly)"
  outf="${outf%.*}" 

  if [ -s $dout/$outf.emapper.annotations ]; then  
    echo "Skipping $outf"
  else
    aa_f="${assembly%assembly/*}annotation/${outf%_*}.faa.gz"
    emapper.py -i $aa_f --output $dout/$outf --cpu 18 --dbmem --tax_scope Lactobacillaceae --num_workers 6 --temp_dir tmp
  fi
done < $1
