#!/bin/bash

threads=24
fin=$1
dout=scarap

scarap core $fin $dout/core -t $threads

dsm=$dout/supermatrix
scarap concat $fin $dout/core/genes.tsv $dsm -t $threads

# trim gaps
trimal \
  -in $dsm/supermatrix_aas.fasta \
  -out $dsm/supermatrix_aas_trimmed.fasta \
  -gt 0.90 -keepheader

iqtree \
  -s $dsm/supermatrix_aas_trimmed.fasta \
  -pre $dout/lab \
  -m LG+F+G4 \
  -alrt 1000 -bb 1000 -nt $threads
