#!/bin/bash

FILELIST=$1
N=$2

for ((i=0;i<$N;i++)); do
  sbatch run_translation.sh $(sed -e 1$'{w/dev/stdout\n;d}' -i~ ${FILELIST})
done
