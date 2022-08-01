#!/bin/bash

KRBCACHE=$1
PROXY=$2
FILELIST=$3
N=$4

for ((i=0;i<$N;i++)); do
  sbatch run_translation.sh $KRBCACHE $PROXY $(sed -e 1$'{w/dev/stdout\n;d}' -i~ ${FILELIST})

  # Caveman way of preventing jobs all trying to copy data into the same disk at once
  if [[ $(( $i + 1 )) -lt 8 ]]; then
    sleep 30
  fi
done
