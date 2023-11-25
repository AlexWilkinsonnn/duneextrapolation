#!/bin/bash

#PBS -l walltime=08:00:00

#PBS -l cput=08:00:00

#PBS -l mem=3500mb

#PBS -l nodes=1

#PBS -m n

infile=$1
nskip=$2
n=$3

cd /unix/dune/awilkinson/extrapolation/larsoft_area
source setup.sh
mkdir nskip${nskip}
cd nskip${nskip}

nend=$(($n+$nskip))
infilelocal=${infile##*/}
lar -c detsim_dune10kt_1x2x6_wirecell_refactored_nooptdet.fcl -s $infile -n $n --nskip $nskip -o ${infilelocal::-5}_${nskip}-${nend}_detsim.root

