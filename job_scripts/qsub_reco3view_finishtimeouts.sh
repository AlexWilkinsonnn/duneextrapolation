#!/bin/bash

#PBS -l walltime=06:00:00

#PBS -l cput=06:00:00

#PBS -l mem=4000mb

#PBS -l nodes=1

#PBS -m n

infile=$1
nskip=$2
n=$3

cd /unix/dune/awilkinson/extrapolation/larsoft_area
source setup.sh
cd nskip${nskip}
ls -lrt *

nend=$(($n+$nskip))
infilelocal=${infile##*/}
lar -c reco3view_dune10kt_1x2x6_networktranslated.fcl -s ${infilelocal::-5}_${nskip}-${nend}_detsim_3viewtrue_3viewalt.root
lar -c reco3view_dune10kt_1x2x6_ndpackets.fcl -s ${infilelocal::-5}_${nskip}-${nend}_detsim_3viewtrue_3viewalt_3viewnetwork.root

