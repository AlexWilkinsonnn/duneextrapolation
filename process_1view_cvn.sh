#!/bin/bash
infile=$1
nskip=$2
n=$3

cd /unix/dune/awilkinson/extrapolation/larsoft_area
source setup.sh
mkdir nskip${nskip}
cd nskip${nskip}

infilelocal=${infile##*/}
nend=$(($n+$nskip))
lar -c cvn1view_dune10kt_1x2x6_truetranslated.fcl -s $infile -n $n --nskip $nskip
lar -c cvn1view_dune10kt_1x2x6_networktranslated.fcl -s ${infilelocal::-5}_cvn1viewtrue.root -n $n
lar -c cvn1view_dune10kt_1x2x6_ndpackets.fcl -s ${infilelocal::-5}_cvn1viewtrue_cvn1viewnetwork.root -n $n -o ${infilelocal::-5}_cvn1viewtrue_cvn1viewnetwork_cvn1viewnd_${nskip}-${nend}.root

