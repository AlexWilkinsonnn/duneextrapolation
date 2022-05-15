#!/bin/bash

#PBS -l walltime=02:00:00

#PBS -l cput=02:00:00

#PBS -l mem=5000mb

#PBS -l nodes=1

#PBS -m n

nskipdir=$1

cd /unix/dune/awilkinson/extrapolation/larsoft_area
source setup.sh
cd $nskipdir

pwd
ls *3viewnd.root

lar -c cvn1view_nosp_dune10kt_1x2x6_truetranslated.fcl -s *3viewnd.root
lar -c cvn1view_nosp_dune10kt_1x2x6_networktranslated.fcl -s *cvn1viewtrue.root

