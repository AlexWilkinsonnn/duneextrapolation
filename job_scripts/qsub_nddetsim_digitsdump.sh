#!/bin/bash

#PBS -l walltime=01:00:00

#PBS -l cput=01:00:00

#PBS -l mem=2000mb

#PBS -l nodes=1

#PBS -m n

infile=$1

cd /unix/dune/awilkinson/extrapolation/larsoft_area
source setup.sh

lar -c run_ExportDigits.fcl -s $infile

