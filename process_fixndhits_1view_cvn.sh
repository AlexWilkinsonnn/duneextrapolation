#!/bin/bash
nskipdir=$1

cd /unix/dune/awilkinson/extrapolation/larsoft_area
source setup.sh
cd $nskipdir
echo $nskipdir

ls -lrt
lar -c run_ShiftHitChannels.fcl -s NDFDTranslation_gen_g4_detsim_trans14latestT10P2train_cvn1viewtrue_cvn1viewnetwork.root
lar -c cvn1view_dune10kt_1x2x6_ndpacketsshifted.fcl -s NDFDTranslation_gen_g4_detsim_trans14latestT10P2train_cvn1viewtrue_cvn1viewnetwork_hitshifted.root

