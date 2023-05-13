#!/bin/bash
################################################################################
# Load depositions in from edep-sim -> detsim -> reco -> dump to flat tree
# Requires list of input file numbers in the tarball that can be used with the
# PROCESS environment variable instead of faffing with sam
################################################################################
# Options

OUTPUT_DIR="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/test/fd"

EDEP_DIR="/pnfs/dune/scratch/users/awilkins/lep_contained_pairs/test/edep"

LARSOFT_LOCAL_DIRNAME="duneextrapolation_larsoft"

################################################################################

echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"
ls

cp -r ${INPUT_TAR_DIR_LOCAL}/file_nums.txt .
line_num=$((PROCESS+1))
file_num=$(sed "${num}q;d" file_nums.txt)
edep_file=FHC.${file_num}.edep_flat.root

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source ${INPUT_TAR_DIR_LOCAL}/${LARSOFT_LOCAL_DIRNAME}/localProducts_larsoft_*/setup_grid
mrbslp

ups active

ifdh cp -D ${EDEP_DIR}/${edep_file}
edep_filepath=$(realpath $edep_file)

sed -e "s#physics.producers.IonAndScint.DepoDataLoc:.*#physics.producers.IonAndScint.DepoDataLoc: \"${edep_filepath}\"#" \
    ${INPUT_TAR_DIR_LOCAL}/${LARSOFT_LOCAL_DIRNAME}/srcs/duneextrapolation/duneextrapolation/MyWork/run_fcls/run_LoadChargeDepositions.fcl > \
    run_LoadChargeDepositions_local.fcl

n_evts=$(echo "std::cout << nd_depos->GetEntries() << std::endl;" | \
         root -l -b "$edep_file" 2>/dev/null | \
         tail -1)

echo "$n_evts in edep-sim file"

# XXX testing
n_evts=5

lar -c ./run_LoadChargeDepositions_local.fcl -n $n_evts

lar -c detsim_dune10kt_1x2x6_wirecell_refactored_nooptdet_dropSC.fcl -s LoadedDepos.root

lar -c reco_dune10kt_1x2x6_loadeddepos.fcl -s LoadedDepos_detsim.root

lar -c run_RecoDumpCVNE.fcl -s LoadedDepos_detsim_reco.fcl

ifdh cp LoadedDepos_detsim_reco_recodump.root ${OUTPUT_DIR}/FHC.${file_num}.LoadedDepos_detsim_reco.root

