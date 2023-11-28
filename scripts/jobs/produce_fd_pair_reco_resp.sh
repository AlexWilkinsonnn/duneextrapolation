#!/bin/bash
################################################################################
# Script to take ND-FD pair HDF5 with truth information for ND and FD and add
# the FD reconstruction and detector response. Charge depositions are loaded
# into larsoft and ran through detector simulation and reconstruction. This
# should be the final step in making ND-FD paired data, leaving ND packets,
# ND packets projected to wires (aligned with FD response), FD wire response,
# FD reco.
################################################################################
# Options

COMPLETE_PAIR_OUTPUT="/pnfs/dune/scratch/users/awilkins/larbath_ndfd_pairs/test_sample_20000/pair_complete_ndfd"
FDRECO_PAIR_OUTPUT="/pnfs/dune/scratch/users/awilkins/larbath_ndfd_pairs/test_sample_20000/fdreco_artroot"
FDRESP_PAIR_OUTPUT="/pnfs/dune/scratch/users/awilkins/larbath_ndfd_pairs/test_sample_20000/fdresp_artroot"

SAVE_FDRECO=false
SAVE_FDRESP=false # this will take a lot of disk
SAVE_COMPLETE_PAIR=true # turn this on if not testing!

INTERACTIVE=false # not running on grid to test (run from inside dir that would be in tarball)

INPUT_PAIR_H5_DIR=$1

################################################################################

if [ "$INTERACTIVE" = true ]; then
  PROCESS=0 # or 1?
  export INPUT_TAR_DIR_LOCAL=${PWD}
else
  echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}. At ${PWD}"
fi

# Setup env
${INPUT_TAR_DIR_LOCAL}/srcs/duneextrapolation/scripts/make_setup_grid.sh localProducts_larsoft_*/setup
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source ${INPUT_TAR_DIR_LOCAL}/localProducts_larsoft_*/setup-grid
setup dunesw v09_78_03d01 -q e20:prof
setup duneprototypes v09_78_03d01 -q e20:prof
setup duneana v09_78_03d01 -q e20:prof
setup dunesim v09_78_03d01 -q e20:prof
setup dunecore v09_78_03d01 -q e20:prof
mrbslp

# Don't try over and over again to copy a file when it isn't going to work
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

input_file=$(ifdh ls $INPUT_PAIR_H5_DIR | head -n $((PROCESS+2)) | tail -n -1)
input_name=${input_file##*/}

echo "input_file is ${input_file}"

ifdh cp $input_file $input_name
input_file_local=$PWD/$input_name

# Prepare fcls
cp ${INPUT_TAR_DIR_LOCAL}/srcs/duneextrapolation/duneextrapolation/NDFDPairs/run_fcls/*.fcl .
sed -i "s#physics.producers.largeant.NDFDH5FileLoc: \"\"#physics.producers.largeant.NDFDH5FileLoc: \"${input_file_local}\"#" run_LoadFDDepos.fcl
sed -i "s#physics.producers.largeant.NDFDH5FileLoc: \"\"#physics.producers.largeant.NDFDH5FileLoc: \"${input_file_local}\"#" run_LoadFDDepos_InsideNDOnly.fcl
sed -i "s#physics.analyzers.addreco.NDFDH5FileLoc: \"\"#physics.analyzers.addreco.NDFDH5FileLoc: \"${input_file_local}\"#" run_AddFDReco.fcl
sed -i "s#physics.analyzers.addresp.NDFDH5FileLoc: \"\"#physics.analyzers.addresp.NDFDH5FileLoc: \"${input_file_local}\"#" run_AddFDResp.fcl

ls -lrth

num_events=$(h5ls-shared $input_name | sed -n "s/^vertices.*Dataset {\([0-9]\+\)}/\1/p")
echo "$input_name has $num_events events"

# Generate FD reco
lar -c ./run_LoadFDDepos.fcl -n $num_events
lar -c ionandscint_dune10kt_1x2x6.fcl -s LoadedFDDeps.root -n -1
lar -c standard_detsim_dune10kt_1x2x6_nooptdet.fcl -s LoadedFDDeps_g4.root -n -1
lar -c standard_reco_dune10kt_nu_1x2x6_nooptreco.fcl -s LoadedFDDeps_g4_detsimnoopt.root -n -1

# Add FD reco to H5 file
lar -c ./run_AddFDReco.fcl -s LoadedFDDeps_g4_detsimnoopt_reconoopt.root -n -1

ls -lrth

# Generate FD detector response
lar -c ./run_LoadFDDepos_NDLAronly.fcl -n $num_events
lar -c ionandscint_dune10kt_1x2x6.fcl -s LoadedFDDepsInsideNDOnly -n -1
lar -c detsim_dune10kt_1x2x6_notpcsigproc_nooptdet.fcl -s LoadedFDDepsInsideNDOnly_g4.root -n -1

# Add FD detector response and wire projected and aligned packets to H5 file
lar -c ./run_AddFDResp.fcl -s LoadedFDDepsInsideNDOnly_g4_detsimnooptnosp.root -n -1

ls -lrth

echo "Copying files to dCache..."
if [ "$SAVE_FDRECO" = true ]; then
  ifdh cp LoadedFDDeps_g4_detsimnoopt_reconoopt.root \
          ${FDRECO_PAIR_OUTPUT}/${input_name%.*}_LoadedFDDeps_g4_detsimnoopt_reconoopt.root
fi
if [ "$SAVE_FDRESP" = true ]; then
  ifdh cp LoadedFDDepsInsideNDOnly_g4_detsimnooptnosp.root \
          ${FDRESP_PAIR_OUTPUT}/${input_name%.*}_LoadedFDDepsInsideNDOnly_g4_detsimnooptnosp.root
fi
if [ "$SAVE_COMPLETE_PAIR" = true ]; then
  ifdh cp ${input_name} ${COMPLETE_PAIR_OUTPUT}/${input_name%.*}_fdreco_fdresp.h5
fi

