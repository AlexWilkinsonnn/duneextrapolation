#!/bin/bash
################################################################################
# Script to take ND-FD pair HDF5 with truth information for ND and FD and add
# the SimChannels for the FD deposition back to the HDF5 file. Charge
# depositions are loaded into larsoft and ionisation + drift modules are
# ran to get SimChannels.
################################################################################
# Options

FD_PAIR_SIMCHANNELS_OUTPUT="/pnfs/dune/scratch/users/awilkins/larbath_ndfd_pairs/tdr_sample/pair_simchannels_ndfd"
FD_SIMCHANNELS_OUTPUT="/pnfs/dune/scratch/users/awilkins/larbath_ndfd_pairs/tdr_sample/fdsimchannels_artroot"

SAVE_FDSC=false
SAVE_PAIR_SC=true # turn this on if not testing!

INTERACTIVE=false # not running on grid to test (run from inside dir that would be in tarball)

INPUT_PAIR_H5_DIR=$1

################################################################################

if [ "$INTERACTIVE" = true ]; then
  PROCESS=0 # or 1?
  export INPUT_TAR_DIR_LOCAL=${PWD}
else
  echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}. At ${PWD}"
fi

echo $INPUT_TAR_DIR_LOCAL

# Setup env
${INPUT_TAR_DIR_LOCAL}/srcs/duneextrapolation/scripts/make_setup_grid.sh ${INPUT_TAR_DIR_LOCAL}/localProducts_larsoft_*/setup \
                                                                         setup-grid
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source setup-grid
setup dunesw v09_78_03d01 -q e20:prof
setup duneprototypes v09_78_03d01 -q e20:prof
setup duneana v09_78_03d01 -q e20:prof
setup dunesim v09_78_03d01 -q e20:prof
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
sed -i "s#physics.producers.largeant.NDFDH5FileLoc: \"\"#physics.producers.largeant.NDFDH5FileLoc: \"${input_file_local}\"#" run_LoadFDDepos_oldgeo.fcl
sed -i "s#physics.analyzers.addsc.NDFDH5FileLoc: \"\"#physics.analyzers.addsc.NDFDH5FileLoc: \"${input_file_local}\"#" run_AddFDSimChannels_oldgeo.fcl

ls -lrth

num_events=$(h5ls-shared $input_name | sed -n "s/^vertices.*Dataset {\([0-9]\+\).*/\1/p")
echo "$input_name has $num_events events"

# Generate FD SimChannels
lar -c ./run_LoadFDDepos_oldgeo.fcl -n $num_events
lar -c ionandscint_elecdrift_dune10kt_1x2x6oldgeo.fcl -s LoadedFDDeps.root -n -1

# Add FD SimChannels to H5 file
lar -c ./run_AddFDSimChannels_oldgeo.fcl -s LoadedFDDeps_g4.root -n -1

ls -lrth

echo "Copying files to dCache..."
if [ "$SAVE_FDSC" = true ]; then
  ifdh cp LoadedFDDeps_g4.root \
          ${FD_SIMCHANNELS_OUTPUT}/${input_name%.*}_LoadedFDDeps_g4.root
fi
if [ "$SAVE_PAIR_SC" = true ]; then
  ifdh cp ${input_name} ${FD_PAIR_SIMCHANNELS_OUTPUT}/${input_name%.*}_fdsimchannels.h5
fi

