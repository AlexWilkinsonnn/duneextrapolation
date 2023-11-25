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
SAVE_COMPLETE_PAIR=false # turn this on if not testing!

INTERACTIVE=true # not running on grid to test (run from inside dir that would be in tarball)

INPUT_PAIR_H5_DIR=$1

################################################################################

if [ "$INTERACTIVE" = true ]; then
  PROCESS=0 # or 1?
else
  echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"
  cd $INPUT_TAR_DIR_LOCAL
fi

source setup.sh

# Don't try over and over again to copy a file when it isn't going to work
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

input_name=$(ifdh ls $INPUT_PAIR_H5_DIR | head -n $(PROCESS+2) | tail -n -1)
input_file=${INPUT_DIR}/${input_name}

echo "Copying files to dCache..."
