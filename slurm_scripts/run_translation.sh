#!/bin/bash
#SBATCH -p gpu_gce
#SBATCH -N1
#SBATCH -c2
#SBATCH --gres=gpu:p100:1
#SBATCH --mail-user=alexander.wilkinson.20@ucl.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/slurm-%j.out

# #SBATCH --nodelist=wcgpu04

###############################################################################
# Slurm submission script that processes ND edep-sim in larnd-sim h5 format to
# art-root files with predicted FD response.
# Made for submission on Wilson Cluster
###############################################################################

indir=/pnfs/dune/persistent/users/awilkins/cafmaker/edep
outdir=/pnfs/dune/persistent/users/awilkins/cafmaker/translated

# Use script like: sbatch run_translation.sh krbcache proxyfile path1 path2 path3 ...
# Paths are relative to indir eg 0m/00/FHC.1000000.edep_dump.h5
RELATIVEPATHS=""
counter=1
for arg in "$@"; do
  if [[ $counter -eq 1 ]]; then
    KRBCACHE="$arg"
  elif [[ $counter -eq 2 ]]; then
    PROXY="$arg"
  else
    RELATIVEPATHS="${RELATIVEPATHS}${arg} "
  fi
  ((counter++))
done

###############################################################################
# Prepare for job

# Need to copy over valid cache files from kinit and setup_fnal_security to worker nodes for I/O
# with /pnfs
echo "Getting kerberos credentials from ${KRBCACHE}"
echo "Ensure kerberos credential cache is reachable (full path somewhere in /wclustre or /work1)"
cp -u $KRBCACHE /tmp/
# Race condition means the cp can fail when the cache file is already there, too lazy to thin
# of a good way around this so just using cp -u and not caring about error handling.
# if [ $? -eq 2 ]; then
#   echo \
#     "kerberos credential cache needs to be reachable (full path somewhere in /wclustre or /work1)"
#   exit 1
# fi

echo "Getting proxy from ${PROXY}"
echo "Ensure proxy is reachable (full path somewhere in /wclustre or /work1)"
cp -u $PROXY /tmp/
# if [ $? -ne 0 ]; then
#   echo "Proxy file needs to be reachable (full path somewhere in /wclustre or /work1)"
#   exit 1
# fi

# Setup environment
source /work1/dune/users/awilkins/extrapolation/setup.sh
source /work1/dune/users/awilkins/extrapolation/prep_libtorch_cu1141_cudnn_sm60sm70.sh
source /work1/dune/users/awilkins/larnd-sim/.venv_v3_9_2_larnd-sim/bin/activate

# Numba uses CUDA_HOME if cudatoolkit not found in a conda path
nvccloc=`which nvcc`
export CUDA_HOME=${nvccloc%/bin/nvcc}

# Cupy needs a valid HOME to write cache stuff to
export HOME=/scratch

# Prepare singularity
module load singularity/3.8.5
export SINGULARITY_BIND="/work1,/wclustre,/cvmfs,/scratch"

cd /scratch

###############################################################################
# Process files and copy results back

for relativepath in $RELATIVEPATHS; do
  echo "Reading ${relativepath} from ${indir}"

  infilename=${relativepath##*/}

  # Use sl7 singularity to do the ifdh cp since required software not on worker nodes
  singularity exec \
    --cleanenv \
    /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest \
    bash -c "
      export IFDH_CP_UNLINK_ON_ERROR=1
      export IFDH_CP_MAXRETRIES=1
      export IFDH_DEBUG=0
      source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh;
      setup ifdhc;
      ifdh cp ${indir}/${relativepath} /scratch/${infilename}
    "

  # ND detsim
  detsimout=${infilename%.edep_dump.h5}.larnd-sim.h5
  python /work1/dune/users/awilkins/larnd-sim/cli/simulate_pixels.py \
    --input_filename $infilename \
    --detector_properties /work1/dune/users/awilkins/larnd-sim/larndsim/detector_properties/ndlar-module.yaml \
    --pixel_layout /work1/dune/users/awilkins/larnd-sim/larndsim/pixel_layouts/multi_tile_layout-3.0.40.yaml \
    --response_file /work1/dune/users/awilkins/larnd-sim/larndsim/response_38.npy \
    --output_filename $detsimout

  detsimdump=${detsimout%.h5}.root
  python /work1/dune/users/awilkins/larpixsoft/export_depos_packets_toroot.py \
    -o $detsimdump \
    --ped 74 \
    --segment_length 0.04 \
    --no_cuts \
    $detsimout

  # Translate to FD response
  detsimloaded=${detsimdump%.root}.tolarsoft.root
  /work1/dune/users/awilkins/extrapolation/srcs/duneextrapolation/duneextrapolation/MyNDFDTranslation/run_fcls/run_LoadNDData_auto.sh \
    $detsimdump \
    -1 \
    $detsimloaded

  lar -c run_NDToFD_inference.fcl -s $detsimloaded
  translated=${detsimloaded%.root}_ndtranslated.root

  # Copy results back, use ifdh to create output dirs if not already there
  singularity exec \
    --cleanenv \
    /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7\:latest \
    bash -c "
      export IFDH_CP_UNLINK_ON_ERROR=1
      export IFDH_CP_MAXRETRIES=1
      export IFDH_DEBUG=0
      source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh;
      setup ifdhc;
      ifdh_mkdir_p() {
        local dir=$1
        local force=$2
        if [ `ifdh ls $dir 0 $force | wc -l` -gt 0 ]; then
          :
        else
          ifdh_mkdir_p `dirname $dir` $force
          ifdh mkdir $dir $force
        fi
      };
      ifdh_mkdir_p ${outdir}/${relativepath%/*}
      ifdh cp $translated ${outdir}/${relativepath%/*}/${translated}
    "

  # Need to clean up scratch area
  rm $translated
  rm $detsimloaded
  rm $detsimdump
  rm $detsimout
  rm $infilename
done
