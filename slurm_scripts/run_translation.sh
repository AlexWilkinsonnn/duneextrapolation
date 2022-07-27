#!/bin/bash
#SBATCH -p gpu_gce
#SBATCH -N1
#SBATCH -c2
#SBATCH --gres=gpu:v100:1
#SBATCH --mail-user=alexander.wilkinson.20@ucl.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/slurm-%j.out

# Paths relative to indir eg 0m/00/FHC.1000000.edep_dump.h5
RELATIVEPATHS="$@"

###############################################################################
# Prepare for job

indir=/pnfs/dune/persistent/users/awilkins/cafmaker/edep
outdir=/pnfs/dune/persistent/users/awilkins/cafmaker/translated

# ifdhc doen't have a mkdir -p equivalent, which is fine
# as long as you always remember to include this convenient function in your scripts
ifdh_mkdir_p() {
  local dir=$1
  local force=$2
  if [ `ifdh ls $dir 0 $force | wc -l` -gt 0 ]
  then
    : # we're done
  else
    ifdh_mkdir_p `dirname $dir` $force
    ifdh mkdir $dir $force
  fi
}

# Don't try over and over again to copy a file when it isn't going to work
export IFDH_CP_UNLINK_ON_ERROR=1
export IFDH_CP_MAXRETRIES=1
export IFDH_DEBUG=0

# Make symlinks to scripts
ln -s /work1/dune/users/awilkins/larnd-sim/cli/simulate_pixels.py /scratch/simulate_pixels.py
ln -s /work1/dune/users/awilkins/larpixsoft/export_depos_packets_toroot.py /scratch/export_depos_packets_toroot.py
ln -s /work1/dune/users/awilkins/extrapolation/srcs/duneextrapolation/duneextrapolation/MyNDFDTranslation/run_fcls/run_LoadNDData_auto.sh /scratch/run_LoadNDData_auto.sh

# Setup environment
source /work1/dune/users/awilkins/extrapolation/setup.sh
source /work1/dune/users/awilkins/extrapolation/prep_libtorch_cu1141_cudnn_sm60sm70.sh
source /work1/dune/users/awilkins/larnd-sim/.venv_v3_9_2_larnd-sim/bin/activate

# Numba uses CUDA_HOME if cudatoolkit not found in a conda path
nvccloc=`which nvcc`
export CUDA_HOME=${nvccloc%/bin/nvcc}

# Cupy needs a valid HOME to write cache stuff to
export HOME=/scratch

cd /scratch

###############################################################################
# Process files and copy results back

for relativepath in $RELATIVEPATHS; do
  echo "Reading ${relativepath} from ${indir}"

  infilename=${relativepath##*/}
  ifdh cp ${indir}/${relativepath} /scratch/${infilename}

  # ND detsim
  detsimout=${infilename%.edep_dump.h5}.larnd-sim.h5
  python simulate_pixels.py \
    --input_filename $infilename \
    --detector_properties /work1/dune/users/awilkins/larnd-sim/larndsim/detector_properties/ndlar-module.yaml \
    --pixel_layout /work1/dune/users/awilkins/larnd-sim/larndsim/pixel_layouts/multi_tile_layout-3.0.40.yaml \
    --response_file /work1/dune/users/awilkins/larnd-sim/larndsim/response_38.npy \
    --n_tracks 100 \
    --output_filename $detsimout

  detsimdump=${detsimout%.h5}.root
  python export_depos_packets_toroot.py \
    -o $detsimdump \
    --ped 74 \
    --segment_length 0.04 \
    --no_cuts \
    $detsimout

  # Translate to FD response
  detsimloaded=${detsimdump%.root}.tolarsoft.root
  ./run_LoadNDData_auto.sh $detsimdump -1 $detsimloaded

  lar -c run_NDToFD_inference.fcl -s $detsimloaded

  # copy results back
  translated=${detsimloaded%.root}_ndtranslated.root
  ifdh_mkdir_p ${outdir}/${relativepath%/*}
  ifdh cp $translated ${outdir}/${relativepath%/*}/${translated}
done
