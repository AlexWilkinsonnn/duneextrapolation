#!/bin/bash
#SBATCH -p gpu_gce
#SBATCH -N1
#SBATCH -c2
#SBATCH --gres=gpu:v100:1
#SBATCH --mail-user=alexander.wilkinson.20@ucl.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH -o ./logs/slurm-%j.out

# Path relative to indir eg 0m/00/FHC.1000000.edep_dump.h5
RELATIVEPATH=$1

indir=/wclustre/dune/awilkins/extrapolation/cafmaker_jobs/data/edep
outdir=/wclustre/dune/awilkins/extrapolation/cafmaker_jobs/data/translated

infilename=${RELATIVEPATH##*/}

echo "Reading ${RELATIVEPATH} from ${indir}"

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

cp ${indir}/${RELATIVEPATH} /scratch/${infilename}
cd /scratch

detsimout=${infilename%.edep_dump.h5}.larnd-sim.h5
python simulate_pixels.py \
  --input_filename $infilename \
  --detector_properties /work1/dune/users/awilkins/larnd-sim/larndsim/detector_properties/ndlar-module.yaml \
  --pixel_layout /work1/dune/users/awilkins/larnd-sim/larndsim/pixel_layouts/multi_tile_layout-3.0.40.yaml \
  --response_file /work1/dune/users/awilkins/larnd-sim/larndsim/response_38.npy \
  --output_filename $detsimout

detsimdump=${detsimout%.h5}.root
python export_depos_packets_toroot.py \
  -o $detsimdump \
  --ped 74 \
  --segment_length 0.04 \
  --no_cuts \
  $detsimout

detsimloaded=${detsimdump%.root}.tolarsoft.root
./run_LoadNDData_auto.sh $detsimdump -1 $detsimloaded

echo
ls
echo

lar -c run_NDToFD_inference.fcl -s $detsimloaded

translated=${detsimloaded%.root}_ndtranslated.root
mkdir -p ${outdir}/${RELATIVEPATH%/*}
cp $translated ${outdir}/${RELATIVEPATH%/*}/${translated}

