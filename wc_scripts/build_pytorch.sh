#!/bin/bash

module load cuda11/11.4.1

unsetup libtorch
setup cmake v3_22_0

# Need the venv with pyyaml
source /work1/dune/users/awilkins/dsim/bin/activate

export CUDNN_LIBRARY=/work1/dune/users/awilkins/dsim/cudnn/lib64
export CUDNN_INCLUDE_DIR=/work1/dune/users/awilkins/dsim/cudnn/include
export CUDNN_LIB_DIR=/work1/dune/users/awilkins/dsim/cudnn/lib64

mkdir -p pytorch-build
mkdir -p pytorch-install
cd pytorch-build

cmake -DBUILD_CUSTOM_PROTOBUF=OFF \
      -DBUILD_PYTHON=OFF \
      -DPYTHON_EXECUTABLE=`which python3` \
      -DUSE_DISTRIBUTED=OFF \
      -DCUDA_HOME=/srv/software/cuda-toolkits/11.4.1 \
      -DCUDNN_LIB_DIR=/work1/dune/users/awilkins/dsim/cudnn/lib64 \
      -DCUDNN_INCLUDE_DIR=/work1/dune/users/awilkins/dsim/cudnn/include \
      -DTORCH_CUDA_ARCH_LIST="6.0;7.0" \
      -DUSE_CUDA=ON \
      -DUSE_CUDNN=ON \
      -DUSE_OPENCV=OFF \
      -DBUILD_TORCH=ON \
      -DCUDA_HOST_COMPILER=cc \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SHARED_LIBS=ON \
      -DCMAKE_INSTALL_PREFIX=../pytorch-install \
      ../pytorch

cmake --build . --target install -j 6

cd ..
