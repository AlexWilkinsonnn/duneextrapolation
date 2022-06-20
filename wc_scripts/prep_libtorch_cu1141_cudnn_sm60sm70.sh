#!/bin/bash

module load cuda11/11.1.4
unsetup libtorch

export LIBTORCH_INC=/work1/dune/users/awilkins/extrapolation/pytorch-install/include
export LIBTORCH_DIR=/work1/dune/users/awilkins/extrapolation/pytorch-install
export LIBTORCH_FQ_DIR=/work1/dune/users/awilkins/extrapolation/pytorch-install
export LIBTORCH_LIB=/work1/dune/users/awilkins/extrapolation/pytorch-install/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work1/dune/users/awilkins/extrapolation/pytorch-install/lib

export CUDNN_LIBRARY=/work1/dune/users/awilkins/dsim/cudnn/lib64
export CUDNN_INCLUDE_DIR=/work1/dune/users/awilkins/dsim/cudnn/include
export CUDNN_LIB_DIR=/work1/dune/users/awilkins/dsim/cudnn/lib64

