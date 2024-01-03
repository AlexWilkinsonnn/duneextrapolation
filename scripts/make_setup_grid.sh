#!/bin/bash

SETUP_SCRIPT=$1
OUTPUT_FILE=$2

local_prods_dirpath=${SETUP_SCRIPT%/*}
local_prods_dirname=${local_prods_dirpath##*/}

sed \
  "s#\(setenv MRB_TOP \).*#\1\"\${INPUT_TAR_DIR_LOCAL}\"#g;
   s#\(setenv MRB_TOP_BUILD \).*#\1\"\${INPUT_TAR_DIR_LOCAL}\"#g;
   s#\(setenv MRB_SOURCE \).*#\1\"\${INPUT_TAR_DIR_LOCAL}/srcs\"#g;
   s#\(setenv MRB_INSTALL \).*#\1\"\${INPUT_TAR_DIR_LOCAL}/${local_prods_dirname}\"#g;
   s#\(setenv CETPKG_INSTALL \).*#\1\"\${INPUT_TAR_DIR_LOCAL}/${local_prods_dirname}\"#g" \
  $SETUP_SCRIPT > $OUTPUT_FILE

