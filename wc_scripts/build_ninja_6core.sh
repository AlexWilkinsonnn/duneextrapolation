#!/bin/bash

cd $MRB_BUILDDIR
mrbsetenv && mrb i --generator=ninja -j6 && mrbslp
cd -
