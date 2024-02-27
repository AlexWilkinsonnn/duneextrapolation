#!/bin/bash
cd $MRB_BUILDDIR
mrbsetenv && mrb i --generator=ninja && mrbslp
cd -

