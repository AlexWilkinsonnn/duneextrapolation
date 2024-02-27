Build with larsoft v09.78.04 e20:prof and a local copy of `dunereco` at the tag `v09_78_03d01` (needed for the CVN thinks). The other dependencies will be setup from the script `scripts/dev/setup.sh`.

This branch also requires a old tag of larsim with a hotfix. Clone `git@github.com:AlexWilkinsonnn/larsim.git` and checkout the branch `v09_38_01-elecdrift_fix`.

NOTE: comment out L241 & L243 in dunereco/CVN/art/PixelMapProducer.cxx for 1view cvn to work.

