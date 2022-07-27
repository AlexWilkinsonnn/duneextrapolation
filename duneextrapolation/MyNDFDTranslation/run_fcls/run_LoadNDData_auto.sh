#!/bin/bash
# Script for running run_LoadNDData.fcl correct
# Sets the number of events correctly and creates a temporary fcl to run over the desired file.

trap '{ rm -f /tmp/run_LoadNDData_ammended.tmp.fcl; }' EXIT

FILE=$1
NEVTS=$2
OUTNAME=$3

if [ ! -f $FILE ]; then
  FILE=${PWD}/${FILE}
  if [ ! -f $FILE ]; then
    echo "File not valid"
    exit 0
  fi
fi

MAXEVTS=$(echo "std::cout << ND_depos_packets->GetEntries() << std::endl;" | root -l -b "$FILE" 2>/dev/null | tail -1)

if [ -z $NEVTS ] || [ $NEVTS -eq -1 ]; then
  NEVTS=$MAXEVTS
elif [ $NEVTS -gt $MAXEVTS ]; then
  echo "Not enough events in $FILE, running over the max $MAXEVTS"
  NEVTS=$MAXEVTS
fi

sed "s/physics\.producers\.IonAndScint\.NDDataLoc:.*/physics\.producers\.IonAndScint\.NDDataLoc: \"${FILE}\"/" ${MRB_TOP}/srcs/duneextrapolation/duneextrapolation/MyNDFDTranslation/run_fcls/run_LoadNDData.fcl > /tmp/run_LoadNDData_ammended.tmp.fcl


if [ ! -z $OUTNAME ]; then
  lar -c /tmp/run_LoadNDData_ammended.tmp.fcl -n $NEVTS -o $OUTNAME
else
  lar -c /tmp/run_LoadNDData_ammended.tmp.fcl -n $NEVTS
fi
