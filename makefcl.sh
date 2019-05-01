#!/bin/bash

# do *not* source env.sh in this script

if [ $# -ne 4 ]; then
  echo Syntax: $(basename $0) analysis_type timestamp gwtimestamp skymap
  echo "Where analysis_type is the XXX in ligojob_XXX.fcl" 
  echo "      timestamp is the RFC time to center the analysis on" 
  echo "      gwtimestamp is the RFC time of the GW event" 
  echo "        (Those are the same time for the signal sample)"
  echo "      skymap is the FITS file"
  exit 1
fi

type=$1
rfctime=$2
realgweventtime=$3
skymap=$4

fcl=ligojob_$type.$rfctime.fcl

if ! [ -e $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl ]; then
  echo Analysis type $type not supported
  exit 1
fi

cat $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl | \
  sed "/this_here_ligoanalysis: @local/a this_here_ligoanalysis.GWEventTime: \"$rfctime\"\
      \nthis_here_ligoanalysis.NeedBGEventTime: \"$realgweventtime\"\
      \nthis_here_ligofilter.GWEventTime: \"$rfctime\"\
      \nthis_here_ligoanalysis.SkyMap: \"$skymap\"" > $SRT_PRIVATE_CONTEXT/job/$fcl

echo $SRT_PRIVATE_CONTEXT/job/$fcl created
