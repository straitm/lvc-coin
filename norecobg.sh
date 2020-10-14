#!/bin/bash

# Use the bg reco files generated for one event to generate others
#
# Takes the full path of a reco file and writes a hist file for
# the event event to the current directory.

if [ $# -ne 2 ]; then
  echo "Syntax: $(basename $0) gweventname infile"
  exit 1
fi

GWNAME=$1
infile=$2

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

if printf $infile | grep -q 'fardet.*t02'; then
  analysistype=minbiasfd_noreco
elif printf $infile | grep -q 'neardet.*ddactivity1'; then
  analysistype=ndactivity_noreco
else
  echo I do not know what to do with $infile
  echo I expected it to be either a fardet-t02 or neardet-ddactivity1 file
  exit 1
fi

if [ ${infile:0:1} != / ]; then
  echo Please give me the full path of the input file
  exit 1
fi

rfctime=$(basename $(dirname $infile) | cut -d- -f 1-5 | sed -e s/-/:/g -e s/:/-/ -e s/:/-/)

out=$rfctime-$(basename $infile .reco.root).hists.root
recoout=$rfctime-$(basename $infile .reco.root).reco.root
log=$rfctime-$(basename $infile .reco.root).log

if ! $SRT_PRIVATE_CONTEXT/ligo/makefcl.sh $analysistype $rfctime \
     $realgweventtime $skymap > $log; then
  exit 2
fi

fcl=$(printf ligojob_$analysistype.$rfctime.$realgweventtime.$(basename $skymap).fcl | sed s/:/-/g)
fqfcl=$SRT_PRIVATE_CONTEXT/job/$fcl

if ! nova -c $fqfcl $infile -o $recoout -T $out | grep -v "Begin processing" >> $log; then
  exit 1
fi

xz $log

rm -f $recoout
