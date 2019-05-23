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

fclout=$SRT_PRIVATE_CONTEXT/job/$fcl

cat $SRT_PRIVATE_CONTEXT/job/ligojob_$type.fcl | \
  sed "/this_here_ligoanalysis: @local/a this_here_ligoanalysis.GWEventTime: \"$rfctime\"\
      \nthis_here_ligoanalysis.NeedBGEventTime: \"$realgweventtime\"\
      \nthis_here_ligofilter.GWEventTime: \"$rfctime\"\
      \nthis_here_ligoanalysis.SkyMap: \"$skymap\"" > $fclout

if grep -q "eliminatebeam.spillfile: " $fclout; then
  spilldir=/pnfs/nova/users/mstrait/spills
  spillbase=spills-*-${rfctime}.txt
  if ! ls $spilldir/$spillbase &> /dev/null; then
    echo Could not find a spill list file like
    echo $spilldir/$spillbase
    echo You cannot pass, flame of Udun, unless you run a job to
    echo make the spill file first.
    exit 1
  fi

  spillfile=$(ls $spilldir/$spillbase | head -n 1)

  if [ $(ls $spilldir/$spillbase | wc -l) -gt 1 ]; then
    echo Found more than one spill file with the same timestamp
    echo I do not know what is up with that. Using the first one.
  fi
  echo Using $spillfile

  # Use the base name only here since the grid will ship to to the PWD
  sed -i '/eliminatebeam.spillfile/s/".*"/"'$(basename $spillfile)'"/' $fclout
fi

echo $fclout created.  Here it is:
echo ========================== BEGIN FCL FILE ==========================
cat $fclout
echo =========================== END FCL FILE ===========================
