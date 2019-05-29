#!/bin/bash

# I don't understand why sometimes novasoft-ligo.tar is written to the PWD
# and not /tmp, but try to make sure that it is somewhere reasonable
cd /nova/ana/users/mstrait/ligometalog/

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

if [ $# -ne 2 ] && [ $# -ne 3 ]; then
  echo "Syntax: $(basename $0) unixtime analysis_type [def]"
  echo "Where unixtime is the time to center the analysis on,"
  echo "      analysis_type is fardet-t02, fardet-ddsnews, fardet-ddenergy"
  echo "                       neardet-ddactivity1 or neardet-ddsnews."
  echo "      def, optional, is the SAM definition to use"
  echo "           (This feature is for the redo definitions)"
  exit 2
fi

unixtime=$1
fracsec=$(cut -d. -f 2 -s <<< $unixtime)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
rfctimesafeforsam=${rfctime//:/-}
analysis_type_key=$2
def=strait-ligo-coincidence-artdaq-${unixtime}-${analysis_type_key}
if [ $3 ]; then
  def=$3
fi

# Don't use noreco for ND because of a bizarre problem in which
# it can't see the slicer product in reco files.
if printf $def | grep fardet | grep -q -- -reco-; then
  typesuffix=_noreco
else
  recoout='--outTier out1:reco'
fi

if [ $analysis_type_key == fardet-t02 ]; then
  type=minbiasfd$typesuffix
  if [ "$typesuffix" == _noreco ]; then
    lifetime=21600
  else
    lifetime=40000
  fi
elif [ $analysis_type_key == fardet-ddsnews ] ||
     [ $analysis_type_key == fardet-ligo ]; then
  disk='--disk 20480'
  if [ "$typesuffix" == _noreco ]; then
    type=minbiasfd$typesuffix
    lifetime=21600
  else
    type=minbiasfd_rawinput
    lifetime=172800
  fi
elif [ $analysis_type_key == fardet-ddenergy ]; then
  type=ddenergy
  lifetime=14400
elif [ $analysis_type_key == neardet-ddactivity1 ]; then
  type=ndactivity$typesuffix
  lifetime=21600
elif [ $analysis_type_key == neardet-ddsnews ] ||
     [ $analysis_type_key == neardet-ligo ]; then
  type=minbiasnd$typesuffix
  lifetime=14400
else
  echo I cannot figure out what analysis type to run for $analysis_type_key
  echo which I got from $def
  exit 2
fi

if [ $GWBLIND ]; then
  echo Doing blind analysis '(livetime report only)'
  if ( [ $analysis_type_key == fardet-ddsnews ] ||
       [ $analysis_type_key == fardet-ligo ] ) &&
       [ $typesuffix != _noreco ]; then
    type=blind_rawinput
  else
    type=blind
  fi
  lifetime=7200
fi

# Tell the fcl we are reading the skymap from the CWD since it will be shipped
# there by the grid.  (*Don't* do that for interactive jobs.)
if ! $SRT_PRIVATE_CONTEXT/ligo/makefcl.sh $type $rfctime \
     $realgweventtime $(basename $skymap); then
  exit 2
fi

fcl=$(printf ligojob_$type.$rfctime.$realgweventtime.$(basename $skymap).fcl | sed s/:/-/g)
fqfcl=$SRT_PRIVATE_CONTEXT/job/$fcl

if grep -q "eliminatebeam.spillfile: " $fqfcl; then
  spilldir=/pnfs/nova/users/mstrait/spills
  spillbase=$(grep eliminatebeam.spillfile: $fqfcl | cut -d\" -f 2)
  spillfile=$spilldir/$spillbase
  if ! [ -e $spilldir/$spillbase ]; then
    echo I am confused.  We just made a fcl assuming
    echo $spillfile existed, but it does not.
    exit 2
  fi
else
  echo Not using EliminateBeamSpills
fi

njobs=$(samweb list-files defname: $def | wc -l)
tag=$(cat $SRT_PRIVATE_CONTEXT/.base_release)
testrel=/nova/app/users/mstrait/novasoft-ligo/

outdir=$outhistdir/$rfctimesafeforsam-$analysis_type_key

mkdir -p $outdir

if [ -e $outdir/bad ]; then
  echo $(basename outdir) is marked bad, will not submit jobs
  exit 2
fi
if [ -e $outdir/complete ]; then
  echo $(basename outdir) is marked complete, will not submit jobs
  exit 2
fi

# This is needed because the way submission works these days is that it
# tars up the whole test release in /tmp, which only has 2GB.
fails=0

# Try to avoid races below
if ! [ $REDOFAST ]; then sleep $((RANDOM % 32 )); fi

while true; do
  kfree=$(df /tmp | tail -n 1 | awk '{print $3}')
  if [ $kfree -lt 1500000 ]; then
    echo /tmp is too full.  Waiting for it to clear out
    df -h /tmp | tail -n 1
    sleep 1m
  elif ls -l /tmp | grep 'mstrait.*tmp......$'; then
    echo waiting for another submit to finish dumping in /tmp
    sleep 1m
  else
    break
  fi
  let fails++
  if [ $fails -gt 10 ]; then
    echo OK, I have had it.  Going to clear my files out of /tmp
    rm -fv /tmp/tmp* /tmp/mstrait.* /tmp/joblist.* /tmp/headtail.*
    fails=0
  fi
done


submit_nova_art.py \
--inputfile $skymap \
$(if [ $spillfile ]; then echo --inputfile $spillfile; fi) \
--opportunistic \
--maxopt \
--logs \
--histTier hists \
$recoout \
$disk \
--copyOut \
--jobname $GWNAME-$def \
--defname $def \
--njobs $njobs \
--files_per_job 1 \
--tag $tag \
--testrel $testrel \
--config $fcl \
--dest $outdir \
--expected_lifetime $lifetime
