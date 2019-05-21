#!/bin/bash

# Cache the files needed for the 8:30 DDsnews trigger of the given date, or the
# specified time, then run jobs to process them all, running redo loops that
# look for failed jobs and resubmit them. Try to do the redo loops gently
# enough that we don't fill up /tmp with overlapping submissions too often or
# cause other problems.

cd /nova/ana/users/mstrait/ligometalog/

if [ $# -ne 5 ] && [ $# -ne 6 ]; then
  echo Syntax: $(basename $0) month day trigger year gwname '[UTC time]'
  echo "For example: $(basename $0) Jan 1 fardet-ddsnews 2019 S190426c"
  echo "  gwname defines the sky map and effective time for pointing"
  echo "If UTC time is not given, the 8:30 SNEWS trigger time is used"
  exit 1
fi

echo $0 running on $HOSTNAME

month=$1
day=$2
trigger=$3
year=$4
export GWNAME=$5
time=$6

# Just (for the moment) to check GWNAME
. $SRT_PRIVATE_CONTEXT/ligo/env.sh

if [ $time ]; then
  unixtime=$(TZ=UTC date +%s -d"$month $day $time $year")
  fracsec=$(cut -d. -f 2 -s <<< $time)
  unixtime+=.$fracsec
else
  unixtime=$(date +%s -d"$month $day 8:29:01 $year")
fi

(
$SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime $trigger
ret=$?
if [ $ret -eq 1 ]; then
  echo getfilelist $unixtime $trigger failed
  exit 1
elif [ $ret -eq 2 ]; then
  # if getfilelist returned 2, it means it found no files to process
  exit 0
fi

defbase=strait-ligo-coincidence
recodef=$defbase-reco-$unixtime-$trigger
rawdef=$defbase-artdaq-$unixtime-$trigger

# getfilelist.sh only allows good definitions to exist
if samweb list-definitions | grep -q $recodef; then
  def=$recodef
else
  def=$rawdef
fi

if $SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh $def; then
  rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
  rfctimesafeforsam=${rfctime//:/-}
  $SRT_PRIVATE_CONTEXT/ligo/combine.sh \
    $outhistdir/$rfctimesafeforsam-$trigger
fi
) 2> /dev/stdout | tee /nova/ana/users/mstrait/ligometalog/$GWNAME-$month-$day-$year-$unixtime-$trigger.$(date "+%Y-%m-%dT%H:%M:%S").log
