#!/bin/bash

# Cache the files needed for the 8:30 DDsnews trigger of the
# given date, then run jobs to process them all, running redo
# loops that look for failed jobs and resubmit them. Try to
# do the redo loops gently enough that we don't fill up /tmp with overlapping
# submissions too often or cause other problems.

cd /nova/ana/users/mstrait/ligometalog/

if [ $# -ne 5 ]; then
  echo Syntax: $(basename $0) month day trigger year gwname
  echo "For example: $(basename $0) Jan 1 fardet-ddsnews 2019 S190426c"
  echo "  gwname defines the sky map and effective time for pointing"
  exit 1
fi

echo $0 running on $HOSTNAME

month=$1
day=$2
trigger=$3
year=$4
export GWNAME=$5

# Just (for the moment) to check GWNAME
. $SRT_PRIVATE_CONTEXT/ligo/env.sh

unixtime=$(date +%s -d"$month $day 8:29:01 $year")

(
if ! $SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime $trigger; then
  echo getfilelist $unixtime $trigger failed
  exit 1
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

$SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh $def
) 2> /dev/stdout | tee $GWNAME-$month-$day-$year-$unixtime-$trigger.log
