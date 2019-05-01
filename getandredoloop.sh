#!/bin/bash

# Cache the files needed for the 8:30 DDsnews trigger of the given month and
# day in 2019 (or another year if specified), then run jobs to process them
# all, running redo loops that look for failed jobs and resubmit them.  Try to
# do the redo loops gently enough that we don't fill up /tmp with overlapping
# submissions too often or cause other problems.

cd /nova/ana/users/mstrait/ligometalog/

if [ $# -ne 5 ]; then
  echo Syntax: $(basename $0) month day trigger year gwname
  echo "For example: $(basename $0) Jan 1 fardet-ddsnews 2019 S190426c"
  exit 1
fi

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

$SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh \
strait-ligo-coincidence-artdaq-$unixtime-$trigger
) 2> /dev/stdout | tee $GWNAME-$month-$day-$year-$unixtime-$trigger.log
