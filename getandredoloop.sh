#!/bin/bash

# Cache the files needed for the 8:30 DDsnews trigger of the given month and day
# in 2019, then run jobs to process them all, running redo loops that look for
# failed jobs and resubmit them.  Try to do the redo loops gently enough that
# we don't fill up /tmp with overlapping submissions too often or cause other
# problems.

md="$1"
trigger=$2

unixtime=$(date +%s -d"$md 8:29:01 2019")

if ! $SRT_PRIVATE_CONTEXT/ligo/getfilelist.sh $unixtime $trigger; then
  echo getfilelist $unixtime $trigger failed
  exit 1
fi

$SRT_PRIVATE_CONTEXT/ligo/redoloop_ligo.sh \
strait-ligo-coincidence-artdaq-$unixtime-$trigger
