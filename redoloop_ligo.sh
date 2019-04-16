#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

# SAM definition that we need to do and redo
realdef=$1
unixtime=$(echo $realdef | cut -d- -f 5)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
stream=$(echo $realdef | cut -d- -f 6-7)

# If empty, the right thing happens (which is that we use a dummy skymap)
skymap=$2

TMP=/tmp/mstrait.redolist.$$

iteration=1

find_redo_list()
{
  samweb list-files defname: $realdef | while read f; do 
    echo $f|cut -d_ -f2-3|sed -e's/r000//' -e's/_s0/ /' -e's/_s/ /'|while read run sr; do
      if ! ls $outhistdir/$rfctime-$stream/*det_r*${run}_*${sr}*_data.hists.root \
           &> /dev/null;then
        echo $f
      fi
    done
  done > $TMP

  if [ $(cat $TMP | wc -l) -eq 0 ]; then
    echo No files need to be redone for $unixtime $stream, exiting
    rm -f $TMP
    exit 0
  elif [ $iteration -gt 1 ]; then
    echo $(cat $TMP | wc -l) files need to be redone:
    cat $TMP
    echo
  fi
  let iteration++
}

resubmitdelay=1

do_a_redo()
{
  def="straitm_$(date +%Y%m%d)_redo_$realdef"
  samweb delete-definition $def 2> /dev/null
  dimensions="$(for f in $(cat $TMP); do
    printf "%s %s or " file_name $(basename $f);
  done | sed 's/or $//')"

  if ! [ "$dimensions" ]; then
    echo Uh oh, I got an empty dimensions statement from this file:
    cat $TMP
    rm -f $TMP
    exit 1
  fi
  samweb create-definition $def "$dimensions"

  N=$(cat $TMP | wc -l)

  rm -f $TMP

  testrel=/nova/app/users/mstrait/novasoft-ligo/
  $testrel/ligo/stage.sh $def
  $testrel/ligo/submit.sh $unixtime $stream "$skymap" $def

  echo Now will watch $rfctime $stream

  # Wait and continue waiting longer between resubmit attempts
  # so that if there is a systemic problem, we don't hammer it
  sleep ${resubmitdelay}m
  let resubmitdelay++

  while true; do
    sleep 3m
    # retry loop for jobsub_q, which is flaky
    while ! jobsub_q --user mstrait > /tmp/joblist.$$; do
      echo jobsub_q failed.  Will try again.
      sleep 15 
    done

    cat /tmp/joblist.$$|grep $def|tee $TMP|grep ' H '|awk '{print $1}'|while read j; do
        echo There is a held $rfctime $stream job.  Killing it.
        jobsub_rm --jobid $j
        sleep 1
    done

    rm -f /tmp/joblist.$$

    if [ $(cat $TMP | wc -l) -eq 0 ]; then
      echo No $rfctime $stream jobs left, stopping watch
      break
    else
      echo At $(date):
      echo $(cat $TMP | grep ' R ' | wc -l) $rfctime $stream 'job(s)' running \
           $(cat $TMP | grep ' I ' | wc -l) idle
      echo
    fi
    rm -f $TMP
  done
}

while true; do
  find_redo_list # exits if nothing to redo
  do_a_redo
done
