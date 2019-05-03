#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

# SAM definition that we need to do and redo
realdef=$1
unixtime=$(echo $realdef | cut -d- -f 5)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
stream=$(echo $realdef | cut -d- -f 6-7)

if ! [ $stream ]; then
  echo Got a null stream from \"$realdef\", cannot continue
  exit 1
fi

# If empty, the right thing happens (which is that we use a dummy skymap)
skymap=$2

TMP=/tmp/mstrait.redolist.$$

iteration=1

vsleep()
{
  echo Sleeping $1 or a little bit more
  sleep $1
  sleep $((RANDOM%10))
}

nsamlistsrunning()
{
  ps f -u mstrait | grep -v grep | grep samweb.py\ list | wc -l
}

# SAM won't let a user do more than 5 queries at a time.  Self-limit
# to three so that I can also do interactive queries while my scripts run.
blocksam()
{
  local try=0
  while true; do
    n=$(nsamlistsrunning)
    if [ $n -le 2 ]; then break; fi
    echo Waiting for $n sam list processes to finish
    sleep $((try + 60 + RANDOM%60))
    let try++
  done
  echo Ok, going ahead
}

find_redo_list()
{
  # This loop allows starting up a new copy of this script if the old copy
  # got killed. It should even protect against multiple running copies all
  # trying to get the same jobs run.
  while true; do
    while ! jobsub_q --user mstrait > /tmp/joblist.$$; do
      echo jobsub_q failed.  Will try again.
      vsleep 45
    done

    deffrag="_redo_$realdef"
    if grep ${GWNAME}.*$deffrag /tmp/joblist.$$; then
      echo Jobs are running already/still for this definition.
      vsleep 4m
    else
      break
    fi
  done

  if [ $stream == neardet-t00 ]; then
    if ! [ $jobid ]; then
      blocksam
      samweb list-files defname: $realdef > $TMP
    else
      jobsub_fetchlog --jobid $jobid
      ngood=$(tar xzf $jobid.tgz -O | \
              grep -c '^Art has completed and will exit with status 0\.$')
      blocksam
      nneed=$(samweb list-files defname: $realdef | wc -l)
      if [ $ngood -eq $nneed ]; then
        tar xzf $jobid.tgz -O | grep '^Spilltime: ' > spills-$unixtime-$rfctime.txt
      else
        blocksam
        samweb list-files defname: $realdef > $TMP
      fi
    fi
  else
    blocksam
    samweb list-files defname: $realdef | while read f; do
      echo $f | cut -d_ -f2-3 | sed -e's/r000//' -e's/_s0/ /' -e's/_s/ /' | \
      while read run sr; do
        if ! ls $outhistdir/$rfctime-$stream/*det_r*${run}_*${sr}*_data.hists.root \
             &> /dev/null;then
          echo $f
        fi
      done
    done > $TMP
  fi

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
  if [ $iteration -gt 10 ]; then
    echo Tried 10 times and still failed, giving up
    exit 1
  fi
}

resubmitdelay=$((60+RANDOM%300))

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

  # Wait and continue waiting longer between resubmit attempts
  # so that if there is a systemic problem, we don't hammer it
  vsleep ${resubmitdelay}
  if [ $resubmitdelay -lt 3600 ]; then
    let resubmitdelay*=2
  fi

  testrel=/nova/app/users/mstrait/novasoft-ligo/
  $testrel/ligo/stage.sh $def
  if [ $stream == neardet-t00 ]; then
    jobid=$($testrel/ligo/spillsubmit.sh $unixtime $def | tee /dev/stderr | \
      grep 'JobsubJobId of first job:' | awk '{print $5}')
  else
    $testrel/ligo/submit.sh $unixtime $stream $def
  fi

  echo Now will watch $rfctime $stream

  while true; do
    vsleep 3m
    # retry loop for jobsub_q, which is flaky
    while ! jobsub_q --user mstrait > /tmp/joblist.$$; do
      echo jobsub_q failed.  Will try again.
      vsleep 45
    done

    cat /tmp/joblist.$$|grep $GWNAME-$def|tee $TMP|grep ' H '|awk '{print $1}'|\
    while read j; do
        echo There is a held $rfctime $stream job.  Killing it.
        jobsub_rm --jobid $j
        vsleep 1
    done

    rm -f /tmp/joblist.$$

    if [ $(cat $TMP | wc -l) -eq 0 ]; then
      echo No $rfctime $stream jobs left, stopping watch
      break
    else
      echo At $(date):
      echo $(cat $TMP | grep ' R ' | wc -l) $rfctime $stream 'job(s)' running \
           $(cat $TMP | grep ' I ' | wc -l) idle
    fi
    rm -f $TMP
  done
}

while true; do
  find_redo_list # exits if nothing to redo
  do_a_redo
done
