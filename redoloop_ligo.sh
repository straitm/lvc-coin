#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

# SAM definition that we need to do and redo
realdef=$1
unixtime=$(echo $realdef | cut -d- -f 5)
fracsec=$(cut -d. -f 2 -s <<< $unixtime)
rfctime=$(TZ=UTC date "+%Y-%m-%dT%H:%M:%S" -d @$unixtime).${fracsec}Z
rfctimesafeforsam=${rfctime//:/-}
stream=$(echo $realdef | cut -d- -f 6-7)

if ! [ $stream ]; then
  echo Got a null stream from \"$realdef\", cannot continue
  exit 1
fi

# If empty, the right thing happens (which is that we use a dummy skymap)
skymap=$2

TMP=/tmp/mstrait.redolist.$$

iteration=0

vsleep()
{
  if ! [ "$2" == really ] && [ $REDOFAST ]; then
    echo Not sleeping since you told me to go fast
    return
  fi

  add=$((RANDOM%10))
  echo Sleeping $1 and $add seconds
  sleep $1
  sleep $add
}

makejoblist()
{
  retrydelay=45
  while ! jobsub_q --user mstrait > /tmp/joblist.$$; do
    echo jobsub_q failed.  Will try again.
    vsleep $retrydelay really
    if [ $retrydelay -lt 3600 ]; then
      let retrydelay*=2
    fi
  done
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
    if [ $n -le 2 ]; then
      break
    fi
    echo Waiting for $n sam list processes to finish
    vsleep $((try + 60 + RANDOM%60)) really
    let try++
  done
  if [ $try -gt 0 ]; then
    echo Ok, going ahead
  fi
}

find_redo_list()
{
  # This loop allows starting up a new copy of this script if the old copy
  # got killed. It should even protect against multiple running copies all
  # trying to get the same jobs run.
  while true; do
    makejoblist
    deffrag="_redo_$realdef"
    if grep ${GWNAME}.*$deffrag -q /tmp/joblist.$$; then
      echo Jobs are running already/still for this definition.
      vsleep 4m really
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
      basename $f | cut -d_ -f2-3 | sed -e's/r000//' -e's/_s/ /' | \
      while read run sr; do
        if ! ls $outhistdir/$rfctimesafeforsam-$stream/*det_r*${run}_*s${sr}*_data*.hists.root \
             &> /dev/null; then
          echo $f
        fi
      done
    done > $TMP
  fi

  if [ $(cat $TMP | wc -l) -eq 0 ]; then
    echo No files need to be redone for $unixtime $stream, exiting
    rm -f $TMP
    exit 0
  elif [ $iteration -gt 0 ]; then
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
  N=$(cat $TMP | wc -l)
  NTOT=$(samweb list-files defname: $realdef | wc -l)

  printf 'Must %sdo %d out of %d file%s\n' \
    "$(if [ $iteration -gt 1 ]; then printf re; fi)" $N $NTOT \
    "$(if [ $NTOT -ne 1 ]; then printf s; fi)"

  if [ $N -eq $NTOT ]; then
    # Avoid super awkward limitations of SAM if we need to process the whole set
    def=$realdef
  else
    def="redo_$(date +%Y%m%d)_$realdef" # defs cannot start with a number!
    samweb delete-definition $def 2> /dev/null

    # Desperate measures. The same file apparently cannot be in two "datasets",
    # as made with sam_add_dataset, which lets you give a list of files in a
    # straightforward way.  "samweb create-definition" can only take a list of
    # files on the command line with lots of "or file_name"s, and has a character
    # limit.  I bet there is a correct way to do this, but it's certainly not
    # easily discoverable.
    for n in `seq $N`; do
      echo Trying to make the definition with $n 'file(s)'
      dimensions="$(for f in $(head -n $n $TMP); do
        printf "file_name %s or " $(basename $f);
      done | sed 's/ or $//')"

      samweb delete-definition $def
      if ! samweb create-definition $def "$dimensions"; then
        if [ $n -le 1 ]; then
          echo Failed with only $n file. Is this because a file name starts
          echo with all numbers? SAM seems to choke on that. Quitting.
          exit 1
        fi
        echo Failed.  Making it with $((n-1)) 'file(s)' and will process rest later
        dimensions="$(for f in $(head -n $((n-1)) $TMP); do
          printf "file_name %s or " $(basename $f);
        done | sed 's/ or $//')"
        samweb create-definition $def "$dimensions"
        break
      else
        echo $n worked
      fi
    done
  fi

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
    vsleep 3m really
    makejoblist

    cat /tmp/joblist.$$|grep $GWNAME-$def|tee $TMP|grep ' H '|awk '{print $1}'|\
    while read j; do
        echo There is a held $rfctime $stream $GWNAME job.  Killing it.
        jobsub_rm --jobid $j
        vsleep 1
    done

    rm -f /tmp/joblist.$$

    if [ $(cat $TMP | wc -l) -eq 0 ]; then
      echo No $rfctime $stream $GWNAME jobs left, stopping watch
      break
    else
      echo At $(date):
      echo $(cat $TMP | grep ' R ' | wc -l) $rfctime $stream $GWNAME 'job(s)' running \
           $(cat $TMP | grep ' I ' | wc -l) idle
    fi
    rm -f $TMP
  done
}

while true; do
  find_redo_list # exits if nothing to redo
  do_a_redo
done
