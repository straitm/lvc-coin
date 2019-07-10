#!/bin/bash

# Make a little spreadsheet of what's processed.

if [ $GWNAME ]; then
  :
else
  GWNAME=$(basename "$PWD")
  if [ $GWNAME == ${GWNAME##*-} ]; then
    echo Please set GWNAME > /dev/stderr
    exit 1
  fi
  export GWNAME=${GWNAME##*-}
  echo I got GWNAME = $GWNAME from the PWD.  I hope that is right. > /dev/stderr
fi

. env.sh

combine=false
commands=false
remove=false

while [ $1 ]; do 
  if [ "$1" == combine ]; then
    combine=true
  fi
  if [ "$1" == commands ]; then
    commands=true
  fi
  if [ "$1" == remove ]; then
    remove=true
  fi
  shift
done

if $combine; then
  printf combining
  combineall.sh 2> /dev/stdout | crap2dot
  echo
else
  ! $commands && echo Not combining.
fi

cd $outhadddir

! $commands && printf '            neardet-ddactivity1   fardet-ddsnews        fardet-t02\n'

nogood()
{
  det=$(echo $trig | cut -d- -f 1)

  # Avoid irritating special case for the first entry below
  if false; then nasal demons

  # t02 has a missing portion outside the ddsnews window, and ddsnews also does,
  # oddly, at a different time.
  elif [ $month == Jan ] && [ $day -eq 1  ] && [ $det ==  fardet ]; then return 0

  elif [ $month == Jan ] && [ $day -eq 11 ] && [ $det == neardet ]; then return 0
  elif [ $month == Jan ] && [ $day -eq 23 ] && [ $det ==  fardet ]; then return 0
  elif [ $month == Jan ] && [ $day -eq 30 ] && [ $det == neardet ]; then return 0
  elif [ $month == Jan ] && [ $day -eq 31 ] && [ $det == neardet ]; then return 0

  # Only one file selected for t02 and ddsnews each
  elif [ $month == Feb ] && [ $day -eq 1  ] && [ $det ==  fardet ]; then return 0

  elif [ $month == Feb ] && [ $day -eq 20 ]; then return 0 #fardet and neardet

  # t02 has some data, but ddsnews doesn't.  Probably t02 therefore bad.
  elif [ $month == Mar ] && [ $day -eq 20 ] && [ $det ==  fardet ]; then return 0

  # These three probably have good data and only appeared blank because of
  # grid problems.  It is expedient to exclude them and I see no way that
  # could bias the sample, since the reasons they failed aren't related
  # to the data in the files.
  elif [ $month == Apr ] && [ $day -eq 2  ] && [ $det ==  fardet ]; then return 0
  elif [ $month == Apr ] && [ $day -eq 4  ] && [ $det ==  fardet ]; then return 0
  elif [ $month == Apr ] && [ $day -eq 12 ] && [ $det ==  fardet ]; then return 0

  elif [ $month == Apr ] && [ $day -eq 7  ]; then return 0 # fardet and neardet
  elif [ $month == Apr ] && [ $day -eq 8  ]; then return 0 # fardet and neardet
  elif [ $month == Apr ] && [ $day -eq 9  ] && [ $det == fardet ]; then return 0
  fi
  return 1
}

TMPJOB=/tmp/mstrait.joblist.$$
TMPLOOP=/tmp/mstrait.looplist.$$

jobsub_q --user mstrait > $TMPJOB

(ps f; ssh novagpvm11 ps f) | tee | cut -d/ -f6- > $TMPLOOP

# Look for whether there's a job running even though there's no redoloop watching it.
# Assume bg event from 8:30.
running()
{
  unixtime=$(TZ=UTC date +%s -d"$month $shortday 8:29:01 2019")

  cat $TMPJOB | grep -q $GWNAME-strait-ligo-coincidence-*-$unixtime-$trig
}

completend=0
completesnews=0
completepulser=0

for month in Jan Feb Mar Apr; do
  for day in {01..31}; do 
    if ( [ $month == Feb ] && [ $day -gt 28 ] ) ||
       ( [ $month == Apr ] && [ $day -gt 30 ] ); then
      continue
    fi

    ! $commands && printf '2019-%s-%s: ' $month $day
    if [ $month == Jan ]; then
      monthnum=01
    elif [ $month == Feb ]; then
      monthnum=02
    elif [ $month == Mar ]; then
      monthnum=03
    elif [ $month == Apr ]; then
      monthnum=04
    else
      echo I only know about 4 months.  I am but a child. >/dev/stderr
      exit 1
    fi

    shortday=$(echo $day | sed s/^0//)

    for trig in neardet-ddactivity1 fardet-ddsnews fardet-t02; do
      if ! nogood; then
        if [ $trig == neardet-ddactivity1 ]; then
        let totalnd++
        elif [ $trig == fardet-ddsnews ]; then
          let totalsnews++
        elif [ $trig == fardet-t02 ]; then
          let totalpulser++
        fi
      fi

      if nogood; then
        ! $commands && printf '              -'
      elif [ -e 2019-${monthnum}-${day}T14:29:01.Z-${trig}.hadded.root ] ||
           [ -e 2019-${monthnum}-${day}T13:29:01.Z-${trig}.hadded.root ]; then
        ! $commands && printf '              Y'
        if [ $trig == neardet-ddactivity1 ]; then
          let completend++
        elif [ $trig == fardet-ddsnews ]; then
          let completesnews++
        elif [ $trig == fardet-t02 ]; then
          let completepulser++
        fi
      elif grep $GWNAME-$month-$shortday-.*-$trig -q $TMPLOOP; then
        ! $commands && printf '              p'
      elif running; then
        ! $commands && printf '              s'
      else
        if $commands; then
          printf "throttle; echo $month $shortday $trig $GWNAME;\n"
          printf "getandredoloop.sh $month $shortday $trig 2019 $GWNAME &> /dev/null & sleep 1\n"
        else
          printf '               '
        fi
      fi
    done
    ! $commands && printf '\n'
  done
  ! $commands && printf -- '--------------------------------------------------------------\n'
done

for trig in neardet-ddactivity1 fardet-ddsnews fardet-t02; do
  for file in *${trig}.{pdf,hadded.root}; do
    # takes care of unexpanded stars
    if ! [ -e $file ]; then continue; fi

    # ignore the 114x-neard-ddactivity1.hadded.root type of files
    if echo ${file:0:5} | grep x -q; then continue; fi

    monthnum=${file:5:2}
    day=${file:8:2}
    if [ $monthnum == 01 ]; then
      month=Jan
    elif [ $monthnum == 02 ]; then
      month=Feb
    elif [ $monthnum == 03 ]; then
      month=Mar
    elif [ $monthnum == 04 ]; then
      month=Apr
    else
      echo I only know about 4 months.  This file is suspicious: > /dev/stderr
      echo $file > /dev/stderr
      if $remove; then rm -v $file; fi
      continue
    fi

    if nogood; then
      echo Extra file $file which is no good > /dev/stderr
      if $remove; then rm -v $file; fi
    fi
  done
done

div()
{
  printf '%5.1f%% (%d/%d)' $(printf "1k $1 100*$2/p"|dc) $1 $2
}

! $commands && echo $GWNAME $(div $completend $totalnd) nd-act \
             $(div $completesnews $totalsnews) fd-ddsnews \
             $(div $completepulser $totalpulser) t02

rm -f $TMPJOB $TMPLOOP
