#!/bin/bash

# Make a little spreadsheet of what's processed.

if [ $GWNAME ]; then
  :
else
  GWNAME=$(basename "$PWD")
  if [ $GWNAME == ${GWNAME##*-} ]; then
    echo Please set GWNAME
    exit 1
  fi
  export GWNAME=${GWNAME##*-}
  echo I got GWNAME = $GWNAME from the PWD.  I hope that is right.
fi

. env.sh

printf combining
combineall.sh 2> /dev/stdout | crap2dot
echo

cd $outhadddir

printf '            neardet-ddactivity1   fardet-ddsnews        fardet-t02\n'

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

  elif [ $month == Apr ] && [ $day -eq 7  ]; then return 0 # fardet and neardet
  elif [ $month == Apr ] && [ $day -eq 8  ]; then return 0 # fardet and neardet
  elif [ $month == Apr ] && [ $day -eq 9  ] && [ $det ==  fardet ]; then return 0
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

for month in Jan Feb Mar Apr; do
  for day in {01..31}; do 
    if ( [ $month == Feb ] && [ $day -gt 28 ] ) ||
       ( [ $month == Apr ] && [ $day -gt 30 ] ); then
      continue
    fi

    printf '2019-%s-%s: ' $month $day
    if [ $month == Jan ]; then
      monthnum=01
    elif [ $month == Feb ]; then
      monthnum=02
    elif [ $month == Mar ]; then
      monthnum=03
    elif [ $month == Apr ]; then
      monthnum=04
    else
      echo I only know about 4 months.  I am but a child.
      exit 1
    fi

    shortday=$(echo $day | sed s/^0//)

    for trig in neardet-ddactivity1 fardet-ddsnews fardet-t02; do
      if nogood; then
        printf '              -'
      elif [ -e 2019-${monthnum}-${day}T14:29:01.Z-${trig}.hadded.root ] ||
           [ -e 2019-${monthnum}-${day}T13:29:01.Z-${trig}.hadded.root ]; then
        printf '              Y'
      elif grep $GWNAME-$month-$shortday-.*-$trig -q $TMPLOOP; then
        printf '              p'
      elif running; then
        printf '              s'
      else
        printf '               '
        if [ $trig == neardet-ddactivity1 ]; then
          incompletend=1
        elif [ $trig == fardet-ddsnews ]; then
          incompletesnews=1
        elif [ $trig == fardet-t02 ]; then
          incompletepulser=1
        fi
      fi
    done
    printf '\n'
  done
  printf -- '--------------------------------------------------------------\n'
done

if ! [ $incompletend ]; then
  echo COMPLETE NEARDET-DDACTIVITY1
fi
if ! [ $incompletesnews ]; then
  echo COMPLETE FARDET-DDSNEWS
fi
if ! [ $incompletepulser ]; then
  echo COMPLETE FARDET-T02
fi

rm -f $TMPJOB $TMPLOOP
