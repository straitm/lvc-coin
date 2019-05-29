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
      elif ps f | grep tee | cut -d/ -f6- | grep $GWNAME-$month-$shortday-.*-$trig -q; then
        printf '              p'
      else
        printf '               '
      fi
    done
    printf '\n'
  done
  printf -- '--------------------------------------------------------------\n'
done
