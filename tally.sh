#!/bin/bash

# Make a little spreadsheet of what's processed.

if ! [ $GWNAME ]; then
  echo Please set GWNAME
  exit 1
fi

. env.sh

combineall.sh

cd $outhadddir

printf '            neardet-ddactivity1   fardet-ddsnews        fardet-t02\n'

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

    for trig in neardet-ddactivity1 fardet-ddsnews fardet-t02; do
      if [ -e 2019-${monthnum}-${day}T14:29:01.Z-${trig}.hadded.root ] ||
         [ -e 2019-${monthnum}-${day}T13:29:01.Z-${trig}.hadded.root ]; then
        printf '        Y'
      else
        printf '        .'
      fi
    done
    printf '\n'
  done
  printf '\n'
done
