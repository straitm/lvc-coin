#!/bin/bash

dir="$(basename "$1")"
rfctime=$(echo $dir | cut -d- -f 1-3)
unixtime=$(date +%s -d "$(echo $rfctime | sed -e 's/T/ /' -e 's/\.//')")
stream=$(echo $dir | cut -d- -f 4-5)

def=strait-ligo-coincidence-artdaq-${unixtime}-${stream}

basedir=/pnfs/nova/scratch/users/mstrait/ligo
status=$(samweb list-files defname: $def | while read f; do 
  echo $f|cut -d_ -f2-3|sed -e's/r000//' -e's/_s0/ /' -e's/_s/ /'| \
    while read run sr; do
    if ! ls $basedir/$dir/*det_r*${run}_*${sr}*_data.hists.root \
       &> /dev/null;then
      echo hist file from $f missing
      let status++
    fi
  done
done | tee /dev/stderr | grep -c missing)

if [ $status -eq 0 ]; then
  echo $dir complete
  exit 0
fi

echo $status files missing
exit $status
