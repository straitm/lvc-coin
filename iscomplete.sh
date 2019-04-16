#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

dir="$(basename "$1")"
def=$($SRT_PRIVATE_CONTEXT/ligo/dir2def.sh "$dir")

status=$(samweb list-files defname: $def | while read f; do 
  echo $f|cut -d_ -f2-3|sed -e's/r000//' -e's/_s0/ /' -e's/_s/ /'| \
    while read run sr; do
    if ! ls $outhistdir/$dir/*det_r*${run}_*${sr}*_data.hists.root \
       &> /dev/null;then
      echo hist file from $f missing
    fi
  done
done | grep -c missing)

if [ $status -eq 0 ]; then
  echo $dir complete
  exit 0
fi

echo $status files missing from $dir
exit $status
