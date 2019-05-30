#!/bin/bash

. $SRT_PRIVATE_CONTEXT/ligo/env.sh

dir="$(basename "$1")"
def=$($SRT_PRIVATE_CONTEXT/ligo/dir2def.sh "$dir")

LIST=/tmp/iscomplete.$$

samweb list-files defname: $def > $LIST

if [ $(cat $LIST | wc -l) -eq 0 ]; then
  echo Empty SAM definition. Cannot tell if this is complete.
  exit 1
fi

missing=$(for f in $(cat $LIST); do
  echo $f|cut -d_ -f2-3|sed -e's/r000//' -e's/_s/ /'| \
    while read run sr; do
    if ! ls $outhistdir/$dir/*det_r*${run}_s${sr}_*_data*.hists.root \
       &> /dev/null;then
      echo hist file from $f missing | tee /dev/stderr
    fi
  done
done | grep -c missing)

rm -f $LIST

if [ $missing -eq 0 ]; then
  echo $dir complete
  exit 0
fi

echo $missing files missing from $dir
exit $missing
