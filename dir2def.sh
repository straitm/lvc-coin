#!/bin/bash

# Convert the output directory name, in my convention, to the SAM definition
# used for it (also my convention)

dir="$(basename "$1")"
rfctime=$(echo $dir | cut -d- -f 1-3)
unixtime=$(date +%s -d "$(echo $rfctime | sed -e 's/T/ /' -e 's/\.//')")
stream=$(echo $dir | cut -d- -f 4-5)

echo strait-ligo-coincidence-artdaq-${unixtime}-${stream}
