#!/bin/bash

# This gets the time the trigger was *issued*, not the time of the GW event

e=$1

set -o pipefail

url=https://gracedb.ligo.org/apiweb/superevents/$e/files/${e}-1-Preliminary.xml,0

if ! wget --no-check-certificate -q -O - $url | grep Date | cut -d\> -f 2 | cut -d\< -f 1; then
  wget $url -O /dev/null
fi
