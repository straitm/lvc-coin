#!/bin/bash

e=$1

set -o pipefail

url=https://gracedb.ligo.org/apiweb/superevents/$e/files/${e}-1-Preliminary.xml,0

if ! wget -q -O - $url | grep Date | cut -d\> -f 2 | cut -d\< -f 1; then
  wget $url -O /dev/null
fi
