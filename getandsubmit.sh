#!/bin/bash

md="$1"

./getfilelist.sh $(date +%s -d"$md 8:29:01 2019")
./submit.sh $(date +%s -d"$md 8:29:01 2019") fardet-ddsnews
./submit.sh $(date +%s -d"$md 8:29:01 2019") fardet-ddenergy
./submit.sh $(date +%s -d"$md 8:29:01 2019") fardet-t02
