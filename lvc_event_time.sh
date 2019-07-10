#!/bin/bash

e=$1

wget -q -O - https://gracedb.ligo.org/apiweb/superevents/S190630ag/files/${e}-1-Preliminary.xml,0 | grep ISOTime | cut -d\> -f 2 | cut -d\< -f 1
