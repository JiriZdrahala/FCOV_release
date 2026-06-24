#!/bin/bash
set -e

if [ $# -lt 2 ] ; then
	echo "1st: OUTFILE"
	echo "2nd: Number of modes wanted in new4"
	exit 1
fi

outfile=$1
NQ=$2

mkdir -p excited
cd excited

geten ../$outfile
gar9 ../$outfile
new1
new2 0 0 y
new4 y 1 $NQ

cd ..


