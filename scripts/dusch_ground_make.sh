#!/bin/bash
set -e

if [ $# -lt 2 ] ; then
	echo "1st: OUTFILE"
	echo "2nd: Number of modes wanted in new4"
	echo "3rd (optional): FREQ OUTFILE"
	exit 1
fi

outfile=$1
NQ=$2
ZPEfile=$3

mkdir -p ground
cd ground

if [ -z "$ZPEfile" ] ; then
	geten ../$outfile
	gar9 ../$outfile
else
	geten_ZPE ../$outfile ../$ZPEfile
	gar9 ../$ZPEfile
fi

new1
new2 0 0 y
new4 y 1 $NQ
cd ..


