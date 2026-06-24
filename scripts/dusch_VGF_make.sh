#!/bin/bash
set -e



if [ $# -lt 2 ] ; then
	echo "1st arg - Gaussian outfile"
	echo "2nd arg - KFAC value"
	echo "3rd arg - IFX value"
	echo "4th arg - IFX6? [t/f]"
	exit 1
fi


outfile=$1

kfac=$2
if [ -z $kfac ] ; then
	kfac="1.0"
fi

ifx=$3
if [ -z $ifx ] ; then
	ifx="-1"
fi

ifx6=$4
if [ -z $ifx6 ] ; then
	ifx6="f"
fi

if (( $(echo "$kfac != 1" |bc -l) )) ; then
	dir="VGFs${kfac}_${ifx}"
else
	dir="VGF${ifx}"
fi

mkdir -p $dir
cp -r ground/ $dir
cp -r excited/ $dir
cp $outfile $dir
cd $dir


cat << EOF > DUSCH.OPT
FIX
$ifx
FIX6
$ifx6
COR
t
SWE
t
KFAC
$kfac
K_DEL
0
WFIX
100
WLIM
200
SSLIM
0.1
VERT_H
t
J_ONE
t
WE_IS_WG
f
NO_REORI
f
SHIFTQ
f
SHIFTI
f
SHIFTE
f
EXP_GEOM
0
EXP_GEOM_FAC
0.2
EOF

dusch_rev1
cd ..



