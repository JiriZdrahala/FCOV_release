#!/bin/bash
set -e



if [ $# -lt 2 ] ; then
	echo "1st arg - Gaussian outfile"
	echo "2nd arg - KFAC value"
	echo "3rd arg - DLT value"
	echo "4th arg - IFX value"
	echo "5th arg - REG? [t/f]"
	echo "6th arg - LSFIX? [t/f]"
	exit 1
fi


outfile=$1

kfac=$2
if [ -z $kfac ] ; then
	kfac=1
fi

dlt=$3
if [ -z $dlt ] ; then
	dlt="0.0"
fi

ifx=$4
if [ -z $ifx ] ; then
	ifx="-1"
fi

reg=$5
reg_str="_REG"
if [ -z $reg ] || [ "$reg" == "f" ] ; then
	reg="f"
	reg_str=""
fi


lsfix=$6
lsfix_str="_LSFIXf"
if [ -z $lsfix ] || [ "$lsfix" == "t" ] ; then
	lsfix="t"
	lsfix_str=""
fi

if [ "$ifx" -eq 2 ] ; then
	if (( $(echo "$kfac == 1" | bc -l) )) ; then
		dir="VH${ifx}_${dlt}${lsfix_str}${reg_str}"
	else
		dir="VHs_${kfac}_${ifx}_${dlt}${lsfix_str}${reg_str}"
	fi
else
	if (( $(echo "$kfac == 1" | bc -l) )) ; then
		dir="VH${ifx}"
	else
		dir="VHs${kfac}_${ifx}"
	fi
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
f
COR
t
SWE
t
KFAC
$kfac
DLT
$dlt
REG
$reg
LSFIX
$lsfix
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
f
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



