#!/bin/bash
set -e

groundDir=$1
outfile=$2
NQ=$3

if [ $# -lt 3 ] ; then
	head -6 $0
	exit 1
fi

cp -r $groundDir .

mkdir -p excited
cd excited
geten ../$outfile
gar9 ../$outfile
new1
new2 0 0 y
new4 y 1 $NQ
cd ..

cat << EOF > DUSCH.OPT
FIX
-1
COR
t
SWE
t
KFAC
1.0
K_DEL
0
WFIX
100
WLIM
200
SSLIM
0.1
VERT_H
f
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
1
EXP_GEOM_FAC
1.0
EOF

dusch_rev1

