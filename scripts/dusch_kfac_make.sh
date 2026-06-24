#!/bin/bash
set -e


outfile=$1

if [ $# -lt 2 ] ; then
	echo "1st arg - Gaussian outfile"
	echo "2nd to last args - KFAC value"
	exit 1
fi

for arg in ${@:2} ; do

argg=${arg//./}
mkdir -p $argg
cp -r ground/ $argg/
cp -r excited/ $argg/
cp $outfile $argg/
if [ -f FCOV_V.OPT ] ; then
	cp FCOV_V.OPT $argg/
elif [ -f ${outfile%.*}.OPT ] ; then
	cp ${outfile%.*}.OPT $argg/
fi

cd $argg


cat << EOF > DUSCH.OPT
FIX
-1
COR
t
SWE
t
KFAC
${arg}
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

dusch_rev1 > /dev/null
echo "KFAC = $arg, <0|0> = $(cat OVERLAP)"
cd ..

done


