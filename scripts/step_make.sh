#!/bin/bash
set -e

echo "${@}"

if [ "$#" -lt 3 ]; then
	echo "ARGUMENTS:"
	echo "1st - Gaussian p-TD-freq input file"
	echo "2nd - FIX option in DUSCH.OPT"
	echo "3rd to last - geometry step factor"
	exit 1
fi


backupMade=0
if [ -f DUSCH.OPT ] ; then
	mv DUSCH.OPT DUSCH.OPT_backup
	backupMade=1
fi

scriptDir=$(dirname $(realpath "$0"))
inputFile=$1
ifx=$2
#for start
for arg in ${@:3} ; do

argg=${arg//./}
echo $argg

cat << EOF > DUSCH.OPT
FIX
${ifx}
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
1
EXP_GEOM_FAC
${arg}

EOF

dusch_rev1

newDir=step/$argg
EGFile=EXTRAP_GEOM_${argg}
mkdir -p $newDir
cp $inputFile $newDir
mv EXTRAP_GEOM $newDir/$EGFile
cd $newDir
$scriptDir/ExtrapGeom2InputFile.py $EGFile $inputFile
#~/subg16-aurum-template-own $inputFile 4:00:00 a36_any 80GB
cd -
rm DUSCH.OPT

done
#for end

if [ $backupMade -gt 0 ] ; then
	mv DUSCH.OPT_backup DUSCH.OPT
fi

