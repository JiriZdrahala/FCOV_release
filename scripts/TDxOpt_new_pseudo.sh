#!/bin/bash
set -e

file=$1
mem=$2
nproc=$3
method=$4
state=$5
nstates=$6
inputFile=$7
comment=$8


if [ $# -lt 8 ] ; then
	echo 'USAGE: '
	head -13 ${0}
	exit -1
fi

name=$(basename $file)

grep -i 'pseudo=read' $inputFile > /dev/null
if [ $? -ne 0 ] ; then
	echo 'No pseudo=read found'
	exit -2
fi

emptyLines=$(grep -n '^$' $inputFile | tr -s ' ' | cut -d':' -f1)
GENStart=$(echo $emptyLines | cut -d" " -f3)
#GENStart=$(($GENStart + 1))
#echo $GENStart
GEN=$(awk "NR > $GENStart && NR < $GENStart + 7" "$inputFile")
PSEUDO=$(awk "NR > $GENStart + 7 && NR < $GENStart + 11" "$inputFile")


if [ -z "$comment" ] ; then
	comment=${file%.*}
fi

cat << EOF
%chk=${name%.*}.chk
%mem=${mem}
%nproc=${nproc}
#p ${method}/GEN iop(2/11=1) td(root=${state},nstates=${nstates}) opt(CalcFC,tight,maxstep=20) scf(xqc) nosymm scrf=(cpcm,solvent=water) pseudo=read

$comment

$(~/scripts/bash/gau_lastCharge.sh $file) $(~/scripts/bash/gau_lastMultiplicity.sh $file)
$(~/scripts/bash/gau_lastInputOrientation.sh $file)

$(echo "$GEN")

$(echo "$PSEUDO")


--link1--
EOF

