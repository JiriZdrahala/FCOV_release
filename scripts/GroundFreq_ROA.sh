#!/bin/bash
set -e

file=$1
mem=$2
nproc=$3
methodANDbasis=$4
comment=$5
excFreq=$6

if [ $# -lt 4 ] ; then
	echo 'USAGE: '
	head -9 ${0}
	exit -1
fi

name=$(basename $file)


if [ -z "$excFreq" ] ; then
	excFreq="532nm"
fi	

if [ -z "$comment" ] ; then
	comment=${file%.*}
fi

cat << EOF
%chk=${name%.*}.chk
%mem=${mem}
%nproc=${nproc}
#p ${methodANDbasis} iop(2/11=1) freq=(ROA) CPHF=InputFreq scf(xqc) scrf=(cpcm,solvent=water)

$comment

$(~/scripts/bash/gau_lastCharge.sh $file) $(~/scripts/bash/gau_lastMultiplicity.sh $file)
$(~/scripts/bash/gau_lastInputOrientation.sh $file)

${excFreq}


EOF

