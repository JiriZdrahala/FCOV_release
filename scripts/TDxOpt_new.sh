#!/bin/bash
set -e

file=$1
mem=$2
nproc=$3
methodANDbasis=$4
state=$5
nstates=$6
comment=$7


if [ $# -lt 7 ] ; then
	echo 'USAGE: '
	head -13 ${0}
	exit -1
fi

name=$(basename $file)

if [ -z "$scrf" ] ; then
	scrf='pcm'
fi

if [ -z "$td" ] ; then
	td='td'
fi

if [ -z "$comment" ] ; then
	comment=${file%.*}
fi


cat << EOF
%chk=${name%.*}.chk
%mem=${mem}
%nproc=${nproc}
#p ${methodANDbasis} iop(2/11=1) td(root=${state},nstates=${nstates}) opt(CalcFC,tight,maxstep=20) scf(xqc) nosymm scrf=(cpcm,solvent=water)

$comment

$(~/scripts/bash/gau_lastCharge.sh $file) $(~/scripts/bash/gau_lastMultiplicity.sh $file)
$(~/scripts/bash/gau_lastInputOrientation.sh $file)




EOF

