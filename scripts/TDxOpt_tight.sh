#!/bin/bash
set -e

file=$1
mem=$2
nproc=$3
methodANDbasis=$4
state=$5
nstates=$6
maxstep=$7
comment=$8


if [ $# -lt 7 ] ; then
	echo 'USAGE: '
	head -13 ${0}
	exit -1
fi

name=$(basename $file)


if [ -z "$comment" ] ; then
	comment=${file%.*}
fi


cat << EOF
%chk=${name%.*}.chk
%mem=${mem}
%nproc=${nproc}
#p ${methodANDbasis} iop(2/11=1) td(root=${state},nstates=${nstates}) opt(tight,maxstep=${maxstep}) scf(xqc) nosymm scrf=(cpcm,solvent=water)

$comment

$(~/scripts/bash/gau_lastCharge.sh $file) $(~/scripts/bash/gau_lastMultiplicity.sh $file)
$(~/scripts/bash/gau_lastInputOrientation.sh $file)



EOF

