#!/bin/bash
set -e

file=$1
mem=$2
nproc=$3
methodANDbasis=$4
state=$5
nstates=$6
comment=$7


if [ $# -lt 6 ] ; then
	echo 'USAGE: '
	head -10 ${0}
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
#p ${methodANDbasis} iop(2/11=1) td(root=${state},nstates=${nstates}) freq(SaveNM) scf(xqc) scrf=(cpcm,solvent=water)

$comment

$(gau_lastCharge.sh $file) $(gau_lastMultiplicity.sh $file)
$(gau_lastInputOrientation.sh $file)



EOF

