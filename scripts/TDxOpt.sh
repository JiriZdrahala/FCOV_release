#!/bin/bash
set -e

file=$1
mem=$2
nproc=$3
methodANDbasis_lower=$4
methodANDbasis=$5
state=$6
nstates=$7
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
# ${methodANDbasis_lower} iop(2/11=1) td(root=${state},nstates=${nstates}) nosymm opt(loose) scf(xqc) scrf=(cpcm,solvent=water)

$comment

$(~/scripts/bash/gau_lastCharge.sh $file) $(~/scripts/bash/gau_lastMultiplicity.sh $file)
$(~/scripts/bash/gau_lastInputOrientation.sh $file)


--link1--
%chk=${name%.*}.chk
%mem=${mem}
%nproc=${nproc}
# ${methodANDbasis} iop(2/11=1) td(root=${state},nstates=${nstates}) geom=allcheck guess=read nosymm opt(maxstep=19) scf(xqc) scrf=(cpcm,solvent=water)




EOF

