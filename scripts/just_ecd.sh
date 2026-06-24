#!/bin/bash

file=$1
mem=$2
nproc=$3
methodANDbasis=$4
comment=$5


if [ $# -lt 4 ] ; then
	echo 'USAGE: '
	head -7 ${0}
	exit 1
fi

name=${file##*/}
name=${name%.*}

cat << EOF
%chk=${name}.chk
%mem=${mem}
%nproc=${nproc}
# ${methodANDbasis} iop(2/11=1) td(nstates=100) scf(novaracc,xqc) scrf=(cpcm,solvent=water)

$comment

$(~/scripts/bash/gau_lastCharge.sh $file) $(~/scripts/bash/gau_lastMultiplicity.sh $file)
$(~/scripts/bash/gau_lastInputOrientation.sh $file)



EOF

