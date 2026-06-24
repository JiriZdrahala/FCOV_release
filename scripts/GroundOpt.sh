#!/bin/bash
set -e

if [ $# -lt 4 ] ; then
	echo 'USAGE: '
	head -15 ${0}
	exit 1
fi

geometry=`cat`
mem=$1
nproc=$2
methodANDbasis=$3
name=$4
comment=$5


cat << EOF
%chk=${name}.chk
%mem=${mem}
%nproc=${nproc}
# ${methodANDbasis} iop(2/11=1) stable=opt scf(xqc) scrf(cpcm,solvent=water)

$comment

$geometry


--link1--
%chk=${name}.chk
%mem=${mem}
%nproc=${nproc}
# ${methodANDbasis} iop(2/11=1) geom=AllCheck Guess=Read opt scf(xqc) scrf=(cpcm,solvent=water)


EOF

