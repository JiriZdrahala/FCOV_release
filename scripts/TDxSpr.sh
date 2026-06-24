#!/bin/bash
set -e


chkGround=$1
chkExcited=$2
nproc=$3
mem=$4
HWHM=$5


if [ $# -lt 5 ] ; then
	head -9 ${0}
	exit 1
fi

if [ -z $HWHM ] ; then
	HWHM=150
fi

cat << EOF
%mem=${mem}
%nproc=${nproc}
%chk=${chkGround}
#p geom=allcheck maxdisk=40GB Frequency=(ReadFC,FCHT,ReadFCHT)

Spectroscopy=CircularDichroism
Method=VerticalHessian
Spectrum=(HWHM=${HWHM},Grain=4)
Print=(matrix=JKABCDE)

${chkExcited}


EOF

