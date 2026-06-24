#!/bin/bash
set -e

nameChk=$1
mem=$2
nproc=$3


if [ $# -lt 3 ] ; then
	echo 'USAGE: '
	head -7 ${0}
	exit -1
fi



cat << EOF
%chk=${nameChk}
%mem=${mem}
%nproc=${nproc}
#p iop(2/11=1) freq(ReadFc) Geom=AllCheck nosymm



EOF

