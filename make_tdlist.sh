#!/bin/bash

DIR=$1
STUFF=$2

if [ -z $DIR ] ; then
	DIR="./"
fi

lendir=${#DIR}
lenk=$((lendir+6))
rm LIST.STEPS.TD > /dev/null 2>&1

if [[ "$STUFF" == *"u"* ]] ; then
	ls ${DIR}FILE_*.OUT | sort -n -k1.${lenk} >> LIST.STEPS.TD
else
	ls ${DIR}FILE_*.out | sort -n -k1.${lenk} >> LIST.STEPS.TD
fi





