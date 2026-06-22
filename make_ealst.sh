#!/bin/bash

activeAt=$1
nat=$2

startt=$((activeAt+1))
cou=$((nat-activeAt))

rm -rf EA.LST

echo "$cou" >> EA.LST
for i in $(seq $startt 1 $nat) ; do 
	echo "$i" >> EA.LST
done


