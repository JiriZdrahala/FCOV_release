#!/bin/bash

pathh=$0
numTrs=$1


input="FILE.TR"

i=0
counter=1
curDir=""
while IFS= read -r line ; do
    if [ $((i % numTrs)) -eq 0 ] ; then
        curDir="TR_$counter"
        counter=$((counter+1))
        if test -d $curDir ; then
            echo "$curDir exists. Save your work and remove it."
            echo "Exiting"
            exit 1
        fi
        mkdir $curDir
        echo "DIR: $curDir"
    fi
    echo $line
    echo $line >> $curDir/FILE.TR
    i=$((i+1))
    
done < $input