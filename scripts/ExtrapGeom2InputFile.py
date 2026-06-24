#!/usr/bin/env python

import re
import sys

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def isMultiplicity(input_text):
    pattern = re.compile(r"-?[0-9]+ [0-9]+", re.IGNORECASE)
    return pattern.match(input_text)

def isAtomDefinition(input_text):
    new_text=input_text.split()
    if(len(new_text)==0): return False
    if(len(new_text)==4 and is_number(new_text[1]) and is_number(new_text[2]) and is_number(new_text[3])):
        return True
    else:
        return False

def isAtomDefinition_regex(input_text):
    pattern = re.compile(r"[A-Za-z0-9]+\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?\s+([+-]?(?=\.\d|\d)(?:\d+)?(?:\.?\d*))(?:[Ee]([+-]?\d+))?", re.IGNORECASE)
    return pattern.match(input_text)

EGFile=sys.argv[1]
InputFile_loc=sys.argv[2]
EGFile=open(EGFile,'r')
InputFile=open(InputFile_loc,'r')

EG=EGFile.readlines()
EGFile.close()

lines=InputFile.readlines()
InputFile.close()
success=False
GeomStart=-1
GeomEnd=-1
for idx,line in enumerate(lines):
    if(isMultiplicity(line)):
        GeomStart=idx+1
        success=True
        continue
    if(success and isAtomDefinition(line)):
        continue
    if(success and not isAtomDefinition(line)):
        GeomEnd=idx-1
        break
        
if(not success or GeomStart==-1 or GeomEnd==-1):
    print("FAIL 1")
    exit(1)
#linesAfter=lines[GeomEnd+1:]

for i in range(0,GeomEnd-GeomStart+1):
    splitt=lines[GeomStart+i].split()
    splitt[1]=EG[i]
    splitt.pop(2)
    splitt.pop(2)
    lines[GeomStart+i]=' '.join(splitt)
    
InputFile=open(InputFile_loc,'w')
InputFile.writelines(lines)

InputFile.close()