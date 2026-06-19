#!/usr/bin/env python

false=False
true=True

import sys
import re
from itertools import combinations
from itertools import permutations
from itertools import product
import math

def MakePerms(arr:list,k:int)->list:
    return [p for p in product(arr, repeat=k)]

def IsTrAllowed(v1:list[int],v1_pos:list[int],v1_class:int,v3:list[int],v3_pos:list[int],v3_class:int)->bool:
    if(abs(v3_class-v1_class)>1):return False
    if(v1_class==0):
        if(v3_class>1):return false
        if(v3[0]>2):return false
        return true
    elif(v1_class==1):
        if(v3_class>2):return false
        idx=-1
        try:
            idx=v3_pos.index(v1_pos[0])
        except:
            idx=-1
        if(idx>-1):
            if(v3_class==1 and v1[0]>=v3[0]):return False
            if(v1[0]-v3[idx]==1 and not (2 in v3)):return true
            if(abs(v1[0]-v3[idx])==1 and v3_class==1):return true
            if(v3_class==2):
                if(v1[0]==1 and (2 in v3)):return false
                if(v1[0]==2 and (2 in v3)):return false
            #return false
        if(not ModeIsSubsetof(v1,v1_pos,v1_class,v3,v3_pos,v3_class)):return False
        return True
    elif(v1_class==2): #¯\_(ツ)_/¯
        return False
    else: #¯\_(ツ)_/¯
        return False

def ModeIsSubsetof(v1:list[int],v1_pos:list[int],v1_class:int,v3:list[int],v3_pos:list[int],v3_class:int)->bool:
    if(v1_class==v3_class):
        for i in range(v1_class):
            idx=-1
            try:
                idx=v3_pos.index(v1_pos[i])
            except:
                idx=-1
            if(idx>-1):
                if(v1[i]<=v3[idx]):continue
            return False
        return True
    elif(v1_class<v3_class):
        for i in range(v1_class):
            idx=-1
            try:
                idx=v3_pos.index(v1_pos[i])
            except:
                idx=-1
            if(idx>-1):
                if(v1[i]<=v3[idx]):continue
            return False
        return True
    elif(v1_class>v3_class):
        for i in range(v3_class):
            idx=-1
            try:
                idx=v1_pos.index(v3_pos[i])
            except:
                idx=-1
            if(idx>-1):
                if(v3[i]<=v1[idx]):continue
            return False
        return True

def Tr2Str(v1:list[int],v1_pos:list[int],v1_class:int,v3:list[int],v3_pos:list[int],v3_class:int)->str:
    strr=""
    if(v1_class==0):
        strr=strr+("0 ")
    else:
        for i in range(v1_class):
            strr=strr + "{0}^{1} ".format(v1_pos[i],v1[i])
    
    strr=strr+"> "
    
    if(v3_class==0):
        strr=strr+"0"
    else:
        for i in range(v3_class):
            strr=strr+"{0}^{1} ".format(v3_pos[i],v3[i])
    return strr

def Eq0(num:int)->bool:
    return num==0


def Read_FINP_Freqs(filepath:str) -> list[float]:
    finp=open(filepath,'r')
    #n3=3*nat
    #nq=n3-6
    line=finp.readline()
    linesplit=line.split()
    linesplit=list(map(int,linesplit))
    nq=linesplit[0]
    n3=linesplit[1]
    nat=linesplit[2]
    for i in range(nat+1):finp.readline()
    
    #L=np.zeros(shape=(n3,nq))
    
    for i in range(nat*nq):
        line=finp.readline()
        #linesplit=line.split()
        #i_r=(int(linesplit[0])-1)*3
        #i_m=int(linesplit[1])
        #vec=list(map(float,linesplit[2:]))
        #L[i_r:i_r+3,nq-i_m]=vec
    finp.readline()
    rows=math.ceil(nq/6.0)
    ws=[]
    for i in range(rows):
        el_L=i*6
       # el_R=min(el_L+6,nq)
        line=finp.readline()
        linesplit=line.split()
        linesplit=list(map(float,linesplit))
        
        ws.extend(linesplit)
    return ws


argc=len(sys.argv)-1

if(argc==0):
    print("USAGE:")
    print('--nq=[number of vib. modes]')
    print('--sel-rules     apply some rudimentary selection rules for the transitions')
    print('--max-v1s       maximum excitations for each class for initial state (v1)')
    print('--max-v3s       maximum excitations for each class for final state (v3)')
    print('--n-ground      number of ground states to include starting from fundamental excitations (different approach from FCOV since frequencies are not read currently)')
    exit(-1)
argv=sys.argv[1:]

nq=0
sel_rules=False
max_v1s=[0]
mc_v1=0
max_v3s=[1]
mc_v3=1
n_ground=1
bf_thr=0
sort=False
wg=[]
for i in range(argc):
    arg=argv[i]
    idx_eq=arg.find('=')+1
    if('--nq=' in arg):
        nq=int(arg[idx_eq:])
    elif('--sel-rules' in arg):
        sel_rules=True
    elif('--max-v1s=' in arg):
        arg=arg[idx_eq:]
        spl=arg.split(',')
        max_v1s=list(map(int,spl))
    elif('--max-v3s=' in arg):
        arg=arg[idx_eq:]
        spl=arg.split(',')
        max_v3s=list(map(int,spl))
    elif('--n-ground=' in arg):
        arg=arg[idx_eq:]
        n_ground=int(arg)
    elif('--bf-tr=' in arg):
        arg=arg[idx_eq:]
        bf_thr=float(arg)
    elif('--sort' in arg):
        sort=True
    else:
        print('Unknown option {0}'.format(arg))
        exit(1)
if(nq==0):
    print('NQ needs to be set')
    exit(2)
if(nq<0):
    print('NQ cannot be smaller than 0')
    exit(3)

if(bf_thr>0 or sort):
    freqfile=open(file='F.INP',mode='r')
    wg=Read_FINP_Freqs('F.INP')


if(all(map(Eq0,max_v1s))):
    mc_v1=0
else:
    mc_v1=len(max_v1s)

if(all(map(Eq0,max_v3s))):
    mc_v3=0
else:
    mc_v3=len(max_v3s)

modes=range(1,nq+1)
modes_c_v1=[[[0]]]
modes_c_v3=[]
for i in range(1,mc_v1+1):
    arr=list(combinations(modes,i))
    arr=list(map(list,arr))
    modes_c_v1.append(arr)

for i in range(1,mc_v3+1):
    arr=list(combinations(modes,i))
    arr=list(map(list,arr))
    modes_c_v3.append(arr)
    # modes_c_v3.append(list(combinations(modes,i)))

filetr=open(file='FILE.TR',mode='w')

i_gr=0
for i in range(mc_v1+1):
    for j in range(len(modes_c_v1[i])):
        v1_pos=modes_c_v1[i][j]
        v1_class=0
        if(v1_pos[0]!=0):v1_class=len(v1_pos)
        v1_list=[[0]]
        if(v1_class!=0):v1_list=MakePerms(range(1,max_v1s[v1_class-1]+1),v1_class)
        for k in range(len(v1_list)):
            v1=list(v1_list[k])
            i_gr+=1
            if(i_gr>n_ground):
                filetr.close()
                print("FILE.TR written")
                exit(0)
            for l in range(mc_v3):
                for m in range(len(modes_c_v3[l])):
                    v3_pos=modes_c_v3[l][m]
                    v3_class=len(v3_pos)
                    v3_list=MakePerms(range(1,max_v3s[v3_class-1]+1),v3_class)
                    for n in range(len(v3_list)):
                        v3=list(v3_list[n])
                        if(v3==v1 and v3_pos==v1_pos):continue
                        if(sel_rules):
                            isallowed=IsTrAllowed(v1,v1_pos,v1_class,v3,v3_pos,v3_class)
                            if(not isallowed):continue
                        strr=Tr2Str(v1,v1_pos,v1_class,v3,v3_pos,v3_class)
                        filetr.write(strr+"\n")
    
filetr.close()
print("FILE.TR written")
exit(0)
