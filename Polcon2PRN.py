#!/usr/bin/env python3

import sys
import numpy as np
import math

def FindIntervalEnd(x:list[float],x_end:float,start:int=0)->int:
    idx=-1
    for i in range(start,len(x),1):
        if(x[i]>x_end):
            idx=i
            return idx

def Conv(x:float,x0:float,y0:float,fwhm:float,prof:str)->float:
    if(prof=='G'):
        return Gauss(x,x0,y0,fwhm)
    elif(prof=='L'):
        return Lorentz(x,x0,y0,fwhm)
    else:
        print('Unknow profile function: '+prof)
        exit(3)

#see Polavarapu2017
def Gauss(x,x0,I,fwhm):
    sk=fwhm/2.3548
    yk0=I/(math.sqrt(2.0)*sk*math.sqrt(math.pi))
    return yk0*math.exp(-0.5*((x-x0)/sk)**2)

def Lorentz(x,x0,I,fwhm):
    dk=fwhm/2.0
    yk0=I/(dk*math.pi)
    return yk0*((dk**2)/((x-x0)**2+dk**2))

filename=sys.argv[1]
file=open(filename,'r')

argc=len(sys.argv)-1
step=2
if(argc>=2):step=float(sys.argv[2])

profile=''
if(argc>=3):profile=sys.argv[3].upper()

fwhm=10
if(argc>=4):fwhm=float(sys.argv[4])

line=file.readline()
line=file.readline()
while(line):
    if('Init.state' in line):
        linesplit_info=line.split()
        
        #class 0 to class 1 transition only
        v1=' '.join(linesplit_info[3:5])
        v1=v1.split('^')
        
        v3=' '.join(linesplit_info[5:7])
        v3=v3.split('^')
        
        line=file.readline()
        x=[]
        y=[]
        while(line != "\n"):
            lineSplit=line.split()
            lenn=len(lineSplit)
            wjn=float(lineSplit[lenn-2])
            wjn=1e7/wjn
            x.append(wjn)
            inten=float(lineSplit[0])
            y.append(inten)       
            line=file.readline()
        
        x_order=[i[0] for i in sorted(enumerate(x),key=lambda y:y[1])]
        x=[x[i] for i in x_order]
        y=[y[i] for i in x_order]
        
        i=1
        # while(i<len(x)):
            # if(x[i]-x[i-1]<x_tol):
                # x_el1=x.pop(i-1)
                # x_el2=x.pop(i-1)
                # avg=(x_el1+x_el2)/2
                # inten1=y.pop(i-1)
                # inten2=y.pop(i-1)
                # intenNew=inten1+inten2
                
                # x.insert(i-1,avg)
                # y.insert(i-1,intenNew)
            # i+=1
            
        
        x_new=[]
        y_new=[]
        n=len(x)
        if(n>0):
            if(step>0):
                xmin=x[0]
                xmax=x[n-1]
                n_steps=round((xmax-xmin)/step)+1
                dx=(xmax-xmin)/n_steps
                idx1=0
                for i in range(n_steps):
                    cur_x=i*dx+xmin
                    actual_x=(i+0.5)*dx+xmin
                    idx2=FindIntervalEnd(x,cur_x,idx1)
                    x_new.append(actual_x)
                    y_new.append(sum(y[idx1:idx2]))
                    idx1=idx2
            else:
                x_new=x
                y_new=y
            n=len(x_new)
            if(profile!=''):
                y_old=y_new
                y_new=np.zeros(n)
                for j in range(n):
                    x0=x_new[j]
                    y0=y_old[j]
                    #x_new.append(x0)
                    for i in range(n):
                        y_new[i]+=Conv(x_new[i],x0,y0,fwhm,profile)
                        
        filename_out='POLCON_{0}.PRN'.format(v3[0])
        fileout=open(filename_out,'w')
        for i in range(len(x_new)):
            fileout.write('{0:>15.8f} {1:>11.4e}\n'.format(x_new[i],y_new[i]))
        fileout.close()
        print('WRITTEN: '+filename_out)
        line=file.readline()
    if(line==""):break
    

file.close()
