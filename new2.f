      PROGRAM INPTRY
      IMPLICIT REAL*4(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      PARAMETER (NTA=10000,MENDELEV=89)
      DIMENSION X(NTA),IX(NTA,NTA),Q(NTA),ITY(NTA),r(NTA,3)
      CHARACTER*4 CODE,st
      CHARACTER*120 s120
      CHARACTER*1 OK
      dimension amas(MENDELEV)
      data amas/1.008,4.003,
     2  6.941, 9.012,   10.810,12.011,14.007,15.999,18.998,20.179,
     3 22.990,24.305,   26.981,28.086,30.974,32.060,35.453,39.948,
     4 39.098,40.080,44.956,47.900,50.941,51.996,54.938,55.847,
     4 58.933,58.700,63.546,65.380,
     4                  69.720,72.590,74.922,78.960,79.904,83.800,
     5 85.468,87.620,88.906,91.220,92.906,95.940,98.906,101.070,
     5 102.906,106.400,107.868,112.410,
     5                 114.82,118.69,121.75,127.600,126.905,131.300,
     6 132.905,137.330,138.906,
     6                 140.120,140.908,144.240,145.000,150.400,
     6 151.960,157.250,158.925,162.500,164.930,167.260,168.934,
     6 173.040,174.970,
     6 178.490,180.948,183.850,186.207,190.200,192.220,195.090,
     6 196.967,207.590,204.370,207.200,208.980,210.000,210.001,
     6 222.02,
     7 223.000,226.025,227.028/     
      logical*4 lauto,lautod
C     
      inquire(file='AUTO',exist=lauto)
      if(.not.lauto)then
       WRITE(*,*)
       WRITE(*,*)' This program makes input for FTRY and NEW4'
       WRITE(*,*)
       WRITE(*,*)'Input:'
       WRITE(*,*)'      FILE.X geometry'
       WRITE(*,*)'        AUTO indicates automatic run (optional)'
       WRITE(*,*)'   SUBST.LST isotopic substitutions'
       WRITE(*,*)'             (with AUTO, optional)'
       WRITE(*,*)
       WRITE(*,*)'Output:'
       WRITE(*,*)'    FTRY.INP INPUT FOR FTRY'
       WRITE(*,*)
      endif
C     ******************************INPUT:****************************
      nlinepar=iargc()
      if(lauto)then
       CODE='AUTO'
      else
       CODE='INTY'
      endif
      IPB=0
      IPF=1
      IPG=0
      IPL=1
      IPP=0
      IPA=1
      IPS=0
      OPEN(7,FILE='FTRY.INP')
      WRITE(7,1111)CODE,IPB,IPF,IPG,IPL,IPP,IPA,IPS
      if(.not.lauto)then
       if(nlinepar.gt.0)then
        call getarg(1,st)
        read(st,*)N
       else
        WRITE(*,*)'  NUMBER OF THE SCALE FACTORS :'
        READ(*,*)N
       endif
      else
       N=0
      endif
      WRITE(7,333)N
333   FORMAT(20I4)
      DO 1 I=1,N
      WRITE(*,222)I
222   FORMAT(I4,'. factor :')
      READ(*,*)X(I)
      WRITE(*,223)
223   FORMAT(' To how many atoms ?')
      READ(*,*)NM
      WRITE(*,224)
224   FORMAT(' List atoms:')
      READ(*,*)(IX(I,J),J=1,NM)
      WRITE(7,4444)X(I),NM,(IX(I,J),J=1,NM)
4444  FORMAT(F9.6,I3,100I3)
1     CONTINUE
1111  FORMAT(A4,7I4)
444   FORMAT(6F12.6)
      I=0

c
      OPEN(4,FILE='FILE.X',STATUS='OLD')
      READ(4,*)
      READ(4,*)NAT
      if(NAT.gt.NTA)then
       write(6,*)'too many atoms'
       close(4)
       stop
      endif
      im=0
      DO 7 I=1,NAT
      if(I.eq.1)then
       READ(4,120)s120
120    format(a120)
       ip=0
       do 1112 j=1,len(s120)
       if(s120(j:j).eq.'.')ip=ip+1
       if(ip.eq.3)then
        do 11113 k=j,len(s120)-1
11113   if(s120(k:k).eq.' '.and.s120(k+1:k+1).ne.' ')im=im+1
        goto 1114
       endif
1112   continue
1114   backspace 4
      endif
      if(im.ge.7)then
       READ(4,*)IZ,(r(i,j),j=1,3),id,id,id,id,id,id,id,Q(I)
      else
       READ(4,*)IZ,(r(i,j),j=1,3)
       Q(I)=0.0d0
      endif
      ITY(I)=IZ
      X(I)=AMAS(IZ)
      IF(IZ.EQ.6)X(I)=12.000000
      IF(IZ.EQ.2)X(I)=15.994915
      IF(IZ.EQ.7)X(I)=14.003074
      IF(IZ.EQ.1)X(I)=1.007825
      IF(IZ.EQ.16)X(I)=31.972000
      IF(IZ.EQ.15)X(I)=30.973700
      IF(IZ.EQ.9)X(I)=18.9984
7     IF(IZ.EQ.17)X(I)=34.968852
      CLOSE(4)
c
2     if(.not.lauto)WRITE(*,*)NAT,' ATOMS FOUND IN FILE.X'
      if(.not.lauto)then
       if(nlinepar.gt.1)then
        call getarg(1,st)
        read(st,*)ISUB
       else
        WRITE(*,*)'HOW MANY ISOTOPIC SUBSTITUTION (neg. for auto-D)?'
        READ(*,*)ISUB
       endif
       if(ISUB.gt.0)then
        DO 5 II=1,ISUB
        WRITE(*,*)' H     1.007825  D    2.014000  T     3.016050'
        WRITE(*,*)' C12  12.000000  C13 13.003355  N14  14.003074'
        WRITE(*,*)' N15  15.000108  O16 15.994915  O17  16.999131'
        WRITE(*,*)' O18  17.999160  F   18.998403  CL35 34.968852'
        WRITE(*,*)' CL37 36.965903  P   30.973762  S    31.972070'
        WRITE(*,*)' GIVE NUMBER OF THE ATOM AND THE NEW MASS'
5       READ(*,*)I,X(I)
       endif
       DO 6 I=1,NAT,3
       IF(NAT-I.GE.2)WRITE(*,6676)I,X(I),I+1,X(I+1),I+2,X(I+2)
       IF(NAT-I.EQ.1)WRITE(*,6676)I,X(I),I+1,X(I+1)
6      IF(NAT-I.EQ.0)WRITE(*,6676)I,X(I)
6676   FORMAT(3(I8,F12.6))
       if(nlinepar.gt.2)then
        call getarg(1,OK)
       else
        WRITE(*,*)'Is it correct (Y/N) ?'
        READ(*,'(A)')OK
       endif
       IF ((OK.EQ.'n').OR.(OK.EQ.'N')) GOTO 2
      else
       inquire(file='SUBST.LST',exist=lautod)
       if(lautod)then
        open(38,file='SUBST.LST')
        read(38,*)ISUB
        if(ISUB.gt.0)then
         DO 51 II=1,ISUB
51       READ(38,*)I,X(I)
        endif
        close(38)
       else
        ISUB=0
       endif
      endif

      if(ISUB.lt.0)then
       ND=0
c
c      Substitute all acidic Hs by Ds
       do 8 i=1,NAT
       a=r(i,1)
       y=r(i,2)
       z=r(i,3)
       if(ITY(i).eq.1)then
        do 3 j=1,NAT
        it=ITY(j)
        dd=sqrt((a-r(j,1))**2+(y-r(j,2))**2+(z-r(j,3))**2)
        if(it.eq.7.or.it.eq.8.or.it.eq.9.or.it.eq.16.or.it.eq.17.or.
     1  it.eq.35.or.it.eq.53)then
         if(dd.lt.1.23)then
          ND=ND+1
          X(i)=2.014000
         endif
        endif
3       continue
       endif
8      continue
       if(.not.lauto)write(6,7009)ND
7009   format(i5,' acidic hydrogens have been substituted by deuteria')
      endif

      WRITE(7,444)(Q(I),I=1,NAT)
      I1=1
      WRITE(7,333)I1
      WRITE(7,444)(X(I),I=1,NAT)
      WRITE(7,*)3*NAT
      A0=0.0
      WRITE(7,444)(A0,I=1,3*NAT)
      CLOSE(7)
      if(.not.lauto)WRITE(*,*)' PROGRAM TERMINATED'
      END

