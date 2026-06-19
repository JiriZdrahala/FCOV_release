program swph
   use iso_fortran_env
   implicit none
   
   
   
   contains
   
      SUBROUTINE SMATOUT(S,AS,NAT3,nq,JOB,ili,nlinepar,IRANGE,JRANGE)
!     dX  =  S dQ   S(nat3xnint)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 JOB
      CHARACTER*6 num
      DIMENSION S(NAT3,NAT3),AS(3,NAT3)
      IRANGE=1
      JRANGE=NAT3
      
!     look at the phase and make the biggest element always positive:
      DO 101 J=IRANGE,JRANGE
      amx=S(1,J)
      DO 102 I=2,NAT3
102   if(abs(S(I,J)).gt.abs(amx))amx=S(I,J)
      if(amx.lt.0)then
       do 103 I=1,NAT3
103    S(I,J)=-S(I,J)
      endif
101   continue
!
      OPEN(34,FILE='F.INP')
      WRITE(34,10)NI,NAT3,NAT3/3
10    FORMAT(3I7)
      OPEN(35,FILE='FILE.X',STATUS='OLD')
      READ(35,*)
      READ(35,*)NAT
      DO 3 I=1,NAT
      READ(35,*)IAT,X,Y,Z
3     WRITE(34,11)IAT,X,Y,Z
11    FORMAT(I7,3F12.6)
      CLOSE(35)
      WRITE(34,14)
14    FORMAT(' Atom Mode    X-disp.    Y-disp.      Z-disp.')
      DO 1 I=1,NAT3/3
      DO 1 J=IRANGE,JRANGE
      if(dabs(AS(1,J)).gt.0.1d0.or.ize.eq.0)then
       IU=3*(I-1)
       K=NAT3-J+1
       WRITE(34,15)I,K,S(IU+1,J),S(IU+2,J),S(IU+3,J)
15     FORMAT(2I7,3F11.6)
      endif
1     continue
      WRITE(34,11)NI
      II=0
      do 4 I=IRANGE,JRANGE
      if(dabs(AS(1,I)).gt.0.1d0.or.ize.eq.0)then
       II=II+1
       WRITE(34,13)AS(1,I)
13     format(F11.3,$)
       if(mod(II,6).eq.0)write(34,*)
      endif
4     continue
      write(34,*)
      CLOSE(34)
      RETURN
      END SUBROUTINE SMATOUT
   
   SUBROUTINE readsi(N3,S,E,NQ,fn,r,ldz,iz)
!     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,iz,fuck1,fuck2
      real*8 CM
      real*8,allocatable :: S(:,:), S_help(:,:),r(:)
      real*8,allocatable :: E(:),E_help(:)
      logical ldz
      character*(*) fn
      CM=219474.630d0
!     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn)
      read(4,*)NQ,nat,nat
      N3=3*nat
      allocate(r(3*nat),S(3*nat,NQ),E(NQ))
      do 1 i=1,NAT
1     read(4,*)r(3*(i-1)+1),(r(3*(i-1)+ix),ix=1,3)
      read(4,*)
      DO 2 I=1,NAT
      DO 2 J=1,NQ
2     read(4,*)fuck1,fuck2,(s(3*(i-1)+ix,NQ-J+1),ix=1,3)
      read(4,*)
      READ(4,4000)(E(NQ-J+1),J=1,NQ)
4000  FORMAT(6F11.6)
      close(4)
!
      write(6,*)NQ,' modes found'
!     delete zero modes if exist:
      if(ldz)then
       iz=0
66     do 6 i=1,NQ
       if(dabs(E(i)).lt.0.1d0)then
        do 7 j=i,NQ-1
        E(j)=E(j+1)
        do 7 ix=1,N3
        s(ix,j)=s(ix,j+1)
7       continue
        iz=iz+1
        NQ=NQ-1
        goto 66
       endif
6      continue
       if(iz.gt.0)then
         write(6,*)iz,' zero modes: deleted'
         allocate(S_help(N3,NQ),E_help(NQ))
         S_help=S(:,1:NQ)
         E_help=E(1:NQ)
         deallocate(S,E)
         S=S_help
         E=E_help
         deallocate(S_help,E_help)
       endif 
      endif
      
      
      write(6,*)NQ,' vibrational modes considered'
      
      DO 3 I=1,NQ
3     E(I)=E(I)/CM
      
      RETURN
      end SUBROUTINE readsi

end program swph