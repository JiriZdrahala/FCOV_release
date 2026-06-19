program getfreqs
   use iso_fortran_env
   implicit none
   integer i,nq,n3,nat
   character(80) finp,nq_ch,line
   logical fex
   double precision,allocatable :: S(:,:),w(:),r(:)
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc>1)then
      call GET_COMMAND_ARGUMENT(1,finp)
   else
      finp='F.INP'
   end if
   inquire(file=finp,exist=fex)
   if(.not. fex)then
      write(output_unit,*)'F.INP not found'
      call exit(1)
   end if
   open(99,file=finp)
   read(99,'(A80)')line
   close(99)
   read(line,*)nq,n3,nat
   
   allocate(S(n3,n3),w(NQ),r(N3)
   call readsi(n3,S,w,nq,line,finp,r,.true.)
   call printFreqs(w,nq)
   
   deallocate(S,w,r)
   
   
   contains
   
   SUBROUTINE printFreqs(w,nq)
      integer nq,J
      double precision w
      
      WRITE(output_unit,4000)(w(NQ-J+1),J=1,NQ)
4000  FORMAT(6F11.6)
   end SUBROUTINE printFreqs
   
!     Ripped from dusch
      SUBROUTINE readsi(N3,S,E,NQ,fn,r,ldz)
!     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,iz,i1,i2
      real*8 S(N3,N3),E(*),CM,r(*)
      logical ldz
      character*(*) fn
      CM=219474.630d0
!     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn)
      read(4,*)NQ,nat,nat
      do 1 i=1,NAT
1     read(4,*)r(3*(i-1)+1),(r(3*(i-1)+ix),ix=1,3)
      read(4,*)
      DO 2 I=1,NAT
      DO 2 J=1,NQ
2     read(4,*)i1,i2,(s(3*(i-1)+ix,NQ-J+1),ix=1,3)
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
       if(iz.gt.0)write(6,*)iz,' zero modes: deleted'
      endif
      
      write(6,*)NQ,' vibrational modes considered'
      
      ! DO 3 I=1,NQ
! 3     E(I)=E(I)/CM
      
      RETURN
      end
   
end program getfreqs