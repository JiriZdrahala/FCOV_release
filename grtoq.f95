program grtoq
   use iso_fortran_env
   implicit none
   integer n3,nq,iz
   double precision,allocatable :: w(:),S(:,:),r(:)
   double precision,allocatable :: gr_car(:),gr_q(:)

   call readsi(N3,S,w,NQ,'F.INP',r,.true.,iz)
   allocate(gr_car(n3))
   call readgrad(n3/3,gr_car,'FILE.GR')
   allocate(gr_q(nq))
   call grad2Q(N3/3,nq,gr_car,gr_q,S)
   call svgrq(gr_q,nq,'FILE.Q.GR')
   write(output_unit,*)'FILE.Q.GR written'
   
   deallocate(gr_car,gr_q,w,s,r)
contains

   subroutine svgrq(g,nq,fille)
      implicit none
      character(*) fille
      real*8 g(*)
      integer*4 nq,i,ix
      open(70,file=fille)
      do 1 i=1,nq
1     write(70,700)g(i)
700   format(f15.9)
      close(70)
      write(6,*)'Gradient written to '//fille
      return
   end
   
   subroutine grad2Q(nat,NQ,gr,gr_q,L)
      implicit none
      integer*4 nat,i,ix,NQ
      real*8 gr(3*nat),gr_q(NQ),L(3*nat,NQ)!,hold(3*nat)
      
      L=L*0.0234280d0
      
      gr_q=matmul(transpose(L),gr)
   end subroutine grad2Q

   subroutine readgrad(nat,gr,filee)
      implicit none
      integer*4 nat,i,ix
      real*8 gr(3*nat)
      character(*) filee
      
      open(77,file=filee)
      do i = 1,nat
         ix=(i-1)*3+1
         read(77,*)gr(ix),gr(ix+1),gr(ix+2)
      end do
      close(77)
   end subroutine readgrad

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
      end

end program grtoq