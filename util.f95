#define TR(arg) trim(adjustl(arg))
module util
   use iso_fortran_env
   use constants
   use strings
   implicit none
   
   type Polar
      double complex :: ap(3,3),G(3,3),Gc(3,3),A(3,3,3),Ac(3,3,3)
      contains
      procedure, pass(pol) :: Polar_assign_d
      generic, public :: assignment(=) => Polar_assign_d
   end type Polar
   
   
   interface operator(+)
      module procedure :: Polars_add
   end interface
   interface operator(*)
      module procedure :: polars_sc_mult,Polars_sc_mult_left
   end interface
   
   type TD_coeff
      integer :: gr=-1,ex=-1
      character(3) :: sym_gr='   ',sym_ex='   '
      double precision :: c=HUGE(1d0)
   end type TD_coeff
   
   type TD_coeff_arr
      type(TD_coeff),allocatable :: arr(:)
   end type TD_coeff_arr
   
   type transition
      type(Polar) polarr
      double precision w1,w3
      integer :: v1_n=HUGE(1),v3_n=HUGE(1)
      integer,allocatable :: v1(:),v3(:)
      integer(int16),allocatable :: v1_pos(:),v3_pos(:)
   end type transition
   
   
   contains
   
   function KD(a,b)result(res)
      integer a,b,res
      
      if(a==b)then
         res=1
      else
         res=0
      end if
   end function KD

   
   function LC(a,b,c)result(res)
      integer a,b,c,res
      
      if(a*b*c/=6)then
         res=0
         return
      end if
      
      if(a==1)then
         if(b==2)then
            res=1
         else
            res=-1
         end if
      else if(a==2)then
         if(b==3)then
            res=1
         else
            res=-1
         end if
      else if(a==3)then
         if(b==1)then
            res=1
         else
            res=-1
         end if
      end if
   end function LC
   
   function VelToLen_vec_cart(v,e_eg)result(res)
      double precision v(3),e_eg,res(3)
      res=-v/e_eg
   end function VelToLen_vec_cart
   
   function VelToLen_vec_cartQ(v,e_eg)result(res)
      double precision v(3,3),e_eg,res(3,3)
      res=-v/e_eg
   end function VelToLen_vec_cartQ
   
   function VelToLen_dervec_cart(v,dv,e_eg,grad_e,grad_g,n)result(res)
      integer n,i
      double precision dv(3,n),v(3),e_eg,grad_e(n),grad_g(n),res(3,n)
      
      do i = 1,n
         res(:,i)=v/e_eg**2 * (grad_e(i)-grad_g(i)) - dv(:,i)/e_eg
      end do
   end function VelToLen_dervec_cart

   function VelToLen_dervec_cartQ(v,dv,e_eg,grad_e,grad_g,n)result(res)
      integer n,i
      double precision dv(3,3,n),v(3,3),e_eg,grad_e(n),grad_g(n),res(3,3,n)
      
      do i = 1,n
         res(:,:,i)=v/e_eg**2 * (grad_e(i)-grad_g(i)) - dv(:,:,i)/e_eg
      end do
   end function VelToLen_dervec_cartQ
   
   function VelToLen_der2(v,dv,dv2,e_eg,grad_ex,grad_gr,ff_ex,ff_gr,n)result(res)
      integer n,k,l
      double precision res(3,n,n),v(3),dv(3,n),dv2(3,n,n),e_eg
      double precision grad_ex(n),grad_gr(n),ff_ex(n,n),ff_gr(n,n)
      double precision buf1(3),buf2(3),buf3(3),buf4(3),buf5(3)
      
      res=0d0
      do k = 1,n
         do l = 1,n
            buf1=-(grad_ex(l)-grad_gr(l))/e_eg**2*dv(:,k)
            buf2=1d0/e_eg*dv2(:,l,k)
            buf3=-(ff_ex(l,k)-ff_gr(l,k))/e_eg**2*v
            buf4=2*(grad_ex(k)-grad_gr(k))*(grad_ex(l)-grad_gr(l))/e_eg**3*v
            buf5=(grad_ex(k)-grad_gr(k))/e_eg**2*dv(:,l)
            res(:,l,k)=-(buf1+buf2+buf3+buf4+buf5)
         end do
      end do
   end function VelToLen_der2
   
   function VelToLen_der2q(v,dv,dv2,e_eg,grad_ex,grad_gr,ff_ex,ff_gr,n)result(res)
      integer n,k,l
      double precision res(3,3,n,n),v(3,3),dv(3,3,n),dv2(3,3,n,n),e_eg
      double precision grad_ex(n),grad_gr(n),ff_ex(n,n),ff_gr(n,n)
      double precision buf1(3,3),buf2(3,3),buf3(3,3),buf4(3,3),buf5(3,3)
      
      res=0d0
      do k = 1,n
         do l = 1,n
            buf1=-(grad_ex(l)-grad_gr(l))/e_eg**2*dv(:,:,k)
            buf2=1d0/e_eg*dv2(:,:,l,k)
            buf3=-(ff_ex(l,k)-ff_gr(l,k))/e_eg**2*v
            buf4=2*(grad_ex(k)-grad_gr(k))*(grad_ex(l)-grad_gr(l))/e_eg**3*v
            buf5=(grad_ex(k)-grad_gr(k))/e_eg**2*dv(:,:,l)
            res(:,:,l,k)=-(buf1+buf2+buf3+buf4+buf5)
         end do
      end do
   end function VelToLen_der2q
   
   
   function RealPart_pol(pol)result(res)
      type(Polar) pol,res
      
      res=0d0
      res%ap=realpart(pol%ap)
      res%G=realpart(pol%G)
      res%Gc=realpart(pol%Gc)
      res%A=realpart(pol%A)
      res%Ac=realpart(pol%Ac)
   end function RealPart_pol
   
   function ImagPart_pol(pol)result(res)
      type(Polar) pol,res
      
      res=0d0
      res%ap=imagpart(pol%ap)
      res%G=imagpart(pol%G)
      res%Gc=imagpart(pol%Gc)
      res%A=imagpart(pol%A)
      res%Ac=imagpart(pol%Ac)
   end function ImagPart_pol
   
   function TD_coeffs_dotproduct(tdc1,tdc2)result(res)
      type(TD_coeff) tdc1(:),tdc2(:)
      type(TD_coeff) tdc_cur1,tdc_cur2
      integer n1,n2,i,j
      double precision res
      
      n1=size(tdc1,dim=1)
      n2=size(tdc2,dim=1)
      res=0d0
      
      outer: do i=1,n1
         tdc_cur1=tdc1(i)
         do j=1,n2
            tdc_cur2=tdc2(j)
            if(tdc_cur1%gr==tdc_cur2%gr .and. tdc_cur1%ex==tdc_cur2%ex)then
               if(tdc_cur1%sym_gr==tdc_cur2%sym_gr .and. tdc_cur1%sym_ex==tdc_cur2%sym_ex)then
                  res=res+tdc_cur1%c*tdc_cur2%c
                  cycle outer
               end if
            end if
         end do 
      end do outer
   end function TD_coeffs_dotproduct
   
   subroutine StateSwitchCheck(tdcs0,tdcs,n,order)
      integer n,idx,bufI
      integer,allocatable :: order(:)
      type(TD_coeff_arr) tdcs0(n),tdcs(n),tdcs_buf
      double precision,allocatable :: S(:)
      
      integer i,j
      
      allocate(S(n)) !overlap matrix
      
      order=seq(n,.false.)
      do i = 1,n
         do j = 1,n
            S(j)=abs(TD_coeffs_dotproduct(tdcs0(i)%arr,tdcs(j)%arr))
         end do
         idx=maxloc(S,dim=1)
         if(idx/=i)then
            bufI=order(i)
            order(i)=order(idx)
            order(idx)=bufI
            
            tdcs_buf=tdcs(i)
            tdcs(i)=tdcs(idx)
            tdcs(idx)=tdcs_buf
         end if
      end do
      
      
      deallocate(S,tdcs_buf%arr)
   end subroutine StateSwitchCheck
   
   
   
   function seq(n,zero)result(res)
      integer n,res(n),i,start,endd
      logical zero
      
      start=1
      endd=n
      if(zero)then
         start=0
         endd=n-1
      end if
      
      do i = start,endd
         res(i)=i
      end do
      
   end function seq
   
   function Polars_add(pol1,pol2)result(polRes)
      type(Polar),intent(in) :: pol1,pol2
      type(Polar) :: polRes
      integer :: i,j,k,l
      
      do j = 1,3
         do i = 1,3
            polRes%Ap(i,j)=pol1%Ap(i,j)+pol2%Ap(i,j)
            polRes%G(i,j)=pol1%G(i,j)+pol2%G(i,j)
            polRes%Gc(i,j)=pol1%Gc(i,j)+pol2%Gc(i,j)
            do k = 1,3
               polRes%A(i,j,k)=pol1%A(i,j,k)+pol2%A(i,j,k)
               polRes%Ac(i,j,k)=pol1%Ac(i,j,k)+pol2%Ac(i,j,k)
            end do
         end do
      end do
   end function Polars_add
   
   Subroutine Polar_assign_c(pol,vall)
      type(Polar), Intent(Out) :: pol
      double complex, Intent(In) :: vall
      
      pol%ap=vall
      pol%G=vall
      pol%Gc=vall
      pol%A=vall
      pol%Ac=vall
   End Subroutine Polar_assign_c
   
   function Polars_sc_mult(pol1,fac)result(pol)
      type(Polar),intent(in) :: pol1
      type(Polar) pol
      integer a,b,c
      double precision,intent(in) :: fac
      
      do a = 1,3
         do b = 1,3
            pol%ap(a,b)=pol1%ap(a,b)*fac
            pol%G(a,b)=pol1%G(a,b)*fac
            pol%Gc(a,b)=pol1%Gc(a,b)*fac
            do c = 1,3
               pol%A(a,b,c)=pol1%A(a,b,c)*fac
               pol%Ac(a,b,c)=pol1%Ac(a,b,c)*fac
            end do
         end do
      end do
   end function Polars_sc_mult
   
   function Polars_sc_mult_left(fac,pol1)result(pol)
      type(Polar),intent(in) :: pol1
      type(Polar) pol
      integer a,b,c
      double precision,intent(in):: fac
      
      do a = 1,3
         do b = 1,3
            pol%ap(a,b)=pol1%ap(a,b)*fac
            pol%G(a,b)=pol1%G(a,b)*fac
            pol%Gc(a,b)=pol1%Gc(a,b)*fac
            do c = 1,3
               pol%A(a,b,c)=pol1%A(a,b,c)*fac
               pol%Ac(a,b,c)=pol1%Ac(a,b,c)*fac
            end do
         end do
      end do
   end function Polars_sc_mult_left
   
   
   
   !Why is it recursive you might ask??? GCC OpenMP was complaining. I checked and did not identify the recursion, not even in the call stack
   pure recursive Subroutine Polar_assign_d(pol,vall)
      integer a,b,c
      class(Polar),intent(inout) :: pol
      double precision, Intent(In) :: vall
      
      do a = 1,3
         do b = 1,3
            pol%ap(a,b)=vall
            pol%G(a,b)=vall
            pol%Gc(a,b)=vall
            do c = 1,3
               pol%A(a,b,c)=vall
               pol%Ac(a,b,c)=vall
            end do
         end do
      end do
   End Subroutine Polar_assign_d 
   
   Subroutine Polar_assign_i(pol,vall)
      Implicit None
      type(Polar), Intent(Out) :: pol
      integer, Intent(In) :: vall
      
      pol%ap=vall
      pol%G=vall
      pol%Gc=vall
      pol%A=vall
      pol%Ac=vall
   End Subroutine Polar_assign_i 
   
   subroutine readDusch(NQ,N3,A,B,C,D,E,J,K,SG,SE,wg,we,filee)
      character(80) c80
      character(*) filee
      integer bufi,i,col,NQ,N3
      double precision,allocatable :: A(:,:),B(:),C(:,:),D(:),E(:,:),J(:,:),K(:),SG(:,:),SE(:,:),wg(:),we(:)
      integer,parameter :: unitt=77
      
      open(unitt,file=filee)
      read(unitt,*)NQ
      ! if(N/=NQ)THEN
         ! write(output_unit,*)'Number of modes in DUSCH.OUT not equal to 3*nat-6'
         ! call exit(2)
      ! end if
      allocate(A(NQ,NQ),B(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ),J(NQ,NQ),K(NQ),SG(N3,NQ),SE(N3,NQ))
      allocate(we(nQ))
      if(.not.allocated(wg))allocate(wg(nq))
      !allocate(A_small(N,N),B_small(N),C_small(N,N),D_small(N),E_small(N,N))
      
      ! A_small=.false.
      ! B_small=.false.
      ! C_small=.false.
      ! D_small=.false.
      ! E_small=.false.
      
      call Forward(unitt,3)
      call readMatrix(unitt,J,NQ)
      
      call Forward(unitt,3)
      call readVector(unitt,K,NQ)
      
      call Forward(unitt,3)
      call readMatrix(unitt,A,NQ)
      
      call Forward(unitt,3)
      call readVector(unitt,B,NQ)
      
      call Forward(unitt,3)
      call readMatrix(unitt,C,NQ)
      
      call Forward(unitt,3)
      call readVector(unitt,D,NQ)
      
      call Forward(unitt,3)
      call readMatrix(unitt,E,NQ)
      
      call Forward(unitt,3)
      call readFreqs(unitt,wg,nQ)
      
      call Forward(unitt,3)
      call readFreqs(unitt,we,nQ)
      
      call Forward(unitt,3)
      call readNonSquareMatrix(unitt,SG,N3,NQ)
      
      call Forward(unitt,3)
      call readNonSquareMatrix(unitt,SE,N3,NQ)
      
      !I swear to god, mankind should start jailing programmers who dont document these conversions
      !I dont think I would say this, but thank Christ Data management plans are starting to become a thing
      !If scientists didn't do sloppy documentation, you would not have to deal with DMPs
      ! SG=SG*0.0234280d0 !this is conversion from amu^-1/2 to au^-1/2
      ! SE=SE*0.0234280d0
      SG=SG*1d0/sqrt(amu_2_au)
      SE=SE*1d0/sqrt(amu_2_au)
      close(unitt)
   end subroutine readDusch
   
   subroutine readFreqs(unitt,w,n)
      integer unitt,n,row,i,idx
      double precision w(n)
      character(80) fn
      
      do row = 1,ceiling(N/6.0),1
         idx=6*(row-1)+1
         read(unitt,'(A80)')fn
         read(fn,*)w(idx:min(idx+5,N))
      end do
   end subroutine readFreqs
   
   subroutine readNonSquareMatrix(unitt,m,n1,n2)
      integer unitt,n1,n2,box,i,idx,bufi
      double precision m(n1,n2)
      character(80) fn
      
      do box = 1,ceiling(n2/5.0),1
         read(unitt,*)
         idx=(box-1)*5+1
         do i = 1,n1
            read(unitt,'(A80)')fn
            !untransposed matrices
            read(fn,*)bufi,m(i,idx:min(idx+4,n2))
         end do
      end do
   end subroutine readNonSquareMatrix
   
   subroutine readMatrix(unitt,m,n)
      integer unitt,n,bufi,i,idx,box
      double precision m(n,n)
      character(80) fn
      
      !the matrices are transposed because Fortran stores 2d-arrays in a column-major order
      !so I am just following the principle of caching
      do box = 1,ceiling(n/5.0),1
         idx=5*(box-1)+1
         read(unitt,*)
         do i = 1,N
            read(unitt,'(A80)')fn
            read(fn,*)bufi,m(idx:min(idx+4,N),i)
         end do
      end do
   end subroutine readMatrix
   
   subroutine readVector(unitt,v,n)
      integer unitt,n,bufi,i
      double precision v(n)
      
      read(unitt,*)
      do i = 1,N
         read(unitt,*)bufi,v(i)
      end do
   end subroutine readVector
   
   
   SUBROUTINE readsi(N3,S,E,NQ,fn,z_at,r,ldz,iz)
!     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,iz,fuck1,fuck2
      integer,allocatable :: z_at(:)
      real*8 CM
      real*8,allocatable :: S(:,:), S_help(:,:),r(:)
      real*8,allocatable :: E(:),E_help(:)
      logical ldz
      character*(*) fn
      CM=219474.630d0
!     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn,status='old')
      read(4,*)NQ,nat,nat
      N3=3*nat
      allocate(z_at(nat),r(3*nat),S(3*nat,NQ),E(NQ))
      do i=1,NAT
      read(4,*)z_at(i),(r(3*(i-1)+ix),ix=1,3)
      end do
      read(4,*)
      DO I=1,NAT
      DO J=1,NQ
      read(4,*)fuck1,fuck2,(s(3*(i-1)+ix,NQ-J+1),ix=1,3)
      end do
      end do
      read(4,*)
      READ(4,4000)(E(NQ-J+1),J=1,NQ)
4000  FORMAT(6F11.6)
      close(4)
!
      !write(6,*)NQ,' modes found'
!     delete zero modes if exist:
      if(ldz)then
       iz=0
66     do i=1,NQ
       if(dabs(E(i)).lt.0.1d0)then
        do j=i,NQ-1
        E(j)=E(j+1)
        do ix=1,N3
        s(ix,j)=s(ix,j+1)
        enddo
        enddo
        iz=iz+1
        NQ=NQ-1
        goto 66
       endif
       end do
       if(iz.gt.0)then
         !write(6,*)iz,' zero modes: deleted'
         allocate(S_help(N3,NQ),E_help(NQ))
         S_help=S(:,1:NQ)
         E_help=E(1:NQ)
         deallocate(S,E)
         S=S_help
         E=E_help
         deallocate(S_help,E_help)
       endif 
      endif
      
      S=S*1d0/sqrt(amu_2_au)
      !write(6,*)NQ,' vibrational modes considered'
      
      ! DO 3 I=1,NQ
! 3     E(I)=E(I)/CM
      
      RETURN
   end

   SUBROUTINE writesi(N3,nat,nq,S,E,fn,z_at,r,convert)
!     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,fuck1,fuck2
      integer :: z_at(nat)
      real*8 CM
      real*8 :: S(n3,nq), r(n3)
      real*8 :: E(nq)
      logical convert
      character*(*) fn
      
      
      CM=219474.630d0
      if(convert)S=S/(1d0/sqrt(amu_2_au))
      
!     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn)
      
      write(4,10)NQ,n3,nat
10    FORMAT(3I7)
      do i=1,NAT
      write(4,11)z_at(i),(r(3*(i-1)+ix),ix=1,3)
11    FORMAT(I7,3F12.6)
      end do
      write(4,14)
14    FORMAT(' Atom Mode    X-disp.    Y-disp.      Z-disp.')
      DO I=1,NAT
      DO J=1,NQ
      write(4,15)I,NQ-J+1,(s(3*(i-1)+ix,NQ-J+1),ix=1,3)
      end do
      end do
15    FORMAT(2I7,3F11.6)
      write(4,11)NQ
      write(4,4000)(E(NQ-J+1),J=1,NQ)
4000  FORMAT(6F11.3)
      close(4)
      
      if(convert)S=S*1d0/sqrt(amu_2_au)
   end subroutine writesi


   function Equal_IArr(arr1,arr2)result(res)
      integer n,m,arr1(:),arr2(:),i
      logical res
      
      n=size(arr1,dim=1)
      m=size(arr2,dim=1)
      if(n==0 .and. m==0)then
         res=.true.
         return
      end if
      res=.false.
      if(m/=n)return
      do i = 1,n
         if(arr1(i)/=arr2(i))then
            return
         end if
      end do
      res=.true.
   end function Equal_IArr
   
   function FC2Str_new(v1,v1_pos,v1_n,gr,printZero)result(str)
      integer v1_n,v1(v1_n),i,ii,leng
      integer(int16) v1_pos(v1_n)
      logical gr,printZero
      integer,parameter :: str_len = 160
      character(str_len) :: str
      character(80) str_help
      
      if(all(v1==0))then
         if(gr)then
            write(str,'(A)')'|0 >'
         else
            write(str,'(A)')'< 0|'
         end if
         return
      end if
      
      write(str,*)' '
      if(gr)THEN
         write(str,'(A)')'|'
      else
         write(str,'(A)')'<'
      end if
      
      ii=2
      do i = 1,v1_n
         if(v1(i)==0 .and. .not. printZero)cycle
         write(str_help,'(I0,A,I0)')v1_pos(i),'^',v1(i)
         leng=Ciphers(int(v1_pos(i),4))+1+Ciphers(v1(i))+1
         str(ii:ii+leng)=TR(str_help)
         ii=ii+leng
      end do
      
      if(gr)then
         str(ii:ii)='>'
      else
         str(ii:ii)='|'
      end if
   end function FC2Str_new
   
   function ReadPolars(fn,n,wexc,e00,isNM)result(polres)
      character(*) fn
      integer n,i
      double precision wexc,e00
      logical endd,isNM
      type(polar),allocatable :: polres(:)
      type(polar) :: curPol
      
      endd=.false.
      open(77,file=fn,status='old')
      read(77,*)n
      allocate(polres(n))
      i=0
      do while(.not. endd)
         i=i+1
         call ReadTenPretty_notr(77,endd,curPol,wexc,isNM)
         if(endd)exit
         polres(i)=curPol
      end do
      close(77)
   end function ReadPolars
   
   subroutine ReadTenPretty(unitt,endd,tr,wexc)
      integer,intent(in) :: unitt
      integer i,j,k
      integer :: v1_n,v3_n
      integer,allocatable :: v1(:),v3(:)
      integer(int16),allocatable :: v3_pos(:),v1_pos(:)
      type(Polar) :: polarr
      type(Transition),intent(out) :: tr
      double precision,intent(out) :: wexc
      double precision w_fi_cm,w1,w3
      logical endd
      character(80) s80
      
      
      endd=.false.
      
      read(unitt,'(1X,F8.2)',end=2002)wexc
      read(unitt,*)s80
      read(s80,*)v1_n
      if(v1_n>0)then
         allocate(v1(v1_n),v1_pos(v1_n))
         read(unitt,*)v1,v1_pos
      else
         v1=[0]
         v1_pos=[0]
         read(unitt,*)
      end if
      
      read(unitt,*)v3_n
      if(v3_n>0)then
         allocate(v3(v3_n),v3_pos(v3_n))
         read(unitt,*)v3,v3_pos
      else
         v3=[0]
         v3_pos=[0]
         read(unitt,*)
      end if
      read(unitt,'(3(1X,F8.3))')w1,w3,w_fi_cm
      
2000  format(3(1X,E15.8,1X,E15.8,2X))
2001  format(9(1X,E15.8,1X,E15.8,2X))      
      !Alpha, electric dipole-electric dipole
      read(unitt,'(A)')
      read(unitt,2000)((polarr%Ap(i,j),j=1,3),i=1,3)
      
      !G, electric dipole-magnetic dipole
      read(unitt,'(A)')
      read(unitt,2000)((polarr%G(i,j),j=1,3),i=1,3)
      
      !G-cursive, magnetic dipole-electric dipole
      read(unitt,'(A)')
      read(unitt,2000)((polarr%Gc(i,j),j=1,3),i=1,3)
      
      !A, electric dipole-electric quadrupole
      read(unitt,'(A)')
      read(unitt,2001)(((polarr%A(i,j,k),k=1,3),j=1,3),i=1,3)
      
      !A-cursive, electric quadrupole-electric dipole
      read(unitt,'(A)')
      read(unitt,2001)(((polarr%Ac(i,j,k),k=1,3),j=1,3),i=1,3)
      read(unitt,*)
      
      if(v1_n>0)then
         tr%v1=v1
         tr%v1_pos=v1_pos
         tr%v1_n=v1_n
      else
         tr%v1=[0]
         tr%v1_pos=[0]
         tr%v1_n=0
      end if
      
      if(v3_n>0)then
         tr%v3=v3
         tr%v3_pos=v3_pos
         tr%v3_n=v3_n
      else
         tr%v3=[0]
         tr%v3_pos=[0]
         tr%v3_n=0
      end if
      tr%w1=w1
      tr%w3=w3
      tr%polarr=polarr
      return
2002  endd=.true.
      
   end subroutine ReadTenPretty
   
   subroutine ReadTenPretty_notr(unitt,endd,polarr,wexc,isNM)
      integer,intent(in) :: unitt
      integer i,j,k
      integer :: v1_n,v3_n
      integer,allocatable :: v1(:),v3(:)
      integer(int16),allocatable :: v3_pos(:),v1_pos(:)
      type(Polar) :: polarr
      !type(Transition),intent(out) :: tr
      double precision,intent(out) :: wexc
      double precision w_fi_cm,w1,w3
      logical endd,isNM
      character(80) s80
      
      
      endd=.false.
      read(unitt,'(1X,F8.2)',end=2002)wexc
      wexc=(1d7/wexc)*cm_2_au
      read(unitt,*)s80
      read(s80,*)v1_n
      if(v1_n>0)then
         allocate(v1(v1_n),v1_pos(v1_n))
         read(unitt,*)v1,v1_pos
      else
         v1=[0]
         v1_pos=[0]
         read(unitt,*)
      end if
      
      read(unitt,*)v3_n
      if(v3_n>0)then
         allocate(v3(v3_n),v3_pos(v3_n))
         read(unitt,*)v3,v3_pos
      else
         v3=[0]
         v3_pos=[0]
         read(unitt,*)
      end if
      
      read(unitt,'(3(1X,F8.3))')w1,w3,w_fi_cm
      isNM=.true.
      if(v1_n==0 .and. v3_n==0)then
         isNM=.false.
      end if 
2000  format(3(1X,E15.8,1X,E15.8,2X))
2001  format(9(1X,E15.8,1X,E15.8,2X))      
      !Alpha, electric dipole-electric dipole
      read(unitt,'(A)')
      read(unitt,2000)((polarr%Ap(i,j),j=1,3),i=1,3)
      
      !G, electric dipole-magnetic dipole
      read(unitt,'(A)')
      read(unitt,2000)((polarr%G(i,j),j=1,3),i=1,3)
      
      !G-cursive, magnetic dipole-electric dipole
      read(unitt,'(A)')
      read(unitt,2000)((polarr%Gc(i,j),j=1,3),i=1,3)
      
      !A, electric dipole-electric quadrupole
      read(unitt,'(A)')
      read(unitt,2001)(((polarr%A(i,j,k),k=1,3),j=1,3),i=1,3)
      
      !A-cursive, electric quadrupole-electric dipole
      read(unitt,'(A)')
      read(unitt,2001)(((polarr%Ac(i,j,k),k=1,3),j=1,3),i=1,3)
      read(unitt,*)
      
      ! if(v1_n>0)then
         ! tr%v1=v1
         ! tr%v1_pos=v1_pos
         ! tr%v1_n=v1_n
      ! else
         ! tr%v1=[0]
         ! tr%v1_pos=[0]
         ! tr%v1_n=0
      ! end if
      
      ! if(v3_n>0)then
         ! tr%v3=v3
         ! tr%v3_pos=v3_pos
         ! tr%v3_n=v3_n
      ! else
         ! tr%v3=[0]
         ! tr%v3_pos=[0]
         ! tr%v3_n=0
      ! end if
      ! tr%w1=w1
      ! tr%w3=w3
      ! tr%polarr=polarr
      return
2002  endd=.true.
      
   end subroutine ReadTenPretty_notr

   subroutine WritePolars(el,nq,e00,unitt,wexc,polars,wg,createFile,isNM)
      type(Polar) polars(nq)
      logical el,createFile,isNM
      integer unitt,i
      integer v1_n,v3_n,nq
      integer,allocatable :: v1(:)
      integer,allocatable :: v3(:)
      integer(int16),allocatable :: v3_pos(:),v1_pos(:)
      double precision w1,w3,wexc_new,wexc,wg(nq),e00
      character(80) filen
      character(5) w_nm_str
      
      if(wexc>1d-40 .and. createFile)then
         wexc_new=1d7/(wexc*au_2_cm)
         write(w_nm_str,'(I5)')NINT(wexc_new)
         if(el)then
            if(isNM)then
               write(filen,'(A,A,A)')'FILE.',TR(w_nm_str),'nm.POLARS.Q.EL'
            else
               write(filen,'(A,A,A)')'FILE.',TR(w_nm_str),'nm.POLARS.EL'
            end if
         else
            if(isNM)then
               write(filen,'(A,A,A)')'FILE.',TR(w_nm_str),'nm.POLARS.Q.EL'
            else
               write(filen,'(A,A,A)')'FILE.',TR(w_nm_str),'nm.POLARS'
            end if
         end if
      end if
      if(createFile)open(unitt,file=filen)
      write(unitt,*)nq,e00*au_2_cm
      v1_n=0
      v3_n=1
      do i = 1,nq
         v1=[0]
         v1_pos=[0]
         w1=0
         
         if(isNM)then
            v3=[1]
            v3_pos=[i]
            w3=wg(i)*au_2_cm
         else
            v3=[0]
            v3_pos=[0]
            w3=0
            v3_n=0
         end if
         
         call WriteTenPretty(unitt,wexc,polars(i),v1,v1_pos,v1_n,w1,v3,v3_pos,v3_n,w3)
      end do
      deallocate(v1,v3,v1_pos,v3_pos)
      if(createFile)close(unitt)
   end subroutine WritePolars
   
   
   subroutine WriteTenPretty(unitt,wexc,polarr,v1,v1_pos,v1_n,w1,v3,v3_pos,v3_n,w3)
      integer unitt,iexc,i,j,k
      integer v1_n,v3_n
      integer v1(v1_n)
      integer v3(v3_n)
      integer(int16) v3_pos(v3_n),v1_pos(v1_n)
      type(Polar) polarr
      double precision w1,w3,wexc
      
      if(wexc>1d-40)then
         write(unitt,'(1X,F8.2)')1d7/(wexc*au_2_cm)
      else
         write(unitt,'(1X,F8.2)')0d0
      end if
      write(unitt,*)v1_n
      write(unitt,*)v1,v1_pos
      write(unitt,*)v3_n
      write(unitt,*)v3,v3_pos
      write(unitt,'(3(1X,F8.3)," cm-1")')w1,w3,w3-w1
      
2000  format(3(1X,E15.8,1X,E15.8,'*i'))
2001  format(9(1X,E15.8,1X,E15.8,'*i'))      
      !Alpha, electric dipole-electric dipole
      write(unitt,'(A)')'Alpha, El.dip.-El.dip.'
      write(unitt,2000)((polarr%Ap(i,j),j=1,3),i=1,3)
      
      !G, electric dipole-magnetic dipole
      write(unitt,'(A)')'G, El.dip.-Mag.dip.'
      write(unitt,2000)((polarr%G(i,j),j=1,3),i=1,3)
      
      !G-cursive, magnetic dipole-electric dipole
      write(unitt,'(A)')'Gc, Mag.dip.-El.dip.'
      write(unitt,2000)((polarr%Gc(i,j),j=1,3),i=1,3)
      
      !A, electric dipole-electric quadrupole
      write(unitt,'(A)')'A, El.dip.-El.Quad.'
      write(unitt,2001)(((polarr%A(i,j,k),k=1,3),j=1,3),i=1,3)
      
      !A-cursive, electric quadrupole-electric dipole
      write(unitt,'(A)')'Ac, El.Quad.-El.dip.'
      write(unitt,2001)(((polarr%Ac(i,j,k),k=1,3),j=1,3),i=1,3)
      write(unitt,*)'--- ---'
      flush(unitt)
   end subroutine WriteTenPretty
   
   subroutine ReadFilettt(fn,polars,n,wexc)
      character(*) fn
      character(80) s80
      integer nat,n,i,aa,bb,ic,cc
      type(Polar),allocatable :: polars(:)
      double precision,allocatable :: ALPHA(:,:,:),G(:,:,:),A(:,:,:,:)
      double precision,allocatable :: ALPHAi(:,:,:),Gi(:,:,:),Ai(:,:,:,:)
      double precision wexc,wexc2
      
      
      open(2,file=fn,status='old')
      READ(2,*)
      READ(2,'(A80)')s80
      read(s80(2:4),*)nat
      read(s80(22:31),*)wexc
      close(2)
      n=3*nat
      
      allocate(polars(n))
      
      allocate(ALPHA(N,3,3),A(N,3,3,3),G(N,3,3))
      ic=0
      CALL READTEN_ttt(fn,N,NAT,ALPHA,G,A,ic)
      call Prim2TrA(N,A)
      allocate(ALPHAi(N,3,3),Ai(N,3,3,3),Gi(N,3,3))
      ic=3
      CALL READTEN_ttt(fn,N,NAT,ALPHAi,Gi,Ai,ic)
      if(ic==0)then
         Alphai=0d0
         Ai=0d0
         gi=0d0
      else
         call Prim2TrA(N,Ai)
      end if
      if(wexc==0d0)then !static limit
         wexc2=-1 
      else
         wexc2=wexc
      end if
      
      do i = 1,n
         do aa = 1,3
            do bb = 1,3
               polars(i)%ap(aa,bb)=alpha(i,aa,bb)+iu*alphai(i,aa,bb)
               if(ic==0)then
                  polars(i)%G(aa,bb)=(iu*G(i,aa,bb))*(-wexc2)
               else
                  polars(i)%G(aa,bb)=(G(i,aa,bb)+iu*Gi(i,aa,bb))*(-wexc2)
               end if
               do cc = 1,3
                  polars(i)%A(aa,bb,cc)=A(i,aa,cc,bb)+iu*Ai(i,aa,cc,bb)
               end do
            end do
         end do
      end do
      
      do i=1,n
         do aa = 1,3
            do bb = 1,3
               polars(i)%Gc(aa,bb)=-polars(i)%G(bb,aa)
            end do
         end do
         polars(i)%Ac=polars(i)%A
      end do
      
      deallocate(ALPHA,A,G,ALPHAi,Ai,Gi)
   end subroutine ReadFilettt
   
   subroutine Prim2TrA(nq,Aten)
      integer nq,i,a,b,c
      double precision,intent(inout) :: Aten(nq,3,3,3)
      double precision Abuf(3,3,3),trace_a
      
      do i = 1,nq
         Abuf=0d0
         do a = 1,3
            trace_a=Aten(i,a,1,1)+Aten(i,a,2,2)+Aten(i,a,3,3)
            do b = 1,3
               do c = 1,3
                  Abuf(a,b,c)=0.5d0*(3*Aten(i,a,b,c)-trace_a*KD(b,c))
               end do
            end do
         end do
         Aten(i,:,:,:)=Abuf
      end do
   end subroutine Prim2TrA
   
   SUBROUTINE READTEN_ttt(fn,N,NAT,P,G,A,ic)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      character(*) fn
      DIMENSION P(N,3,3),G(N,3,3),A(N,3,3,3)
      
      OPEN(2,FILE=fn,STATUS='OLD')
      READ(2,*)
      READ(2,*)NAT
      if(NAT.gt.N/3)stop 666
      READ(2,*)
      READ(2,*)
      DO 1 I=1,3
      READ(2,*)
      DO 1 L=1,NAT
      DO 1 IX=1,3
      M=3*(L-1)+IX
1     READ(2,*,err=2000)(P(M,I,J),J=1,2),(P(M,I,J),J=1,ic),(P(M,I,J),J=1,3)
      READ(2,*)
      READ(2,*)
      DO 2 I=1,3
      READ(2,*)
      DO 2 L=1,NAT
      DO 2 IX=1,3
      M=3*(L-1)+IX
!     last index magnetic, inner
2     READ(2,*)(G(M,I,J),J=1,2),(G(M,I,J),J=1,ic),(G(M,I,J),J=1,3)
      READ(2,*)
      READ(2,*)
      DO 3 I=1,3
      DO 3 J=1,3
      READ(2,*)
      DO 3 L=1,NAT
      DO 3 IX=1,3
      M=3*(L-1)+IX
3     READ(2,*)(A(M,I,J,K),K=1,2),(A(M,I,J,K),K=1,ic), &
      (A(M,I,J,K),K=1,3)
      CLOSE(2)
      WRITE(6,*)NAT,' atoms, tensors read in from FILE.TTT'
      RETURN
2000  ic=0
      close(2)
   END SUBROUTINE READTEN_ttt
   
   subroutine WriteFilettt(fn,polars,n3,nat,wexc)
      character(*) fn
      integer nat,n3,i,aa,bb,cc,idx,idxA
      type(Polar) :: polars(n3)
      double precision,allocatable :: ALPHA(:,:,:),G(:,:,:),A(:,:,:)
      double precision wexc,wexc2
      
      allocate(Alpha(1,9*n3,2),G(1,9*n3,2),A(1,27*n3,2))
      
      if(wexc==0d0)then
         wexc2=-1
      else
         wexc2=wexc
      end if
      do i = 1,n3
         do aa = 1,3
            do bb = 1,3
               idx=aa+3*(bb-1)+9*(i-1) !TODO check, taken from gar9_cheese.f, this hurts me physically
               Alpha(1,idx,1)=realpart(polars(i)%ap(aa,bb))
               Alpha(1,idx,2)=imagpart(polars(i)%ap(aa,bb))
               G(1,idx,1)=realpart(polars(i)%G(aa,bb))/(-wexc2)
               G(1,idx,2)=imagpart(polars(i)%G(aa,bb))/(-wexc2)
               
               do cc = 1,3
                  idxA=3*(bb-1)+9*(aa-1)+27*(i-1)+cc
                  A(1,idxA,1)=realpart(polars(i)%A(aa,cc,bb))
                  A(1,idxA,2)=imagpart(polars(i)%A(aa,cc,bb))
               end do
            end do
         end do
      end do
      
      call writettt09(fn,nat,alpha,G,A,n3,1,.true.,wexc)
      
      deallocate(Alpha,G,A)
   end subroutine WriteFilettt
   
   SUBROUTINE WRITETTT09(filename,NAT,ALPHA09,GTENS09,ATENS09,MX3,MFR,lcompl,wexc)
      IMPLICIT none
      INTEGER*4 MX3,ifr,MFR,I,J,K,L,IX,IIND,IM,IA,N,NQ,NAT,IQ,N7,ii
      real*8 ALPHA09(MFR,9*MX3,2),GTENS09(MFR,9*MX3,2), &
      ATENS09(MFR,27*MX3,2),amu,a,A0(3,3),G0(3,3),AT0(3,3,3),fr,p,wexc
      character(*) filename
      character*10 fchar
      real*8,allocatable::e(:),s(:,:),ALPHAQ(:,:,:),GTENSQ(:,:,:), &
      ATENSQ(:,:,:,:),ALPHA(:,:,:),GTENS(:,:,:),ATENS(:,:,:,:)
      integer*4,allocatable::nml(:)
      logical lex,lcompl

      fr=wexc
      !p=3.0d0/2.0d0
      p=1d0
      ifr=1
      !write(fchar,'(i10)')ifr
      !filename='FILE.TTT.f'//fchar(10:10)
      !if(ifr.gt.9)filename='FILE.TTT.f'//fchar(9:10)
      !if(ifr.gt.99)filename='FILE.TTT.f'//fchar(8:10)
      !if(ifr.gt.999)filename='FILE.TTT.f'//fchar(7:10)
      write(*,*)
      write(*,*)filename
      OPEN(2,FILE=filename)
!     for the first frequency, write the file twice:
      if(ifr.eq.1)OPEN(21,FILE=filename)
      WRITE(2,2000)NAT,ifr,fr
2000  FORMAT(' ROA tensors, cartesian derivatives',/, &
      I4,' atoms, freq. ',i2,f11.6/, &
      ' The electric-dipolar electric-dipolar polarizability:',/, &
      ' Atom/x    jx           jy           jz')
      if(ifr.eq.1)WRITE(21,2000)NAT,ifr,fr
      DO 1 I=1,3
      WRITE(2,2002)I
      if(ifr.eq.1)WRITE(21,2002)I
2002  FORMAT(' Alpha(',I1,',J):')
      DO 1 L=1,NAT
      DO 1 IX=1,3
      IIND=3*(L-1)+IX
      if(lcompl)then
       WRITE(2 ,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3), &
                          (ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),2),J=1,3)
       if(ifr.eq.1) &
       WRITE(21,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3), &
                          (ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),2),J=1,3)
      else
       WRITE(2 ,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3)
       if(ifr.eq.1) &
       WRITE(21,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3)
      endif
2001  FORMAT(I5,1H ,I1,6g15.7)
1     continue
      WRITE(2,2004)
      if(ifr.eq.1)WRITE(21,2004)
2004  FORMAT(' The electric dipole magnetic dipole polarizability:',/, &
             ' Atom/x    jx(Bx)       jy(By)       jz(Bz)')
      DO 2 I=1,3
      WRITE(2,2003)I
      if(ifr.eq.1)WRITE(21,2003)I
2003  FORMAT(' G(',I1,',J):')
      DO 2 L=1,NAT
      DO 2 IX=1,3
      IIND=3*(L-1)+IX
      if(lcompl)then
       WRITE(2 ,2001)L,IX,(GTENS09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3), &
                          (GTENS09(ifr,I+3*(J-1)+9*(IIND-1),2),J=1,3)
       if(ifr.eq.1) &
       WRITE(21,2001)L,IX,(GTENS09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3), &
                          (GTENS09(ifr,I+3*(J-1)+9*(IIND-1),2),J=1,3)
      else
       WRITE(2 ,2001)L,IX,(GTENS09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3)
       if(ifr.eq.1) &
       WRITE(21,2001)L,IX,(GTENS09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3)
      endif
2     continue
      WRITE(2,2005)
      if(ifr.eq.1)WRITE(21,2005)
2005  FORMAT(' The electric dipole electric quadrupole polarizability:', &
      /,     ' Atom/x    kx           ky           kz')
      DO 3 I=1,3
      DO 3 J=1,3
      WRITE(2,2006)I,J
      if(ifr.eq.1)WRITE(21,2006)I,J
2006  FORMAT(' A(',I1,',',I1,',K):')
      DO 3 L=1,NAT
      DO 3 IX=1,3
      ii=3*(J-1)+9*(I-1)+27*(3*(L-1)+IX-1)
      if(lcompl)then
       WRITE(2,20077)L,IX, &
       (ATENS09(ifr,K+ii,1)*p,K=1,3), &
       (ATENS09(ifr,K+ii,2)*p,K=1,3),L,IX,I,J
       if(ifr.eq.1)WRITE(21,20077)L,IX, &
       (ATENS09(ifr,K+ii,1)*p,K=1,3), &
       (ATENS09(ifr,K+ii,2)*p,K=1,3),L,IX,I,J
20077  FORMAT(I5,1H ,I1,6g15.7,' ',4i3)
      else
       WRITE(2,2007)L,IX, &
       (ATENS09(ifr,K+ii,1)*p,K=1,3),L,IX,I,J
       if(ifr.eq.1)WRITE(21,2007)L,IX, &
       (ATENS09(ifr,K+ii,1)*p,K=1,3),L,IX,I,J
2007   FORMAT(I5,1H ,I1,3g15.7,' ',4i3)
      endif
3     continue
      write(2,*)
      if(ifr.eq.1)write(21,*)
      write(2,*)'dummy alpha v:'
      if(ifr.eq.1)write(21,*)'dummy alpha v:'
      DO 4 I=1,3
      WRITE(2,2002)I
      if(ifr.eq.1)WRITE(21,2002)I
      DO 4 L=1,NAT
      DO 4 IX=1,3
      IIND=3*(L-1)+IX
      if(lcompl)then
       WRITE(2 ,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3), &
                          (ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),2),J=1,3)
       if(ifr.eq.1) &
       WRITE(21,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3), &
                          (ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),2),J=1,3)
      else
       WRITE(2 ,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3)
       if(ifr.eq.1) &
       WRITE(21,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1),1),J=1,3)
      endif
4     continue
      CLOSE(2)
      if(ifr.eq.1)CLOSE(21)

      RETURN
   END
   
   subroutine Readfileqttt(fn,polars,velocity,E,nq,wexc)
      integer nq,unitt,im
      integer i,ii,j,idx,aa,bb,K
      logical velocity
      character(*) fn
      type(Polar),allocatable :: polars(:)
      double precision,allocatable :: E(:)
      double precision alpha(3,3),alphai(3,3)
      double precision G(3,3),Gi(3,3)
      double precision A(3,3,3),Ai(3,3,3),wexc
      double precision :: alpha_vel(3,3),alpha_veli(3,3)
      
      unitt=22
      open(unitt,file=fn,status='old')
      read(unitt,*)
      read(unitt,2003,err=1111)nq,wexc
      
      goto 2222
1111  continue
      backspace(unitt)
      read(unitt,1999)nq
2222  continue
      
      read(unitt,*)
      read(unitt,*)
      allocate(polars(nq))
      do i=1,nq
         polars(i)%ap=(0d0,0d0)
         polars(i)%G=(0d0,0d0)
         polars(i)%Gc=(0d0,0d0)
         polars(i)%A=(0d0,0d0)
         polars(i)%Ac=(0d0,0d0)
      end do
2003  FORMAT(I4,7X,F15.8)
1999  FORMAT(I4,10X)
      
      allocate(E(nq))
      do ii = 1,nq
       read(unitt,221)im,E(nq-ii+1)
221    format(i5,f12.2)
       DO I=1,3
         read(unitt,2001)im,(ALPHA(I,J),J=1,3),(ALPHAi(I,J),J=1,3)
       end do
2001   FORMAT(I3,6E14.6)
       polars(nq-ii+1)%Ap=Alpha+Alphai*iu
      end do
      
      read(unitt,*)
      read(unitt,*)
      
      do ii = 1,nq
       read(unitt,221)im,E(nq-ii+1)
       DO I=1,3
          read(unitt,2001)im,(G(I,J),J=1,3),(Gi(I,J),J=1,3)
       end do
       polars(nq-ii+1)%G=(G*iu+Gi)*(-wexc)
      end do
      
      read(unitt,*)
      read(unitt,*)
      
      do ii = 1,nq
       read(unitt,221)im,E(nq-ii+1)
       DO I=1,3
          DO J=1,3
             read(unitt,2002)im,idx,(A(I,J,K),K=1,3),(Ai(I,J,K),K=1,3)
          end do
       end do
2002   FORMAT(2I3,6F14.6)
       polars(nq-ii+1)%A=A+Ai*iu
      end do
      
      do i=1,nq
         do aa = 1,3
            do bb = 1,3
               polars(i)%Gc(aa,bb)=-polars(i)%G(bb,aa)
            end do
         end do
         polars(i)%Ac=polars(i)%A
      end do
      
      if(.not.velocity)then
         close(unitt)
         return
      end if
      read(unitt,*)
      read(unitt,*)
      im=0
      DO II=1,nq
       !im=im+1
       read(unitt,221)im,E(nq-ii+1)
       DO I=1,3
         read(unitt,2001)im,(alpha_vel(I,J),J=1,3),(alpha_veli(I,J),J=1,3)
       end do
       polars(nq-ii+1)%ap=alpha_vel+alpha_veli*iu
      end do
      close(unitt)
   end subroutine Readfileqttt
   
   subroutine writefileqttt(fn,polars,nq,wg,wexc)
      character(*) fn
      integer nq,i,aa,bb,cc
      double precision wg(nq),wg_new(nq),wexc
      logical lcompl
      type(Polar) polars(nq)
      double precision, allocatable :: Alpha(:,:,:),Alphai(:,:,:)
      double precision, allocatable :: G(:,:,:),Gi(:,:,:)
      double precision, allocatable :: A(:,:,:,:),Ai(:,:,:,:)
      double precision, allocatable :: av(:,:,:),avi(:,:,:)
      
      allocate(Alpha(nq,3,3),Alphai(nq,3,3),G(nq,3,3),Gi(nq,3,3),A(nq,3,3,3),Ai(nq,3,3,3))
      allocate(av(nq,3,3),avi(nq,3,3))
      
      do i = 1,nq
         wg_new(nq-i+1)=wg(i)
         do aa=1,3
            do bb=1,3
               Alpha(nq-i+1,aa,bb)=realpart(polars(i)%ap(aa,bb))
               Alphai(nq-i+1,aa,bb)=imagpart(polars(i)%ap(aa,bb))
               
               G(nq-i+1,aa,bb)=imagpart(polars(i)%G(aa,bb))/(-wexc)
               Gi(nq-i+1,aa,bb)=realpart(polars(i)%G(aa,bb))/(-wexc)
               do cc=1,3
                  A(nq-i+1,aa,bb,cc)=realpart(polars(i)%A(aa,bb,cc))
                  Ai(nq-i+1,aa,bb,cc)=imagpart(polars(i)%A(aa,bb,cc))
               end do
            end do
         end do
      end do
      av=0d0
      avi=0d0
      lcompl=.true.
      call wtqi(fn,nq,1,nq,alpha,a,g,alphai,ai,gi,av,wg_new,0d0,40000d0,lcompl,wexc)
      deallocate(Alpha,Alphai,G,Gi,A,Ai,av,avi)
   end subroutine writefileqttt
   
   SUBROUTINE WTQi(fn,N,N7,NINT,ALPHA,A,G,ALPHAi,Ai,Gi,av,E,WMIN,WMAX,lcompl,wexc)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      character(*) fn
      logical lcompl
      real*8 wexc
      DIMENSION ALPHA(N,3,3),G(N,3,3),A(N,3,3,3), &
      ALPHAi(N,3,3),Gi(N,3,3),Ai(N,3,3,3),av(N,3,3),E(*)
!
      CM=219470.0d0
      im=0
      DO 4 II=N7,NINT
4     if(E(II)*CM.GE.WMIN.AND.E(II).LE.WMAX)im=im+1
      open(22,file=fn)
      WRITE(22,2003)im,wexc
2003  FORMAT(' ROA tensors, normal modes derivatives',/,I4,' modes ',F15.8,/, &
      ' The electric-dipolar electric-dipolar polarizability:',/, &
              ' mode e(cm-1)   jx           jy           jz')
      im=0
      DO 1 II=N7,NINT
      if(E(II)*CM.GE.WMIN.AND.E(II).LE.WMAX)then
       im=im+1
       write(22,221)im,E(II)*CM
221    format(i5,f12.2)
       DO I=1,3
         if(lcompl)then
            WRITE(22,2001)I,(ALPHA(II,I,J),J=1,3),(ALPHAi(II,I,J),J=1,3)
         else
            WRITE(22,2000)I,(ALPHA(II,I,J),J=1,3)
         end if
       end do
2001   FORMAT(I3,6E14.6)
2000   FORMAT(I3,3E14.6)
      endif
1     continue

      WRITE(22,2004)
2004  FORMAT(' The electric dipole magnetic dipole polarizability:',/, &
            ' mode  e(cm-1)  jx(Bx)       jy(By)       jz(Bz)')
      im=0
      DO 2 II=N7,NINT
      if(E(II)*CM.GE.WMIN.AND.E(II).LE.WMAX)then
       im=im+1
       write(22,221)im,E(II)*CM
       DO I=1,3
         if(lcompl)then
            WRITE(22,2001)I,(G(II,I,J),J=1,3),(Gi(II,I,J),J=1,3)
         else
            WRITE(22,2000)I,(G(II,I,J),J=1,3)
         end if
       end do
!      last inner index magnetic
      endif
2     continue

      WRITE(22,2005)
2005  FORMAT(' The electric dipole electric quadrupole polarizability:', &
      /,     ' mode  e(cm-1)  kx           ky           kz')
      im=0
      DO 3 II=N7,NINT
      if(E(II)*CM.GE.WMIN.AND.E(II).LE.WMAX)then
       im=im+1
       write(22,221)im,E(II)*CM
       DO I=1,3
       DO J=1,3
         if(lcompl)then
            WRITE(22,2002)I,J,(A(II,I,J,K),K=1,3),(Ai(II,I,J,K),K=1,3)
         else
            WRITE(22,1999)I,J,(A(II,I,J,K),K=1,3)
         end if
       end do
       end do
2002   FORMAT(2I3,6F14.6)
1999   FORMAT(2I3,3F14.6)
      endif
3     continue

      WRITE(22,2009)
2009  FORMAT('The velocity form of alpha',/, &
      ' mode e(cm-1)    vx           vy           vz')
      im=0
      DO 5 II=N7,NINT
      if(E(II)*CM.GE.WMIN.AND.E(II).LE.WMAX)then
       im=im+1
       write(22,221)im,E(II)*CM
       DO I=1,3
         if(lcompl)then
            WRITE(22,2001)I,(av(II,I,J),J=1,3),(av(II,I,J),J=1,3)
         else
            WRITE(22,2000)I,(av(II,I,J),J=1,3)
         end if
       end do
      endif
5     continue

      close(22)
      RETURN
   END
   
   
   function Equal_I2Arr(arr1,arr2,n)result(res)
      integer n,i
      integer(int16)arr1(n),arr2(n)
      logical res
      
      
      if(n==0)then
         res=.true.
         return
      end if
      
      do i = 1,n
         if(arr1(i)/=arr2(i))then
            res=.false.
            return
         end if
      end do
      res=.true.
   end function Equal_I2Arr
   
   function is_orca_output(filename) result(is_orca)
       implicit none
       character(len=*), intent(in) :: filename
       logical :: is_orca
       
       integer :: iunit, ierr, i
       character(len=256) :: line

       is_orca = .false.

       open(newunit=iunit, file=filename, status='old', action='read')

       ! Limit the search to the first 50 lines to avoid I/O bottlenecks 
       ! on large files that are not ORCA outputs.
       do i = 1, 50
           read(iunit, '(A)') line
           
           ! Check for the signature ORCA banner
           if (index(line, '* O   R   C   A *') > 0) then
               is_orca = .true.
               exit
           end if
       end do

       close(iunit)
   end function is_orca_output   
   
   
   subroutine rd_casscf_orca(iunit,filename, n_trans,nat, energies,r,z, elec_dipole, mag_dipole)
        ! ---------------------------------------------------------------------
        ! Extracts transitions from the ground state from ORCA 6 CASSCF output.
        ! ---------------------------------------------------------------------
        character(len=*), intent(in) :: filename
        integer, intent(inout) :: n_trans
        
        ! Allocatable arrays. Allocation happens dynamically inside the routine.
        real(8), allocatable, intent(out) :: energies(:)      ! Units: eV
        real(8), allocatable, intent(out) :: elec_dipole(:,:) ! Units: a.u.
        real(8), allocatable, intent(out) :: mag_dipole(:,:)  ! Units: a.u.
        real(8), allocatable, intent(out) :: r(:)  ! Units: angstrom
        integer,allocatable :: z(:)
        
        integer :: iunit, ierr, i,nat
        character(len=512) :: line
        character(len=32)  :: dummy(7)
        character(2) :: sym
        logical :: found_abs, found_cd
        
        
        found_abs = .false.
        found_cd = .false.
        
        open(iunit, file=filename, status='old', action='read')
        
        
        ! Pass 2: Rewind to precisely extract the numerical values
        
        do 
         read(iunit,'(a)')line
         if(index(line,'Number of atoms')>0)then
            read(line(48:54),*)nat
            exit
         end if
        end do
        allocate(r(3*nat),z(nat))
        
        rewind(iunit)
        do 
         read(iunit,'(a)')line
         if(index(line,'CARTESIAN COORDINATES (ANGSTROEM)')>0)then
            read(iunit,*)
            do i = 1,nat
               read(iunit,*)sym,r((i-1)*3+1),r((i-1)*3+2),r((i-1)*3+3)
               call To_upper(sym)
               z(i)=findloc(elements,sym,dim=1)
            end do
            exit
         end if
        end do
        
        
        if(n_trans<=0)then
           n_trans = 0
           ! Pass 1: Locate the electric dipole block to determine the number of transitions
           do
               read(iunit, '(A)', iostat=ierr) line
               if (ierr /= 0) exit
               
               ! The absorption spectrum table defaults to transitions from the ground state
               if (index(line, 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS') > 0) then
                   found_abs = .true.
                   ! Skip the 4 table header lines
                   read(iunit, '(A)', iostat=ierr) line
                   read(iunit, '(A)', iostat=ierr) line
                   read(iunit, '(A)', iostat=ierr) line
                   read(iunit, '(A)', iostat=ierr) line
                   
                   ! Count transitions until a separator or blank line is encountered
                   do
                       read(iunit, '(A)', iostat=ierr) line
                       if (ierr /= 0) exit
                       if (len_trim(line) == 0 .or. index(line, '----') > 0) exit
                       n_trans = n_trans + 1
                   end do
                   exit
               end if
           end do
           
           if (n_trans == 0 .or. .not. found_abs) then
               print *, "Warning: No transitions found or missing absorption spectrum block."
               close(iunit)
               return
           end if
        end if
        ! n_trans=n_trans-1
        ! Allocate the arrays based on the detected root count
        allocate(energies(n_trans))
        allocate(elec_dipole(3, n_trans))
        allocate(mag_dipole(3, n_trans))
        
        energies = 0.0d0
        elec_dipole = 0.0d0
        mag_dipole = 0.0d0
        
        rewind(iunit)
        ! Extract Energies and Electric Dipole Moments
        do
            read(iunit, '(A)',end=30) line
            
            if (index(line, 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS') > 0) then
                ! Skip the 4 header lines
                read(iunit, '(A)') line
                read(iunit, '(A)') line
                read(iunit, '(A)') line
                read(iunit, '(A)') line
                
                do i = 1, n_trans
                    ! Target format: State_Init -> State_Final Energy_eV Energy_cm-1 Wavelength_nm fosc D2 DX DY DZ
                    read(iunit, *) dummy(1:3), energies(i), dummy(4:7), elec_dipole(1, i), elec_dipole(2, i), elec_dipole(3, i)
                    energies(i)=energies(i)/27.2110d0
                end do
                
            end if
        
           ! Extract Magnetic Dipole Moments (continues reading from the current file position)
            
            if (index(line, 'CD SPECTRUM') > 0 .and. index(line, 'MOMENTS') > 0) then
                found_cd = .true.
                ! Skip the 4 header lines
                read(iunit, '(A)') line
                read(iunit, '(A)') line
                read(iunit, '(A)') line
                read(iunit, '(A)') line
                
                do i = 1, n_trans
                    ! Target format: State_Init -> State_Final Energy_eV Energy_cm-1 Wavelength_nm R MX MY MZ
                    read(iunit, *) dummy(1:7), mag_dipole(1, i), mag_dipole(2, i), mag_dipole(3, i)
                end do
                
            end if
        end do
        
30      if (.not. found_cd) then
            print *, "Warning: CD SPECTRUM block not found. Magnetic dipoles set to zero."
        end if
        
        close(iunit)
    end subroutine rd_casscf_orca
   
   
   subroutine rd_td_new(unitt,filename,nstates,z,r,nat,u,v,m,q,e_gr,ens,mult_tm,dnst,primq)
      character(*) filename
      character(80),allocatable :: splitt(:)
      integer unitt,nstates,nat,i,ii,sizee,j,k,l,idx,istate,dnst,nstates_r,td,bufi
      integer,parameter :: maxat=5000
      integer,parameter :: maxtd=2000
      integer,parameter :: max_tdc=100
      double precision,allocatable :: u(:,:),v(:,:),m(:,:),q(:,:,:)
      double precision,allocatable :: r(:),r_h(:),ens(:),ens_h(:)
      double precision bufr
      logical primq
      
      integer tdc_cur_n,orb_gr_idx,orb_ex_idx,tdcs_n,n,num
      integer,allocatable :: order(:)
      !type(TD_coeff),allocatable :: tdc_cur(:)
      ! type(TD_coeff_arr),allocatable :: tdcs_h(:),tdcs(:),tdcs_buf(:)
      character(7) orb_gr,orb_ex
      character(3) orb_gr_sym,orb_ex_sym
      
      
      
      double precision e,qq(6),qqq(6),mult_tm(5)
      integer,allocatable :: z(:),z_h(:)
      double precision e_gr,bl(5,16),vall
      character(100) s80,s80_2,s_arr(9),buf
      character(8) bufc
      logical endd
      
      allocate(z_h(maxat),r_h(3*maxat),ens_h(maxtd))
      ! allocate(tdc_cur(max_tdc),tdcs_h(maxtd))
      z_h=0
      r_h=0
      istate=1
      nstates=0
      nstates_r=0
      open(unitt,file=filename,status='old')
30    read(unitt,'(A100)',end=40)s80
      
      if(s80(27:44)=='Input orientation:')then
         read(unitt,*)
         read(unitt,*)
         read(unitt,*)
         read(unitt,*)
         i=1
         read(unitt,'(A100)')s80_2
         do while(s80_2(2:4)/='---')
            idx=(i-1)*3+1
            read(s80_2,*)bufc,z_h(i),k,r_h(idx:idx+2)
            i=i+1
            read(unitt,'(A100)')s80_2
         end do
         nat=i-1
      elseif(s80(2:10)=='SCF Done:')then
         splitt=SplitString(s80,80,' ')
         read(splitt(5),*)e_gr
      elseif(s80(2:14)=='Excited State')then
         idx=index(s80,' nm')
         read(s80(idx-8:idx-1),*)e
         e=(1d7/e)*cm_2_au
         nstates_r=nstates_r+1
         ens_h(nstates_r)=e
      elseif(s80(2:31)=='Electronic transition elements')then
         if(dnst<=0)dnst=nstates_r
         nstates=dnst
         if(allocated(ens))deallocate(ens)
         ens=ens_h(1:dnst)
         ! if(allocated(tdcs))deallocate(tdcs)
         ! tdcs=tdcs_h(1:dnst)
         endd=.false.

         allocate(u(3,dnst),v(3,dnst),m(3,dnst),q(3,3,dnst))
         ii=0
         bl=ReadBlock(unitt,endd,sizee)
         outer : do while(.not.endd)
            ii=ii+1
            do i = 1,sizee
               td=(ii-1)*5+i
               if(td>dnst)exit outer
               u(:,td)=bl(i,2:4)
               v(:,td)=bl(i,5:7)
               m(:,td)=bl(i,8:10)
               qqq=bl(i,11:16)
               !qq(1),qq(4),qq(6),qq(2),qq(3),qq(5)
               qq(1)=qqq(1)
               qq(4)=qqq(2)
               qq(6)=qqq(3)
               qq(2)=qqq(4)
               qq(3)=qqq(5)
               qq(5)=qqq(6)
               q(:,:,td)=TracelessQ(qq)
               
            end do
            bl=ReadBlock(unitt,endd,sizee)
         end do outer
      end if
      
      goto 30
40    close(unitt)
      
      !the #p option was not specified, read the common td output instead
      if(.not.allocated(u))then
         open(unitt,file=filename,status='old')
50       read(unitt,'(A100)',end=60)s80
         if(s80(2:64)=='Ground to excited state transition electric dipole moments (Au)')then
            read(unitt,*)
            if(dnst<=0)dnst=nstates_r
            nstates=dnst
            allocate(u(3,dnst),v(3,dnst),m(3,dnst),q(3,3,dnst))
            if(allocated(ens))deallocate(ens)
            ens=ens_h(1:dnst)
            
            !electric dipole
            do i = 1,dnst
               read(unitt,*)bufI,u(1,i),u(2,i),u(3,i),bufr,bufr
            end do
            !skip unwanted TDMs
            do i = dnst+1,nstates_r
               read(unitt,*)
            end do
            read(unitt,*)
            read(unitt,*)
            
            !velocity dipole
            do i = 1,dnst
               read(unitt,*)bufI,v(1,i),v(2,i),v(3,i),bufr,bufr
            end do
            !skip unwanted TDMs
            do i = dnst+1,nstates_r
               read(unitt,*)
            end do
            read(unitt,*)
            read(unitt,*)
            
            !magnetic dipole
            do i = 1,dnst
               read(unitt,*)bufI,m(1,i),m(2,i),m(3,i)
            end do
            !skip unwanted TDMs
            do i = dnst+1,nstates_r
               read(unitt,*)
            end do
            read(unitt,*)
            read(unitt,*)
            
            !electric (mixed) quadrupole
            do i = 1,dnst
               read(unitt,*)bufI,qqq(1),qqq(2),qqq(3),qqq(4),qqq(5),qqq(6)
               qq(1)=qqq(1)
               qq(4)=qqq(2)
               qq(6)=qqq(3)
               qq(2)=qqq(4)
               qq(3)=qqq(5)
               qq(5)=qqq(6)
               if(primq)then
                  q(:,:,i)=PrimitiveQ(qq)
               else
                  q(:,:,i)=TracelessQ(qq)
               end if
            end do
            goto 60
         end if
         goto 50
60       close(unitt)
      end if
      
      z=z_h(1:nat)
      r=r_h(1:3*nat)
      deallocate(z_h,r_h)
      u=u*mult_tm(1)
      v=v*mult_tm(2)
      m=m*mult_tm(3)
      
      !deallocate(tdc_cur,tdcs_h)
      deallocate(ens_h)
   end subroutine rd_td_new
   
   subroutine rd_td_derivatives(unitt,filename,z,r,nat,u,v,m,q,du,dv,dm,dq,e_gr,e_tr,iroot,mult_tm)
      character(*) filename
      character(80),allocatable :: splitt(:)
      integer unitt,nat,n3,i,ii,sizee,j,k,l,idx,istate,td,el_c_a,el_c_b,el_c
      integer,parameter :: maxat=5000
      integer,parameter :: maxtd=2000
      integer,parameter :: max_tdc=100
      integer :: iroot
      double precision :: u(3),v(3),m(3),q(3,3),e_tr
      double precision,allocatable :: r(:),r_h(:)
      double precision,allocatable :: du(:,:),dv(:,:),dm(:,:),dq(:,:,:)
      logical rootFound
      
      integer tdc_cur_n,orb_gr_idx,orb_ex_idx,tdcs_n,n,num
      integer,allocatable :: order(:)

      character(7) orb_gr,orb_ex
      character(3) orb_gr_sym,orb_ex_sym
      
      
      
      double precision e,qq(6),qqq(6),mult_tm(5)
      integer,allocatable :: z(:),z_h(:)
      double precision e_gr,bl(5,16),vall
      character(100) s80,s80_2,s_arr(9),buf
      character(8) bufc
      logical endd
      
      allocate(z_h(maxat),r_h(3*maxat))
      rootFound=.false.
      z_h=0
      r_h=0
      istate=1
      iroot=0
      open(unitt,file=filename,status='old')
30    read(unitt,'(A100)',end=40)s80
      
      if(s80(27:44)=='Input orientation:')then
         read(unitt,*)
         read(unitt,*)
         read(unitt,*)
         read(unitt,*)
         i=1
         read(unitt,'(A100)')s80_2
         do while(s80_2(2:4)/='---')
            idx=(i-1)*3+1
            read(s80_2,*)bufc,z_h(i),k,r_h(idx:idx+2)
            i=i+1
            read(unitt,'(A100)')s80_2
         end do
         nat=i-1
         n3=3*nat
      elseif(s80(2:10)=='SCF Done:')then
         splitt=SplitString(s80,80,' ')
         read(splitt(5),*)e_gr
      elseif(s80(2:14)=='Excited State' .and. .not.rootFound)then
10       iroot=iroot+1
         idx=index(s80,' nm')
         read(s80(idx-8:idx-1),*)e
         e=(1d7/e)*cm_2_au
         
         do
            read(unitt,'(A100)')s80
            if(s80(2:59)=='This state for optimization and/or second-order correction')then
               rootFound=.true.
               exit
            else if(s80(2:14)=='Excited State')then
               goto 10
            end if
         end do
         e_tr=e
      elseif(s80(2:34)=='Electronic Transition Derivatives')then
         endd=.false.
         if(.not.allocated(du))allocate(du(3,n3),dv(3,n3),dm(3,n3),dq(3,3,n3))
         ii=0
         bl=ReadBlock(unitt,endd,sizee)
         outer : do while(.not.endd)
            ii=ii+1
            do i = 1,sizee
               if(ii==1 .and. i<=3)cycle !skip first 3 columns
               td=(ii-1)*5+i-3
               du(:,td)=bl(i,2:4)
               dv(:,td)=bl(i,5:7)
               dm(:,td)=bl(i,8:10)
               qqq=bl(i,11:16)
               !qq(1),qq(4),qq(6),qq(2),qq(3),qq(5)
               qq(1)=qqq(1)
               qq(4)=qqq(2)
               qq(6)=qqq(3)
               qq(2)=qqq(4)
               qq(3)=qqq(5)
               qq(5)=qqq(6)
               dq(:,:,td)=TracelessQ(qq)
               
            end do
            bl=ReadBlock(unitt,endd,sizee)
         end do outer
         goto 40
      end if
      
      goto 30
40    close(unitt)
      
      z=z_h(1:nat)
      r=r_h(1:3*nat)
      deallocate(z_h,r_h)
      u=u*mult_tm(1)
      v=v*mult_tm(2)
      m=m*mult_tm(3)
      
   end subroutine rd_td_derivatives
   
   subroutine TDOrbitalStrExtract(str,num,sym)
      character(*) str
      integer num,n,i
      character(3) sym
      
      n=len(str)
      sym='   '
      do i = n,1,-1
         if(IsStrNumber_I(sym(i:i)))exit
      end do
      read(str(1:i),'(I3)')num
      if(i==n)return
      read(str(i:n),'(A3)')sym
   end subroutine TDOrbitalStrExtract
   
   function IsStrNumber_I(str)result(res)
      character(*) str
      integer i,n
      logical res
      
      n=len(str)
      res=.false.
      if(index(str,'E')>0 .or. index(str,'e')>0)return
      if(index(str,'D')>0 .or. index(str,'d')>0)return
      if(index(str,'.')>0)return
      
      read(str,'(I8)',err=111)n
      res=.true.
      return
111   continue
   end function IsStrNumber_I
   
   function ReadBlock(unitt,endd,sizee)result(res)
      integer unitt,i,j,sizee
      logical endd
      double precision res(5,16)
      character(80) s80
      integer header(5)
      
      res=0d0
      sizee=0
      header=0
      read(unitt,'(A80)')s80
      read(s80,*,end=2000,err=2000)header
2000  continue
      sizee=count(header/=0)
      if(sizee==0)then
         endd=.true.
         return
      end if
      do i = 1,16
         read(unitt,*)j,res(1:sizee,i)
      end do
   end function ReadBlock
   
   subroutine Forward(unitt,n)
      integer unitt,n,i
      do i = 1,n
         read(unitt,*)
      end do
   end subroutine Forward
   
   !Taken from guvcde by Petr Bour
   subroutine rdd(f,d,nat,u0,m0,a,v,v0,q,q0,Line,nroot,eau,mult_tm,mult_dtm)
!     read transition dipole derivatives from excited state freq calc
      implicit none
      integer nat,ln,ia,xa,i,xar,iar,ibr,ib,ic,ii,ix,iwr,nd, &
                ntr,j,ir,idx, &
                iLeft,iRight,iTD,iFREQ,iNSTATES
      integer,intent(inout) :: Line !Given, skips number of lines; when returned gives last line read
      integer,intent(out) :: nroot
      double precision d(9*nat),u(6),u0(3),step,sign,stepau,m0(3),a(9*nat), &
             v(9*nat),v0(3),q(18*nat),q0(6),ev,eau
      logical lnumerical,isAH
      character*(*) f
      character*70 s80 !gaussian might change the output file row length for some stupid reason, modify this if they do
      character*200 s200
      character*70,dimension(10) :: s80_arr !Gaussian input parameters
      character*800 :: s800
      double precision,allocatable::v16(:,:)
      double precision mult_tm(5),mult_dtm(5)
!
!     total number of roots (exc. states) calculated:
      ntr=0
      iwr=0
      lnumerical=.true.
      open(9,file=f,status='old',action='READ')
      
      if(Line>0)then
         do i= 1,line
            read(9,*)
         end do
         ln=line
      else
         ln=0
      end if
      
      iar=0
      xar=0
      ibr=0 

!     where it starts:
1     read(9,900,end=88,err=88)s80
900   format(1x,a70)
      ln=ln+1
      if(index(s80,'#')>0 .and. index(s80,'\#')==0)then
         ic=0
         s80_arr(1)=s80
         idx=1
         read(9,900)s80
         ln=ln+1
         do while (index(s80,'-------')==0)
           s80_arr(idx+1)=s80
           idx=idx+1
           read(9,900)s80
           ln=ln+1
         end do
         s800=''
         do i=1,idx
           s800(len(s80)*(i-1)+1:len(s80)*i)=s80_arr(i)
         end do
         
         call To_Lower(s800)
         iFreq=index(s800,'freq')
         iTD=index(s800,'td')
         if(iTD > 0 .and. iFreq > 0)then
            iLeft=index(s800(iTD:),'(')+itd-1
            iRight=index(s800(iTD:),')')+itd-1
            s200=TR(s800(ileft:iright))
            
            iNSTATES=index(s200,'nstates=')
            iRight=findAnyChar(s200,[')',','],2)
            read(s200(iNSTATES+8:iright-1),*)ntr
            iNSTATES=index(s200,'root=')
            iRight=findAnyChar(s200,[')',','],2)
            read(s200(iNSTATES+5:iright-1),*)nroot
            goto 2
         end if
      end if
      goto 1
      
88    close(9)
      write(output_unit,*)'TD freq not found'
      stop

2     ia=1
      xa=1
      ib=0
      ii=0
      step=0.001d0
      a=0.0d0
      v=0.0d0
      q=0.0d0
      d=0.0d0
      nd=0
      
4     read(9,901,end=78)s200
901   format(A200)
      if(s200(2:14).eq.'Nuclear step=')then
       read(s200(15:23),*)step
      endif
      if(s200(2:24).eq.'Re-enter D2Numr: IAtom=')then
       read(s200(25:27),*)iar
       read(s200(34:34),*)xar
       read(s200(42:43),*)ibr
      endif

      if(s200(2:65).eq.'Ground to excited state transition electric dipole moments (Au):')then
       if(ii.eq.0)then
        do i=1,nroot
           read(9,*)
        enddo
        read(9,*)u0(1),(u0(ix),ix=1,3)
       else
        ib=ib+1
        if(ib.gt.2)then
         ib=1
         xa=xa+1
         if(xa.gt.3)then
          xa=1
          ia=ia+1
         endif
        endif
        if(iwr.gt.1)write(6,*)' d ',ia,xa,ib,iar,xar,ibr
        do i=1,nroot
           read(9,*)
        end do
        read(9,*)u(1),(u(ix),ix=1,3)
        if(ib.eq.1)then
         sign=1.0d0
        else
         sign=-1.0d0
        endif
        do ix=1,3
           d(ix+3*(xa-1)+9*(ia-1))=d(ix+3*(xa-1)+9*(ia-1))+sign*u(ix)
        enddo
       endif
      endif

      if(s200(2:65).eq.'Ground to excited state transition velocity dipole moments (Au):')then
       if(ii.eq.0)then
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)v0(1),(v0(ix),ix=1,3)
       else
        if(iwr.gt.1)write(6,*)' v ',ia,xa,ib,iar,xar,ibr
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)u(1),(u(ix),ix=1,3)
        do ix=1,3
        v(ix+3*(xa-1)+9*(ia-1))=v(ix+3*(xa-1)+9*(ia-1))+sign*u(ix)
        end do
       endif
      endif

      if(s200(2:65).eq.'Ground to excited state transition magnetic dipole moments (Au):')then
       if(ii.eq.0)then
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)m0(1),(m0(ix),ix=1,3)
       else
        if(iwr.gt.1)write(6,*)' m ',ia,xa,ib,iar,xar,ibr
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)u(1),(u(ix),ix=1,3)
        do ix=1,3
        a(ix+3*(xa-1)+9*(ia-1))=a(ix+3*(xa-1)+9*(ia-1))+sign*u(ix)
        end do
       endif
      endif

      if(s200(2:69).eq.'Ground to excited state transition velocity quadrupole moments (Au):')then
       if(ii.eq.0)then
        do i=1,nroot
        read(9,*)
        end do
!       our order         1    2      3     4     5     6
!       our order        xx    xy    xz    yy    yz    zz
!       gaussian   #     xx    yy    zz    xy    xz    yz
        read(9,*)q0(1),q0(1),q0(4),q0(6),q0(2),q0(3),q0(5)
       else
        if(iwr.gt.1)write(6,*)' m ',ia,xa,ib,iar,xar,ibr
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)u(1),u(1),u(4),u(6),u(2),u(3),u(5)
        do ix=1,3
        q(ix+6*(xa-1)+18*(ia-1))=q(ix+6*(xa-1)+18*(ia-1))+sign*q(ix)
        end do
       endif
      endif

      if(s200(2:15).eq.'Excited State '.and.s200(16:16).ne.'s')then
       if(ii.eq.0)then
        do i=1,len(s200)-1
        if(s200(i:i)  .eq.':') read(s200(15:i-1),*)nd
        if(s200(i:i+1).eq.'eV')read(s200(i-10:i-2),*)ev
        end do
        if(nd.eq.nroot)then
         eau=ev/27.211384205943d0
         if(ia.eq.nat.and.ix.eq. 3 .and.ib .eq. 2)goto 78
         ii=ii+1
        endif
       endif
      endif

!     analytical frequency output:
!     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      if(s200(2:31).eq.'Electronic transition elements')then
!      read nroot moments
       allocate(v16(16,nroot))
       call ru(9,v16,16,nroot)
       do ix=1,3
       u0(ix)=v16(ix+ 1,nroot)
       v0(ix)=v16(ix+ 4,nroot)
       m0(ix)=v16(ix+ 7,nroot)
       end do
       q0(1 )=v16(11,nroot)
       q0(2 )=v16(14,nroot)
       q0(3 )=v16(15,nroot)
       q0(4 )=v16(12,nroot)
       q0(5 )=v16(16,nroot)
       q0(6 )=v16(13,nroot)
       deallocate(v16)
      endif

      if(s200(2:34).eq.'Electronic Transition Derivatives')then
       allocate(v16(16,3 + 3*nat))
       call ru(9,v16,16,3 + 3*nat) !the first three are derivatives of the applied electrical field
       do ia=1,nat
       do xa=1,3
       ii=xa+3*(ia-1)
       do ix=1,3
       d(ix  +3*(ii-1))=v16(ix+ 1,3+ii)
       v(ix  +3*(ii-1))=v16(ix+ 4,3+ii)
       a(ix  +3*(ii-1))=v16(ix+ 7,3+ii)
       end do
       !           1  2  3  4  5  6
       !Gaussian: xx yy zz xy xz yz
       !  guvcde: xx xy xz yy yz zz
       q(1+6*(ii-1))=v16(11,3+ii)
       q(2+6*(ii-1))=v16(14,3+ii)
       q(3+6*(ii-1))=v16(15,3+ii)
       q(4+6*(ii-1))=v16(12,3+ii)
       q(5+6*(ii-1))=v16(16,3+ii)
       q(6+6*(ii-1))=v16(13,3+ii)
       end do
       end do
       deallocate(v16)
       lnumerical=.false.
      endif
!     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
!     checkpoint style
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(s200(1:18).eq.'ETran state values')then
       allocate(v16(1,16*ntr+48+16*3*nat))
       read(9,*)v16
       ii=0
       do ir=1,ntr
       ii=ii+1
       do ix=1,3
       ii=ii+1
       if(ir.eq.nroot)u0(ix)=v16(1,ii)
       end do
       do ix=1,3
       ii=ii+1
       if(ir.eq.nroot)v0(ix)=v16(1,ii)
       end do
       do ix=1,3
       ii=ii+1
       if(ir.eq.nroot)m0(ix)=v16(1,ii)
       end do
       ii=ii+1
       if(ir.eq.nroot)q0(1 )=v16(1,ii)
       ii=ii+1
       if(ir.eq.nroot)q0(4 )=v16(1,ii)
       ii=ii+1
       if(ir.eq.nroot)q0(6 )=v16(1,ii)
       ii=ii+1
       if(ir.eq.nroot)q0(2 )=v16(1,ii)
       ii=ii+1
       if(ir.eq.nroot)q0(3 )=v16(1,ii)
       ii=ii+1
       if(ir.eq.nroot)q0(5 )=v16(1,ii)
       end do
       ii=ii+48

       do j=1,3*nat
       ii=ii+1
       do ix=1,3
       ii=ii+1
       d(ix+3*(j-1))=v16(1,ii)
       end do
       do ix=1,3
       ii=ii+1
       v(ix+3*(j-1))=v16(1,ii)
       end do
       do ix=1,3
       ii=ii+1
       a(ix+3*(j-1))=v16(1,ii)
       end do
       ii=ii+1
       q(1+6*(j-1))=v16(1,ii)
       ii=ii+1
       q(4+6*(j-1))=v16(1,ii)
       ii=ii+1
       q(6+6*(j-1))=v16(1,ii)
       ii=ii+1
       q(2+6*(j-1))=v16(1,ii)
       ii=ii+1
       q(3+6*(j-1))=v16(1,ii)
       ii=ii+1
       q(5+6*(j-1))=v16(1,ii)
       end do
       deallocate(v16)
       lnumerical=.false.
      endif
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      goto 4

78    close(9)


      if(lnumerical)then
       stepau=2.0d0*step/0.529177d0
       do ix=1,9*nat
       d(ix      )=d(ix      )/stepau
       v(ix      )=v(ix      )/stepau
       a(ix      )=a(ix      )/stepau
       q(      ix)=q(      ix)/stepau
       q(9*nat+ix)=q(9*nat+ix)/stepau
       end do
      endif

      !v  = -v/eau !transformation from velocity, not to length, but just same units Im told
      !q  = -q/eau !transformation from velocity
      !a  =  a!/2.0d0 !transformation cause gradient
      ! v0 = v0/eau
      ! q0 = q0/eau
      u0 = u0*mult_tm(1)
      v0 = v0*mult_tm(2)
      m0 = m0*mult_tm(3)
      
      q0(1) = q0(1)*mult_tm(4)
      q0(4) = q0(4)*mult_tm(4)
      q0(6) = q0(6)*mult_tm(4)
      
      q0(2) = q0(2)*mult_tm(5)
      q0(3) = q0(3)*mult_tm(5)
      q0(5) = q0(5)*mult_tm(5)
      
      d=d*mult_dtm(1)
      v=v*mult_dtm(2)
      a=a*mult_dtm(3)
      
      do i = 1,3*nat
         ii=(i-1)*6+1
         q(ii)=q(ii)*mult_dtm(4)
         q(ii+3)=q(ii+3)*mult_dtm(4)
         q(ii+5)=q(ii+5)*mult_dtm(4)
         
         q(ii+1)=q(ii+1)*mult_dtm(5)
         q(ii+2)=q(ii+2)*mult_dtm(5)
         q(ii+4)=q(ii+4)*mult_dtm(5)
      end do
      
!     dipole derivatives from gaussian:
      open(40,file='DEG.TEN')
      write(40,*)'atom coord dipx        dipy        dipz'
      do ia=1,nat
      do xa=1,3
      write(40,400)ia,xa,(d(ix+3*(xa-1)+9*(ia-1)),ix=1,3)
      end do
      end do
      write(40,*)'equilibrium transition electric dipole:'
      write(40,400)0,0,(u0(ix),ix=1,3)
      write(40,*)'atom coordmdipx       mdipy       mdipz'
      do ia=1,nat
      do xa=1,3
      write(40,400)ia,xa,(a(ix+3*(xa-1)+9*(ia-1)),ix=1,3)
      end do
      end do
      write(40,*)'equilibrium transition magnetic dipole:'
      write(40,400)0,0,(m0(ix),ix=1,3)
      write(40,*)'atom coord Qv xx xy xz yy yz zz:'
      do ia=1,nat
      do xa=1,3
      write(40,400)ia,xa,(q(ix+6*(xa-1)+18*(ia-1)),ix=1,6)
      end do
      end do
      write(40,*)'equilibrium transition quadrupole velocity:'
      write(40,400)0,0,(q0(ix),ix=1,6)
      write(40,*)'atom coord grad x y z:'
      do ia=1,nat
      do xa=1,3
      write(40,400)ia,xa,(v(ix+3*(xa-1)+9*(ia-1)),ix=1,3)
      end do
      end do
      write(40,*)'equilibrium transition velocity:'
      write(40,400)0,0,(v0(ix),ix=1,3)
400   format(i5,i2,6g12.4)
      close(40)
      
      line=ln
      return
      
      contains
      
      subroutine ru(io,A,n,m)
!     matrix n x m
      implicit none
      DOUBLE PRECISION A(n,m)
      INTEGER io,n,m,N1,N3,LN,J
      N1=1
1     N3=MIN(N1+4,m)
      read(io,*)
      DO LN=1,N
         read(io,*)A(LN,N1),(A(LN,J),J=N1,N3)
      end do
      N1=N1+5
      IF(N3.LT.m)GOTO 1
      return
      end subroutine ru

      
   end subroutine rdd
   
   SUBROUTINE READFF(N,FCAR,file_fc)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FCAR(N,N)
      CHARACTER(*) file_fc
      OPEN(20,FILE=file_fc,STATUS='OLD')
      N1=1
1     N3=N1+4
      IF(N3.GT.N)N3=N
      DO 130 LN=N1,N
130   READ(20,17)(FCAR(LN,J),J=N1,MIN(LN,N3))
      N1=N1+5
      IF(N3.LT.N)GOTO 1
17    FORMAT(4X,5D14.6)
!
!      CONST=4.359828d0/0.5291772d0/0.5291772d0 ! conversion to whatever
      CONST=1
      DO 3 I=1,N
      DO 3 J=I,N
3     FCAR(J,I)=FCAR(J,I)*CONST
      DO 31 I=1,N
      DO 31 J=I+1,N
   31 FCAR(I,J)=FCAR(J,I)
      CLOSE(20)
      !WRITE(*,*)' Cartesian FF read in ... '
      RETURN
   END SUBROUTINE READFF
   
   subroutine WriteGrad(unitt,fn,grad,n3)
      integer unitt,n3,ix,i
      character(*) fn
      double precision grad(n3)
      
      open(unitt,file=fn)
      !write(unitt,*)n3
      do i = 1,n3/3
         ix=(i-1)*3+1
         write(unitt,'(3F15.8)')grad(ix),grad(ix+1),grad(ix+2)
      end do
      close(unitt)
   end subroutine WriteGrad
   
   subroutine readgrad(nat,gr,filee)
      integer*4 nat,i,ix
      real*8 gr(3*nat)
      character(*) filee
      
      open(77,file=filee,status='old')
      do i = 1,nat
         ix=(i-1)*3+1
         read(77,*)gr(ix),gr(ix+1),gr(ix+2)
      end do
      close(77)
   end subroutine readgrad

   subroutine readgradQ(NQ,gr,filee)
      integer*4 NQ,i
      real*8 gr(NQ)
      character(*) filee
      
      open(77,file=filee,status='old')
      do i = 1,NQ
         read(77,*)gr(i)
      end do
      close(77)
   end subroutine readgradQ
   
   subroutine trafN_new(N,dd,ddi,nat,NQ,s)
      integer N,nat,NQ,ix,iq,ia,xa,iqi,ii
      double PRECISION dd(*),ddi(*),s(:,:)
      do ix=1,N
      do iq=1,NQ
      iqi=ix+N*(iq-1)
      ddi(iqi)=0.0d0
      do ia=1,nat
      do xa=1,3
      ii=xa+3*(ia-1)
!     s-opposite order of normal modes:
      ddi(iqi)=ddi(iqi)+dd(ix+N*(ii-1))*s(ii,iq)
      enddo
      enddo
      enddo
      enddo
      return
   end subroutine trafN_new
   
   function Car2NM_Pol1(n3,nq,iq,dpol_c,smat)result(dpol_q)
      integer n3,nq,a,b,c,i,iq
      double precision smat(n3,nq)
      type(Polar) dpol_c(n3),dpol_q
      
      dpol_q=0d0
      do i = 1,n3
         do a = 1,3
            do b = 1,3
               dpol_q%ap(a,b)=dpol_q%ap(a,b)+smat(i,iq)*dpol_c(i)%ap(a,b)
               dpol_q%G(a,b)=dpol_q%G(a,b)+smat(i,iq)*dpol_c(i)%G(a,b)
               dpol_q%Gc(a,b)=dpol_q%Gc(a,b)+smat(i,iq)*dpol_c(i)%Gc(a,b)
               do c = 1,3
                  dpol_q%A(a,b,c)=dpol_q%A(a,b,c)+smat(i,iq)*dpol_c(i)%A(a,b,c)
                  dpol_q%Ac(a,b,c)=dpol_q%Ac(a,b,c)+smat(i,iq)*dpol_c(i)%Ac(a,b,c)
               end do
            end do
         end do
      end do
      
   end function Car2NM_Pol1
   
   function NM2Car_Pol1(n3,nq,ic,dpol_q,smat)result(dpol_c)
      integer n3,nq,a,b,c,i,ic
      double precision smat(n3,nq)
      type(Polar) dpol_c,dpol_q(nq)
      
      dpol_c=0d0
      do i = 1,nq
         do a = 1,3
            do b = 1,3
               dpol_c%ap(a,b)=dpol_c%ap(a,b)+smat(ic,i)*dpol_q(i)%ap(a,b)
               dpol_c%G(a,b)= dpol_c%G(a,b)+smat(ic,i)*dpol_q(i)%G(a,b)
               dpol_c%Gc(a,b)=dpol_c%Gc(a,b)+smat(ic,i)*dpol_q(i)%Gc(a,b)
               do c = 1,3
                  dpol_c%A(a,b,c) =dpol_c%A(a,b,c)+smat(ic,i)*dpol_q(i)%A(a,b,c)
                  dpol_c%Ac(a,b,c)=dpol_c%Ac(a,b,c)+smat(ic,i)*dpol_q(i)%Ac(a,b,c)
               end do
            end do
         end do
      end do
      
   end function NM2Car_Pol1
   
   pure function Car2NM_FF(n3,nq,ff_c,smat)result(ff_q)
      integer,intent(in) :: n3,nq
      integer i,j,k,l
      double precision,intent(in) :: ff_c(n3,n3),smat(n3,nq)
      double precision :: ff_q(nq,nq)
      
      do i = 1,nq
         do j = 1,nq
            ff_q(j,i)=0d0
            do k = 1,n3
               do l = 1,n3
                  ff_q(j,i)=ff_q(j,i)+ff_c(k,l)*smat(k,i)*smat(l,j)
               end do
            end do
         end do
      end do
   end function Car2NM_FF
   
   
   function Car2NM_flat(n3,nq,dc,smat)result(dq)
      integer n3,nq,i,j
      double precision dc(n3),smat(n3,nq),dq(nq)
      
      do i = 1,nq
         dq(i)=0d0
         do j = 1,n3
            dq(i)=dq(i)+dc(j)*smat(j,i)
         end do
      end do
   end function Car2NM_flat
   
   function Car2NM(n3,nq,dc,smat)result(dq)
      integer n3,nq,i,j
      double precision dc(3,n3),smat(n3,nq),dq(3,nq)
      
      do i = 1,nq
         dq(:,i)=0d0
         do j = 1,n3
            dq(1,i)=dq(1,i)+dc(1,j)*smat(j,i)
            dq(2,i)=dq(2,i)+dc(2,j)*smat(j,i)
            dq(3,i)=dq(3,i)+dc(3,j)*smat(j,i)
         end do
      end do
   end function Car2NM
   
   function Car2NM_2(n3,nq,dc,smat)result(dq)
      integer n3,nq,i,j,k,l
      double precision dc(3,n3,n3),smat(n3,nq),dq(3,nq,nq)
      
      do i = 1,nq
         do j = 1,nq
            dq(:,i,j)=0d0
            do k = 1,n3
               do l = 1,n3
                  dq(1,i,j)=dq(1,i,j)+dc(1,k,l)*smat(k,i)*smat(l,j)
                  dq(2,i,j)=dq(2,i,j)+dc(2,k,l)*smat(k,i)*smat(l,j)
                  dq(3,i,j)=dq(3,i,j)+dc(3,k,l)*smat(k,i)*smat(l,j)
               end do
            end do
         end do
      end do
   end function Car2NM_2
   
   function Car2NM_Q(n3,nq,dc,smat)result(dq)
      integer n3,nq,i,j,ii
      double precision dc(3,3,n3),smat(n3,nq),dq(3,3,nq)
      
      do i = 1,nq
         dq(:,:,i)=0d0
         do j = 1,n3
            do ii = 1,3
               dq(1,ii,i)=dq(1,ii,i)+dc(1,ii,j)*smat(j,i)
               dq(2,ii,i)=dq(2,ii,i)+dc(2,ii,j)*smat(j,i)
               dq(3,ii,i)=dq(3,ii,i)+dc(3,ii,j)*smat(j,i)
            end do
         end do
      end do
   end function Car2NM_Q
   
   pure function Car2NM_Q_2(n3,nq,dc,smat)result(dq)
      integer,intent(in) :: n3,nq
      integer i,j,k,l,a
      double precision,intent(in) :: dc(3,3,n3,n3),smat(n3,nq)
      double precision :: dq(3,3,nq,nq)
      
      do i = 1,nq
         do j = 1,nq
            dq(:,:,i,j)=0d0
            do k = 1,n3
               do l = 1,n3
                  do a = 1,3
                     dq(1,a,i,j)=dq(1,a,i,j)+dc(1,a,k,l)*smat(k,i)*smat(l,j)
                     dq(2,a,i,j)=dq(2,a,i,j)+dc(2,a,k,l)*smat(k,i)*smat(l,j)
                     dq(3,a,i,j)=dq(3,a,i,j)+dc(3,a,k,l)*smat(k,i)*smat(l,j)
                  end do
               end do
            end do
         end do
      end do
   end function Car2NM_Q_2
   
   subroutine WriteDerVec_normalmode(unitt,namee,vecs,nq)
      integer unitt,n,i,ix,nq
      character(*) namee
      double precision vecs(3,nq),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nq
      do i = 1,nq
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,vecs(:,i)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVec_normalmode
   
   subroutine WriteDerVecQ_normalmode(unitt,namee,vecs,nq)
      integer unitt,n,i,ix,nq
      character(*) namee
      double precision vecs(6,nq),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nq
      do i = 1,nq
         q=NewQ(vecs(:,i))
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')  q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')  q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVecQ_normalmode
   
   subroutine WriteDerVecQtr_normalmode(unitt,namee,vecs,nq)
      integer unitt,n,i,ix,nq
      character(*) namee
      double precision vecs(6,nq),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nq
      do i = 1,nq
         q=TracelessQ(vecs(:,i))
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVecQtr_normalmode
   
   
   subroutine TDD2_write(unitt,fn,d2,n3)
      integer n3
      double precision d2(3,n3,n3)
      character(*) fn
      integer unitt,a,i,j
      
      open(unitt,file=fn)
      write(unitt,*)n3
      do i = 1,n3
         do j = 1,n3
            write(unitt,'(1X,I4,1X,I4,1X,3(G15.8,1X))')j,i,d2(:,j,i)
         end do
      end do
      close(unitt)
      write(output_unit,*)'WRITTEN: ',trim(adjustl(fn))
   end subroutine TDD2_write
   
   subroutine TDD3_sparse_write(unitt,fn,d3,n3)
      integer n3
      double precision d3(3,n3,n3,n3)
      character(*) fn
      integer unitt,a,i,j
      
      open(unitt,file=fn)
      write(unitt,*)n3
      do i = 1,n3
         do j = 1,n3
            write(unitt,'(1X,I4,1X,I4,1X,I4,1X,3(G15.8,1X))')j,i,i,d3(:,j,i,i)
            write(unitt,'(1X,I4,1X,I4,1X,I4,1X,3(G15.8,1X))')i,i,j,d3(:,i,i,j)
         end do
      end do
      close(unitt)
      write(output_unit,*)'WRITTEN: ',trim(adjustl(fn))
   end subroutine TDD3_sparse_write
   
   subroutine TDD2_writeq(unitt,fn,d2,n3)
      integer n3
      double precision d2(3,3,n3,n3)
      character(*) fn
      integer unitt,a,b,i,j
      
      open(unitt,file=fn)
      write(unitt,*)n3
      do i = 1,n3
         do j = 1,n3
            do a = 1,3
               write(unitt,'(1X,I1,1X,I4,1X,I4,1X,3(G15.8,1X))')a,j,i,d2(:,a,j,i)
            end do
         end do
      end do
      close(unitt)
      write(output_unit,*)'WRITTEN: ',trim(adjustl(fn))
   end subroutine TDD2_writeq
   
   subroutine TDD2_read(unitt,fn,d2,n3)
      integer n3
      double precision,allocatable :: d2(:,:,:)
      character(*) fn
      integer unitt,a,i,j
      
      open(unitt,file=fn)
      read(unitt,*)n3
      allocate(d2(3,n3,n3))
      do i = 1,n3
         do j = 1,n3
            read(unitt,'(1X,I4,1X,I4,1X,3(G15.8,1X))')a,a,d2(:,j,i)
         end do
      end do
      close(unitt)
      write(output_unit,*)'READ: ',trim(adjustl(fn))
   end subroutine TDD2_read
   
   subroutine TDD2_readq(unitt,fn,d2,n3)
      integer n3
      double precision,allocatable :: d2(:,:,:,:)
      character(*) fn
      integer unitt,a,b,i,j
      
      open(unitt,file=fn)
      read(unitt,*)n3
      allocate(d2(3,3,n3,n3))
      do i = 1,n3
         do j = 1,n3
            do a = 1,3
               read(unitt,'(1X,I1,1X,I4,1X,I4,1X,3(G15.8,1X))')b,b,b,d2(:,a,j,i)
            end do
         end do
      end do
      close(unitt)
      write(output_unit,*)'READ: ',trim(adjustl(fn))
   end subroutine TDD2_readq
   
   subroutine WriteDerVec(unitt,namee,vecs,n)
      integer unitt,n,i,ix
      character(*) namee
      double precision vecs(3,n),curvec(3)
      
      open(unitt,file=namee)
      write(unitt,*)n
      do i = 1,n
         curvec=vecs(:,i)
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,curvec
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVec
   
   subroutine WriteDerQ(unitt,namee,qtr,n)
      integer unitt,n,i,ix
      character(*) namee
      double precision q(3,3),qtr(3,3,n)
      
      open(unitt,file=namee)
      write(unitt,*)n
      do i = 1,n
         q=qtr(:,:,i)
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerQ
   
   subroutine ReadDerVec(unitt,namee,vecs,n3)
      integer unitt,n3,i,ix,nq,n
      character(*) namee
      double precision,allocatable :: vecs(:,:)
      
      open(unitt,file=namee)
      read(unitt,*)n3
      allocate(vecs(3,n3))
      do i = 1,n3
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,vecs(:,i)
      end do
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadDerVec
   
   subroutine ReadDerQ(unitt,namee,dq,n3)
      integer unitt,n,i,ix,n3
      character(*) namee
      double precision,allocatable :: dq(:,:,:)
      
      open(unitt,file=namee)
      read(unitt,*)n3
      allocate(dq(3,3,n3))
      do i = 1,n3
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,dq(1,:,i)
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,dq(2,:,i)
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,dq(3,:,i)
      end do
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadDerQ
   
   subroutine ReadDerVec_normalmode(unitt,namee,vecs,nq)
      integer unitt,n,i,ix,nq
      character(*) namee
      double precision,allocatable :: vecs(:,:)
      
      open(unitt,file=namee)
      read(unitt,*)nq
      allocate(vecs(3,nq))
      do i = 1,nq
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,vecs(:,i)
      end do
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadDerVec_normalmode
   
   subroutine ReadDerVecQ_normalmode(unitt,namee,q,nq)
      integer unitt,n,i,ix,nq
      character(*) namee
      double precision,allocatable :: q(:,:,:)
      
      open(unitt,file=namee)
      read(unitt,*)nq
      allocate(q(3,3,nq))
      do i = 1,nq
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,q(1,:,i)
         read(unitt,'(1X,6X,1X,3(E15.8E2,1X))')  q(2,:,i)
         read(unitt,'(1X,6X,1X,3(E15.8E2,1X))')  q(3,:,i)
      end do
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadDerVecQ_normalmode
   
   subroutine ReadDerVecQtr_normalmode(unitt,namee,q,nq)
      integer unitt,n,i,ix,nq
      character(*) namee
      double precision,allocatable :: q(:,:,:)
      
      open(unitt,file=namee)
      read(unitt,*)nq
      allocate(q(3,3,nq))
      do i = 1,nq
         read(unitt,'(1X,I6,1X,3(E15.8E2,1X))')n,q(1,:,i)
         read(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(2,:,i)
         read(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(3,:,i)
      end do
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadDerVecQtr_normalmode
   
   subroutine WriteVec(unitt,namee,vec)
      integer unitt
      character(*) namee
      double precision vec(3)
      
      open(unitt,file=namee)
      write(unitt,'(3(1X,F15.8))')vec
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteVec
   
   subroutine WriteTen(unitt,namee,ten)
      integer unitt
      character(*) namee
      double precision ten(3,3)
      
      open(unitt,file=namee)
      write(unitt,'(3(1X,F15.8))')ten(1,:)
      write(unitt,'(3(1X,F15.8))')ten(2,:)
      write(unitt,'(3(1X,F15.8))')ten(3,:)
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteTen
   
   subroutine ReadVec(unitt,namee,vec)
      integer unitt
      character(*) namee
      double precision vec(3)
      
      open(unitt,file=namee)
      read(unitt,'(3(1X,F15.8))')vec
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadVec
   
   subroutine ReadTen(unitt,namee,ten)
      integer unitt
      character(*) namee
      double precision ten(3,3)
      
      open(unitt,file=namee)
      read(unitt,'(3(1X,F15.8))')ten(1,:)
      read(unitt,'(3(1X,F15.8))')ten(2,:)
      read(unitt,'(3(1X,F15.8))')ten(3,:)
      close(unitt)
      write(output_unit,*)'Read: '//TR(namee)
   end subroutine ReadTen
   
   function PrimitiveQ(q_6)result(q)
      double precision q_6(6),q(3,3)
      !1  2  3  4  5  6
      !xx xy xz yy yz zz
      q(1,1)=q_6(1)
      q(2,2)=q_6(4)
      q(3,3)=q_6(6)
      
      q(1,2)=q_6(2)
      q(2,1)=q(1,2)
      
      q(1,3)=q_6(3)
      q(3,1)=q(1,3)
      
      q(2,3)=q_6(5)
      q(3,2)=q(2,3)
   end function PrimitiveQ
   
   function TracelessQ(q_6)result(q_tr)
      double precision q_6(6),q_tr(3,3)
      double precision trace
      !Q
      !1  2  3  4  5  6
      !xx xy xz yy yz zz
      !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa dont trust old Peter's code
      trace=q_6(1)+q_6(4)+q_6(6)
      q_tr(1,1)=(3d0*q_6(1)-trace)/2d0
      q_tr(2,2)=(3d0*q_6(4)-trace)/2d0
      q_tr(3,3)=(3d0*q_6(6)-trace)/2d0
      q_tr(1,2)=1.5d0*q_6(2)
      q_tr(1,3)=1.5d0*q_6(3)
      q_tr(2,1)=q_tr(1,2)
      q_tr(2,3)=1.5d0*q_6(5)
      q_tr(3,1)=q_tr(1,3)
      q_tr(3,2)=q_tr(2,3)
   end function TracelessQ
   
   function NewQ(q_6)result(q_tr)
      double precision q_6(6),q_tr(3,3)
      
      !Q
      !1  2  3  4  5  6
      !xx xy xz yy yz zz
      q_tr(1,1)=q_6(1)
      q_tr(2,2)=q_6(4)
      q_tr(3,3)=q_6(6)
      q_tr(1,2)=q_6(2)
      q_tr(1,3)=q_6(3)
      q_tr(2,1)=q_tr(1,2)
      q_tr(2,3)=q_6(5)
      q_tr(3,1)=q_tr(1,3)
      q_tr(3,2)=q_tr(2,3)
   end function NewQ
   
end module util