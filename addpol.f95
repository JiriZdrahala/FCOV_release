#define TR(arg) trim(adjustl(arg))
program addpol
   use iso_fortran_env
   use constants
   use util
   use strings
   implicit none
   
   integer argc,i,n3,n,polars_n,ii,j,ix,a,b,c,d,endd,st
   character(500) s500,filepol
   character(11) ft1,ft2,cur_sign
   character(80) s80,s80_l
   integer,parameter :: Polars_max=50
   type(Polar),allocatable :: PolarsRes(:),Polars(:,:),polars_cur(:)
   double precision, allocatable :: wexcs(:)
   double precision wexc,sign2,sign1
   double precision signs(Polars_max)
   double complex :: sum_ap(3,3,3),sum_G(3,3,3),sum_A(3,3,3,3)
   double complex :: sum_G_an(3,3,3),sum_A_an(3,3,3,3)
   double precision :: sum_sq
   logical do_conjg(Polars_max),mult_i(Polars_max)
   logical :: output_qttt=.false.,output_ttt=.false.
   logical :: output_polars_c=.false.,output_polars_q=.false.
   logical :: dummy=.false.
   
   
   if(dummy)then
      call printA(sum_ap)
   end if
   
   do_conjg=.false.
   mult_i=.false.
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc==0)then
      write(output_unit,*)'USAGE:'
      write(output_unit,*)''
      call exit(1)
   end if
   
   call GET_COMMAND_ARGUMENT(1,s80)
   s80_l=s80
   call to_lower(s80_l)
   if(index(s80_l,'qttt')>0)then
      output_qttt=.true.
      write(output_unit,*)'Warning: cannot turn A to primitive quadrupole form'
   elseif(index(s80_l,'ttt')>0)then
      output_ttt=.true.
      write(output_unit,*)'Warning: cannot turn A to primitive quadrupole form'
   elseif(index(s80_l,'polq')>0)then
      output_polars_q=.true.
   elseif(index(s80_l,'polc')>0)then
      output_polars_c=.true.
   elseif(index(s80_l,'zero')>0)then
      write(output_unit,*)'Setting polarizabilities to 0 in files'
      call GET_COMMAND_ARGUMENT(2,s80)
      read(s80,*)endd
      write(output_unit,*)'Starting with atom: ',endd
      endd=endd*3
      
      do i = 3,argc
         call GET_COMMAND_ARGUMENT(i,filepol)
         polars_cur=readPol(filepol,ft1,wexc,.false.)
         if(TR(ft1)/='FILE.TTT')then
            write(output_unit,*)'Zeroing only implemented for FILE.TTT'
            stop 66
         end if
         
         n3=size(polars_cur,dim=1)
         do j = endd+1,n3
            call ZeroPol(polars_cur(j))
         end do
         call writepol(polars_cur,n3,ft1,wexc)
         deallocate(polars_cur)
      end do
      
      return
   end if
   
   polars_n=0
   do i = 2,argc,2
      call GET_COMMAND_ARGUMENT(i,cur_sign)
      j=index(cur_sign,'*')
      if(j>0)then
         do_conjg(i/2)=.true.
         cur_sign(j:j)=' '
      end if
      j=index(cur_sign,'i')
      if(j>0)then
         mult_i(i/2)=.true.
         cur_sign(j:j)=' '
      end if
      
      read(cur_sign,*)signs(i/2)
      
      call GET_COMMAND_ARGUMENT(i+1,filepol)
      polars_cur=readPol(filepol,ft1,wexc,.true.)
      n=size(polars_cur,dim=1)
      if(.not.allocated(polars))then
         allocate(polars(n,polars_max),wexcs(polars_max))
         
         do ii = 1,n
            do j = 1,polars_max
               polars(ii,j)=0d0
               wexcs(j)=0d0
            end do
         end do
      end if
      wexcs(i/2)=wexc
      Polars(:,i/2)=polars_cur
      polars_n=polars_n+1
   end do
   wexc=wexcs(1)
   allocate(polarsres(n))
   do i = 1,n
      polarsres(i)=0d0
   end do
   do i = 1,polars_n
      do ii = 1,n
         if(do_conjg(i))call ConjPol(Polars(ii,i))
         if(mult_i(i))call MultIPol(Polars(ii,i))
         PolarsRes(ii)=PolarsRes(ii)+signs(i)*Polars(ii,i)
      end do
   end do
   
   if(output_ttt .or. output_polars_c)then
      sum_ap=(0d0,0d0)
      sum_G=(0d0,0d0)
      sum_a=(0d0,0d0)
      do i = 1,n
         ix=mod(i-1,3)+1
         do a = 1,3
            do b = 1,3
               sum_ap(a,b,ix)=sum_ap(a,b,ix)+PolarsRes(i)%ap(a,b)
               sum_G(a,b,ix)=sum_G(a,b,ix)+PolarsRes(i)%G(a,b)
               sum_G_an(a,b,ix)=0d0
               do d = 1,3
                  sum_G_an(a,b,ix)=sum_G_an(a,b,ix)+LC(b,d,ix)*polarsRes(i)%ap(a,d)
               end do
               sum_G_an(a,b,ix)=sum_G_an(a,b,ix)*wexc/2d0
               
               do c = 1,3
                  sum_a(a,b,c,ix)=sum_a(a,b,c,ix)+PolarsRes(i)%A(a,b,c)
                  sum_a_an(a,b,c,ix)=0.5d0*(3d0*(KD(b,ix)*polarsres(i)%ap(a,c)+KD(c,ix)*polarsres(i)%ap(a,b))-2*polarsres(i)%ap(a,ix)*KD(b,c))
               end do
            end do
         end do
      end do
      
      sum_sq=0d0
      do ix = 1,3
         do a = 1,3
            do b = 1,3
               sum_sq=sum_sq+abs(sum_ap(a,b,ix))**2
            end do 
         end do 
      end do
      
      write(output_unit,*)sum_sq
      
      sum_sq=0d0
      do ix = 1,3
         do a = 1,3
            do b = 1,3
               sum_sq=sum_sq+abs(sum_G(a,b,ix)-sum_G_an(a,b,ix))**2
            end do 
         end do 
      end do
      
      write(output_unit,*)sum_sq
      
      sum_sq=0d0
      do ix = 1,3
         do a = 1,3
            do b = 1,3
               do c = 1,3
                  sum_sq=sum_sq+abs(sum_A(a,b,c,ix)-sum_A_an(a,b,c,ix))**2
               end do
            end do 
         end do 
      end do
      
      write(output_unit,*)sum_sq
      
      
      
   end if
   
   if(output_qttt)then
      call writePol(PolarsRes,n,'FILE.Q.TTT',wexc)
   end if
   if(output_ttt)then
      call writePol(PolarsRes,n,'FILE.TTT',wexc)
   end if
   if(output_polars_q)then
      call writePol(PolarsRes,n,'FILE.POLARS.Q',wexc)
   end if
   if(output_polars_c)then
      call writePol(PolarsRes,n,'FILE.POLARS',wexc)
   end if
   
   deallocate(Polars,polarsres)
   contains
   
   subroutine ZeroPol(pol)
      integer a,b,c
      type(Polar) pol
      
      do a = 1,3
         do b = 1,3
            pol%ap(a,b)=0d0
            pol%G(a,b)=0d0
            pol%Gc(a,b)=0d0
            do c = 1,3
               pol%A(a,b,c)=0d0
               pol%Ac(a,b,c)=0d0
            end do
         end do
      end do
   end subroutine ZeroPol
   
   subroutine ConjPol(pol)
      integer a,b,c
      type(Polar) pol
      
      do a = 1,3
         do b = 1,3
            pol%ap(a,b)=conjg(pol%ap(a,b))
            pol%G(a,b)=conjg(pol%G(a,b))
            pol%Gc(a,b)=conjg(pol%Gc(a,b))
            do c = 1,3
               pol%A(a,b,c)=conjg(pol%A(a,b,c))
               pol%Ac(a,b,c)=conjg(pol%Ac(a,b,c))
            end do
         end do
      end do
   end subroutine ConjPol
   
   subroutine MultIPol(pol)
      integer a,b,c
      type(Polar) pol
      
      do a = 1,3
         do b = 1,3
            pol%ap(a,b)=(pol%ap(a,b))*iu
            pol%G(a,b)=(pol%G(a,b))*iu
            pol%Gc(a,b)=(pol%Gc(a,b))*iu
            do c = 1,3
               pol%A(a,b,c)=(pol%A(a,b,c))*iu
               pol%Ac(a,b,c)=(pol%Ac(a,b,c))*iu
            end do
         end do
      end do
   end subroutine MultIPol
   
   subroutine writePol(pol,n,ft,wexc)
      character(*) ft
      integer n,nq,n3,iz,i
      type(Polar) pol(n)
      double precision :: wg_fake(n),wexc
      double precision,allocatable :: wg(:),smat(:,:),r(:)
      integer,allocatable :: z(:)
      type(Polar),allocatable :: pol_q(:)
      
      if(ft=='FILE.TTT')then
         call WriteFilettt('FILE.TTT',pol,n,n/3,wexc)
      elseif(ft=='FILE.Q.TTT')then
         call readsi(n3,smat,wg,nq,'F.INP',z,r,.true.,iz)
         wg=wg*cm_2_au
         do i = 1,nq
            pol(i)=sqrt(2d0*wg(i))*pol(i)
         end do
         call writefileqttt('FILE.Q.TTT',pol,n,wg,wexc)
         
         deallocate(wg,smat,r,z)
      elseif(ft=='FILE.POLARS.Q')then
         call readsi(n3,smat,wg,nq,'F.INP',z,r,.true.,iz)
         wg=wg*cm_2_au
         open(77,file='FILE.POLARS.Q')
         call WritePolars(.false.,nq,0d0,77,wexc,pol,wg,.false.,.true.)
         close(77)
         
         deallocate(wg,smat,r,z)
      elseif(ft=='FILE.POLARS')then
         wg_fake=0d0
         open(77,file='FILE.POLARS')
         call WritePolars(.false.,n,0d0,77,wexc,pol,wg_fake,.false.,.false.)
         close(77)
      end if
      
   end subroutine writePol
   
   subroutine printA(A)
      double complex A(3,3,3)
      integer i
      character(1) axes(3)
      axes=['X','Y','Z']
      
      do i = 1,3
         write(output_unit,*)'i','j',axes(i)
         write(output_unit,*)realpart(A(:,1,i))
         write(output_unit,*)realpart(A(:,2,i))
         write(output_unit,*)realpart(A(:,3,i))
      end do
            
   end subroutine printA
   
   function readPol(fn,ft,wexc,q2tr)result(polres)
      character(*) fn
      type(Polar),allocatable :: polres(:),polbuf(:)
      character(11) ft
      double precision wexc,e00
      double precision,allocatable :: wg(:),smat(:,:)
      integer n
      logical isNM,q2tr
      
      ft=filetype(fn)
      if(TR(ft)=='FILE.TTT')then
         call ReadFilettt(fn,polres,n,wexc,q2tr)
         isNM=.false.
         if(output_qttt.or.output_polars_q)then
            call Transform2NM(polres,n,.true.)
         end if
      elseif(TR(ft)=='FILE.Q.TTT')then
         wexc=0d0
         call Readfileqttt(fn,polres,.false.,wg,n,wexc,q2tr)
         isNM=.true.
         if(output_ttt .or. output_polars_c)then
            call Transform2Cart(polres,n,.false.)
         end if
      elseif(TR(ft)=='FILE.POLARS')then
         polres=ReadPolars(fn,n,wexc,e00,isNM)
         if(isNM .and. (output_ttt .or. output_polars_c))then
            call Transform2Cart(polres,n,.true.)
         elseif(.not. isNM .and. (output_qttt .or. output_polars_q))then
            call Transform2NM(polres,n,.true.)
         end if
      end if
   end function readPol
   
   subroutine Transform2NM(pol,n3,sqrt_w)
      logical sqrt_w
      integer n3,nq,i,iz
      type(Polar),allocatable :: pol(:),polbuf(:)
      double precision,allocatable :: smat(:,:),wg(:),r(:)
      integer,allocatable :: z(:)
      
      call readsi(n3,smat,wg,nq,'F.INP',z,r,.true.,iz)
      allocate(polbuf(nq))
      
      do i = 1,nq
         if(sqrt_w)then
            polbuf(i)=(1d0/sqrt(2d0*wg(i)*cm_2_au))*Car2NM_Pol1(n3,nq,i,pol,smat)
         else
            polbuf(i)=Car2NM_Pol1(n3,nq,i,pol,smat)
         end if
      end do
      deallocate(pol)
      pol=polbuf
      deallocate(polbuf)
      deallocate(smat,wg,r,z)
   end subroutine Transform2NM
   
   subroutine Transform2Cart(pol,nq,sqrt_w)
      integer nq,n3,iz,i,ix
      type(Polar),allocatable :: pol(:),polbuf(:)
      double precision,allocatable :: smat(:,:),wg(:),r(:)
      double precision m_cur
      integer,allocatable :: z(:)
      logical sqrt_w
      
      call readsi(n3,smat,wg,nq,'F.INP',z,r,.true.,iz)
      
      allocate(polbuf(n3))
     ! m=0d0
      do i = 1,n3/3
         ix=(i-1)*3+1
         m_cur=(amas(z(i))*amu_2_au)
         ! m(ix,ix)=m_cur
         ! m(ix+1,ix+1)=m_cur
         ! m(ix+2,ix+2)=m_cur
         do j = 1,nq
            smat(ix,j)=smat(ix,j)*m_cur
            smat(ix+1,j)=smat(ix+1,j)*m_cur
            smat(ix+2,j)=smat(ix+2,j)*m_cur
         end do
      end do
      
      if(sqrt_w)then
         do i = 1,nq
            pol(i)=pol(i)*sqrt(2d0*wg(i)*cm_2_au)
         end do
      end if
      
      do i = 1,n3
         polbuf(i)=NM2Car_Pol1(n3,nq,i,pol,smat)
      end do
      deallocate(pol)
      pol=polbuf
      deallocate(polbuf)
      deallocate(smat,wg,r,z)
   end subroutine Transform2Cart
   
   function filetype(fn)result(res)
      character(*) fn
      character(80) s80
      character(11) res
      
      open(77,file=fn,status='old')
      read(77,'(A80)')s80
      close(77)
      if(s80(1:35)==' ROA tensors, cartesian derivatives')then
         write(res,'(A)')'FILE.TTT'
      else if(s80(1:38)==' ROA tensors, normal modes derivatives')then
         write(res,'(A)')'FILE.Q.TTT'
      else
         write(res,'(A)')'FILE.POLARS'
      end if
      
   end function filetype
   
end program addpol