program sumrules
   use iso_fortran_env
   use util
   implicit none
   integer :: nq,nat,n3,i,ii,ix
   integer :: a,b,c,argc
   double precision :: u(3),v(3),m(3),q(3,3)
   double precision,allocatable :: du(:,:),dv(:,:),dm(:,:),dq(:,:,:)
   double precision :: sum_u(3),sum_q(3,3,3),sum_q_an(3,3,3)
   double precision :: sum_m(3),sum_m_an(3),w_n0
   double precision :: vec(3),step(3)
   character(80) s80
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc>=1)then
      call GET_COMMAND_ARGUMENT(1,s80)
      read(s80,*)step
   end if
   
   
   call ReadVec(77,'TM.U',u)
   call ReadVec(77,'TM.M',m)
   call ReadVec(77,'TM.V',v)
   call ReadTen(77,'TM.Qtr',q)
   
   call ReadDerVec(77,'DTM.U',du,n3)
   call ReadDerVec(77,'DTM.M',dm,n3)
   call ReadDerVec(77,'DTM.V',dv,n3)
   call ReadDerQ(77,'DTM.Qtr',dq,n3)
   nat=n3/3
   nq=n3-6
   
   open(77,file='ENERGY.VERTICAL.nm',status='old')
   read(77,*)w_n0
   close(77)
   w_n0=(1d7/w_n0)*cm_2_au
   
   write(output_unit,*)
   
   sum_u=0d0
   do i = 1,nat
      ii=(i-1)*3+1
      sum_u=sum_u+du(:,ii)
   end do
   write(output_unit,10)'Sum rule du/dx:',sum_u
   
   sum_u=0d0
   do i = 1,nat
      ii=(i-1)*3+2
      sum_u=sum_u+du(:,ii)
   end do
   write(output_unit,10)'Sum rule du/dy:',sum_u
   
   sum_u=0d0
   do i = 1,nat
      ii=(i-1)*3+3
      sum_u=sum_u+du(:,ii)
   end do
   write(output_unit,10)'Sum rule du/dz:',sum_u

10 format(A,3(1X,G12.5))
   write(output_unit,*)
   
   sum_q=0d0
   do i = 1,nat
      do ix = 1,3
         ii=(i-1)*3+ix
         sum_q(:,:,ix)=sum_q(:,:,ix)+dq(:,:,ii)
      end do
   end do
   
   sum_q_an=0d0
   do a = 1,3
      do b = 1,3
         do c = 1,3
            sum_q_an(a,b,c)=0.5*(3*(KD(a,c)*v(b)+KD(b,c)*v(a))-2*v(c)*KD(a,b))
         end do
      end do
   end do
   
   
   write(output_unit,20)'Sum(dQ/dX_i):',sum_q(:,:,1)
   write(output_unit,20)'Sum rule X:  ',sum_q_an(:,:,1)
   write(output_unit,*)

   write(output_unit,20)'Sum(dQ/dY_i):',sum_q(:,:,2)
   write(output_unit,20)'Sum rule Y:  ',sum_q_an(:,:,2)
   write(output_unit,*)
   
   write(output_unit,20)'Sum(dQ/dZ_i):',sum_q(:,:,3)
   write(output_unit,20)'Sum rule Z:  ',sum_q_an(:,:,3)
   write(output_unit,*)

20 format(A,3(3(1X,G12.5)))   

   write(output_unit,*)'Sum rules (dm/dr):'

   sum_m=0d0
   do i = 1,nat
      sum_m(1)=sum_m(1)+dm(1,(i-1)*3+2)
      sum_m(2)=sum_m(2)+dm(2,(i-1)*3+3)
      sum_m(3)=sum_m(3)+dm(3,(i-1)*3+1)
   end do

   sum_m_an=0d0
   sum_m_an(1)=-v(3)*w_n0
   sum_m_an(2)=-v(1)*w_n0
   sum_m_an(3)=-v(2)*w_n0
   sum_m_an=sum_m_an/(2d0)
   
   print *,sum_m
   print *,sum_m_an
   
   ! write(output_unit,*)
   
   ! vec=0.5d0*cross(step,-v*w_n0)
   ! print *,vec
   ! print *,sum_m-sum_m_an
   ! write(output_unit,*)
   sum_m=0d0
   do i = 1,nat
      sum_m(1)=sum_m(1)+dm(1,(i-1)*3+3)
      sum_m(2)=sum_m(2)+dm(2,(i-1)*3+1)
      sum_m(3)=sum_m(3)+dm(3,(i-1)*3+2)
   end do

   sum_m_an=0d0
   sum_m_an(1)=v(2)*w_n0
   sum_m_an(2)=v(3)*w_n0
   sum_m_an(3)=v(1)*w_n0
   sum_m_an=sum_m_an/(2d0)
   
   print *,sum_m
   print *,sum_m_an
   
   
   
   
   contains

   function cross(a,b)result(res)
      double precision a(3),b(3),res(3)
       
      res(1)=a(2)*b(3)-a(3)*b(2)
      res(2)=a(3)*b(1)-a(1)*b(3)
      res(3)=a(1)*b(2)-a(2)*b(1)
   end function cross

end program sumrules