#define TR(arg) trim(adjustl(arg))
program elpol
   use iso_fortran_env
   use util
   use constants
   use strings
   implicit none
   
   integer i,argc,idx_eq
   integer n3,nq,iz,nat,nq_buf
   integer,allocatable :: z_at(:)
   double precision u(3),v(3),m(3),q_tr(3,3)
   double precision,allocatable :: du(:,:),dv(:,:),dm(:,:),dq_tr(:,:,:),wg_fake(:)
   double precision,allocatable :: wg(:),smat(:,:),r(:),grad_gr(:),grad_ex(:)
   double precision :: wexc(1),gamma=500,theta=250,e00,e_gr,e_ex
   logical :: finp_found=.false.,fex=.false.,grad=.false.,e00_found=.false.,vel=.false.,st=.true.,cart=.false.,vert=.true.
   character(500) s80,s80_l,s80_2
   type(Polar),allocatable :: polars(:),polars_new(:)
   
   wexc(1)=532
   argc=COMMAND_ARGUMENT_COUNT()
   do i = 1,argc
      call GET_COMMAND_ARGUMENT(i,s80)
      call To_Lower(s80_l)
      
      if(index(s80,'vel')>0)then
         vel=.true.
      elseif(index(s80,'grad')>0)then
         grad=.true.
      elseif(index(s80,'nost')>0)then
         st=.false.
      elseif(index(s80,'cart')>0)then
         cart=.true.
      elseif(index(s80,'novert')>0)then
         vert=.false.
      elseif(index(s80,'gamma=')>0)then
         read(s80(7:),*)gamma
      elseif(index(s80,'wexc=')>0)then
         read(s80(6:),*)wexc(1)
      end if
      
      ! idx_eq=index(s80_l,'=')
      ! if(idx_eq==0)idx_eq=index(s80_l,' ')-1
      
      
      ! select case(s80_l(1:idx_eq))
         ! case('wexc=')
            ! read(s80(idx_eq+1:),*)wexc(1)
         ! case('grad')
            ! grad=.true.
         ! case('vel')
            ! vel=.true.
         
      ! end select
   end do
   
   call readsi(n3,smat,wg,nq,'ground/F.INP',z_at,r,.true.,iz)
   finp_found=.true.
   wg=wg*cm_2_au
   if(.not.e00_found)then
      if(vert)then
         open(77,file='ENERGY.VERTICAL')
         read(77,*)e00
         close(77)
      else
         open(77,file='ground/ENERGY')
         read(77,*)e_gr
         close(77)
         open(77,file='excited/ENERGY')
         read(77,*)e_ex
         close(77)
         e00=e_ex-e_gr
      end if
   end if
   
   nat=size(z_at,dim=1)
   allocate(polars(nq))
   gamma=gamma*cm_2_au
   theta=theta*cm_2_au
   wexc(1)=(1d7/wexc(1))*cm_2_au
   
   call ReadVec(77,'TM.U',u)
   call ReadVec(77,'TM.M',m)
   call ReadVec(77,'TM.V',v)
   call ReadTen(77,'TM.Qtr',q_tr)
   
   call ReadDerVec_normalmode(77,'DTM.Q.U',du,nq_buf)
   call ReadDerVec_normalmode(77,'DTM.Q.M',dm,nq_buf)
   call ReadDerVec_normalmode(77,'DTM.Q.V',dv,nq_buf)
   call ReadDerVecQtr_normalmode(77,'DTM.Q.Qtr',dq_tr,nq_buf)
   
   if(grad)then
      allocate(grad_gr(nq),grad_ex(nq))
      call readgradq(nq,grad_gr,'FILE.Q.GR.ground')
      call readgradq(nq,grad_ex,'FILE.Q.GR.excited.ground')
   end if
   
   call ElPoldo(wg,nat,nq,e00,wexc(1),gamma,theta,u,v,m,q_tr,du,dv,dm,dq_tr,polars,grad,grad_gr,grad_ex,vel,st)
   if(.not.cart)then !NM
      call WritePolars(.true.,nq,e00,77,wexc(1),polars,wg,.true.,.true.)
   else !Cartesian
      n3=nq+6
      allocate(polars_new(n3))
      do i = 1,n3
         polars_new(i)=NM2Car_Pol1(n3,nq,i,polars,smat)
      end do
      allocate(wg_fake(n3))
      wg_fake=0d0
      call WritePolars(.true.,n3,e00,77,wexc(1),polars_new,wg_fake,.true.,.false.)
      deallocate(wg_fake)
   end if
   
   
   contains

   subroutine ElPoldo(wg,nat,nq,w_jn,wexc,gamma,theta,u,v,m,q,du,dv,dm,dq,polars,grad,grad_gr,grad_ex,vel,st)
      logical grad,vel,st
      integer nat,dl_idx,nq_orig,n3,iz,q_idx
      integer a,b,c,i,iexc
      integer,value :: nq
      double precision w_jn,wexc,gamma,sqrt_w,theta
      double precision u(3),v(3),m(3),q(3,3)
      double precision grad_ex(:),grad_gr(:)
      double precision :: wg(nq)
      double complex :: f,f2
      double precision :: du(3,nq),dm(3,nq),dv(3,nq),dq(3,3,nq)
      type(Polar),allocatable,target :: polars(:)
      type(Polar),pointer :: pol
      double complex :: df,df2
            
      
      
      if(vel)then
         u=v
         du=dv
      end if

      
      df=0
      df2=0
      do i = 1,nq
         sqrt_w=1d0/sqrt(2d0*wg(i))
         pol=>polars(i)
         f=((w_jn-wexc-iu*gamma))
         f2=((w_jn+(wexc-wg(i))-iu*gamma))
         if(grad)then
            df=1d0/(f**2)*(grad_ex(i)-grad_gr(i))*sqrt_w
            df2=1d0/(f2**2)*(grad_ex(i)-grad_gr(i))*sqrt_w
         end if
         f=1d0/f*sqrt_w
         f2=1d0/f2*sqrt_w
         if(.not.st)then
            f2=0d0
            df2=0d0
         end if
         do a = 1,3
            do b = 1,3
               pol%ap(a,b)=du(a,i)*u(b)*f + u(a)*du(b,i)*f + u(a)*u(b)*df + &
                                du(b,i)*u(a)*f2 + u(b)*du(a,i)*f2 + u(b)*u(a)*df2
               pol%G(a,b)=(du(a,i)*m(b)*f + u(a)*dm(b,i)*f + u(a)*m(b)*df - &
                               dm(b,i)*u(a)*f2 - m(b)*du(a,i)*f2 - m(b)*u(a)*df2)*iu
               pol%Gc(a,b)=(-dm(a,i)*u(b)*f - m(a)*du(b,i)*f - m(a)*u(b)*df + &
                                du(b,i)*m(a)*f2 + u(b)*dm(a,i)*f2 + u(b)*m(a)*df2)*iu
               do c = 1,3
                  pol%A(a,b,c)=du(a,i)*q(b,c)*f + u(a)*dq(b,c,i)*f + u(a)*q(b,c)*df + &
                                    dq(b,c,i)*u(a)*f2 + q(b,c)*du(a,i)*f2 + q(b,c)*u(a)*df2
                  pol%Ac(a,b,c)=u(a)*dq(b,c,i)*f + du(a,i)*q(b,c)*f + u(a)*q(b,c)*df + &
                                    q(b,c)*du(a,i)*f2 + dq(b,c,i)*u(a)*f2 + q(b,c)*u(a)*df2 !TODO check
               end do
            end do
         end do
      end do
      
   end subroutine ElPoldo


end program elpol