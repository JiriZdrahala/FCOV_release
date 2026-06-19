#define TR(arg) trim(adjustl(arg))
program tdd2_num
   use iso_fortran_env
   use util
   implicit none
   integer argc,nat,desired_n_files,i,j,n3,a,b,iroot,nfiles,ix,iy
   double precision :: u0(3),v0(3),m0(3),q0(3,3)
   double precision :: u_t(3),v_t(3),m_t(3),q_t(3,3)
   double precision :: e0,e_tr0,e_t,e_tr_t
   double precision,allocatable :: u(:,:),v(:,:),m(:,:),q(:,:,:)
   double precision,allocatable :: du(:,:,:),dv(:,:,:),dm(:,:,:),dq(:,:,:,:)
   double precision,allocatable :: du0(:,:),dv0(:,:),dm0(:,:),dq0(:,:,:)
   double precision,allocatable :: du_t(:,:),dv_t(:,:),dm_t(:,:),dq_t(:,:,:)
   double precision,allocatable :: du2(:,:,:),dv2(:,:,:),dm2(:,:,:),dq2(:,:,:,:)
   double precision,allocatable :: du3(:,:,:,:),dv3(:,:,:,:),dm3(:,:,:,:),dq3(:,:,:,:,:)
   double precision,allocatable :: e_trs(:),e_gr(:)
   double precision,allocatable :: grad_gr(:),grad_ex(:),ff_ex(:,:),ff_gr(:,:)
   integer,allocatable :: z(:),z_t(:)
   double precision,allocatable :: r(:),r_t(:)
   double precision :: mult_tm(5),step
   character(500) fn,fn_cur
   character(80) s80
   character(4) plm
   logical :: cen=.false.,for=.false.,bac=.false.,doubl=.false.,der3=.false.
   double precision :: sr_u(3)
   
   
   mult_tm=[1d0,1d0,0.5d0,1d0,1d0]
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc<2)then
      write(output_unit,*)'USAGE: '
      write(output_unit,*)'tdd2_num [td-outfile-0.out] [+-] [options]'
      write(output_unit,*)'LIST.STEPS.TD is the list of displaced geometries'
      write(output_unit,*)'In case of central-difference, the files in the list are organized as'
      write(output_unit,*)'1+'
      write(output_unit,*)'2+'
      write(output_unit,*)'3+'
      write(output_unit,*)'...'
      write(output_unit,*)'1-'
      write(output_unit,*)'2-'
      write(output_unit,*)'3-'
      write(output_unit,*)'etc.'
      call exit(1)
   end if
   
   call GET_COMMAND_ARGUMENT(1,fn)
   call GET_COMMAND_ARGUMENT(2,plm)
   
   for=index(plm,'+')>0
   bac=index(plm,'-')>0
   cen=for.and.bac
   
   if(cen)then
      a=CountSubstring(plm,'+')
      b=CountSubstring(plm,'-')
      if(a==2 .and. b==2)doubl=.true.
   elseif(for)then
      a=CountSubstring(plm,'+')
      if(a==2)doubl=.true.
      stop 333
   elseif(bac)then
      b=CountSubstring(plm,'-')
      if(b==2)doubl=.true.
      stop 334
   end if
   
   do i = 3,argc
      call GET_COMMAND_ARGUMENT(i,s80)
      call To_lower(s80)
      select case(TR(s80))
         case('der3')
            der3=.true.
         case default
            write(output_unit,*)'Unknown option: ',TR(s80)
      end select
   end do
   
   call rd_td_derivatives(77,fn,z,r,nat,u0,v0,m0,q0,du0,dv0,dm0,dq0,e0,e_tr0,iroot,mult_tm)  
   n3=3*nat
   if(cen)then
      desired_n_files=3*nat*2
   else
      desired_n_files=3*nat
   end if
   if(doubl)desired_n_files=desired_n_files*2
   
   allocate(du(3,n3,desired_n_files),dv(3,n3,desired_n_files),dm(3,n3,desired_n_files),dq(3,3,n3,desired_n_files),u(3,desired_n_files),v(3,desired_n_files),m(3,desired_n_files),q(3,3,desired_n_files),e_trs(desired_n_files),e_gr(desired_n_files))
   
   step=0d0
   open(77,file='LIST.STEPS.TD')
   i=0
   do
      read(77,'(a500)',end=20)fn_cur
      i=i+1
      call rd_td_derivatives(78,fn_cur,z_t,r_t,nat,u_t,v_t,m_t,q_t,du_t,dv_t,dm_t,dq_t,e_t,e_tr_t,iroot,mult_tm)
      if(step==0d0)then
         step=sum(r_t-r)
         write(output_unit,'(A,G11.4,A)')'STEP: ',step,' (Angstrom)'
         step=step/BohrR
         write(output_unit,'(A,G11.4,A)')'      ',step,' (a.u.)'
      end if
      
      u(:,i)=u_t
      v(:,i)=v_t
      m(:,i)=m_t
      q(:,:,i)=q_t
      
      du(:,:,i)=du_t
      dv(:,:,i)=dv_t
      dm(:,:,i)=dm_t
      dq(:,:,:,i)=dq_t
      
      e_trs(i)=e_tr_t
      e_gr(i)=e_t
   end do
20 continue
   nfiles=i
   close(77)
   if(nfiles/=desired_n_files)then
      write(output_unit,*)'Wrong number of finite-difference files!'
      write(output_unit,*)'Expected: ',desired_n_files
      write(output_unit,*)'     Got: ',nfiles
      call exit(6)
   end if
   
   allocate(grad_gr(n3),grad_ex(n3))
   call readgrad(nat,grad_gr,'FILE.GR.GROUND')
   call readgrad(nat,grad_ex,'FILE.GR.EXCITED')
   ! do i = 1,n3
      ! if(cen)then
         ! grad_gr(i)=cendiff_R(e_gr(i),e_gr(i+n3),step)
         ! grad_ex(i)=cendiff_R(e_trs(i),e_trs(i+n3),step)
      ! else
         ! stop 2 !TODO
      ! end if
   ! end do
   
   
   allocate(du2(3,n3,n3),dv2(3,n3,n3),dm2(3,n3,n3),dq2(3,3,n3,n3))
   if(der3)then 
      allocate(du3(3,n3,n3,n3),dv3(3,n3,n3,n3),dm3(3,n3,n3,n3),dq3(3,3,n3,n3,n3))
      du3=0d0
      dv3=0d0
      dm3=0d0
      dq3=0d0
   end if
   do i = 1,n3
      do j = 1,n3
         if(cen.and.doubl)then
            do a = 1,3
               du2(a,j,i)=cendiff_4p_R(du(a,j,i+2*n3),du(a,j,i),du(a,j,i+n3),du(a,j,i+3*n3),step) 
               dm2(a,j,i)=cendiff_4p_R(dm(a,j,i+2*n3),dm(a,j,i),dm(a,j,i+n3),dm(a,j,i+3*n3),step) 
               dv2(a,j,i)=cendiff_4p_R(dv(a,j,i+2*n3),dv(a,j,i),dv(a,j,i+n3),dv(a,j,i+3*n3),step) 
               if(der3)stop 55
               do b = 1,3
                  dq2(b,a,j,i)=cendiff_4p_R(dq(b,a,j,i+2*n3),dq(b,a,j,i),dq(b,a,j,i+n3),dq(b,a,j,i+3*n3),step)
               end do
               
            end do
         elseif(cen)then
            do a = 1,3
               du2(a,j,i)=cendiff_R(du(a,j,i),du(a,j,i+n3),step) 
               dm2(a,j,i)=cendiff_R(dm(a,j,i),dm(a,j,i+n3),step) 
               dv2(a,j,i)=cendiff_R(dv(a,j,i),dv(a,j,i+n3),step)
               if(der3)then
                  du3(a,j,i,i)=cendiff2_R(du(a,j,i),du0(a,j),du(a,j,i+n3),step)
                  du3(a,i,i,j)=du3(a,j,i,i)
                  dm3(a,j,i,i)=cendiff2_R(dm(a,j,i),dm0(a,j),dm(a,j,i+n3),step)
                  dm3(a,i,i,j)=dm3(a,j,i,i)
                  dv3(a,j,i,i)=cendiff2_R(dv(a,j,i),dv0(a,j),dv(a,j,i+n3),step)
                  dv3(a,i,i,j)=dv3(a,j,i,i)
               end if
               do b = 1,3
                  dq2(b,a,j,i)=cendiff_R(dq(b,a,j,i),dq(b,a,j,i+n3),step)
                  if(der3)then
                     dq3(b,a,j,i,i)=cendiff2_R(dq(b,a,j,i),dq0(b,a,j),dq(b,a,j,i+n3),step)
                     dq3(b,a,i,i,j)=dq3(b,a,j,i,i)
                  end if
               end do
            end do
         else
            stop 3
         end if
      end do
   end do
   
   allocate(ff_ex(n3,n3),ff_gr(n3,n3))
   call readff(n3,ff_gr,'FILE.FC.GROUND')
   call readff(n3,ff_ex,'FILE.FC.EXCITED')
   dv2=VelToLen_der2(v,dv,dv2,e_tr0,grad_ex,grad_gr,ff_ex,ff_gr,n3)
   dq2=VelToLen_der2q(q,dq,dq2,e_tr0,grad_ex,grad_gr,ff_ex,ff_gr,n3)
   
   call TDD2_write(77,'DTM2.U',du2,n3)
   call TDD2_write(77,'DTM2.M',dm2,n3)
   call TDD2_write(77,'DTM2.V',dv2,n3)
   call TDD2_writeq(77,'DTM2.Qtr',dq2,n3)
   
   if(der3)then
      call TDD3_sparse_write(77,'DTM3.U',du3,n3)
      call TDD3_sparse_write(77,'DTM3.M',dm3,n3)
   end if
   
   sr_u=0d0
   do i = 1,nat
      ix=(i-1)*3+1 !x
      iy=(i-1)*3+1 !x
      sr_u=sr_u+dm2(:,ix,iy)
   end do
   print *,sr_u
   
   sr_u=0d0
   do i = 1,nat
      ix=(i-1)*3+1 !x
      iy=(i-1)*3+2 !x
      sr_u=sr_u+dm2(:,ix,iy)
   end do
   print *,sr_u
   
   sr_u=0d0
   do i = 1,nat
      ix=(i-1)*3+1 !x
      iy=(i-1)*3+3 !x
      sr_u=sr_u+dm2(:,ix,iy)
   end do
   print *,sr_u
   
   
   contains
   
   
   function diff_R(fx,fxh,h)result(res)
      double precision fx,fxh,res
      double precision h
      res=(fxh-fx)/h
   end function diff_R
   
   function cendiff_R(fxhp,fxhm,h)result(res)
      double precision h
      double precision fxhm,fxhp,res
      res=(fxhp-fxhm)/(2d0*h)
   end function cendiff_r
   
   function cendiff2_R(fxhp,fx,fxhm,h)result(res)
      double precision h
      double precision fxhm,fxhp,res,fx
      res=(fxhp-2*fx+fxhm)/h**2
   end function cendiff2_R
   
   function cendiff_4p_R(fx2hp,fxhp,fxhm,fx2hm,h)result(res)
      double precision h
      double precision fx2hp,fxhm,fxhp,fx2hm,res
      res=(-fx2hp+8d0*fxhp-8*fxhm+fx2hm)/(12d0*h)
   end function cendiff_4p_r
   
   function cendiff2_4p_R(fx2hp,fxhp,fx,fxhm,fx2hm,h)result(res)
      double precision h
      double precision fx2hp,fxhm,fx,fxhp,fx2hm,res
      res=(-fx2hp+16*fxhp-30*fx+16*fxhm-fx2hm)/(12*h**2)
   end function cendiff2_4p_R
   
end program tdd2_num