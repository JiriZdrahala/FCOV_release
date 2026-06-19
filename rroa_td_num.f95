#define TR(arg) trim(adjustl(arg))

program rroa_td_num
   use iso_fortran_env
   use strings
   use util
   use constants
   use wrram
   implicit none
      
   character(500) filename,s80,s80_2
   character(80) str_v1,str_v3,istr,derm_file,dummy_str
   character(80),allocatable :: strs(:)
   character(10) wexc_str,buf
   integer :: wst_i=0,iii
   double precision wst_w
   integer nstates,nat,nstates2,ifile,argc,n3,nq,dnst,output_deriv
   integer i,ii,j,k,nexc,npx,steps,a,b
   integer,allocatable :: z(:),z_t(:),igns(:),nograds(:)
   double precision e0_gr,e,step,h_t,gamma,theta,rroa_inv(13),fwhm,e_tr,kbt,temp,w1_au,w3_au,wexc_au,dfac,max_t,sqrt_w
   double precision,allocatable :: r0(:),u0(:,:),v0(:,:),m0(:,:),q0(:,:,:),ens0(:),ens_buf(:),wg(:),wg_fake(:)
   
   double precision,allocatable :: e_gr(:)
   double precision,allocatable :: r(:,:),u(:,:,:),v(:,:,:),m(:,:,:),q(:,:,:,:),ens(:,:)
   double precision,allocatable :: r_t(:),u_t(:,:),v_t(:,:),m_t(:,:),q_t(:,:,:),ens_t(:)
   integer :: n_files,desired_n_files,dummy_int
   integer :: tabOutputs(ns0),gene_end=0
   logical :: overs=.false.,combs=.false.,fex=.false.
   logical :: for=.false.,bac=.false.,cen=.false.,vel=.false.,grad=.true.
   logical :: rroa_spr_do(ns0),spectrum_temp=.true.,use_gauss=.false.,derm=.false.
   logical :: usea=.true.,useac=.true.,useg=.true.,usegc=.true.,st=.true.,wrpol=.true.,wrpolq=.false.
   logical :: stateSwitchChec=.false.,wexc_adapt=.false.,isOrca=.false.,analy=.false.,gene=.false.,corren=.false.
   type(TD_coeff_arr),allocatable :: tdcs0(:),tdcs_t(:),tdcs(:,:)
   logical,allocatable :: ign_state(:)
   
   double precision,allocatable :: du_cur(:,:),dv_cur(:,:),dm_cur(:,:),dq_cur(:,:,:)
   double precision,allocatable :: du(:,:,:),dv(:,:,:),dm(:,:,:),dq(:,:,:,:)
   ! double precision :: du(3),dv(3),dm(3),dq(3,3)
   double precision,allocatable :: du_nm(:,:,:),dv_nm(:,:,:),dm_nm(:,:,:),dq_nm(:,:,:,:)
   double precision,allocatable :: grad_gr(:),grad_ex(:)
   double precision,allocatable :: grad_gr_nm(:),grad_ex_nm(:,:)
   double precision,allocatable :: smat(:,:),wexc(:),wexc_nm(:),wexc_h(:)
   double precision,allocatable :: rroa_spr_temp(:,:,:),rroa_spr(:,:,:)
   double complex,allocatable :: ap_sum(:,:,:),G_sum(:,:,:),A_sum(:,:,:,:)
   double complex,allocatable :: ap_sum_an(:,:,:),G_sum_an(:,:,:),A_sum_an(:,:,:,:)
   double precision :: sum_sq(3,3)
   
   double precision :: mult_tm(5)
   
   integer,allocatable :: v1(:),v3(:)
   integer(int16),allocatable :: v1_pos(:),v3_pos(:)
   integer v1_n,v3_n,ix,iexc,c,d
   
   type(Polar),allocatable :: d_pol(:,:,:),d2_pol(:,:,:),d11_pol(:,:,:)
   type(Polar),allocatable :: pols0(:,:),polars(:,:,:),polar_buf(:),polars_buf(:,:),polars_buf2(:)
   type(Polar),allocatable :: polars_nm(:),polars_cart(:)
   type(Polar) :: curpolar,d_pol_nm1,pol
   
   double precision :: w,wmin=0d0,wmax=4000d0
   integer,allocatable :: rroa_tabOutputs(:)
   
   allocate(rroa_tabOutputs(ns0))
   do i = 1,ns0
      rroa_tabOutputs(i)=99+i!starts at 100, 1st wexc ends at 118
   end do
   rroa_spr_do=.false.
   rroa_spr_do(12)=.true. !SCP180
   mult_tm=[1d0,1d0,0.5d0,1d0,1d0]
   gamma=1000
   theta=0
   dfac=1
   max_t=100 !I dunno, this is just an arbitrary value
   steps=2**8
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc<2)then
      write(output_unit,*)'USAGE: '
      write(output_unit,*)'rroa_td_num [td-outfile-0.out] [+-] [532,542,552...] [options]'
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
   
   do i = 4,argc
      call GET_COMMAND_ARGUMENT(i,s80)
      if(index(s80,'dnst=')>0)then
         read(s80(6:),*)dnst
      end if
   end do
   
   call GET_COMMAND_ARGUMENT(1,s80)
   isOrca=is_orca_output(s80)
   if(isOrca)then
      call rd_casscf_orca(77,s80,dnst,nat,ens0,r0,z,u0,m0)
      nstates=dnst
      allocate(q0(3,3,nstates),v0(3,nstates))
      q0=0d0
      v0=0d0
   else
      call rd_td_new(77,s80,nstates,z,r0,nat,u0,v0,m0,q0,e0_gr,ens0,mult_tm,dnst)
   end if
   
   
   output_deriv=0
   do i = 4,argc
      call GET_COMMAND_ARGUMENT(i,s80)
      if(index(s80,'vel')>0)then
         vel=.true.
      elseif(index(s80,'noaten')>0)then
         usea=.false.
         useac=.false.
      elseif(index(s80,'nogten')>0)then
         useg=.false.
         usegc=.false.
      elseif(index(s80,'overtones')>0)then
         overs=.true.
      elseif(index(s80,'combinations')>0)then
         combs=.true.
      elseif(index(s80,'dnst')>0)then
         continue
      elseif(index(s80,'analyitical')>0)then
         analy=.true.
      elseif(index(s80,'nograd')>0)then
         if(index(s80,'=')>0)then
            strs=splitString(s80(8:),80,',')
            allocate(nograds(size(strs,dim=1)))
            do ii = 1,size(strs,dim=1)
               read(strs(ii),*)nograds(ii)
            end do
            grad=.false.
         else
            grad=.false.
            allocate(nograds(nstates))
            do ii = 1,nstates
               nograds(ii)=ii
            end do
         end if
      elseif(index(s80,'nost')>0)then
         st=.false.
      elseif(index(s80,'statesw')>0)then
         stateSwitchChec=.true.
      elseif(index(s80,'wrpolq')>0)then
         wrpolq=.true.
      elseif(index(s80,'encorr')>0)then
         corren=.true.
      elseif(index(s80,'gene=')>0)then
         gene=.true.
         read(s80(6:),*)gene_end
      elseif(index(s80,'adapt')>0)then
         wexc_adapt=.true.
      elseif(index(s80,'deriv=')>0)then
         read(s80(7:),*)output_deriv
      elseif(index(s80,'theta=')>0)then
         read(s80(7:),*)theta
      elseif(index(s80,'gamma=')>0)then
         read(s80(7:),*)gamma
      elseif(index(s80,'steps=')>0)then
         read(s80(7:),*)steps
      elseif(index(s80,'max_t=')>0)then
         read(s80(7:),*)max_t
      elseif(index(s80,'fwhm=')>0)then
         read(s80(6:),*)fwhm
      elseif(index(s80,'wst=')>0)then
         read(s80(5:),*)buf
         ii=index(s80,':')
         read(s80(5:ii-1),*)wst_i
         read(s80(ii+1:),*)wst_w
      elseif(index(s80,'igno=')>0)then
         strs=splitString(s80(6:),80,',')
         allocate(igns(size(strs,dim=1)))
         do ii = 1,size(strs,dim=1)
            read(strs(ii),*)igns(ii)
         end do
      elseif(index(s80,'dfac=')>0)then
         read(s80(6:),*)dfac
      else
         write(output_unit,*)'Unknown option ',TR(s80)
         call exit(1)
      end if
   end do
   
   if(wst_i/=0)then
      ens0(wst_i)=1d7/(1d7/(ens0(wst_i)*au_2_cm)+wst_w)*cm_2_au
   end if
   n3=3*nat
   call GET_COMMAND_ARGUMENT(2,s80)
   ens0=ens0
   fwhm=5
   npx=4001
   write(output_unit,*)'GAMMA   = ',gamma
   !write(output_unit,*)'THETA   = ',theta
   write(output_unit,*)'NSTATES = ',nstates
   if(allocated(igns))then
      write(output_unit,*)'IGNO:'
      do i = 1,size(igns,dim=1)
         write(output_unit,'(1X,I4)',advance='no')igns(i)
      end do
      write(output_unit,*)
   end if
   
   ! if(allocated(nograds))then
      ! write(output_unit,*)'nograds:'
      ! do i = 1,size(nograds,dim=1)
         ! write(output_unit,'(1X,I4)',advance='no')nograds(i)
      ! end do
      ! write(output_unit,*)
   ! end if
   
   
   if(allocated(igns))then
      allocate(ign_state(nstates))
      ign_state=.false.
      do i = 1,size(igns,dim=1)
         ign_state(igns(i))=.true.
      end do
   end if
   
   
   allocate(rroa_spr_temp(npx,2,ns0),rroa_spr(npx,2,ns0))
   
   for=index(s80,'+')>0
   bac=index(s80,'-')>0
   cen=for.and.bac
   
   call GET_COMMAND_ARGUMENT(3,s80)
   wexc=ReadWexc(s80)
   wexc_nm=wexc
   wexc=(1d7/wexc)*cm_2_au
   nexc=size(wexc,dim=1)
   do i = 1,ns0
      tabOutputs(i)=99+i
   end do
   
   if(cen)then
      if(combs)then
         desired_n_files=2*n3+n3*(n3-1)
      else
         desired_n_files=2*n3
      end if
   else
      if(combs)then
         desired_n_files=n3+n3*(n3-1)/2
      else
         desired_n_files=n3
      end if
   end if
   
   call readsi(n3,smat,wg,nq,'F.INP',z_t,r_t,.true.,i)
   deallocate(z_t,r_t)
   
   allocate(grad_gr(n3))
   open(77,file='FILE.GR',status='old')
   do i = 1,nat
      ii=(i-1)*3+1
      read(77,*)grad_gr(ii:ii+2)
   end do
   close(77)
   
   allocate(e_gr(desired_n_files),r(n3,desired_n_files),u(3,nstates,desired_n_files), &
      v(3,nstates,desired_n_files),m(3,nstates,desired_n_files),q(3,3,nstates,desired_n_files), &
      ens(nstates,desired_n_files))
   n_files=0
   step=0
   open(77,file='LIST.STEPS.TD',status='old')
! 40 read(77,'(a500)',end=50)filename
   ! n_files=n_files+1
   ! goto 40
! 50 rewind(77)
   ifile=1
   write(output_unit,*)'Reading files'

   if(isOrca)then
      allocate(q_t(3,3,dnst),v_t(3,dnst))
      q_t=0d0
      v_t=0d0
   end if

60 read(77,'(a500)',end=70)filename
   !ens_t=ens0
   !tdcs_t=tdcs0
   if(isOrca)then
      call rd_casscf_orca(78,filename,dnst,nat,ens_t,r_t,z_t,u_t,m_t)
      nstates2=dnst
   else
      call rd_td_new(78,filename,nstates2,z_t,r_t,nat,u_t,v_t,m_t,q_t,e,ens_t,mult_tm,dnst)
   end if
   if(wst_i/=0)then
      ens_t(wst_i)=1d7/(1d7/(ens_t(wst_i)*au_2_cm)+wst_w)*cm_2_au
   end if
   if(nstates2/=nstates)then
      write(output_unit,*)'ERROR: Different TD nstates'
      call exit(2)
   end if
   h_t=(sum(r_t-r0))
   if(step==0)then
      step=h_t
      write(output_unit,*)'STEP: ',step
   else 
      if(abs(step/h_t)-1>0.0001d0)then
         write(output_unit,*)'WARNING: Different step in file ',ifile
         stop 666
      end if
   end if
   e_gr(ifile)=e
   r(:,ifile)=r_t
   u(:,:,ifile)=u_t
   v(:,:,ifile)=v_t
   m(:,:,ifile)=m_t
   q(:,:,:,ifile)=q_t
   !tdcs(:,ifile)=tdcs_t
   ens(:,ifile)=ens_t
   if(.not.isOrca)deallocate(v_t,q_t)
   deallocate(r_t,u_t,m_t,ens_t,z_t)
   ifile=ifile+1
   if(desired_n_files==ifile-1)goto 70
   goto 60
70 close(77)  
   step=step/BohrR
   ifile=ifile-1
   if(desired_n_files/=ifile)then
      write(output_unit,*)'ERROR: File count in LIST.STEPS.TD wrong'
      write(output_unit,*)'Expected: ',desired_n_files,'Got: ',ifile
      call exit(3)
   end if
   
   if(gene)then
      open(77,file='EN.CORR')
      write(77,*)gene_end,ifile,' au'
      
      !equilibrium
      do i = 1,gene_end
         write(77,*)0,i,ens0(i)
      end do
      
      !finite differences
      do j = 1,ifile
         do i = 1,gene_end
            write(77,*)j,i,ens(i,j)
         end do
      end do
      
      close(77)
      write(output_unit,*)'Written EN.CORR'
      
      
      !dipole moments
      open(77,file='U.CORR')
      open(78,file='M.CORR')
      write(77,*)gene_end,ifile,' au'
      write(78,*)gene_end,ifile,' au'
      
      do i = 1,gene_end
         write(77,*)0,i,u0(:,i)
         write(78,*)0,i,m0(:,i)
      end do
      
      do j = 1,ifile
         do i = 1,gene_end
            write(77,*)j,i,u(:,i,j)
            write(78,*)j,i,m(:,i,j)
         end do
      end do
      
      close(77)
      write(output_unit,*)'Written U.CORR'
      close(78)
      write(output_unit,*)'Written M.CORR'
      call exit(0)
   end if
   
   if(corren)then
      open(77,file='EN.CORR',status='old')
      read(77,*)gene_end,dummy_int,dummy_str
      
      !equilibrium
      do i = 1,gene_end
         read(77,*)dummy_int,dummy_int,ens0(i)
      end do
      
      !finite differences
      do j = 1,ifile
         do i = 1,gene_end
            read(77,*)dummy_int,dummy_int,ens(i,j)
         end do
      end do
      
      close(77)
      write(output_unit,*)'Read EN.CORR ',gene_end
      
      inquire(file='U.CORR',exist=fex)
      if(fex)then
         open(77,file='U.CORR')
         read(77,*)gene_end,dummy_int,dummy_str
         
         do i = 1,gene_end
            read(77,*)dummy_int,dummy_int,u0(:,i)
         end do
         
         do j = 1,ifile
            do i = 1,gene_end
               read(77,*)dummy_int,dummy_int,u(:,i,j)
            end do
         end do
         
         close(77)
         write(output_unit,*)'Read U.CORR'
      end if
      
      inquire(file='M.CORR',exist=fex)
      if(fex)then
         open(77,file='M.CORR')
         read(77,*)gene_end,dummy_int,dummy_str
         
         do i = 1,gene_end
            read(77,*)dummy_int,dummy_int,m0(:,i)
         end do
         
         do j = 1,ifile
            do i = 1,gene_end
               read(77,*)dummy_int,dummy_int,m(:,i,j)
            end do
         end do
         
         close(77)
         write(output_unit,*)'Read M.CORR'
      end if
      
   end if
   
   write(output_unit,*)'First 10 excited states:'
   do i = 1,min(10,nstates)
      write(output_unit,'(1X,F9.2)',advance='no')1d7/((ens0(i))*au_2_cm)
   end do
   write(output_unit,*)
   
   do i=1,nstates
      e_tr=ens0(i)
      ens0(i)=e_tr
      v0(:,i)=VelToLen_vec_cart(v0(:,i),e_tr)
      q0(:,:,i)=VelToLen_vec_cartQ(q0(:,:,i),e_tr)
      do ii = 1,ifile
         e_tr=ens(i,ii)
         v(:,i,ii)=VelToLen_vec_cart(v(:,i,ii),e_tr)
         q(:,:,i,ii)=VelToLen_vec_cartQ(q(:,:,i,ii),e_tr)
      end do
   end do
   
   do ii = 1,n3
      if(cen)then
         grad_gr(ii)=cendiff_R(e_gr(ii),e_gr(ii+n3),step)
      else
         grad_gr(ii)=diff_R(e0_gr,e_gr(ii),step)
      end if
   end do
   
   if(analy)allocate(du(3,n3,nstates),dm(3,n3,nstates),dv(3,n3,nstates),dq(3,3,n3,nstates))
   if(output_deriv /= 0 .or. analy)then
      allocate(du_cur(3,n3),dv_cur(3,n3),dm_cur(3,n3),dq_cur(3,3,n3),grad_ex(n3))
      
      if(analy)i=0
      
111   continue
      if(analy)then
         i=i+1
      else
         i=output_deriv
      end if
      
      do j=1,n3
         if(cen)then
            grad_ex(j)=cendiff_r(ens(i,j),ens(i,j+n3),step)
         else
            grad_ex(j)=cendiff_r(ens0(i),ens(i,j),step)
         end if
         do a = 1,3
            if(cen)then
               du_cur(a,j)=cendiff_R(u(a,i,j),u(a,i,j+n3),step)
               dv_cur(a,j)=cendiff_R(v(a,i,j),v(a,i,j+n3),step)
               dm_cur(a,j)=cendiff_R(m(a,i,j),m(a,i,j+n3),step)
               do b = 1,3
                  dq_cur(a,b,j)=cendiff_R(q(a,b,i,j),q(a,b,i,j+n3),step)
               end do
            else
               du_cur(a,j)=diff_R(u0(a,i),u(a,i,j),step)
               dv_cur(a,j)=diff_R(v0(a,i),v(a,i,j),step)
               dm_cur(a,j)=diff_R(m0(a,i),m(a,i,j),step)
               do b = 1,3
                  dq_cur(a,b,j)=diff_R(q0(a,b,i),q(a,b,i,j),step)
               end do
            end if
         end do
      end do
      if(analy)then
         du(:,:,i)=du_cur
         dm(:,:,i)=dm_cur
         dv(:,:,i)=dv_cur
         dq(:,:,:,i)=dq_cur
      end if
      
      
      
      if(analy .and. i<nstates)goto 111
      if(output_deriv/=0)then
         call WriteVec(77,'TM.U',u0(:,i))
         call WriteVec(77,'TM.M',m0(:,i))
         call WriteVec(77,'TM.V',v0(:,i))
         call WriteTen(77,'TM.Qtr',q0(:,:,i))
         
         call WriteDerVec(77,'DTM.U',du_cur,n3)
         call WriteDerVec(77,'DTM.M',dm_cur,n3)
         call WriteDerVec(77,'DTM.V',dv_cur,n3)
         call WriteDerQ(77,'DTM.Qtr',dq_cur,n3)
         
         call WriteGrad(77,'FILE.GR.excited',grad_ex,n3)
      end if
      deallocate(du_cur,dv_cur,dm_cur,dq_cur,grad_ex)
      if(output_deriv/=0) call exit(0)
   end if
   !if(derm)deallocate(du,dv,dm,dq)
   
   if(vel)then
      u=v
   end if
   
   if(.not.grad)then
      allocate(polars(nexc,1,nq))
   else
      allocate(polars(nexc,ifile,nq))
   end if
   gamma=gamma*cm_2_au
   allocate(d_pol(nexc,n3,nq))
   if(overs)allocate(d2_pol(nexc,n3,nq))
   if(combs)allocate(d11_pol(nexc,n3*(n3-1)/2,nq))
   
   write(output_unit,*)'Doing geometry polarizabilities'
   !if(.not.cen .or. overs .or. combs)then
      allocate(pols0(nexc,nq))
   !end if
   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP PRIVATE(polars_buf,polars_buf2) &
   !$OMP PRIVATE(i,j,k,sqrt_w) &
   !$OMP SHARED(n3,step,nq,ifile,u,m,q,ens,nstates,wexc,nexc,gamma,wg,polars,d_pol,wexc_adapt) &
   !$OMP SHARED(st,pols0,ens0,u0,m0,v0,q0,ign_state,grad,smat,cen,overs,d2_pol,combs)
   allocate(polars_buf(nexc,ifile),polars_buf2(nexc))
   
   
   !$OMP MASTER
   write(output_unit,*)'Fundamentals...'
   !$OMP END MASTER
   !$OMP DO
   do i = 1,nq
      if(.not.grad)then
         call TDPolar_nograd(i,ifile,nq,u,m,q,u0,m0,q0,ens0,nstates,ign_state,wexc,nexc,gamma,wg(i)*cm_2_au,polars_buf2,st,wexc_adapt,smat)
         polars(:,1,i)=polars_buf2
      else
         !if(.not.cen .or. overs)then
         call TDPolar(u0,m0,q0,ens0,nstates,ign_state,wexc,nexc,gamma,wg(i)*cm_2_au,pols0(:,i),st,wexc_adapt)
         !end if
         
         do k = 1,ifile
            call TDPolar(u(:,:,k),m(:,:,k),q(:,:,:,k),ens(:,k),nstates,ign_state,wexc,nexc,gamma,wg(i)*cm_2_au,polars_buf2,st,wexc_adapt)
            polars_buf(:,k)=polars_buf2
         end do
         polars(:,:,i)=polars_buf
      end if
   end do
   !$OMP END DO
   deallocate(polars_buf,polars_buf2)
   !$OMP MASTER
   if(grad)write(output_unit,*)'Doing polarizability derivatives'
   !$OMP END MASTER
   if(grad)then
      if(cen)then
         !$OMP DO
         do i = 1,nq
            d_pol(:,:,i)=CenDiff_pols(polars(:,:,i),n3,nexc,step)
         end do
         !$OMP END DO
         
         if(overs)then
            !$OMP DO
            do i = 1,nq
               d2_pol(:,:,i)=CenDiff_pols(polars(:,:,i),n3,nexc,step)
            end do
            !$OMP END DO
         end if
      else 
         !$OMP DO
         do i = 1,nq
            d_pol(:,:,i)=Diff_pols(pols0,polars(:,:,i),n3,nexc,step)
         end do
         !$OMP END DO
      end if
   end if
   
   if(combs)then
      
   end if
   
   !$OMP END PARALLEL
   
   do ii = 1,nexc
      
      write(wexc_str,'(F6.1,"nm")')dble(NINT(1d7/(wexc(ii)*au_2_cm)*10))/10d0
      wexc_str=TR(wexc_str)
         
      wexc_au=wexc(ii)
      rroa_spr=0
      do i = 1,ns0
         if(rroa_spr_do(i))then
            open(rroa_tabOutputs(i),file=TR(rroa_exp(i))//'_'//TR(wexc_str)//'.TAB')
         end if
      end do
      allocate(polars_nm(nq))
      if(grad)then
         do i =1,nq
            polar_buf=d_pol(ii,:,i)
            d_pol_nm1=Car2NM_Pol1(n3,nq,i,polar_buf,smat)
            sqrt_w=1d0/sqrt(2d0*wg(i)*cm_2_au)
            polars_nm(i)=d_pol_nm1*sqrt_w
         end do
         allocate(polars_cart(n3))
         do i = 1,n3
            polars_cart(i)=d_pol(ii,i,max(1,nq/2))
         end do
         
         if(wrpol)then
            allocate(wg_fake(n3))
            wg_fake=0d0
            call WritePolars(.true.,n3,0d0,55,wexc_au,polars_cart,wg_fake,.true.,.false.)
            write(output_unit,'(A,A)')'Written Cartesian polarizability derivative for ',wexc_str
            deallocate(wg_fake)
         end if
         if(wrpolq)then
            call WritePolars(.true.,nq,0d0,55,wexc_au,polars_nm,wg*cm_2_au,.true.,.true.)
            write(output_unit,'(A,A)')'Written transition polarizability for ',wexc_str
         end if
         deallocate(polars_cart)
      else
         polars_nm=polars(ii,1,:)
         if(wrpol.or.wrpolq)then
            call WritePolars(.true.,nq,0d0,55,wexc_au,polars_nm,wg*cm_2_au,.true.,.true.)
            write(output_unit,'(A,A)')'Written transition polarizability for ',wexc_str
         end if
      end if
      
      do i=1,nq
         d_pol_nm1=polars_nm(i)
         v1=[0]
         v1_pos=[0]
         v1_n=0
         v3=[1]
         v3_pos=[i]
         v3_n=1
         str_v1=FC2str_new(v1,v1_pos,v1_n,.true.,.true.)
         str_v3=FC2str_new(v3,v3_pos,v3_n,.true.,.true.)
         
         w1_au=0d0
         w3_au=wg(i)*cm_2_au
         temp=300
         spectrum_temp=.true.
         use_gauss=.false.
         !curPolar=polars(i)
         kbt=1d0
         
         rroa_spr_temp=0
         call wrram3(wexc_au,temp,spectrum_temp,w1_au,w3_au,TINY(1d0),i,wexc_nm(ii),0d0,4000d0,npx,use_gauss,FWHM,.true.,rroa_inv,&
                      realpart(d_pol_nm1%Ap),imagpart(d_pol_nm1%Ap), &
                      realpart(d_pol_nm1%G),imagpart(d_pol_nm1%G), &
                      realpart(d_pol_nm1%Gc),imagpart(d_pol_nm1%Gc), &
                      realpart(d_pol_nm1%A),imagpart(d_pol_nm1%A), &
                      realpart(d_pol_nm1%Ac),imagpart(d_pol_nm1%Ac), &
                      rroa_spr_temp,kBT,useA.or.useAc,useG.or.useGc,rroa_spr_do,tabOutputs,str_v1,str_v3, &
                      realpart(d_pol_nm1%Ap),imagpart(d_pol_nm1%Ap))
                      
         rroa_spr(:,1,:)=rroa_spr(:,1,:)+rroa_spr_temp(:,1,:)
         rroa_spr(:,2,:)=rroa_spr(:,2,:)+rroa_spr_temp(:,2,:)
      end do
      deallocate(polars_nm)
      do i = 1,ns0
         if(rroa_spr_do(i))then
            close(rroa_tabOutputs(i))
         end if
      end do
      do i = 1,ns0
         if(.not.rroa_spr_do(i))cycle
         write(s80,3000)TR(rroa_exp(i)),TR(wexc_str)
         write(s80_2,3001)TR(rroa_exp(i)),TR(wexc_str)
   3000  format('RAM',A,'_',A,'.PRN')
   3001  format('ROA',A,'_',A,'.PRN')

         open(77,file=s80)
         open(78,file=s80_2)
         do k = 1,npx
            w=dble(wmax-wmin)/npx*(k-1)+wmin
            write(77,'(f14.4,g25.12)')w,rroa_spr(k,1,i)
            write(78,'(f14.4,g25.12)')w,rroa_spr(k,2,i)
         end do
         close(77)
         write(output_unit,*)'Written file: '//TR(s80)
         close(78)
         write(output_unit,*)'Written file: '//TR(s80_2)
      end do
      
      do i = 1,ns0
         if(rroa_spr_do(i))then
            close(rroa_tabOutputs(i))
         end if
      end do
   end do
   deallocate(polars)
   contains
   
   subroutine NoGradCorrection(n3,nq,polars,wexc,nexc,w_jn,wg,gamma,grad_gr,grad_ex,u,v,m,q,st)
      integer n3,nq,nexc
      integer i,j,ii,jj,k,a,b,c,iexc
      double precision wexc(nexc),grad_gr(n3),grad_ex(n3),w_jn,gamma
      double precision u(3),v(3),m(3),q(3,3),wg(nq)
      double complex f2,g2
      type(Polar) polars(nexc,n3,nq)
      logical st
      
      g2=(0d0,0d0)
      do j = 1,nq
         do i = 1,n3
            do iexc=1,nexc
               f2=1d0/(w_jn-wexc(iexc)-iu*gamma)**2
               if(st)g2=1d0/(w_jn+(wexc(iexc)-wg(j))+iu*gamma)**2
               do a = 1,3
                  do b = 1,3
                     polars(iexc,i,j)%ap(a,b)=polars(iexc,i,j)%ap(a,b)+conjg(u(a)*u(b)*(f2+g2)*(grad_ex(i)-grad_gr(i)))
                     polars(iexc,i,j)%G(a,b)=polars(iexc,i,j)%G(a,b)+conjg((-u(a)*m(b)*f2+m(b)*u(a)*g2)*(grad_ex(i)-grad_gr(i)))
                     polars(iexc,i,j)%Gc(a,b)=polars(iexc,i,j)%Gc(a,b)+conjg((m(a)*u(b)*f2-u(b)*m(a)*g2)*(grad_ex(i)-grad_gr(i)))
                     do c = 1,3
                        polars(iexc,i,j)%A(a,b,c)=polars(iexc,i,j)%A(a,b,c)+conjg(u(a)*q(b,c)*(f2+g2)*(grad_ex(i)-grad_gr(i)))
                        polars(iexc,i,j)%Ac(a,b,c)=polars(iexc,i,j)%Ac(a,b,c)+conjg(u(a)*q(b,c)*(f2+g2)*(grad_ex(i)-grad_gr(i)))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine NoGradCorrection
   
   
   
   function ReadWexc(s80)result(wexcs)
      character(*) s80
      double precision,allocatable :: wexcs(:)
      integer i,j,co,start,endd,step,lenn
      i=index(s80,':')
      lenn=len(s80)
      if(i>0)then
         j=index(s80,':',back=.true.)
         read(s80(1:i-1),*)start
         read(s80(i+1:j-1),*)step
         read(s80(j+1:),*)endd
         co=(endd-start)/step+1
         allocate(wexcs(co))
         do i = 1,co
            wexcs(i)=start+(i-1)*step
         end do
      else
         co=1
         i=index(s80,',')
         do while(i>0)
            s80(i:i)=' '
            co=co+1
            i=index(s80,',')
         end do
         
         allocate(wexcs(co))
         read(s80,*)wexcs
      end if
   end function ReadWexc
   
   function Voigt(wjn,wexc,gamma,theta,t_max,steps,st,f1_arr,wg)result(res)
      double precision wjn,wexc,gamma,theta,t_max,dt,t_cur,wg
      double complex res,f1_arr(steps)
      integer i,steps
      logical st
      
      dt=t_max/dble(steps)
      if(st)then
         res=(fun_st_faster(wjn,wexc,gamma,theta,0d0,f1_arr(1),wg)+fun_st_faster(wjn,wexc,gamma,theta,t_max,f1_arr(steps),wg))/2d0
      else
         f1_arr(1)=fun(wjn,wexc,gamma,theta,0d0)
         f1_arr(steps)=fun(wjn,wexc,gamma,theta,t_max)
         res=(f1_arr(1)+f1_arr(steps))/2d0
      end if
      do i = 2,steps-1
         t_cur=(i-1)*dt
         if(st)then
            res=res+fun_st_faster(wjn,wexc,gamma,theta,t_cur,f1_arr(i),wg)
         else
            f1_arr(i)=fun(wjn,wexc,gamma,theta,t_cur)
            res=res+f1_arr(i)
         end if
      end do
      
      res=res*dt
   end function Voigt
   
   function fun(wjn,wexc,gamma,theta,t)
      double precision t,wjn,wexc,gamma,theta
      double complex fun
      fun=exp(-iu*t*(wjn-wexc-iu*gamma-iu*t*theta**2/2d0))
   end function fun
   
   function fun_st(wjn,wexc,gamma,theta,t)
      double precision t,wjn,wexc,gamma,theta
      double complex fun_st
      fun_st=exp(-iu*t*(-wjn-wexc-iu*gamma-iu*t*theta**2/2d0))
   end function fun_st
      
   function fun_st_faster(wjn,wexc,gamma,theta,t,f1,wg)
      double precision t,wjn,wexc,gamma,theta,wg
      double complex fun_st_faster,f1
      fun_st_faster=f1*exp(-iu*t*(2*wjn-wg))
   end function fun_st_faster
      
   
   subroutine TDPolar(u,m,q,ens,nst,ign_state,wexc,nexc,gamma,wg,res,st,wexc_adapt)
      integer nst,i,nexc,a,b,c,iexc
      double precision u(3,nst),m(3,nst),q(3,3,nst),ens(nst),wexc(nexc),wexc_cur,gamma,wg
      double precision wr
      double complex f,f2
      type(Polar) res(nexc)
      logical st,ffr,wexc_adapt
      logical,allocatable :: ign_state(:)
      
      do iexc=1,nexc
         res(iexc)=0d0
      end do
      f2=0d0
      ffr=.false.
      do i = 1,nst
         if(allocated(ign_state))then
            if(ign_state(i))cycle
         end if
         do iexc = 1,nexc
            wexc_cur=wexc(iexc)
            if(wexc_adapt)wexc_cur=wexc_cur-wg/2d0
            if(gamma==0d0)then !FFR
               f=2*ens(i)/(ens(i)**2-wexc_cur**2)
               f2=-2*wexc_cur/(ens(i)**2-wexc_cur**2)
               ffr=.true.
            else !in resonance, Lorentzian, currently accepted approach (2026)
               f=1d0/(ens(i)-wexc_cur-iu*gamma)
               wr=(wexc_cur-wg)
               if(st)f2=1d0/(ens(i)+(wexc_cur-wg)+iu*gamma)
            end if
            if(ffr)then
               do a = 1,3
                  do b = 1,3
                     res(iexc)%ap(a,b)=res(iexc)%ap(a,b) + u(a,i)*u(b,i)*f
                     res(iexc)%G(a,b)=res(iexc)%G(a,b) + u(a,i)*m(b,i)*iu*f2
                     res(iexc)%Gc(a,b)=res(iexc)%Gc(a,b) - m(a,i)*u(b,i)*iu*f2
                     do c = 1,3
                        res(iexc)%A(a,b,c)=res(iexc)%A(a,b,c)+u(a,i)*q(b,c,i)*f
                        res(iexc)%Ac(a,b,c)=res(iexc)%Ac(a,b,c)+u(a,i)*q(b,c,i)*f
                     end do
                  end do
               end do
            else
               do a = 1,3
                  do b = 1,3
                     res(iexc)%ap(a,b)=res(iexc)%ap(a,b) + conjg(u(a,i)*u(b,i)*f + u(b,i)*u(a,i)*f2)
                     ! res(iexc)%G(a,b)=res(iexc)%G(a,b) + u(a,i)*m(b,i)*iu*f - m(b,i)*u(a,i)*iu*f2
                     ! res(iexc)%Gc(a,b)=res(iexc)%Gc(a,b) - m(a,i)*u(b,i)*iu*f + u(b,i)*m(a,i)*iu*f2
                     res(iexc)%G(a,b)=res(iexc)%G(a,b) + conjg(u(a,i)*m(b,i)*iu*f - m(b,i)*u(a,i)*iu*f2)
                     res(iexc)%Gc(a,b)=res(iexc)%Gc(a,b) - conjg(m(a,i)*u(b,i)*iu*f + u(b,i)*m(a,i)*iu*f2)
                     do c = 1,3
                        res(iexc)%A(a,b,c)=res(iexc)%A(a,b,c)+conjg(u(a,i)*q(b,c,i)*f + q(b,c,i)*u(a,i)*f2)
                        res(iexc)%Ac(a,b,c)=res(iexc)%Ac(a,b,c)+conjg(u(a,i)*q(b,c,i)*f + q(b,c,i)*u(a,i)*f2)
                     end do
                  end do
               end do
            end if
         end do
      end do
   end subroutine TDPolar
   
   subroutine TDPolar_nograd(iq,nfiles,nq,u,m,q,u0,m0,q0,ens,nst,ign_state,wexc,nexc,gamma,wg,res,st,wexc_adapt,smat)
      integer nst,i,nexc,a,b,c,iexc,steps,nfiles,iq,nq
      double precision u(3,nst,nfiles),m(3,nst,nfiles),q(3,3,nst,nfiles),ens(nst),wexc(nexc),wexc_cur,gamma,wg
      double precision u0(3,nst),m0(3,nst),q0(3,3,nst)
      double precision wr,smat(n3,nq),sqrt_w
      double precision duu_1,duu_2,dum_1,dum_2,dmu_1,dmu_2,duq_1,duq_2,dqu_1,dqu_2
      
      double complex f,f2
      type(Polar) res(nexc)
      logical st,ffr,wexc_adapt
      logical,allocatable :: ign_state(:)
      
      do iexc=1,nexc
         res(iexc)=0d0
      end do
      f2=0d0
      ffr=.false.
      sqrt_w=1d0/sqrt(2d0*wg)
      do i = 1,nst
         if(allocated(ign_state))then
            if(ign_state(i))cycle
         end if
         do iexc = 1,nexc
            wexc_cur=wexc(iexc)
            if(wexc_adapt)wexc_cur=wexc_cur-wg/2d0
            if(gamma==0d0)then !FFR
               f=2*ens(i)/(ens(i)**2-wexc_cur**2)
               f2=-2*wexc_cur/(ens(i)**2-wexc_cur**2)
               ffr=.true.
            else !in resonance, Lorentzian, currently accepted approach (2026)
               f=1d0/(ens(i)-wexc_cur-iu*gamma)
               wr=(wexc_cur-wg)
               if(st)f2=1d0/(ens(i)+(wexc_cur-wg)+iu*gamma)
            end if
            if(ffr)then
               do a = 1,3
                  do b = 1,3
                     call TDPolar_derivs(iq,i,nq,n3,nst,a,b,duu_1,duu_2,dum_1,dum_2,dmu_1,dmu_2,u,m,u0,m0,step,cen,nfiles,smat)
                     res(iexc)%ap(a,b)=res(iexc)%ap(a,b) + duu_1*f
                     res(iexc)%G(a,b)=res(iexc)%G(a,b) + dum_1*iu*f2
                     res(iexc)%Gc(a,b)=res(iexc)%Gc(a,b) + dmu_1*iu*f2
                     do c = 1,3
                        call TDPolar_derivsQ(iq,i,nq,n3,nst,a,b,c,duq_1,duq_2,dqu_1,dqu_2,u,q,u0,q0,step,cen,nfiles,smat)
                        res(iexc)%A(a,b,c)=res(iexc)%A(a,b,c)+duq_1*f
                        res(iexc)%Ac(a,b,c)=res(iexc)%Ac(a,b,c)+dqu_1*f
                     end do
                  end do
               end do
            else
               do a = 1,3
                  do b = 1,3
                     call TDPolar_derivs(iq,i,nq,n3,nst,a,b,duu_1,duu_2,dum_1,dum_2,dmu_1,dmu_2,u,m,u0,m0,step,cen,nfiles,smat)
                  
                     res(iexc)%ap(a,b)=res(iexc)%ap(a,b) + conjg(duu_1*f + duu_2*f2)
                     ! res(iexc)%G(a,b)=res(iexc)%G(a,b) + u(a,i)*m(b,i)*iu*f - m(b,i)*u(a,i)*iu*f2
                     ! res(iexc)%Gc(a,b)=res(iexc)%Gc(a,b) - m(a,i)*u(b,i)*iu*f + u(b,i)*m(a,i)*iu*f2
                     res(iexc)%G(a,b)=res(iexc)%G(a,b) + conjg(dum_1*iu*f + dmu_2*iu*f2)
                     res(iexc)%Gc(a,b)=res(iexc)%Gc(a,b) + conjg(dmu_1*iu*f + dum_2*iu*f2)
                     do c = 1,3
                        call TDPolar_derivsQ(iq,i,nq,n3,nst,a,b,c,duq_1,duq_2,dqu_1,dqu_2,u,q,u0,q0,step,cen,nfiles,smat)
                        res(iexc)%A(a,b,c)=res(iexc)%A(a,b,c)+conjg(duq_1*f + dqu_2*f2)
                        res(iexc)%Ac(a,b,c)=res(iexc)%Ac(a,b,c)+conjg(dqu_1*f + duq_2*f2)
                     end do
                  end do
               end do
            end if
         end do
      end do
      do iexc=1,nexc
         res(iexc)%ap=res(iexc)%ap*sqrt_w
         res(iexc)%G=res(iexc)%G*sqrt_w
         res(iexc)%Gc=res(iexc)%Gc*sqrt_w
         res(iexc)%A=res(iexc)%A*sqrt_w
         res(iexc)%Ac=res(iexc)%Ac*sqrt_w
      end do
   end subroutine TDPolar_nograd
   
   subroutine TDPolar_derivs(iq,ist,nq,n3,nst,a,b,duu_1,duu_2,dum_1,dum_2,dmu_1,dmu_2,u_all,m_all,u0,m0,step,cen,nfiles,smat)
      double precision,intent(in) :: u_all(3,nst,nfiles),m_all(3,nst,nfiles),u0(3,nst),m0(3,nst)
      double precision,intent(in) :: step,smat(n3,nq)
      integer,intent(in) :: nfiles,ist,a,b,n3,iq,nq,nst
      logical,intent(in) :: cen
      double precision,intent(out) :: duu_1,duu_2,dum_1,dum_2,dmu_1,dmu_2
      double precision :: dsdri !d(something)/dr
      integer i
      
      duu_1=0d0
      duu_2=0d0
      dum_1=0d0
      dum_2=0d0
      dmu_1=0d0
      dmu_2=0d0
      if(cen)then
         do i = 1,n3
            !<n|u_a|j>*<j|u_b|n>
            dsdri=cendiff_r(u_all(a,ist,i)*u_all(b,ist,i),u_all(a,ist,i+n3)*u_all(b,ist,i+n3),step)
            duu_1=duu_1+smat(i,iq)*dsdri
            
            !<n|u_a|j>*<j|m_b|n>
            dsdri=cendiff_r(u_all(a,ist,i)*m_all(b,ist,i),u_all(a,ist,i+n3)*m_all(b,ist,i+n3),step)
            dum_1=dum_1+smat(i,iq)*dsdri
            
            dsdri=cendiff_r(-m_all(a,ist,i)*u_all(b,ist,i),-m_all(a,ist,i+n3)*u_all(b,ist,i+n3),step)
            dmu_1=dmu_1+smat(i,iq)*dsdri
         end do
      else
         do i = 1,n3
            !<n|u_a|j>*<j|u_b|n>
            dsdri=diff_r(u0(a,ist)*u0(b,ist),u_all(a,ist,i)*u_all(b,ist,i),step)
            duu_1=duu_1+smat(i,iq)*dsdri
            
            !<n|u_a|j>*<j|m_b|n>
            dsdri=diff_r(u0(a,ist)*m0(b,ist),u_all(a,ist,i)*m_all(b,ist,i),step)
            dum_1=dum_1+smat(i,iq)*dsdri
            
            dsdri=diff_r(-m0(a,ist)*u0(b,ist),-m_all(a,ist,i)*u_all(b,ist,i),step)
            dmu_1=dmu_1+smat(i,iq)*dsdri
         end do
      end if
      duu_2=duu_1
      dmu_2=-dum_1
      dum_2=-dmu_1
   end subroutine TDPolar_derivs
   
   subroutine TDPolar_derivsQ(iq,ist,nq,n3,nst,a,b,c,duq_1,duq_2,dqu_1,dqu_2,u_all,q_all,u0,q0,step,cen,nfiles,smat)
      double precision,intent(in) :: u_all(3,nst,nfiles),q_all(3,3,nst,nfiles),u0(3,nst),q0(3,3,nst)
      double precision,intent(in) :: step,smat(n3,nq)
      integer,intent(in) :: nfiles,ist,a,b,c,n3,nst,nq,iq
      logical,intent(in) :: cen
      double precision,intent(out) :: duq_1,duq_2,dqu_1,dqu_2
      double precision :: dsdri !d(something)/dr
      integer i
      
      duq_1=0d0
      duq_2=0d0
      dqu_1=0d0
      dqu_2=0d0
      if(cen)then
         do i = 1,n3
            dsdri=cendiff_r(u_all(a,ist,i)*q_all(b,c,ist,i),u_all(a,ist,i+n3)*q_all(b,c,ist,i+n3),step)
            duq_1=duq_1+smat(i,iq)*dsdri
            
            dsdri=cendiff_r(q_all(b,c,ist,i)*u_all(a,ist,i),q_all(b,c,ist,i+n3)*u_all(a,ist,i+n3),step)
            dqu_1=dqu_1+smat(i,iq)*dsdri
         end do
      else
         do i = 1,n3
            !<n|u_a|j>*<j|u_b|n>
            dsdri=diff_r(u0(a,ist)*q0(b,c,ist),u_all(a,ist,i)*q_all(b,c,ist,i),step)
            duq_1=duq_1+smat(i,iq)*dsdri
            
            dsdri=diff_r(q0(b,c,ist)*u0(a,ist),q_all(b,c,ist,i)*u_all(a,ist,i),step)
            dqu_1=dqu_1+smat(i,iq)*dsdri
         end do
      end if
      duq_2=duq_1
      dqu_2=dqu_1
   end subroutine TDPolar_derivsQ
   

   function ElPoldo(wg,nat,w_jn,wexc,gamma,u,v,m,q,du,dv,dm,dq,grad,grad_gr,grad_ex,vel)result(polars)
      logical grad,vel
      integer nat,dl_idx,nq_orig,n3,iz,q_idx
      integer a,b,c,i,iexc
      double precision w_jn,wexc,gamma,sqrt_w
      double precision u(3),v(3),m(3),q(3,3)
      double precision grad_ex(3*nat),grad_gr(3*nat)
      double precision :: wg,wgg
      double complex :: f,f2
      double precision :: du(3,3*nat),dm(3,3*nat),dv(3,3*nat),dq(3,3,3*nat)
      type(Polar),target :: polars(3*nat)
      type(Polar),pointer :: pol
      double complex :: df,df2
            
      
      
      if(vel)then
         u=v
         du=dv
      end if

      df=0
      df2=0
      do i = 1,3*nat
         call Polar_assign_d(polars(i),0d0)
         sqrt_w=1d0/sqrt(2d0*wg)
         pol=>polars(i)
         f=((w_jn-wexc-iu*gamma))
         f2=((w_jn+(wexc-wg)+iu*gamma))
         if(grad)then
            df=-1d0/(f**2)*(grad_ex(i)-grad_gr(i))
            df2=-1d0/(f2**2)*(grad_ex(i)-grad_gr(i))
         end if
         f=1d0/f
         f2=1d0/f2
         
         do a = 1,3
            do b = 1,3
               pol%ap(a,b)=(du(a,i)*u(b)*f + u(a)*du(b,i)*f + u(a)*u(b)*df + &
                                du(b,i)*u(a)*f2 + u(b)*du(a,i)*f2 + u(b)*u(a)*df2)*sqrt_w
               pol%G(a,b)=(-du(a,i)*m(b)*f - u(a)*dm(b,i)*f - u(a)*m(b)*df + &
                               dm(b,i)*u(a)*f2 + m(b)*du(a,i)*f2 + m(b)*u(a)*df2)*iu*sqrt_w
               pol%Gc(a,b)=(dm(a,i)*u(b)*f + m(a)*du(b,i)*f + m(a)*u(b)*df - &
                                du(b,i)*m(a)*f2 - u(b)*dm(a,i)*f2 - u(b)*m(a)*df2)*iu*sqrt_w
               do c = 1,3
                  pol%A(a,b,c)=(du(a,i)*q(b,c)*f + u(a)*dq(b,c,i)*f + u(a)*q(b,c)*df + &
                                    dq(b,c,i)*u(a)*f2 + q(b,c)*du(a,i)*f2 + q(b,c)*u(a)*df2)*sqrt_w
                  pol%Ac(a,b,c)=(u(a)*dq(b,c,i)*f + du(a,i)*q(b,c)*f + u(a)*q(b,c)*df + &
                                    q(b,c)*du(a,i)*f2 + dq(b,c,i)*u(a)*f2 + q(b,c)*u(a)*df2)*sqrt_w !TODO check
               end do
            end do
         end do
      end do
      
   end function ElPoldo
   
   function CenDiff_pols(pols,n3,nexc,step)result(dpol)
      integer n3,nexc,iexc,a,b,c,i
      double precision :: step
      type(Polar) pols(nexc,2*n3),dpol(nexc,n3)
      
      
      do iexc=1,nexc
         do i = 1,n3
            do a = 1,3
               do b = 1,3
                  dpol(iexc,i)%ap(a,b)=cendiff(pols(iexc,i)%ap(a,b),pols(iexc,i+n3)%ap(a,b),step)
                  dpol(iexc,i)%G(a,b)=cendiff(pols(iexc,i)%G(a,b),pols(iexc,i+n3)%G(a,b),step)
                  dpol(iexc,i)%Gc(a,b)=cendiff(pols(iexc,i)%Gc(a,b),pols(iexc,i+n3)%Gc(a,b),step)
                  do c = 1,3
                     dpol(iexc,i)%A(a,b,c)=cendiff(pols(iexc,i)%A(a,b,c),pols(iexc,i+n3)%A(a,b,c),step)
                     dpol(iexc,i)%Ac(a,b,c)=cendiff(pols(iexc,i)%Ac(a,b,c),pols(iexc,i+n3)%Ac(a,b,c),step)
                  end do
               end do
            end do
         end do
      end do
   end function CenDiff_pols
   
   function CenDiff2_pols(pols0,pols,n3,nexc,step)result(d2pol)
      integer n3,nexc,iexc,a,b,c,i
      double precision :: step
      type(Polar) pols(nexc,2*n3),d2pol(nexc,n3),pols0(nexc)
      
      
      do iexc=1,nexc
         do i = 1,n3
            do a = 1,3
               do b = 1,3
                  d2pol(iexc,i)%ap(a,b)=cendiff2(pols0(iexc)%ap(a,b),pols(iexc,i)%ap(a,b),pols(iexc,i+n3)%ap(a,b),step)
                  d2pol(iexc,i)%G(a,b)=cendiff2(pols0(iexc)%G(a,b),pols(iexc,i)%G(a,b),pols(iexc,i+n3)%G(a,b),step)
                  d2pol(iexc,i)%Gc(a,b)=cendiff2(pols0(iexc)%Gc(a,b),pols(iexc,i)%Gc(a,b),pols(iexc,i+n3)%Gc(a,b),step)
                  do c = 1,3
                     d2pol(iexc,i)%A(a,b,c)=cendiff2(pols0(iexc)%A(a,b,c),pols(iexc,i)%A(a,b,c),pols(iexc,i+n3)%A(a,b,c),step)
                     d2pol(iexc,i)%Ac(a,b,c)=cendiff2(pols0(iexc)%Ac(a,b,c),pols(iexc,i)%Ac(a,b,c),pols(iexc,i+n3)%Ac(a,b,c),step)
                  end do
               end do
            end do
         end do
      end do
   end function CenDiff2_pols
   
   function CenDiff11_pols(pols0,pols,d2pol,n3,nexc,step)result(d11pol)
      integer n3,nexc,iexc,a,b,c,i,ix,j,add
      double precision :: step
      type(Polar) pols(nexc,2*n3),d11pol(nexc,n3*(n3-1)/2),pols0(nexc),d2pol(nexc,n3)
      
      
      add=n3*(n3-1)/2
      do iexc=1,nexc
         do i = 1,n3
            do j = i+1,n3
               ix=(i-1)*n3+j
               do a = 1,3
                  do b = 1,3
                     d11pol(iexc,ix)%ap(a,b)=cendiff_mixed(pols0(iexc)%ap(a,b),pols(iexc,ix)%ap(a,b),pols(iexc,ix+add)%ap(a,b),d2pol(iexc,i)%ap(a,b),d2pol(iexc,j)%ap(a,b),step,step)
                     d11pol(iexc,ix)%G(a,b)=cendiff_mixed(pols0(iexc)%G(a,b),pols(iexc,ix)%G(a,b),pols(iexc,ix+add)%G(a,b),d2pol(iexc,i)%G(a,b),d2pol(iexc,j)%G(a,b),step,step)
                     d11pol(iexc,ix)%Gc(a,b)=cendiff_mixed(pols0(iexc)%Gc(a,b),pols(iexc,ix)%Gc(a,b),pols(iexc,ix+add)%Gc(a,b),d2pol(iexc,i)%Gc(a,b),d2pol(iexc,j)%Gc(a,b),step,step)
                     do c = 1,3
                        d11pol(iexc,ix)%A(a,b,c)=cendiff_mixed(pols0(iexc)%A(a,b,c),pols(iexc,ix)%A(a,b,c),pols(iexc,ix+add)%A(a,b,c),d2pol(iexc,i)%A(a,b,c),d2pol(iexc,j)%A(a,b,c),step,step)
                        d11pol(iexc,ix)%Ac(a,b,c)=cendiff_mixed(pols0(iexc)%Ac(a,b,c),pols(iexc,ix)%Ac(a,b,c),pols(iexc,ix+add)%Ac(a,b,c),d2pol(iexc,i)%Ac(a,b,c),d2pol(iexc,j)%Ac(a,b,c),step,step)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end function CenDiff11_pols
   
   function Diff_pols(pols0,pols,n3,nexc,step)result(dpol)
      integer n3,nexc,iexc,a,b,c,i
      double precision :: step
      type(Polar) pols(nexc,n3),pols0(nexc),dpol(nexc,n3)
      
      
      do iexc=1,nexc
         do i = 1,n3
            do a = 1,3
               do b = 1,3
                  dpol(iexc,i)%ap(a,b)=diff(pols0(iexc)%ap(a,b),pols(iexc,i)%ap(a,b),step)
                  dpol(iexc,i)%G(a,b)=diff(pols0(iexc)%G(a,b),pols(iexc,i)%G(a,b),step)
                  dpol(iexc,i)%Gc(a,b)=diff(pols0(iexc)%Gc(a,b),pols(iexc,i)%Gc(a,b),step)
                  do c = 1,3
                     dpol(iexc,i)%A(a,b,c)=diff(pols0(iexc)%A(a,b,c),pols(iexc,i)%A(a,b,c),step)
                     dpol(iexc,i)%Ac(a,b,c)=diff(pols0(iexc)%Ac(a,b,c),pols(iexc,i)%Ac(a,b,c),step)
                  end do
               end do
            end do
         end do
      end do
   end function Diff_pols
   
   function diff(fx,fxh,h)result(res)
      double complex fx,fxh,res
      double precision h
      res=(fxh-fx)/h
   end function diff
   
   function diff_R(fx,fxh,h)result(res)
      double precision fx,fxh,res
      double precision h
      res=(fxh-fx)/h
   end function diff_R
   
   function cendiff(fxhp,fxhm,h)result(res)
      double precision h
      double complex fxhm,fxhp,res
      res=(fxhp-fxhm)/(2d0*h)
   end function cendiff
   
   function cendiff2(fx,fxhp,fxhm,h)result(res)
      double precision h
      double complex fx,fxhm,fxhp,res
      res=(fxhp-2*fx+fxhm)/(h**2)
   end function cendiff2
   
   function cendiff_mixed(f0,f_xp_yp,f_xm_ym,d2f_x2,d2f_y2,dx,dy)result(res)
      double precision dx,dy
      double complex f0,f_xp_yp,f_xm_ym,d2f_x2,d2f_y2,res
      res=(f_xp_yp+f_xm_ym-2*f0-d2f_x2*dx**2-d2f_y2*dy**2)/(2*dx*dy)
   end function cendiff_mixed
   
   function cendiff_R(fxhp,fxhm,h)result(res)
      double precision h
      double precision fxhm,fxhp,res
      res=(fxhp-fxhm)/(2d0*h)
   end function cendiff_r
   
end program rroa_td_num