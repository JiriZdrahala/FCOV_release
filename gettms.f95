#define TR(arg) trim(adjustl(arg))
program gettms
   use iso_fortran_env
   use constants
   use strings
   implicit none
   integer,parameter :: char_len=200
   character(char_len) filee,finp,duschout,s80,s80_l
   integer n3,nat,nq,line,nroot,argc,iz,del_modes_c,optend,i,ncar
   
   logical :: file_found=.false.
   logical :: islinear = .false.,transform_vel=.false.
   logical :: isAH=.false.,AHAS=.false.,vel_grad=.true.,coordinates=.true.
   logical :: finp_found=.false.,dusch_found=.false.,fex,ground=.true.
   
   double precision :: v(3),u(3),m(3),q(6),q_tr(3,3),q_new(3,3)
   double precision :: eau,e_00,e_gr,e_ex
   double precision :: mult_tm(5),mult_dtm(5)
   double precision,allocatable :: smat(:,:),r(:),wg(:)
   double precision,allocatable :: aa(:),vv(:),dd(:),qq(:)
   double precision,allocatable :: aa_q(:),vv_q(:),dd_q(:),qq_q(:)
   double precision,allocatable :: du(:,:),dv(:,:),dm(:,:),dq(:,:)
   
   logical :: grad_e_found=.false.,grad_g_found=.false.
   character(char_len) :: grad_e_file,grad_g_file
   double precision,allocatable :: grad_e(:),grad_g(:)
   
   
   mult_tm=[1d0,1d0,0.5d0,1d0,1d0]
   mult_dtm=[1d0,1d0,0.5d0,1d0,1d0]
   e_gr=0
   e_ex=0
   e_00=0
   finp='F.INP'
   duschout='DUSCH.OUT'
   !call ReadArgs(argc,optend)
   argc=COMMAND_ARGUMENT_COUNT()
   do i = 1,argc
      call GET_COMMAND_ARGUMENT(i,s80)
      if(s80(1:1)=='-')then !OPTION
         s80_l=s80
         call To_Lower(s80_l)
         if(TR(s80_l(2:))=='ground')then
            ground=.true.
         else if(TR(s80_l(2:))=='excited')then
            ground=.false.
         else if(TR(s80_l(2:))=='linear')then
            islinear=.true.
         else if(TR(s80_l(2:5))=='e_eg=')then
            read(s80(6:),*)e_00
         else if(TR(s80_l(2:6))=='trvel')then
            transform_vel=.true.
         else
            write(output_unit,*)'ERROR: Unknown option '//TR(s80)
            call exit(1)
         end if
         cycle
      end if
      
      if(index(s80,'F.INP')>0)then
         finp=s80
         finp_found=.true.
      else if(index(s80,'DUSCH.OUT')>0)then
         duschout=s80
         dusch_found=.true.
      else
         if(file_found)then
            write(output_unit,*)'ERROR: Only 1 "#p TD-freq" file can be specified'
            call exit(3)
         end if
         inquire(file=s80,exist=fex)
         if(.not.fex)then
            write(output_unit,*)'ERROR: Unknown option or file does not exist '//TR(s80)
            call exit(98)
         else
            filee=s80
            file_found=.true.
         end if
      end if
   end do
   if(.not.file_found)then
      write(output_unit,*)"USAGE:"
      write(output_unit,*)"gettms [[-ground -excited -linear]] (DUSCH.OUT/F.INP filepath) (#p-TD-freq file)"
      call exit(99)
   end if
   if(dusch_found.and.finp_found)then
      write(output_unit,*)'ERROR: Specify F.INP or DUSCH.OUT, not both.'
      call exit(2)
   end if
   
   open(77,file='ground/ENERGY',status='old')
   read(77,*)e_gr
   close(77)
   write(output_unit,*)'Read ground/ENERGY'
   open(77,file='excited/ENERGY',status='old')
   read(77,*)e_ex
   close(77)
   write(output_unit,*)'Read excited/ENERGY'

   if(e_00==0d0)then
      e_00=e_ex-e_gr
   end if
   
   nat=readnat(77,filee)
   ncar=9*nat
   n3=3*nat
   
   grad_g_file='ground/FILE.GR'
   grad_e_file='excited/FILE.GR'
   if(transform_vel)grad_g=readfilegr(77,nat,grad_g_file)
   if(transform_vel)grad_e=readfilegr(77,nat,grad_e_file)
      
   if(finp_found)then
      call readsi(n3,smat,wg,nq,finp,r,.true.,iz)
   else if(dusch_found)then
      smat=readDusch(nq,n3,ground,duschout)
   end if
   allocate(aa(ncar),vv(ncar),dd(ncar),qq(2*ncar))
   line=0
   call rdd(filee,dd,nat,u,m,aa,vv,v,qq,q,line,nroot,eau,mult_tm,mult_dtm)
   
   
   
   if(transform_vel)then
      if(e_00==0d0)then
         write(output_unit,*)'ERROR: Velocity transformation to length was requested but transition energy is 0'
         call exit(55)
      end if
      vv=VelToLen_dervec_cart(3,v,vv,e_00,grad_e,grad_g,nat,ncar)
      v=VelToLen_vec_cart(v,e_00)
      qq=VelToLen_dervec_cartQ(6,q,qq,e_00,grad_e,grad_g,nat,2*ncar)
      q=VelToLen_vec_cartq(q,e_00)
   end if
   
   q_tr=TracelessQ(q)
   q_new=NewQ(q)
   call WriteVec(77,'TM.U',u)
   call WriteVec(77,'TM.M',m)
   call WriteVec(77,'TM.V',v)
   call WriteTen(77,'TM.Q',q_new)
   call WriteTen(77,'TM.Qtr',q_tr)
   
   call WriteDerVec(77,'DTM.U',dd,nat,ncar)
   call WriteDerVec(77,'DTM.M',aa,nat,ncar)
   call WriteDerVec(77,'DTM.V',vv,nat,ncar)
   call WriteDerQ(77,'DTM.Q',qq,nat,2*ncar)
   call WriteDerQtr(77,'DTM.Qtr',qq,nat,2*ncar)
   
   if(.not.allocated(smat))call exit(0)
   allocate(aa_q(ncar),vv_q(ncar),dd_q(ncar),qq_q(2*ncar))
   call trafN_new(3,dd,dd_q,nat,nq,smat)
   call trafN_new(3,aa,aa_q,nat,nq,smat)
   call trafN_new(6,qq,qq_q,nat,nq,smat)
   call trafN_new(3,vv,vv_q,nat,nq,smat)
   deallocate(aa,vv,dd,qq)
   
   ! allocate(du(3,nq),dv(3,nq),dm(3,nq),dq(6,nq))
   call TidyDerivatives(dd_q,ncar,du,3,NQ)
   call TidyDerivatives(aa_q,ncar,dm,3,NQ)
   call TidyDerivativesQ(qq_q,2*ncar,dq,6,NQ)
   call TidyDerivatives(vv_q,ncar,dv,3,NQ)
   deallocate(dd_q,aa_q,qq_q,vv_q)
   
   ! if(transform_vel)then
      ! call VelToLen_dip(v,dv,e_00,nq,isAH,AHAS,vel_grad,coordinates)
      ! call VelToLen_quad(q_tr,dq,e_00,nq,isAH,AHAS,vel_grad,coordinates)
   ! end if
   
   call WriteDerVec_normalmode(77,'DTM.Q.U',du,nat,nq)
   call WriteDerVec_normalmode(77,'DTM.Q.M',dm,nat,nq)
   call WriteDerVec_normalmode(77,'DTM.Q.V',dv,nat,nq)
   call WriteDerVecQ_normalmode(77,'DTM.Q.Q',dq,nat,nq)
   call WriteDerVecQtr_normalmode(77,'DTM.Q.Qtr',dq,nat,nq)
   
   
   deallocate(du,dv,dm,dq)
   
   contains
   
   function readfilegr(unitt,nat,filee)result(gr)
      integer nat,unitt,ix,i
      character(*) filee
      double precision gr(nat*3)
      
      open(unitt,file=filee,status='old')
      do i = 1,nat
         ix=(i-1)*3+1
         read(unitt,*)gr(ix:ix+2)
      end do
      close(unitt)
   end function readfilegr
   
   function readnat(unitt,filee)result(nat)
      character(*) filee
      integer nat,unitt
      character(80) s80
      
      open(unitt,file=filee,status='old')
100   read(unitt,'(A80)',end=200)s80
      if(s80(2:8)=='NAtoms=')then
         if(s80(17:17)=='N')goto 100
         read(s80(9:13),*)nat
         goto 200
      end if
      goto 100
200   close(unitt)
   end function readnat
   
   subroutine WriteDerVec(unitt,namee,vecs,nat,n)
      integer unitt,n,i,nat,ix
      character(*) namee
      double precision vecs(n)
      
      open(unitt,file=namee)
      write(unitt,*)nat*3
      do i = 1,nat*3
         ix=(i-1)*3+1
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,vecs(ix:ix+2)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVec
   
   
   subroutine WriteDerVec_normalmode(unitt,namee,vecs,nat,nq)
      integer unitt,n,i,nat,ix,nq
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
   
   subroutine WriteDerVecQ_normalmode(unitt,namee,vecs,nat,nq)
      integer unitt,n,i,nat,ix,nq
      character(*) namee
      double precision vecs(6,nq),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nq
      do i = 1,nq
         q=NewQ(vecs(:,i))
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')  q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')  q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVecQ_normalmode
   
   subroutine WriteDerVecQtr_normalmode(unitt,namee,vecs,nat,nq)
      integer unitt,n,i,nat,ix,nq
      character(*) namee
      double precision vecs(6,nq),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nq
      do i = 1,nq
         q=TracelessQ(vecs(:,i))
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerVecQtr_normalmode
   
   subroutine WriteDerQ(unitt,namee,vecs,nat,n)
      integer unitt,n,i,nat,ix
      character(*) namee
      double precision vecs(n),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nat*3
      do i = 1,nat*3
         ix=(i-1)*6+1
         q=NewQ(vecs(ix:ix+5))
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerQ
   
   subroutine WriteDerQtr(unitt,namee,vecs,nat,n)
      integer unitt,n,i,nat,ix
      character(*) namee
      double precision vecs(n),q(3,3)
      
      open(unitt,file=namee)
      write(unitt,*)nat*3
      do i = 1,nat*3
         ix=(i-1)*6+1
         q=TracelessQ(vecs(ix:ix+5))
         write(unitt,'(1X,I6,1X,3(E15.8E2,1X))')i,q(1,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(2,:)
         write(unitt,'(1X,6X,1X,3(E15.8E2,1X))')q(3,:)
      end do
      close(unitt)
      write(output_unit,*)'Written: '//TR(namee)
   end subroutine WriteDerQtr
      
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
   
   subroutine ReadArgs(argc,optEnd)
      character(80) s80,s80_2
      integer i,argc,optEnd,idx_eq
      
      
      do i = 1,argc
         call GET_COMMAND_ARGUMENT(i,s80)
         call To_lower(s80)
         
         if(s80(1:1)=='-')then
            idx_eq=index(s80,'=')+1
            
            if(idx_eq==1)then
               s80_2=s80(1:80)
            else
               s80_2=s80(1:idx_eq-1)
            end if
            select case(TR(s80))
               case('')
               
               case default
                  write(output_unit,*)'Unknown option '//s80
            end select
         else
            optEnd=i-1
            return
         end if
      end do
      optEnd=argc
   end subroutine ReadArgs
   
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
   end
   
   subroutine report(s)
   character*(*) s
   write(6,*)s
   stop
   end
   
   SUBROUTINE readsi(N3,S,E,NQ,fn,r,ldz,iz)
!     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,iz,fuck1,fuck2
      real*8 CM
      real*8,allocatable :: S(:,:), S_help(:,:)
      real*8,allocatable :: E(:),E_help(:)
      double precision,allocatable :: r(:)
      logical ldz
      character*(*) fn
      
      CM=219474.630d0
!     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn,status='old')
      read(4,*)NQ,nat,nat
      N3=3*nat
      allocate(S(n3,n3),E(nq),r(n3))
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
      !write(6,*)NQ,' modes found'
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
         !write(6,*)iz,' zero modes: deleted'
         allocate(S_help(N3,NQ),E_help(NQ))
         S_help=S(:,1:NQ)
         E_help=E(1:NQ)
         deallocate(S,E)
         S=S_help/sqrt(amu_2_au)
         E=E_help
         deallocate(S_help,E_help)
       endif 
      endif
      
      
      !write(6,*)NQ,' vibrational modes considered'
      
      DO 3 I=1,NQ
3     E(I)=E(I)/CM
      
      RETURN
   end SUBROUTINE readsi
   
   function TracelessQ(q_6)result(q_tr)
      double precision q_6(6),q_tr(3,3)
      
      !Q
      !1  2  3  4  5  6
      !xx xy xz yy yz zz
      q_tr(1,1)=q_6(1)-0.5d0*(q_6(4)+q_6(6))
      q_tr(2,2)=q_6(4)-0.5d0*(q_6(1)+q_6(6))
      q_tr(3,3)=q_6(6)-0.5d0*(q_6(1)+q_6(4))
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
   
   subroutine TidyDerivatives(du,n,du_new,a,nq)
      integer n,nq,a,idx,i,i_new
      double precision du(n)
      double precision, allocatable :: du_new(:,:)
      
      allocate(du_new(a,nq))
      i_new=1
      do i = 1,nq
         !if(findloc(del_modes,i,dim=1)>0)cycle
         idx=(i_new-1)*a+1
         du_new(:,i_new)=du(idx:idx+a-1)
         i_new=i_new+1
      end do
   end subroutine TidyDerivatives
   
   subroutine TidyDerivativesQ(dq,n,dq_new,a,nq)
      integer n,nq,a,idx,i,i_new
      double precision dq(n)
      double precision, allocatable :: dq_new(:,:)
      
      allocate(dq_new(a,nq))
      i_new=1
      do i = 1,nq
         idx=(i_new-1)*a+1
         dq_new(:,i_new)=dq(idx:idx+a-1)
         i_new=i_new+1
      end do
   end subroutine TidyDerivativesQ
   
   function VelToLen_vec_cart(v,e_eg)result(res)
      double precision v(3),e_eg,res(3)
      res=-v/e_eg
   end function VelToLen_vec_cart
   
   function VelToLen_dervec_cart(dimm,v,dv,e_eg,grad_e,grad_g,nat,n)result(res)
      integer nat,n,ix,dimm
      double precision dv(n),v(dimm),e_eg,grad_e(3*nat),grad_g(3*nat),res(n)
      
      do i = 1,nat*3
         ix=(i-1)*dimm+1
         res(ix:ix+dimm-1)=v/e_eg**2 * (grad_e(i)-grad_g(i)) - dv(ix:ix+dimm-1)/e_eg
      end do
   end function VelToLen_dervec_cart
   
   function VelToLen_vec_cartQ(Q,e_eg)result(res)
      double precision q(6),e_eg,res(6)
      res=-q/e_eg
   end function VelToLen_vec_cartQ
   
   function VelToLen_dervec_cartQ(dimm,q,dq,e_eg,grad_e,grad_g,nat,n)result(res)
      integer nat,n,ix,dimm
      double precision dq(n),q(dimm),e_eg,grad_e(3*nat),grad_g(3*nat),res(n)
      
      do i = 1,nat*3
         ix=(i-1)*dimm+1
         res(ix:ix+dimm-1)=q/e_eg**2 * (grad_e(i)-grad_g(i)) - dq(ix:ix+dimm-1)/e_eg
      end do
   end function VelToLen_dervec_cartQ
   
   subroutine VelToLen_dip(v,dv,e00,nq,isAH,AHAS,vel_grad,coordinates)
      integer nq,i
      double precision v(3),dv(3,nq),e00,e00_copy
      double precision gr_q(nq),gr_q2(nq)
      logical isAH,AHAS,vel_grad,coordinates
      
      e00_copy=e00
      if(isAH)then !adiabatic
         open(77,file='ENERGY.VERTICAL')
         read(77,*)e00_copy
         close(77)
         if(AHAS)then
            call readgradQ(nq,gr_q,'FILE.Q.GR.excited')
            call readgradQ(nq,gr_q2,'FILE.Q.GR.ground.excited')
            if(.not.vel_grad)then
               gr_q=0
               gr_q2=0
            end if
            do i = 1,nq
               dv(:,i)=(-v(:)/e00_copy**2)*(gr_q(i)-gr_q2(i))+dv(:,i)/e00_copy
            end do
            v=v/e00_copy
            return
         end if
         if(coordinates)then
            call readgradq(nq,gr_q,'FILE.Q.GR.ground.excited')
         else
            call readgradQ(nq,gr_q,'FILE.Q.GR.EXTRAP')
         end if
         if(.not.vel_grad)gr_q=0
         gr_q=-gr_q
      else !vertical
         if(coordinates)then
            call readgradq(nq,gr_q,'FILE.Q.GR.excited.ground')
         else
            call readgradQ(nq,gr_q,'FILE.Q.GR.excited')
         end if
      end if
      !u=v/e00_copy
      if(.not.vel_grad)gr_q=0
      do i = 1,nq
         dv(:,i)=(-v(:)/e00_copy**2)*gr_q(i)+dv(:,i)/e00_copy
      end do
      v=v/e00_copy
   end subroutine VelToLen_dip
   
   subroutine VelToLen_quad(q,dq,e00,nq,isAH,AHAS,vel_grad,coordinates)
      integer nq,i
      double precision q(3,3),dq(3,3,nq),e00,e00_copy
      double precision gr_q(nq),gr_q2(nq)
      logical isAH,AHAS,vel_grad,coordinates
      
      !dq_help=dv
      e00_copy=e00
      if(isAH)then
         open(77,file='ENERGY.VERTICAL')
         read(77,*)e00_copy
         close(77)
         if(AHAS)then
            call readgradQ(nq,gr_q,'FILE.Q.GR.excited')
            call readgradQ(nq,gr_q2,'FILE.Q.GR.excited.ground')
            if(.not.vel_grad)then
               gr_q=0
               gr_q2=0
            end if
            do i = 1,nq
               dq(:,:,i)=(-q(:,:)/e00_copy**2)*(gr_q(i)-gr_q2(i))+dq(:,:,i)/e00_copy
            end do
            q=q/e00_copy
            return
         end if
         if(coordinates)then
            call readgradq(nq,gr_q,'FILE.Q.GR.ground.excited')
         else
            call readgradQ(nq,gr_q,'FILE.Q.GR.EXTRAP')
         end if
         if(.not.vel_grad)gr_q=0
         gr_q=-gr_q
      else
         if(coordinates)then
            call readgradq(nq,gr_q,'FILE.Q.GR.excited.ground')
         else
            call readgradQ(nq,gr_q,'FILE.Q.GR.excited')
         end if
      end if
      if(.not.vel_grad)gr_q=0
      do i = 1,nq
         dq(:,:,i)=(-q(:,:)/e00_copy**2)*gr_q(i)+dq(:,:,i)/e00_copy
      end do
      q=q/e00_copy
   end subroutine VelToLen_quad
   
   subroutine rdd(f,d,nat,u0,m0,a,v,v0,q,q0,Line,nroot,eau,mult_tm,mult_dtm)
!     read transition dipole derivatives from excited state freq calc
      implicit none
      integer nat,ln,ia,xa,i,xar,iar,ibr,ib,ic,ii,ix,iwr,nd, &
                ntr,j,ir,idx, &
                iLeft,iRight,iTD,iFREQ,iNSTATES
      integer,intent(inout) :: Line !Given, skips number of lines; when returned gives last line read
      integer,intent(out) :: nroot
      double precision d(9*nat),u(6),u0(3),step,sign,stepau,m0(3),a(9*nat), &
             v(9*nat),v0(3),q(18*nat),q0(6),ev,eau,mult_tm(5),mult_dtm(5)
      logical lnumerical,isAH
      character*(*) f
      character*70 s80 !gaussian might change the output file row length for some stupid reason, modify this if they do
      character*200 s200
      character*70,dimension(10) :: s80_arr !Gaussian input parameters
      character*800 :: s800
      double precision,allocatable::v16(:,:)
!
      ! write(6,*)
      ! write(6,*)' Rdd: reading '//f
      ! write(6,*)
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
            ! write(output_unit,*)'TD freq job starts at line',ln
            
            iNSTATES=index(s200,'nstates=')
            iRight=findAnyChar(s200,[')',','],2)
            read(s200(iNSTATES+8:iright-1),*)ntr
            iNSTATES=index(s200,'root=')
            iRight=findAnyChar(s200,[')',','],2)
            read(s200(iNSTATES+5:iright-1),*)nroot
            ! write(output_unit,*)'TD Root is: ',nroot
            goto 2
         end if
      end if
      goto 1
      
88    close(9)
      ! write(output_unit,*)'TD freq not found'
      stop

2     ia=1
      xa=1
      ib=0
      ii=0
      step=0.001d0
      ! write(6,*)'iwr:',iwr
      a=0.0d0
      v=0.0d0
      q=0.0d0
      d=0.0d0
      nd=0
      
4     read(9,901,end=78)s200
901   format(A200)
      if(s200(2:14).eq.'Nuclear step=')then
       ! write(6,*)s200(1:23)
       read(s200(15:23),*)step
      endif
      if(s200(2:24).eq.'Re-enter D2Numr: IAtom=')then
       read(s200(25:27),*)iar
       read(s200(34:34),*)xar
       read(s200(42:43),*)ibr
       if(ia.ne.iar)call report('ia <> iar')
       if(xa.ne.xar)call report('xa <> xar')
       if(ib.ne.ibr)call report('ib <> ibr')
      endif

      if(s200(2:65).eq.'Ground to excited state transition electric dipole moments (Au):')then
       if(ii.eq.0)then
        ! write(6,*)'Zero point electric'
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
        ! if(iwr.gt.1)write(6,*)' d ',ia,xa,ib,iar,xar,ibr
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
        ! write(6,*)'Zero point velocity'
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)v0(1),(v0(ix),ix=1,3)
       else
        ! if(iwr.gt.1)write(6,*)' v ',ia,xa,ib,iar,xar,ibr
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
        ! write(6,*)'Zero point magnetic'
        do i=1,nroot
        read(9,*)
        end do
        read(9,*)m0(1),(m0(ix),ix=1,3)
       else
        ! if(iwr.gt.1)write(6,*)' m ',ia,xa,ib,iar,xar,ibr
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
        ! write(6,*)'Zero point quadrupole'
        do i=1,nroot
        read(9,*)
        end do
!       our order         1    2      3     4     5     6
!       our order        xx    xy    xz    yy    yz    zz
!       gaussian   #     xx    yy    zz    xy    xz    yz
        read(9,*)q0(1),q0(1),q0(4),q0(6),q0(2),q0(3),q0(5)
       else
        ! if(iwr.gt.1)write(6,*)' m ',ia,xa,ib,iar,xar,ibr
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
         ! write(output_unit,*)'Its energy is: ',ev,' eV'
         eau=ev/27.211384205943d0
         if(ia.eq.nat.and.ix.eq. 3 .and.ib .eq. 2)goto 78
         ii=ii+1
         ! write(6,*)'energy read, go to the next point'
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
       ! write(6,*)'Analytical transitions read'
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
       ! write(6,*)'Analytical transition derivatives read'
       lnumerical=.false.
      endif
!     AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
!     checkpoint style
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(s200(1:18).eq.'ETran state values')then
       if(ntr.eq.0)call report('nstates not found')
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
       ! write(6,*)'Analytical transition derivatives read from chk'
       lnumerical=.false.
      endif
!     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      goto 4

78    close(9)

      if(nd.eq.0)call report('Energy not found')

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
      ! open(40,file='DEG.TEN')
      ! write(40,*)'atom coord dipx        dipy        dipz'
      ! do ia=1,nat
      ! do xa=1,3
      ! write(40,400)ia,xa,(d(ix+3*(xa-1)+9*(ia-1)),ix=1,3)
      ! end do
      ! end do
      ! write(40,*)'equilibrium transition electric dipole:'
      ! write(40,400)0,0,(u0(ix),ix=1,3)
      ! write(40,*)'atom coordmdipx       mdipy       mdipz'
      ! do ia=1,nat
      ! do xa=1,3
      ! write(40,400)ia,xa,(a(ix+3*(xa-1)+9*(ia-1)),ix=1,3)
      ! end do
      ! end do
      ! write(40,*)'equilibrium transition magnetic dipole:'
      ! write(40,400)0,0,(m0(ix),ix=1,3)
      ! write(40,*)'atom coord Qv xx xy xz yy yz zz:'
      ! do ia=1,nat
      ! do xa=1,3
      ! write(40,400)ia,xa,(q(ix+6*(xa-1)+18*(ia-1)),ix=1,6)
      ! end do
      ! end do
      ! write(40,*)'equilibrium transition quadrupole velocity:'
      ! write(40,400)0,0,(q0(ix),ix=1,6)
      ! write(40,*)'atom coord grad x y z:'
      ! do ia=1,nat
      ! do xa=1,3
      ! write(40,400)ia,xa,(v(ix+3*(xa-1)+9*(ia-1)),ix=1,3)
      ! end do
      ! end do
      ! write(40,*)'equilibrium transition velocity:'
      ! write(40,400)0,0,(v0(ix),ix=1,3)
! 400   format(i5,i2,6g12.4)
      ! close(40)
      ! write(6,*)'DEG.TEN written'
      
      line=ln
      return
   end subroutine rdd
   
   function readDusch(nq,N3,ground,filee)result(S)
      character(80) c80
      character(*) filee
      integer bufi,i,col,NQ,N3
      double precision,allocatable :: A(:,:),B(:),C(:,:),D(:),E(:,:),J(:,:),K(:),SG(:,:),SE(:,:),wg(:),we(:),S(:,:)
      integer,parameter :: unitt=77
      logical ground
      
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
      
      deallocate(A,B,C,D,E,wg,we,J,K)
      
      call Forward(unitt,3)
      call readNonSquareMatrix(unitt,SG,N3,NQ)
      
      call Forward(unitt,3)
      call readNonSquareMatrix(unitt,SE,N3,NQ)
      
      SG=SG/sqrt(amu_2_au)
      SE=SE/sqrt(amu_2_au)
      if(ground)then
         S=SG
      else
         S=SE
      end if
      close(unitt)
   end function readDusch
   
   subroutine Forward(unitt,n)
      integer unitt,n,i
      
      do i = 1,N
         read(unitt,*)
      end do
   end subroutine Forward
   
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
   
end program gettms