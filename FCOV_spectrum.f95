#define TR(arg) (trim(adjustl(arg)))

program FCOV_spectrum
   use iso_fortran_env
   use constants
   use wrram
   use strings
   use util
   implicit none
   integer argc,i,unitt,polar_i,ii,idx,k,idx2
   
   double precision :: wmin=0d0,wmax=4000d0,wexc=0d0,wexc_tmp=0d0,wexc_tmp2=0d0,fwhm=5d0,temp=300,kBT=1d0
   double precision :: wmin_fqttt=0d0,wmax_fqttt=1d6
   integer,allocatable :: modes_fqttt(:)
   logical isModeFqttt,conjug(3)
   double precision :: mult_el=1d0,mult_mag=1d0,mult_quad=1d0,w_tr=huge(1d0)
   double precision :: mult_el_fqttt=1d0,mult_mag_fqttt=1d0,mult_quad_fqttt=1d0
   integer :: npx=4001,nq_fileqttt=0
   logical :: use_gauss=.false.,spectrum_temp=.true.
   integer,parameter :: array_temp_c=1000
   logical :: rroa_spr_do(ns0)
   logical :: useA=.true.,useG=.true.
   logical :: useAc=.true.,useGc=.true.
   logical :: dcp1=.false.,dcp2=.false.
   logical :: aboveFundamentals=.false.
   double precision,allocatable :: rroa_spr(:,:,:),rroa_spr_temp(:,:,:),rroa_inv(:)
   integer,allocatable :: rroa_tabOutputs(:)
   double precision,allocatable :: x(:)
   
   
   type(Polar),allocatable :: polars(:)
   type(Polar) polar_sum,curPolar
   character(200) filee
   character(80) str_v1,str_v3
   character(80) c80,c80_2
   character(20) wexc_str
   
   
   type(polar),allocatable :: polars_fileqttt(:)
   double precision :: polar_fileqttt_coeff=1d0,polars_el_coeff=1d0,spr_coeff=1d0,vibronic_coeff=1d0
   double precision :: coeff_cur=1d0,vib_min=0,vib_max=1d6
   
   type(transition),allocatable :: trs(:),trs_elpol(:)
   type(transition),allocatable :: trs_anharm(:)
   type(transition) :: tr_fake
   integer :: trs_anharm_n,trs_c_elpol=0
   type(transition) :: cur_tr
   integer,parameter :: max_trs = 10000
   integer :: trs_c=0,polarfile_trs_c
   integer :: v1_class,v3_class
   integer :: optEnd
   integer :: fileqttt_c,filepolars_c,polarfile_el_trs_c
      
   integer,allocatable :: v1_pos(:),v1(:),v3_pos(:),v3(:)   
   integer :: v1_n,v3_n
   double precision :: w1,w3,w13,w1_au,w3_au,wexc_au,wr,w,cl,cr
   double precision,allocatable :: wexcs(:),fileqttt_E(:),finp(:,:),finp_qttt(:,:)
   logical endd
   logical,allocatable :: switchPh(:)
   
   allocate(trs(max_trs),trs_elpol(max_trs))
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc==0)then 
10    write(output_unit,*)'OPTIONS:'
      write(output_unit,*)'--npx='
      write(output_unit,*)'--x-min='
      write(output_unit,*)'--x-max='
      write(output_unit,*)'--fwhm='
      write(output_unit,*)'--temp='
      write(output_unit,*)'--lineshape='
      write(output_unit,*)'--wexc='
      write(output_unit,*)'--spectrum-temp'
      write(output_unit,*)'--vibronic-coeff='
      write(output_unit,*)'--qttt-coeff='
      write(output_unit,*)'--polars-el-coeff='
      write(output_unit,*)'--spr-coeff='
      write(output_unit,*)'--el-coeff-ttt='
      write(output_unit,*)'--el-coeff='
      write(output_unit,*)'--mag-coeff='
      write(output_unit,*)'--quad-coeff='
      write(output_unit,*)'--no-G'
      write(output_unit,*)'--no-Gc'
      write(output_unit,*)'--no-A'
      write(output_unit,*)'--no-Ac'
      write(output_unit,*)'--dcp1'
      write(output_unit,*)'--dcp2'
      call exit(1)
   end if
   
   fileqttt_c=0
   filepolars_c=0
   conjug=[.false.,.false.,.false.]
   call ReadArgs(optEnd,trs,fileqttt_c,polar_fileqttt_coeff,trs_c,finp,finp_qttt,switchPh,wexc)
   !if(optEnd==argc)goto 10
   !allocate(wexcs(argc-optEnd))
   !wexcs=0
   do i = optEnd+1,argc
      call GET_COMMAND_ARGUMENT(i,filee)
      
      open(77,file=filee,status='old')
      if(index(TR(filee),'FILE.AHHARM.FREQ')>0)then
         call Readfileanharmfreq(77,trs_anharm,trs_anharm_n)
         close(77)
         cycle
      end if
      
      if(index(TR(filee),'nm.POLARS.EL')>0)then
         coeff_cur=polars_el_coeff
      else
         coeff_cur=vibronic_coeff
      end if
      read(77,*)polarfile_trs_c
      endd=.false.
      do while(.not.endd)
         call ReadTenPretty(77,endd,cur_tr,wexc_tmp)
         if(wexc==0d0)then
            wexc=wexc_tmp
         end if
         if(fileqttt_c>0 .and. ((cur_tr%w3-cur_tr%w1)<vib_min .or. (cur_tr%w3-cur_tr%w1)>vib_max))cycle
         if(.not.endd)then
            idx=FindTransition(cur_tr,trs,max_trs,trs_c)
            ! if(switchPh(idx))coeff_cur=-coeff_cur
            cur_tr%polarr=coeff_cur*cur_tr%polarr
            cur_tr%polarr%ap=cur_tr%polarr%ap*mult_el**2
            cur_tr%polarr%G=cur_tr%polarr%G*mult_mag*mult_el
            cur_tr%polarr%Gc=cur_tr%polarr%Gc*mult_mag*mult_el
            if(conjug(2))then
               cur_tr%polarr%G=conjg(cur_tr%polarr%G)
               cur_tr%polarr%Gc=conjg(cur_tr%polarr%Gc)
            end if
            cur_tr%polarr%A=cur_tr%polarr%A*mult_quad*mult_el
            cur_tr%polarr%Ac=cur_tr%polarr%Ac*mult_quad*mult_el
            if(idx>0)then
               call AddTransition(trs(idx),cur_tr,.true.)
            else
               trs_c=trs_c+1
               ii=trs_c
               trs(ii)=cur_tr
            end if
         end if
      end do
      filepolars_c=filepolars_c+1
      !wexcs(i)=wexc
      close(77)
   end do
   
   write(output_unit,*)'FILE.Q.TTT count: ',fileqttt_c
   write(output_unit,*)'FILE...POLARS count: ',filepolars_c
   allocate(rroa_spr(npx,2,ns0),rroa_inv(ni0),rroa_tabOutputs(ns0),x(npx))
   do i = 1,ns0
      rroa_tabOutputs(i)=99+i!starts at 100, 1st wexc ends at 118
   end do
   rroa_spr=0d0
   rroa_inv=0d0
   do i = 1,npx
      x(i)=dble(wmax-wmin)/(npx-1)*(i-1)+wmin
   end do
   allocate(rroa_spr_temp(npx,2,ns0))
   rroa_spr_do=.false.
   rroa_spr_do(12)=.true. !SCP180
   if(dcp1)rroa_spr_do(15)=.true. !DCPI180
   if(dcp2)rroa_spr_do(18)=.true. !DCPII180
   
   if(wexc==0d0)then
      write(output_unit,*)'WARNING: "--wexc" is 0'
     ! call exit(123)
   end if
   write(wexc_str,'(F6.1,"nm")')dble(NINT(wexc*10))/10d0
   wexc_str=TR(wexc_str)
   
   do i = 1,ns0
      if(rroa_spr_do(i))then
         open(rroa_tabOutputs(i),file=TR(rroa_exp(i))//'_'//TR(wexc_str)//'.TAB')
      end if
   end do
   
   
   do i = 1,trs_c
      w1_au=trs(i)%w1*cm_2_au
      w3_au=trs(i)%w3*cm_2_au
      
      
      curPolar=trs(i)%polarr
      v1_class=count(trs(i)%v1>0)
      v3_class=count(trs(i)%v3>0)
      
      if(.not.IsFundamentalTransition(trs(i)%v1,v1_class,trs(i)%v3,v3_class).and..not.aboveFundamentals)cycle
      
      str_v1=FC2Str_new(trs(i)%v1,trs(i)%v1_pos,trs(i)%v1_n,.true.,.true.)
      str_v3=FC2Str_new(trs(i)%v3,trs(i)%v3_pos,trs(i)%v3_n,.true.,.true.)
      
      
      if(allocated(trs_anharm))then
         idx=FindTransition(trs(i),trs_anharm,trs_anharm_n,trs_anharm_n)
         if(idx>0)then
            w3_au=trs_anharm(idx)%w3*cm_2_au
         end if
         if(v1_class/=0)then
            tr_fake%v1=[0]
            tr_fake%v1_pos=[0]
            tr_fake%v1_n=0
            tr_fake%v3=trs(i)%v1
            tr_fake%v3_pos=trs(i)%v1_pos
            tr_fake%v3_n=trs(i)%v1_n
            idx=FindTransition(tr_fake,trs_anharm,trs_anharm_n,trs_anharm_n)
            if(idx>0)then
               w1_au=trs_anharm(idx)%w3*cm_2_au
            end if
         end if
      end if
      if(wexc/=0d0)then
         wexc_au=1d7/wexc*cm_2_au
      else
         wexc_au=0d0
      end if
      wr=wexc_au-w3_au
      kbt=exp(-h*cc*trs(i)%w1*100/(kb*temp))
      rroa_spr_temp=0d0
      
      
      
      if(.not.useG)curPolar%G=(0d0,0d0)
      if(.not.useGc)curPolar%Gc=(0d0,0d0)
      if(.not.useA)curPolar%A=(0d0,0d0)
      if(.not.useAc)curPolar%Ac=(0d0,0d0)
      
      call wrram3(wexc_au,temp,spectrum_temp,w1_au,w3_au,TINY(1d0),i,wexc,dble(wmin),dble(wmax),npx,use_gauss,FWHM,.true.,rroa_inv,&
                   realpart(curPolar%Ap),imagpart(curPolar%Ap), &
                   realpart(curPolar%G),imagpart(curPolar%G), &
                   realpart(curPolar%Gc),imagpart(curPolar%Gc), &
                   realpart(curPolar%A),imagpart(curPolar%A), &
                   realpart(curPolar%Ac),imagpart(curPolar%Ac), &
                   rroa_spr_temp,kBT,useA.or.useAc,useG.or.useGc,rroa_spr_do,rroa_tabOutputs,str_v1,str_v3, &
                   realpart(curPolar%Ap),imagpart(curPolar%Ap))
      rroa_spr(:,1,:)=rroa_spr(:,1,:)+rroa_spr_temp(:,1,:)
      rroa_spr(:,2,:)=rroa_spr(:,2,:)+rroa_spr_temp(:,2,:)
   end do
   
   do i = 1,ns0
      if(.not.rroa_spr_do(i))cycle
      write(c80,3000)TR(rroa_exp(i)),TR(wexc_str)
      write(c80_2,3001)TR(rroa_exp(i)),TR(wexc_str)
3000  format('RAM',A,'_',A,'.PRN')
3001  format('ROA',A,'_',A,'.PRN')
      !c80='RAM'//TR(rroa_exp(i))//'_'//TR(wexc_nm_str(j))//'.PRN'
      !c80_2='ROA'//TR(rroa_exp(i))//'_'//TR(wexc_nm_str(j))//'.PRN'
      open(77,file=c80)
      open(78,file=c80_2)
      do k = 1,npx
         w=dble(wmax-wmin)/npx*(k-1)+wmin
         write(77,'(f14.4,g25.12)')w,rroa_spr(k,1,i)
         write(78,'(f14.4,g25.12)')w,rroa_spr(k,2,i)
      end do
      close(77)
      write(output_unit,*)'Written file: '//c80
      close(78)
      write(output_unit,*)'Written file: '//c80_2
   end do
   
   do i = 1,ns0
      if(rroa_spr_do(i))then
         close(rroa_tabOutputs(i))
      end if
   end do
   
   contains
   
   function IsFundamentalTransition(v1,v1_class,v3,v3_class)result(res)
      integer v1_class,v3_class
      integer v1(v1_class),v3(v3_class)
      logical res
      
      res=.false.
      if(v1_class > 0 .or. v3_class > 1)return
      if(v3(1)>1)return
      res=.true.
   end function IsFundamentalTransition
   
   subroutine Readfileanharmfreq(unitt,trs,n_tr)
      integer unitt
      type(Transition),allocatable :: trs(:)
      integer :: n_tr,i,v_n
      
      n_tr=0
      do while(.true.)
         read(unitt,*,end=1000)
         n_tr=n_tr+1
      end do
1000  rewind(unitt)
      n_tr=n_tr/2
      allocate(trs(n_tr))
      do i = 1,n_tr
         read(unitt,*)v_n
         trs(i)%v1_n=0
         trs(i)%v3_n=v_n
         trs(i)%v1=[0]
         trs(i)%v1_pos=[0]
         allocate(trs(i)%v3(v_n),trs(i)%v3_pos(v_n))
         read(unitt,*)trs(i)%v3,trs(i)%v3_pos,trs(i)%w1,trs(i)%w3
      end do
   end subroutine Readfileanharmfreq
   
   function readsi_wrap(fn,z_at,E)result(S)
      character(*) fn
      double precision,allocatable :: S(:,:),E(:),r(:)
      integer,allocatable :: z_at(:)
      integer n3,nq,iz
      logical ldz
      
      ldz=.true.
      call readsi(n3,S,E,NQ,fn,z_at,r,ldz,iz)
      S=S*0.0234280d0 !this is conversion from amu^-1/2 to au^-1/2
      deallocate(r)
   end function readsi_wrap
   
   
   
   
   
   subroutine ReadArgs(optend,trs,fileqttt_c,polar_fileqttt_coeff,trs_c,finp,finp_qttt,switchPh,wexc)
      integer argc,i,idx_eq,optend,buf,fileqttt_c,trs_c,ii,jj,n3,ix
      double precision polar_fileqttt_coeff,vec3(3),wexc,wexc_tmp
      double precision,allocatable :: finp(:,:),finp_qttt(:,:),vec(:),finp_copy(:,:),m(:),m_qttt(:),E(:)
      double precision,allocatable :: E_finp(:),E_finp_qttt(:)
      double complex bufc
      double complex,allocatable :: vecC(:),pols(:)
      type(transition) cur_tr,trs(max_trs)
      type(polar),allocatable :: polars_fileqttt(:),polars_fileqttt_copy(:)
      type(polar) polar3(3)
      integer, parameter :: str_l=500
      character(str_l) s80,s80_2,ccc,s80_c
      character(1) ls
      integer perm(3),a,b,c,cur_nat,bufI,nq_fileqttt
      integer,allocatable :: perm_at(:),z_at(:),z_at_qttt(:),vecI(:),order(:)
      logical isFilettt,finp_found,finp_qttt_found
      logical,allocatable :: switchPh(:)
      
      isfilettt=.false.
      finp_found=.false.
      finp_qttt_found=.false.
      argc=COMMAND_ARGUMENT_COUNT()
      perm=[1,2,3]
      do i = 1,argc
         call GET_COMMAND_ARGUMENT(i,s80)
         s80_c=s80
         call To_lower(s80)
         if(s80(1:2)=='--')then
            idx_eq=index(s80,'=')+1
            
            if(idx_eq==1)then
               s80_2=s80(1:str_l)
            else
               s80_2=s80(1:idx_eq-1)
            end if
            select case(TR(s80_2))
               case('--abovefund')
                  aboveFundamentals=.true.
               case('--vib-min-cm=')
                  read(s80(idx_eq:),*)vib_min
               case('--vib-max-cm=')
                  read(s80(idx_eq:),*)vib_max
               case('--dcp1')
                  dcp1=.true.
               case('--dcp2')
                  dcp2=.true.
               case('--npx=')
                  read(s80(idx_eq:),*)npx
               case('--wmin=')
                  read(s80(idx_eq:),*)wmin
               case('--wmax=')
                  read(s80(idx_eq:),*)wmax
               case('--fwhm=')
                  read(s80(idx_eq:),*)fwhm
               case('--temp=')
                  read(s80(idx_eq:),*)temp
               case('--lineshape=')
                  ls=s80(idx_eq:idx_eq)
                  if(ls=='g')then
                     use_gauss=.true.
                  else if(ls=='l')then
                     use_gauss=.false.
                  else
                     write(output_unit,*)'Unknown lineshape: '//ls
                     call exit(23)
                  end if
               case('--finp-qttt=')
                  finp_qttt=readsi_wrap(s80_c(idx_eq:),z_at_qttt,E_finp_qttt)
                  finp_qttt_found=.true.
               case('--finp=')
                  finp=readsi_wrap(s80_c(idx_eq:),z_at,E_finp)
                  finp_found=.true.
               case('--fileqttt=')
                  isfilettt=index(s80_c(idx_eq:),'.Q.TTT')==0
                  if(.not.isfilettt)then
                     !stop 3
                     !open(77,file=s80_c(idx_eq:),status='old',err=555)
                     call Readfileqttt('FILE.Q.TTT',polars_fileqttt,.false.,E,nq_fileqttt,wexc_tmp)
                     n3=nq_fileqttt+6
                  else
                     call ReadFilettt(s80_c(idx_eq:),polars_fileqttt,n3,wexc_tmp)
                     nq_fileqttt=n3-6
                  end if
                  fileqttt_c=fileqttt_c+1
                  if(wexc==0d0)wexc=1d7/(wexc_tmp*au_2_cm)
               case('--fileqttt-modes=')
                  read(s80(idx_eq:),*)ccc
                  
                  allocate(modes_fqttt(array_temp_c))
                  modes_fqttt=0
                  read(ccc,*)modes_fqttt
                  buf=count(modes_fqttt/=0)
                  deallocate(modes_fqttt)
                  allocate(modes_fqttt(buf))
                  read(ccc,*)modes_fqttt
               case('--fileqttt-permute=')
                  read(s80(idx_eq:idx_eq),*)perm(1)
                  read(s80(idx_eq+1:idx_eq+1),*)perm(2)
                  read(s80(idx_eq+2:idx_eq+2),*)perm(3)
               case('--filettt-permute-atoms=')
                  open(777,file=s80(idx_eq:),status='old')
                  read(777,*)cur_nat
                  allocate(perm_at(cur_nat))
                  ii=1
                  do while(.true.)
                     read(777,*,end=573)bufI
                     perm_at(ii)=bufI
                     ii=ii+1
                  end do
                  
573               close(777)
                  call CheckPermAt(perm_at,cur_nat)
               case('--wexc=')
                  read(s80(idx_eq:),*)wexc
               case('--no-spectrum-temp')
                  spectrum_temp=.false.
               case('--vibronic-coeff=')
                  read(s80(idx_eq:),*)vibronic_coeff
               case('--fileqttt-coeff=')
                  read(s80(idx_eq:),*)polar_fileqttt_coeff
                  if(fileqttt_c>0)then
                     do ii = 1,nq_fileqttt
                        trs(ii)%polarr=polar_fileqttt_coeff*trs(ii)%polarr
                     end do
                  end if
               case('--polar-el-coeff=')
                  read(s80(idx_eq:),*)polars_el_coeff
               case('--fileqttt-wmin')
                  read(s80(idx_eq:),*)wmin_fqttt
               case('--fileqttt-wmax')
                  read(s80(idx_eq:),*)wmax_fqttt
               case('--mag-coeff=')
                  read(s80(idx_eq:),*)mult_mag
               case('--el-coeff-ttt=')
                  read(s80(idx_eq:),*)mult_el_fqttt
               case('--el-coeff=')
                  read(s80(idx_eq:),*)mult_el
               case('--quad-coeff=')
                  read(s80(idx_eq:),*)mult_quad
               case('--no-a')
                  useA=.false.
               case('--no-ac')
                  useAc=.false.
               case('--no-g')
                  useG=.false.
               case('--no-gc')
                  useGc=.false.
               case('--spr-coeff=')
                  read(s80(idx_eq:),*)spr_coeff
               case default
                  write(output_unit,*)'Unknown option: '//TR(s80)
            end select
         else
            optEnd=i-1
            exit
            !return
         end if
      end do
      if(i>argc)optEnd=argc
      
      if(.not.fileqttt_c>0)return
      if(fileqttt_c>0 .and. isfilettt .and. (.not. finp_found))then
         write(output_unit,*)'ERROR: With FILE.TTT you need F.INP'
         call exit(6)
      end if
      
      
      
      do i=1,nq_fileqttt
         do a = 1,3
            do b = 1,3
               polars_fileqttt(i)%Gc(a,b)=-polars_fileqttt(i)%G(b,a)
            end do
         end do
         polars_fileqttt(i)%Ac=polars_fileqttt(i)%A
      end do
      
      
      call FILEQTTT_2_Transitions(nq_fileqttt,E,polars_fileqttt,polar_fileqttt_coeff,trs,trs_c,max_trs)
      deallocate(polars_fileqttt)
      return
      
555   write(output_unit,*)'ERROR FILE ',TR(s80(idx_eq:)),'NOT FOUND!'
      call exit(5)
   end subroutine ReadArgs
   
   subroutine FixPhases(nq,n3,finp,finpqttt,switchPh)
      integer nq,n3
      double precision finp(n3,nq),finpqttt(n3,nq)
      logical switchPh(nq)
      integer i,j,k
      double precision buf1,buf2
      
      switchPh=.false.
      do i = 1,nq
         buf1=sum((finp(:,i)-finpqttt(:,i))**2)
         buf2=sum((finp(:,i)+finpqttt(:,i))**2)
         if(buf2<buf1)then
            switchPh(i)=.true.
         end if
      end do
      
   end subroutine FixPhases
   
   
   subroutine CheckFINP(nq,n3,finp,finp_q)
      integer nq,n3
      double precision finp(n3,nq),finp_q(n3,nq),vec(nq)
      integer i,j
      
      
      do i = 1,nq
         vec=0d0
         do j = 1,nq
            vec(j)=sum(abs(finp(:,i))*abs(finp_q(:,j)))
         end do
         j=0
      end do
   end subroutine CheckFINP
   
   subroutine FixPermutation(perm,order,n)
      integer n, perm(n),perm_c(n), order(n)
      integer i,buf,place
      
      do i = 1,n
         order(i)=i
      end do
      
      i=1
      perm_c=perm
      do while(i<=n)
         if(perm_c(i)==i)then
            i=i+1
            cycle
         end if
         place=perm_c(i)
         buf=perm_c(place)
         perm_c(place)=place
         perm_c(i)=buf
         
         buf=order(place)
         order(place)=order(i)
         order(i)=buf
      end do
      
   end subroutine FixPermutation
   
   subroutine CheckPermAt(perm,n)
      integer n,perm(n),i
      logical check(n)
      
      check=.false.
      do i = 1,n
         check(perm(i))=.not.check(perm(i))
      end do
      
      if(.not.all(check))then
         write(output_unit,*)'ERROR: Permutation file does not contain correct atom numbers.'
         call exit(25)
      end if
   end subroutine CheckPermAt
   
   subroutine FILEQTTT_2_Transitions(nq_fileqttt,fileqttt_e,polars,polar_fileqttt_coeff,trs,trs_c,trs_max)
      integer trs_c,nq_fileqttt,trs_max,i
      double precision :: fileqttt_e(nq_fileqttt),polar_fileqttt_coeff
      type(Polar) polars(nq_fileqttt)
      type(transition) trs(trs_max),cur_tr
      
      do i = 1,nq_fileqttt
         cur_tr%w1=0
         cur_tr%w3=fileqttt_e(i)
         cur_tr%v1_n=0
         cur_tr%v3_n=1
         cur_tr%v1=[0]
         cur_tr%v3=[1]
         cur_tr%v1_pos=[0]
         cur_tr%v3_pos=[i]
         cur_tr%polarr=1d0/sqrt(2d0*cur_tr%w3*cm_2_au)*polar_fileqttt_coeff*polars(i)
         trs(i)=cur_tr
      end do
      trs_c=nq_fileqttt
   end subroutine FILEQTTT_2_Transitions
   
   function TrPol(iq,nq,n3,finp,pols)
      integer nq,n3,iq,i,j
      double precision finp(n3,nq)
      double complex TrPol,pols(n3)
      
      TrPol=(0d0,0d0)
      do i = 1,n3
         TrPol=TrPol+finp(i,iq)*pols(i)
      end do
      
   end function TrPol
   
   
   subroutine AddTransition(tr,tr2,check)
      type(transition) tr,tr2
      logical check
      
      if(check)then
         if(.not.TransitionsAreEqual(tr,tr2))then
            write(output_unit,*)'Transitions are not equal'
            stop 33
         end if
      end if
      
      tr%polarr=tr%polarr+tr2%polarr
   end subroutine AddTransition
   
   function FindTransition(tr,trs,trs_n,trs_c)result(idx)
      integer trs_n
      type(transition) tr,trs(trs_n)
      integer i,idx,trs_c
      
      idx=0
      do i = 1,trs_c
         if(TransitionsAreEqual(tr,trs(i)))then
            idx=i
            return
         end if
      end do
   end function FindTransition
   
   function TransitionsAreEqual(tr1,tr2)result(res)
      type(transition) tr1,tr2
      integer n1,n3
      logical res
      
      res=.false.
      if(tr1%v1_n/=tr2%v1_n .or. tr1%v3_n/=tr2%v3_n)return
      n1=tr1%v1_n
      n3=tr1%v3_n
      if(Equal_I2Arr(tr1%v1_pos,tr2%v1_pos,n1) .and. Equal_I2Arr(tr1%v3_pos,tr2%v3_pos,n3))then
         if(Equal_IArr(tr1%v1,tr2%v1) .and. Equal_IArr(tr1%v3,tr2%v3))then
            res=.true.
         end if
      end if
   end function TransitionsAreEqual
   
   
   subroutine Permute2D(arr,perm)
      double complex arr(3,3),arr_buf(3,3)
      integer perm(3),i,j
      
      arr_buf=arr
      
      do i = 1,3
         do j = 1,3
            arr(i,j)=arr_buf(perm(i),perm(j))
         end do
      end do
   end subroutine Permute2D
   
   subroutine Permute3D(arr,perm)
      double complex arr(3,3,3),arr_buf(3,3,3)
      integer perm(3),i,j,k
      
      arr_buf=arr
      
      do i = 1,3
         do j = 1,3
            do k = 1,3
               arr(i,j,k)=arr_buf(perm(i),perm(j),perm(k))
            end do
         end do
      end do
   end subroutine Permute3D
   
   
   
   subroutine PolarZero(polarr)
      type(Polar) :: polarr
      polarr%ap=0d0
      polarr%G=0d0
      polarr%Gc=0d0
      polarr%A=0d0
      polarr%Ac=0d0
   end subroutine PolarZero
   
   
   
   
end program FCOV_spectrum