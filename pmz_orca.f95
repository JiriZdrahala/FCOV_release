#define TR(arg) trim(adjustl(arg))

program pmz_orca
   use iso_fortran_env
   use strings
   implicit none
   
   integer idiff,cc
   logical no_zero,q_dimless,only_double,lvstep
   logical :: isOrca=.false.
   double precision step

   integer nat,n3,nq,ch,m,inat,i,ic,ix
   integer,allocatable :: z(:)
   double precision, allocatable :: r(:,:),r_cur(:,:)
   character(80) filename
   character(1) signn,coord,axis
   character(4) num
   
   call ReadPar()
   call ReadCM(77,ch,m)
   call ReadGeom(r,z,nat)
   n3=3*nat
   nq=n3-6
   cc=0
   select case(idiff)
      case(2)
         write(output_unit,'(A)')'Central finite difference, 2 displacement types: dx, -dx'
         write(output_unit,'(A,I8,A)')'Will write ',n3*2,' files'
         call Diff_Central(0,nat,r,z,ch,m,step,isOrca,1,cc)
      case(22)
         write(output_unit,'(A)')'Central finite difference, 4 displacement types: 2*dx, dx, -dx, -2*dx'
         write(output_unit,'(A,I8,A)')'Will write ',n3*2*2,' files'
         call Diff_Central(0,nat,r,z,ch,m,step,isOrca,1,cc)
         call Diff_Central(0,nat,r,z,ch,m,step,isOrca,2,cc)
      case(4)
         write(output_unit,'(A)')'Central finite difference, 6 displacement types: dx, dy, dxdy, -dx, -dy, -dx-dy,'
         write(output_unit,'(A,I8,A)')'Will write ',n3*2+n3*(n3-1),' files'
         call Diff_Central(0,nat,r,z,ch,m,step,isOrca,1,cc)
         ic=2*3*nat
         call Diff_CentralMixed(ic,nat,r,z,ch,m,step,isOrca)
      case default
         stop 2
   end select
   
   contains
   
   subroutine Diff_CentralMixed(startIndex,nat,r,z,ch,m,step,isOrca)
      character(1) signn,coord,axis,axisj
      character(4) num,numj
      character(80) filename
      double precision r(3,nat),step
      double precision, allocatable :: r_cur(:,:)
      integer startIndex,z(nat),idx,missI
      logical isOrca
      integer nat,ch,m,inat,i,j,ix,jnat,jx
      
      coord='r'
      signn='+'
      allocate(r_cur(3,nat))
      r_cur=r
      do i = 1,3*nat
         inat=(i-1)/3+1
         ix=mod(i-1,3)+1
         write(num,'(I4)')inat
         call GetAxis(ix,axis)
         do j = i+1,3*nat
            jnat=(j-1)/3+1
            jx=mod(j-1,3)+1
            write(numj,'(I4)')jnat
            call GetAxis(jx,axisj)
            missI=i*(i+1)/2
            idx=(i-1)*3*nat+j-missI+startIndex
            call CreateFilename(filename,idx)
            open(77,file=filename)
            call WriteInputHead(77,78)
            write(77,*)

            !r_cur=r
            r_cur(ix,inat)=r_cur(ix,inat)+step
            r_cur(jx,jnat)=r_cur(jx,jnat)+step
            call WriteComment(77,coord//TR(num)//axis//signn//' '//coord//TR(numj)//axisj//signn,isOrca)
            write(77,*)
            call WriteGeomHeader(77,ch,m,isOrca)
            call WriteGeom(77,z,r_cur,nat,isOrca)
            close(77)
            r_cur(ix,inat)=r(ix,inat)
            r_cur(jx,jnat)=r(jx,jnat)
            if(mod(idx,100)==0)write(output_unit,*)filename
         end do
      end do
      signn='-'
      do i = 1,3*nat
         inat=(i-1)/3+1
         ix=mod(i-1,3)+1
         write(num,'(I4)')inat
         call GetAxis(ix,axis)
         do j = i+1,3*nat
            jnat=(j-1)/3+1
            jx=mod(j-1,3)+1
            write(numj,'(I4)')jnat
            call GetAxis(jx,axisj)
            missI=i*(i+1)/2
            idx=(i-1)*3*nat+j-missI+startIndex+n3*(n3-1)/2
            call CreateFilename(filename,idx)
            open(77,file=filename)
            call WriteInputHead(77,78)
            write(77,*)
            
            !r_cur=r
            r_cur(ix,inat)=r_cur(ix,inat)-step
            r_cur(jx,jnat)=r_cur(jx,jnat)-step
            call WriteComment(77,coord//TR(num)//axis//signn//' '//coord//TR(numj)//axisj//signn,isOrca)
            write(77,*)
            call WriteGeomHeader(77,ch,m,isOrca)
            call WriteGeom(77,z,r_cur,nat,isOrca)
            close(77)
            r_cur(ix,inat)=r(ix,inat)
            r_cur(jx,jnat)=r(jx,jnat)
            if(mod(idx,100)==0)write(output_unit,*)filename
         end do
      end do
      
      deallocate(r_cur)
   end subroutine Diff_CentralMixed
   
   subroutine Diff_Central(startIndex,nat,r,z,ch,m,step,isOrca,mult,counter)
      character(1) signn,coord,axis
      character(3) mult_str
      character(4) num
      character(80) filename
      logical isOrca
      double precision r(3,nat)
      double precision,value :: step
      double precision, allocatable :: r_cur(:,:)
      integer :: counter
      integer startIndex,z(nat),mult
      integer nat,ch,m,inat,i,ix
      
      step=step*mult
      coord='r'
      signn='+'
      mult_str='   '
      if(mult>1)mult_str=' *'//I2Str(mult,1)
      do i = 1,3*nat
         inat=(i-1)/3+1
         ix=mod(i-1,3)+1
         write(num,'(I4)')inat
         call GetAxis(ix,axis)
         call CreateFilename(filename,i+startIndex+counter)
         open(77,file=filename)
         call WriteInputHead(77,78)
         write(77,*)
         r_cur=r
         
         !f(x+h)
         r_cur(ix,inat)=r_cur(ix,inat)+step
         
         call WriteComment(77,coord//TR(num)//axis//signn//mult_str,isOrca)
         write(77,*)
         call WriteGeomHeader(77,ch,m,isOrca)
         call WriteGeom(77,z,r_cur,nat,isOrca)
         close(77)
         write(output_unit,*)filename
      end do
      
      signn='-'
      do i = 1,3*nat
         inat=(i-1)/3+1
         ix=mod(i-1,3)+1
         write(num,'(I4)')inat
         call GetAxis(ix,axis)
         call CreateFilename(filename,i+3*nat+startIndex+counter)
         open(77,file=filename)
         call WriteInputHead(77,78)
         write(77,*)
         r_cur=r
         
         !f(x-h)
         r_cur(ix,inat)=r_cur(ix,inat)-step
         
         call WriteComment(77,coord//TR(num)//axis//signn//mult_str,isOrca)
         write(77,*)
         call WriteGeomHeader(77,ch,m,isOrca)
         call WriteGeom(77,z,r_cur,nat,isOrca)
         close(77)
         write(output_unit,*)filename
      end do
      counter=counter+2*3*nat
      deallocate(r_cur)
   end subroutine Diff_Central
   
   subroutine GetAxis(i,str)
      integer i
      character(1) str
      
      select case(i)
         case(1)
            str='x'
         case(2)
            str='y'
         case(3)
            str='z'
         case default
         
         stop 66
      end select
   end subroutine GetAxis
   
   subroutine WriteComment(unitt,str,isOrca)
      integer unitt
      logical isorca
      character(*) str
      
      if(isOrca)then
         write(unitt,'(A,A)')'# ',TR(str)
      else
         write(unitt,'(A)')TR(str)
      end if
   end subroutine WriteComment
   
   subroutine WriteGeomHeader(unitt,ch,m,isorca)
      integer unitt,ch,m
      logical isOrca
      
      if(isOrca)then
         write(unitt,'(A,I2,I2)')'* xyz ',ch,m
      else
         write(unitt,'(1X,I2,I2)')ch,m
      end if
   end subroutine WriteGeomHeader
   
   subroutine WriteGeom(unitt,z,r,nat,isOrca)
      logical isOrca
      integer nat,unitt,z(nat),i
      double precision r(3,nat)
      
      do i = 1,nat
         write(unitt,'(1X,I3,3(1X,F15.8))')z(i),r(:,i)
      end do
      
      if(isOrca)write(unitt,'(A)')'*'
      write(unitt,*)
   end subroutine WriteGeom
   
   subroutine CreateFilename(filename,i)
      integer i
      character(80) filename
      
      if(i<10)then
         write(filename,'(A,I1,A)')'FILE_',i,'.inp'
      elseif(i<100)then
         write(filename,'(A,I2,A)')'FILE_',i,'.inp'
      elseif(i<1000)then
         write(filename,'(A,I3,A)')'FILE_',i,'.inp'
      elseif(i<10000)then
         write(filename,'(A,I4,A)')'FILE_',i,'.inp'
      elseif(i<100000)then
         write(filename,'(A,I5,A)')'FILE_',i,'.inp'
      elseif(i<1000000)then
         write(filename,'(A,I6,A)')'FILE_',i,'.inp'
      else
         stop 3
      end if
   end subroutine CreateFilename
   
   subroutine WriteInputHead(unitt,unittHead)
      integer unitt,unittHead
      character(500) s500
      
      open(unittHead,file='G.TXT')
      
      do
         read(unittHead,'(A)',end=20)s500
         write(unitt,'(A)')trim(s500)
      end do
      
20    close(unittHead)
   end subroutine WriteInputHead
   
   subroutine ReadCM(unitt,ch,m)
      integer unitt,ch,m
      
      open(unitt,file='CM.TXT')
      read(unitt,*)ch,m
      close(unitt)
   end subroutine ReadCM
   
   subroutine ReadGeom(r,z,nat)
      integer nat
      integer,allocatable :: z(:)
      double precision,allocatable :: r(:,:)
      
      open(77,file='FILE.X',status='old')
      
      read(77,*)
      read(77,*)nat
      
      allocate(r(3,nat),z(nat))
      do i = 1,nat
         read(77,*)z(i),r(1,i),r(2,i),r(3,i)
      end do
      
      close(77)
   end subroutine ReadGeom
   
   subroutine ReadPar()
      logical fex
      character(4) key
      
      inquire(file='PMZ.PAR',exist=fex)
      if(.not.fex)stop 1
      
      open(7,file='PMZ.PAR')

10    read(7,*,end=20)key      
      
      if(key.eq.'NO00')read(7,*)no_zero
      if(key.eq.'ORCA')isOrca=.true.
      if(key.eq.'STEP')read(7,*)step
      if(key.eq.'QDIM')read(7,*)q_dimless
      if(key.eq.'DONL')read(7,*)only_double !Do only double steps d/dq_i^2, no cross terms like d/dq_i.dq_j
      if(key.eq.'IDIF')read(7,*)idiff
      if(key(1:2).eq.'IC')read(7,*)ic
      goto 10
      
20    close(7)
      
   end subroutine ReadPar
   
end program pmz_orca