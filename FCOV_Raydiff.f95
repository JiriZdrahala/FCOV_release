#define TR(arg) trim(adjustl(arg))
program FCOV_Raydiff
   use iso_fortran_env
   use util
   use constants
   implicit none
   character(300) fn,dir
   character(80) s80
   character(10) wexc_str
   integer n,n3,nq,nat,steps
   integer a,b,c,i
   double precision wexc,e00,step
   double precision,allocatable :: smat(:,:),wg(:),r(:)
   !double precision,allocatable :: smats(:,:,:),wgs(:,:),rs(:,:)
   integer,allocatable :: z_at(:)
   logical isnm
   type(polar),allocatable :: curpolar(:),polars(:,:),polars_R(:),polars_Q(:)
   
   call Get_Command_Argument(1,s80)
   read(s80,*)step
   step=step/BohrR
   call readsi(n3,smat,wg,nq,'F.INP',z_at,r,.true.,n)
   
   steps=2*n3
   
   !allocate(smats(n3,nq,steps),wgs(nq,steps),rs(n3,steps))
   allocate(polars(1,steps),polars_R(n3),polars_Q(nq))
   
   open(77,file='LIST.STEPS.RRay',status='old')
   i=0
90 read(77,'(A300)',end=99)fn
   i=i+1
   !dir=dirname(fn)
   curpolar=ReadPolars(666,fn,n,wexc,e00,isnm)
   polars(1,i)=curpolar(1)
   deallocate(curpolar) !slow in a loop, todo fix
   goto 90
   
99 close(77)
   if(i/=steps)then
      write(output_unit,*)'No'
      call exit(2)
   end if
   
   do i = 1,n3
      do a = 1,3
         do b = 1,3
            polars_R(i)%ap(a,b)=diffC(polars(1,i)%ap(a,b),polars(1,i+n3)%ap(a,b),step)
            polars_R(i)%G(a,b)=diffC(polars(1,i)%G(a,b),polars(1,i+n3)%G(a,b),step)
            polars_R(i)%Gc(a,b)=diffC(polars(1,i)%Gc(a,b),polars(1,i+n3)%Gc(a,b),step)
            do c = 1,3
               polars_R(i)%A(a,b,c)=diffC(polars(1,i)%A(a,b,c),polars(1,i+n3)%A(a,b,c),step)
               polars_R(i)%Ac(a,b,c)=diffC(polars(1,i)%Ac(a,b,c),polars(1,i+n3)%Ac(a,b,c),step)
            end do
         end do
      end do
   end do
   
   do i = 1,nq
      polars_q(i)=Car2NM_Pol1(n3,nq,i,polars_R,smat)*(1d0/sqrt(2d0*wg(i)*cm_2_au))
   end do
   
   call WritePolars(.true.,nq,0d0,55,wexc,polars_q,wg*cm_2_au,.true.,.true.)
   write(wexc_str,'(F6.1,"nm")')dble(NINT(1d7/(wexc*au_2_cm)*10))/10d0
   wexc_str=TR(wexc_str)
   write(output_unit,'(A,A)')'Written transition polarizabilities for ',wexc_str
   
   contains
   
   function diffC(fxhp,fxhm,h)result(res)
      double complex fxhp,fxhm,res
      double precision h
      res=(fxhp-fxhm)/(2*h)
   end function diffC
   
   function dirname(fn)
      character(*) fn
      character(len(fn)) dirname
      integer i
      do i = len(fn), 1, -1
         if (fn(i:i) == '/' .or. fn(i:i) == '\') then
            dirname=fn(1:i)
            return
         end if
      end do
   end function dirname
   
end program FCOV_Raydiff