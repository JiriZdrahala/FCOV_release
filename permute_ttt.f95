program permute_ttt
   use iso_fortran_env
   use util
   implicit none
   character(500) s500
   character(3) permCart_c,invert_c
   integer permCart(3),perm_n
   integer,allocatable :: perm(:)
   logical invert(3)
   type(Polar),allocatable :: polars(:)
   integer n,i
   double precision wexc
   
   call GET_COMMAND_ARGUMENT(1,s500)
   call GET_COMMAND_ARGUMENT(2,permCart_c)
   do i =1,3
      read(permCart_c(i:i),*)permCart(i)
   end do
   
   call GET_COMMAND_ARGUMENT(3,invert_c)
   do i=1,3
      if(invert_c(i:i)=='+')then
         invert(i)=.false.
      elseif(invert_c(i:i)=='-')then
         invert(i)=.true.
      else
         stop 11
      end if
   end do
   
   
   call ReadFilettt(s500,polars,n,wexc)
   call ReadPerm(perm,n)
   !call PermuteAtoms(polars,n,perm)
   call PermuteCartesianAxes(polars,n,permCart,invert)
   call WriteFilettt('FILE.TTT',polars,n,n/3,wexc)
   
   contains
   
   subroutine ReadPerm(perm,n)
      integer n,i
      integer,allocatable :: perm(:)
      logical ex
      
      
      inquire(file='PERMUTATION.TXT',exist=ex)
      if(.not.ex)then
         perm=seq(n,.false.)
         return
      end if
      
      
      open(77,file='PERMUTATION.TXT')
      allocate(perm(n))
      do i = 1,n
         read(77,*)perm(i)
      end do
      close(77)
   end subroutine ReadPerm
   
   subroutine PermuteAtoms(polars,n,perm)
      integer n,perm(n),summ,i,ia,ib,idx
      type(Polar) :: polars(n)
      type(Polar),allocatable :: polars_buf(:)
      
      summ=0
      do i =1,n
         summ=summ+i
      end do
      if(sum(perm)/=summ)then
         stop 2
      end if
      allocate(polars_buf(n))
      do i = 1,n
         polars_buf(i)=polars(i)
      end do
      
      do i = 1,n/3
         ia=(perm(i)-1)*3+1
         ib=ia+2
         idx=(i-1)*3+1
         polars(idx:idx+2)=polars_buf(ia:ib)
      end do
      deallocate(polars_buf)
   end subroutine PermuteAtoms
   
   subroutine PermuteCartesianAxes(polars,n,perm,invert)
      integer n,perm(3),summ,i,a,b,c,idx
      type(Polar) :: polars(n)
      type(Polar),allocatable :: polars_buf(:)
      logical invert(3)
      
      if(sum(perm)/=6 .or. product(perm)/=6)then
         stop 3
      end if
      
      allocate(polars_buf(n))
      do i = 1,n
         polars_buf(i)=polars(i)
      end do
      
      !permute derivative coordinates
      do i = 1,n/3
         idx=(i-1)*3
         polars(idx+1)=polars_buf(idx+perm(1))   !dx_i
         polars(idx+2)=polars_buf(idx+perm(2)) !dy_i
         polars(idx+3)=polars_buf(idx+perm(3)) !dz_i
      end do
      
      !permute polarizability coordinates
      ! do i = 1,n
         ! do a = 1,3
            ! do b = 1,3
               ! polars(i)%ap(a,b)=polars_buf(i)%ap(perm(a),perm(b))
               ! polars(i)%G(a,b)=polars_buf(i)%G(perm(a),perm(b))
               ! polars(i)%Gc(a,b)=polars_buf(i)%Gc(perm(a),perm(b))
               ! do c = 1,3
                  ! polars(i)%a(a,b,c)=polars_buf(i)%a(perm(a),perm(b),perm(c))
                  ! polars(i)%ac(a,b,c)=polars_buf(i)%ac(perm(a),perm(b),perm(c))
               ! end do
            ! end do
         ! end do
      ! end do
      
      
      
      deallocate(polars_buf)
      
   end subroutine PermuteCartesianAxes
   
end program permute_ttt