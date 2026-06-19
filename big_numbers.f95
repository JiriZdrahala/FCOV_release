module ownmath
   use iso_fortran_env
   implicit none
   
   double complex, parameter :: iu = (0d0,1d0)
   double precision,parameter :: pi = 4*atan(1d0)
   
   private iu,pi
   
   type,public :: big_double
      integer p
      double precision num
   end type big_double
   
   type,public :: big_complex
      type(big_double) r
      double precision phi
   end type big_complex
   
   ! interface operator(+)
      ! module procedure bc_add_complex,complex_add_bc
   ! end interface
   
   interface operator(*)
      module procedure bd_mult
      module procedure bd_mult_double,double_mult_bd
      module procedure bd_mult_int,int_mult_bd
      
      module procedure bc_mult
      module procedure bc_mult_bd,bd_mult_bc
      module procedure bc_mult_complex,complex_mult_bc
      module procedure bc_mult_double,double_mult_bc
      module procedure bc_mult_int,int_mult_bc
      
   end interface
   
   interface operator(/)
      module procedure bd_div,bd_div_double,double_div_bd
      module procedure bc_div
      module procedure bc_div_complex
   end interface
   
   interface BD
      module procedure BD_direct,BD_derive
   end interface BD
   
   interface BC
      module procedure BC_derive
   end interface BC
   
   interface co_sqrt
      module procedure co_sqrt_cmplx,co_sqrt_big
   end interface co_sqrt
   
   contains
   
   
   ! pure function floor_own(a)result(res)
      ! double precision,intent(in) :: a
      ! integer res
      ! if(a<0d0)then
         ! res=ceiling(a)
      ! else
         ! res=floor(a)
      ! end if
   ! end function floor_own
   
   pure function BC_derive(z)result(res)
      double complex,intent(in) :: z
      double precision sz,x,y
      type(big_complex) res
      
      x=realpart(z)
      y=imagpart(z)
      sz=sqrt(x**2+y**2)
      res%r=BD_derive(sz)
      res%phi=atan2(y,x)
   end function BC_derive
   
   pure function BD_direct(num,p)result(res)
      integer,intent(in) :: p
      double precision,intent(in) :: num
      type(big_double) res
      
      res%num=num
      res%p=p
      call BD_Normalize(res)
   end function BD_direct
   
   pure function BD_derive(a)result(res)
      double precision,intent(in) :: a
      type(big_double) res
      
      res%p=floor(log10(a))
      res%num=a/(10d0**res%p)
   end function BD_derive
   
   pure subroutine BD_Normalize(bd)
      type(big_double),intent(inout) :: bd
      double precision l
      
      l=log10(bd%num)
      if(l>=1d0 .or. l<=-1d0)then
         bd%p=bd%p+floor(l)
         bd%num=bd%num*(10d0**(-floor(l)))
      end if
      
      ! if(abs(bd%num)>1d1)then
         ! bd%p=bd%p+1
         ! bd%num=bd%num*1d-1
      ! else if(abs(bd%num)<1d-1)then
         ! bd%p=bd%p-1
         ! bd%num=bd%num*1d1
      ! end if
   end subroutine BD_Normalize
   
   pure subroutine BC_Normalize(bc)
      type(big_complex),intent(inout) :: bc
      call BD_Normalize(bc%r)
      ! if(bc%r%num<0d0)stop 567
      return
      
      bc%phi=mod(bc%phi,2*pi)
      if(bc%phi>pi)then
         bc%phi=bc%phi-2*pi
      elseif(bc%phi<=-pi)then
         bc%phi=bc%phi+2*pi
      end if
   end subroutine BC_Normalize
   
   pure function BD_to_BC(bd)result(res)
      type(big_double),intent(in) :: bd
      type(big_complex) res
      
      res%r%p=bd%p
      res%r%num=abs(bd%num)
      if(bd%num<0d0)then
         res%phi=pi
      else
         res%phi=0d0
      end if
      
   end function BD_to_BC
   
   pure function bd_mult(l,r)result(res)
      type(big_double),intent(in) :: l,r
      type(big_double) res
      
      res%p=l%p+r%p
      res%num=l%num*r%num
      call BD_Normalize(res)
   end function bd_mult
   
   pure function bd_mult_double(l,r)result(res)
      type(big_double),intent(in) :: l
      double precision,intent(in) :: r
      type(big_double) res
      
      res%num=l%num*r
      res%p=l%p
      call BD_Normalize(res)
   end function bd_mult_double
   
   pure function double_mult_bd(l,r)result(res)
      double precision,intent(in) :: l
      type(big_double),intent(in) :: r
      type(big_double) res
      
      res=bd_mult_double(r,l)
   end function double_mult_bd
   
   pure function bd_mult_int(l,r)result(res)
      type(big_double),intent(in) :: l
      integer,intent(in) :: r
      type(big_double) res
      
      res%num=l%num*r
      res%p=l%p
      call BD_Normalize(res)
   end function bd_mult_int
   
   pure function int_mult_bd(l,r)result(res)
      integer,intent(in) :: l
      type(big_double),intent(in) :: r
      type(big_double) res
      
      res=bd_mult_int(r,l)
   end function int_mult_bd
   
   pure function bc_mult(l,r)result(res)
      type(big_complex),intent(in) :: l,r
      type(big_complex) res
      
      res%r=l%r*r%r
      res%phi=l%phi+r%phi
      call BC_Normalize(res)
   end function bc_mult
   
   pure function bc_mult_complex(l,r)result(res)
      type(big_complex),intent(in) :: l
      double complex,intent(in) :: r
      type(big_complex) res
      double precision sz,phi
      
      call co_polar(r,sz,phi)
      
      res%r=l%r*sz
      res%phi=l%phi+phi
      call BC_Normalize(res)
   end function bc_mult_complex
   
   pure function bc_div_complex(l,r)result(res)
      type(big_complex),intent(in) :: l
      double complex,intent(in) :: r
      type(big_complex) res
      double precision sz,phi
      
      call co_polar(r,sz,phi)
      res%r=l%r/sz
      res%phi=l%phi-phi
      call BC_Normalize(res)
   end function bc_div_complex
   
   pure function complex_mult_bc(l,r)result(res)
      double complex,intent(in) :: l
      type(big_complex),intent(in) :: r
      type(big_complex) res
      
      res=bc_mult_complex(r,l)
   end function complex_mult_bc
   
   pure function bc_mult_double(l,r)result(res)
      type(big_complex),intent(in) :: l
      double precision,intent(in) :: r
      double complex :: r_c
      type(big_complex) res
      
      r_c=dble(r)
      res=l*r_c
      ! res%r=l%r*r
      ! res%phi=l%phi
      call BC_Normalize(res)
   end function bc_mult_double
   
   pure function double_mult_bc(l,r)result(res)
      double precision,intent(in) :: l
      type(big_complex),intent(in) :: r
      type(big_complex) res
      
      res=bc_mult_double(r,l)
   end function double_mult_bc
   
   pure function bc_mult_int(l,r)result(res)
      type(big_complex),intent(in) :: l
      integer,intent(in) :: r
      double complex r_c
      type(big_complex) res
      
      
      r_c=dble(r)
      res=l*r_c
      ! res%r=l%r
      ! res%phi=l%phi
      call BC_Normalize(res)
   end function bc_mult_int
   
   pure function int_mult_bc(l,r)result(res)
      integer,intent(in) :: l
      type(big_complex),intent(in) :: r
      type(big_complex) res
      
      res=bc_mult_int(r,l)
   end function int_mult_bc
   
   pure function bc_div(l,r)result(res)
      type(big_complex),intent(in) :: l,r
      type(big_complex) res
      
      res%r=l%r/r%r
      res%phi=l%phi-r%phi
      call BC_Normalize(res)
   end function bc_div
   
   pure function bc_mult_bd(l,r)result(res)
      type(big_complex),intent(in) :: l
      type(big_double),intent(in) :: r
      type(big_complex) res,r_c
      
      r_c=BD_to_BC(r)
      res=l*r_c
      call BC_Normalize(res)
   end function bc_mult_bd
   
   pure function bd_mult_bc(l,r)result(res)
      type(big_double),intent(in) :: l
      type(big_complex),intent(in) :: r
      type(big_complex) res
      
      res=bc_mult_bd(r,l)
   end function bd_mult_bc
   
   pure function bd_div(l,r)result(res)
      type(big_double),intent(in) :: l,r
      type(big_double) res
      
      res%p=l%p-r%p
      res%num=l%num/r%num
      call BD_Normalize(res)
   end function bd_div
   
   pure function bd_div_double(l,r)result(res)
      type(big_double),intent(in) :: l
      double precision,intent(in) :: r
      type(big_double) res
      
      res%num=l%num/r
      res%p=l%p
      call BD_Normalize(res)
   end function bd_div_double
   
   pure function double_div_bd(l,r)result(res)
      double precision,intent(in) :: l
      type(big_double),intent(in) :: r
      type(big_double) res
      
      res=bd_div_double(r,l)
   end function double_div_bd
   
   pure function co_size(z)result(res)
      double complex,intent(in) :: z
      double precision res,x,y
      
      x=realpart(z)
      y=imagpart(z)
      res=sqrt(x**2+y**2)
   end function co_size
   
   pure function co_sqrt_cmplx(z,principal)result(z_sqrt)
      logical,intent(in) :: principal
      double complex,intent(in) :: z
      
      double complex z_sqrt
      double precision x,y,z_sz
      
      x=realpart(z)
      y=imagpart(z)
      z_sz=sqrt(x**2+y**2)
      if(principal)then
         z_sqrt=sqrt(0.5d0*(z_sz+x))+iu*sgn(y)*sqrt(0.5d0*(z_sz-x))
      else
         z_sqrt=-(sqrt(0.5d0*(z_sz+x))+iu*sgn(y)*sqrt(0.5d0*(z_sz-x)))
      end if
      
   end function co_sqrt_cmplx
   
   pure function bd_sqrt(a)result(res)
      type(big_double),intent(in) :: a
      type(big_double) res
      double precision pw,num
      integer pw_i
      
      pw=dble(a%p)/2d0
      if(pw<0d0)then
         pw_i=ceiling(pw)
      else
         pw_i=floor(pw)
      end if
      pw=pw-pw_i
      res%num=sqrt(a%num)*10d0**pw
      res%p=pw_i
      !res=BD(num*10d0**pw)
      call BD_Normalize(res)
   end function bd_sqrt
   
   pure function co_sqrt_big(z,principal)result(z_sqrt)
      logical,intent(in) :: principal
      type(big_complex),intent(in) :: z
      
      type(big_complex) z_sqrt,z_intm
      
      z_intm=z
      call BC_Normalize(z_intm)
      if(principal)then
         z_sqrt%r=bd_sqrt(z_intm%r)
         z_sqrt%phi=z_intm%phi/2d0
      else
         z_sqrt%r=bd_sqrt(z_intm%r)
         z_sqrt%phi=z_intm%phi/2d0
      end if
      call BC_Normalize(z_sqrt)
   end function co_sqrt_big
   
   pure function bc_extract(bc)result(z)
      type(big_complex),intent(in) :: bc
      double complex z
      
      z=bd_extract(bc%r)*(cos(bc%phi)+iu*sin(bc%phi))
   end function bc_extract
   
   pure function bd_extract(bd)result(a)
      type(big_double),intent(in) :: bd
      double precision a
      
      a=bd%num*(10d0**bd%p)
   end function bd_extract
   
   pure subroutine co_conjg(z)
      type(big_complex),intent(inout) :: z
      z%phi=-z%phi
   end subroutine co_conjg
   
   pure subroutine co_polar(z,r,phi)
      double complex,intent(in) :: z
      double precision, intent(out) :: r,phi
      
      double precision x,y
      
      x=realpart(z)
      y=imagpart(z)
      r=sqrt(x**2+y**2)
      phi=atan2(y,x)
   end subroutine co_polar
   
   pure function sgn(a)
      double precision,intent(in) :: a
      double precision sgn
      
      if(a>0)then
         sgn=1
      else if(a<0)then
         sgn=-1
      else
         sgn=0
      end if
   end function sgn
   
end module ownmath
