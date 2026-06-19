module strings
#define TR(arg) (trim(adjustl(arg)))

   implicit none
   !tried making my own String type, does not work well with gfortran (parametrized allocatable types are no bueno)
   
   integer,private :: ColumnWrite_curColumn=1
   
   contains
   
   
   function Ciphers(num)result(cif)
      integer num,cif
      if(num==0)then
         cif=1
         return
      end if
      cif=floor(log10(real(abs(num),8)))+1
   end function Ciphers

   function ShaveFront(str,charr)result(newStr)
      character(*) str
      character(*) charr
      
      integer idx
      character(200)newstr
      
      idx=index(str,charr,back=.true.)
      newStr=TR(str(idx+1:))
   end function ShaveFront
   
   subroutine ReplaceSuffix(filename,newSuffix,suffix)
      character(*) filename,newSuffix,suffix
      integer period,space,leng
      
      period=index(filename,suffix,.true.)
      leng=len(newSuffix)
      if(period==0)then
         space=index(filename,' ',.true.)
         filename(space:space+leng+1)=EmptyString(leng)
         filename(space:space+leng+1-1)=suffix//newSuffix
      else
         
         filename(period+1:period+1+leng-1)=newSuffix
         filename(period+1+leng:)=EmptyString(max(len(newSuffix),(period+1+leng)))
      end if
   end subroutine ReplaceSuffix
   
   function EmptyString(n)result(str)
      integer :: n,i
      character(n) :: str
      
      do i = 1,n
         str(i:i)=' '
      end do
   end function EmptyString
   
   function getOrder(num)
      integer num,getOrder
      getOrder=floor(log10(real(num,8)))+1
   end function getOrder
   
   function D2str(num,n,m)
      integer n,m
      double precision num
      character(n) D2str
      write(D2str,'('//'F0.'//I2str(m,getOrder(m))//')')num
   end function D2Str
      
   !shamelessly borrowed from: http://rosettacode.org/wiki/String_case#Fortran
   subroutine To_upper(str)
     character(*), intent(inout) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("a":"z")
           str(i:i) = achar(iachar(str(i:i))-32)
       end select
     end do 
   end subroutine To_upper
 
   subroutine To_lower(str)
     character(*), intent(inout) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do  
   end subroutine To_Lower   
   
   
   function StrAreEqual(str1,str2,caseSensitive)result(res)
      character(*),intent(in) :: str1,str2
      logical,optional,intent(in) :: caseSensitive
      
      logical res
      character(len(str1)) :: one,two
      
      ! if(len(str1)/=len(str2))then
         ! res=.false.
         ! return
      ! end if
      
      one=str1
      two=str2
      call to_upper(one)
      call to_upper(two)
      if(present(caseSensitive))then
         if(.not.caseSensitive)then
            one = (adjustl(trim(one)))
            two = (adjustl(trim(two)))
            res = one == two
         else
            one = adjustl(trim(one))
            two = adjustl(trim(two))
            res = one == two
         end if
      else
         one = (adjustl(trim(one)))
         two = (adjustl(trim(two)))
         res = one == two
      end if
      
   end function StrAreEqual
   
   function splitString(str,n,sep)result(arr)
      character(*),intent(in) :: str,sep
      integer,intent(in) :: n
      character(n),allocatable :: arr(:)
      character(len(str)) str2
      
      integer :: m
      integer :: arr_help(100,2)
      integer :: idx1,idx2,i,idx3,leng,idx
      
      str2=trim(adjustl(str))
      idx1=1
      idx2=index(str2,sep)
      leng=len(str2)
      if(idx2==0)then
         allocate(arr(1))
         arr(1)=str2
         !m=1
         return
      end if
      i=1
      do while(idx2>0)
         arr_help(i,:)=[idx1,idx2]
         idx1=idx2+1
         idx3=idx2
         idx2=index(str2(idx1:),sep)
         i=i+1
         if(idx2==0)then
            arr_help(i,:)=[idx1,leng]
            exit
         end if
         idx2=idx3+idx2
      end do
      m=i
      allocate(arr(m))
      idx=1
      do i = 1,m
         if(IsOnlyWhitespace(str2(arr_help(i,1):arr_help(i,2))))cycle
         arr(idx)=str2(arr_help(i,1):arr_help(i,2))
         idx=idx+1
      end do
   end function splitString
   
   function IsOnlyWhitespace(str)result(res)
      character(*) str
      logical res
      integer i,n
      
      n = len(str)
      
      res=.true.
      do i = 1,n
         if(str(i:i)/=' ')then
            res=.false.
            return
         end if
      end do
      
   end function IsOnlyWhitespace
   
   !will not work with allocatable chars, cause INDEX intrinsic is lame
   function findAnyChar(str,sep,n)result(idx)
      integer :: n,i,idx
      character(*) str,sep(n)
      
      idx=0
      do i = 1,n
         idx = index(str,sep(i))
         if(idx>0)return
      end do
      
   end function findAnyChar
   
   function findString(arr,n,str)result(res)
      integer(4) res,i,n
      character(2) arr(n),str
      
      res=0
      do i = 1,n
         if(str==arr(i))then
            res=i
            exit
         end if
      end do
   end function findString
   
   recursive function I2str(num,n)result(res)
      integer num,n
      character(n) res
      write(res,'('//'I0'//')')num
   end function I2str
   
   function f2str(num,n,m)
      integer n,m
      double precision num
      character(n) f2str
      write(f2str,'('//'F'//TR(i2str(n,2))//'.'//i2str(m,1)//')')num
   end function f2str
   
   function CountSubstring(str,sub)result(num)
      character(*) str,sub
      integer i,sub_len,str_len,num
      
      sub_len=len(sub)
      str_len=len(str)
      num=0
      do i = 1,str_len-sub_len+1
         if(str(i:i+sub_len-1)==sub)then
            num=num+1
         end if
      end do
   end function CountSubstring
   
   subroutine DelFirstChar(str,cha)
      character(*) str,cha
      integer i,idx,leng
      
      leng=len(cha)
      idx=index(str,cha)
      
      do i = 1,leng
         str(idx:)=str(idx+leng:)
      end do
   end subroutine DelFirstChar
   
   function DFC(str,cha)result(str_new)
      character(*) str,cha
      character(len(str)) str_new
      str_new=str
      call DelFirstChar(str_new,cha)
   end function DFC
   
end module strings
