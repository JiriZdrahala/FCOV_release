#define TR(arg) (trim(adjustl(arg)))

module stuff
   use iso_fortran_env
   implicit none
   logical :: nck_arr_allocated=.false.
   integer(int64),allocatable :: nck_arr(:,:)
   
   private nck_arr_allocated
   private nck_arr
   
   contains
   
   
   !rewritten in free-format to not throw warnings
   subroutine INV(a,ai,n,e,IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      dimension ai(n,n),a(n,n),e(n,2*n)
      TOL=1.0d-10
10000 IERR=0
      DO ii=1,n
      DO jj=1,n
      e(ii,jj)=a(ii,jj)
      e(ii,jj+n)=0.0D0
      end do
      end do
      do ii=1,n
      e(ii,ii+n)=1.0D0
      end do
!
      DO ii=1,n-1
      if (ABS(e(ii,ii)).LE.TOL) then
       DO io=ii+1,n
       if (ABS(e(io,ii)).GT.TOL) goto 11
       end do
       IERR=1
       write(6,*)ii
       write(6,*)'tol = ',tol
       tol=tol*0.50d0
       if(tol.gt.1.0d-20)goto 10000
       RETURN
!
11     CONTINUE
       DO kk=1,2*n
       w=e(ii,kk)
       e(ii,kk)=e(io,kk)
       e(io,kk)=w
       end do
      ENDIF
      eii=e(ii,ii)
      DO jj=ii+1,n
      e1=e(jj,ii)/eii
      DO kk=ii+1, 2*n
      e(jj,kk)=e(jj,kk)-e(ii,kk)*e1
      end do
      e(jj,ii)=0.0D0
      end do
      end do
!
      DO i2=n,2,-1
      eii=e(i2,i2)
      DO j2=i2-1,1,-1
      e1=e(j2,i2)/eii
      DO kk=1, n
      e(j2,kk+n)=e(j2,kk+n)-e(i2,kk+n)*e1
      end do
      e(j2,i2)=0.0d0
      end do
      end do
!
      DO ii=1,n
      ei=1.0d0/e(ii,ii)
      DO jj=1,n
      ai(ii,jj)=e(ii,jj+n)*ei
      end do
      end do
!
      RETURN
   END
   
         
   pure subroutine DetPermMat(Perm,N,det)
      implicit none
      integer,intent(in) :: N,Perm(N)
      integer, intent(out) :: det
      integer*4 i,j,buf
      
      det=1
      do i = 1,N
         if(Perm(i)/=i)det=-det
      end do
   end subroutine DetPermMat
   
   subroutine INV_C_lapack(mat,piv,lwork,n)
      integer n,lwork,piv(n)
      double complex mat(n,n)
      !lwork=INV_C_lapack_getLWork(mat,n)
      call INV_C_lapack_do(mat,piv,n,lwork)
   end subroutine INV_C_lapack
   
   subroutine INV_C_lapack_do(mat,piv,n,lwork)
      integer n,ierr,lwork,piv(n)
      double complex mat(n,n),work(lwork)
      call zgetri(n,mat,n,piv,work,lwork,ierr)
      if(ierr/=0)stop 3
   end subroutine INV_C_lapack_do
   
   function INV_C_lapack_getLWork(mat,n)result(lwork)
      integer n,lwork,piv(n),info
      double complex mat(n,n),work(1)
      lwork=-1
      call zgetri(n,mat,n,piv,work,lwork,info)
      lwork=work(1)
   end function INV_C_lapack_getLWork
      
   
   function eps(alpha,beta,gamma)result(res)
      integer alpha,beta,gamma,res
      integer i,j,k
      
      res=0
      if(alpha*beta*gamma /= 6)return
      
      if(alpha==1)then
         if(beta==2)then
            res=1
         else
            res=-1
         end if
      else if(beta==1)then
         if(gamma==2)then
            res=1
         else
            res=-1
         end if
      else if(gamma==1)then
         if(alpha==2)then
            res=1
         else
            res=-1
         end if
      end if
   end function eps
   
   
   subroutine linsort_D_noord(arr,n,descending)
      double precision arr(n)
      integer n,i,j,buf
      logical Descending
      
      if(descending)then
         do i = 1,n
            do j = i,n
               if(arr(j)>arr(i))then
                  buf=arr(i)
                  arr(i)=arr(j)
                  arr(j)=buf
               end if
            end do
         end do
      else
         do i = 1,n
            do j = i,n
               if(arr(j)<arr(i))then
                  buf=arr(i)
                  arr(i)=arr(j)
                  arr(j)=buf
               end if
            end do
         end do
      end if
   end subroutine linsort_D_noord
   
   subroutine linsort_I2(arr,n,order,descending)
      integer(int16) arr(n)
      integer n,i,j,buf,order(n)
      logical Descending
      
      if(descending)then
         do i = 1,n
            do j = i,n
               if(arr(j)>arr(i))then
                  buf=arr(i)
                  arr(i)=arr(j)
                  arr(j)=buf
                  
                  buf=order(i)
                  order(i)=order(j)
                  order(j)=buf
               end if
            end do
         end do
      else
         do i = 1,n
            do j = i,n
               if(arr(j)<arr(i))then
                  buf=arr(i)
                  arr(i)=arr(j)
                  arr(j)=buf
                  
                  buf=order(i)
                  order(i)=order(j)
                  order(j)=buf
               end if
            end do
         end do
      end if
   end subroutine linsort_I2
   
   subroutine linsort_D_big(arr,n,order)
      integer(int64) n
      double precision arr(n),buf
      integer i,j,bufI,order(n)
      
      do i = 1,n
         order(i)=i
      end do
      do i = 1,n
         do j = i,n
            if(arr(j)<arr(i))then
               buf=arr(i)
               arr(i)=arr(j)
               arr(j)=buf
               
               bufI=order(i)
               order(i)=order(j)
               order(j)=bufI
            end if
         end do
      end do
   end subroutine linsort_D_big
   
   function GetSeq(n)result(seq)
      integer i,n,seq(n)
      do i = 1,n
         seq(i)=i
      end do
   end function getseq
   
   function GetSeq_zero(n)result(seq)
      integer i,n,seq(n)
      do i = 0,n-1
         seq(i+1)=i
      end do
   end function GetSeq_zero

   function GetSizeStr(bytes)result(SizeStr)
      integer(int64) bytes
      character(20) SizeStr
      
      if(bytes<1024)then !B
         write(SizeStr,'(I8," B")')bytes
      else if(bytes<1024**2)then !kB
         write(SizeStr,'(F8.2," kB")')dble(bytes)/1024.0d0
      else if(bytes<1024**3)then !MB
         write(SizeStr,'(F8.2," MB")')dble(bytes)/1024.0d0**2
      else if(bytes<1024_8**4)then !GB
         write(SizeStr,'(F8.2," GB")')dble(bytes)/1024.0d0**3
      else !TB
         write(SizeStr,'(F8.2," TB")')dble(bytes)/1024.0d0**4
      end if
   end function GetSizeStr
   
   function GetSizeStr_D(bytes)result(SizeStr)
      double precision bytes
      character(20) SizeStr
      
      if(bytes<1024)then !B
         write(SizeStr,'(F8.2," B")')bytes
      else if(bytes<1024**2)then !kB
         write(SizeStr,'(F8.2," kB")')(bytes)/1024.0d0
      else if(bytes<1024**3)then !MB
         write(SizeStr,'(F8.2," MB")')(bytes)/1024.0d0**2
      else if(bytes<1024_8**4)then !GB
         write(SizeStr,'(F8.2," GB")')(bytes)/1024.0d0**3
      else !TB
         write(SizeStr,'(F8.2," TB")')(bytes)/1024.0d0**4
      end if
   end function GetSizeStr_D
   
   function GetTimeStr(seconds)result(TimeStr)
      double precision seconds,hours,mins
      character(120) TimeStr
      
      
      if(seconds<60.0)then
         write(TimeStr,'(F8.2," seconds")')seconds
      else if(seconds<3600.0)then
         mins=seconds/60.0
         write(TimeStr,'(F8.2," minutes")')mins
      else
         hours=seconds/3600.0
         write(TimeStr,'(F8.2," hours")')hours
      end if
      TimeStr=TR(TimeStr)
   end function GetTimeStr
   
   function TimeStr_ZeroFront(vall)result(str)
      integer vall
      character(10) str
      
      if(vall<10)then
         write(str,'(A1,I1)')"0",vall
      else if(vall<100)then
         write(str,'(I2)')vall
      else
         write(str,'(I10)')vall
      end if
      str=TR(str)
   end function TimeStr_ZeroFront
   
   function GetTimeStr_butBig(seconds,short)result(TimeStr)
      integer hours,mins
      double precision seconds,seconds_sub
      logical short
      character(2) secs_str,mins_str
      character(80) TimeStr
      
      hours=FLOOR(seconds)/3600 !integer division baybee
      mins=FLOOR(seconds-3600*hours)/60
      seconds_sub=(seconds-3600*hours-60*mins)
      
      if(short)then
         secs_str=TimeStr_ZeroFront(floor(seconds_sub))
         mins_str=TimeStr_ZeroFront(mins)
         write(TimeStr,'(I8,":",A2,":",A2)')hours,mins_str,secs_str
      else
         write(TimeStr,'(I8," hours",1X,I2," minutes",1X,F4.1," seconds")')hours,mins,seconds_sub
      end if
      TimeStr=TR(TimeStr)
   end function GetTimeStr_butBig
   
   
         
   function Ciphers(num)result(cif)
      integer num,cif
      if(num==0)then
         cif=1
         return
      end if
      cif=floor(log10(real(abs(num),8)))+1
   end function Ciphers
   
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
      
   subroutine Prepare_NCK_arr(max_n,max_k)
      integer(int32) n,k,max_n,max_k
      allocate(nck_arr(0:max_n,0:max_k))
      
      do n = 0,max_n
         do k = 0,max_k
            nck_arr(n,k)=nck(n,k)
         end do
      end do
      nck_arr_allocated=.true.
   end subroutine Prepare_NCK_arr
   
   subroutine Dispose_NCK_arr()
      deallocate(nck_arr)
      nck_arr_allocated=.false.
   end subroutine Dispose_NCK_arr
   
   
   function nck(n,k)result(res)
      integer(int32),intent(in) :: n,k
      integer(int64) res,i
      if(nck_arr_allocated)then
         res=nck_arr(n,k)
         return
      end if
      
      if(N<k)then
         res=0
         return
      else if(N==k)then
         res=1
         return
      end if
      
      res=n
      do i = n-1,n-k+1,-1
         res=res*i
      end do
      res=res/fac(k)
   end function nck
   
   function fac(n)result(res)
      integer(int32),intent(in) :: n
      integer(int64) res
      integer(int32) i
      
      res=1
      if(n<=1)return
      
      do i = 2,n
         res=res*i
      end do
   end function fac
   
   function minutes(secs)result(mins)
      double precision secs
      integer mins
      mins=floor(secs/60.0d0)
   end function minutes
   
   function hours(s)result(hou)
      double precision s
      integer hou
      hou=floor(s/3600.0d0)
   end function hours
   
   
end module stuff


module ANALYTICAL
   use stuff
   implicit none
   contains
   
   function FC_v1_0_v2_1(FC_00,D1)result(fc)
      double precision FC_00,D1,fc
      fc=(sqrt(2d0)*D1*FC_00)/2d0
   end function FC_v1_0_v2_1
   
   function FC_v1_0_v2_2(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(sqrt(2d0)*FC_00*(D1**2/2d0 + C1))/2d0
   end function FC_v1_0_v2_2
   
   function FC_v1_0_v2_3(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(sqrt(3d0)*FC_00*(D1**3/6d0 + C1*D1))/2d0
   end function FC_v1_0_v2_3
   
   function FC_v1_0_v2_4(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(sqrt(2d0)*sqrt(3d0)*FC_00*(C1**2/2d0 + (C1*D1**2)/2d0 + D1**4/24d0))/2d0
   end function FC_v1_0_v2_4
   
   function FC_v1_0_v2_5(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(sqrt(15d0)*FC_00*((C1**2*D1)/2d0 + (C1*D1**3)/6d0 + D1**5/120d0))/2d0
   end function FC_v1_0_v2_5
   
   function FC_v1_0_v2_6(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(3d0*sqrt(5d0)*FC_00*(C1**3/6d0 + (C1**2*D1**2)/4d0 + (C1*D1**4)/24d0 + D1**6/720d0))/2d0
   end function FC_v1_0_v2_6
   
   function FC_v1_0_v2_7(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(3d0*sqrt(2d0)*sqrt(35d0)*FC_00*((C1**3*D1)/6d0 + (C1**2*D1**3)/12d0 + (C1*D1**5)/120d0 + D1**7/5040d0))/4d0
   end function FC_v1_0_v2_7
   
   function FC_v1_0_v2_8(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(3d0*sqrt(2d0)*sqrt(35d0)*FC_00*(C1**4/24d0 + (C1**3*D1**2)/12d0 + (C1**2*D1**4)/48d0 + (C1*D1**6)/720d0 + D1**8/40320d0))/2d0
   end function FC_v1_0_v2_8
   
   function FC_v1_0_v2_9(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(9d0*sqrt(35d0)*FC_00*((C1**4*D1)/24d0 + (C1**3*D1**3)/36d0 + (C1**2*D1**5)/240d0 + (C1*D1**7)/5040d0 + D1**9/362880d0))/2d0
   end function FC_v1_0_v2_9
   
   function FC_v1_0_v2_10(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=(45d0*sqrt(7d0)*FC_00*(C1**5/120d0 + (C1**4*D1**2)/48d0 + (C1**3*D1**4)/144d0 + (C1**2*D1**6)/1440d0 + (C1*D1**8)/40320d0 + D1**10/3628800d0))/2d0
   end function FC_v1_0_v2_10
   
   
   function FC_v1_0_v2_11(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=FC_00*((144115188075855872d0*C1**5*D1)/123873416550785415d0 + (72057594037927936d0*C1**4*D1**3)/74324049930471249d0 + (72057594037927936d0*C1**3*D1**5)/371620249652356245d0 + &
       (36028797018963968d0*C1**2*D1**7)/2601341747566493715d0 + (9007199254740992d0*C1*D1**9)/23412075728098443435d0 + (4503599627370496d0*D1**11)/1287664165045414388925d0)
   end function FC_v1_0_v2_11
   
   function FC_v1_0_v2_12(FC_00,C1,D1)result(fc)
      double precision FC_00,C1,D1,fc
      fc=FC_00*((18014398509481984d0*C1**6)/37928332905582105d0 + (18014398509481984d0*C1**5*D1**2)/12642777635194035d0 + (4503599627370496d0*C1**4*D1**4)/7585666581116421d0 + &
       (9007199254740992d0*C1**3*D1**6)/113784998716746315d0 + (1125899906842624d0*C1**2*D1**8)/265498330339074735d0 + (1125899906842624d0*C1*D1**10)/11947424865258363075d0 + &
       (281474976710656d0*D1**12)/394265020553525981475d0)
   end function FC_v1_0_v2_12
   
   function FC_v1_0_v2_13(FC_00,C_11,D1)result(fc)
      double precision C_11,D1,FC_00,fc
      fc=FC_00*((36028797018963968d0*C_11**6*D1)/29753432246708685d0 + (36028797018963968d0*C_11**5*D1**3)/29753432246708685d0 + &
      (9007199254740992d0*C_11**4*D1**5)/29753432246708685d0 + (18014398509481984d0*C_11**3*D1**7)/624822077180882385d0 + &
      (2251799813685248d0*C_11**2*D1**9)/1874466231542647155d0 + (2251799813685248d0*C_11*D1**11)/103095642734845593525d0 &
      + (562949953421312d0*D1**13)/4020730066658978147475d0)
   end Function FC_v1_0_v2_13
   
   function FC_v1_0_v2_14(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((576460752303423488d0*C1**7)/1259522918006415615d0 + (288230376151711744d0*C1**6*D1**2)/179931845429487945d0 + &
      (144115188075855872d0*C1**5*D1**4)/179931845429487945d0 + (72057594037927936d0*C1**4*D1**6)/539795536288463835d0 + &
      (36028797018963968d0*C1**3*D1**8)/3778568754019246845d0 + (18014398509481984d0*C1**2*D1**10)/56678531310288702675d0 + &
      (9007199254740992d0*C1*D1**12)/1870391533239527188275d0 + (4503599627370496d0*D1**14)/170205629524796974133025d0)
   end Function FC_v1_0_v2_14
   
   function FC_v1_0_v2_15(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((1152921504606846976d0*C1**7*D1)/919825485182458155d0 + (576460752303423488d0*C1**6*D1**3)/394210922221053495d0 + &
      (288230376151711744d0*C1**5*D1**5)/657018203701755825d0 + (144115188075855872d0*C1**4*D1**7)/2759476455547374465d0 + &
      (72057594037927936d0*C1**3*D1**9)/24835288099926370185d0 + (36028797018963968d0*C1**2*D1**11)/455313615165316786725d0 + &
      (18014398509481984d0*C1*D1**13)/17757230991447354682275d0 + (9007199254740992d0*D1**15)/1864509254101972241638875d0)
   end Function FC_v1_0_v2_15
   
   function FC_v1_0_v2_16(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((288230376151711744d0*C1**8)/650414838080722365d0 + (1152921504606846976d0*C1**7*D1**2)/650414838080722365d0 + &
      (288230376151711744d0*C1**6*D1**4)/278749216320309585d0 + (288230376151711744d0*C1**5*D1**6)/1393746081601547925d0 + &
      (36028797018963968d0*C1**4*D1**8)/1951244514242167095d0 + (72057594037927936d0*C1**3*D1**10)/87806003140897519275d0 + &
      (18014398509481984d0*C1**2*D1**12)/965866034549872712025d0 + (18014398509481984d0*C1*D1**14)/87893809144038416794275d0 + &
      (1125899906842624d0*D1**16)/1318407137160576251914125d0)
   end Function FC_v1_0_v2_16
   
   function FC_v1_0_v2_17(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((1152921504606846976d0*C1**8*D1)/892361795892302925d0 + (4611686018427387904d0*C1**7*D1**3)/2677085387676908775d0 + &
      (1152921504606846976d0*C1**6*D1**5)/1912203848340649125d0 + (1152921504606846976d0*C1**5*D1**7)/13385426938384543875d0 + &
      (144115188075855872d0*C1**4*D1**9)/24093768489092178975d0 + (288230376151711744d0*C1**3*D1**11)/1325157266900069843625d0 + &
      (72057594037927936d0*C1**2*D1**13)/17227044469700907967125d0 + (72057594037927936d0*C1*D1**15)/1808839669318595336548125d0 + &
      (4503599627370496d0*D1**17)/30750274378416120721318125d0)
   end Function FC_v1_0_v2_17
   
   function FC_v1_0_v2_18(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((9223372036854775808d0*C1**9)/21416683101415269255d0 + (4611686018427387904d0*C1**8*D1**2)/2379631455712807695d0 + &
      (9223372036854775808d0*C1**7*D1**4)/7138894367138423085d0 + (4611686018427387904d0*C1**6*D1**6)/15297630786725192325d0 + &
      (1152921504606846976d0*C1**5*D1**8)/35694471835692115425d0 + (576460752303423488d0*C1**4*D1**10)/321250246521229038825d0 + &
      (576460752303423488d0*C1**3*D1**12)/10601258135200558281225d0 + (288230376151711744d0*C1**2*D1**14)/321571496767750267863825d0 + &
      (36028797018963968d0*C1*D1**16)/4823572451516254017957375d0 + (18014398509481984d0*D1**18)/738006585081986864747478375d0)
   end Function FC_v1_0_v2_18
   
   function FC_v1_0_v2_19(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((18446744073709551616d0*C1**9*D1)/13896979074475086555d0 + (9223372036854775808d0*C1**8*D1**3)/4632326358158362185d0 + &
      (18446744073709551616d0*C1**7*D1**5)/23161631790791810925d0 + (9223372036854775808d0*C1**6*D1**7)/69484895372375432775d0 + &
      (2305843009213693952d0*C1**5*D1**9)/208454686117126298325d0 + (1152921504606846976d0*C1**4*D1**11)/2293001547288389281575d0 + &
      (1152921504606846976d0*C1**3*D1**13)/89427060344247181981425d0 + (576460752303423488d0*C1**2*D1**15)/3129947112048651369349875d0 + &
      (72057594037927936d0*C1*D1**17)/53209100904827073278947875d0 + (36028797018963968d0*D1**19)/9098756254725429530700086625d0)
   end function FC_v1_0_v2_19
   
   function FC_v1_0_v2_20(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((36893488147419103232d0*C1**10)/87892212942080024625d0 + (36893488147419103232d0*C1**9*D1**2)/17578442588416004925d0 + &
      (9223372036854775808d0*C1**8*D1**4)/5859480862805334975d0 + (36893488147419103232d0*C1**7*D1**6)/87892212942080024625d0 + &
      (4611686018427387904d0*C1**6*D1**8)/87892212942080024625d0 + (4611686018427387904d0*C1**5*D1**10)/1318383194131200369375d0 + &
      (1152921504606846976d0*C1**4*D1**12)/8701329081265922437875d0 + (2305843009213693952d0*C1**3*D1**14)/791820946395198941846625d0 + &
      (144115188075855872d0*C1**2*D1**16)/3959104731975994709233125d0 + (144115188075855872d0*C1*D1**18)/605743023992327190512668125d0 + &
      (36028797018963968d0*D1**20)/57545587279271083098703471875d0)
   end Function FC_v1_0_v2_20
   
   function FC_v1_0_v2_21(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((147573952589676412928d0*C1**10*D1)/108496503140355821925d0 + (147573952589676412928d0*C1**9*D1**3)/65097901884213493155d0 + &
      (36893488147419103232d0*C1**8*D1**5)/36165501046785273975d0 + (147573952589676412928d0*C1**7*D1**7)/759475521982490753475d0 + &
      (18446744073709551616d0*C1**6*D1**9)/976468528263202397325d0 + (18446744073709551616d0*C1**5*D1**11)/17901923018158710617625d0 + &
      (4611686018427387904d0*C1**4*D1**13)/139634999541637942817475d0 + (9223372036854775808d0*C1**3*D1**15)/14661674951871983995834875d0 + &
      (576460752303423488d0*C1**2*D1**17)/83082824727274575976397625d0 + (576460752303423488d0*C1*D1**19)/14207163028363952491963993875d0 + &
      (144115188075855872d0*D1**21)/1491752117978215011656219356875d0)
   end Function FC_v1_0_v2_21
   
   function FC_v1_0_v2_22(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((295147905179352825856d0*C1**11)/719684383964353263225d0 + (147573952589676412928d0*C1**10*D1**2)/65425853087668478475d0 + &
      (73786976294838206464d0*C1**9*D1**4)/39255511852601087085d0 + (36893488147419103232d0*C1**8*D1**6)/65425853087668478475d0 + &
      (36893488147419103232d0*C1**7*D1**8)/457980971613679349325d0 + (18446744073709551616d0*C1**6*D1**10)/2944163388945081531375d0 + &
      (9223372036854775808d0*C1**5*D1**12)/32385797278395896845125d0 + (4611686018427387904d0*C1**4*D1**14)/589421510466805322581275d0 + &
      (1152921504606846976d0*C1**3*D1**16)/8841322657002079838719125d0 + (576460752303423488d0*C1**2*D1**18)/450907455507106071774675375d0 + &
      (288230376151711744d0*C1*D1**20)/42836208273175076818594160625d0 + (144115188075855872d0*D1**22)/9895164111103442745095251104375d0)
   end Function FC_v1_0_v2_22
   
   function FC_v1_0_v2_23(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((1180591620717411303424d0*C1**11*D1)/848893387086690578325d0 + (590295810358705651712d0*C1**10*D1**3)/231516378296370157725d0 + &
      (295147905179352825856d0*C1**9*D1**5)/231516378296370157725d0 + (147573952589676412928d0*C1**8*D1**7)/540204882691530368025d0 + &
      (147573952589676412928d0*C1**7*D1**9)/4861843944223773312225d0 + (73786976294838206464d0*C1**6*D1**11)/38200202418901076024625d0 + &
      (36893488147419103232d0*C1**5*D1**13)/496602631445713988320125d0 + (18446744073709551616d0*C1**4*D1**15)/10428655260359993754722625d0 + &
      (4611686018427387904d0*C1**3*D1**17)/177287139426119893830284625d0 + (2305843009213693952d0*C1**2*D1**19)/10105366947288833948326223625d0 + &
      (1152921504606846976d0*C1*D1**21)/1061063529465327564574253480625d0 + (36028797018963968d0*D1**23)/16778067059670493135625701929447d0)
   end Function FC_v1_0_v2_23
   
   function FC_v1_0_v2_24(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((295147905179352825856d0*C1**12)/735163238321691015375d0 + (590295810358705651712d0*C1**11*D1**2)/245054412773897005125d0 + &
      (147573952589676412928d0*C1**10*D1**4)/66833021665608274125d0 + (147573952589676412928d0*C1**9*D1**6)/200499064996824822375d0 + &
      (18446744073709551616d0*C1**8*D1**8)/155943717219752639625d0 + (73786976294838206464d0*C1**7*D1**10)/7017467274888868783125d0 + &
      (18446744073709551616d0*C1**6*D1**12)/33082345724476095691875d0 + (18446744073709551616d0*C1**5*D1**14)/1003497820309108235986875d0 + &
      (1152921504606846976d0*C1**4*D1**16)/3010493460927324707960625d0 + (2305843009213693952d0*C1**3*D1**18)/460605499521880680317975625d0 + &
      (576460752303423488d0*C1**2*D1**20)/14585840818192888210069228125d0 + (828662331436171264d0*C1*D1**22)/4843410766691176235934339103095d0 + &
      (1125899906842624d0*D1**24)/3632558075018381784046452942205d0)
   end Function FC_v1_0_v2_24
   
   function FC_v1_0_v2_25(FC_00,C1,D1)result(fc)
      double precision C1,D1,FC_00,fc
      fc=FC_00*((4722366482869645213696d0*C1**12*D1)/3326969031016509734625d0 + (9444732965739290427392d0*C1**11*D1**3)/3326969031016509734625d0 + &
      (2361183241434822606848d0*C1**10*D1**5)/1512258650462049879375d0 + (2361183241434822606848d0*C1**9*D1**7)/6351486331940609493375d0 + &
      (295147905179352825856d0*C1**8*D1**9)/6351486331940609493375d0 + (1180591620717411303424d0*C1**7*D1**11)/349331748256733522135625d0 + &
      (295147905179352825856d0*C1**6*D1**13)/1946276883144658194755625d0 + (295147905179352825856d0*C1**5*D1**15)/68119690910063036816446875d0 + &
      (18446744073709551616d0*C1**4*D1**17)/231606949094214325175919375d0 + (36893488147419103232d0*C1**3*D1**19)/39604788295110649605082213125d0 + &
      (145844570332766142464d0*C1**2*D1**21)/21918775022075301474374107486905d0 + (432345564227567616d0*C1*D1**23)/16439081266556474327698687067795d0 + &
      (2251799813685248d0*D1**25)/51372128957988984941181237407935d0)
   end Function FC_v1_0_v2_25
   
   
   function FC_c1(v,fc_00,C_ii,D_i)result(fc)
      integer v
      double precision C_ii,D_i,fc,fc_00
      
      select case(v)
         case(0)
            fc=0 !intentionally made 0, I dont want to return <0|0>
         case(1)
            fc=FC_v1_0_v2_1(Fc_00,D_i)
         case(2)
            fc=FC_v1_0_v2_2(Fc_00,C_ii,D_i)
         case(3)
            fc=FC_v1_0_v2_3(Fc_00,C_ii,D_i)
         case(4)
            fc=FC_v1_0_v2_4(Fc_00,C_ii,D_i)
         case(5)
            fc=FC_v1_0_v2_5(Fc_00,C_ii,D_i)
         case(6)
            fc=FC_v1_0_v2_6(Fc_00,C_ii,D_i)
         case(7)
            fc=FC_v1_0_v2_7(Fc_00,C_ii,D_i)
         case(8)
            fc=FC_v1_0_v2_8(Fc_00,C_ii,D_i)
         case(9)
            fc=FC_v1_0_v2_9(Fc_00,C_ii,D_i)
         case(10)
            fc=FC_v1_0_v2_10(Fc_00,C_ii,D_i)
         case(11)
            fc=FC_v1_0_v2_11(Fc_00,C_ii,D_i)
         case(12)
            fc=FC_v1_0_v2_12(Fc_00,C_ii,D_i)
         case(13)
            fc=FC_v1_0_v2_13(Fc_00,C_ii,D_i)
         case(14)
            fc=FC_v1_0_v2_14(Fc_00,C_ii,D_i)
         case(15)
            fc=FC_v1_0_v2_15(Fc_00,C_ii,D_i)
         case(16)
            fc=FC_v1_0_v2_16(Fc_00,C_ii,D_i)
         case(17)
            fc=FC_v1_0_v2_17(Fc_00,C_ii,D_i)
         case(18)
            fc=FC_v1_0_v2_18(Fc_00,C_ii,D_i)
         case(19)
            fc=FC_v1_0_v2_19(Fc_00,C_ii,D_i)
         case(20)
            fc=FC_v1_0_v2_20(Fc_00,C_ii,D_i)
         case(21)
            fc=FC_v1_0_v2_21(Fc_00,C_ii,D_i)
         case(22)
            fc=FC_v1_0_v2_22(Fc_00,C_ii,D_i)
         case(23)
            fc=FC_v1_0_v2_23(Fc_00,C_ii,D_i)
         case(24)
            fc=FC_v1_0_v2_24(Fc_00,C_ii,D_i)
         case(25)
            fc=FC_v1_0_v2_25(Fc_00,C_ii,D_i)
         case default
            write(output_unit,*)'Analytical form not available for class 1, v = ',v
            call exit(112)
      end select
   end function FC_c1
   
   
   
   
   !These were done in Matlab
   !It is pretty slow to obtain them past class 3
   !I tried parallelizing the code, but failed, since I dont have that experience with Matlab
   !Maybe if I knew before-hand that I do not need some terms and deleted them?
   
   
   !<0|1,1,1>
   function FC3_v1_0_v2_1_1_1(FC_00,C2_1,C3_1,C3_2,D1,D2,D3)result(fc)
      double precision fc,fc_00,C2_1,C3_1,C3_2,D1,D2,D3
      fc=(2d0**(0.5d0)*FC_00*(2d0*C2_1*D3 + 2d0*C3_1*D2 + 2d0*C3_2*D1 + D1*D2*D3))/4d0
   end function FC3_v1_0_v2_1_1_1
   
   !<0|2,1,1>
   function FC3_v1_0_v2_2_1_1(FC_00,C1_1,C2_1,C3_1,C3_2,D1,D2,D3)result(fc)
      double precision fc,fc_00,C1_1,C2_1,C3_1,C3_2,D1,D2,D3
      fc=(2d0**(0.5d0)*FC_00*(C3_2*D1**2 + 2*C1_1*C3_2 + 4*C2_1*C3_1 + C1_1*D2*D3 + 2*C2_1*D1*D3 + 2*C3_1*D1*D2 + (D1**2*D2*D3)/2d0))/4d0
   end function FC3_v1_0_v2_2_1_1
   
   !<0|3,1,1>
   function FC3_v1_0_v2_3_1_1(FC_00,C1_1,C2_1,C3_1,C3_2,D1,D2,D3)result(fc)
      double precision fc,fc_00,C1_1,C2_1,C3_1,C3_2,D1,D2,D3
      fc=(3d0**(0.5d0)*FC_00*((C3_2*D1**3)/3d0 + 2d0*C1_1*C2_1*D3 + 2d0*C1_1*C3_1*D2 + 2d0*C1_1*C3_2*D1 + 4d0*C2_1*C3_1*D1 + C2_1*D1**2*D3 + C3_1*D1**2*D2 + (D1**3*D2*D3)/6d0 + C1_1*D1*D2*D3))/4d0
   end function FC3_v1_0_v2_3_1_1
   
   !<0|1,1,1,1>
   function FC4_v1_0_v2_1_1_1_1(FC_00,C2_1,C3_1,C3_2,C4_1,C4_2,C4_3,D1,D2,D3,D4)result(fc)
      double precision fc,fc_00
      double precision C2_1,C3_1,C3_2,C4_3,C4_2,C4_1
      double precision D1,D2,D3,D4
      fc=FC_00*(C2_1*C4_3 + C3_1*C4_2 + C3_2*C4_1 + (C2_1*D3*D4)/2d0 + (C3_1*D2*D4)/2d0 + (C3_2*D1*D4)/2d0 + (C4_1*D2*D3)/2d0 + (C4_2*D1*D3)/2d0+ (C4_3*D1*D2)/2d0 + (D1*D2*D3*D4)/4d0)
   end function FC4_v1_0_v2_1_1_1_1
   
   !<1,0,0|1,1,1>
   function FC3_v1_1_0_0_v2_1_1_1(fc_00,B1,C2_1,C3_1,C3_2,D1,D2,D3,E1_1,E2_1,E3_1)result(fc)
      double precision fc,fc_00
      double precision B1
      double precision C2_1,C3_1,C3_2
      double precision D1,D2,D3
      double precision E1_1,E2_1,E3_1
      fc=FC_00*((C3_2*E1_1)/2d0 + (C2_1*E3_1)/2d0 + (C3_1*E2_1)/2d0 + (B1*C2_1*D3)/2d0 + (B1*C3_1*D2)/2d0 + (B1*C3_2*D1)/2d0 + (D2*D3*E1_1)/4d0 + (D1*D3*E2_1)/4d0 + (D1*D2*E3_1)/4d0 + (B1*D1*D2*D3)/4d0)
   end function FC3_v1_1_0_0_v2_1_1_1
   
   
   
   !<2,0|1,1>
   function FC2_v1_2_0_v2_1_1(fc_00,A1_1,B1,C2_1,D1,D2,E1_1,E2_1)result(fc)
      double precision fc,fc_00
      double precision A1_1
      double precision B1
      double precision C2_1
      double precision D1,D2
      double precision E1_1,E2_1
      fc=(sqrt(2.0)*FC_00*(B1**2*C2_1 + 2*A1_1*C2_1 + E1_1*E2_1 + A1_1*D1*D2 + B1*D2*E1_1 + B1*D1*E2_1 + (B1**2*D1*D2)/2d0))/4d0
   end function FC2_v1_2_0_v2_1_1
   
   !<2,0,0|1,1,1>
   function FC3_v1_2_0_0_v2_1_1_1(fc_00,A1_1,B1,C2_1,C3_1,C3_2,D1,D2,D3,E1_1,E2_1,E3_1)result(fc)
      double precision fc,fc_00
      double precision A1_1
      double precision B1
      double precision C2_1,C3_1,C3_2
      double precision D1,D2,D3
      double precision E1_1,E2_1,E3_1
      fc=FC_00*((A1_1*C2_1*D3)/2d0 + (A1_1*C3_1*D2)/2d0 + (A1_1*C3_2*D1)/2d0 + (B1*C3_2*E1_1)/2d0 + (B1*C2_1*E3_1)/2d0 + (B1*C3_1*E2_1)/2d0 + &
      (D3*E1_1*E2_1)/4d0 + (D2*E1_1*E3_1)/4d0 + (D1*E2_1*E3_1)/4d0 + (B1**2*C2_1*D3)/4d0 + (B1**2*C3_1*D2)/4d0 + (B1**2*C3_2*D1)/4d0 + (A1_1*D1*D2*D3)/4d0 + &
      (B1*D2*D3*E1_1)/4d0 + (B1*D1*D3*E2_1)/4d0 + (B1*D1*D2*E3_1)/4d0 + (B1**2*D1*D2*D3)/8d0)
   end function FC3_v1_2_0_0_v2_1_1_1
   
   !<1,1,0|1,1,1>
   function FC3_v1_1_1_0_v2_1_1_1(fc_00,A2_1,B1,B2,C2_1,C3_1,C3_2,D1,D2,D3,E1_1,E1_2,E2_1,E2_2,E3_1,E3_2)result(fc)
      double precision fc,fc_00
      double precision A2_1
      double precision B1,B2
      double precision C2_1,C3_1,C3_2
      double precision D1,D2,D3
      double precision E1_1,E1_2,E2_1,E3_1,E3_2,E2_2
      fc=(sqrt(2.0)*FC_00*(4*A2_1*C2_1*D3 + 4*A2_1*C3_1*D2 + 4*A2_1*C3_2*D1 + 2*B1*C3_2*E1_2 + 2*B2*C3_2*E1_1 + 2*B1*C2_1*E3_2 + 2*B1*C3_1*E2_2 + 2*B2*C2_1*E3_1 + 2*B2*C3_1*E2_1 &
      + D3*E1_1*E2_2 + D3*E1_2*E2_1 + D2*E1_1*E3_2 + D2*E1_2*E3_1 + D1*E2_1*E3_2 + D1*E2_2*E3_1 + 2*B1*B2*C2_1*D3 + 2*B1*B2*C3_1*D2 + 2*B1*B2*C3_2*D1 + 2*A2_1*D1*D2*D3 + B1*D2*D3*E1_2 &
      + B2*D2*D3*E1_1 + B1*D1*D3*E2_2 + B2*D1*D3*E2_1 + B1*D1*D2*E3_2 + B2*D1*D2*E3_1 + B1*B2*D1*D2*D3))/8d0
   end function FC3_v1_1_1_0_v2_1_1_1
   
   
   
end module ANALYTICAL

module Fcov_primitive
   use constants
   use stuff
   implicit none
   
   type, public :: Mode
      integer,allocatable :: v(:)
      integer(int16),allocatable :: v_pos(:)
      double precision :: w = HUGE(1d0)
   end type Mode
   
   type, public :: Mode_col
      type(Mode),allocatable :: vs(:)
   end type Mode_col
   
   type, public :: Transition_R
      integer,allocatable :: vi(:),vf(:)
      integer(int16),allocatable :: vi_pos(:),vf_pos(:)
      double precision :: wif=HUGE(1d0)
      double precision :: ram=huge(1d0),roa=huge(1d0)
   end type Transition_R
   
   double precision,allocatable, public :: sqrt_arr(:)
   
   
   !public A,B,C,D,E,J,K,SG,SE,wg,we
   
   contains
   
   
   
   
   !TODO check, not sure about exc_m13
   subroutine v_union(v1,v1_pos,v1_n,v3,v3_pos,v3_n,v13,v13_pos,v13_n,exc_m13,mc_gfs,nq)
      integer,intent(in) :: v1_n,v3_n
      integer,intent(in) :: v1(v1_n),v3(v3_n)
      integer(int16),intent(in) :: v1_pos(v1_n),v3_pos(v3_n)
      integer,intent(out) :: v13_n
      integer(int16),allocatable,intent(out) :: v13_pos(:)
      integer,allocatable,intent(out) :: v13(:),exc_m13(:,:)
      integer v13_pos_help(v1_n+v3_n),v13_help(v1_n+v3_n),whereFrom_help(v1_n+v3_n)
      integer v1_nz(v1_n),v3_nz(v3_n),nq
      
      integer,allocatable :: order(:),whereFrom(:)
      integer i,j,idx,idx2,v1_clas,v3_clas
      integer c_gfs,mc_gfs
      
      
      v1_clas=FCClass(v1,v1_n,v1_nz,v1_n)
      v3_clas=FCClass(v3,v3_n,v3_nz,v3_n)
      v13_pos_help=HUGE(int(1,kind=kind(v1_pos)))
      idx=1
      do i = 1,v1_n
         v13_pos_help(idx)=v1_pos(i)
         v13_help(idx)=v1(i)
         whereFrom_help(idx)=1
         idx=idx+1
      end do
      
      do i = 1,v3_n
         idx2=findloc(v13_pos_help,v3_pos(i),dim=1)
         if(idx2==0)then
            v13_pos_help(idx)=v3_pos(i)
            v13_help(idx)=v3(i)
            whereFrom_help(idx)=3
            idx=idx+1
         else
            v13_help(idx2)=max(v13_help(idx2),v3(i))
            whereFrom_help(idx2)=3
         end if
      end do
      
      idx=idx-1
      v13_n=idx
      allocate(v13(v13_n),v13_pos(v13_n),whereFrom(v13_n),order(v13_n))
      do i = 1,v13_n
         v13_pos(i)=v13_pos_help(i)
         v13(i)=v13_help(i)
         whereFrom(i)=whereFrom_help(i)
      end do
      order=getseq(v13_n)
      call linsort_I2(v13_pos,v13_n,order,.false.)
      v13_help(1:v13_n)=v13
      whereFrom_help(1:v13_n)=whereFrom
      do i = 1,v13_n
         v13(i)=v13_help(order(i))
         whereFrom(i)=whereFrom_help(order(i))
      end do
      
      mc_gfs=max(v1_clas,v3_clas)
      allocate(exc_m13(nq,mc_gfs))
      exc_m13=0
      ! do i = 1,v1_n
         ! idx=findloc(v13_pos,v1_pos(i),dim=1)
         ! do c_gfs = 1,v1_clas
            ! exc_m13(idx,c_gfs)=max(exc_m13(idx,c_gfs),v1(i))
         ! end do
      ! end do
      
      ! do i = 1,v3_n
         ! idx=findloc(v13_pos,v3_pos(i),dim=1)
         ! do c_gfs = 1,v3_clas
            ! exc_m13(idx,c_gfs)=max(exc_m13(idx,c_gfs),v3(i))
         ! end do
      ! end do
      
      do i = 1,v1_n
         idx=v1_pos(i)
         do c_gfs = 1,v1_clas
            exc_m13(idx,c_gfs)=max(exc_m13(idx,c_gfs),v1(i))
         end do
      end do
      
      do i = 1,v3_n
         idx=v3_pos(i)
         do c_gfs = 1,v3_clas
            exc_m13(idx,c_gfs)=max(exc_m13(idx,c_gfs),v3(i))
         end do
      end do
      
      
      deallocate(order)
   end subroutine v_union
   
   function MakeMode(v,v_pos,n)result(m)
      integer(int16) v_pos(n)
      integer v(n),n
      type(Mode) m
      m%v=v
      m%v_pos=v_pos
   end function MakeMode
   
   
   
   
   subroutine MakeSqrtArr(N)
      integer N,i
      
      allocate(sqrt_arr(0:N))
      do i = 0,N
         sqrt_arr(i)=sqrt(dble(i))
      end do
   end subroutine MakeSqrtArr
   
   !incase someone uses valgrind on this
   subroutine UnmakeSqrtArr()
      deallocate(sqrt_arr)
   end subroutine UnmakeSqrtArr
      
   
   
   function LOMC_Wrap(n,k,comb_n)result(storage)
      integer n,k
      integer(int16),allocatable :: storage(:,:)
      integer(int64),intent(in) :: comb_n
      
      integer i,ii,stor_idx
      integer(int16), allocatable :: seq(:),arr(:)
      
      allocate(seq(k),storage(k,comb_n))
      arr=GetSeq(n)
      i=1
      ii=1
      stor_idx=1
      call LOMC(seq,arr,n,k,i,ii,storage,stor_idx)
      deallocate(seq,arr)
   end function LOMC_Wrap

   
   !Loop Over Mode Combinations
   recursive subroutine LOMC(sequence,arr,n,k,i,ii,stor,stor_idx)
      integer n,k,j
      integer,value :: i,ii
      integer :: idx
      integer(int16) :: sequence(k),arr(n),stor(:,:)
      integer :: stor_idx
      
      if(i>k)then
         stor(:,stor_idx)=sequence
         stor_idx=stor_idx+1
         return
      end if
      
      do j = ii,n
         sequence(i)=arr(j)
         call LOMC(sequence,arr,n,k,i+1,j+1,stor,stor_idx)
      end do
   end subroutine LOMC
      
   function FCClass(v,n,nonzero,nonzero_n)result(clas)
      integer n,v(n),i,v_pos_help(n),clas,nonzero_n,nonzero(nonzero_n)
      clas=0
      do i = 1,N
         if(v(i)/=0)then
            clas=clas+1
            nonzero(clas)=i
         end if
      end do
   end function FCClass
   
   function FCClass_fromIdx(idx,clas_starts)result(clas)
      integer clas,i
      integer(int64) idx,clas_starts(:)
      
      clas=0
      if(idx==1)return
      clas=1
      do while(idx>clas_starts(clas+1)) !off by one error?
         clas=clas+1
      end do
      clas=clas
   end function FCClass_fromIdx
   
   pure subroutine LowerMode_inplace(v,i,by,N,neg)
      logical,intent(out) :: neg
      integer,intent(in) :: N,by,i
      integer,intent(inout) :: v(N)
      neg=.false.
      v(i)=v(i)-by
      if(v(i)<0)then
         neg=.true.
      end if
   end subroutine LowerMode_inplace
   
   pure subroutine IncrementMode_inplace(v,i,by,N)
      integer,intent(in) :: N,by,i
      integer,intent(inout) :: v(N)
      v(i)=v(i)+by
   end subroutine IncrementMode_inplace
   
   function getFreq_short(v,v_pos,w,n)result(res)
      integer n,v(n),i
      integer(int16) v_pos(n)
      double precision :: w(:),res
      
      res = 0
      do i = 1,n
         res=res+w(v_pos(i))*v(i)
      end do
   end function getFreq_short
   
end module Fcov_primitive

module FCOV_storage
   use Fcov_primitive
   use ANALYTICAL
   use stuff
   use util
   implicit none
   
   type,public :: v_col
      integer,allocatable :: arr(:)
   end type v_col
   
   type,public :: Modes_arr_t
      integer :: clas
      integer(int16), allocatable :: arr(:,:)
   end type Modes_arr_t
   
   type, public :: FC_CombsIndices
      integer :: clas = -1
      integer(int64),allocatable :: idxs(:)
   end type FC_CombsIndices
   
   type, public :: FC_storSys_Dusch
      double precision,allocatable :: A(:,:),B(:),C(:,:),D(:),E(:,:)
   end type FC_storSys_Dusch
   
   type, public :: FC_storSys_derivs
      double precision,allocatable :: du(:,:),dm(:,:),dq(:,:,:)
   end type FC_storSys_derivs
   
   !storage system
   type, public :: FC_storSys
      integer :: nq=0
      integer,allocatable :: red_dims(:,:),exc_m(:,:),v_map(:)
      integer,allocatable :: v_map_rd(:)
      double precision, allocatable :: w(:),sqrt_w(:)
      integer :: red_n=0
      integer(int64) :: fc_arr_size=1
      
      integer :: mc_v
      integer(int64),allocatable :: cs(:)
      type(FC_CombsIndices),allocatable :: fc_ci(:)
   end type FC_storSys
   
   type Polar_exc
      double complex,allocatable :: ap_fcht(:,:,:)
      double complex,allocatable :: ap(:,:,:),G(:,:,:),Gc(:,:,:),A(:,:,:,:),Ac(:,:,:,:)
   end type Polar_exc
   
   integer,private :: Polar_nexc = -1
   
   interface operator(+)
      module procedure :: Polars_exc_add
   end interface
      
   interface operator(*)
      module procedure :: Polars_mult_sc
   end interface
      
   !$OMP DECLARE REDUCTION(+:Polar_exc:omp_out=omp_out + omp_in) INITIALIZER ( omp_priv = OMP_POLAR_NEW() )
   
   
   !Vibronic system definition
   !One ground state and one excited state
   type, public :: ExcState
      integer :: root
      integer :: NQ=-1,nq_gr=-1,nq_ex=-1
      double precision,allocatable :: wg_gr(:)
      integer :: n_gr
      double precision,allocatable :: A(:,:),B(:),C(:,:),D(:),E(:,:),J(:,:),K(:),SG(:,:),SE(:,:),J_i(:,:),K_i(:)
      double precision,allocatable :: A_gr(:,:),B_gr(:),C_gr(:,:),D_gr(:),E_gr(:,:),J_gr(:,:),K_gr(:),SG_gr(:,:),SE_gr(:,:),J_gr_i(:,:),K_gr_i(:)
      double precision,allocatable :: wg(:),we(:),we_gr(:)
      integer,allocatable :: wg_gr_idx(:),wg_ex_idx(:)
      double precision,allocatable :: sqrt_w(:),sqrt_wg(:)
      integer,allocatable :: w_map(:)
      integer,allocatable :: del_modes(:)
      character(2) :: TMExpand
      logical :: ExpandInExc1,ExpandinExc2
      
      !transition (multi)poles
      !u = dipole, length formalism
      !v = dipole, velocity formalism
      !q = quadrupole
      !m = magnetic dipole
      double precision :: m(3),q(3,3),v(3),u(3)
      double precision :: q_raw(6)
      
      double precision :: m_gr(3),q_gr(3,3),v_gr(3),u_gr(3) !<ex_el|M(Q')|gr_el>
      
      double precision :: m_gr_tr(3),q_gr_tr(3,3),v_gr_tr(3),u_gr_tr(3) !<ex_el|M^tr(Q')|gr_el>
      
      double precision :: m_ex(3),q_ex(3,3),v_ex(3),u_ex(3) !<ex_el|M(Q'')|gr_el>
      
      double precision :: m_ex_tr(3),q_ex_tr(3,3),v_ex_tr(3),u_ex_tr(3) !<ex_el|M^tr(Q'')|gr_el>
      
      !transition (multi)pole first derivatives
      !du = derivative dipole - length formalism, it will have dimension(9*nat)
      !dv = derivative dipole - velocity formalism, dimension(9*nat)
      !dq = derivative quadrupole - dimension(18*nat)
      !dm = derivative magnetic dipole - dimension(9*nat)
      double precision,allocatable :: du(:,:),dv(:,:),dq(:,:,:),dm(:,:)
      double precision,allocatable :: du_raw(:),dv_raw(:),dq_raw(:),dm_raw(:)
      
      double precision,allocatable :: du_ex(:,:),dv_ex(:,:),dq_ex(:,:,:),dm_ex(:,:)
      double precision,allocatable :: du_ex_tr(:,:),dv_ex_tr(:,:),dq_ex_tr(:,:,:),dm_ex_tr(:,:)
      
      double precision,allocatable :: du_gr(:,:),dv_gr(:,:),dq_gr(:,:,:),dm_gr(:,:)
      double precision,allocatable :: du_gr_tr(:,:),dv_gr_tr(:,:),dq_gr_tr(:,:,:),dm_gr_tr(:,:)
      
      !transition (multi)pole second derivatives
      double precision,allocatable :: du2_ex(:,:,:),dv2_ex(:,:,:),dq2_ex(:,:,:,:),dm2_ex(:,:,:)
      double precision,allocatable :: du2_ex_tr(:,:,:),dv2_ex_tr(:,:,:),dq2_ex_tr(:,:,:,:),dm2_ex_tr(:,:,:)
      
      double precision,allocatable :: du2_gr(:,:,:),dv2_gr(:,:,:),dq2_gr(:,:,:,:),dm2_gr(:,:,:)
      double precision,allocatable :: du2_gr_tr(:,:,:),dv2_gr_tr(:,:,:),dq2_gr_tr(:,:,:,:),dm2_gr_tr(:,:,:)
      
      type(Polar_exc),allocatable :: Polars(:,:)
      type(Polar_exc),allocatable :: elpolars(:)
      logical :: dusch_gauss = .false.,antiStokes=.false.
      double precision :: uncoup_lim
      
      double precision :: FC_00
      double precision :: e_00,e_00_gr,e00_vertical,e00_vertical_exc
      
      !FC storage
      integer :: mc_gs,mc_fs,mc_ms,mc_gfs
      integer :: mc_v1,mc_v2
      integer,allocatable :: max_v1s(:),max_v3s(:),max_v2s(:),min_v1s(:),min_v2s(:)
      integer,allocatable :: max_v2s_ps(:),min_v2s_ps(:)
      integer :: max_osc_v2_algo=0
      double precision :: max_osc_v2_coeff=0.49
      integer,allocatable :: max_osc_v2(:)
      
      integer,allocatable :: min_v2s_excm(:)
      integer :: mc_ms_ps
      double precision :: min_freq_mode=-HUGE(1d0),max_freq_mode=HUGE(1d0)
      
      
      double precision,allocatable :: fc_arr(:,:)
      integer,allocatable :: fc_arr_check(:,:)
      integer,allocatable :: red_dims1(:,:),red_dims2(:,:)
      integer,allocatable :: v1_map(:),v2_map(:),v1_map_w(:),v2_map_w(:)
      logical,allocatable :: im_modes_gr(:),im_modes_ex_gr(:),im_modes_ex_ex(:)
      integer :: red_n1,red_n2,red_n2_max,red_n2_max_algo,fc_sum_min_ps_its=10
      
      integer,allocatable :: exc_m1(:,:),exc_m2(:,:)
      double precision :: thr_v1,gamma,theta,eps,fc_sum_min_ps=20d0
      double precision :: thr_v2,thr_v2_comb,thr_v2_over,thr_v2_other,thr_v2_0
      
   end type ExcState
      
   contains
   
   function PolarExc2Polar(polar_exc_o,iexc)result(polarr)
      integer iexc
      type(Polar) polarr
      type(polar_exc) polar_exc_o
      
      polarr%ap=polar_exc_o%ap(:,:,iexc)
      polarr%G=polar_exc_o%G(:,:,iexc)
      polarr%Gc=polar_exc_o%Gc(:,:,iexc)
      polarr%A=polar_exc_o%A(:,:,:,iexc)
      polarr%Ac=polar_exc_o%Ac(:,:,:,iexc)
      
   end function PolarExc2Polar
   
   subroutine AllocatePolar(polar1,nexc)
      integer nexc
      type(Polar_exc) polar1
      allocate(polar1%Ap(3,3,nexc),polar1%G(3,3,nexc),polar1%Gc(3,3,nexc),polar1%A(3,3,3,nexc),polar1%Ac(3,3,3,nexc))
      allocate(polar1%Ap_fcht(3,3,nexc))
   end subroutine AllocatePolar
   
   subroutine DeallocatePolar(polar1)
      type(Polar_exc) polar1
      deallocate(polar1%Ap,polar1%G,polar1%Gc,polar1%A,polar1%Ac,polar1%Ap_fcht)
   end subroutine DeallocatePolar
   
   
   
   function MakeExcMaxUnion_1ExcSt(exc_m_a,mc_a,exc_m_b,mc_b,nq)result(exc_m_c)
      integer nq,mc_a,mc_b,exc_m_a(nq,mc_a),exc_m_b(nq,mc_b)
      integer,allocatable :: exc_m_c(:,:)
      
      integer ii,clas
      
      allocate(exc_m_c(nq,max(mc_a,mc_b)))
      do clas = 1,min(mc_a,mc_b)
         do ii = 1,nq
            exc_m_c(ii,clas)=max(exc_m_a(ii,clas),exc_m_b(ii,clas))
         end do
      end do
      clas=clas-1
      
      if(mc_a<mc_b)then
         exc_m_c(:,clas:)=exc_m_b(:,clas:)
      else if(mc_b<mc_a)then
         exc_m_c(:,clas:)=exc_m_a(:,clas:)
      end if
   end function MakeExcMaxUnion_1ExcSt
   
   function FC_storSys_size(fc_ss)result(sz)
      type(fc_storSys) fc_ss
      integer(int64) sz
      integer i
      
      sz=0
      sz=kind(fc_ss%nq)+kind(fc_ss%red_n)+kind(fc_ss%fc_arr_size)
      sz=sz+kind(fc_ss%mc_v)
      sz=sz+kind(fc_ss%exc_m)*size(fc_ss%exc_m,dim=1)*size(fc_ss%exc_m,dim=2)
      if(allocated(fc_ss%red_dims))sz=sz+kind(fc_ss%red_dims)*size(fc_ss%red_dims,dim=1)*size(fc_ss%red_dims,dim=2)
      if(allocated(fc_ss%v_map))sz=sz+kind(fc_ss%v_map)*size(fc_ss%v_map,dim=1)
      if(allocated(fc_ss%w))sz=sz+kind(fc_ss%w)*size(fc_ss%w,dim=1)
      sz=sz+kind(fc_ss%cs)*size(fc_ss%cs,dim=1)
      if(allocated(fc_ss%fc_ci))then
         do i = 2,fc_ss%mc_v
            sz=sz+int(kind(fc_ss%fc_ci(i)%idxs),kind=8)*size(fc_ss%fc_ci(i)%idxs,dim=1,kind=8)
            sz=sz+kind(fc_ss%fc_ci(i)%clas)
         end do
      end if
      if(allocated(fc_ss%sqrt_w))sz=sz+kind(fc_ss%sqrt_w)*size(fc_ss%sqrt_w,dim=1)
   end function FC_storSys_size
   
   function FC_storSys_make(exc_m,exc_m_mc,nq,true_red_n,mc,w)result(fc_ss)
      integer nq,mc,exc_m_mc,exc_m(nq,exc_m_mc),i
      double precision,optional :: w(nq)
      type(FC_storSys) fc_ss
      logical true_red_n
      
      fc_ss%exc_m=exc_m
      fc_ss%mc_v=mc
      fc_ss%nq=nq
      if(true_red_n)then
         fc_ss%red_n=Count_red_n(exc_m,exc_m_mc,nq,mc,fc_ss%v_map)
      else 
         fc_ss%red_n=nq
         fc_ss%v_map=GetSeq(nq)
      end if
      if(mc>0)then
         call Make_red_dims(exc_m,exc_m_mc,fc_ss%red_n,nq,mc,fc_ss%red_dims,fc_ss%v_map)
         call PrepareIdxs(mc,fc_ss%red_n,fc_ss%red_dims,fc_ss%fc_arr_size,fc_ss%fc_ci,fc_ss%cs)
         if(present(w))then
            fc_ss%w=MakeW(w,nq,fc_ss%v_map,fc_ss%red_n)
            allocate(fc_ss%sqrt_w(fc_ss%red_n))
            do i = 1,fc_ss%red_n
               fc_ss%sqrt_w(i)=dsqrt(1d0/(2d0*fc_ss%w(i)))
            end do
         end if
      end if
      
   end function FC_storSys_make
   
   function FC_storSys_make_0(nq)result(fc_ss)
      integer nq
      type(FC_storSys) fc_ss
      
      fc_ss%mc_v=0
      fc_ss%nq=nq
      fc_ss%red_n=0
   end function FC_storSys_make_0
   
   function FC_storSys_make_uncoup(exc_m,exc_m_mc,nq,true_red_n,mc,uncoup_m,w)result(fc_ss)
      integer nq,mc,exc_m_mc,exc_m(nq,exc_m_mc),i,j
      double precision,optional :: w(nq)
      type(FC_storSys) fc_ss
      logical true_red_n,uncoup_m(nq)
      
      allocate(fc_ss%exc_m(nq,exc_m_mc))
      do i = 1,nq
         if(uncoup_m(i))then
            fc_ss%exc_m(i,:)=0
         else
            do j = 1,exc_m_mc
               fc_ss%exc_m(i,j)=exc_m(i,j)
            end do
         end if
      end do
      fc_ss%mc_v=mc
      fc_ss%nq=nq
      if(true_red_n)then
         fc_ss%red_n=Count_red_n(exc_m,exc_m_mc,nq,mc,fc_ss%v_map)
      else 
         fc_ss%red_n=nq
         fc_ss%v_map=GetSeq(nq)
      end if
      if(mc>0)then
         call Make_red_dims(exc_m,exc_m_mc,fc_ss%red_n,nq,mc,fc_ss%red_dims,fc_ss%v_map)
         call PrepareIdxs(mc,fc_ss%red_n,fc_ss%red_dims,fc_ss%fc_arr_size,fc_ss%fc_ci,fc_ss%cs)
         if(present(w))then
            fc_ss%w=MakeW(w,nq,fc_ss%v_map,fc_ss%red_n)
         end if
      end if
   end function FC_storSys_make_uncoup
   
   function IsDecreasingSequence(arr,n)result(res)
      integer n,arr(n),buf,i
      logical res
      
      res=.false.
      do i = 2,n
         if(arr(i)>arr(i-1))return
      end do
      
      res=.true.
   end function IsDecreasingSequence
   
   function FC_storSys_make_minmaxVs(exc_m,exc_m_mc,nq,min_vs,max_vs,true_red_n,full_ht,w2sqrt,mc,w, &
                                     max_osc_v2,max_osc_v2_algo,max_osc_coeff,fc_sharp_vec)result(fc_ss)
      integer nq,mc,exc_m_mc,exc_m(nq,exc_m_mc),min_vs(mc),max_vs(mc)
      double precision :: w(nq),w2sqrt(nq),max_osc_coeff,fc_sharp_vec(nq)
      double precision, allocatable :: fc_sharp_vec_red(:)
      type(FC_storSys) fc_ss
      logical true_red_n,full_ht
      integer i,j
      
      integer :: max_osc_v2_algo,osc_max,mc_nonzero
      integer,allocatable :: max_osc_v2(:),exc_m_1mode(:)
      integer,allocatable :: red_dims_slice(:)
      
      
      do i = 1,nq
         exc_m_1mode=exc_m(i,:)
         if(.not.IsDecreasingSequence(exc_m_1mode,exc_m_mc))then
            write(Output_unit,*)'Exc M is not a decreasing sequence for mode ',i
            call exit(22)
         end if
      end do
      deallocate(exc_m_1mode)
      fc_ss%exc_m=exc_m
      fc_ss%mc_v=mc
      fc_ss%nq=nq
      if(true_red_n)then
         fc_ss%red_n=Count_red_n(exc_m,exc_m_mc,nq,mc,fc_ss%v_map)
      else 
         fc_ss%red_n=nq
         fc_ss%v_map=GetSeq(nq)
      end if
      call Make_red_dims(exc_m,exc_m_mc,fc_ss%red_n,nq,mc,fc_ss%red_dims,fc_ss%v_map)
      
      
      ! if(allocated(max_osc_v2))then
         ! fc_ss%red_n=min(fc_ss%red_n,maxval(max_osc_v2)) !oh no, if I will ever make the red_n(:) scheme, Im effed
      ! end if
      do i = 1,fc_ss%red_n
         exc_m_1mode=fc_ss%red_dims(i,:)
         if(.not.IsDecreasingSequence(exc_m_1mode,mc))then
            write(Output_unit,*)'ERROR 24: Red Dims is not a decreasing sequence for mode ',i
            call exit(24)
         end if
      end do
      
      do i = 1,mc
         do j = 1,fc_ss%red_n
            if(max_vs(i)<0)then
               fc_ss%red_dims(j,i)=-max_vs(i)
            else
               fc_ss%red_dims(j,i)=min(fc_ss%red_dims(j,i),max_vs(i))
            end if
            fc_ss%red_dims(j,i)=max(fc_ss%red_dims(j,i),min_vs(i))
         end do
         
         if(allocated(max_osc_v2))then
            osc_max=max_osc_v2(i)
         else
            osc_max=999
         end if
         if(count(fc_ss%red_dims(:,i)>0)>osc_max)then
            allocate(red_dims_slice(fc_ss%red_n),fc_sharp_vec_red(fc_ss%red_n))
            
            do j = 1,fc_ss%red_n
               fc_sharp_vec_red(j)=fc_sharp_vec(fc_ss%v_map(j))
            end do
            
            red_dims_slice=fc_ss%red_dims(:,i)
            
            
            call FC_storSys_smaller_RedDims(red_dims_slice,i,fc_ss%red_n,osc_max,max_osc_v2_algo,max_osc_coeff,fc_sharp_vec_red)
            
            fc_ss%red_dims(:,i)=red_dims_slice
            
            deallocate(red_dims_slice,fc_sharp_vec_red)
         end if
      end do
      
      do i = 1,fc_ss%red_n
         exc_m_1mode=fc_ss%red_dims(i,:)
         if(.not.IsDecreasingSequence(exc_m_1mode,mc))then
            write(Output_unit,*)'ERROR 23: Red Dims is not a decreasing sequence for mode ',i
            call exit(23)
         end if
      end do
      
      if(mc==0)then
         fc_ss%fc_arr_size=1
         allocate(fc_ss%cs(1))
         fc_ss%cs(1)=1
      else
         call PrepareIdxs(mc,fc_ss%red_n,fc_ss%red_dims,fc_ss%fc_arr_size,fc_ss%fc_ci,fc_ss%cs)
      end if
      
      if(mc==0)then
         
      else
         fc_ss%w=MakeW(w,nq,fc_ss%v_map,fc_ss%red_n)
      end if
      
      if(full_HT)then
         allocate(fc_ss%v_map_rd(nq))
         fc_ss%v_map_rd=0
         do i = 1,fc_ss%red_n
            fc_ss%v_map_rd(fc_ss%v_map(i))=i
         end do
      else
         allocate(fc_ss%sqrt_w(fc_ss%red_n))
         do i = 1,fc_ss%red_n
            fc_ss%sqrt_w(i)=dsqrt(1d0/(2d0*w2sqrt(fc_ss%v_map(i))))
         end do
      end if
   end function FC_storSys_make_minmaxVs
   
   subroutine FC_storSys_smaller_RedDims(red_dims_1c,clas,red_n,max_red_n,algo_type,coeff,FC_Sharp_vec)
      integer algo_type,red_n,max_red_n,cur_red_n
      integer red_dims_1c(red_n),clas,idx,vi
      character(3) cl_str
      double precision :: coeff,FC_Sharp_vec(red_n)
      logical :: maskk(red_n)
      
      
      select case(algo_type)
         case default
         case(0) ! * coeff 
            red_dims_1c=FLOOR(red_dims_1c*coeff)
            cur_red_n=count(red_dims_1c>0)
            if(cur_red_n>max_red_n)then
               cl_str=ClassString(clas)
               write(output_unit,'("WARNING: Class: red_n > max_osc_v2","",A,I4," > ",I4)')TR(cl_str),red_n,max_red_n
            end if
         case(1) ! remove smallest excitations
            cur_red_n=count(red_dims_1c>0)
            do while(cur_red_n>max_red_n)
               idx=minloc(red_dims_1c,dim=1,back=.true.)
               red_dims_1c(idx)=0
               cur_red_n=count(red_dims_1c>0)
            end do
         case(2) ! remove excitations with smaller B (ground states) or D (excited states) vector
            maskk=.true.
            cur_red_n=count(red_dims_1c>0)
            do while(cur_red_n>max_red_n)
               idx=minloc(FC_Sharp_vec,dim=1,mask=maskk)
               red_dims_1c(idx)=0
               maskk(idx)=.false.
               cur_red_n=count(red_dims_1c>0)
            end do
         case(3) ! remove excitations with lowest excitations from highest frequencies first
            cur_red_n=count(red_dims_1c>0)
            vi=1
            do while(cur_red_n>max_red_n)
               idx=findloc(red_dims_1c==vi,.true.,dim=1,back=.true.) !OPTimize this
               if(idx==0)then
                  vi=vi+1
                  cycle
               end if
               red_dims_1c(idx)=0
               cur_red_n=count(red_dims_1c>0)
            end do
         case(4) ! remove excitations from highest frequencies first
            cur_red_n=count(red_dims_1c>0)
            vi=1
            do while(cur_red_n>max_red_n)
               idx=findloc(red_dims_1c>0,.true.,dim=1,back=.true.) !OPTimize this
               if(idx==0)then
                  exit
               end if
               red_dims_1c(idx)=0
               cur_red_n=count(red_dims_1c>0)
            end do
      end select
   end subroutine FC_storSys_smaller_RedDims
   
   function ClassString(clas)result(clas_str)
      integer clas
      character(3) clas_str
      
      if(clas<=9)then
         write(clas_str,'("C",I1)')clas
      else
         write(clas_str,'("C",I2)')clas
      end if

   end function ClassString
   
   subroutine FC_storSys_dispose(fc_ss)
      type(fc_storSys) fc_ss
      if(allocated(fc_ss%red_dims))deallocate(fc_ss%red_dims)
      deallocate(fc_ss%exc_m)
      if(allocated(fc_ss%v_map))deallocate(fc_ss%v_map)
      if(allocated(fc_ss%fc_ci))deallocate(fc_ss%fc_ci) !this one might not be allocated if mc_v <= 2
      deallocate(fc_ss%cs)
      if(allocated(fc_ss%w))deallocate(fc_ss%w) !this one is not used so often
      if(allocated(fc_ss%v_map_rd))deallocate(fc_ss%v_map_rd)
      if(allocated(fc_ss%sqrt_w))deallocate(fc_ss%sqrt_w)
   end subroutine FC_storSys_dispose
   
   subroutine FC_storSys_dusch_dispose(fc_ss_dusch)
      type(fc_storSys_dusch) fc_ss_dusch
      deallocate(fc_ss_dusch%A)
      deallocate(fc_ss_dusch%B)
      deallocate(fc_ss_dusch%C)
      deallocate(fc_ss_dusch%D)
      deallocate(fc_ss_dusch%E)
   end subroutine FC_storSys_dusch_dispose
   
   function Modes_arr_t_size(modes_arr_t_obj)result(siz)
      type(modes_arr_t),allocatable :: modes_arr_t_obj(:)
      integer(int64) siz
      integer i
      
      siz=0
      do i = 1,size(modes_arr_t_obj,dim=1)
         siz=siz+size(modes_arr_t_obj(i)%arr)*kind(modes_arr_t_obj(i)%arr)
         siz=siz+kind(modes_arr_t_obj(i)%clas)
      end do
   end function Modes_arr_t_size
   
   subroutine Modes_arr_t_dispose(modes_arr_t_obj,mc)
      integer mc,i
      type(modes_arr_t),allocatable :: modes_arr_t_obj(:)
      
      do i = 1,mc
         deallocate(modes_arr_t_obj(i)%arr)
      end do
      deallocate(modes_arr_t_obj)
   end subroutine Modes_arr_t_dispose
   

   function FC_storSys_Dusch_make(fc_ss1,fc_ss2,A,B,C,D,E,nq)result(fc_ss_dusch)
      integer nq
      double precision A(nq,nq),B(nq),C(nq,nq),D(nq),E(nq,nq)
      type(FC_storSys) fc_ss1,fc_ss2
      type(FC_storSys_Dusch) FC_ss_dusch
      integer i,j
      
      allocate(fc_ss_dusch%A(fc_ss1%red_n,fc_ss1%red_n), &
      fc_ss_dusch%B(fc_ss1%red_n), &
      fc_ss_dusch%C(fc_ss2%red_n,fc_ss2%red_n), &
      fc_ss_dusch%D(fc_ss2%red_n), &
      fc_ss_dusch%E(fc_ss1%red_n,fc_ss2%red_n))
      
      do i = 1,fc_ss1%red_n
         fc_ss_dusch%B(i)=B(fc_ss1%v_map(i))
         do j = 1,fc_ss1%red_n
            fc_ss_dusch%A(i,j)=A(fc_ss1%v_map(i),fc_ss1%v_map(j))
         end do
         do j = 1,fc_ss2%red_n
            fc_ss_dusch%E(i,j)=E(fc_ss1%v_map(i),fc_ss2%v_map(j))
         end do
      end do
      
      do i = 1,fc_ss2%red_n
         fc_ss_dusch%D(i)=D(fc_ss2%v_map(i))
         do j = 1,fc_ss2%red_n
            fc_ss_dusch%C(i,j)=C(fc_ss2%v_map(i),fc_ss2%v_map(j)) !maybe? the indices might be switched here
         end do
      end do
   end function FC_storSys_Dusch_make
   
   
   subroutine Make_red_dims(exc_m,exc_m_mc,red_n,nq,mc,red_dims,v_map)
      integer exc_m_mc,nq,mc,exc_m(nq,exc_m_mc),red_n,v_map(red_n)
      integer,allocatable :: red_dims(:,:)
      integer clas,i,idx
      
      if(mc==0)return
      allocate(red_dims(red_n,mc))
      red_dims=0
      if(red_n==nq)then
         do clas=1,min(exc_m_mc,mc)
            red_dims(:,clas)=exc_m(:,clas)
         end do
         do clas=min(exc_m_mc,mc)+1,mc
            red_dims(:,clas)=red_dims(:,exc_m_mc)
         end do
      else
         do clas=1,min(exc_m_mc,mc)
            !idx=1
            do i = 1,red_n
               !if(exc_m(i,clas)>0)then
                  red_dims(i,clas)=exc_m(v_map(i),clas)
                  !idx=idx+1
               !end if
            end do
         end do
         do clas=min(exc_m_mc,mc)+1,mc
            red_dims(:,clas)=red_dims(:,exc_m_mc)
         end do
      end if
   end subroutine Make_red_dims
   
   function Count_red_n(exc_m,exc_m_mc,nq,mc,v_map)result(red_n)
      integer nq,mc,exc_m(nq,exc_m_mc),red_n,red_n_idx
      integer ii,counts(exc_m_mc),clas,idx,exc_m_mc
      integer,allocatable,intent(out) :: v_map(:)
      
      if(mc==0 .or. nq==0)then
         red_n=0
         return
      end if
      
      do clas = 1,exc_m_mc
         counts(clas)=count(exc_m(:,clas)/=0)
      end do
      red_n_idx=maxloc(counts,dim=1)
      red_n=counts(red_n_idx)
      allocate(v_map(red_n))
      
      idx=1
      do ii = 1,nq
         if(exc_m(ii,red_n_idx)==0)cycle
         v_map(idx)=ii
         idx=idx+1
      end do
   end function Count_red_n
   
   function MakeW(w,nq,v_map,red_n)result(w_red)
      integer nq,red_n,v_map(red_n)
      double precision w(nq),w_red(red_n)
      
      integer i
      
      do i = 1,red_n
         w_red(i)=w(v_map(i))
      end do
   end function MakeW
   
   subroutine Polar_set_nexc(nexc)  
      integer nexc
      Polar_nexc=nexc
   end subroutine Polar_set_nexc
   
   
   subroutine SymmetrizePolarAAc(polarr)
      type(Polar_exc) polarr
      integer aa,iexc
      
      
      do aa = 1,3
         do iexc = 1,size(polarr%ap,dim=3)
            polarr%A(aa,1,3,iexc)=polarr%A(aa,3,1,iexc)
            polarr%A(aa,1,2,iexc)=polarr%A(aa,2,1,iexc)
            polarr%A(aa,2,3,iexc)=polarr%A(aa,3,2,iexc)
            
            polarr%Ac(aa,1,3,iexc)=polarr%Ac(aa,3,1,iexc)
            polarr%Ac(aa,1,2,iexc)=polarr%Ac(aa,2,1,iexc)
            polarr%Ac(aa,2,3,iexc)=polarr%Ac(aa,3,2,iexc)
         end do
      end do
   end subroutine SymmetrizePolarAAc
   
   function OMP_POLAR_NEW()result(polarr)
      integer i,j,k,l
      type(Polar_exc) polarr
      
      allocate(polarr%Ap(3,3,Polar_nexc))
      allocate(polarr%Ap_fcht(3,3,Polar_nexc))
      allocate(polarr%G(3,3,Polar_nexc))
      allocate(polarr%Gc(3,3,Polar_nexc))
      allocate(polarr%A(3,3,3,Polar_nexc))
      allocate(polarr%Ac(3,3,3,Polar_nexc))
      
      do l = 1,Polar_nexc
         do j = 1,3
            do i = 1,3
               polarr%Ap(i,j,l)=(0d0,0d0)
               polarr%Ap_fcht(i,j,l)=(0d0,0d0)
               polarr%G(i,j,l)=(0d0,0d0)
               polarr%Gc(i,j,l)=(0d0,0d0)
               do k = 1,3
                  polarr%A(i,j,k,l)=(0d0,0d0)
                  polarr%Ac(i,j,k,l)=(0d0,0d0)
               end do
            end do
         end do
      end do
      ! polarr%Ap=0
      ! polarr%G =0
      ! polarr%Gc=0
      ! polarr%A =0
      ! polarr%Ac=0
   end function OMP_POLAR_NEW
   
   subroutine Polar_allocate(polarr,nexc)
      integer nexc
      type(Polar_exc) polarr
      
      allocate(polarr%Ap(3,3,nexc))
      allocate(polarr%Ap_fcht(3,3,nexc))
      allocate(polarr%G(3,3,nexc))
      allocate(polarr%Gc(3,3,nexc))
      allocate(polarr%A(3,3,3,nexc))
      allocate(polarr%Ac(3,3,3,nexc))
   end subroutine Polar_allocate
   
   subroutine Polar_new(pol)
      type(Polar_exc) pol
      integer i,j,k,l
      
      do l = 1,Polar_nexc
         do j = 1,3
            do i = 1,3
               pol%Ap(i,j,l)=(0d0,0d0)
               pol%ap_fcht(i,j,l)=(0d0,0d0)
               pol%G(i,j,l)=(0d0,0d0)
               pol%Gc(i,j,l)=(0d0,0d0)
               do k = 1,3
                  pol%A(i,j,k,l)=(0d0,0d0)
                  pol%Ac(i,j,k,l)=(0d0,0d0)
               end do
            end do
         end do
      end do
   end subroutine Polar_new
   
   subroutine Polar_new_nexc(pol,nexc)
      type(Polar_exc) pol
      integer i,j,k,l,nexc
      
      do l = 1,nexc
         do j = 1,3
            do i = 1,3
               pol%Ap(i,j,l)=(0d0,0d0)
               pol%ap_fcht(i,j,l)=(0d0,0d0)
               pol%G(i,j,l)=(0d0,0d0)
               pol%Gc(i,j,l)=(0d0,0d0)
               do k = 1,3
                  pol%A(i,j,k,l)=(0d0,0d0)
                  pol%Ac(i,j,k,l)=(0d0,0d0)
               end do
            end do
         end do
      end do
   end subroutine Polar_new_nexc
   
   
   
   function Polars_exc_add(pol1,pol2)result(polRes)
      type(Polar_exc),intent(in) :: pol1,pol2
      type(Polar_exc) :: polRes
      integer :: i,j,k,l
      
      call Polar_allocate(polres,Polar_nexc)
      do l = 1,Polar_nexc
         do j = 1,3
            do i = 1,3
               polRes%Ap(i,j,l)=pol1%Ap(i,j,l)+pol2%Ap(i,j,l)
               polRes%Ap_fcht(i,j,l)=pol1%Ap_fcht(i,j,l)+pol2%Ap_fcht(i,j,l)
               polRes%G(i,j,l)=pol1%G(i,j,l)+pol2%G(i,j,l)
               polRes%Gc(i,j,l)=pol1%Gc(i,j,l)+pol2%Gc(i,j,l)
               do k = 1,3
                  polRes%A(i,j,k,l)=pol1%A(i,j,k,l)+pol2%A(i,j,k,l)
                  polRes%Ac(i,j,k,l)=pol1%Ac(i,j,k,l)+pol2%Ac(i,j,k,l)
               end do
            end do
         end do
      end do
      ! polres%ap=pol1%ap+pol2%ap
      ! polres%G=pol1%G+pol2%G
      ! polres%Gc=pol1%Gc+pol2%Gc
      ! polres%A=pol1%A+pol2%A
      ! polres%Ac=pol1%Ac+pol2%Ac
   end function Polars_exc_add
   
   
   function Polars_mult_sc(pol1,sc)result(polRes)
      type(Polar_exc),intent(in) :: pol1
      double precision,intent(in) :: sc
      type(Polar_exc) :: polRes
      integer :: i,j,k,l
      
      call Polar_allocate(polres,Polar_nexc)
      do l = 1,Polar_nexc
         do j = 1,3
            do i = 1,3
               polRes%Ap(i,j,l)=pol1%Ap(i,j,l)*sc
               polRes%Ap_fcht(i,j,l)=pol1%Ap_fcht(i,j,l)*sc
               polRes%G(i,j,l)=pol1%G(i,j,l)*sc
               polRes%Gc(i,j,l)=pol1%Gc(i,j,l)*sc
               do k = 1,3
                  polRes%A(i,j,k,l)=pol1%A(i,j,k,l)*sc
                  polRes%Ac(i,j,k,l)=pol1%Ac(i,j,k,l)*sc
               end do
            end do
         end do
      end do
   end function Polars_mult_sc
   
   
   
   
   function MapVPos(v_pos,v_n,v_map,red_n)result(v_pos_real)
      integer(int16) v_pos(v_n),v_pos_real(v_n)
      integer v_n,red_n,v_map(red_n)
      integer i
      
      do i = 1,v_n
         v_pos_real(i)=v_map(v_pos(i))
      end do
   end function MapVPos
   
   function DemapVPos(v_pos_real,v_n,v_map,red_n)result(v_pos)
      integer(int16) v_pos(v_n),v_pos_real(v_n)
      integer v_n,red_n,v_map(red_n)
      integer i
      
      do i = 1,v_n
         v_pos(i)=findloc(v_map,v_pos_real(i),dim=1)
      end do
   end function DemapVPos
   
   
   subroutine PrepareIdxs(max_clas,red_n,red_dims,max_storage,FC_CI,clas_starts)
      integer max_clas,red_n,red_dims(red_n,max_clas),MI_st
      integer(int64),intent(out) :: max_storage
      integer(int64),allocatable :: clas_starts(:)
     ! integer(int64),intent(in) :: MI
      type(FC_CombsIndices),allocatable :: FC_CI(:)
      
      integer,allocatable :: seq(:),seq2(:),prev_seq(:),v1_pos(:)
      integer i,j,k
      integer(int64) :: l
      logical first
      
      if(max_clas==0)then
         allocate(clas_starts(1))
         clas_starts(1)=1
         max_storage=1
         return
      else if(max_clas==1)then
         allocate(clas_starts(2))
         clas_starts(1)=1
         clas_starts(2)=sum(red_dims(:,1))+clas_starts(1)
         max_storage=clas_starts(2)
         return
      end if
      allocate(clas_starts(max_clas+1),FC_CI(2:max_clas))
      clas_starts=1
      v1_pos=GetSeq_zero(red_n)
      
      !Im not touching this part
      !It was made over a 3 days coding binge and worked flawlessly
      do i = 2,2
         allocate(seq(i-1),seq2(i),FC_CI(i)%idxs(nck(red_n,i)),prev_seq(i))
         FC_CI(i)%idxs=-1
         FC_CI(i)%clas=i
         prev_seq=0
         ! prev_seq(i)=1
         j=1
         k=1
         clas_starts(i)=clas_starts(i-1)
         call GetClassEnd(seq,v1_pos+1,red_n,i-1,j,k,red_dims(:,i-1),red_n,clas_starts(i))
         j=1
         k=1
         l=0
         first=.true.
         call PopulatePreviousKCombs(seq2,prev_seq,v1_pos,red_n,i,j,k,red_dims(:,i),red_N,l,FC_CI(i)%idxs,first)
         FC_CI(i)%idxs(1)=0
         deallocate(seq,seq2,prev_seq)
      end do
      
      do i = 3,max_clas
         allocate(seq(i-1))
         clas_starts(i)=clas_starts(i-1)
         j=1
         k=1
         call GetClassEnd(seq,v1_pos+1,red_n,i-1,j,k,red_dims(:,i-1),red_n,clas_starts(i))
         deallocate(seq)
      end do
      
      do i = 3,max_clas
         allocate(seq(i-1),seq2(i),FC_CI(i)%idxs(nck(red_n,i)),prev_seq(i))
         FC_CI(i)%idxs=-1
         FC_CI(i)%clas=i
         prev_seq=0
         ! prev_seq(i)=1
         ! j=1
         ! k=1
         ! clas_starts(i)=clas_starts(i-1)
         ! call GetClassEnd(seq,v1_pos+1,red_n,i-1,j,k,red_dims(:,i-1),red_n,clas_starts(i))
         j=1
         k=1
         l=0
         first=.true.
         call PopulatePreviousKCombs(seq2,prev_seq,v1_pos,red_n,i,j,k,red_dims(:,i),red_N,l,FC_CI(i)%idxs,first)
         FC_CI(i)%idxs(1)=0
         deallocate(seq,seq2,prev_seq)
      end do
      
      if(max_clas==2)then
         j=1
         k=1
         clas_starts(i)=clas_starts(i-1)
         allocate(seq(i-1))
         call GetClassEnd(seq,v1_pos+1,red_n,i-1,j,k,red_dims(:,max_clas),red_n,clas_starts(i))
         deallocate(seq)
         max_storage=clas_starts(i)
      else
         j=1
         k=1
         clas_starts(i)=clas_starts(i-1)
         allocate(seq(i-1))
         call GetClassEnd(seq,v1_pos+1,red_n,i-1,j,k,red_dims(:,max_clas),red_n,clas_starts(i))
         deallocate(seq)
         max_storage=clas_starts(i)
      end if
      deallocate(v1_pos)
   end subroutine PrepareIdxs
   
   function dimprod_vpos(alpha,n,v_pos,full_dims,nq)result(res)
      integer(2) v_pos(n)
      integer alpha,n,nq,full_dims(nq),i,res
      select case(alpha)
         case(2)
            res=full_dims(v_pos(1))
         case(3)
            res=full_dims(v_pos(1))*full_dims(v_pos(2))
         case(4)
            res=full_dims(v_pos(1))*full_dims(v_pos(2))*full_dims(v_pos(3))
         case(5)
            res=full_dims(v_pos(1))*full_dims(v_pos(2))*full_dims(v_pos(3))*full_dims(v_pos(4))
         case default
            res=1
            do i = alpha,2,-1
               res=res*full_dims(v_pos(i-1))
            end do
      end select
   end function dimprod_vpos

   function dimprod_vpos_zero(alpha,n,v_pos,full_dims,nq,nonzeros,clas,mc)result(res)
      integer(int16) v_pos(n)
      integer alpha,n,nq,mc,full_dims(nq,mc),i,res,nonzeros(n),clas
      ! res=1
      ! i=1
      ! do while(nonzeros(i)>0 .and. i<=n-alpha+1) !i = alpha,2,-1
         ! res=res*full_dims(v_pos(nonzeros(i)))
         ! i=i+1
      ! end do
      select case(alpha)
         case(2)
            res=full_dims(v_pos(nonzeros(1)),clas)
         case(3)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)
         case(4)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)
         case(5)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas)
         case(6)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)
         case(7)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)*full_dims(v_pos(nonzeros(6)),clas)
         case(8)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)*full_dims(v_pos(nonzeros(6)),clas)*full_dims(v_pos(nonzeros(7)),clas)
         case(9)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)*full_dims(v_pos(nonzeros(6)),clas)*full_dims(v_pos(nonzeros(7)),clas)*full_dims(v_pos(nonzeros(8)),clas)
         case(10)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)*full_dims(v_pos(nonzeros(6)),clas)*full_dims(v_pos(nonzeros(7)),clas)*full_dims(v_pos(nonzeros(8)),clas) &
            *full_dims(v_pos(nonzeros(9)),clas)
         case(11)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)*full_dims(v_pos(nonzeros(6)),clas)*full_dims(v_pos(nonzeros(7)),clas)*full_dims(v_pos(nonzeros(8)),clas) &
            *full_dims(v_pos(nonzeros(9)),clas)*full_dims(v_pos(nonzeros(10)),clas)
         case(12)
            res=full_dims(v_pos(nonzeros(1)),clas)*full_dims(v_pos(nonzeros(2)),clas)*full_dims(v_pos(nonzeros(3)),clas)*full_dims(v_pos(nonzeros(4)),clas) &
            *full_dims(v_pos(nonzeros(5)),clas)*full_dims(v_pos(nonzeros(6)),clas)*full_dims(v_pos(nonzeros(7)),clas)*full_dims(v_pos(nonzeros(8)),clas) &
            *full_dims(v_pos(nonzeros(9)),clas)*full_dims(v_pos(nonzeros(10)),clas)*full_dims(v_pos(nonzeros(11)),clas)
         case default
            res=1
            do i = alpha,2,-1
               res=res*full_dims(v_pos(nonzeros(i-1)),clas)
            end do
      end select
   end function dimprod_vpos_zero

   
   function ReduceArr_idx_vpos(idxs,n,v_pos,full_dims,nq)result(idx)
      integer n,idxs(n),nq,full_dims(nq)
      integer(int16) v_pos(n)
      integer idx,i,j,prod
      
      idx=idxs(1)
      do i = 2,n
         prod=dimprod_vpos(i,n,v_pos,full_dims,nq)
         idx=idx+(idxs(i)-1)*prod
      end do
   end function ReduceArr_idx_vpos
   
   function UnreduceArr_idx_vpos(idx0,n,v_pos,full_dims,nq)result(idxs)
      integer,value :: idx0
      integer(int16) v_pos(n)
      integer nq,n,full_dims(nq),idxs(n),dp
      integer i
      
      do i=n,2,-1
         dp=dimprod_vpos(i,n,v_pos,full_dims,nq)
         idxs(i)=(idx0-1)/dp+1
         idx0=idx0-dp*(idxs(i)-1)
      end do
      idxs(1)=idx0
   end function UnreduceArr_idx_vpos
   
   function ReduceArr_idx_vpos_zeros(idxs,n,v_pos,full_dims,nq,nonzeros,clas,mc)result(idx)
      integer(int16) v_pos(n)
      integer n,idxs(n),mc,nq,full_dims(nq,mc),nonzeros(n)
      integer idx,i,j,prod,clas

      idx=idxs(nonzeros(1))
      do i = 2,clas
         prod=dimprod_vpos_zero(i,n,v_pos,full_dims,nq,nonzeros,clas,mc)
         idx=idx+(idxs(nonzeros(i))-1)*prod
      end do
   end function ReduceArr_idx_vpos_zeros
   
   function ReduceFC_idx_fcss(v,v_pos,v_n,fc_ss)result(idx)
      integer v_n,v(v_n)
      integer clas,nonzero(v_n)
      integer(int16) v_pos(v_n)
      integer(int64) :: c
      type(FC_storSys) fc_ss
      integer(int64) idx
      
      integer :: i,ii
      integer,allocatable :: red_maxes(:),red_v(:),red_vpos(:)
      
      clas=FCClass(v,v_n,nonzero,v_n)
      select case(clas)
         case(0)
            idx=1
         case(1)
            idx=1
            if(v_pos(nonzero(1))>1)idx=idx+sum(fc_ss%red_dims(1:v_pos(nonzero(1))-1,1))
            idx=idx+v(nonzero(1))
         case default
            ! allocate(red_maxes(clas),red_v(clas))
            ! do i = 1,clas
               ! red_maxes(i)=red_dims(v_pos(nonzero(i)))
               ! red_v(i)=v(nonzero(i))
            ! end do
            c=0
            
            i=KCombs_idx_zeros(v_pos,v_n,nonzero,clas) !1)
            c=fc_ss%fc_ci(clas)%idxs(i) !1)
            ! call GetPreviousKCombs(seq,getseq(red_n),red_n,clas,i,ii,red_vPos,red_dims,red_n,c,sig) !2)
            
            idx=c+ReduceArr_idx_vpos_zeros(v,v_n,v_pos,fc_ss%red_dims,fc_ss%red_n,nonzero,clas,fc_ss%mc_v)+fc_ss%cs(clas)
            !idx=c+ReduceArr_idx(red_v,clas,red_maxes)+clas_starts(clas)
            ! deallocate(red_maxes,red_v)
      end select
   end function ReduceFC_idx_fcss
   
   
   !it works
   !its just kinda slow btw :(
   !So im just going to iterate through combinations
   subroutine UnReduceFC_idx(idx,v,v_pos,v_n,maxes,mc_v,NQ,clas_starts,clas_idxs)
      integer NQ,mc_v,maxes(NQ,mc_v)
      integer clas,v_n
      integer c
      integer(int64) idx,clas_starts(:)
      type(FC_CombsIndices) clas_idxs(2:)
      
      integer :: i,ii,jj,prev_i,prev_clases,prev_prev_clases,new_idx,cur,cur_next
      integer :: min_dif2,min_dif,dif
      integer,allocatable :: v(:)
      integer(int16),allocatable :: v_pos(:)
      
      clas=FCClass_fromIdx(idx,clas_starts)
      select case(clas)
         case(0)
            allocate(v(0),v_pos(0))
            v_n=0
         case(1)
            allocate(v(1),v_pos(1))
            v_n=clas
            prev_clases=1
            prev_prev_clases=prev_clases
            ii=0
            do while(prev_clases<idx)
               ii=ii+1
               prev_prev_clases=prev_clases
               prev_clases=prev_clases+maxes(ii,1)
            end do
            v_pos(1)=ii
            v(1)=idx-prev_prev_clases
         case default
            new_idx=idx-clas_starts(clas)
            v_n=clas
            ii=0
            jj=0
            min_dif=HUGE(1)
            min_dif2=min_dif
            do i = 1,size(clas_idxs(clas)%idxs,dim=1)
               cur=clas_idxs(clas)%idxs(i)
               dif=new_idx-cur
               if(cur<=new_idx .and. dif<min_dif)then
                  min_dif=dif
                  ii=i
               end if
               if(cur>=new_idx .and. -dif < min_dif2)then
                  min_dif2=-dif
                  jj=i
               end if
            end do
            if(ii==jj)then !excitations are maximal for given modes
               v_pos=KCombs_idx_rev(ii,clas)+1
               call PrevKComb(v_pos,clas,NQ)
               allocate(v(clas))
               do i = 1,clas
                  v(i)=maxes(v_pos(i),clas)
               end do
            else !excitations are not maximal
               v_pos=KCombs_idx_rev(ii,clas)+1
               new_idx=new_idx-clas_idxs(clas)%idxs(ii)
               v=UnreduceArr_idx_vpos(new_idx,clas,v_pos,maxes(:,clas),nq)
               !v_pos_high=KCombs_idx_rev(jj,clas)+1
            end if
            !new_idx=KCombs_idx([3,4],clas)
      end select
   end subroutine UnReduceFC_idx
   
   function KCombs_idx(arr,k)result(res)
      integer,intent(in) :: k,arr(k)
      integer(int64) res
      integer i
      
      res=1
      do i = 0,k-1
         res=res+nck(arr(k-i),k-i)
      end do
   end function KCombs_idx
   
   function KCombs_idx_zeros(arr,k,nonzeros,clas)result(res)
      integer(int16),intent(in) :: arr(k)
      integer,intent(in) :: k,nonzeros(k),clas
      integer(int64) res
      integer i
      
      res=1
      do i = 0,clas-1
         res=res+nck(arr(nonzeros(clas-i))-1,clas-i)
      end do
   end function KCombs_idx_zeros
   
   function KCombs_idx_rev(idx,k)result(arr)
      integer,intent(in) :: k,idx
      integer arr(k),i,ck
      integer(int64) c_k,c_k_prev,summ
      
      summ=0
      do i = k,1,-1
         ck=-1
         c_k=0
         do while(.not. c_k > idx-1-summ)
            ck=ck+1
            c_k_prev=c_k
            c_k=nck(ck,i)
         end do
         arr(i)=ck-1
         summ=summ+nck(arr(i),i)
      end do
   end function KCombs_idx_rev
   
   !gives previous k-combination
   !assuming some ordering rules (according to makeCombinations)
   subroutine PrevKComb(arr,k,max_n)
      integer(int16) arr(k)
      integer k,max_n
      integer idx
      logical check
      
      if(arr(k)==k)then
         arr=0
         return
      end if
      arr(k)=arr(k)-1
      idx=k
      do while(.not.IsSeq(arr,k))
         arr(idx)=max_n-(k-idx)
         idx=idx-1
         if(idx==0)return
         arr(idx)=arr(idx)-1
      end do
   end subroutine PrevKComb
   
   function IsSeq(arr,n)result(res)
      integer(int16) arr(n)
      integer n,i
      logical res
      
      do i = n-1,1,-1
         if(.not.arr(i)<arr(i+1))then
            res=.false.
            return
         end if
      end do
      res=.true.
   end function IsSeq
   
   recursive subroutine PopulatePreviousKCombs(sequence,prev_seq,arr,n,k,i,ii,maxes,NQ,amount,stor_arr,first)
      integer n,arr(n),k,j,jj
      integer :: i
      integer :: sequence(k)
      integer :: ii
      integer(int64) :: amount,amount2
      integer :: NQ,maxes(NQ)
      integer(int64) :: stor_arr(:)
      integer :: prev_seq(k)
      logical first
      
      if(i>k)then
         ! if(all(prev_seq==0))then
            ! goto 20 !maybe put prev_seq=[0,0,1] when calling this subroutine?
         ! end if
         if(first)goto 20
         amount2=1
         do j = 1,k
            amount2=amount2*maxes(prev_seq(j)+1)
         end do
         amount=amount+amount2
         stor_arr(KCombs_idx(sequence,k))=amount
20       first=.false.
         prev_seq=sequence
         return
      end if
      
      do j = ii,n
         sequence(i) = arr(j)
         call PopulatePreviousKCombs(sequence,prev_seq,arr,n,k,i+1,j+1,maxes,NQ,amount,stor_arr,first)
      end do
   end subroutine PopulatePreviousKCombs
   
   
   recursive subroutine GetClassEnd(sequence,arr,n,k,i,ii,maxes,NQ,amount)
      integer n,arr(n),k,j
      integer :: i
      integer :: sequence(k)
      integer(int64) :: amount
      integer :: ii,amount2
      integer :: NQ,maxes(NQ)
      logical :: ind(k)
      
      if(i>k)then
         amount2=1
         do j = 1,k
            amount2=amount2*maxes(sequence(j))
         end do
         amount=amount+amount2
         return
      end if
      
      do j = ii,n
         sequence(i) = arr(j)
         call GetClassEnd(sequence,arr,n,k,i+1,j+1,maxes,NQ,amount)
      end do
   end subroutine GetClassEnd
   
end module FCOV_storage

module FCOV_FCFuns
   use Fcov_primitive
   use FCOV_storage

   implicit none
   
   contains
   
   
   function FCOV1_stored_HT_high_generic(v1,v1_pos,v1_n,v2,v2_pos,v2_n,k,k_mapped,td,fc_ss1,fc_ss2,gr_idx,exc_idx,fc_arr)result(res)
      type(ExcState) td
      type(fc_storSys) fc_ss1,fc_ss2
      logical neg1,neg2,check
      integer :: k,k_mapped
      integer(int16) v1_pos(v1_n),v2_pos(v2_n)
      integer j,v1(v1_n),v2(v2_n),v1_n,v2_n
      integer(int64) gr_idx,gr_idx2,exc_idx,exc_idx2
      double precision sum1,sum2,intg1,intg2,res,fc_arr(:,:)
      
      sum1=0d0
      sum2=0d0
      intg1=0d0
      intg2=0d0
      !you can check that k == fc_ss1%v_map(v1_pos(k_mapped))
      intg1=td%B(k)*FC_arr(gr_idx,exc_idx)
      
      if(k_mapped==0)then
         !nothing
      else
         call LowerMode_inplace(v1,k_mapped,1,v1_n,neg1)
         if(.not.neg1)then
            gr_idx2=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
            ! 2*(v1+1) because it is decremented
            intg2=sqrt_arr(2*(v1(k_mapped)+1))*td%A(k,k)*fc_arr(gr_idx2,exc_idx)
         end if
         call incrementMode_inplace(v1,k_mapped,1,v1_n)
      end if
      
      do j = 1,v2_n
         call LowerMode_inplace(v2,j,1,v2_n,neg2)
         if(.not.neg2)then
            exc_idx2=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
            sum2=sum2+sqrt((v2(j)+1)*0.5)*td%E(k,fc_ss2%v_map(v2_pos(j)))*FC_arr(gr_idx,exc_idx2)
         end if
         call IncrementMode_inplace(v2,j,1,v2_n)
      end do
      
      
      
      do j = 1,v1_n
         if(k_mapped==j)cycle
         call LowerMode_inplace(v1,j,1,v1_n,neg2)
         if(.not.(neg2))then
            gr_idx2=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
            sum1=sum1+sqrt_arr(2*(v1(j)+1))*td%A(fc_ss1%v_map(v1_pos(j)),k)*FC_arr(gr_idx2,exc_idx)
         end if
         call IncrementMode_inplace(v1,j,1,v1_n)
      end do
        
      if(k_mapped==0)then
         res=1d0/sqrt_arr(2)*(intg1+sum1+sum2)
      else
         res=1d0/sqrt_arr(2*(v1(k_mapped)+1))*(intg1+intg2+sum1+sum2)
      end if
   end function FCOV1_stored_HT_high_generic
   
   function FCOV2_stored_HT_high_generic(v1,v1_pos,v1_n,v2,v2_pos,v2_n,k,k_mapped,td,fc_ss1,fc_ss2,gr_idx,exc_idx,fc_arr)result(res)
      type(ExcState) td
      type(fc_storSys) fc_ss1,fc_ss2
      logical neg1,neg2,check
      integer :: k,k_mapped
      integer(int16) v1_pos(v1_n),v2_pos(v2_n)
      integer j,v1(v1_n),v2(v2_n),v1_n,v2_n
      integer(int64) gr_idx,gr_idx2,exc_idx,exc_idx2
      double precision sum1,sum2,intg1,intg2,res,fc_arr(:,:)
      
      sum1=0d0
      sum2=0d0
      intg1=0d0
      intg2=0d0
      !you can check that k == fc_ss2%v_map(v2_pos(k_mapped))
      intg1=td%D(k)*FC_arr(gr_idx,exc_idx)
      
      if(k_mapped==0)then
         !nothing
      else
         call LowerMode_inplace(v2,k_mapped,1,v2_n,neg1)
         if(.not.neg1)then
            exc_idx2=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
            ! 2*(v1+1) because it is decremented
            intg2=sqrt_arr(2*(v2(k_mapped)+1))*td%C(k,k)*fc_arr(gr_idx,exc_idx2)
         end if
         call incrementMode_inplace(v2,k_mapped,1,v2_n)
      end if
      
      do j = 1,v1_n
         call LowerMode_inplace(v1,j,1,v1_n,neg2)
         if(.not.neg2)then
            gr_idx2=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
            sum2=sum2+sqrt((v1(j)+1)*0.5)*td%E(fc_ss1%v_map(v1_pos(j)),k)*FC_arr(gr_idx2,exc_idx)
         end if
         call IncrementMode_inplace(v1,j,1,v1_n)
      end do
      
      
      
      do j = 1,v2_n
         if(k_mapped==j)cycle !optimize sometime?
         call LowerMode_inplace(v2,j,1,v2_n,neg2)
         if(.not.(neg2))then
            exc_idx2=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
            sum1=sum1+sqrt_arr(2*(v2(j)+1))*td%C(fc_ss2%v_map(v2_pos(j)),k)*FC_arr(gr_idx,exc_idx2)
         end if
         call IncrementMode_inplace(v2,j,1,v2_n)
      end do
            
      if(k_mapped==0)then
         res=1d0/sqrt_arr(2)*(intg1+sum1+sum2)
      else
         res=1d0/sqrt_arr(2*(v2(k_mapped)+1))*(intg1+intg2+sum1+sum2)
      end if
   end function FCOV2_stored_HT_high_generic
   
   
   
   !constructs an FC integral array with given parameters
   function MakeFCArr_recalculate(FC_ss1,FC_ss2,FC_ss_dusch,fc_00,modes_arr_v2,n_thr,fc_check,doOutput,write_fcarr)result(fc_arr)
      double precision,allocatable :: fc_arr(:,:)
      double precision fc_00
      integer,allocatable :: fc_arr_check(:,:),seq(:)
      integer(int16),allocatable :: modes_arr_gs(:,:)
      integer :: n_thr
      type(FC_storSys) FC_ss1,FC_ss2
      type(FC_storSys_Dusch) FC_ss_dusch
      logical fc_check,doOutput,write_fcarr
      type(modes_arr_t),allocatable,intent(out) :: modes_arr_v2(:)
      type(v_col),allocatable :: v2_arr(:)
      integer(int16),allocatable :: v2_pos(:)

      integer c_ms,ii,iii
      integer(int64) comb_n,i_64
      integer(int64) summ,multi
      integer c_gs

      if(write_fcarr)open(999,file='FC_ARR')
      allocate(fc_arr(FC_ss1%fc_arr_size,FC_ss2%fc_arr_size))
      if(fc_check)then
         fc_arr=HUGE(1d0)
         allocate(fc_arr_check(FC_ss1%fc_arr_size,FC_ss2%fc_arr_size)) !TODO, switch the indices, because the excited state is iterated over a lot more than ground state. Right now this might create a lot of cache misses.
         fc_arr_check=0
         fc_arr_check(1,1)=1
      end if
      !call PopulateFCStorage_Class1_fcss(.true.,FC_ss_dusch,fc_ss1,fc_ss2,fc_check,fc_arr,fc_arr_check)
      !call PopulateFCStorage_Class1_fcss(.false.,FC_ss_dusch,fc_ss1,fc_ss2,fc_check,fc_arr,fc_arr_check)
#ifdef DEBUG
      if(write_fcarr)then
         write(999,'(A,G10.4)')'<0|0>: ',fc_00
         write(999,'(A)')'Gr.state: < 0|'
      end if
#endif
      fc_arr(1,1)=FC_00
      !<0|v>
      allocate(modes_arr_v2(fc_ss2%mc_v),v2_arr(fc_ss2%mc_v))
      !write(output_unit,*)'55555'
      do c_ms=1,fc_ss2%mc_v
         allocate(v2_arr(c_ms)%arr(c_ms),v2_pos(c_ms))
         comb_n=nck(FC_ss2%red_n,c_ms)
         allocate(modes_arr_v2(c_ms)%arr(c_ms,comb_n),seq(c_ms))
         
         modes_arr_v2(c_ms)%arr=LOMC_Wrap(FC_ss2%red_n,c_ms,comb_n)
         if(doOutput)then
            if(c_ms<=9)then
               write(output_unit,'("C",I1)',advance='no')c_ms
            else
               write(output_unit,'("C",I2)',advance='no')c_ms
            end if
         end if
         !Calculate the <0|v> integrals
         !if(c_ms>1)then
         !$OMP PARALLEL DO &
         !$OMP SCHEDULE(DYNAMIC,1) DEFAULT(NONE) &
         !$OMP SHARED(comb_n,modes_arr_v2,c_ms,fc_check,fc_arr,fc_arr_check) &
         !$OMP SHARED(fc_ss2,FC_ss_dusch,write_fcarr) &
         !$OMP PRIVATE(i_64,ii,seq,v2_pos)
         do i_64 = 1,comb_n
            ii=1
            v2_pos=modes_arr_v2(c_ms)%arr(:,i_64)
            call ExcM_v2_FC_0_inner(ii,seq,v2_pos,c_ms,fc_ss2,FC_ss_dusch,fc_arr,fc_arr_check,fc_check,write_fcarr)
         end do
         !$OMP END PARALLEL DO
         !end if
         deallocate(seq,v2_pos)
      end do
      
      !write(output_unit,*)'66666'
      !<v|v>
      if(fc_ss1%mc_v>0)then
      do c_gs=1,fc_ss1%mc_v
         comb_n=nck(fc_ss1%red_n,c_gs)
         allocate(seq(c_gs))
         modes_arr_gs=LOMC_Wrap(fc_ss1%red_n,c_gs,comb_n)
         
         !!$OMP PARALLEL DO &
         !!$OMP DEFAULT(NONE) &
         !!$OMP SHARED(comb_n,modes_arr_gs,c_gs,modes_arr_v2,fc_arr,fc_arr_check,fc_check) &
         !!$OMP SHARED(fc_ss1,fc_ss2,fc_ss_dusch) &
         !!$OMP PRIVATE(i_64,ii,seq,v2_arr)
         do i_64 = 1,comb_n
            ii=1
            call Excm_v1_FC_inner(ii,seq,modes_arr_gs(:,i_64),c_gs,fc_ss1,modes_arr_v2,fc_ss2%mc_v,v2_arr,fc_ss2,fc_arr,fc_arr_check,n_thr,fc_check,doOutput,write_fcarr)
         end do
         !!$OMP END PARALLEL DO
         
         deallocate(modes_arr_gs,seq)
      end do
      end if
      !write(output_unit,*)'77777'
      
      !allocate(modes_arr_v2(fc_ss2%mc_v,))
      if(fc_check)then
         summ=0
         do ii=1,fc_ss1%fc_arr_size
            do iii=1,fc_ss2%fc_arr_size
               summ=summ+fc_arr_check(ii,iii)
            end do
         end do
         multi=FC_ss1%fc_arr_size*FC_ss2%fc_arr_size
         if(summ /= multi)then
            write(output_unit,*)'ERROR on FC_ARR_CHECK'
            call exit(8)
         end if
         deallocate(fc_arr_check)
      end if
#ifdef DEBUG
      if(write_fcarr)then
         close(999)
         write_fcarr=.false.
      end if
#endif   
      contains
      
      recursive subroutine ExcM_v2_FC_0_inner(i,v2,v2_pos,v2_clas,fc_ss2,FC_ss_dusch,fcarr,fcarr_check,check,write_fcarr)
         integer(int16) v2_pos(v2_clas)
         integer i,j,v2_clas,v2(v2_clas)
         integer,allocatable :: fcarr_check(:,:)
         type(FC_storSys) fc_ss2
         type(FC_storSys_Dusch) FC_ss_dusch
         double precision fcarr(:,:),res,res2
         logical check,write_fcarr
         !integer idx1,idx2
         
         if(i>v2_clas)then
            res=FCOV2_stored_fcss(v2,v2_pos,v2_clas,1,fc_ss2,FC_ss_dusch,fcarr,fcarr_check,check)
#ifdef DEBUG
            if(v2_clas==3 .and. check)then
               if(all(v2==1))then
                  res2=FC3_v1_0_v2_1_1_1(fcarr(1,1),fc_ss_dusch%C(v2_pos(2),v2_pos(1)),fc_ss_dusch%C(v2_pos(3),v2_pos(1)),fc_ss_dusch%C(v2_pos(3),v2_pos(2)), &
                  FC_ss_dusch%D(v2_pos(1)),FC_ss_dusch%D(v2_pos(2)),FC_ss_dusch%D(v2_pos(3)))
                  res2=0d0 !just something I can put a breakpoint on
               else if(v2(1)==2 .and. v2(2)==1 .and. v2(3)==1)then
                  res2=FC3_v1_0_v2_2_1_1(fcarr(1,1),fc_ss_dusch%C(v2_pos(1),v2_pos(1)),fc_ss_dusch%C(v2_pos(2),v2_pos(1)),fc_ss_dusch%C(v2_pos(3),v2_pos(1)),fc_ss_dusch%C(v2_pos(3),v2_pos(2)), &
                  FC_ss_dusch%D(v2_pos(1)),FC_ss_dusch%D(v2_pos(2)),FC_ss_dusch%D(v2_pos(3)))
                  res2=0d0
               end if
            end if
            if(write_fcarr)then
               !$OMP CRITICAL
               write(999,'(A,1X,G10.4)')TR(FC2Str_new(v2,MapVPos(v2_pos,v2_clas,fc_ss2%v_map,fc_ss2%red_n),v2_clas,.false.,.false.)),res
               !$OMP END CRITICAL
            end if
#endif
            
            return
         end if
        
         do j = 1,fc_ss2%red_dims(v2_pos(v2_clas-i+1),v2_clas)
            v2(v2_clas-i+1)=j
            call ExcM_v2_FC_0_inner(i+1,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fcarr,fcarr_check,check,write_fcarr)
         end do
      end subroutine ExcM_v2_FC_0_inner
      
      recursive subroutine Excm_v1_FC_inner(i,v1,v1_pos,v1_clas,fc_ss1,modes_arr_v2,mc_v2,v2_arr,fc_ss2,fc_arr,fc_arr_check,n_thr,check,doOutput,write_fcarr)
         integer mc_v2,c_v2,n_thr
         integer i,ii,j,v1_clas,v1(v1_clas)
         type(v_col) v2_arr(mc_v2) !just some memory to pass
         integer, allocatable :: v2(:)
         integer :: v1_idx
         integer,allocatable :: fc_arr_check(:,:)
         type(FC_storSys) fc_ss1,fc_ss2
         type(modes_arr_t) :: modes_arr_v2(mc_v2)
         logical check,doOutput,write_fcarr
         double precision fc_arr(:,:),res
         integer(int16) :: v1_pos(v1_clas)
         integer(int16),allocatable :: v2_pos(:)
         
         
         if(i>v1_clas)then !Set final state mode excitations
#ifdef DEBUG
            if(write_fcarr)write(999,'(A,A)')'Gr.state: ',TR(FC2Str_new(v1,MapVPos(v1_pos,v1_clas,fc_ss1%v_map,fc_ss1%red_n),v1_clas,.true.,.false.))
#endif            
            v1_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_clas,fc_ss1)
            !v1_modes(v1_idx)=MakeMode(v1,v1_pos,v1_clas)
            ! print *,v1
            ! print *,v1_pos
            if(v1_clas>=1)then
               res=FCOV1_stored_fcss(v1,v1_pos,v1_clas,1,fc_ss1,[0],[0_2],0,0,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check) !<v|0>
#ifdef DEBUG
               if(write_fcarr)write(999,'(A,1X,G10.4)')TR(FC2Str_new([0],[0_2],0,.false.,.true.)),res
#endif            
               !fc_arr_checK(v1_idx,1)=fc_arr_checK(v1_idx,1)+1
            end if
            do c_v2 = 1,mc_v2
               if(doOutput)then
                  if(c_v2<=9)then
                     write(output_unit,'("C",I1)',advance='no')c_v2
                  else
                     write(output_unit,'("C",I2)',advance='no')c_v2
                  end if
               end if
               v2=v2_arr(c_v2)%arr
               allocate(v2_pos(c_v2))
               !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC,2) &
               !$OMP SHARED(modes_arr_v2,c_v2,fc_ss_dusch) &
               !$OMP SHARED(v1_pos,v1_clas,fc_ss1,fc_ss2,fc_arr,FC_arr_check,check,write_fcarr) &
               !$OMP FIRSTPRIVATE(v1) &
               !$OMP PRIVATE(v2,ii,j,v2_pos)
               do j = 1,size(modes_arr_v2(c_v2)%arr,dim=2)
                  ii=1
                  v2_pos=modes_arr_v2(c_v2)%arr(:,j)
                  call ExcM_v2_FC_inner(ii,v1,v1_pos,v1_clas,fc_ss1,v2,v2_pos,c_v2,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,write_fcarr)
               end do
               !$OMP END PARALLEL DO
               deallocate(v2_pos)
            end do
            return
         end if
         
         do j = 1,fc_ss1%red_dims(v1_pos(v1_clas-i+1),v1_clas)
            v1(v1_clas-i+1)=j
            call Excm_v1_FC_inner(i+1,v1,v1_pos,v1_clas,fc_ss1,modes_arr_v2,mc_v2,v2_arr,fc_ss2,fc_arr,fc_arr_check,n_thr,check,doOutput,write_fcarr)
         end do
      end subroutine Excm_v1_FC_inner
      
      
      recursive subroutine ExcM_v2_FC_inner(i,v1,v1_pos,v1_clas,fc_ss1,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,write_fcarr)
         integer(int16) v2_pos(v2_clas),v1_pos(v1_clas)
         integer i,j,v2_clas,v2(v2_clas),v1_clas,v1(v1_clas)
         integer,allocatable :: fc_arr_check(:,:)
         logical check,write_fcarr
         type(FC_storSys) fc_ss1,fc_ss2
         type(fc_storSys_dusch) fc_ss_dusch
         double precision fc_arr(:,:),res
         
         if(i>v2_clas)then
            res = FCOV1_stored_fcss(v1,v1_pos,v1_clas,1,fc_ss1,v2,v2_pos,v2_clas,1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
#ifdef DEBUG
            if(write_fcarr)then
               !$OMP CRITICAL
               write(999,'(A,1X,G10.4)')TR(FC2Str_new(v2,MapVPos(v2_pos,v2_clas,fc_ss2%v_map,fc_ss2%red_n),v2_clas,.false.,.true.)),res
               !$OMP END CRITICAL
            end if
#endif
            return
         end if
         
         do j = 1,fc_ss2%red_dims(v2_pos(v2_clas-i+1),v2_clas)
            v2(v2_clas-i+1)=j
            call ExcM_v2_FC_inner(i+1,v1,v1_pos,v1_clas,fc_ss1,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,write_fcarr)
         end do
      end subroutine ExcM_v2_FC_inner
      
      
   end function MakeFCArr_recalculate
   
   function FCOV1_stored_fcss(v1,v1_pos,v1_n,v1_i,fc_ss1,v2,v2_pos,v2_n,v2_i,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)result(res)
      logical neg1,neg2,check
      integer :: v1_i
      integer(int16) v2_pos(v2_n),v1_pos(v1_n)
      integer v2_i,j,v1(v1_n),v2(v2_n),v1_n,v2_n
      integer,allocatable :: fc_arr_check(:,:)
      integer(int64) this1_idx,this2_idx,gr_idx,gr_idx2,exc_idx
      double precision sum1,sum2,intg1,intg2,res,fc_arr(:,:)
      type(FC_storSys) fc_ss1,fc_ss2
      type(FC_storSys_Dusch) fc_ss_dusch
      
      if(v1_n==0)then
         res=FCOV2_stored_fcss(v2,v2_pos,v2_n,v2_i,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
         return
      end if
      sum1=0d0
      sum2=0d0
      !intg1=0d0
      intg2=0d0
      
      this1_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
      this2_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
      CALL LowerMode_inplace(v1,v1_i,1,v1_n,neg1)
      gr_idx2=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
      do j = 1,v2_n
         CALL LowerMode_inplace(v2,j,1,v2_n,neg2)
         if(.not.(neg2))then
            exc_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
            sum2=sum2+sqrt(dble(v2(j)+1)/2d0)*fc_ss_dusch%E(v1_pos(v1_i),v2_pos(j))*FC_arr(gr_idx2,exc_idx)
         end if
         CALL incrementMode_inplace(v2,j,1,v2_n)
      end do

      do j = 1,v1_i-1
         call LowerMode_inplace(v1,j,1,v1_n,neg2)
         if(.not.(neg2))then
            gr_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
            sum1=sum1+sqrt(2d0*(v1(j)+1))*fc_ss_dusch%A(v1_pos(j),v1_pos(v1_i))*FC_arr(gr_idx,this2_idx)
         end if
         call IncrementMode_inplace(v1,j,1,v1_n)
      end do

      do j = v1_i+1,v1_n
         call LowerMode_inplace(v1,j,1,v1_n,neg2)
         if(.not.(neg2))then
            gr_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
            sum1=sum1+sqrt(2d0*(v1(j)+1))*fc_ss_dusch%A(v1_pos(j),v1_pos(v1_i))*FC_arr(gr_idx,this2_idx)
         end if
         call IncrementMode_inplace(v1,j,1,v1_n)
      end do
      
      intg1=fc_ss_dusch%B(v1_pos(v1_i))*FC_arr(gr_idx2,this2_idx)
      
      call LowerMode_inplace(v1,v1_i,1,v1_n,neg1)
      if(.not.neg1)then
         gr_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_n,fc_ss1)
         intg2=sqrt(2d0*(v1(v1_i)-1+2))*fc_ss_dusch%A(v1_pos(v1_i),v1_pos(v1_i))*FC_arr(gr_idx,this2_idx)
      end if
      call IncrementMode_inplace(v1,v1_i,2,v1_n)
      
      res=1d0/sqrt(2d0*v1(v1_i))*(intg1+intg2+sum1+sum2)
      FC_arr(this1_idx,this2_idx)=res
      if(check)FC_arr_check(this1_idx,this2_idx)=FC_arr_check(this1_idx,this2_idx)+1
   end function FCOV1_stored_fcss
   

   function FCOV2_stored_fcss(v2,v2_pos,v2_n,v2_i,fc_ss2,fc_ss_dusch,fc_arr,FC_arr_check,check)result(res)
      logical neg1,neg2,check
      integer(int16) v2_pos(v2_n)
      integer :: v2_i
      integer j,j2,v2_n,v2(v2_n),idx2,N
      integer,allocatable :: fc_arr_check(:,:)
      integer(int64) this2_idx,exc_idx
      type(FC_storSys) fc_ss2
      type(FC_storSys_dusch) fc_ss_dusch
      double precision sum1,intg1,intg2,fc_arr(:,:),res
      
      sum1=0d0
      intg1=0d0
      intg2=0d0
      
      this2_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
      
      !write(99,*)this2_idx
      call LowerMode_inplace(v2,v2_i,1,v2_n,neg1)
      
      do j = 1,v2_i-1
         call LowerMode_inplace(v2,j,1,v2_n,neg2)
         if(.not.(neg2))then
            exc_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
            sum1=sum1+sqrt(2d0*(v2(j)+1))*fc_ss_dusch%C(v2_pos(j),v2_pos(v2_i))*FC_arr(1,exc_idx)
         end if
         call IncrementMode_inplace(v2,j,1,v2_n)
      end do
      
      do j = v2_i+1,v2_n
         call LowerMode_inplace(v2,j,1,v2_n,neg2)
         if(.not.(neg2))then
            exc_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
            sum1=sum1+sqrt(2d0*(v2(j)+1))*fc_ss_dusch%C(v2_pos(j),v2_pos(v2_i))*FC_arr(1,exc_idx)
         end if
         call IncrementMode_inplace(v2,j,1,v2_n)
      end do
      
      exc_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
      intg1=fc_ss_dusch%D(v2_pos(v2_i))*FC_arr(1,exc_idx)
      
      call LowerMode_inplace(v2,v2_i,1,v2_n,neg2)
      if(.not.neg2)then
         exc_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_n,fc_ss2)
         intg2=sqrt(2d0*(v2(v2_i)-1+2))*fc_ss_dusch%C(v2_pos(v2_i),v2_pos(v2_i))*FC_arr(1,exc_idx)
      end if
      call IncrementMode_inplace(v2,v2_i,2,v2_n)
      res=1d0/sqrt(2d0*v2(v2_i))*(intg1+intg2+sum1)
      FC_Arr(1,this2_idx)=res
      if(check)FC_arr_check(1,this2_idx)=FC_arr_check(1,this2_idx)+1
   end function FCOV2_stored_fcss
   
   
      
   
   
end module FCOV_FCFuns

module OMP_functions
   !$ use omp_lib
   implicit none
   contains
   
   function getTime()result(t)
      double precision t
      call CPU_TIME(t)
      !$ t=OMP_get_wtime()
   end function getTime

end module OMP_functions

module spectra
   use iso_fortran_env
   use constants
   implicit none
   
   contains
   
   !w are HO frequencies in cm-1
   function qpar_vib(w,n,temp)result(q)
      integer :: n,i
      double precision :: w(n),q,temp
      
      q=1d0
      do i = 1,n
         if(w(i)<0d0)cycle
         q=q*(1d0/(1d0-exp(-h*cc*w(i)*100/(kb*temp))))
      end do
   end function qpar_vib
   
   
   subroutine ap3(is,ns0,ECM,a,s,wrin,wrax,npx,fwhh,lglg,temp,temp_spr)
   integer*4 ns0,npx,i,is
   real*8 ECM,wrin,wrax,fwhh,w,dw,a(2),s(npx,2,ns0),t,temp,fac
   logical lglg,temp_spr
   dw=(wrax-wrin)/(npx-1)
   w=wrin-dw
   if(temp_spr)then
      fac=1d0/(ECM*cm_2_au*(1d0-exp(-ECM*100*h*cc/(kb*temp))))
   else
      fac=1
   end if
   do i=1,npx
      w=w+dw
      t=sf(fwhh,lglg,w,ECM)
      s(i,1,is)=s(i,1,is)+t*a(1)*fac
      s(i,2,is)=s(i,2,is)+t*a(2)*fac
   end do
   return
   end
   
   elemental function sf(d,g,x,x0)
      logical,intent(in) :: g
      real*8,intent(in) :: d,x,x0
      real*8 sf,dd
      real*8,parameter :: pi=3.14159265358979d0
      real*8,parameter :: spi=1.77245385090552d0

      dd=((x-x0)/d)**2
      if(g)then
       if(dd.lt.32.0d0)then
        sf=exp(-dd)/d/spi
       else
        sf=0.0d0
       endif
      else
       if(dd.lt.1.0d14)then
        sf=1.0d0/d/(dd+1.0d0)/pi
       else
        sf=0.0d0
       endif
      endif
      return
   end function sf
   
end module spectra

module FFT
   use iso_fortran_env
   use constants
   implicit none
   
   contains
      
   subroutine fft_ct_br(a_,br,k)
      integer n,k
      integer s,m,j,kk
      integer br(2**k)
      !the pseudo code expects arrays starting at 0
      double complex a_(0:(2**k-1)),w,wm,u,t
      
      n=2**k
      call ReorderArr_C_inplace(a_,br,n)
      
      do s = 1,k
         m=2**s
         wm=exp(-2*pi*iu/m)
         do kk = 0,n-1,m
            w=1
            do j=0,m/2-1
               t=w*a_(kk+j+m/2)
               u=a_(kk+j)
               a_(kk+j)=u+t
               a_(kk+j+m/2)=u-t
               w=w*wm
            end do
         end do
      end do
      
   end subroutine fft_ct_br
   
   function ReorderArr_C(arr,order,n)result(newarr)
      integer n,order(n),i
      double complex arr(n),newarr(n)
      
      do i = 1,n
         newarr(i)=arr(order(i)+1)
      end do
   end function ReorderArr_C
   
   subroutine ReorderArr_C_inplace(arr,order,n)
      integer n,order(n),i,ii,bufi
      double complex arr(n),buf
      
      i=1
      do while(i<n)
         if(order(i)+1==i)then
            i=i+1
            cycle
         end if
         ii=order(i)+1
         buf=arr(ii)
         arr(ii)=arr(i)
         arr(i)=buf
         
         bufi=order(ii)
         order(ii)=order(i)
         order(i)=bufi
         
      end do
   end subroutine ReorderArr_C_inplace
   
   function BRMake(k)result(order)
      integer i,k
      integer :: n
      integer order(2**k)
      logical(1),allocatable :: order_bits(:,:)
      
      n=2**k
      allocate(order_bits(k,2**k))
      order=GetSeq0(n)
      do i = 1,n
         order_bits(:,i)=Dec2Bin(order(i),k)
         call BitReversal(order_bits(:,i),k)
         order(i)=Bin2Dec(order_bits(:,i),k)
      end do
      deallocate(order_bits)
   end function BRMake
     
   function Bin2Dec(bits,n)result(num)
      integer n,num,i
      logical(1) bits(n)
      
      num=0
      do i = 1,n
         if(bits(n-i+1))num=num+2**(i-1)
      end do
   end function Bin2Dec
   
   function Dec2Bin(num,n)result(bits)
      integer num,n,i,buf
      logical(1) bits(n)
      
      bits=.false.
      
      if(2**n<num)then
         write(output_unit,*)'ERROR: small amount of bits to hold ',num,', got',n
         call exit(2)
      end if
      buf=num
      do i = n-1,0,-1
         if(2**i<=buf)then
            bits(n-i)=.true.
            buf=buf-2**i
         end if
      end do
   end function Dec2Bin
   
   
   subroutine BitReversal(bits,n)
      integer n,i,i_rev
      logical(1) bits(n),buf
      
      do i = 1,n/2,1
         buf=bits(i)
         i_rev=n-i+1
         bits(i)=bits(i_rev)
         bits(i_rev)=buf
      end do
   end subroutine BitReversal
   
   function GetSeq0(n)result(seq)
      integer n,seq(n),i
      do i = 1,n
         seq(i)=i-1
      end do
   end function GetSeq0
end module FFT


module RROA_TD
   use constants
   use iso_fortran_env
   use stuff
   use FCOV_storage
   use FFT
   use strings
#ifdef _OPENACC
   use openacc
#endif
#ifdef _OPENMP
   use omp_lib
#endif
   use ownmath
   implicit none
   
   integer :: z_lwork = -2
   integer :: d_lwork = -2
   integer :: k_inc = 0
   logical :: debug=.false.
   logical :: firstpass=.true.
   double precision :: t_off
   double precision,parameter :: t_off_fac=0.1d0
   integer(int64) :: nq_8=huge(1_8),bigdim=huge(1_8),dp1=huge(1_8),dp2=huge(1_8),dp3=huge(1_8),dp4=huge(1_8)
   !I have no shame
   
   type list_arr_int32
      integer,allocatable :: arr(:)
      integer :: nz=huge(1)
   end type list_arr_int32
   
   type RROA_Options_TD
      !TODO
   end type RROA_Options_TD
   
   type TMS
      double precision :: u(3),m(3),q(3,3)
      double precision,allocatable :: du(:,:),dm(:,:),dq(:,:,:)
      double precision,allocatable :: du2(:,:,:),dm2(:,:,:),dq2(:,:,:,:)
   end type TMS
   
   type TD_SplitP_Mats_T1T1 !Usable for alpha, G and Gc ROA tensors
      double complex :: a_eta(3),b_eta(3)
      double complex,allocatable :: j_zeta_a(:,:),j_zeta_p_b(:,:) !FC-HT,HT-FC
      double complex :: a_zeta_p_b(3,3) !HT-HT
      
      double complex :: eta_a2_eta(3),eta_b2_eta(3),a2_zeta_f(3),b2_zeta_f(3) !FC-HT2,HT2-FC
      double complex,allocatable :: j_zeta_p_b2_eta(:,:),j_zeta_a2_eta(:,:) !FC-HT2,HT2-FC
      
      double complex :: a_zeta_p_b2_eta(3,3),b_zeta_p_a2_eta(3,3)
      double complex,allocatable :: j_zeta_p_b2_zeta_p_a(:,:,:),j_zeta_a2_zeta_p_b(:,:,:) !HT-HT2,HT2-HT
      
      double complex :: zeta_p_a2_zeta_p_b2_f(3,3) !HT2-HT2
      double complex,allocatable :: j_zeta_a2_zeta_p_b2_eta(:,:,:),j_zeta_p_b2_zeta_p_a2_eta(:,:,:)
      double complex :: eta_a2_zeta_p_b2_eta(3,3)
   end type TD_SplitP_Mats_T1T1
   
   type TD_SplitP_Mats_T1T2 !Usable for A, ROA tensor
      double complex :: a_eta(3),b_eta(3,3)
      double complex,allocatable :: j_zeta_p_b(:,:,:),j_zeta_a(:,:) !FC-HT,HT-FC
      double complex :: a_zeta_p_b(3,3,3) !HT-HT
      
      double complex :: eta_a2_eta(3),eta_b2_eta(3,3),a2_zeta_f(3),b2_zeta_f(3,3) !FC-HT2,HT2-FC
      double complex,allocatable :: j_zeta_p_b2_eta(:,:,:),j_zeta_a2_eta(:,:) !FC-HT2,HT2-FC
      
      double complex :: a_zeta_p_b2_eta(3,3,3),b_zeta_p_a2_eta(3,3,3)
      double complex,allocatable :: j_zeta_p_b2_zeta_p_a(:,:,:,:),j_zeta_a2_zeta_p_b(:,:,:,:) !HT-HT2,HT2-HT
      
      double complex :: zeta_p_a2_zeta_p_b2_f(3,3,3) !HT2-HT2
      double complex,allocatable :: j_zeta_a2_zeta_p_b2_eta(:,:,:,:),j_zeta_p_b2_zeta_p_a2_eta(:,:,:,:)
      double complex :: eta_a2_zeta_p_b2_eta(3,3,3)
   end type TD_SplitP_Mats_T1T2
   
   type TD_SplitP_Mats_T2T1 !Ac, ROA tensor
      double complex :: a_eta(3,3),b_eta(3)
      double complex,allocatable :: j_zeta_p_b(:,:),j_zeta_a(:,:,:) !FC-HT,HT-FC
      double complex :: a_zeta_p_b(3,3,3) !HT-HT
      
      double complex :: eta_a2_eta(3,3),eta_b2_eta(3),a2_zeta_f(3,3),b2_zeta_f(3) !FC-HT2,HT2-FC
      double complex,allocatable :: j_zeta_p_b2_eta(:,:),j_zeta_a2_eta(:,:,:) !FC-HT2,HT2-FC
      
      double complex :: a_zeta_p_b2_eta(3,3,3),b_zeta_p_a2_eta(3,3,3)
      double complex,allocatable :: j_zeta_p_b2_zeta_p_a(:,:,:,:),j_zeta_a2_zeta_p_b(:,:,:,:) !HT-HT2,HT2-HT
      
      double complex :: zeta_p_a2_zeta_p_b2_f(3,3,3) !HT2-HT2
      double complex,allocatable :: j_zeta_a2_zeta_p_b2_eta(:,:,:,:),j_zeta_p_b2_zeta_p_a2_eta(:,:,:,:)
      double complex :: eta_a2_zeta_p_b2_eta(3,3,3)
   end type TD_SplitP_Mats_T2T1
   
   contains
   
   function ITG_TR(t_max,dt,points_c,X,wexc)result(res)
      double precision t_max,wexc,wexc2,dt,t
      integer(int64) points_C,i
      double complex X(points_C),res,buf
      
      wexc2=wexc
      buf=(X(1)*exp(iu*wexc2*dt*t_off)+X(points_c)*exp(iu*wexc2*(t_max+dt*t_off)))*0.5
      ! buf=(X(1)*exp(iu*wexc2*dt*t_off)+X(points_c)*exp(iu*wexc2*(t_max+dt*t_off)))
      
      res=0
      !!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
      !!$OMP SHARED(points_c,X,dt,wexc) &
      !!$OMP PRIVATE(i,t) &
      !!$OMP REDUCTION(+:res)
      do i = 2,points_c-1
         t=dt*(i-1_8+t_off)
         res=res+X(i)*exp(iu*wexc2*t)
      end do
      !!$OMP END PARALLEL DO
      res=dt*(buf+res)
   end function ITG_TR
   
   function DoRROA_TD_k(trType,nq,k_batch,batch_n,wg,we,u_gr,m_gr,q_gr,du_gr,dm_gr,dq_gr,u_ex,m_ex,q_ex,du_ex,dm_ex,dq_ex,du2_ex,dm2_ex,dq2_ex, &
      v,gamma_gr,kgk,jgj,J_dusch,K_dusch,J_i,K_i,gamma,theta,eps,w_ad,N_points,t_max,fs,wmax,sparse,td_alt,J_nonzero,nz, &
      ht,ht2,fixphase,tmexp,st,contr,n_thr,wexc,nexc,bro,X_tol,write_corrf,write_fft_cf,num_integ,norm_fft,correctPhaseX,correctPhaseX_abs,interpolateFFT,wexc_adapt,w_ad_zero)result(polars)
      
      integer nq,kk,N_points,n_thr,nexc,batch_n,i,j
      integer(int64) wexc_idx(nexc),ii,points_c
      integer k_batch(batch_n),kk_idx
      integer a,b,c,iexc,w_idx,X_tol,nz,trType
      type(list_arr_int32) J_nonzero(batch_n)
      type(big_double) detG
      integer :: bro(2_8**N_points) !bit reversal order for fft
      integer,allocatable :: bro_copy(:)
      type(Polar_exc),target :: polars(batch_n)
      type(Polar_exc),pointer :: polarr
      double precision wg(nq),we(nq),wexc(nexc),wr_exc(nexc),maxw,eps
      double precision u_gr(3),m_gr(3),q_gr(3,3),du_gr(3,nq),dm_gr(3,nq),dq_gr(3,3,nq)
      double precision u_ex(3),m_ex(3),q_ex(3,3),du_ex(3,nq),dm_ex(3,nq),dq_ex(3,3,nq)
      double precision,allocatable :: du2_ex(:,:,:),dm2_ex(:,:,:),dq2_ex(:,:,:,:)
      
      double complex, allocatable :: test_mat(:,:)
      double complex integ,bufC(nexc)
      double precision v(nq),gamma_gr(nq,nq),kgk,jgj(nq,nq),w_ad
      double precision J_dusch(nq,nq),K_dusch(nq)
      double precision J_i(nq,nq),K_i(nq)
      double precision t_max,fs,wmax,dw,gamma,theta
      double precision X_max,X_last,dt,w_ad_new
      double complex,allocatable :: X_0(:)
      double complex, allocatable :: X_ap(:,:,:,:),X_G(:,:,:,:),X_Gc(:,:,:,:),X_A(:,:,:,:,:),X_Ac(:,:,:,:,:)
      double complex,allocatable :: X(:)
      double precision,allocatable :: X_copy(:),X_copy2(:)
      logical ht,ht2,tmexp(2),sparse,write_actual,write_now,write_corrf,write_fft_cf,num_integ,norm_fft,st,correctPhaseX,correctPhaseX_abs
      logical :: interpolateFFT,st_actual,wexc_adapt,w_ad_correct,contr(9),td_alt,w_ad_zero,fixphase
      integer, parameter :: X_last_c=500
      
      type(tms) tmss
      
      character(80) file_zz
      
      points_c=2_8**N_points
      dw=2*pi/t_max
      allocate(test_mat(nq,nq))
      z_lwork=INV_C_lapack_getLWork(test_mat,nq)
      deallocate(test_mat)
     
#ifdef DEBUG
      debug=.true.
#endif
      nq_8=int(nq,kind=int64)
      bigdim=nq_8**5
      dp1=nq_8
      dp2=nq_8**2
      dp3=nq_8**3
      dp4=nq_8**4
      bro_copy=bro
      maxw=dw*(n_points/2-1)
      dt=t_max/dble(points_c)
      t_off=t_off_fac !GLOBAL
      write_actual=write_corrf .and. trType==0
      st_actual=st .and. trType==0 !second term only for fundamentals right now
      
      do i = 1,batch_n
         call AllocatePolar(polars(i),nexc)
         call Polar_new_nexc(polars(i),nexc)
      end do
      ! do iexc=1,nexc
         ! wexc_idx(iexc)=points_c-NINT(wexc(iexc)/(dw),kind=8)+1_8
      ! end do
      flush(output_unit)
      detG=DiagProd_big_vec(wg,nq)
      !if(w_ad_correct)then
      w_ad_new=w_ad-sum(we)/2d0
         !w_ad_new=0d0
      !else
      !   w_ad_new=w_ad
      !end if
      if(w_ad_zero)w_ad_new=0d0
      
      tmss%u=u_ex
      tmss%m=m_ex
      tmss%q=q_ex
      if(ht)then
         if(td_alt)then
            tmss%du=transpose(du_ex)
            tmss%dm=transpose(dm_ex)
            allocate(tmss%dq(nq,3,3))
            do i = 1,nq
               do a = 1,3
                  do b = 1,3
                     tmss%dq(i,a,b)=dq_ex(a,b,i)
                  end do
               end do
            end do
         else
            tmss%du=du_ex
            tmss%dm=dm_ex
            tmss%dq=dq_ex
         end if
      end if
      
      if(ht2)then
         if(td_alt)then
            allocate(tmss%dq2(nq,nq,3,3),tmss%du2(nq,nq,3),tmss%dm2(nq,nq,3))
            do i = 1,nq
               do j = 1,nq
                  do a = 1,3
                     tmss%du2(i,j,a)=du2_ex(a,i,j)
                     tmss%dm2(i,j,a)=dm2_ex(a,i,j)
                     do b = 1,3
                        tmss%dq2(i,j,a,b)=dq2_ex(a,b,i,j)
                     end do
                  end do
               end do
            end do
         else
            tmss%du2=du2_ex
            tmss%dm2=dm2_ex
            tmss%dq2=dq2_ex
         end if
      end if
      if(td_alt)then
      call Make_Corrf_SplitPropagator(nq,wg,we,tmss,v,detg,kgk,jgj,J_dusch,K_dusch,gamma,theta,eps,w_ad_new,N_points,t_max,ht,ht2,fixphase,contr,n_thr,X_ap,X_g,x_gc,x_a,x_ac,X_0)
      else
      call Make_corrf_k(trType,nq,k_batch,batch_n,wg,we,tmss, &
         v,gamma_gr,detG,kgk,J_dusch,K_dusch,gamma,theta,w_ad_new,N_points,t_max,sparse,J_nonzero,nz,ht,ht2,tmexp,fixphase,contr,n_thr,10d0**(-X_tol), &
         X_ap,X_g,x_gc,x_a,x_ac,X_0)
      end if
      if(firstPass)then
         open(777,file='FILE.X0')
         do ii = 1,2**N_points
            write(777,*)realpart(X_0(ii)),imagpart(X_0(ii))
         end do
         close(777)
         if(write_fft_cf)then
            bro_copy=bro
            call fft_ct_br(X_0,bro_copy,N_points)
            X_0=X_0/dble(points_c)*t_max
            call ShiftFFT(X_0,points_c,t_max,dt,t_off)
            open(777,file='FILE.X0_FFT')
            do ii = 1,2**N_points
               write(777,*)realpart(X_0(ii)),imagpart(X_0(ii))
            end do
            close(777)
         end if
         firstpass=.false.
      end if
      allocate(X(2_8**N_points),X_copy(2_8**N_points),X_copy2(min(X_last_c+1,points_c)))
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP FIRSTPRIVATE(X,X_copy,X_copy2) &
      !$OMP PRIVATE(kk_idx,kk,wr_exc,polarr,bufC,x_max,x_last,write_now) &
      !$OMP SHARED(X_ap,X_G,X_Gc,X_A,X_Ac,dt) &
      !$OMP SHARED(batch_n,X_tol,polars,k_batch,wg,wexc,points_c,dw,n_points,write_actual) &
      !$OMP SHARED(norm_fft,num_integ,interpolatefft,st_actual,bro,nexc,t_max,t_off,trType,wexc_adapt,td_alt)
      do kk_idx=1,batch_n
         kk=k_batch(kk_idx) 
         wr_exc=wexc-wg(kk)
         if(TD_alt)then
            X=X_ap(kk_idx,2,2,:)
         else
            X=X_ap(2,2,kk_idx,:)
         end if
         X_copy=abs(X)
         if(sum(X_copy)>=1d-14)then
            ! X_max=0
            ! do i = 1,points_c
               ! if(X_copy(i)>X_max)then
                  ! X_max_idx=i
               ! end if
            ! end do
            X_max=maxval(X_copy)
            X_last=maxval(X_copy(points_c-NINT(points_c*0.2):points_c))
            !deallocate(X_copy,X_copy2)
            if(X_last/X_max>10d0**(-X_tol))then
               !$OMP CRITICAL
               write(output_unit,*)'Warning: X not converged for mode: ',kk,trType
               write(output_unit,*)'Ratio last/max = ',X_last/X_max
               !$OMP END CRITICAL
            end if
         end if
         polarr=>polars(kk_idx)
         do a =1,3  
            do b =1,3
               if(a==2 .and. b==2)then
                  write_now=write_actual
               else
                  write_now=.false.
               end if
               
               call GetXabk(2_8**N_points,X,X_ap,a,b,kk_idx,td_alt)
               bufC=ExtractPol(X,wg(kk),kk,1,bro,norm_fft,num_integ,wexc,dw,nexc,t_max,dt,t_off,N_points,points_c,interpolateFFT,write_now,.false.,wexc_adapt)
               polarr%Ap(a,b,:)=bufC
               
               call GetXabk(2_8**N_points,X,X_G,a,b,kk_idx,td_alt)
               bufC=ExtractPol(X,wg(kk),kk,2,bro,norm_fft,num_integ,wexc,dw,nexc,t_max,dt,t_off,N_points,points_c,interpolateFFT,.false.,.false.,wexc_adapt)
               polarr%G(a,b,:)=bufC
               
               call GetXabk(2_8**N_points,X,X_Gc,a,b,kk_idx,td_alt)
               bufC=ExtractPol(X,wg(kk),kk,3,bro,norm_fft,num_integ,wexc,dw,nexc,t_max,dt,t_off,N_points,points_c,interpolateFFT,.false.,.false.,wexc_adapt)
               polarr%Gc(a,b,:)=bufC
               
               
               do c = 1,3
                  call GetXabkQ(2_8**N_points,X,X_A,a,b,c,kk_idx,td_alt)
                  bufC=ExtractPol(X,wg(kk),kk,4,bro,norm_fft,num_integ,wexc,dw,nexc,t_max,dt,t_off,N_points,points_c,interpolateFFT,.false.,.false.,wexc_adapt)
                  polarr%A(a,b,c,:)=bufC
                  
                  call GetXabkQ(2_8**N_points,X,X_Ac,a,b,c,kk_idx,td_alt)
                  bufC=ExtractPol(X,wg(kk),kk,5,bro,norm_fft,num_integ,wexc,dw,nexc,t_max,dt,t_off,N_points,points_c,interpolateFFT,.false.,.false.,wexc_adapt)
                  polarr%Ac(a,b,c,:)=bufC
                  
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO
      deallocate(X_ap,X_g,x_gc,x_a,x_ac,x_0)
      deallocate(X,X_copy,X_copy2,bro_copy)
   end function DoRROA_TD_k
   
   subroutine GetXabk(n,Xcur,X,a,b,k,tdalt)
      integer(int64) n
      integer a,b,k
      double complex :: X(:,:,:,:)
      double complex,allocatable :: Xcur(:)
      logical tdalt
      
      if(.not.allocated(Xcur))allocate(Xcur(n))
      if(tdalt)then
         Xcur=X(k,a,b,:) !splicing the wrong dimension but oh well
      else
         Xcur=X(a,b,k,:)
      end if
      
   end subroutine GetXabk
   
   subroutine GetXabkQ(n,Xcur,X,a,b,c,k,tdalt)
      integer(int64) n
      integer a,b,k,c
      double complex :: X(:,:,:,:,:)
      double complex,allocatable :: Xcur(:)
      logical tdalt
      
      if(.not.allocated(Xcur))allocate(Xcur(n))
      if(tdalt)then
         Xcur=X(k,a,b,c,:)
      else
         Xcur=X(a,b,c,k,:)
      end if
      
   end subroutine GetXabkQ
   
   subroutine ShiftFFT(X,n,tmax,dt,t_off)
      integer(int64) n,ii
      double precision tmax,t_off,cur_w,dw,dt
      double complex X(n),shift_factor
      
      dw=2*pi/tmax
      !dt=tmax/n
      do ii = 1, n !thank you Gemini
          ! Determine the frequency for this index (0-based math for physics)
          ! For the first half (positive freqs): w = (k-1) * dw
          ! For the second half (negative freqs): w = (k-1-N) * dw
          
          ! Assuming you only care about the Positive Spectrum (indices 1 to N/2 + 1)
          ! which is typical for vibronic spectra:
          cur_w = dble(ii-1) * dw 
          
          ! Calculate the Phase Shift Factor
          ! This rotates the result from the "FFT time frame" (starting at 0)
          ! to the "Physical time frame" (starting at t_off)
          shift_factor = exp(-iu * cur_w * t_off*dt)
          
          ! Apply correction
          X(ii) = X(ii) * shift_factor
      end do   
   end subroutine ShiftFFT
   
   function ExtractPol(X,wg,kk,which,bro,norm_fft,num_integ,wexc,dw,nexc,t_max,dt,t_off,N_points,points_c,interpolateFFT,output,st,wexc_adapt)result(res)
      logical output,norm_fft,num_integ,interpolateFFT,st,wexc_adapt
      integer(int64) points_c,ii
      integer bro(points_c),N_points
      integer,allocatable :: bro_copy(:)
      double precision wexc(nexc),wexc_c(nexc),t_max,dw,diff,dw_conv,wg,dt,t_off
      double complex X(2**N_points),res(nexc),res1,res2
      integer which,iexc,nexc,w_idxs(nexc),w_idx,w_idx2,kk,bufi
      character(80) file_zz
      character(:),allocatable :: numm
      
      if(wexc_adapt)then
         wexc_c=wexc-wg/2d0
      else
         wexc_c=wexc
      end if
      bufi=floor(log10(dble(kk)))+1
      numm=i2STR(kk,bufi)
      select case(which)
         case(1) 
            file_zz=TR('_ap_YY_'//numm)
         case(2)
            file_zz=TR('_G_YY_'//numm)
         case(3)
            file_zz=TR('_Gc_YY_'//numm)
         case(4)
            file_zz=TR('_A_YY_'//numm)
         case(5)
            file_zz=TR('_Ac_YY_'//numm)
      end select
      if(output)then
         !$OMP CRITICAL
         if(st)then
            open(777,file='X2'//file_zz)
         else
            open(777,file='X'//file_zz)
         end if
         do ii = 1,points_c
            write(777,*)realpart(X(ii)),imagpart(X(ii))
         end do
         close(777)
         !$OMP END CRITICAL
      end if
      
      if(.not.num_integ)then
         bro_copy=bro
         call fft_ct_br(X,bro_copy,N_points)
         call ShiftFFT(X,points_c,t_max,dt,t_off)
         X=X/dble(points_c)*t_max
      end if
      ! dw_conv=dw*2*pi
      dw_conv=dw
      do iexc=1,nexc
         if(num_integ)then
            res(iexc)=iu*(ITG_TR(t_max,dt,points_c,X,wexc_c(iexc)))
         else
            if(interpolateFFT)then
               w_idx=points_c-floor(wexc_c(iexc)/(dw_conv))
               w_idx2=points_c-ceiling(wexc_c(iexc)/(dw_conv))
               diff=(wexc_c(iexc)-(floor(wexc_c(iexc)/dw_conv))*dw_conv)/dw_conv
               res(iexc)=iu*Lerp_C((floor(wexc_c(iexc)/dw_conv))*dw_conv,X(w_idx),(ceiling(wexc_c(iexc)/dw_conv))*dw_conv,X(w_idx2),diff)
            else
               w_idx=points_c-NINT(wexc_c(iexc)/dw_conv)
               res(iexc)=iu*X(w_idx)
            end if
         end if
      end do
      
      if(output)then
         if(num_integ)then
            bro_copy=bro
            call fft_ct_br(X,bro_copy,N_points)
            X=X/dble(points_c)*t_max
            call ShiftFFT(X,points_c,t_max,dt,t_off)
         end if
         !$OMP CRITICAL
         if(st)then
            open(777,file='W2'//file_zz)
         else
            open(777,file='W'//file_zz)
         end if
         do ii = 1,points_c
            write(777,*)realpart(X(ii)),imagpart(X(ii))
         end do
         close(777)
         !$OMP END CRITICAL
      end if
      if(.not.num_integ)deallocate(bro_copy)
   end function ExtractPol
      
   function Lerp_C(x1,y1,x2,y2,alpha)result(res)
      double complex y1,y2,res
      double precision x1,x2,alpha
      double complex k,q
      
      k=(y1-y2)/(x1-x2)
      q=y2-k*x2
      res=k*(x1+(x2-x1)*alpha)+q
   end function Lerp_C
   
   subroutine PrintCDiagonal(mat,n)
      integer n,i,start,endd,ii
      integer,parameter :: step=5
      double complex mat(n,n)
      
      i=1
      do while(i<n)
         start=(i-1)*step+1
         endd=start+step-1
         do ii = start,endd
            write(output_unit,'(E15.8,1X,E15.8,"*i",1X)',advance='no')mat(ii,ii)
         end do
         write(output_unit,*)
      end do
   end subroutine PrintCDiagonal
   
   subroutine Make_Corrf_SplitPropagator(nq,wg,we,tmss,v,detG,kgk,JGJ,J,K,gamma,theta,eps,w_ad,N_points,t_max,ht,ht2,fixphase,contrs,n_thr,X_ap,X_G,X_Gc,X_A,X_Ac,X_0)
      integer nq,N_points,n_thr,i_thr,a,b,c
      integer(int64) points_c,chunk,i,c_arr(9)
      logical ht,ht2,fixphase,contrs(9)
      type(tms) tmss
      type(big_double) detG
      type(TD_SplitP_Mats_T1T1) :: mats_ap,mats_G,mats_Gc
      type(TD_SplitP_Mats_T1T2) mats_A
      type(TD_SplitP_Mats_T2T1) mats_Ac
      double precision J(nq,nq),JGJ(nq,nq),K(nq),kgk,t_max,wg(nq),we(nq),v(nq),damp,damp2,w_ad,t
      double precision gamma,theta,cur_phase,cur_sizee,dt,phase,eps
      double precision,allocatable :: JT(:,:),phases(:),C_G(:)
      double complex :: shift,tau,x0
      double complex, allocatable :: eta(:),zeta(:,:),zeta_p(:,:)
      double complex, allocatable :: Bmat(:,:),Fi(:,:),D(:,:),aa(:),bb(:),abia(:)
      double complex, allocatable :: x_fc_fc(:),x_fc_ht(:),x_ht_fc(:),x_ht_ht(:)
      double complex, allocatable :: x_fc_ht2(:),x_ht2_fc(:),x_ht_ht2(:),x_ht2_ht(:),x_ht2_ht2(:)
      
      double complex :: dip_eta(3),mag_eta(3),quad_eta(3,3)
      
      double complex,allocatable :: Jt_eta(:)
      double complex,allocatable :: j_zeta_p(:,:),j_zeta(:,:)
      
      double complex,allocatable,intent(out) :: X_ap(:,:,:,:),X_G(:,:,:,:),X_Gc(:,:,:,:),X_A(:,:,:,:,:),X_Ac(:,:,:,:,:),X_0(:)
      
      integer ii,kk
      
      points_c=2_8**N_points
      allocate(X_ap(nq,3,3,points_c),X_G(nq,3,3,points_c),X_Gc(nq,3,3,points_c))
      allocate(X_A(nq,3,3,3,points_c),X_Ac(nq,3,3,3,points_c))
      
      if(ht)write(output_unit,'(A)',advance='no')'HT'
      if(ht2)write(output_unit,'(A)',advance='no')'HT2'
      flush(output_unit)
      
      chunk=points_c/n_thr
      dt=t_max/dble(points_c-1_8) !Unsure about the last element of x0, similar problem with x of the spectrum, todo check
      allocate(C_G(nq))
      do i=1,nq
         C_G(i)=2*sqrt(wg(i)*0.5d0)
      end do
      c_arr=0
      do i =1,9
         if(contrs(i))c_arr(i)=1
      end do
      write(output_unit,*)
      write(output_unit,'(9(1X,I1))')c_arr
      allocate(X_0(points_c))
      allocate(phases(points_c))
      allocate(JT(nq,nq))
      JT=transpose(J)
      i_thr=1
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(X_ap,X_G,X_Gc,X_A,X_Ac,X_0) &
      !$OMP SHARED(gamma,theta,J,JT,K,jgj,kgk,tmss,v,wg,we,nq) &
      !$OMP SHARED(w_ad,points_c,t_max,ht,ht2,n_thr,t_off,dt) &
      !$OMP SHARED(chunk,detG,phases,c_g,fixphase,eps,c_arr) &
      !$OMP PRIVATE(i_thr,i,t,aa,bb,abia,Bmat,D,Fi) &
      !$OMP PRIVATE(eta,zeta,zeta_p,phase,a,b,c) &
      !$OMP PRIVATE(j_zeta_p,j_zeta) &
      !$OMP PRIVATE(x_fc_fc,x0,tau,shift,cur_sizee,cur_phase,damp) &
      !$OMP PRIVATE(X_fc_ht,X_ht_fc,X_ht_ht,jt_eta) &
      !$OMP PRIVATE(x_fc_ht2,x_ht2_fc,x_ht_ht2,x_ht2_ht,x_ht2_ht2) &
      !$OMP PRIVATE(mats_ap,mats_G,mats_Gc,mats_A,mats_Ac)
      !$ i_thr=omp_get_thread_num()+1
      
      allocate(eta(nq),zeta(nq,nq),zeta_p(nq,nq))
      allocate(Fi(nq,nq),D(nq,nq))
      allocate(aa(nq),bb(nq),abia(nq),Bmat(nq,nq))
      
      allocate(x_fc_fc(nq))
      allocate(jt_eta(nq))
      
      if(ht)then
         allocate(x_fc_ht(nq),x_ht_fc(nq),x_ht_ht(nq))
         allocate(j_zeta(nq,nq),j_zeta_p(nq,nq))
      end if
      if(ht2)then
         allocate(x_fc_ht2(nq),x_ht2_fc(nq),x_ht_ht2(nq),x_ht2_ht(nq),x_ht2_ht2(nq))
      end if
      call AllocateMats_T1T1(nq,mats_ap,ht,ht2)
      call AllocateMats_T1T1(nq,mats_G,ht,ht2)
      call AllocateMats_T1T1(nq,mats_Gc,ht,ht2)
      call AllocateMats_T1T2(nq,mats_A,ht,ht2)
      call AllocateMats_T2T1(nq,mats_Ac,ht,ht2)
      
      
      if(fixphase)then
         !$OMP MASTER
         write(output_unit,*)'Will artifically fix X0 phase'
         !$OMP END MASTER
         phase=0
         !$OMP DO SCHEDULE(STATIC,chunk)
         do i = 1,points_c
            t=dt*(i-1_8+t_off) !workaround against t=0 (division by zero)
            call x0_ab_SplitPropagator(nq,J,JGJ,wg,we,t,aa,bb,Bmat,abia,eps)
            x0=make_x0_SplitPropagator(nq,detG,v,kgk,aa,bb,Bmat,abia,D,Fi,eta,zeta,zeta_p,phase) !Here C_i and D_i are inverted
            X_0(i)=x0
         end do
         !$OMP END DO
      
         !$OMP MASTER
         call FixPhase_X0(X_0,points_c,phases)
         !$OMP END MASTER
      end if
      phase=0
      !$OMP DO SCHEDULE(STATIC,chunk)
      do i = 1,points_c
         t=dt*(i-1_8+t_off) !workaround against t=0 (division by zero)
         damp=exp(-gamma*t - theta**2*t**2*0.5d0) 
         shift=exp(-iu*t*w_ad)
         
         call x0_ab_SplitPropagator(nq,J,JGJ,wg,we,t,aa,bb,Bmat,abia,eps) !Here C_i and D_i are not inverted yet
         x0=make_x0_SplitPropagator(nq,detG,v,kgk,aa,bb,Bmat,abia,D,Fi,eta,zeta,zeta_p,phase) !Here C_i and D_i are inverted
         if(fixphase)then
            cur_sizee=sqrt(realpart(x0)**2+imagpart(x0)**2)
            cur_phase=phases(i)
            x0=cur_sizee*(cos(cur_phase)+iu*sin(cur_phase))
         end if
         x0=x0*shift*damp
         X_0(i)=x0
         Jt_eta=matmul(J,eta)
            
         call MakeMatrices(ht,ht2,nq,J,eta,zeta,zeta_p,tmss,mats_ap,mats_G,mats_Gc,mats_A,mats_Ac,j_zeta,j_zeta_p)
         
         x_fc_fc=(K+Jt_eta)
         
         do b = 1,3
            do a = 1,3
               X_ap(:,a,b,i)=tmss%u(a)*tmss%u(b)*x_fc_fc*c_arr(1)
               X_G(:,a,b,i)=tmss%u(a)*tmss%m(b)*x_fc_fc*c_arr(1)
               X_Gc(:,a,b,i)=-tmss%m(a)*tmss%u(b)*x_fc_fc*c_arr(1)
               do c = 1,3
                  X_a(:,a,b,c,i)=tmss%u(a)*tmss%q(b,c)*x_fc_fc*c_arr(1)
                  X_ac(:,c,a,b,i)=tmss%q(a,b)*tmss%u(c)*x_fc_fc*c_arr(1)
                  if(ht)then
                     x_fc_ht=Calc_FC_HT(nq,X_fc_fc,tmss%u(a),mats_a%b_eta(b,c),mats_a%J_zeta_p_b(:,b,c))
                     x_ht_fc=Calc_HT_FC(nq,X_fc_fc,tmss%q(b,c),mats_a%a_eta(a),mats_a%J_zeta_a(:,a))
                     x_ht_ht=Calc_HT_HT(nq,X_fc_fc,mats_a%a_eta(a),mats_a%b_eta(b,c),mats_a%a_zeta_p_b(a,b,c),mats_a%j_zeta_a(:,a),mats_a%j_zeta_p_b(:,b,c))
                     X_A(:,a,b,c,i)=X_A(:,a,b,c,i)+x_fc_ht*c_arr(2)+x_ht_fc*c_arr(3)+x_ht_ht*c_arr(4)
                     
                     x_fc_ht=Calc_FC_HT(nq,X_fc_fc,tmss%q(a,b),mats_Ac%b_eta(c),mats_Ac%J_zeta_p_b(:,c))
                     x_ht_fc=Calc_HT_FC(nq,X_fc_fc,tmss%u(c),mats_Ac%a_eta(a,b),mats_Ac%J_zeta_a(:,a,b))
                     x_ht_ht=Calc_HT_HT(nq,X_fc_fc,mats_Ac%a_eta(a,b),mats_Ac%b_eta(c),mats_Ac%a_zeta_p_b(a,b,c),mats_Ac%j_zeta_a(:,a,b),mats_Ac%j_zeta_p_b(:,c))
                     X_Ac(:,c,a,b,i)=X_Ac(:,c,a,b,i)+x_fc_ht*c_arr(2)+x_ht_fc*c_arr(3)+x_ht_ht*c_arr(4)
                  end if
                  if(ht2)then
                     x_fc_ht2=Calc_FC_HT2(nq,X_fc_fc,tmss%u(a),mats_A%eta_b2_eta(b,c),mats_A%b2_zeta_f(b,c),mats_A%j_zeta_p_b2_eta(:,b,c))
                     x_ht2_fc=Calc_HT2_FC(nq,X_fc_fc,tmss%q(b,c),mats_A%eta_a2_eta(a),mats_A%a2_zeta_f(a),mats_A%j_zeta_a2_eta(:,a))
                     x_ht_ht2=Calc_HT_HT2(nq,X_fc_fc,mats_A%a_eta(a),mats_A%eta_b2_eta(b,c),mats_A%b2_zeta_f(b,c),mats_A%a_zeta_p_b2_eta(a,b,c),mats_A%j_zeta_a(:,a),mats_A%j_zeta_p_b2_eta(:,b,c),mats_A%j_zeta_p_b2_zeta_p_a(:,a,b,c))
                     x_ht2_ht=Calc_HT2_HT(nq,X_fc_fc,mats_A%b_eta(b,c),mats_A%eta_a2_eta(a),mats_A%a2_zeta_f(a),mats_A%b_zeta_p_a2_eta(a,b,c),mats_A%j_zeta_p_b(:,b,c),mats_A%j_zeta_a2_eta(:,a),mats_A%j_zeta_a2_zeta_p_b(:,a,b,c))
                     x_ht2_ht2=Calc_HT2_HT2(nq,X_fc_fc,mats_A%j_zeta_a2_eta(:,a),mats_A%b2_zeta_f(b,c),mats_A%eta_b2_eta(b,c),mats_A%j_zeta_a2_zeta_p_b2_eta(:,a,b,c),mats_A%j_zeta_p_b2_eta(:,b,c),mats_A%a2_zeta_f(a),mats_A%eta_a2_eta(a),mats_A%j_zeta_p_b2_zeta_p_a2_eta(:,a,b,c),mats_A%eta_a2_zeta_p_b2_eta(a,b,c),mats_A%zeta_p_a2_zeta_p_b2_f(a,b,c))
                     X_A(:,a,b,c,i)=X_A(:,a,b,c,i)+x_fc_ht2*c_arr(5)+x_ht2_fc*c_arr(6)+x_ht_ht2*c_arr(7)+x_ht2_ht*c_arr(8)+x_ht2_ht2*c_arr(9)
                     
                     x_fc_ht2=Calc_FC_HT2(nq,X_fc_fc,tmss%q(a,b),mats_Ac%eta_b2_eta(c),mats_Ac%b2_zeta_f(c),mats_Ac%j_zeta_p_b2_eta(:,c))
                     x_ht2_fc=Calc_HT2_FC(nq,X_fc_fc,tmss%u(c),mats_Ac%eta_a2_eta(a,b),mats_Ac%a2_zeta_f(a,b),mats_Ac%j_zeta_a2_eta(:,a,b))
                     x_ht_ht2=Calc_HT_HT2(nq,X_fc_fc,mats_Ac%a_eta(a,b),mats_Ac%eta_b2_eta(c),mats_Ac%b2_zeta_f(c),mats_Ac%a_zeta_p_b2_eta(a,b,c),mats_Ac%j_zeta_a(:,a,b),mats_Ac%j_zeta_p_b2_eta(:,c),mats_Ac%j_zeta_p_b2_zeta_p_a(:,a,b,c))
                     x_ht2_ht=Calc_HT2_HT(nq,X_fc_fc,mats_Ac%b_eta(c),mats_Ac%eta_a2_eta(a,b),mats_Ac%a2_zeta_f(a,b),mats_Ac%b_zeta_p_a2_eta(a,b,c),mats_Ac%j_zeta_p_b(:,c),mats_Ac%j_zeta_a2_eta(:,a,b),mats_Ac%j_zeta_a2_zeta_p_b(:,a,b,c))
                     x_ht2_ht2=Calc_HT2_HT2(nq,X_fc_fc,mats_Ac%j_zeta_a2_eta(:,a,b),mats_Ac%b2_zeta_f(c),mats_Ac%eta_b2_eta(c),mats_Ac%j_zeta_a2_zeta_p_b2_eta(:,a,b,c),mats_Ac%j_zeta_p_b2_eta(:,c),mats_Ac%a2_zeta_f(a,b),mats_Ac%eta_a2_eta(a,b),mats_Ac%j_zeta_p_b2_zeta_p_a2_eta(:,a,b,c),mats_Ac%eta_a2_zeta_p_b2_eta(a,b,c),mats_Ac%zeta_p_a2_zeta_p_b2_f(a,b,c))
                     X_Ac(:,c,a,b,i)=X_AC(:,c,a,b,i)+x_fc_ht2*c_arr(5)+x_ht2_fc*c_arr(6)+x_ht_ht2*c_arr(7)+x_ht2_ht*c_arr(8)+x_ht2_ht2*c_arr(9)
                  end if
               end do
               if(ht)then
                  !alpha
                  x_fc_ht=Calc_FC_HT(nq,X_fc_fc,tmss%u(a),mats_ap%b_eta(b),mats_ap%J_zeta_p_b(:,b))
                  x_ht_fc=Calc_HT_FC(nq,X_fc_fc,tmss%u(b),mats_ap%a_eta(a),mats_ap%J_zeta_a(:,a))
                  x_ht_ht=Calc_HT_HT(nq,X_fc_fc,mats_ap%a_eta(a),mats_ap%b_eta(b),mats_ap%a_zeta_p_b(a,b),mats_ap%j_zeta_a(:,a),mats_ap%j_zeta_p_b(:,b))
                  X_ap(:,a,b,i)=X_ap(:,a,b,i)+x_fc_ht*c_arr(2)+x_ht_fc*c_arr(3)+x_ht_ht*c_arr(4)
                  
                  x_fc_ht=Calc_FC_HT(nq,X_fc_fc,tmss%u(a),mats_G%b_eta(b),mats_G%J_zeta_p_b(:,b))
                  x_ht_fc=Calc_HT_FC(nq,X_fc_fc,tmss%m(b),mats_G%a_eta(a),mats_G%J_zeta_a(:,a))
                  x_ht_ht=Calc_HT_HT(nq,X_fc_fc,mats_G%a_eta(a),mats_G%b_eta(b),mats_G%a_zeta_p_b(a,b),mats_G%j_zeta_a(:,a),mats_G%j_zeta_p_b(:,b))
                  X_G(:,a,b,i)=X_G(:,a,b,i)+x_fc_ht*c_arr(2)+x_ht_fc*c_arr(3)+x_ht_ht*c_arr(4)
                  
                  x_fc_ht=Calc_FC_HT(nq,X_fc_fc,-tmss%m(a),mats_Gc%b_eta(b),mats_Gc%J_zeta_p_b(:,b))
                  x_ht_fc=Calc_HT_FC(nq,X_fc_fc,tmss%u(b),mats_Gc%a_eta(a),mats_Gc%J_zeta_a(:,a))
                  x_ht_ht=Calc_HT_HT(nq,X_fc_fc,mats_Gc%a_eta(a),mats_Gc%b_eta(b),mats_Gc%a_zeta_p_b(a,b),mats_Gc%j_zeta_a(:,a),mats_Gc%j_zeta_p_b(:,b))
                  X_Gc(:,a,b,i)=X_Gc(:,a,b,i)+x_fc_ht*c_arr(2)+x_ht_fc*c_arr(3)+x_ht_ht*c_arr(4)
                  
               end if
               if(ht2)then
                  x_fc_ht2=Calc_FC_HT2(nq,X_fc_fc,tmss%u(a),mats_ap%eta_b2_eta(b),mats_ap%b2_zeta_f(b),mats_ap%j_zeta_p_b2_eta(:,b))
                  x_ht2_fc=Calc_HT2_FC(nq,X_fc_fc,tmss%u(b),mats_ap%eta_a2_eta(a),mats_ap%a2_zeta_f(a),mats_ap%j_zeta_a2_eta(:,a))
                  x_ht_ht2=Calc_HT_HT2(nq,X_fc_fc,mats_ap%a_eta(a),mats_ap%eta_b2_eta(b),mats_ap%b2_zeta_f(b),mats_ap%a_zeta_p_b2_eta(a,b),mats_ap%j_zeta_a(:,a),mats_ap%j_zeta_p_b2_eta(:,b),mats_ap%j_zeta_p_b2_zeta_p_a(:,a,b))
                  x_ht2_ht=Calc_HT2_HT(nq,X_fc_fc,mats_ap%b_eta(b),mats_ap%eta_a2_eta(a),mats_ap%a2_zeta_f(a),mats_ap%b_zeta_p_a2_eta(a,b),mats_ap%j_zeta_p_b(:,b),mats_ap%j_zeta_a2_eta(:,a),mats_ap%j_zeta_a2_zeta_p_b(:,a,b))
                  x_ht2_ht2=Calc_HT2_HT2(nq,X_fc_fc,mats_ap%j_zeta_a2_eta(:,a),mats_ap%b2_zeta_f(b),mats_ap%eta_b2_eta(b),mats_ap%j_zeta_a2_zeta_p_b2_eta(:,a,b),mats_ap%j_zeta_p_b2_eta(:,b),mats_ap%a2_zeta_f(a),mats_ap%eta_a2_eta(a),mats_ap%j_zeta_p_b2_zeta_p_a2_eta(:,a,b),mats_ap%eta_a2_zeta_p_b2_eta(a,b),mats_ap%zeta_p_a2_zeta_p_b2_f(a,b))
                  X_ap(:,a,b,i)=X_ap(:,a,b,i)+x_fc_ht2*c_arr(5)+x_ht2_fc*c_arr(6)+x_ht_ht2*c_arr(7)+x_ht2_ht*c_arr(8)+x_ht2_ht2*c_arr(9)
                  
                  x_fc_ht2=Calc_FC_HT2(nq,X_fc_fc,tmss%u(a),mats_G%eta_b2_eta(b),mats_G%b2_zeta_f(b),mats_G%j_zeta_p_b2_eta(:,b))
                  x_ht2_fc=Calc_HT2_FC(nq,X_fc_fc,tmss%m(b),mats_G%eta_a2_eta(a),mats_G%a2_zeta_f(a),mats_G%j_zeta_a2_eta(:,a))
                  x_ht_ht2=Calc_HT_HT2(nq,X_fc_fc,mats_G%a_eta(a),mats_G%eta_b2_eta(b),mats_G%b2_zeta_f(b),mats_G%a_zeta_p_b2_eta(a,b),mats_G%j_zeta_a(:,a),mats_G%j_zeta_p_b2_eta(:,b),mats_G%j_zeta_p_b2_zeta_p_a(:,a,b))
                  x_ht2_ht=Calc_HT2_HT(nq,X_fc_fc,mats_G%b_eta(b),mats_G%eta_a2_eta(a),mats_G%a2_zeta_f(a),mats_G%b_zeta_p_a2_eta(a,b),mats_G%j_zeta_p_b(:,b),mats_G%j_zeta_a2_eta(:,a),mats_G%j_zeta_a2_zeta_p_b(:,a,b))
                  x_ht2_ht2=Calc_HT2_HT2(nq,X_fc_fc,mats_G%j_zeta_a2_eta(:,a),mats_G%b2_zeta_f(b),mats_G%eta_b2_eta(b),mats_G%j_zeta_a2_zeta_p_b2_eta(:,a,b),mats_G%j_zeta_p_b2_eta(:,b),mats_G%a2_zeta_f(a),mats_G%eta_a2_eta(a),mats_G%j_zeta_p_b2_zeta_p_a2_eta(:,a,b),mats_G%eta_a2_zeta_p_b2_eta(a,b),mats_G%zeta_p_a2_zeta_p_b2_f(a,b))
                  X_G(:,a,b,i)=X_G(:,a,b,i)+x_fc_ht2*c_arr(5)+x_ht2_fc*c_arr(6)+x_ht_ht2*c_arr(7)+x_ht2_ht*c_arr(8)+x_ht2_ht2*c_arr(9)
                  
                  x_fc_ht2=Calc_FC_HT2(nq,X_fc_fc,-tmss%m(a),mats_Gc%eta_b2_eta(b),mats_Gc%b2_zeta_f(b),mats_Gc%j_zeta_p_b2_eta(:,b))
                  x_ht2_fc=Calc_HT2_FC(nq,X_fc_fc,tmss%u(b),mats_Gc%eta_a2_eta(a),mats_Gc%a2_zeta_f(a),mats_Gc%j_zeta_a2_eta(:,a))
                  x_ht_ht2=Calc_HT_HT2(nq,X_fc_fc,mats_Gc%a_eta(a),mats_Gc%eta_b2_eta(b),mats_Gc%b2_zeta_f(b),mats_Gc%a_zeta_p_b2_eta(a,b),mats_Gc%j_zeta_a(:,a),mats_Gc%j_zeta_p_b2_eta(:,b),mats_Gc%j_zeta_p_b2_zeta_p_a(:,a,b))
                  x_ht2_ht=Calc_HT2_HT(nq,X_fc_fc,mats_Gc%b_eta(b),mats_Gc%eta_a2_eta(a),mats_Gc%a2_zeta_f(a),mats_Gc%b_zeta_p_a2_eta(a,b),mats_Gc%j_zeta_p_b(:,b),mats_Gc%j_zeta_a2_eta(:,a),mats_Gc%j_zeta_a2_zeta_p_b(:,a,b))
                  x_ht2_ht2=Calc_HT2_HT2(nq,X_fc_fc,mats_Gc%j_zeta_a2_eta(:,a),mats_Gc%b2_zeta_f(b),mats_Gc%eta_b2_eta(b),mats_Gc%j_zeta_a2_zeta_p_b2_eta(:,a,b),mats_Gc%j_zeta_p_b2_eta(:,b),mats_Gc%a2_zeta_f(a),mats_Gc%eta_a2_eta(a),mats_Gc%j_zeta_p_b2_zeta_p_a2_eta(:,a,b),mats_Gc%eta_a2_zeta_p_b2_eta(a,b),mats_Gc%zeta_p_a2_zeta_p_b2_f(a,b))
                  X_Gc(:,a,b,i)=X_Gc(:,a,b,i)+x_fc_ht2*c_arr(5)+x_ht2_fc*c_arr(6)+x_ht_ht2*c_arr(7)+x_ht2_ht*c_arr(8)+x_ht2_ht2*c_arr(9)
               end if
            end do
         end do
         if(i_thr==1)then
            if(i==NINT(0.1*chunk))then
               write(output_unit,'(A4)',advance='no')'10%'
               flush(output_unit)
            elseif(i==NINT(0.2*chunk))then
               write(output_unit,'(A4)',advance='no')'20%'
               flush(output_unit)
            elseif(i==NINT(0.3*chunk))then
               write(output_unit,'(A4)',advance='no')'30%'
               flush(output_unit)
            elseif(i==NINT(0.4*chunk))then
               write(output_unit,'(A4)',advance='no')'40%'
               flush(output_unit)
            elseif(i==NINT(0.5*chunk))then
               write(output_unit,'(A4)',advance='no')'50%'
               flush(output_unit)
            elseif(i==NINT(0.6*chunk))then
               write(output_unit,'(A4)',advance='no')'60%'
               flush(output_unit)
            elseif(i==NINT(0.7*chunk))then
               write(output_unit,'(A4)',advance='no')'70%'
               flush(output_unit)
            elseif(i==NINT(0.8*chunk))then
               write(output_unit,'(A4)',advance='no')'80%'
               flush(output_unit)
            elseif(i==NINT(0.9*chunk))then
               write(output_unit,'(A4)',advance='no')'90%'
               flush(output_unit)
            end if
         end if
         do b = 1,3
            do a = 1,3
               X_ap(:,a,b,i)=X_ap(:,a,b,i)*x0*C_G
               X_G(:,a,b,i)=X_G(:,a,b,i)*x0*C_G*iu
               X_Gc(:,a,b,i)=X_Gc(:,a,b,i)*x0*C_G*iu
               do c = 1,3
                  X_A(:,c,a,b,i)=X_A(:,c,a,b,i)*x0*C_g
                  X_Ac(:,c,a,b,i)=X_Ac(:,c,a,b,i)*x0*C_g
               end do
            end do
         end do
      end do
      !$OMP END DO
      deallocate(eta,zeta,zeta_p,Fi,Bmat,aa,bb,abia)
      deallocate(x_fc_fc,jt_eta)
      if(ht)then
         deallocate(x_fc_ht,x_ht_fc,x_ht_ht)
         deallocate(j_zeta,j_zeta_p)
         ! deallocate(jt_zeta_p_dip,dip_zeta_p_dip,jt_zeta_dip)
         ! deallocate(jt_zeta_p_mag,dip_zeta_p_mag,mag_zeta_p_dip,jt_zeta_mag)
         ! deallocate(jt_zeta_p_quad,dip_zeta_p_quad,quad_zeta_p_dip,jt_zeta_quad)
      end if
      if(ht2)then
         deallocate(x_fc_ht2,x_ht2_fc,x_ht_ht2,x_ht2_ht,x_ht2_ht2)
      end if
      call deAllocateMats_T1T1(nq,mats_ap,ht,ht2)
      call deAllocateMats_T1T1(nq,mats_G,ht,ht2)
      call deAllocateMats_T1T1(nq,mats_Gc,ht,ht2)
      call deAllocateMats_T1T2(nq,mats_A,ht,ht2)
      call deAllocateMats_T2T1(nq,mats_Ac,ht,ht2)
      
      !$OMP END PARALLEL
      write(output_unit,'(A5)')' 100%'
      flush(output_unit)
      deallocate(c_g)
   end subroutine Make_Corrf_SplitPropagator
   
   pure function Calc_FC_HT(nq,X_fc_fc,A_a,B_eta_b,J_zeta_p_B_b)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      
      double precision,intent(in) :: A_a
      double complex,intent(in) :: B_eta_b,J_zeta_p_B_b(nq)
      
      double complex :: res(nq)
      res=A_a*(X_fc_fc*B_eta_b+J_zeta_p_B_b)
   end function Calc_FC_HT
   
   pure function Calc_HT_FC(nq,X_fc_fc,B_b,A_eta_a,J_zeta_A_a)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      
      double precision,intent(in) :: B_b
      double complex,intent(in) :: A_eta_a,J_zeta_A_a(nq)
      
      double complex :: res(nq)
      res=B_b*(X_fc_fc*A_eta_a+J_zeta_A_a)
   end function Calc_HT_FC
   
   pure function Calc_HT_HT(nq,X_fc_fc,A_eta_a,B_eta_b,A_a_zeta_p_B_b,j_zeta_A_a,j_zeta_p_B_b)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      
      double complex,intent(in) :: A_eta_a,B_eta_b,A_a_zeta_p_B_b,j_zeta_A_a(nq),j_zeta_p_B_b(nq)
      
      double complex :: res(nq)
      res=(A_eta_a*B_eta_b+A_a_zeta_p_B_b)*X_fc_fc+j_zeta_p_B_b*A_eta_a+j_zeta_A_a*B_eta_b
   end function Calc_HT_HT
   
   pure function Calc_FC_HT2(nq,X_fc_fc,A_a,eta_b2_eta,b2_zeta_f,j_zeta_p_b2_eta)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      double precision,intent(in) :: A_a
      double complex,intent(in) :: eta_b2_eta,b2_zeta_f,j_zeta_p_b2_eta(nq)
      
      double complex :: res(nq)
      res=A_a*(X_fc_fc*(eta_b2_eta+b2_zeta_f)+2*j_zeta_p_b2_eta)
   end function Calc_FC_HT2
   
   pure function Calc_HT2_FC(nq,X_fc_fc,B_b,eta_a2_eta,a2_zeta_f,j_zeta_a2_eta)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      double precision,intent(in) :: B_b
      double complex,intent(in) :: eta_a2_eta,a2_zeta_f,j_zeta_a2_eta(nq)
      
      double complex :: res(nq)
      res=B_b*(X_fc_fc*(eta_a2_eta+a2_zeta_f)+2*j_zeta_a2_eta)
   end function Calc_HT2_FC
   
   pure function Calc_HT_HT2(nq,X_fc_fc,eta_A,eta_b2_eta,b2_zeta_f,a_zeta_p_b2_eta,j_zeta_a,j_zeta_p_b2_eta,j_zeta_p_b2_zeta_p_a)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      
      double complex,intent(in) :: eta_A,eta_b2_eta,b2_zeta_f,a_zeta_p_b2_eta,j_zeta_a(nq),j_zeta_p_b2_eta(nq),j_zeta_p_b2_zeta_p_a(nq)
      
      double complex :: res(nq)
      res=X_fc_fc*(eta_a*(eta_b2_eta+b2_zeta_f)+2*a_zeta_p_b2_eta)+j_zeta_a*(eta_b2_eta+b2_zeta_f)+2*eta_a*j_zeta_p_b2_eta+2*J_zeta_p_b2_zeta_p_a
   end function Calc_HT_HT2
   
   pure function Calc_HT2_HT(nq,X_fc_fc,eta_b,eta_a2_eta,a2_zeta_f,b_zeta_p_a2_eta,j_zeta_p_b,j_zeta_a2_eta,j_zeta_a2_zeta_p_b)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      
      double complex,intent(in) :: eta_b,eta_a2_eta,a2_zeta_f,b_zeta_p_a2_eta,j_zeta_p_b(nq),j_zeta_a2_eta(nq),j_zeta_a2_zeta_p_b(nq)
      
      double complex :: res(nq)
      res=X_fc_fc*(eta_b*(eta_a2_eta+a2_zeta_f)+2*b_zeta_p_a2_eta)+j_zeta_p_b*(eta_a2_eta+a2_zeta_f)+2*eta_b*j_zeta_a2_eta+2*j_zeta_a2_zeta_p_b
   end function Calc_HT2_HT
   
   pure function Calc_HT2_HT2(nq,X_fc_fc,j_zeta_a2_eta,b2_zeta_f,eta_b2_eta,j_zeta_a2_zeta_p_b2_eta,j_zeta_p_b2_eta,a2_zeta_f,eta_a2_eta,j_zeta_p_b2_zeta_p_a2_eta,eta_a2_zeta_p_b2_eta,zeta_p_a2_zeta_p_b2_f)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: X_fc_fc(nq)
      
      double complex,intent(in) :: j_zeta_a2_eta(nq),b2_zeta_f,eta_b2_eta,j_zeta_a2_zeta_p_b2_eta(nq),j_zeta_p_b2_eta(nq),a2_zeta_f,eta_a2_eta,j_zeta_p_b2_zeta_p_a2_eta(nq),eta_a2_zeta_p_b2_eta,zeta_p_a2_zeta_p_b2_f
      
      double complex :: res(nq)
      res=2*j_zeta_a2_eta*(b2_zeta_f+eta_b2_eta)+4*j_zeta_a2_zeta_p_b2_eta+2*j_zeta_p_b2_eta*(a2_zeta_f+eta_a2_eta)+4*j_zeta_p_b2_zeta_p_a2_eta &
       +X_fc_fc*((eta_a2_eta+a2_zeta_f)*(eta_b2_eta+b2_zeta_f)+4*eta_a2_zeta_p_b2_eta+2*zeta_p_a2_zeta_p_b2_f)
   end function Calc_HT2_HT2
   
   pure subroutine AllocateMats_T1T1(nq,mats,ht,ht2)
      integer,intent(in) :: nq
      type(TD_SplitP_Mats_T1T1),intent(inout) :: mats
      logical,intent(in) :: ht,ht2
      
      if(ht)allocate(mats%j_zeta_p_b(nq,3),mats%j_zeta_a(nq,3))
      if(.not.ht2)return
      allocate(mats%j_zeta_p_b2_eta(nq,3),mats%j_zeta_a2_eta(nq,3))
      allocate(mats%j_zeta_p_b2_zeta_p_a(nq,3,3),mats%j_zeta_a2_zeta_p_b(nq,3,3))
      allocate(mats%j_zeta_a2_zeta_p_b2_eta(nq,3,3),mats%j_zeta_p_b2_zeta_p_a2_eta(nq,3,3))
   end subroutine AllocateMats_T1T1
   
   pure subroutine DeAllocateMats_T1T1(nq,mats,ht,ht2)
      integer,intent(in) :: nq
      type(TD_SplitP_Mats_T1T1),intent(inout) :: mats
      logical,intent(in) :: ht,ht2
      
      if(ht)then
         deallocate(mats%j_zeta_p_b)
         deallocate(mats%j_zeta_a)
      end if
      if(.not.ht2)return
      deallocate(mats%j_zeta_p_b2_eta,mats%j_zeta_a2_eta)
      deallocate(mats%j_zeta_p_b2_zeta_p_a,mats%j_zeta_a2_zeta_p_b)
      deallocate(mats%j_zeta_a2_zeta_p_b2_eta,mats%j_zeta_p_b2_zeta_p_a2_eta)
   end subroutine DeAllocateMats_T1T1
   
   pure subroutine AllocateMats_T1T2(nq,mats,ht,ht2)
      integer,intent(in) :: nq
      type(TD_SplitP_Mats_T1T2),intent(inout) :: mats
      logical,intent(in) :: ht,ht2
      
      if(ht)allocate(mats%j_zeta_p_b(nq,3,3),mats%j_zeta_a(nq,3))
      if(.not.ht2)return
      allocate(mats%j_zeta_p_b2_eta(nq,3,3),mats%j_zeta_a2_eta(nq,3))
      allocate(mats%j_zeta_p_b2_zeta_p_a(nq,3,3,3),mats%j_zeta_a2_zeta_p_b(nq,3,3,3))
      allocate(mats%j_zeta_a2_zeta_p_b2_eta(nq,3,3,3),mats%j_zeta_p_b2_zeta_p_a2_eta(nq,3,3,3))
   end subroutine AllocateMats_T1T2
   
   pure subroutine DeAllocateMats_T1T2(nq,mats,ht,ht2)
      integer,intent(in) :: nq
      type(TD_SplitP_Mats_T1T2),intent(inout) :: mats
      logical,intent(in) :: ht,ht2
      
      if(ht)deallocate(mats%j_zeta_p_b,mats%j_zeta_a)
      if(.not.ht2)return
      deallocate(mats%j_zeta_p_b2_eta,mats%j_zeta_a2_eta)
      deallocate(mats%j_zeta_p_b2_zeta_p_a,mats%j_zeta_a2_zeta_p_b)
      deallocate(mats%j_zeta_a2_zeta_p_b2_eta,mats%j_zeta_p_b2_zeta_p_a2_eta)
   end subroutine DeAllocateMats_T1T2
   
   pure subroutine AllocateMats_T2T1(nq,mats,ht,ht2)
      integer,intent(in) :: nq
      type(TD_SplitP_Mats_T2T1),intent(inout) :: mats
      logical,intent(in) :: ht,ht2
      
      if(ht)allocate(mats%j_zeta_p_b(nq,3),mats%j_zeta_a(nq,3,3))
      if(.not.ht2)return
      allocate(mats%j_zeta_p_b2_eta(nq,3),mats%j_zeta_a2_eta(nq,3,3))
      allocate(mats%j_zeta_p_b2_zeta_p_a(nq,3,3,3),mats%j_zeta_a2_zeta_p_b(nq,3,3,3))
      allocate(mats%j_zeta_a2_zeta_p_b2_eta(nq,3,3,3),mats%j_zeta_p_b2_zeta_p_a2_eta(nq,3,3,3))
   end subroutine AllocateMats_T2T1
   
   pure subroutine DeAllocateMats_T2T1(nq,mats,ht,ht2)
      integer,intent(in) :: nq
      type(TD_SplitP_Mats_T2T1),intent(inout) :: mats
      logical,intent(in) :: ht,ht2
      
      if(ht)deallocate(mats%j_zeta_p_b,mats%j_zeta_a)
      if(.not.ht2)return
      deallocate(mats%j_zeta_p_b2_eta,mats%j_zeta_a2_eta)
      deallocate(mats%j_zeta_p_b2_zeta_p_a,mats%j_zeta_a2_zeta_p_b)
      deallocate(mats%j_zeta_a2_zeta_p_b2_eta,mats%j_zeta_p_b2_zeta_p_a2_eta)
   end subroutine DeAllocateMats_T2T1
   
   !Workhorse
   subroutine MakeMatrices(ht,ht2,nq,J,eta,zeta,zeta_p,tmss, &
      m_ap,m_G,m_Gc,m_A,m_Ac,j_zeta,j_zeta_p)
      integer,intent(in) :: nq
      type(TMS),intent(in) :: tmss
      logical,intent(in) :: ht2,ht
      double complex,intent(inout) :: j_zeta(nq,nq),j_zeta_p(nq,nq) !passed in work arrays
      type(TD_SplitP_Mats_T1T1),intent(inout) :: m_ap,m_G,m_Gc
      type(TD_SplitP_Mats_T1T2),intent(inout) :: m_A
      type(TD_SplitP_Mats_T2T1),intent(inout) :: m_Ac
      double precision,intent(in) :: J(nq,nq)
      double complex,intent(in) :: eta(nq),zeta(nq,nq),zeta_p(nq,nq)
      double complex :: dip_eta(3),mag_eta(3),quad_eta(3,3)
      double complex :: j_zeta_dip(nq,3),j_zeta_p_dip(nq,3)
      double complex :: j_zeta_mag(nq,3),j_zeta_p_mag(nq,3)
      double complex :: j_zeta_quad(nq,3,3),j_zeta_p_quad(nq,3,3)
      double complex :: zeta_p_dip(nq,3),zeta_p_mag(nq,3),zeta_p_quad(nq,3,3)
      
      double complex :: tmp(nq,nq)
      !HT2 contribution
      double complex :: d2_eta(nq,3),m2_eta(nq,3),q2_eta(nq,3,3)
      double complex :: eta_d2(nq,3),eta_m2(nq,3),eta_q2(nq,3,3)
      double complex :: eta_d2_eta(3),eta_m2_eta(3),eta_q2_eta(3,3)
      double complex :: d2_zeta_f(3),m2_zeta_f(3),q2_zeta_f(3,3)
      double complex :: dip_zeta_p(nq,3),mag_zeta_p(nq,3),quad_zeta_p(nq,3,3)
      double complex :: j_zeta_d2(nq,nq,3),j_zeta_p_d2(nq,nq,3),j_zeta_m2(nq,nq,3),j_zeta_p_m2(nq,nq,3),j_zeta_q2(nq,nq,3,3),j_zeta_p_q2(nq,nq,3,3)
      double complex :: zeta_p_d2(nq,nq,3),zeta_p_m2(nq,nq,3),zeta_p_q2(nq,nq,3,3)
      double complex :: j_zeta_d2_eta(nq,3),j_zeta_p_d2_eta(nq,3),j_zeta_m2_eta(nq,3),j_zeta_p_m2_eta(nq,3),j_zeta_q2_eta(nq,3,3),j_zeta_p_q2_eta(nq,3,3)
      
      
      integer a,b,c
      
      
      if(.not.(ht.or.ht2))return
      J_zeta=matmul(J,zeta)
      j_zeta_p=matmul(J,zeta_p)
      
      do a=1,3 !dipole and magnet non-cross-terms
         dip_eta(a)=zd_dot(nq,eta,tmss%du(:,a))
         mag_eta(a)=zd_dot(nq,eta,tmss%dm(:,a))
         call zd_mv_mult_sym_2times2(nq,j_zeta,j_zeta_p,tmss%du(:,a),tmss%dm(:,a),j_zeta_dip(:,a),j_zeta_p_dip(:,a),j_zeta_mag(:,a),j_zeta_p_mag(:,a))
         zeta_p_dip(:,a)=matmul(zeta_p,tmss%du(:,a))
         zeta_p_mag(:,a)=matmul(zeta_p,tmss%dm(:,a))
         
         if(ht2)then
            d2_eta(:,a)=matmul(tmss%du2(:,:,a),eta)
            eta_d2(:,a)=matmul(eta,tmss%du2(:,:,a))
            eta_d2_eta(a)=zz_dot(nq,eta,d2_eta(:,a))
            m2_eta(:,a)=matmul(tmss%dm2(:,:,a),eta)
            eta_m2(:,a)=matmul(eta,tmss%dm2(:,:,a))
            eta_m2_eta(a)=zz_dot(nq,eta,m2_eta(:,a))
            
            call dz_frob_sym(nq,tmss%du2(:,:,a),zeta,d2_zeta_f(a))
            call dz_frob_sym(nq,tmss%dm2(:,:,a),zeta,m2_zeta_f(a))
            
            dip_zeta_p(:,a)=matmul(tmss%du(:,a),zeta_p)
            mag_zeta_p(:,a)=matmul(tmss%dm(:,a),zeta_p)
            
            j_zeta_d2_eta(:,a)=matmul(J_zeta,d2_eta(:,a))
            j_zeta_p_d2_eta(:,a)=matmul(J_zeta_p,d2_eta(:,a))
            j_zeta_m2_eta(:,a)=matmul(J_zeta,m2_eta(:,a))
            j_zeta_p_m2_eta(:,a)=matmul(J_zeta_p,m2_eta(:,a))
            
            !TODO, check whether my own matrix multiplication would be faster here
            !I have a feeling that matmul implicitly casts double precision arrays to double complex
            j_zeta_d2(:,:,a)=matmul(J_zeta,tmss%du2(:,:,a))
            zeta_p_d2(:,:,a)=matmul(zeta_p,tmss%du2(:,:,a))
            j_zeta_p_d2(:,:,a)=matmul(J,zeta_p_d2(:,:,a))
            
            j_zeta_m2(:,:,a)=matmul(J_zeta,tmss%dm2(:,:,a))
            zeta_p_m2(:,:,a)=matmul(zeta_p,tmss%dm2(:,:,a))
            j_zeta_p_m2(:,:,a)=matmul(J,zeta_p_m2(:,:,a))
            
         end if
         
         
         do b=1,3 !quadrupole non-cross-terms
            zeta_p_quad(:,b,a)=matmul(zeta_p,tmss%dq(:,b,a))
            if(ht2)then
               q2_eta(:,b,a)=matmul(tmss%dq2(:,:,b,a),eta)
               eta_q2(:,b,a)=matmul(eta,tmss%dq2(:,:,b,a))
               eta_q2_eta(b,a)=zz_dot(nq,eta,q2_eta(:,b,a))
               
               call dz_frob_sym(nq,tmss%dq2(:,:,b,a),zeta,q2_zeta_f(b,a))
               
               quad_zeta_p(:,b,a)=matmul(tmss%dq(:,b,a),zeta_p)
               
               j_zeta_q2_eta(:,b,a)=matmul(J_zeta,q2_eta(:,b,a))
               j_zeta_p_q2_eta(:,b,a)=matmul(J_zeta_p,q2_eta(:,b,a))
               
               j_zeta_q2(:,:,b,a)=matmul(J_zeta,tmss%dq2(:,:,b,a))
               zeta_p_q2(:,:,b,a)=matmul(zeta_p,tmss%dq2(:,:,b,a))
               j_zeta_p_q2(:,:,b,a)=matmul(J,zeta_p_q2(:,:,b,a))
            end if
            quad_eta(b,a)=zd_dot(nq,eta,tmss%dq(:,b,a))
            j_zeta_p_quad(:,b,a)=matmul(j_zeta_P,tmss%dq(:,b,a))
            j_zeta_quad(:,b,a)=matmul(j_zeta_P,tmss%dq(:,b,a))
         end do
      end do
      
      !Terms that are not cross-terms
      m_ap%a_eta=dip_eta
      m_ap%b_eta=dip_eta
      m_ap%j_zeta_p_b=j_zeta_p_dip
      m_ap%j_zeta_a=j_zeta_dip
      
      m_G%a_eta=dip_eta
      m_G%b_eta=mag_eta
      m_G%j_zeta_p_b=j_zeta_p_mag
      m_G%j_zeta_a=j_zeta_dip
      
      m_Gc%b_eta=dip_eta
      m_Gc%a_eta=-mag_eta
      m_Gc%j_zeta_p_b=j_zeta_p_dip
      m_Gc%j_zeta_a=-j_zeta_mag
      
      m_A%a_eta=dip_eta
      m_A%b_eta=quad_eta
      m_A%j_zeta_a=j_zeta_dip
      m_A%j_zeta_p_b=j_zeta_p_quad
      
      m_Ac%a_eta=quad_eta
      m_Ac%b_eta=dip_eta
      m_Ac%j_zeta_a=j_zeta_quad
      m_Ac%j_zeta_p_b=j_zeta_p_dip
      
      if(ht2)then
         m_ap%eta_a2_eta=eta_d2_eta
         m_ap%eta_b2_eta=eta_d2_eta
         m_ap%a2_zeta_f=d2_zeta_f
         m_ap%b2_zeta_f=d2_zeta_f
         m_ap%j_zeta_a2_eta=j_zeta_d2_eta
         m_ap%j_zeta_p_b2_eta=j_zeta_p_d2_eta
         
         m_G%eta_a2_eta=eta_d2_eta
         m_G%eta_b2_eta=eta_m2_eta
         m_G%a2_zeta_f=d2_zeta_f
         m_G%b2_zeta_f=m2_zeta_f
         m_G%j_zeta_a2_eta=j_zeta_d2_eta
         m_G%j_zeta_p_b2_eta=j_zeta_p_m2_eta
         
         m_Gc%eta_a2_eta=-eta_m2_eta
         m_Gc%eta_b2_eta=eta_d2_eta
         m_Gc%a2_zeta_f=-m2_zeta_f
         m_Gc%b2_zeta_f=d2_zeta_f
         m_Gc%j_zeta_a2_eta=-j_zeta_m2_eta
         m_GC%j_zeta_p_b2_eta=j_zeta_p_d2_eta
         
         m_A%eta_a2_eta=eta_d2_eta
         m_a%eta_b2_eta=eta_q2_eta
         m_A%a2_zeta_f=d2_zeta_f
         m_A%b2_zeta_f=q2_zeta_f
         m_A%j_zeta_a2_eta=j_zeta_d2_eta
         m_A%j_zeta_p_b2_eta=j_zeta_p_q2_eta
         
         m_Ac%eta_a2_eta=eta_q2_eta
         m_Ac%eta_b2_eta=eta_d2_eta
         m_Ac%a2_zeta_f=q2_zeta_f
         m_Ac%b2_zeta_f=d2_zeta_f
         m_Ac%j_zeta_a2_eta=j_zeta_q2_eta
         m_Ac%j_zeta_p_b2_eta=j_zeta_p_d2_eta
      end if
      
      do a = 1,3
         do b = 1,3
            m_ap%a_zeta_p_b(a,b)=dz_dot(nq,tmss%du(:,a),zeta_p_dip(:,b))
            m_G%a_zeta_p_b(a,b)=dz_dot(nq,tmss%du(:,a),zeta_p_mag(:,b))
            m_Gc%a_zeta_p_b(a,b)=-dz_dot(nq,tmss%dm(:,a),zeta_p_dip(:,b))
            
            if(ht2)then
               m_ap%a_zeta_p_b2_eta(a,b)=zz_dot(nq,dip_zeta_p(:,a),d2_eta(:,b))
               m_ap%b_zeta_p_a2_eta(a,b)=zz_dot(nq,dip_zeta_p(:,b),d2_eta(:,a))
               m_ap%j_zeta_p_b2_zeta_p_a(:,a,b)=matmul(j_zeta_p_d2(:,:,b),zeta_p_dip(:,a))
               m_ap%j_zeta_a2_zeta_p_b(:,a,b)=matmul(j_zeta_d2(:,:,a),zeta_p_dip(:,b))
               m_ap%j_zeta_a2_zeta_p_b2_eta(:,a,b)=matmul(j_zeta_d2(:,:,a),matmul(zeta_p,d2_eta(:,b)))
               m_ap%j_zeta_p_b2_zeta_p_a2_eta(:,a,b)=matmul(j_zeta_p_d2(:,:,b),matmul(zeta_p,d2_eta(:,a)))
               m_ap%eta_a2_zeta_p_b2_eta(a,b)=zz_dot(nq,eta_d2(:,a),matmul(zeta_p,d2_eta(:,b)))
               tmp=matmul(tmss%du2(:,:,a),zeta_p_d2(:,:,b))
               call zz_frob_sym(nq,zeta_p,tmp,m_ap%zeta_p_a2_zeta_p_b2_f(a,b))
               
               m_G%a_zeta_p_b2_eta(a,b)=zz_dot(nq,dip_zeta_p(:,a),m2_eta(:,b))
               m_G%b_zeta_p_a2_eta(a,b)=zz_dot(nq,mag_zeta_p(:,b),d2_eta(:,a))
               m_G%j_zeta_p_b2_zeta_p_a(:,a,b)=matmul(j_zeta_p_m2(:,:,b),zeta_p_dip(:,a))
               m_G%j_zeta_a2_zeta_p_b(:,a,b)=matmul(j_zeta_d2(:,:,a),zeta_p_mag(:,b))
               m_G%j_zeta_a2_zeta_p_b2_eta(:,a,b)=matmul(j_zeta_d2(:,:,a),matmul(zeta_p,m2_eta(:,b)))
               m_G%j_zeta_p_b2_zeta_p_a2_eta(:,a,b)=matmul(j_zeta_p_m2(:,:,b),matmul(zeta_p,d2_eta(:,a)))
               m_G%eta_a2_zeta_p_b2_eta(a,b)=zz_dot(nq,eta_d2(:,a),matmul(zeta_p,m2_eta(:,b)))
               tmp=matmul(tmss%du2(:,:,a),zeta_p_m2(:,:,b))
               call zz_frob_sym(nq,zeta_p,tmp,m_G%zeta_p_a2_zeta_p_b2_f(a,b))
               
               m_Gc%a_zeta_p_b2_eta(a,b)=-zz_dot(nq, mag_zeta_p(:,a), d2_eta(:,b))
               m_Gc%b_zeta_p_a2_eta(a,b)=-zz_dot(nq, dip_zeta_p(:,b), m2_eta(:,a))
               m_Gc%j_zeta_p_b2_zeta_p_a(:,a,b)=-matmul(j_zeta_p_d2(:,:,b),zeta_p_mag(:,a))
               m_Gc%j_zeta_a2_zeta_p_b(:,a,b)=-matmul(j_zeta_m2(:,:,a),zeta_p_dip(:,b))
               m_Gc%j_zeta_a2_zeta_p_b2_eta(:,a,b)=-matmul(j_zeta_m2(:,:,a),matmul(zeta_p,d2_eta(:,b)))
               m_Gc%j_zeta_p_b2_zeta_p_a2_eta(:,a,b)=-matmul(j_zeta_p_d2(:,:,b),matmul(zeta_p,m2_eta(:,a)))
               m_Gc%eta_a2_zeta_p_b2_eta(a,b)=-zz_dot(nq,eta_m2(:,a),matmul(zeta_p,d2_eta(:,b)))
               tmp=-matmul(tmss%dm2(:,:,a),zeta_p_d2(:,:,b))
               call zz_frob_sym(nq,zeta_p,tmp,m_Gc%zeta_p_a2_zeta_p_b2_f(a,b))
            end if
            
            do c = 1,3
               m_A%a_zeta_p_b(a,b,c)=dz_dot(nq,tmss%du(:,a),zeta_p_quad(:,b,c))
               
               m_Ac%a_zeta_p_b(a,b,c)=dz_dot(nq,tmss%dq(:,a,b),zeta_p_dip(:,c)) !TODO, not certain about the indexing of Ac
               
               if(ht2)then
                  m_A%a_zeta_p_b2_eta(a,b,c)=zz_dot(nq,dip_zeta_p(:,a),q2_eta(:,b,c))
                  m_A%b_zeta_p_a2_eta(a,b,c)=zz_dot(nq,quad_zeta_p(:,b,c),d2_eta(:,a))
                  m_A%j_zeta_p_b2_zeta_p_a(:,a,b,c)=matmul(j_zeta_p_q2(:,:,b,c),zeta_p_dip(:,a))
                  m_A%j_zeta_a2_zeta_p_b(:,a,b,c)=matmul(j_zeta_d2(:,:,a),zeta_p_quad(:,b,c))
                  m_A%j_zeta_a2_zeta_p_b2_eta(:,a,b,c)=matmul(j_zeta_d2(:,:,a),matmul(zeta_p,q2_eta(:,b,c)))
                  m_A%j_zeta_p_b2_zeta_p_a2_eta(:,a,b,c)=matmul(j_zeta_p_q2(:,:,b,c),matmul(zeta_p,d2_eta(:,a)))
                  m_A%eta_a2_zeta_p_b2_eta(a,b,c)=zz_dot(nq,eta_d2(:,a),matmul(zeta_p,q2_eta(:,b,c)))
                  tmp=matmul(tmss%du2(:,:,a),zeta_p_q2(:,:,b,c))
                  call zz_frob_sym(nq,zeta_p,tmp,m_A%zeta_p_a2_zeta_p_b2_f(a,b,c))
                  
                  m_Ac%a_zeta_p_b2_eta(a,b,c)=zz_dot(nq,quad_zeta_p(:,a,b),d2_eta(:,c))
                  m_Ac%b_zeta_p_a2_eta(a,b,c)=zz_dot(nq,dip_zeta_p(:,c),q2_eta(:,a,b))
                  m_Ac%j_zeta_p_b2_zeta_p_a(:,a,b,c)=matmul(j_zeta_p_d2(:,:,c),zeta_p_quad(:,a,b))
                  m_Ac%j_zeta_a2_zeta_p_b(:,a,b,c)=matmul(j_zeta_q2(:,:,a,b),zeta_p_dip(:,c))
                  m_Ac%j_zeta_a2_zeta_p_b2_eta(:,a,b,c)=matmul(j_zeta_q2(:,:,a,b),matmul(zeta_p,d2_eta(:,c)))
                  m_Ac%j_zeta_p_b2_zeta_p_a2_eta(:,a,b,c)=matmul(j_zeta_p_d2(:,:,c),matmul(zeta_p,q2_eta(:,a,b)))
                  m_Ac%eta_a2_zeta_p_b2_eta(a,b,c)=zz_dot(nq,eta_q2(:,a,b),matmul(zeta_p,d2_eta(:,c)))
                  tmp=matmul(tmss%dq2(:,:,a,b),zeta_p_d2(:,:,c))
                  call zz_frob_sym(nq,zeta_p,tmp,m_Ac%zeta_p_a2_zeta_p_b2_f(a,b,c))
               end if
            end do
         end do
      end do
   end subroutine MakeMatrices
   
   subroutine dz_frob_sym_twice(nq,a1,a2,b,res1,res2)
      integer,intent(in) :: nq
      double precision,intent(in) :: a1(nq,nq),a2(nq,nq)
      double complex,intent(in) :: b(nq,nq)
      integer i,j
      double complex,intent(out) :: res1,res2
      
      res1=0d0
      res2=0d0
      do i = 1,nq
         res1=res1+a1(i,i)*b(i,i)
         res2=res2+a2(i,i)*b(i,i)
         !$OMP SIMD REDUCTION(+:res1,res2)
         do j = i+1,nq
            res1=res1+2*a1(j,i)*b(j,i)
            res2=res2+2*a2(j,i)*b(j,i)
         end do
         !$OMP END SIMD
      end do
   end subroutine dz_frob_sym_twice
   
   subroutine dz_frob_sym(nq,a1,b,res1)
      integer,intent(in) :: nq
      double precision,intent(in) :: a1(nq,nq)
      double complex,intent(in) :: b(nq,nq)
      integer i,j
      double complex,intent(out) :: res1
      
      res1=0d0
      do i = 1,nq
         res1=res1+a1(i,i)*b(i,i)
         !$OMP SIMD REDUCTION(+:res1)
         do j = i+1,nq
            res1=res1+2*a1(j,i)*b(j,i)
         end do
         !$OMP END SIMD
      end do
   end subroutine dz_frob_sym
   
   subroutine zz_frob_sym(nq,a1,b,res1)
      integer,intent(in) :: nq
      double complex,intent(in) :: a1(nq,nq)
      double complex,intent(in) :: b(nq,nq)
      integer i,j
      double complex,intent(out) :: res1
      
      res1=0d0
      do i = 1,nq
         res1=res1+a1(i,i)*b(i,i)
         !$OMP SIMD REDUCTION(+:res1)
         do j = i+1,nq
            res1=res1+2*a1(j,i)*b(j,i)
         end do
         !$OMP END SIMD
      end do
   end subroutine zz_frob_sym
   
   function dzd_mmm_mult_sym(nq,a,b,c)result(res)
      integer,intent(in) :: nq
      double precision,intent(in) :: a(nq,nq),c(nq,nq)
      double complex,intent(in) :: b(nq,nq)
      double complex res(nq,nq),M(nq,nq)
      
      integer i,j,k,l
      
      res=0D0
      M=0d0
      do j = 1, nq
        do k = 1, nq
            do i = 1, nq
               M(i,j) = M(i,j) + b(i,k) * c(k,j)
            end do
        end do
      end do
      
      do j = 1, nq
        do k = 1, nq
            do i = 1, nq
                res(i,j) = res(i,j) + a(i,k) * M(k,j)
            end do
        end do
      end do
      
   end function dzd_mmm_mult_sym
   
   function zdz_vmv_mult_sym(nq,u,a,v)result(res)
      integer,intent(in) :: nq
      double precision,intent(in) :: a(nq,nq)
      double complex,intent(in) :: v(nq),u(nq)
      double complex res,summ
      
      integer i,j
      
      res=(0d0,0d0)
      do i = 1,nq
         summ=(0d0,0d0)
         do j = 1,nq
            summ=summ+a(j,i)*v(j)
         end do
         res=res+summ*u(i)
      end do
   end function zdz_vmv_mult_sym
   
   function dzd_vmv_mult_sym(nq,u,a,v)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq,nq)
      double precision,intent(in) :: v(nq),u(nq)
      double complex res,summ
      integer i,j
      
      res=(0d0,0d0)
      do i = 1,nq
         summ=0d0
         !$OMP SIMD &
         !$OMP REDUCTION(+:summ)
         do j = 1,nq
            summ=summ+a(j,i)*v(j) !matrix a is symmetric
         end do
         !$OMP END SIMD
         res=res+summ*u(i)
      end do
   end function dzd_vmv_mult_sym
   
   subroutine dzd_vmv_mult_sym_thrice(nq,a,u,v,res1,res2,res3)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq,nq)
      double precision,intent(in) :: v(nq),u(nq)
      double complex summ1,summ2,summ3,a_ji
      double complex,intent(out) :: res1,res2,res3
      integer i,j
      
      res1=(0d0,0d0)
      res2=(0d0,0d0)
      res3=(0d0,0d0)
      do i = 1,nq
         summ1=0d0
         summ2=0d0
         summ3=0d0
         !$OMP SIMD &
         !$OMP REDUCTION(+:summ1,summ2,summ3) PRIVATE(a_ji)
         do j = 1,nq
            a_ji=a(j,i)
            summ1=summ1+a_ji*u(j)
            summ2=summ2+a_ji*v(j)
            summ3=summ3+a_ji*u(j)
         end do
         !$OMP END SIMD
         res1=res1+summ1*u(i)
         res2=res2+summ2*u(i)
         res3=res3+summ3*v(i)
      end do
   end subroutine dzd_vmv_mult_sym_thrice
   
   function zd_mv_mult_sym(nq,a,v)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq,nq)
      double precision,intent(in) :: v(nq)
      double complex res(nq),summ
      integer i,j
      
      do i = 1,nq
         summ=0d0
         !$OMP SIMD &
         !$OMP REDUCTION(+:SUMM)
         do j = 1,nq
            summ=summ+a(j,i)*v(j)
         end do
         !$OMP END SIMD
         res(i)=summ
      end do
   end function zd_mv_mult_sym
   
   subroutine zd_mv_mult_sym_2times2(nq,a,b,v1,v2,resa1,resa2,resb1,resb2)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq,nq),b(nq,nq)
      double precision,intent(in) :: v1(nq),v2(nq)
      double complex summa1,summa2,summb1,summb2
      double complex,intent(out) :: resa1(nq),resa2(nq),resb1(nq),resb2(nq)
      integer i,j
      
      do i = 1,nq
         summa1=0d0
         summa2=0d0
         summb1=0d0
         summb2=0d0
         !$OMP SIMD &
         !$OMP REDUCTION(+:SUMMa1,summa2,summb1,summb2)
         do j = 1,nq
            summa1=summa1+a(j,i)*v1(j)
            summa2=summa2+a(j,i)*v2(j)
            summb1=summb1+b(j,i)*v1(j)
            summb2=summb2+b(j,i)*v2(j)
         end do
         !$OMP END SIMD
         resa1(i)=summa1
         resa2(i)=summa2
         resb1(i)=summb1
         resb2(i)=summb2
      end do
   end subroutine zd_mv_mult_sym_2times2
   
   subroutine dz_mv_mult_sym_twice(nq,a1,a2,v,res1,res2)
      integer,intent(in) :: nq
      double precision,intent(in) :: a1(nq,nq),a2(nq,nq)
      double complex,intent(in) :: v(nq)
      double complex,intent(out) :: res1(nq),res2(nq)
      double complex summ1,summ2
      integer i,j
      
      do i = 1,nq
         summ1=0d0
         summ2=0d0
         !$OMP SIMD &
         !$OMP REDUCTION(+:SUMM1,SUMM2)
         do j = 1,nq
            summ1=summ1+a1(j,i)*v(j)
            summ2=summ2+a2(j,i)*v(j)
         end do
         !$OMP END SIMD
         res1(i)=summ1
         res2(i)=summ2
      end do
   end subroutine dz_mv_mult_sym_twice
   
   function zd_dot(nq,a,b)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq)
      double precision,intent(in) :: b(nq)
      double complex res
      integer i
      
      res=(0d0,0d0)
      !$OMP SIMD REDUCTION(+:res)
      do i =1,nq
         res=res+a(i)*b(i)
      end do
      !$END OMP SIMD
   end function zd_dot
   
   function dz_dot(nq,a,b)result(res)
      integer,intent(in) :: nq
      double precision,intent(in) :: a(nq)
      double complex,intent(in) :: b(nq)
      double complex res
      integer i
      
      res=(0d0,0d0)
      !$OMP SIMD REDUCTION(+:res)
      do i =1,nq
         res=res+a(i)*b(i)
      end do
      !$END OMP SIMD
   end function dz_dot
   
   function zz_dot(nq,a,b)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq)
      double complex,intent(in) :: b(nq)
      double complex res
      integer i
      
      res=(0d0,0d0)
      !$OMP SIMD REDUCTION(+:res)
      do i =1,nq
         res=res+a(i)*b(i)
      end do
      !$END OMP SIMD
   end function zz_dot
   
   subroutine zd_dot_twice(nq,a,b1,b2,res1,res2)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq)
      double precision,intent(in) :: b1(nq),b2(nq)
      double complex,intent(out) :: res1,res2
      double complex :: summi1,summi2
      integer i,j
      
      res1=(0d0,0d0)
      res2=(0d0,0d0)
      !$OMP SIMD REDUCTION(+:res1,res2)
      do i =1,nq
         res1=res1+a(i)*b1(i)
         res2=res2+a(i)*b2(i)
      end do
      !$OMP END SIMD
   end subroutine zd_dot_twice
   
   pure function matmul_aBc(n,Avecs,Bmat,Cvecs)result(res)
      integer,intent(in) :: n
      double precision,intent(in) :: Avecs(3,n),Cvecs(3,n)
      double complex,intent(in) :: Bmat(n,n)
      double complex :: res(3,3)
      integer a,b,i,j
      
      do b = 1,3
         do a = 1,3
            res(a,b)=0d0
            do j = 1,n
               do i = 1,n
                  res(a,b)=res(a,b)+Avecs(a,i)*Bmat(i,j)*Cvecs(b,j)
               end do
            end do
         end do
      end do
   end function matmul_aBc
   
   pure function matmul_aBc_ten2(n,Avecs,Bmat,Cten2)result(res)
      integer,intent(in) :: n
      double precision,intent(in) :: Avecs(3,n),Cten2(3,3,n)
      double complex,intent(in) :: Bmat(n,n)
      double complex :: res(3,3,3)
      integer a,b,c,i,j
      
      do c = 1,3
         do b = 1,3
            do a = 1,3
               res(a,b,c)=0d0
               do j = 1,n
                  do i = 1,n
                     res(a,b,c)=res(a,b,c)+Avecs(a,i)*Bmat(i,j)*Cten2(b,c,j)
                  end do
               end do
            end do
         end do
      end do
   end function matmul_aBc_ten2
   
   pure function matmul_ten2_aBc(n,Aten2,Bmat,Cvecs)result(res)
      integer,intent(in) :: n
      double precision,intent(in) :: Aten2(3,3,n),Cvecs(3,n)
      double complex,intent(in) :: Bmat(n,n)
      double complex :: res(3,3,3)
      integer a,b,c,i,j
      
      do c = 1,3
         do b = 1,3
            do a = 1,3
               res(a,b,c)=0d0
               do j = 1,n
                  do i = 1,n
                     res(a,b,c)=res(a,b,c)+Aten2(b,c,i)*Bmat(i,j)*Cvecs(a,j)
                  end do
               end do
            end do
         end do
      end do
   end function matmul_ten2_aBc
   
   pure function matmul_J_Zeta_v(n,J,zeta,v)result(res)
      integer,intent(in) :: n
      double precision,intent(in) :: J(n,n),v(3,n)
      double complex,intent(in) :: zeta(n,n)
      double complex :: res(n,3)
      integer a,k,i,jj
      
      do a = 1,3
         do k = 1,n
            res(k,a)=0d0
            do jj = 1,n
               do i = 1,n
                  res(k,a)=res(k,a)+J(k,i)*zeta(i,jj)*v(a,jj)
               end do
            end do
         end do
      end do
   end function matmul_J_zeta_v
   
   !t2 = tensor 2nd order (i.e. quadrupole)
   pure function matmul_J_Zeta_t2(n,J,zeta,ten)result(res)
      integer,intent(in) :: n
      double precision,intent(in) :: J(n,n),ten(3,3,n)
      double complex,intent(in) :: zeta(n,n)
      double complex :: res(3,3,n)
      integer a,b,k,i,jj
      
      do b = 1,3
         do a = 1,3
            do k = 1,n
               res(a,b,k)=0d0
               do jj = 1,n
                  do i = 1,n
                     res(a,b,k)=res(a,b,k)+J(k,i)*zeta(i,jj)*ten(a,b,jj)
                  end do
               end do
            end do
         end do
      end do
   end function matmul_J_zeta_t2
   
   pure subroutine x0_ab_SplitPropagator(nq,J,JGJ,wg,we,t,aa,bb,B,abia,eps)
      integer,intent(in) :: nq
      double precision,intent(in) :: wg(nq),we(nq),J(nq,nq),JGJ(nq,nq)
      double precision,intent(in) :: t,eps
      double complex,intent(out) :: aa(nq),bb(nq),abia(nq),B(nq,nq)
      integer i,jj,k
      
      do i = 1,nq
         aa(i)=iu*we(i)/sin(t/2d0*we(i)+iu*eps)
         bb(i)=iu*we(i)/tan(t/2d0*we(i)+iu*eps)
         abia(i)=aa(i)*aa(i)/bb(i)
         B(i,i)=JGJ(i,i)-bb(i)
         do jj = i+1,nq
            B(i,jj)=JGJ(i,jj)
            B(jj,i)=JGJ(jj,i)
         end do
      end do
   end subroutine x0_ab_SplitPropagator
   
   
   function make_x0_SplitPropagator(nq,detG,v,kgk,aa,bb,B,abia,D,Fi,eta,zeta,zeta_p,phasePrev)result(x0)
      integer :: nq,i,j,detPB,detPD,ierr,Piv(nq)
      double complex :: aa(nq),bb(nq),abia(nq),B(nq,nq),bc_ex
      double complex :: D(nq,nq),Fi(nq,nq),eta(nq),zeta(nq,nq),zeta_p(nq,nq),x0,vtbiv,vtzpv
      double precision :: phasePrev,kgk,v(nq)
      double precision :: lndetA,lndetD,lndetB,lndetG,lndetprod,ln2n
      double precision :: phaseA,phaseD,phaseB,phaseG,phaseprod,iu_phase
      type(big_double) detG
      
      lndetG=(log(detG%num)+detG%p*log(10d0))
      phaseG=0d0
      
      call DiagProd_c_big_vec_ln(aa,nq,lndetA,phaseA)
      ln2n=dble(nq)/2d0*log(2d0)
      iu_phase=-pi*nq !TODO, this line is under scrutiny, I got something else in the equations but Gemini said THIS particular phase is the correct one
      
      Fi=B
      do i = 1,nq
         Fi(i,i)=Fi(i,i)+abia(i)
      end do
      
      call zgetrf(nq,nq,B,nq,Piv,ierr) !LU form
      call DetFromLU_C_big_ln(B,nq,piv,lndetB,phaseB,detPB)
      call INV_C_lapack(B,piv,z_lwork,nq) !inverse
      
      D=0d0
      do i = 1,nq
         D(i,i)=bb(i)
         do j = 1,nq
            D(i,j)=D(i,j)+aa(i)*B(i,j)*aa(j)
            D(i,j)=-2d0*D(i,j)
         end do
      end do
      
      call zgetrf(nq,nq,D,nq,Piv,ierr) !LU form
      call DetFromLU_C_big_ln(D,nq,piv,lndetD,phaseD,detPD)
      call INV_C_lapack(D,piv,z_lwork,nq) !inverse
      
      
      call zgetrf(nq,nq,Fi,nq,Piv,ierr) !LU form
      call INV_C_lapack(Fi,piv,z_lwork,nq) !inverse
      
      
      ! do i = 1,nq
         ! eta(i)=0d0
         ! do j = 1,nq
            ! eta(i)=eta(i)-Fi(i,j)*v(j)
            ! zeta(i,j)=0.5d0*(Fi(i,j)+B(i,j))
            ! zeta_p(i,j)=0.5d0*(Fi(i,j)-B(i,j))
         ! end do
      ! end do
      eta=-matmul(Fi,v)
      zeta=0.5d0*(Fi+B)
      zeta_p=0.5d0*(Fi-B)
      
      lndetProd=lndetG/2d0+lndetA-lndetB-lndetD/2d0+ln2n
      phaseProd=phaseG/2d0+phaseA-phaseB-phaseD/2d0+iu_phase
      phaseProd=PhaseNorm(phaseProd)
      ! if(abs(cos(phaseprod)-cos(phaseprev))>0.1d0 .or. abs(sin(phaseprod)-sin(phaseprev))>0.1d0)then
         ! phaseProd=phaseProd+pi
      ! end if
      phasePrev=phaseProd
      
      vTzPv=0d0
      vtbiv=0d0
      do i = 1,nq
         do j = 1,nq
            vtbiv=vtbiv+v(i)*B(j,i)*v(j)
            vtzpv=vtzpv+v(i)*zeta_p(j,i)*v(j)
         end do
      end do
      
      
      bc_ex=exp(lndetProd)*(cos(phaseProd)+iu*sin(phaseProd))
      x0=bc_ex*exp(-kgk+2*vTzPv+vTBiv)
   end function make_x0_SplitPropagator
   
   
   subroutine Make_corrf_k(trType,nq,k_batch,batch_n,wg,we,tmss,v,gamma_gr,detG,kgk,J,K_dusch,gamma,theta,w_ad,N_points,t_max,sparse,J_nz,nz,ht,ht2,tmexp,fixphase, &
   contr,n_thr,tol,X_ap,X_G,X_Gc,X_A,X_Ac,X_0)
      integer kk,kkk,nq,N_points,n_thr,i_thr,batch_n
      integer jj,mm,ii,ll,nz,pp,nn
      integer a,b,c,d,e
      integer ai,bi,ci,di,ei
      integer k,l,m,n,p,q
      integer d_idx,e_idx
      type(list_arr_int32) :: J_nz(batch_n)
      integer(int64) points_c,i,chunk
      integer k_batch(batch_n),kk_idx,trType
      double precision J(nq,nq),K_dusch(nq),t_max,wg(nq),we(nq),v(nq),gamma_gr(nq,nq),damp,damp2,w_ad,phase
      double precision,allocatable :: X_holdr(:),X_holdi(:),JT(:,:),phases(:)
      double precision kgk,wg_k,dt,t,gamma,theta,sqrt_wg_k,sqrt_wg_k_arr(batch_n),tol,cur_sizee,cur_phase
      type(big_double) detG
      type(tms) tmss
      double complex, allocatable :: eta(:)
      double complex :: dx0_dv_,d2x0_dvdv_,d3x0_dvdvdv_,d4x0_dvdvdvdv_,d5x0_dvdvdvdvdv_
      double complex :: dx0_dc_,d2x0_dcdc_
      double complex :: d2x0_dvdc_1,d2x0_dvdc_2,d2x0_dvdc_3
      double complex :: d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6
      double complex :: d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5,d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10
      double complex :: d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5
      
      
      double complex tau,shift,shift2,kgkp,vdv,test,test2,test3,test4
      double complex :: x0,x0_st,x_fcht1,x_fcht2,x_ht
      
      double complex :: qx_l,xq_l,qxq_lm,qqx_lm,xqq_lm,qqxq_mln,qxqq_mln,qqqx_mln,qqqxq_lmnp,qqxqq_lmnp,qqqxqq_lmnpq
      double complex,allocatable :: x_fc(:),x_fc_ht(:),x_ht_fc(:),x_ht_ht(:),x_fc_ht2(:),x_ht2_fc(:),x_ht_ht2(:),x_ht2_ht(:),x_ht2_ht2(:)
      double complex,allocatable :: qqxqq(:,:,:,:),qqqx(:,:,:),qqqxq(:,:,:,:),qqqxqq(:,:,:,:,:),qxqq(:,:,:)
      double complex,allocatable :: qjx(:),xqj(:),qxqj(:,:),qqjx(:,:),qjxqj(:,:),qqjxqj(:,:,:)
      
      double complex,allocatable :: aa_ex(:),C_i(:,:),D_i(:,:),vp(:)
      double complex,allocatable,intent(out) :: X_ap(:,:,:,:),X_G(:,:,:,:),X_Gc(:,:,:,:),X_A(:,:,:,:,:),X_Ac(:,:,:,:,:),X_0(:)
      double complex,allocatable :: w1(:),w2(:,:),w2p(:,:),X_hold(:)
      double complex buf1,buf
      logical ht,ht2,tmexp(2),groundleft,groundright,sparse,gg,ee,st,correctPhaseX,ovcomb
      logical contr(9),fixphase
      
      groundleft=tmexp(1)
      groundright=tmexp(2)
      gg=groundleft.and.groundright
      ee=.not.groundleft.and..not.groundright
      
      points_c=2_8**int(N_points,kind=8)
      chunk=points_c/n_thr
      dt=t_max/dble(points_c-1_8) !Unsure about the last element of x0, similar problem with x of the spectrum, todo check
      allocate(X_ap(3,3,batch_n,points_c),X_G(3,3,batch_n,points_c),X_Gc(3,3,batch_n,points_c))
      allocate(X_A(3,3,3,batch_n,points_c),X_Ac(3,3,3,batch_n,points_c))
      
      X_Ap=0d0
      X_G=0d0
      X_Gc=0d0
      X_A=0d0
      X_Ac=0d0
      
      do kk=1,batch_n
         k=k_batch(kk)
         wg_k=wg(k)
         sqrt_wg_k_arr(kk)=2*sqrt(wg_k*0.5d0)
      end do
      
      allocate(X_0(points_c))
      allocate(phases(points_c))
      allocate(JT(nq,nq))
      JT=transpose(J)
      i_thr=1
      ovcomb=trtype==1 .or. trtype==2
      !$OMP PARALLEL DEFAULT(NONE) &
      !$OMP SHARED(points_c,dt,wg,we,J,JT,J_nz,K_dusch,gamma_gr,sqrt_wg_k_arr) &
      !$OMP SHARED(tmss) &
      !$OMP SHARED(nq,v,kgk,groundleft,groundright,ht,ht2,gamma,theta,z_lwork,tmexp,ovcomb) &
      !$OMP SHARED(X_ap,X_G,X_Gc,X_A,X_Ac,X_0,n_thr,detG,chunk,w_ad) &
      !$OMP SHARED(phases,fixphase) &
      !$OMP SHARED(batch_n,k_batch,nz,sparse,gg,ee,st,correctPhaseX,t_off,trType,contr) &
      !$OMP PRIVATE(c_i,d_i,aa_ex,vp,kgkp,vdv) &
      !$OMP PRIVATE(kk,wg_k,sqrt_wg_k) &
      !$OMP PRIVATE(eta,dx0_dv_,d2x0_dvdv_,dx0_dc_,d3x0_dvdvdv_,d4x0_dvdvdvdv_) &
      !$OMP PRIVATE(d2x0_dvdc_1,d2x0_dvdc_2,d2x0_dvdc_3,d2x0_dcdc_) &
      !$OMP PRIVATE(d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6) &
      !$OMP PRIVATE(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5) &
      !$OMP PRIVATE(d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10) &
      !$OMP PRIVATE(d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5) &
      !$OMP PRIVATE(d5x0_dvdvdvdvdv_) &
      !$OMP PRIVATE(qx_l,xq_l,qxq_lm,qqx_lm,xqq_lm,qqxq_mln,qxqq_mln,qqqx_mln,qqqxq_lmnp,qqxqq_lmnp,qqqxqq_lmnpq) &
      !$OMP PRIVATE(kk_idx) &
      !$OMP PRIVATE(X0,X_fc,X_fc_ht,X_ht_fc,X_ht_ht,x_fc_ht2,x_ht2_fc,x_ht_ht2,x_ht2_ht,x_ht2_ht2) &
      !$OMP PRIVATE(i,jj,t,tau,damp,shift,i_thr,phase) &
      !$OMP PRIVATE(cur_phase,cur_sizee) &
      !$OMP PRIVATE(k,l,m,n,p,q)
      !$ i_thr=omp_get_thread_num()+1
      
      allocate(aa_ex(nq),C_i(nq,nq),D_i(nq,nq))
      allocate(eta(nq))
      allocate(vp(nq))
      
      allocate(x_fc(batch_n))
      
      if(ht)then
         allocate(x_fc_ht(batch_n),x_ht_fc(batch_n),x_ht_ht(batch_n))
         if(ht2)then
            allocate(x_fc_ht2(batch_n),x_ht2_fc(batch_n),x_ht_ht2(batch_n),x_ht2_ht(batch_n),x_ht2_ht2(batch_n))
         end if
      end if
      
      x_fc=0d0
      vp=v
      kgkp=kgk
      vdv=0
      phase=0
      !$OMP DO SCHEDULE(STATIC,chunk)
      do i = 1,points_c
         t=dt*(i-1_8+t_off) !workaround against t=0 (division by zero)
         call x0_acd(nq,J,wg,we,t,aa_ex,C_i,D_i,vp,kgkp) !Here C_i and D_i are not inverted yet
         x0=make_x0(i,nq,wg,aa_ex,vp,kgkp,C_i,D_i,detG,phase) !Here C_i and D_i are inverted
         X_0(i)=x0
      end do
      !$OMP END DO
      
      if(fixphase)then
         !$OMP MASTER
         call FixPhase_X0(X_0,points_c,phases)
         !$OMP END MASTER
      end if
      phase=0
      !$OMP DO SCHEDULE(STATIC,chunk)
      do i = 1,points_c
         t=dt*(i-1_8+t_off) !workaround against t=0 (division by zero)
         tau=iu*t
         damp=exp(-gamma*t - theta**2*t**2*0.5d0) 
         shift=exp(-tau*w_ad)
         
         call x0_acd(nq,J,wg,we,t,aa_ex,C_i,D_i,vp,kgkp) !Here C_i and D_i are not inverted yet
         x0=make_x0(i,nq,wg,aa_ex,vp,kgkp,C_i,D_i,detG,phase) !Here C_i and D_i are inverted
         if(fixphase)then
            cur_sizee=sqrt(realpart(x0)**2+imagpart(x0)**2)
            cur_phase=phases(i)
            x0=cur_sizee*(cos(cur_phase)+iu*sin(cur_phase))
         end if
         X_0(i)=x0*shift*damp
         do l = 1,nq
            eta(l)=0d0
            do jj = 1,nq
               eta(l)=eta(l)+2*D_i(l,jj)*vp(jj)
            end do
         end do
         
         x_fc=0d0
         do l = 1,nq
            qx_l=-0.5d0*x0*dx0_dv(nq,eta,l)
            xq_l=qx_l
            do kk = 1,batch_n
               k=k_batch(kk)
               x_fc(kk)=x_fc(kk)+J(k,l)*qx_l
               if(.not.ht)cycle
               x_fc_ht(kk)=K_dusch(k)*xq_l
               x_ht_fc(kk)=K_dusch(k)*qx_l
            end do
            
            
            if(.not.ht)cycle
            do m = 1,nq
               d2x0_dvdv_=d2x0_dvdv(nq,eta,D_i,m,l)
               dx0_dc_=dx0_dc(nq,C_i,m,l)
               
               qxq_lm=x0*(0.25d0*d2x0_dvdv_+dx0_dc_)
               qqx_lm=x0*(0.25d0*d2x0_dvdv_-dx0_dc_)
               xqq_lm=x0*(0.25d0*d2x0_dvdv_-dx0_dc_)
               do kk = 1,batch_n
                  k=k_batch(kk)
                  x_fc_ht(kk)=x_fc_ht(kk)+J(k,m)*qxq_lm
                  x_ht_fc(kk)=x_ht_fc(kk)+J(k,m)*qqx_lm
                  x_ht_ht(kk)=K_dusch(k)*qxq_lm
                  if(.not.ht2)cycle
                  x_fc_ht2(kk)=K_dusch(k)*xqq_lm
                  x_ht2_fc(kk)=K_dusch(k)*qqx_lm
               end do
               do n = 1,nq
                  d3x0_dvdvdv_=d3x0_dvdvdv(nq,eta,d_i,n,m,l)
                  
                  d2x0_dvdc_1=d2x0_dvdc(nq,C_i,eta,n,m,l)
                  d2x0_dvdc_2=d2x0_dvdc(nq,C_i,eta,l,n,m)
                  d2x0_dvdc_3=d2x0_dvdc(nq,C_i,eta,m,l,n)
                  
                  qqxq_mln=x0*(-0.125d0*d3x0_dvdvdv_-0.5d0*d2x0_dvdc_1-0.5d0*d2x0_dvdc_3+0.5d0*d2x0_dvdc_2)
                  qxqq_mln=x0*(-0.125d0*d3x0_dvdvdv_+0.5d0*d2x0_dvdc_1-0.5d0*d2x0_dvdc_3-0.5d0*d2x0_dvdc_2)
                  qqqx_mln=x0*(-0.125d0*d3x0_dvdvdv_+0.5d0*d2x0_dvdc_1+0.5d0*d2x0_dvdc_2+0.5d0*d2x0_dvdc_3)
                  do kk = 1,batch_n
                     k=k_batch(kk)
                     x_ht_ht(kk)=x_ht_ht(kk)+J(k,n)*qqxq_mln
                     if(.not.ht2)cycle
                     x_fc_ht2(kk)=x_fc_ht2(kk)+J(k,n)*qxqq_mln
                     x_ht2_fc(kk)=x_ht2_fc(kk)+J(k,n)*qqqx_mln
                     X_ht_ht2(kk)=K_dusch(kk)*qxqq_mln
                     X_ht2_ht(kk)=K_dusch(kk)*qqxq_mln
                  end do
                  if(.not.ht2)cycle
                  do p = 1,nq
                                                            !l,m,n,p  old
                                                            !p,n,m,l  new
                     d4x0_dvdvdvdv_=d4x0_dvdvdvdv(nq,eta,D_i,p,n,m,l)
                                                             !k,l,m,n
                   ! d3x0_dvdvdc_1=d3x0_dvdvdc(nq,C_i,eta,D_i,l,m,n,p)
                     d3x0_dvdvdc_1=d3x0_dvdvdc(nq,C_i,eta,D_i,p,n,m,l)
                     ! d3x0_dvdvdc_2=d3x0_dvdvdc(nq,C_i,eta,D_i,l,n,m,p)
                     d3x0_dvdvdc_2=d3x0_dvdvdc(nq,C_i,eta,D_i,p,m,n,l)
                     ! d3x0_dvdvdc_3=d3x0_dvdvdc(nq,C_i,eta,D_i,l,p,m,n)
                     d3x0_dvdvdc_3=d3x0_dvdvdc(nq,C_i,eta,D_i,p,l,n,m)
                     ! d3x0_dvdvdc_4=d3x0_dvdvdc(nq,C_i,eta,D_i,m,n,l,p)
                     d3x0_dvdvdc_4=d3x0_dvdvdc(nq,C_i,eta,D_i,n,m,p,l)
                     ! d3x0_dvdvdc_5=d3x0_dvdvdc(nq,C_i,eta,D_i,m,p,l,n)
                     d3x0_dvdvdc_5=d3x0_dvdvdc(nq,C_i,eta,D_i,n,l,p,m)
                     d3x0_dvdvdc_6=d3x0_dvdvdc(nq,C_i,eta,D_i,m,l,p,n)
                     
                     d2x0_dcdc_=d2x0_dcdc(nq,C_i,p,n,m,l)
                     
                     qqxqq_lmnp=Calc_QQxQQ(x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_)
                     qqqxq_lmnp=Calc_QQQxQ(x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_)
                     do kk = 1,batch_n
                        k=k_batch(kk)
                        X_ht2_ht2(kk)=K_dusch(k)*qqxqq_lmnp
                        X_ht_ht2(kk)=X_ht_ht2(kk)+J(k,p)*qqxqq_lmnp
                        X_ht2_ht(kk)=X_ht2_ht(kk)+J(k,p)*qqqxq_lmnp
                     end do
                     do q = 1,nq
                        ! d5x0_dvdvdvdvdv_=d5x0_dvdvdvdvdv(nq,eta,D_i,l,m,n,p,q)
                        d5x0_dvdvdvdvdv_=d5x0_dvdvdvdvdv(nq,eta,D_i,q,p,n,m,l)
                        
                        call PopDer(nq,C_i,D_i,eta,d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
                         d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
                         d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,q,p,n,m,l)
                         
                        qqqxqq_lmnpq=Calc_QQQxQQ(x0,d5x0_dvdvdvdvdv_,d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
                         d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
                         d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5)
                        
                        do kk = 1,batch_n
                           k=k_batch(kk)
                           X_ht2_ht2(kk)=X_ht2_ht2(kk)+J(k,q)*qqqxqq_lmnpq
                        end do
                     end do
                     do kk = 1,batch_n
                        k=k_batch(kk)
                        call X_tensors_add(x_ht2_ht2(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,9,p,n,m,l)
                     end do
                  end do
                  do kk = 1,batch_n
                     call X_tensors_add(x_ht_ht2(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,7,n,m,l,0)
                     call X_tensors_add(x_ht2_ht(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,8,n,m,l,0)
                  end do
               end do
               do kk = 1,batch_n
                  k=k_batch(kk)
                  call X_tensors_add(x_ht_ht(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,4,m,l,0,0)
                  if(.not.ht2)cycle
                  call X_tensors_add(x_fc_ht2(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,5,m,l,0,0)
                  call X_tensors_add(x_ht2_fc(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,6,m,l,0,0)
               end do
            end do
            do kk = 1,batch_n
               k=k_batch(kk)
               call X_tensors_add(X_fc_ht(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,2,l,0,0,0)
               call X_tensors_add(X_ht_fc(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,3,l,0,0,0)
            end do
         end do
         do kk = 1,batch_n
            k=k_batch(kk)
            x_fc(kk)=x_fc(kk)+K_dusch(k)*x0
            call X_tensors_add(X_fc(kk),i,kk,k,batch_n,points_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,1,0,0,0,0)
         end do
         
            ! qx(k)=-0.5d0*x0*dx0_dv(nq,eta,k)
            ! xq(k)=qx(k)
            ! do l=k,nq
               ! d2x0_dvdv_=d2x0_dvdv(nq,eta,D_i,k,l)!(eta(l)*eta(k)+2*D_i(l,k))
               ! dx0_dc_=dx0_dc(nq,C_i,k,l)
               
               ! qxq(k,l)=x0*(0.25d0*d2x0_dvdv_+dx0_dc_)
               ! qxq(l,k)=qxq(k,l)

               ! qqx(k,l)=x0*(0.25d0*d2x0_dvdv_-dx0_dc_)
               ! qqx(l,k)=qqx(k,l)
               
               ! xqq(k,l)=x0*(0.25d0*d2x0_dvdv_-dx0_dc_)
               ! xqq(l,k)=xqq(k,l)
               ! do m=l,nq
                  ! d3x0_dvdvdv_=d3x0_dvdvdv(nq,eta,d_i,k,l,m)
                  
                  ! d2x0_dvdc_1=d2x0_dvdc(nq,C_i,eta,k,l,m)
                  ! d2x0_dvdc_2=d2x0_dvdc(nq,C_i,eta,m,k,l)
                  ! d2x0_dvdc_3=d2x0_dvdc(nq,C_i,eta,l,m,k)
                  
                  ! qqxq(k,l,m)=x0*(-0.125d0*d3x0_dvdvdv_-0.5d0*d2x0_dvdc_1-0.5d0*d2x0_dvdc_3+0.5d0*d2x0_dvdc_2)
                  ! qqxq(m,k,l)=x0*(-0.125d0*d3x0_dvdvdv_-0.5d0*d2x0_dvdc_2-0.5d0*d2x0_dvdc_1+0.5d0*d2x0_dvdc_3)
                  ! qqxq(l,m,k)=x0*(-0.125d0*d3x0_dvdvdv_-0.5d0*d2x0_dvdc_3-0.5d0*d2x0_dvdc_2+0.5d0*d2x0_dvdc_1)
                  ! qqxq(k,m,l)=qqxq(m,k,l)
                  ! qqxq(m,l,k)=qqxq(l,m,k)
                  ! qqxq(l,k,m)=qqxq(k,l,m)
                  ! if(.not. ht2)cycle
                  
                  ! qxqq(k,l,m)=x0*(-0.125d0*d3x0_dvdvdv_+0.5d0*d2x0_dvdc_1-0.5d0*d2x0_dvdc_3-0.5d0*d2x0_dvdc_2)
                  ! qxqq(m,k,l)=x0*(-0.125d0*d3x0_dvdvdv_+0.5d0*d2x0_dvdc_2-0.5d0*d2x0_dvdc_1-0.5d0*d2x0_dvdc_3)
                  ! qxqq(l,m,k)=x0*(-0.125d0*d3x0_dvdvdv_+0.5d0*d2x0_dvdc_3-0.5d0*d2x0_dvdc_2-0.5d0*d2x0_dvdc_1)
                  ! qxqq(k,m,l)=qxqq(k,l,m)
                  ! qxqq(m,l,k)=qxqq(m,k,l)
                  ! qxqq(l,k,m)=qxqq(l,m,k)
                  
                  ! qqqx(k,l,m)=x0*(-0.125d0*d3x0_dvdvdv_+0.5d0*d2x0_dvdc_1+0.5d0*d2x0_dvdc_2+0.5d0*d2x0_dvdc_3)
                  ! qqqx(m,k,l)=qqqx(k,l,m)
                  ! qqqx(l,m,k)=qqqx(k,l,m)
                  ! qqqx(k,m,l)=qqqx(k,l,m)
                  ! qqqx(m,l,k)=qqqx(k,l,m)
                  ! qqqx(l,k,m)=qqqx(k,l,m)
                  
                  ! do n = m,nq
                     ! d4x0_dvdvdvdv_=d4x0_dvdvdvdv(nq,eta,D_i,k,l,m,n)
                     
                     ! d3x0_dvdvdc_1=d3x0_dvdvdc(nq,C_i,eta,D_i,k,l,m,n)
                     ! d3x0_dvdvdc_2=d3x0_dvdvdc(nq,C_i,eta,D_i,k,m,l,n)
                     ! d3x0_dvdvdc_3=d3x0_dvdvdc(nq,C_i,eta,D_i,k,n,l,m)
                     ! d3x0_dvdvdc_4=d3x0_dvdvdc(nq,C_i,eta,D_i,l,m,k,n)
                     ! d3x0_dvdvdc_5=d3x0_dvdvdc(nq,C_i,eta,D_i,l,n,k,m)
                     ! d3x0_dvdvdc_6=d3x0_dvdvdc(nq,C_i,eta,D_i,m,n,k,l)
                     
                     ! d2x0_dcdc_=d2x0_dcdc(nq,C_i,k,l,m,n)
                     
                     ! call Populate_QQxQQ_and_QQQxQ(x0,nq,qqxqq,qqqxq,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_,k,l,m,n)
                     ! do p = n,nq
                        ! call Populate_QQQxQQ(x0,nq,qqqxqq,C_i,eta,D_i,k,l,m,n,p)
                     ! end do
                  ! end do
               ! end do
            ! end do
         
         if(i_thr==1)then
            if(i==NINT(0.1*chunk))then
               write(output_unit,'(A4)',advance='no')'10%'
               flush(output_unit)
            elseif(i==NINT(0.2*chunk))then
               write(output_unit,'(A4)',advance='no')'20%'
               flush(output_unit)
            elseif(i==NINT(0.3*chunk))then
               write(output_unit,'(A4)',advance='no')'30%'
               flush(output_unit)
            elseif(i==NINT(0.4*chunk))then
               write(output_unit,'(A4)',advance='no')'40%'
               flush(output_unit)
            elseif(i==NINT(0.5*chunk))then
               write(output_unit,'(A4)',advance='no')'50%'
               flush(output_unit)
            elseif(i==NINT(0.6*chunk))then
               write(output_unit,'(A4)',advance='no')'60%'
               flush(output_unit)
            elseif(i==NINT(0.7*chunk))then
               write(output_unit,'(A4)',advance='no')'70%'
               flush(output_unit)
            elseif(i==NINT(0.8*chunk))then
               write(output_unit,'(A4)',advance='no')'80%'
               flush(output_unit)
            elseif(i==NINT(0.9*chunk))then
               write(output_unit,'(A4)',advance='no')'90%'
               flush(output_unit)
            end if
         end if
         do kk = 1,batch_n
            k=k_batch(kk)
            wg_k=sqrt_wg_k_arr(kk)
            X_ap(:,:,kk,i)=X_ap(:,:,kk,i)*wg_k*damp*shift
            X_G(:,:,kk,i)=X_G(:,:,kk,i)*wg_k*iu*damp*shift
            X_Gc(:,:,kk,i)=X_Gc(:,:,kk,i)*wg_k*iu*damp*shift
            X_A(:,:,:,kk,i)=X_A(:,:,:,kk,i)*wg_k*damp*shift
            X_Ac(:,:,:,kk,i)=X_Ac(:,:,:,kk,i)*wg_k*damp*shift
         end do
      end do
      !$OMP END DO
      deallocate(aa_ex,C_i,D_i)
      deallocate(eta)
      deallocate(vp)
      deallocate(x_fc)
      if(ht)then
         deallocate(x_fc_ht,x_ht_fc,x_ht_ht)
         if(ht2)then
            deallocate(x_fc_ht2,x_ht2_fc,x_ht_ht2,x_ht2_ht,x_ht2_ht2)
         end if
      end if
      
      !$OMP END PARALLEL
      deallocate(phases)
      write(output_unit,*)'100%'
   end subroutine Make_corrf_k
   
   pure subroutine X_tensors_add(X_inc,i,k_idx,kk,batch_n,p_c,X_ap,X_G,X_Gc,X_a,X_Ac,tmss,typee,k,l,m,n)
      integer,intent(in) :: typee,l,k,m,n
      integer,intent(in) :: k_idx,kk,batch_n
      integer(int64),intent(in) :: p_c,i
      type(tms),intent(in) :: tmss
      double complex,intent(inout) :: X_ap(3,3,batch_n,p_c),X_G(3,3,batch_n,p_c),X_Gc(3,3,batch_n,p_c),X_A(3,3,3,batch_n,p_c),X_Ac(3,3,3,batch_n,p_c)
      double complex,intent(in) :: X_inc
      integer :: a,b,c
      
      select case(typee)
         case(1)!FC
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%u(a),tmss%u(b),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%u(a),tmss%m(b),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%m(a),tmss%u(b),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%u(a),tmss%q(b,c),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%q(b,c),tmss%u(a),X_inc)
                  end do
               end do
            end do
         case(2)!FC_HT
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%u(a),tmss%du(b,k),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%u(a),tmss%dm(b,k),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%m(a),tmss%du(b,k),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%u(a),tmss%dq(b,c,k),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%q(b,c),tmss%du(a,k),X_inc)
                  end do
               end do
            end do
         case(3)!HT_FC
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%du(a,k),tmss%u(b),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%du(a,k),tmss%m(b),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%dm(a,k),tmss%u(b),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%du(a,k),tmss%q(b,c),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%dq(b,c,k),tmss%u(a),X_inc)
                  end do
               end do
            end do
         case(4)!HT_HT
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%du(a,k),tmss%du(b,l),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%du(a,k),tmss%dm(b,l),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%dm(a,k),tmss%du(b,l),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%du(a,k),tmss%dq(b,c,l),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%dq(b,c,k),tmss%du(a,l),X_inc)
                  end do
               end do
            end do
         case(5)!FC_HT2
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%u(a),tmss%du2(b,k,l),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%u(a),tmss%dm2(b,k,l),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%m(a),tmss%du2(b,k,l),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%u(a),tmss%dq2(b,c,k,l),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%q(b,c),tmss%du2(a,k,l),X_inc)
                  end do
               end do
            end do
         case(6)!HT2_FC
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%du2(a,k,l),tmss%u(b),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%du2(a,k,l),tmss%m(b),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%dm2(a,k,l),tmss%u(b),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%du2(a,k,l),tmss%q(b,c),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%dq2(b,c,k,l),tmss%u(a),X_inc)
                  end do
               end do
            end do
         case(7)!HT_HT2
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%du(a,k),tmss%du2(b,l,m),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%du(a,k),tmss%dm2(b,l,m),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%dm(a,k),tmss%du2(b,l,m),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%du(a,k),tmss%dq2(b,c,l,m),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%dq(b,c,k),tmss%du2(a,l,m),X_inc)
                  end do
               end do
            end do
         case(8)!HT2_HT
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%du2(a,k,l),tmss%du(b,m),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%du2(a,k,l),tmss%dm(b,m),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%dm2(a,k,l),tmss%du(b,m),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%du2(a,k,l),tmss%dq(b,c,m),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%dq2(b,c,k,l),tmss%du(a,m),X_inc)
                  end do
               end do
            end do
         case(9)!HT2_HT2
            do a = 1,3
               do b = 1,3
                  call X_FC_add(X_ap(a,b,k_idx,i),tmss%du2(a,k,l),tmss%du2(b,m,n),X_inc)
                  call X_FC_add(X_G(a,b,k_idx,i),tmss%du2(a,k,l),tmss%dm2(b,m,n),X_inc)
                  call X_FC_add(X_Gc(a,b,k_idx,i),-tmss%dm2(a,k,l),tmss%du2(b,m,n),X_inc)
                  do c = 1,3
                     call X_FC_add(X_a(a,b,c,k_idx,i),tmss%du2(a,k,l),tmss%dq2(b,c,m,n),X_inc)
                     call X_FC_add(X_ac(a,b,c,k_idx,i),tmss%dq2(b,c,k,l),tmss%du2(a,m,n),X_inc)
                  end do
               end do
            end do
      end select
   end subroutine X_tensors_add
   
   pure subroutine X_FC_add(X_el,tm_A,tm_B,X_fc)
      double complex,intent(inout) :: X_el
      double precision,intent(in) :: tm_A,tm_B
      double complex,intent(in) :: X_FC
      
      X_el=X_el+tm_A*tm_B*X_fc
   end subroutine X_FC_add
   
   pure function Calc_QQxQQ(x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_)result(res)
      double complex,intent(in) :: x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_
      double complex :: res
      res=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_1-d3x0_dvdvdc_2-d3x0_dvdvdc_3-d3x0_dvdvdc_4-d3x0_dvdvdc_5+d3x0_dvdvdc_6)+d2x0_dcdc_)
   end function Calc_QQxQQ
   
   pure function Calc_QQQxQ(x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_)result(res)
      double complex,intent(in) :: x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_
      double complex :: res
      res=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_2+d3x0_dvdvdc_4-d3x0_dvdvdc_6+d3x0_dvdvdc_1-d3x0_dvdvdc_3-d3x0_dvdvdc_5)-d2x0_dcdc_)
   end function Calc_QQQxQ
   
   pure subroutine Populate_QQxQQ_and_QQQxQ(x0,nq,qqxqq,qqqxq,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_,a,b,c,d)
      integer,intent(in) :: nq,a,b,c,d
      double complex,intent(in) :: x0,d4x0_dvdvdvdv_,d3x0_dvdvdc_1,d3x0_dvdvdc_2,d3x0_dvdvdc_3,d3x0_dvdvdc_4,d3x0_dvdvdc_5,d3x0_dvdvdc_6,d2x0_dcdc_
      double complex,intent(inout) :: qqxqq(nq,nq,nq,nq),qqqxq(nq,nq,nq,nq)
      integer :: k,l,m,n
      double complex vall
      
      k=a
      l=b
      m=c
      n=d
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_1-d3x0_dvdvdc_2-d3x0_dvdvdc_3-d3x0_dvdvdc_4-d3x0_dvdvdc_5+d3x0_dvdvdc_6)+d2x0_dcdc_)
      !Switching either on the left or on the right
      qqxqq(k,l,m,n)=vall
      qqxqq(l,k,m,n)=vall
      qqxqq(k,l,n,m)=vall
      qqxqq(l,k,n,m)=vall
      
      !Switching left and right Qs
      qqxqq(m,n,k,l)=vall
      qqxqq(n,m,k,l)=vall
      qqxqq(m,n,l,k)=vall
      qqxqq(n,m,l,k)=vall
      
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_4-d3x0_dvdvdc_5-d3x0_dvdvdc_1-d3x0_dvdvdc_6-d3x0_dvdvdc_2+d3x0_dvdvdc_3)+d2x0_dcdc_)
      qqxqq(l,m,n,k)=vall
      qqxqq(m,l,n,k)=vall
      qqxqq(l,m,k,n)=vall
      qqxqq(m,l,k,n)=vall
      
      qqxqq(n,k,l,m)=vall
      qqxqq(k,n,l,m)=vall
      qqxqq(n,k,m,l)=vall
      qqxqq(k,n,m,l)=vall
      
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_2-d3x0_dvdvdc_1-d3x0_dvdvdc_3-d3x0_dvdvdc_4-d3x0_dvdvdc_6+d3x0_dvdvdc_5)+d2x0_dcdc_)
      qqxqq(k,m,l,n)=vall
      qqxqq(k,m,n,l)=vall
      qqxqq(m,k,l,n)=vall
      qqxqq(m,k,n,l)=vall
      
      qqxqq(l,n,k,m)=vall
      qqxqq(n,l,k,m)=vall
      qqxqq(l,n,m,k)=vall
      qqxqq(n,l,m,k)=vall
      
      
      
      
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_2+d3x0_dvdvdc_4-d3x0_dvdvdc_6+d3x0_dvdvdc_1-d3x0_dvdvdc_3-d3x0_dvdvdc_5)-d2x0_dcdc_)
      qqqxq(k,l,m,n)=vall
      qqqxq(k,m,l,n)=vall
      qqqxq(l,k,m,n)=vall
      qqqxq(l,m,k,n)=vall
      qqqxq(m,k,l,n)=vall
      qqqxq(m,l,k,n)=vall
      
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_3+d3x0_dvdvdc_5-d3x0_dvdvdc_6+d3x0_dvdvdc_1-d3x0_dvdvdc_2-d3x0_dvdvdc_4)-d2x0_dcdc_)
      qqqxq(k,l,n,m)=vall
      qqqxq(k,n,l,m)=vall
      qqqxq(n,k,l,m)=vall
      qqqxq(n,l,k,m)=vall
      qqqxq(l,n,k,m)=vall
      qqqxq(l,k,n,m)=vall
      
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_2+d3x0_dvdvdc_6-d3x0_dvdvdc_4+d3x0_dvdvdc_3-d3x0_dvdvdc_1-d3x0_dvdvdc_5)-d2x0_dcdc_)
      qqqxq(k,n,m,l)=vall
      qqqxq(k,m,n,l)=vall
      qqqxq(m,k,n,l)=vall
      qqqxq(m,n,k,l)=vall
      qqqxq(n,k,m,l)=vall
      qqqxq(n,m,k,l)=vall
      
      vall=x0*(1d0/16d0*d4x0_dvdvdvdv_+0.25d0*(d3x0_dvdvdc_6+d3x0_dvdvdc_4-d3x0_dvdvdc_2+d3x0_dvdvdc_5-d3x0_dvdvdc_3-d3x0_dvdvdc_1)-d2x0_dcdc_)
      qqqxq(n,l,m,k)=vall
      qqqxq(n,m,l,k)=vall
      qqqxq(m,l,n,k)=vall
      qqqxq(m,n,l,k)=vall
      qqqxq(l,m,n,k)=vall
      qqqxq(l,n,m,k)=vall
   end subroutine Populate_QQxQQ_and_QQQxQ
   
   pure function Calc_QQQxQQ(x0,d5x0_dvdvdvdvdv_,d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,&
    d4x0_dvdvdvdc_5,d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10,d3x0_dvdcdc_1, &
    d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5)result(res)
      double complex,intent(in) :: x0,d5x0_dvdvdvdvdv_
      double complex,intent(in) :: d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5,d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10
      double complex,intent(in) :: d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5
      double complex res
      res=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5 &
       +d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+ d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10) &
       +0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
   end function Calc_QQQxQQ
   
   pure subroutine Populate_QQQxQQ(x0,nq,qqqxqq,Ci,eta,Di,k,l,m,n,p)
      integer,intent(in) :: nq,k,l,m,n,p
      double complex,intent(in) :: x0,Ci(nq,nq),Di(nq,nq),eta(nq)
      double complex,intent(inout) :: qqqxqq(nq,nq,nq,nq,nq)
      double complex :: d5x0_dvdvdvdvdv_
      double complex :: d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5,d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10
      double complex :: d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5
      double complex vall
      
      ! d5x0_dvdvdvdvdv_=d5x0_dvdvdvdvdv(nq,eta,Di,k,l,m,n,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,l,m,n,p)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,k,l,m,n,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,l,n,m,p)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,k,l,n,m,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,l,p,m,n)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,k,l,p,m,n)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,m,n,l,p)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,k,m,n,l,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,m,p,l,n)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,k,m,p,l,n)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,n,p,l,m)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,k,n,p,l,m)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,l,m,n,k,p)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,l,m,n,k,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,l,m,p,k,n)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,l,m,n,k,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,l,n,p,k,m)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,l,m,n,k,p)
      
      ! call PopDer(d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
       ! d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
       ! d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,m,n,p,k,l)
      ! vall=x0*(-1d0/32d0*d5x0_dvdvdvdvdv_+1d0/8d0*(d4x0_dvdvdvdc_1-d4x0_dvdvdvdc_2-d4x0_dvdvdvdc_3-d4x0_dvdvdvdc_4-d4x0_dvdvdvdc_5+d4x0_dvdvdvdc_6-d4x0_dvdvdvdc_7-d4x0_dvdvdvdc_8+d4x0_dvdvdvdc_9+d4x0_dvdvdvdc_10)+0.5d0*(-d3x0_dvdcdc_1-d3x0_dvdcdc_2-d3x0_dvdcdc_3+d3x0_dvdcdc_4+d3x0_dvdcdc_5))
      ! call symmetrize(nq,qqqxqq,vall,l,m,n,k,p)
      
      contains
      
      
      pure subroutine symmetrize(nq,qqqxqq,vall,k,l,m,n,p)
         integer,intent(in) :: k,l,m,n,p,nq
         double complex,intent(in) :: vall
         double complex,intent(inout) :: qqqxqq(nq,nq,nq,nq,nq)
         
         qqqxqq(k,l,m,n,p)=vall
         qqqxqq(k,m,l,n,p)=vall
         qqqxqq(l,k,m,n,p)=vall
         qqqxqq(l,m,k,n,p)=vall
         qqqxqq(m,l,k,n,p)=vall
         qqqxqq(m,k,l,n,p)=vall
         
         qqqxqq(l,k,m,p,n)=vall
         qqqxqq(l,m,k,p,n)=vall
         qqqxqq(k,m,l,p,n)=vall
         qqqxqq(k,l,m,p,n)=vall
         qqqxqq(m,l,k,p,n)=vall
         qqqxqq(m,k,l,p,n)=vall
      end subroutine symmetrize
      
   end subroutine Populate_QQQxQQ
   
   pure subroutine PopDer(nq,ci,di,eta,d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5, &
    d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10, &
    d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5,k,l,m,n,p)
      integer,intent(in) :: k,l,m,n,p,nq
      double complex,intent(in) :: ci(nq,nq),di(nq,nq),eta(nq)
      double complex,intent(out) :: d4x0_dvdvdvdc_1,d4x0_dvdvdvdc_2,d4x0_dvdvdvdc_3,d4x0_dvdvdvdc_4,d4x0_dvdvdvdc_5
      double complex,intent(out) :: d4x0_dvdvdvdc_6,d4x0_dvdvdvdc_7,d4x0_dvdvdvdc_8,d4x0_dvdvdvdc_9,d4x0_dvdvdvdc_10
      double complex,intent(out) :: d3x0_dvdcdc_1,d3x0_dvdcdc_2,d3x0_dvdcdc_3,d3x0_dvdcdc_4,d3x0_dvdcdc_5
      
      d4x0_dvdvdvdc_1=d4x0_dvdvdvdc(nq,Ci,eta,Di,k,l,m,n,p)
      d4x0_dvdvdvdc_2=d4x0_dvdvdvdc(nq,Ci,eta,Di,k,m,n,l,p)
      d4x0_dvdvdvdc_3=d4x0_dvdvdvdc(nq,Ci,eta,Di,k,m,p,l,n)
      d4x0_dvdvdvdc_4=d4x0_dvdvdvdc(nq,Ci,eta,Di,l,m,n,k,p)
      d4x0_dvdvdvdc_5=d4x0_dvdvdvdc(nq,Ci,eta,Di,l,m,p,k,n)
      d4x0_dvdvdvdc_6=d4x0_dvdvdvdc(nq,Ci,eta,Di,m,n,p,k,l)
      d4x0_dvdvdvdc_7=d4x0_dvdvdvdc(nq,Ci,eta,Di,k,l,n,m,p)
      d4x0_dvdvdvdc_8=d4x0_dvdvdvdc(nq,Ci,eta,Di,k,l,p,m,n)
      d4x0_dvdvdvdc_9=d4x0_dvdvdvdc(nq,Ci,eta,Di,k,n,p,l,m)
      d4x0_dvdvdvdc_10=d4x0_dvdvdvdc(nq,Ci,eta,Di,l,n,p,k,m)
      
      d3x0_dvdcdc_1=d3x0_dvdcdc(nq,Ci,eta,Di,m,k,l,n,p)
      d3x0_dvdcdc_2=d3x0_dvdcdc(nq,Ci,eta,Di,k,l,m,n,p)
      d3x0_dvdcdc_3=d3x0_dvdcdc(nq,Ci,eta,Di,l,k,m,n,p)
      d3x0_dvdcdc_4=d3x0_dvdcdc(nq,Ci,eta,Di,n,k,l,m,p)
      d3x0_dvdcdc_5=d3x0_dvdcdc(nq,Ci,eta,Di,p,k,l,m,n)
   end subroutine PopDer
   
   pure function dx0_dv_an(vk,wgk,wek,t)result(res)
      double complex,intent(in) :: vk
      double precision,intent(in) :: wgk,wek,t
      double complex res
      
      res=(2d0)*vk/(wek*tanh(wek*iu*t)+wgk)
   end function dx0_dv_an
   
   pure function d2x0_dvdv_an(k,l,v,Di,nq)result(res)
      integer,intent(in) :: k,l,nq
      double complex,intent(in) :: v(nq),Di(nq,nq)
      double complex res
      
      res=Di(l,l)*(conjg(v(l))+v(l))*Di(k,k)*(conjg(v(k))+v(k))+Di(k,l)+Di(l,k)
   end function d2x0_dvdv_an
   
   pure function d3x0_dvdvdv_an(k,l,j,v,Di,nq)result(res)
      integer,intent(in) :: k,l,j,nq
      double complex,intent(in) :: v(nq),Di(nq,nq)
      double complex res
      
      res=Di(k,k)*(Di(l,j)+Di(j,l))*(conjg(v(k))+v(k))+ &
          Di(l,l)*(Di(k,j)+Di(j,k))*(conjg(v(l))+v(l))+ &
          Di(j,j)*(Di(k,l)+Di(l,k))*(conjg(v(j))+v(j))+ &
          Di(l,l)*Di(k,k)*Di(j,j)*(conjg(v(l))+v(l))*(conjg(v(k))+v(k))*(conjg(v(j))+v(j))
   end function d3x0_dvdvdv_an
   
   pure function d2x0_dCdv_an(vj,Ci_kl,Di_jj)result(res)
      double complex,intent(in) :: vj,Ci_kl,Di_jj
      double complex res
      
      res=(-0.5d0*Ci_kl*(conjg(vj)*Di_jj+vj*Di_jj))
   end function d2x0_dCdv_an
   
   pure function dx0_dc(nq,Ci,i,j)result(res)
      integer,intent(in) :: nq,i,j
      double complex,intent(in) :: Ci(nq,nq)
      double complex res
      res=-1d0/2d0*Ci(i,j)
   end function dx0_dc
   
   pure function d2x0_dcdc(nq,Ci,i,j,k,l)result(res)
      double complex,intent(in) :: Ci(nq,nq)
      integer,intent(in) :: i,j,k,l,nq
      double complex :: res
      
      res=1d0/4d0*(Ci(i,j)*Ci(k,l)+Ci(i,l)*Ci(j,k)+Ci(i,k)*Ci(j,l))
   end function d2x0_dcdc
   
   pure function d2x0_dvdc(nq,Ci,eta,i,j,k)result(res)
      integer,intent(in) :: nq,i,j,k
      double complex,intent(in) :: Ci(nq,nq),eta(nq)
      double complex res
      res=dx0_dv(nq,eta,i)*dx0_dc(nq,Ci,j,k)
   end function d2x0_dvdc
   
   pure function d3x0_dvdvdc(nq,Ci,eta,Di,i,j,k,l)result(res)
      integer,intent(in) :: nq,i,j,k,l
      double complex,intent(in) :: Ci(nq,nq),eta(nq),Di(nq,nq)
      double complex res
      res=d2x0_dvdv(nq,eta,Di,i,j)*dx0_dc(nq,Ci,k,l)
   end function d3x0_dvdvdc
   
   pure function d3x0_dvdcdc(nq,Ci,eta,Di,i,j,k,l,m)result(res)
      integer,intent(in) :: nq,i,j,k,l,m
      double complex,intent(in) :: Ci(nq,nq),eta(nq),Di(nq,nq)
      double complex res
      res=dx0_dv(nq,eta,i)*d2x0_dcdc(nq,Ci,j,k,l,m)
   end function d3x0_dvdcdc
   
   pure function d4x0_dvdvdvdc(nq,Ci,eta,Di,i,j,k,l,m)result(res)
      integer,intent(in) :: nq,i,j,k,l,m
      double complex,intent(in) :: Ci(nq,nq),eta(nq),Di(nq,nq)
      double complex res
      res=d3x0_dvdvdv(nq,eta,Di,i,j,k)*dx0_dc(nq,Ci,l,m)
   end function d4x0_dvdvdvdc
   
   pure function dx0_dv(nq,eta,i)result(res)
      integer,intent(in) :: i,nq
      double complex,intent(in) :: eta(nq)
      double complex res
      
      res=eta(i)
   end function dx0_dv
   
   pure function d2x0_dvdv(nq,eta,Di,i,j)result(res)
      integer,intent(in) :: i,j,nq
      double complex,intent(in) :: Di(nq,nq),eta(nq)
      double complex res
      
      res=eta(i)*eta(j)+2d0*Di(i,j)
   end function d2x0_dvdv
   
   pure function d3x0_dvdvdv(nq,eta,Di,i,j,k)result(res)
      integer,intent(in) :: i,j,k,nq
      double complex,intent(in) :: Di(nq,nq),eta(nq)
      double complex res
      
      res=eta(i)*eta(j)*eta(k)+2d0*(Di(i,j)*eta(k)+Di(i,k)*eta(j)+Di(j,k)*eta(i))
   end function d3x0_dvdvdv
   
   pure function d4x0_dvdvdvdv(nq,eta,Di,i,j,k,l)result(res)
      double complex,intent(in) :: Di(nq,nq),eta(nq)
      integer,intent(in) :: i,j,k,l,nq
      double complex :: res
      
      res=eta2(i,j)*eta2(k,l)+2d0*(Di(i,j)*eta2(k,l)+Di(i,l)*eta2(k,j)+Di(i,k)*eta2(j,l)+Di(j,k)*eta2(i,l)+Di(j,l)*eta2(k,i)+Di(k,l)*eta2(i,j)) &
       + 4d0*(Di2(i,j,k,l)+Di2(i,k,l,j)+Di2(i,l,j,k))
      contains
      
      pure function Di2(i,j,k,l)result(res2)
         integer,intent(in) :: i,j,k,l
         double complex :: res2
         res2=Di(i,j)*Di(k,l)
      end function Di2
      
      pure function eta2(i,j)result(res2)
         integer,intent(in) :: i,j
         double complex :: res2
         res2=eta(i)*eta(j)
      end function eta2
      
   end function d4x0_dvdvdvdv
   
   pure function d5x0_dvdvdvdvdv(nq,eta,Di,i,j,k,l,m)result(res)
      double complex,intent(in) :: Di(nq,nq),eta(nq)
      integer,intent(in) :: i,j,k,l,m,nq
      double complex :: res
      
      res=eta(i)*eta(j)*eta(k)*eta(l)*eta(m) &
       +2d0*(Di(i,j)*eta3(k,l,m)+Di(i,k)*eta3(j,l,m)+Di(i,l)*eta3(j,k,m) &
       +Di(i,m)*eta3(j,k,l)+Di(j,k)*eta3(i,l,m)+Di(j,l)*eta3(i,k,m)+Di(j,m)*eta3(i,k,l) &
       +Di(k,l)*eta3(i,j,m)+Di(k,m)*eta3(i,j,l)+Di(l,m)*eta3(i,j,k)) &
       +4d0*eta(m)*(Di2(i,j,k,l)+Di2(i,k,j,l)+Di2(i,l,j,k)) &
       +4d0*eta(l)*(Di2(i,j,k,m)+Di2(i,k,j,m)+Di2(i,m,j,k)) &
       +4d0*eta(k)*(Di2(i,j,l,m)+Di2(i,l,j,m)+Di2(i,m,j,l)) &
       +4d0*eta(j)*(Di2(i,k,l,m)+Di2(i,l,k,m)+Di2(i,m,k,l)) &
       +4d0*eta(i)*(Di2(j,k,l,m)+Di2(j,l,k,m)+Di2(j,m,k,l))
      contains
       
      pure function Di2(i,j,k,l)result(res2)
         integer,intent(in) :: i,j,k,l
         double complex :: res2
         res2=Di(i,j)*Di(k,l)
      end function Di2
       
      pure function eta3(i,j,k)result(res2)
         integer,intent(in) :: i,j,k
         double complex :: res2
         res2=eta(i)*eta(j)*eta(k)
      end function eta3
   end function d5x0_dvdvdvdvdv
   
   subroutine FixPhase_X0(X0_arr,points_c,phases)
      integer(int64) points_c,i
      double complex X0_arr(points_c)
      double precision phase,phasePrev,sizee,sizeePrev,phases(points_c)
      double precision im,re,tolsin,tolcos,cosres,sinres
      
      re=realpart(X0_arr(1))
      im=imagpart(X0_arr(1))
      phasePrev=atan2(im,re)
      phases(1)=phasePrev
      sizeePrev=sqrt(re**2+im**2)
      
      re=realpart(X0_arr(2))
      im=imagpart(X0_arr(2))
      phasePrev=atan2(im,re)
      phases(2)=phasePrev
      sizeePrev=sqrt(re**2+im**2)
      
      ! tolsin=(sin(phases(2))-sin(phases(1)))*1.5d0
      ! tolcos=(cos(phases(2))-cos(phases(1)))*1.5d0
      tolsin=0.1
      tolcos=0.1
      do i = 3,points_c
         re=realpart(X0_arr(i))
         im=imagpart(X0_arr(i))
         phase=atan2(im,re)
         sizee=sqrt(re**2+im**2)
         cosres=cos(phase)-cos(phaseprev)
         sinres=sin(phase)-sin(phaseprev)
         if(abs(cosres)>tolcos .or. abs(sinres)>tolsin)then
            phase=PhaseNorm(phase+pi)
            X0_arr(i)=sizee*(cos(phase)+iu*sin(phase))
         end if
         sizeePrev=sizee
         phasePrev=phase
         phases(i)=phase
      end do
   end subroutine FixPhase_X0
   
   pure subroutine PermSym_3(vall,arr,n,i,j,k)
      integer,intent(in) :: n,i,j,k
      double complex,intent(inout) :: arr(n,n,n)
      double complex,intent(in) :: vall
      
      arr(i,j,k)=vall
      arr(j,k,i)=vall
      arr(k,i,j)=vall
      
      arr(j,i,k)=vall
      arr(i,k,j)=vall
      arr(k,j,i)=vall
   end subroutine PermSym_3
   
   pure subroutine PermSym_4(vall,arr,n,i,j,k,l)
      integer,intent(in) :: n,i,j,k,l
      double complex,intent(inout) :: arr(n,n,n,n)
      double complex,intent(in) :: vall
      
      arr(i,j,k,l)=vall
      arr(j,k,l,i)=vall
      arr(k,l,i,j)=vall
      arr(l,i,j,k)=vall
      
      arr(j,i,k,l)=vall
      arr(i,k,l,j)=vall
      ! arr()=vall
      ! arr()=vall
      ! arr()=vall
      ! arr()=vall
      ! arr()=vall
   end subroutine PermSym_4
   
   pure subroutine RepairPhaseWrap(X,points_c)
      integer(int64),intent(in) :: points_c
      double complex,intent(inout) :: X(points_c)
      
      
      call PhaseUnwrapping(X,points_c)
      ! X_r=realpart(X)
      ! coeff_r=0
      ! if(.not.all(X_r==0d0))then
         ! call NormalizeCorrF(X_r,points_c,coeff_r)
         ! if(abss)then
            ! call RepairPhaseCorrF_abs(X_r,points_c,dt,t_off,thr)
         ! else
            ! call RepairPhaseCorrF(X_r,points_c,dt,t_off,thr)
         ! end if
      ! end if
      
      
      ! X_i=imagpart(X)
      ! coeff_i=0
      ! if(.not.all(X_i==0d0))then
         ! call NormalizeCorrF(X_i,points_c,coeff_i)
         ! if(abss)then
            ! call RepairPhaseCorrF_abs(X_i,points_c,dt,t_off,thr)
         ! else
            ! call RepairPhaseCorrF(X_i,points_c,dt,t_off,thr)
         ! end if
      ! end if
      ! X=X_r*coeff_r + iu*X_i*coeff_i
   end subroutine RepairPhaseWrap
   
   pure subroutine PhaseUnwrapping(X,n)
      integer(int64),intent(in) :: n
      integer(int64) i
      double complex,intent(inout) :: X(n)
      double precision r,phi,phi_prev,dphi,dphi_cur
      
      call co_polar(X(1),r,phi_prev)
      call co_polar(X(2),r,phi)
      dphi=phi-phi_prev
      if(dphi>pi)dphi=dphi-2*pi
      phi_prev=phi
      do i = 3,n
         call co_polar(X(i),r,phi)
         dphi_cur=phi-phi_prev
         if(dphi_cur>pi)then
            dphi_cur=dphi_cur-2*pi
         end if
         if(abs(dphi_cur-dphi)>1d-4)then
            phi_prev=PhaseNorm(phi_prev+dphi)
            X(i)=r*(cos(phi_prev)+iu*sin(phi_prev))
         else
            phi_prev=phi
         end if
      end do
   end subroutine PhaseUnwrapping
   
   pure subroutine NormalizeCorrF(X_dbl,points_c,coeff)
      integer(int64),intent(in) :: points_c
      double precision,intent(inout) :: X_dbl(points_c)
      double precision,intent(out) :: coeff
      double precision maxx,minn
      
      maxx=maxval(X_dbl,dim=1)
      minn=minval(X_dbl,dim=1)
      coeff=max(abs(maxx),abs(minn))
      X_dbl=X_dbl/coeff
   end subroutine NormalizeCorrF
      
   pure subroutine RepairPhaseCorrF(X,points_c,dt,t_off,thr)
      integer(int64),intent(in) :: points_c
      double precision,intent(inout) :: X(points_c)
      double precision,intent(in) :: dt,t_off,thr
      
      double precision t,X_abs,dx_dt,dx_dt_prev
      integer(int64) :: i
      
      do i = 5,points_c
         dx_dt=Der3AtPoint_forward(X,points_c,i,dt)
         dx_dt_prev=Der3AtPoint_forward(X,points_c,i-1,dt)
         if(abs(abs(dx_dt)-abs(dx_dt_prev))>thr)then
            X(i)=-X(i)
         end if
      end do
   end subroutine RepairPhaseCorrF
   
   pure subroutine RepairPhaseCorrF_Abs(X,points_c,dt,t_off,thr)
      integer(int64),intent(in) :: points_c
      double precision,intent(inout) :: X(points_c)
      double precision,intent(in) :: dt,t_off,thr
      
      double precision,allocatable :: X_buf(:)
      double precision t,X_abs,dx_dt,dx_dt_prev
      integer(int64) :: i
      
      
      X_buf=abs(X)
      do i = 1,points_c
         if(X_buf(i)-X(i)>thr/1d6)then
            X(i)=-X_buf(i)
         else
            X(i)=X_buf(i)
         end if
      end do
      deallocate(X_buf)
   end subroutine RepairPhaseCorrF_Abs
   
   pure function Der3AtPoint_forward(X,n,i,h)result(res)
      integer(int64),intent(in) :: n,i
      double precision,intent(in) :: X(n)
      double precision,intent(in) :: h
      double precision res
      
      res=(X(i)-3*X(i-1)+3*X(i-2)-X(i-3))/h**3
   end function Der3AtPoint_forward
   
   pure function Der2AtPoint_forward(X,n,i,h)result(res)
      integer(int64),intent(in) :: n,i
      double precision,intent(in) :: X(n)
      double precision,intent(in) :: h
      double precision res
      
      res=(X(i)-2*X(i-1)+X(i-2))/h**2
   end function Der2AtPoint_forward
   
   pure function Der1AtPoint_forward(X,n,i,h)result(res)
      integer(int64),intent(in) :: n,i
      double precision,intent(in) :: X(n)
      double precision,intent(in) :: h
      double precision res
      
      res=(X(i)-X(i-1))/h
   end function Der1AtPoint_forward
   
   pure subroutine Make_X_fcht12HT(ht2,nq,i,points_c,sqrt_wg_k,kk,J,JT,K,qxqq,qqxqq,qqqx,qqqxq,qqqxqq,qqxq,xqq,qqx,qxq,qx,xq,qqjxqj,qjxqj,qqjx,qxqj,qjx,xqj,x0,&
      tmexp,sparse,nz,J_nz,u,m,q,du,dm,dq,u_ex,m_ex,q_ex,du_ex,dm_ex,dq_ex,du2_ex,dm2_ex,dq2_ex,batch_n,kk_idx,contr,X_Ap,X_G,X_Gc,X_A,X_Ac)
      integer,intent(in) :: nq,kk,nz,kk_idx,batch_n
      type(list_arr_int32),intent(in) :: J_nz(nq)
      integer(int64),intent(in) :: points_c,i
      integer l,lp,lpp,jj,mm,kkk,nn,pp,ii,ll
      integer li,lpi,lppi
      integer a,b,c,d,e,f
      integer(int64) idx,newidx
      logical,intent(in) :: tmexp(2),sparse,contr(9),ht2
      double precision,intent(in) :: sqrt_wg_k,J(nq,nq),JT(nq,nq),K(nq)
      double complex,intent(in) :: qqxq(nq,nq,nq),qqx(nq,nq),qxq(nq,nq),qx(nq),xq(nq),x0
      double complex,intent(in) :: qxqq(:,:,:),qqxqq(:,:,:,:),qqqx(:,:,:),qqqxq(:,:,:,:),qqqxqq(:,:,:,:,:),xqq(:,:)
      double complex,intent(in) :: qqjxqj(:,:,:),qjxqj(:,:),qqjx(:,:),qxqj(:,:),qjx(:),xqj(:)
      double complex :: X_fcht2_inc,X_fcht1_inc,Ksum,X_ht_inc
      double complex :: X_fc_ht2,X_ht_ht2,X_ht2_fc,X_ht2_ht,X_ht2_ht2
      
      
      double precision,intent(in) :: u(3),m(3),q(3,3)
      double precision,intent(in) :: du(3,nq),dm(3,nq),dq(3,3,nq)
      
      double precision,intent(in) :: u_ex(3),m_ex(3),q_ex(3,3)
      double precision,intent(in) :: du_ex(3,nq),dm_ex(3,nq),dq_ex(3,3,nq)
      double precision,intent(in) :: du2_ex(3,nq,nq),dm2_ex(3,nq,nq),dq2_ex(3,3,nq,nq)
      
      double complex,intent(inout) :: X_ap(3,3,batch_n,points_c),X_G(3,3,batch_n,points_c),X_Gc(3,3,batch_n,points_c),X_A(3,3,3,batch_n,points_c),X_Ac(3,3,3,batch_n,points_c)
      
      
      
      if(tmexp(1).and.tmexp(2))then !GG
         do jj = 1,nq
            X_fcht1_inc=qjxqj(kk,jj)
            X_fcht2_inc=K(kk)*qjx(jj)
            do mm = 1,nq
               !X_fcht1_inc=X_fcht1_inc+J(kk,mm)*qxqj(mm,jj)
               X_fcht2_inc=X_fcht2_inc+JT(mm,kk)*qqjx(mm,jj)
               X_ht_inc=K(kk)*qjxqj(mm,jj)
               do nn = 1,nq
                  X_ht_inc=X_ht_inc+JT(nn,kk)*qqjxqj(nn,mm,jj)
               end do
               do a = 1,3
                  do b = 1,3
                     X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+X_ht_inc*du(a,jj)*du(b,mm)
                     X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+X_ht_inc*du(a,jj)*dm(b,mm)*iu
                     X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-X_ht_inc*(dm(a,jj))*du(b,mm)*iu
                     do c = 1,3
                        X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+X_ht_inc*du(a,jj)*dq(b,c,mm)
                        X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+X_ht_inc*dq(b,c,jj)*du(a,mm)
                     end do
                  end do
               end do
            end do
            !X_fcht1=X_fcht1+dtm_B(jj)*X_fcht1_inc
            !X_fcht2=X_fcht2+dtm_A(jj)*X_fcht2_inc
            !TODO check
            do a = 1,3
               do b = 1,3
                  X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+(X_fcht1_inc*u(a)*du(b,jj)+du(a,jj)*u(b)*X_fcht2_inc)
                  X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+(X_fcht1_inc*u(a)*dm(b,jj)+du(a,jj)*m(b)*X_fcht2_inc)*iu
                  X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-(X_fcht1_inc*(m(a))*du(b,jj)+(dm(a,jj))*u(b)*X_fcht2_inc)*iu
                  do c = 1,3
                     X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+(X_fcht1_inc*u(a)*dq(b,c,jj)+X_fcht2_inc*du(a,jj)*q(b,c))
                     X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+(X_fcht1_inc*q(b,c)*du(a,jj)+X_fcht2_inc*dq(b,c,jj)*u(a))
                  end do
               end do
            end do
         end do
      else if(.not.tmexp(1) .and. .not.tmexp(2))then!EE
         do ll=1,nq
            X_fcht1_inc=K(kk)*xq(ll)
            X_fcht2_inc=K(kk)*qx(ll)
            do mm=1,nq
               X_fcht1_inc=X_fcht1_inc+JT(mm,kk)*qxq(mm,ll)
               x_fcht2_inc=X_fcht2_inc+JT(mm,kk)*qqx(mm,ll)
               x_ht_inc=K(kk)*qxq(mm,ll)
               
               if(ht2)X_fc_ht2=K(kk)*xqq(mm,ll)
               if(ht2)X_ht2_fc=K(kk)*qqx(mm,ll)
               
               do nn=1,nq
                  x_ht_inc=x_ht_inc+JT(nn,kk)*qqxq(nn,mm,ll)
                  
                  if(.not.ht2)cycle
                  x_fc_ht2=x_fc_ht2+JT(nn,kk)*qxqq(nn,mm,ll)
                  x_ht2_fc=x_ht2_fc+JT(nn,kk)*qqqx(nn,mm,ll)
                  X_ht_ht2=K(kk)*qxqq(nn,mm,ll)
                  X_ht2_ht=K(kk)*qqxq(nn,mm,ll)
                  
                  do pp=1,nq
                     X_ht_ht2=X_ht_ht2+JT(pp,kk)*qqxqq(pp,nn,mm,ll)
                     X_ht2_ht=X_ht2_ht+JT(pp,kk)*qqqxq(pp,nn,mm,ll)
                  
                     X_ht2_ht2=K(kk)*qqxqq(ll,mm,nn,pp)
                     do ii = 1,nq
                        X_ht2_ht2=X_ht2_ht2+JT(ii,kk)*qqqxqq(ii,pp,nn,mm,ll)
                     end do
                  end do
               end do
               
               do a = 1,3
                  do b = 1,3
                     X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+X_ht_inc*du_ex(a,ll)*du_ex(b,mm)
                     X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+X_ht_inc*du_ex(a,ll)*dm_ex(b,mm)*iu
                     X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-X_ht_inc*(dm_ex(a,ll))*du_ex(b,mm)*iu
                     
                     if(ht2)X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+X_fc_ht2*u_ex(a)*du2_ex(b,mm,ll)
                     if(ht2)X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+X_fc_ht2*u_ex(a)*dm2_ex(b,mm,ll)*iu
                     if(ht2)X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-X_fc_ht2*m_ex(a)*du2_ex(b,mm,ll)*iu
                     
                     if(ht2)X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+X_ht2_fc*du2_ex(a,mm,ll)*u_ex(b)
                     if(ht2)X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+X_ht2_fc*du2_ex(a,mm,ll)*m_ex(b)*iu
                     if(ht2)X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-X_ht2_fc*dm2_ex(a,mm,ll)*u_ex(b)*iu
                     
                     do c = 1,3
                        X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+X_ht_inc*du_ex(a,ll)*dq_ex(b,c,mm)
                        X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+X_ht_inc*dq_ex(b,c,ll)*du_ex(a,mm)
                        
                        if(ht2)X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+X_fc_ht2*u_ex(a)*dq2_ex(b,c,mm,ll)
                        if(ht2)X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+X_fc_ht2*q_ex(b,c)*du2_ex(a,mm,ll)
                        
                        if(ht2)X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+X_ht2_fc*du2_ex(a,mm,ll)*q_ex(b,c)
                        if(ht2)X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+X_ht2_fc*dq2_ex(b,c,mm,ll)*u_ex(a)
                     end do
                  end do
               end do
            end do

            do a = 1,3
               do b = 1,3
                  X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+(X_fcht1_inc*u_ex(a)*du_ex(b,ll))
                  X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+(X_fcht1_inc*u_ex(a)*dm_ex(b,ll))*iu
                  X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-(X_fcht1_inc*(m_ex(a))*du_ex(b,ll))*iu
                  do c = 1,3
                     X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+(X_fcht1_inc*u_ex(a)*dq_ex(b,c,ll))
                     X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+(X_fcht1_inc*q_ex(b,c)*du_ex(a,ll))
                  end do
               end do
            end do
            do a = 1,3
               do b = 1,3
                  X_ap(a,b,kk_idx,i)=X_ap(a,b,kk_idx,i)+(du_ex(a,ll)*u_ex(b)*X_fcht2_inc)
                  X_G(a,b,kk_idx,i)=X_G(a,b,kk_idx,i)+(du_ex(a,ll)*m_ex(b)*X_fcht2_inc)*iu
                  X_Gc(a,b,kk_idx,i)=X_Gc(a,b,kk_idx,i)-(dm_ex(a,ll)*u_ex(b)*X_fcht2_inc)*iu
                  do c = 1,3
                     X_A(a,b,c,kk_idx,i)=X_A(a,b,c,kk_idx,i)+(X_fcht2_inc*du_ex(a,ll)*q_ex(b,c))
                     X_Ac(a,b,c,kk_idx,i)=X_Ac(a,b,c,kk_idx,i)+(X_fcht2_inc*dq_ex(b,c,ll)*u_ex(a))
                  end do
               end do
            end do
         end do
      end if
   end subroutine Make_X_fcht12HT   
   
   pure subroutine x0_acd(nq,J,wg,we,t,aa_ex,C,D,v,kdk)
      integer,intent(in) :: nq
      double precision,intent(in) :: wg(nq),we(nq),J(nq,nq)
      double precision,intent(in) :: t
      double complex,intent(inout) :: v(nq),kdk
      double complex :: ce,de,tanh_res,bb_ex(nq)
      double complex,intent(out) :: aa_ex(nq),C(nq,nq),D(nq,nq)
      double precision,parameter :: eps = 0.000005
      
      
      integer i,jj,k
      
      C=0
      D=0
      do i = 1,nq
         jj=i
         aa_ex(i)=we(i)/sin(t*we(i)+iu*eps)
         ! bb_ex(i)=we(i)/tan(t*we(i))
         tanh_res=iu*tan(we(i)*t/2d0+iu*eps)

         ce=we(i)/tanh_res
         de=we(i)*tanh_res
         ! ce=-iu*(bb_ex(i)+aa_ex(i))
         ! de=-iu*(bb_ex(i)-aa_ex(i))
         C(i,i)=ce
         D(i,i)=de
         do jj = 1,nq
            do k=1,nq
               C(i,jj)=C(i,jj)+J(k,i)*wg(k)*J(k,jj)
               D(i,jj)=D(i,jj)+J(k,i)*wg(k)*J(k,jj)
            end do
         end do
         
      end do
   end subroutine x0_acd
   
   pure function DetFromLU_C_big(mat,N,piv)result(det)
      integer,intent(in) :: N,piv(n)
      integer detP
      double complex,intent(in) :: mat(N,N)
      type(big_complex) det
      
      ! call zgetrf(N,N,mat,N,Piv,info)
      ! if(info/=0)then
         ! stop 2
      ! end if
      
      call DetPermMat(Piv,N,detP)
      det=DiagProd_C_big(mat,n)*detP
   end function DetFromLU_C_big
   
   pure subroutine DetFromLU_C_big_ln(mat,N,piv,lndet,phi,detP)
      integer,intent(in) :: N,piv(n)
      integer,intent(out) :: detP
      integer i
      double complex,intent(in) :: mat(N,N)
      double precision,intent(out) :: lndet,phi
      
      ! call zgetrf(N,N,mat,N,Piv,info)
      ! if(info/=0)then
         ! stop 2
      ! end if
      
      call DetPermMat(Piv,N,detP)
      
      lndet=0d0
      phi=0d0
      do i = 1,n
         lndet=lndet+log(abs(mat(i,i)))
         phi=phi+atan2(imagpart(mat(i,i)),realpart(mat(i,i)))
      end do
      
      if(detP==-1)phi=phi+pi
      !phi=PhaseNorm(phi)
   end subroutine DetFromLU_C_big_ln
   
   pure function DiagProd_C_big(mat,N)result(res)
      integer,intent(in) :: N
      double complex,intent(in) :: mat(N,N)
      type(big_complex) res
      
      integer i
      
      res=BC((1d0,0d0))
      
      do i = 1,N
         res=res*mat(i,i)
      end do
   end function DiagProd_C_big
   
   pure function DetFromLU_big(mat,N,piv)result(det)
      integer,intent(in) :: N,piv(n)
      integer detP
      double precision,intent(in) :: mat(N,N)
      type(big_double) det
      
      ! call dgetrf(N,N,mat,N,Piv,info)
      ! if(info/=0)then
         ! stop 2
      ! end if
      
      call DetPermMat(Piv,N,detP)
      det=DiagProd_big(mat,n)*detP
   end function DetFromLU_big
   
   pure function DiagProd_big(mat,N)result(res)
      integer,intent(in) :: N
      double precision,intent(in) :: mat(N,N)
      type(big_double) res
      
      integer i
      
      res=BD(1d0)
      
      do i = 1,N
         res=res*mat(i,i)
      end do
   end function DiagProd_big
   
   pure function DiagProd_big_vec(vec,N)result(res)
      integer,intent(in) :: N
      double precision,intent(in) :: vec(n)
      type(big_double) res
      
      integer i
      
      res=bd(1d0)
      
      do i = 1,N
         res=res*vec(i)
      end do
   end function DiagProd_big_vec
   
   pure function DiagProd_C_big_vec(vec,N)result(res)
      integer,intent(in) :: N
      double complex,intent(in) :: vec(n)
      type(big_complex) res
      
      integer i
      
      res=bc((1d0,0d0))
      
      do i = 1,N
         res=res*vec(i)
      end do
   
   end function DiagProd_C_big_vec
   
   pure subroutine DiagProd_C_big_vec_ln(vec,N,ln_r,phi)
      integer,intent(in) :: N
      double complex,intent(in) :: vec(n)
      double precision,intent(out) :: ln_r,phi
      
      integer i
      
      !res=bc((1d0,0d0))
      ln_r=0d0
      phi=0d0
      do i = 1,N
         ln_r=ln_r+log(abs(vec(i)))
         phi=phi+atan2(imagpart(vec(i)),realpart(vec(i)))
      end do
      !phi=PhaseNorm(phi)
   end subroutine DiagProd_C_big_vec_ln
   
   pure function PhaseNorm(phi)result(res)
      double precision,intent(in) :: phi
      double precision res
      
      ! res=atan2(cos(phi),sin(phi))
      
      res=mod(phi+pi,2d0*pi)
      if(res<0)res=res+2*pi
      res=res-pi
      ! res=phi
   end function PhaseNorm
   
   function Make_x0(i,nq,wg,aa_ex,v,kgk,C,D,detG,phasePrev)result(x0)
      integer,intent(in) :: nq
      integer(int64),intent(in) :: i
      integer ierr,piv(nq),ii,jj
      double complex,intent(inout) :: C(nq,nq),D(nq,nq)
      double complex :: bc_ex,vdv
      double complex,intent(in) :: kgk,v(nq)
      double precision,intent(in) :: wg(nq)
      double precision phase
      double complex x0
      double complex aa_ex(nq),mat(nq,nq)
      type(big_double) detG
      type(big_complex) detA,detA_new,detC,detD,detProd
      
      
      double precision lndetG,lndetA,lndetC,lndetD,lndetProd
      double precision phaseG,phaseA,phaseC,phaseD,phaseProd
      double precision phasePrev,dphi,phasesDiff(4)
      integer detPc,detPd
      
      logical,parameter :: pv=.true.
      
      lndetG=(log(detG%num)+detG%p*log(10d0))
      phaseG=0d0 !negative frequencies do not exist in FCOV :))), we will ignore that case completely
      call DiagProd_c_big_vec_ln(aa_ex,nq,lndetA,phaseA)
      lndetA=lndetA+nq*log(2d0)
      ! phaseA=PhaseNorm(phaseA-pi/2d0*nq)
      phaseA=phaseA-pi/2d0*nq
      
      
      call zgetrf(nq,nq,C,nq,Piv,ierr) !LU form
      call DetFromLU_C_big_ln(C,nq,piv,lndetC,phaseC,detPc)
      call INV_C_lapack(C,piv,z_lwork,nq) !inverse
      
      lndetC=lndetC
      phaseC=phaseC
      
      call zgetrf(nq,nq,D,nq,Piv,ierr) !LU form
      call DetFromLU_C_big_ln(D,nq,piv,lndetD,phaseD,detPd)
      call INV_C_lapack(D,piv,z_lwork,nq) !inverse
      
      lndetD=lndetD
      phaseD=phaseD
      
      vdv=(0d0,0d0)
      do ii = 1,nq
         do jj = 1,nq
            vdv=vdv+conjg(v(ii))*D(ii,jj)*v(jj)
         end do
      end do
      
      lndetProd=lndetG/2d0+lndetA/2d0-lndetC/2d0-lndetD/2d0
      phaseProd=phaseG/2d0+phaseA/2d0-phaseC/2d0-phaseD/2d0
      phaseProd=PhaseNorm(phaseProd)
      ! if(abs(cos(phaseprod)-cos(phaseprev))>0.1d0 .or. abs(sin(phaseprod)-sin(phaseprev))>0.1d0)phaseProd=phaseProd+pi
      !phaseProd=PhaseNorm(phaseProd)
      phasePrev=phaseProd
      bc_ex=exp(lndetProd)*(cos(phaseProd)+iu*sin(phaseProd))
      x0=bc_ex*exp(-kgk+vdv)
   end function Make_x0
   
   
   pure function Make_x0_matmul(nq,sc,a,wg,c_i,d_i)result(res)
      integer,intent(in) :: nq
      double complex,intent(in) :: a(nq),c_i(nq,nq),d_i(nq,nq),sc
      double precision,intent(in) :: wg(nq)
      double complex res(nq,nq)
      
      integer i,j,m
      
      do i = 1,nq
         do j = 1,nq
            res(j,i)=(0d0,0d0)
            do m = 1,nq
               res(j,i)=res(j,i)+c_i(j,m)*d_i(m,i)
            end do
            res(j,i)=res(j,i)*a(j)*wg(j)*sc
         end do
      end do
      
   end function Make_x0_matmul
   
end module RROA_TD

!Handles the vibronic RROA calculation, mostly the time-independent part
module RROA
   use constants
   use iso_fortran_env
   use FCOV_storage
   use FCOV_FCFuns
   use OMP_functions
   use stuff
   use strings
   use spectra
   use RROA_TD !handles the time-dependent RROA calculation
   implicit none
   
   type RROA_Options_TI
      !TODO
   end type RROA_Options_TI
   
   integer,parameter :: op_polcontrs_unit=500
   
   contains
   
   !ripped from the EVCD module
   subroutine DoDipRot(dip_12,mag_21,w0,spr,x,npx,dw,dw_which,kbt,lg, &
                        v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,i_thr,thr_D,thr_R,dd)
      logical lg
      integer npx,i,dw_which
      double precision dip_12(3),mag_21(3),diprot(2),w_nm,w_cm
      double precision x(npx),spr(npx,2),kbt,dw,w0,D,R
      integer v1_clas,v1(v1_clas)
      integer v2_clas,v2(v2_clas)
      integer(int16) v1_pos(v1_clas),v2_pos(v2_clas)
      double precision cm(npx),dd,dd_arr(npx),sfc(npx),wm,xx
      double precision,parameter :: D_lim=1d-20,R_lim=1d-22
      
      integer i_thr
      double precision thr_D,thr_R
      
      D=dip_12(1)**2+dip_12(2)**2+dip_12(3)**2 &
       /3d0*AU2SI_debye**2 &
       *108.7d0*kbt
      R=DOT_PRODUCT(dip_12,mag_21) &
       /3d0*2*2.541765*9.2740154d-21*1e-18/1e-36 &
       *435.0d0*kbt
      !R=dip_12(1)*mag_21(1)+dip_12(2)*mag_21(2)+dip_12(3)*mag_21(3)
      w_nm=10d0**7/(w0*au_2_cm)
      
      !D=D/3d0*AU2SI_debye**2
      !taken from guvcd4, R [cgs/10**-36]
      !R=R/3d0*2*2.541765*9.2740154d-21*1e-18/1e-36 
      
      ! if(D>thr_D .or. abs(R)>thr_R)then
         ! curTrans=Trans_Make(w0,D,R,v1,v1_pos,v1_clas,v2,v2_pos,v2_clas)
         ! call TransList_AddTrans(tr_list(i_thr+1),curTrans)
      ! end if
      
      
      !D=108.7d0*D*kbt
      !R=435.0d0*R*kbt
      select case(dw_which) !nm
         case(1)
            wm=w_nm
            !do i = 1,npx
               !xx=x(i)
               ! sfc=sf(dd,LG,x(i),wm)
               ! spr(i,1)=spr(i,1)+sfc*D*wm
               ! spr(i,2)=spr(i,2)+sfc*R*wm
               sfc=sf(dd,LG,x,wm)*wm
               if(D>D_lim)spr(:,1)=spr(:,1)+sfc*D
               if(R>R_lim)spr(:,2)=spr(:,2)+sfc*R
            !end do
         case(2) !cm-1
            wm=w0*au_2_cm
            !do i = 1,npx
               !xx=x(i)
               !cm=1.0d7/x
               sfc=sf(dd,LG,cm,wm)*wm
               if(D>D_lim)spr(:,1)=spr(:,1)+sfc*D
               if(R>R_lim)spr(:,2)=spr(:,2)+sfc*R
            !end do
         case(3) !x in nm, but dw is in energy constant units in cm
            wm=w_nm
            !do i = 1,npx
               !xx=x(i)
               !if(la)then
               dd_arr=dd*x**2/1.0d7
               !end if
               sfc=sf(dd_arr,LG,x,wm)*wm
               if(D>D_lim)spr(:,1)=spr(:,1)+sfc*D
               if(R>R_lim)spr(:,2)=spr(:,2)+sfc*R
            !end do
         end select
   end subroutine DoDipRot
   
   
   !w - excitation frequency in au
   !wr - scattered frequency in au
   !ei - initial state energy in au
   !ef - final state energy in au
   !THRC - if the Raman intensity is below this threshold, do not add ROA
   !nt - transition index (just for some output)
   !si - initial state, mode vector?
   !ni - mode vector length?
   !nf,sf - same but for final state
   !EXCNM -excitation frequency in nm
   !wrin - minimal Raman freq in cm-1
   !wrax - maximal Raman freq in cm-1
   !npx - number of spectral points
   !lglg - Gaussian line shape (else Lorentzian)
   !fwhh - FWHH in cm-1
   !ltab - make RROA.TAB,ECDI.TAB
   !inv - the ROA invariants
   !alr..atci - polarizabilities in order as they usually appear (see Hecht 1991)
   !r - real part, i - imaginary part
   !sr - spectrum
   !bf - Boltzmann factor
   !lusea - use the A,Ac polarizabilities
   !luseg - use the G,Gc polarizabilities
   !ldo - vector of different experimental setups
   !sr_ -
   !si_ - factors of polarizabilities (ie. sr_*<a|u|b>*<b|u|c>)
   subroutine wrram3(w,wR,temp,spectrum_temp,ei,ef,THRC,nt,wrin,wrax,npx,lglg,fwhh,ltab, &
            inv,alr,ali,gtr,gti,gtcr,gtci,atr,ati,atcr,atci,sr,bf,lusea,luseg,ldo,&
            tabOutputs,str_v1,str_v3,alr_fcht,ali_fcht)
      
      
      
!     bf .. Boltzmann factor
      implicit none
      integer*4 ni,nf,nt,npx,ns0,ni0, &
      ia,b,id,e,ii,is
      parameter (ns0=19,ni0=13)
      integer tabOutputs(ns0)
!     ns0 : number of experimental setups
!     ni0 : number of invariants
      real*8 ef,ei,fwhh,THRC,CM,AMU,BOHR,ECM,wrax,wrin,w,clight, &
      YDY,YDX,sr(npx,2,ns0),gpisvejc,roa1,ram1,tr,ti,a(2),bf, &
      inv(ni0),co(ni0+2,ns0),aaar,aaai,gtaar,gtaai,gtcaar,gtcaai, &
      alr(3,3),ali(3,3),gtr(3,3),gti(3,3),gtcr(3,3),gtci(3,3), &
      atr(3,3,3),ati(3,3,3),atcr(3,3,3),atci(3,3,3), &
      alsr(3,3),alsi(3,3),gtsr(3,3),gtsi(3,3),gtcsr(3,3),gtcsi(3,3), &
      alar(3,3),alai(3,3),gtar(3,3),gtai(3,3),gtcar(3,3),gtcai(3,3), &
      ear(3,3),eapr(3,3),easr(3,3),eapsr(3,3),eaar(3,3),eapar(3,3), &
      eai(3,3),eapi(3,3),easi(3,3),eapsi(3,3),eaai(3,3),eapai(3,3), &
      eacr(3,3),eapcr(3,3),eacsr(3,3),eapcsr(3,3),eacar(3,3), &
      eapcar(3,3),sr_,si_, &
      eaci(3,3),eapci(3,3),eacsi(3,3),eapcsi(3,3),eacai(3,3),eapcai(3,3), &
      alr_fcht(3,3),ali_fcht(3,3),alsr_fcht(3,3), alsi_fcht(3,3), alar_fcht(3,3), & 
      alai_fcht(3,3),aaar_fcht,aaai_fcht
      double precision wR,temp,Kcon,e0
      character(20) str_v1,str_v3
      logical ltab,lglg,lusea,luseg,ldo(ns0)
      logical opl,spectrum_temp
!             prefix RAM, prefix ROA
!                 a2  bs(a)2  ba(a)2      aG  bs(G)2 ba(G)2  bs(A)2
!             ba(A)2      aG  bs(G)2  ba(G)2  bs(A)2 ba(A)2
!     script:              ^     ^       ^       ^      ^
!
!    co: Raman prefactor, ROA prefactor,3 Raman invariants,
!        10 ROA invariants, all for ns0=18 spectral kinds
!    1 1 2 ICP(0o): Nafie
      data co/ 4.0d0,  8.0d0, &
             45.0d0,  7.0d0,  5.0d0, 45.0d0,  7.0d0,  5.0d0, 1.0d0, &
             -1.0d0,-45.0d0,  5.0d0, -5.0d0, -3.0d0,  1.0d0,&
!    2 3 4 ICPx(90o): Nafie
              2.0d0,  4.0d0,&
             45.0d0,  7.0d0,  5.0d0, 45.0d0,  7.0d0,  5.0d0,  1.0d0,&
             -1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    3 5 6 ICPz(90o): Nafie
              4.0d0,  8.0d0,&
              0.0d0,  3.0d0,  5.0d0,  0.0d0,  3.0d0,  5.0d0, -1.0d0,&
              1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    4 7 8 ICP*(90o): (magic) Nafie
            6.666d0,  13.3333d0,&
              9.0d0,  2.0d0,  2.0d0,  9.0d0,  2.0d0,  2.0d0,  0.0d0,&
              0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    5 9 10 ICPu(90o): Nafie
              4.0d0,  4.0d0,&
             45.0d0, 13.0d0, 15.0d0, 45.0d0, 13.0d0, 15.0d0, -1.0d0,&
              1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    6 11 12 ICP(180o): Nafie
              4.0d0,  8.0d0, &
             45.0d0,  7.0d0,  5.0d0, 45.0d0,  7.0d0,  5.0d0,  1.0d0, &
             -1.0d0, 45.0d0, -5.0d0,  5.0d0,  3.0d0, -1.0d0, &
!    7 13 14 SCP(0o): nafie
              4.0d0,  8.0d0,&
             45.0d0,  7.0d0,  5.0d0, 45.0d0, -5.0d0,  5.0d0, -3.0d0,&
             -1.0d0,-45.0d0, -7.0d0, -5.0d0,  1.0d0,  1.0d0,&
!    8 15 16 SCPx(90o): nafie
              2.0d0,  4.0d0, &
             45.0d0,  7.0d0,  5.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
              0.0d0,-45.0d0, -7.0d0, -5.0d0,  1.0d0,  1.0d0,&
!    9 17 18 SCPz(90o): Nafie
              4.0d0,  8.0d0,&
              0.0d0,  3.0d0,  5.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
              0.0d0,  0.0d0, -3.0d0, -5.0d0, -1.0d0, -1.0d0,&
!    10 19 20 SCP*(90o): Nafie
            6.666d0,  13.3333d0,&
              9.0d0,  2.0d0,  2.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
              0.0d0, -9.0d0, -2.0d0, -2.0d0,  0.0d0,  0.0d0,&
!    11 21 22 SCPu(90o): Nafie
              4.0d0,  4.0d0, &
             45.0d0, 13.0d0, 15.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
              0.0d0,-45.0d0,-13.0d0,-15.0d0, -1.0d0, -1.0d0,&
!    12 23 24 SCP(180o):Nafie
              4.0d0,  8.0d0, &
             45.0d0,  7.0d0,  5.0d0,-45.0d0,  5.0d0, -5.0d0,  3.0d0,&
              1.0d0,-45.0d0, -7.0d0, -5.0d0,  1.0d0,  1.0d0,& 
!    13 25 26 DCPI(0o): Nafie
              4.0d0,  8.0d0,&
             45.0d0,  1.0d0,  5.0d0, 45.0d0,  1.0d0,  5.0d0, -1.0d0,&
             -1.0d0,-45.0d0, -1.0d0, -5.0d0, -1.0d0, +1.0d0,&
!    14 27 28 DCPI(90o): Nafie
              2.0d0,  2.0d0, &
             45.0d0, 13.0d0, 15.0d0, 45.0d0, 13.0d0, 15.0d0, -1.0d0,&
             +1.0d0,-45.0d0,-13.0d0,-15.0d0, -1.0d0, -1.0d0,&
!    15 29 30 DCPI(180o): Nafie
             24.0d0, 16.0d0,&
              0.0d0,  1.0d0,  0.0d0,  0.0d0,  3.0d0,  0.0d0,  1.0d0,&
              0.0d0,  0.0d0, -3.0d0,  0.0d0,  1.0d0,  0.0d0,&
!    16 31 32 DCPII(0o): Nafie
              24.0d0, 16.0d0,&
               0.0d0,  1.0d0,  0.0d0,  0.0d0,  3.0d0,  0.0d0,  1.0d0,&
               0.0d0,  0.0d0,  3.0d0,  0.0d0, -1.0d0,  0.0d0,&
!    17 33 34 DCPII(90o): 
               2.0d0,  2.0d0,&
              45.0d0, 13.0d0, 15.0d0, 45.0d0, 13.0d0, 15.0d0, -1.0d0,&
               1.0d0, 45.0d0, 13.0d0, 15.0d0,  1.0d0, 1.0d0,&
!    18 35 36 DCPII(180o): 
               4.0d0,  8.0d0,&
              45.0d0,  1.0d0,  5.0d0, 45.0d0,  1.0d0,  5.0d0, -1.0d0,&
              -1.0d0, 45.0d0,  1.0d0,  5.0d0,  1.0d0, -1.0d0,&
!    19 special
               1.0d0,  0.0d0,&
               1.0d0,   .0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
               0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0/
!             prefix      a2  bs(a)2  ba(a)2      aG  bs(G)2 ba(G)2
!             bs(A)2  ba(A)2      aG  bs(G)2  ba(G)2  bs(A)2 ba(A)2
!     script:            ^          ^    ^       ^        ^     ^
!     invariants:
!     1 a2     =(1/9)Re (as_aa as*_bb)
!     2 bs(a)2 =(1/2)Re(3as_ab as*_ab-as_aa as*_bb)
!     3 ba(a)2 =(1/2)Re(3as_ab as*_ab-as_aa as*_bb)
!     4 aG     =(1/9)Im(as_aa Gs*_bb)
!     5 bs(G)2 =(1/2)Im(3as_ab Gs*_ab-as_aa Gs*_bb)
!     6 ba(G)2 =(3/2)Im(3aa_ab Ga*_ab)
!     7 bs(A)2 =(w/2)Im(i as_ab eas*_ab)
!                ea_ab=eps_adg A_d,gb
!     8 ba(A)2 =(w/2)Im(i aa_ab eaa*_ab)+i aa_ab eap_ab
!                eap_ab=eps_abg A_d,gd
!       script tensors:
!     9 aG     =(1/9)Im(as_aa Gs*_bb)
!    10 bs(G)2 =(1/2)Im(3as_ab Gs*_ab-as_aa Gs*_bb)
!    11 ba(G)2 =(3/2)Im(3aa_ab Ga*_ab)
!    12 bs(A)2 =(w/2)Im(i as_ab eas*_ab)
!    13 ba(A)2 =(w/2)Im(i aa_ab eaa*_ab)+i aa_ab eap_ab
      
      CM=219474.63d0
      AMU=1822.0d0
      BOHR=0.529177d0
      ECM=(ef-ei)*CM
      !gpisvejc=(AMU*BOHR**5)*1.0d4*2.0d0*4.0d0*atan(1.0d0)/EXCA
      clight=137.03599d0
      e0=1d0/(4d0*pi)
      
      Kcon=((pi/e0)**2)*((w-ef/(2*pi*cc_AU))**4)/90
      ! Kcon=((pi/e0)**2)*((w-ei/(2*pi*cc_AU))**4)/90
      !if(.not.spectrum_temp)Kcon=Kcon*bf
      Kcon=Kcon*bf
      
!     eap_ab = eps_abc A_d,cd
      call calcep(eapr ,atr )
      call calcep(eapi ,ati )
      call calcep(eapci,atci)
      call calcep(eapcr,atcr)
!      
!     ea_ab = eps_adc A_d,cb
      ear =0.0d0
      eai =0.0d0
      eacr=0.0d0
      eaci=0.0d0
      do ia=1,3
      id=ia+1
      if(id.gt.3)id=1
      e=id+1
      if(e.gt.3)e=1
      do b=1,3
         ear( ia,b)=atr (id,e,b)-atr (e,id,b)
         eai( ia,b)=ati (id,e,b)-ati (e,id,b)
         eacr(ia,b)=atcr(id,e,b)-atcr(e,id,b)
         eaci(ia,b)=atci(id,e,b)-atci(e,id,b)
      end do
      end do

!     symmetric and antisymmetric tensor combinations
      call dsa(alr , ali, alsr, alsi, alar, alai)
      call dsa(alr_fcht , ali_fcht, alsr_fcht, alsi_fcht, alar_fcht, alai_fcht)
      call dsa(gtr , gti, gtsr, gtsi, gtar, gtai)
      call dsa(gtcr,gtci,gtcsr,gtcsi,gtcar,gtcai)
      call dsa(ear  ,eai  ,easr  ,easi  ,eaar  ,eaai  )
      call dsa(eapr ,eapi ,eapsr ,eapsi ,eapar ,eapai )
      call dsa(eacr ,eaci ,eacsr ,eacsi ,eacar ,eacai )
      call dsa(eapcr,eapci,eapcsr,eapcsi,eapcar,eapcai)
      inv=0.0d0
!     a2: (1/9) Re (als_aa als*_bb)
      aaar=alsr(1,1)+alsr(2,2)+alsr(3,3)
      aaai=alsi(1,1)+alsi(2,2)+alsi(3,3)
      
      aaar_fcht=alsr_fcht(1,1)+alsr_fcht(2,2)+alsr_fcht(3,3)
      aaai_fcht=alsi_fcht(1,1)+alsi_fcht(2,2)+alsi_fcht(3,3)
      inv(1)=(aaar_fcht*aaar_fcht+aaai_fcht*aaai_fcht)/9.0d0
!     beta_s(alpha)2:(1/2)Re(3als_ab als*_ab-als_aa als*_bb)
      call abab(tr,ti,alsr_fcht,alsi_fcht,alsr_fcht,alsi_fcht)
      inv(2)=1.5d0*tr-4.5d0*inv(1)
!     beta_a(alpha)2:(3/2)Re(als_ab als*_ab)
      call abab(tr,ti,alar_fcht,alai_fcht,alar_fcht,alai_fcht)
      inv(3)=1.5d0*tr

      if(luseg)then
!      aG: (1/9) Im(als_aa gs*_bb)
       gtaar=gtsr(1,1)+gtsr(2,2)+gtsr(3,3)
       gtaai=gtsi(1,1)+gtsi(2,2)+gtsi(3,3)
       inv(4)=(aaai*gtaar-aaar*gtaai)/9.0d0
!      beta_s(G)2:(1/2)Im(3als_ab Gs*_ab-als_aa Gs*_bb)
       call abab(tr,ti,alsr,alsi,gtsr,gtsi)
       inv(5)=1.5d0*ti-4.5d0*inv(4)
!      beta_a(G)2: (3/2)Im(ala_ab Ga*_ab)
       call abab(tr,ti,alar,alai,gtar,gtai)
       inv(6)=1.5d0*ti
!      aGc:
       gtcaar=gtcsr(1,1)+gtcsr(2,2)+gtcsr(3,3)
       gtcaai=gtcsi(1,1)+gtcsi(2,2)+gtcsi(3,3)
       inv(9)=(aaai*gtcaar-aaar*gtcaai)/9.0d0
!      beta_s(Gc)2:
       call abab(tr,ti,alsr,alsi,gtcsr,gtcsi)
       inv(10)=1.5d0*ti-4.5d0*inv(9)
!      beta_a(Gc)2: (3/2)Im(ala_ab Gca*_ab)
       call abab(tr,ti,alar,alai,gtcar,gtcai)
       inv(11)=1.5d0*ti
      endif

      if(lusea)then
!      betasA2: (w/2)Im(i als_ab e_agd As*g,db)
       call abab(tr,ti,alsr,alsi,easr,easi)
       inv(7)=0.5d0*w*tr/clight
!      betaaA2: (w/2)Im(i[ala_ab[e_agd Aa*g,db+ e_abg Aa*d,gd)
       call abab(tr,ti,alar,alai,eaar,eaai)
       inv(8)=tr
       call abab(tr,ti,alar,alai,eapar,eapai)
       inv(8)=0.5d0*w*(inv(8)+tr)/clight
!      betasAc2: (w/2)Im(i als_ab e_agd Acs*g,db)
       call abab(tr,ti,alsr,alsi,eacsr,eacsi)
       inv(12)=0.5d0*wR*tr/clight
!      betaaAc2:
       call abab(tr,ti,alar,alai,eacar,eacai)
       inv(13)=tr
       call abab(tr,ti,alar,alai,eapcar,eapcai)
       inv(13)=0.5d0*wR*(inv(13)+tr)/clight
      endif
      
      
!     spectral intensities:
      do is=1,ns0
      if(ldo(is))then
!      Raman:
       a(1)= &
       Kcon*co(1,is)*(co(3,is)*inv(1)+co(4,is)*inv(2)+co(5,is)*inv(3))
!      ROA:
       a(2)=0.0d0
       do ii=4,ni0
         a(2)=a(2)+inv(ii)*co(2+ii,is)
       end do
       a(2)=co(2,is)*a(2)*kcon/cc_au
       
       if(dabs(a(1)).gt.THRC)then
        if(is.le.18)call ap3(is,ns0,ECM,a,sr,wrin,wrax,npx,fwhh,lglg,temp,spectrum_temp)
        if(ltab .and. .false.)then
         YDX=0.0d0
         YDY=0.0d0
         ! if(is.eq.19)then
          ! ram1=sr_*sr_+si_*si_
          ! roa1=0.0d0
         ! else
          ! ram1=a(1)*gpisvejc
          ! roa1=a(2)*gpisvejc
          ram1=a(1)
          roa1=a(2)
         ! endif
         WRITE(tabOutputs(is),3001)nt,ECM,YDX,YDY,ram1,roa1,bf,TR(str_v1),TR(str_v3)
3001     FORMAT(I7,f9.2,2f3.0,g12.4,' 0 0 0 0',2g12.4,A,' --> ',A)
!        initial and final states:
         ! if(is.eq.12)then
          ! call wrs(50+is,si,ni)
          ! call wrs(50+is,sf,nf)
         ! endif
         ! write(50+is,*)
        endif
       endif
       end if
      end do
      
      return
   end subroutine wrram3
   
   subroutine wrram3_noSpectra(w,wR,ei,ef,THRC,nt,EXCNM, &
            inv,alr,ali,gtr,gti,gtcr,gtci,atr,ati,atcr,atci,bf,lusea,luseg,ldo,&
            alr_fcht,ali_fcht,tabOutputs,str_v1,str_v3,a_arr)
      
      
      
!     EXCNM ... excitation frequency in nm
!     bf .. Boltzmann factor
      implicit none
      integer*4 ni,nf,nt,npx,ns0,ni0, &
      ia,b,id,e,ii,is
      parameter (ns0=19,ni0=13)
      integer tabOutputs(ns0)
!     ns0 : number of experimental setups
!     ni0 : number of invariants
      real*8 ef,ei,fwhh,THRC,CM,AMU,BOHR,ECM,EXCNM,wrax,wrin,w,clight, &
      YDY,YDX,gpisvejc,EXCA,roa1,ram1,tr,ti,a(2),bf, &
      inv(ni0),a_arr(2,ns0),co(ni0+2,ns0),aaar,aaai,gtaar,gtaai,gtcaar,gtcaai, &
      alr(3,3),ali(3,3),gtr(3,3),gti(3,3),gtcr(3,3),gtci(3,3), &
      atr(3,3,3),ati(3,3,3),atcr(3,3,3),atci(3,3,3), &
      alsr(3,3),alsi(3,3),gtsr(3,3),gtsi(3,3),gtcsr(3,3),gtcsi(3,3), &
      alar(3,3),alai(3,3),gtar(3,3),gtai(3,3),gtcar(3,3),gtcai(3,3), &
      ear(3,3),eapr(3,3),easr(3,3),eapsr(3,3),eaar(3,3),eapar(3,3), &
      eai(3,3),eapi(3,3),easi(3,3),eapsi(3,3),eaai(3,3),eapai(3,3), &
      eacr(3,3),eapcr(3,3),eacsr(3,3),eapcsr(3,3),eacar(3,3), &
      eapcar(3,3),sr_,si_, &
      eaci(3,3),eapci(3,3),eacsi(3,3),eapcsi(3,3),eacai(3,3),eapcai(3,3), &
      alr_fcht(3,3),ali_fcht(3,3),alsr_fcht(3,3), alsi_fcht(3,3), alar_fcht(3,3), & 
      alai_fcht(3,3),aaar_fcht,aaai_fcht
      double precision wR,temp,Kcon,e0
      character(20) str_v1,str_v3
      logical ltab,lglg,lusea,luseg,ldo(ns0)
      logical opl,spectrum_temp
!             prefix RAM, prefix ROA
!                 a2  bs(a)2  ba(a)2      aG  bs(G)2 ba(G)2  bs(A)2
!             ba(A)2      aG  bs(G)2  ba(G)2  bs(A)2 ba(A)2
!     script:              ^     ^       ^       ^      ^
!
!    co: Raman prefactor, ROA prefactor,3 Raman invariants,
!        10 ROA invariants, all for ns0=18 spectral kinds
!    1 1 2 ICP(0o): Nafie
      data co/ 4.0d0,  8.0d0, &
             45.0d0,  7.0d0,  5.0d0, 45.0d0,  7.0d0,  5.0d0, 1.0d0, &
             -1.0d0,-45.0d0,  5.0d0, -5.0d0, -3.0d0,  1.0d0,&
!    2 3 4 ICPx(90o): Nafie
              2.0d0,  4.0d0,&
             45.0d0,  7.0d0,  5.0d0, 45.0d0,  7.0d0,  5.0d0,  1.0d0,&
             -1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    3 5 6 ICPz(90o): Nafie
              4.0d0,  8.0d0,&
              0.0d0,  3.0d0,  5.0d0,  0.0d0,  3.0d0,  5.0d0, -1.0d0,&
              1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    4 7 8 ICP*(90o): (magic) Nafie
            6.666d0,  13.3333d0,&
              9.0d0,  2.0d0,  2.0d0,  9.0d0,  2.0d0,  2.0d0,  0.0d0,&
              0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    5 9 10 ICPu(90o): Nafie
              4.0d0,  4.0d0,&
             45.0d0, 13.0d0, 15.0d0, 45.0d0, 13.0d0, 15.0d0, -1.0d0,&
              1.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
!    6 11 12 ICP(180o): Nafie
              4.0d0,  8.0d0, &
             45.0d0,  7.0d0,  5.0d0, 45.0d0,  7.0d0,  5.0d0,  1.0d0, &
             -1.0d0, 45.0d0, -5.0d0,  5.0d0,  3.0d0, -1.0d0, &
!    7 13 14 SCP(0o): nafie
              4.0d0,  8.0d0,&
             45.0d0,  7.0d0,  5.0d0, 45.0d0, -5.0d0,  5.0d0, -3.0d0,&
             -1.0d0,-45.0d0, -7.0d0, -5.0d0,  1.0d0,  1.0d0,&
!    8 15 16 SCPx(90o): nafie
              2.0d0,  4.0d0, &
             45.0d0,  7.0d0,  5.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0, &
              0.0d0,-45.0d0, -7.0d0, -5.0d0,  1.0d0,  1.0d0,&
!    9 17 18 SCPz(90o): Nafie
              4.0d0,  8.0d0,&
              0.0d0,  3.0d0,  5.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
              0.0d0,  0.0d0, -3.0d0, -5.0d0, -1.0d0, -1.0d0,&
!    10 19 20 SCP*(90o): Nafie
            6.666d0,  13.3333d0,&
              9.0d0,  2.0d0,  2.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
              0.0d0, -9.0d0, -2.0d0, -2.0d0,  0.0d0,  0.0d0,&
!    11 21 22 SCPu(90o): Nafie
              4.0d0,  4.0d0, &
             45.0d0, 13.0d0, 15.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
              0.0d0,-45.0d0,-13.0d0,-15.0d0, -1.0d0, -1.0d0,&
!    12 23 24 SCP(180o):Nafie
              4.0d0,  8.0d0, &
             45.0d0,  7.0d0,  5.0d0,-45.0d0,  5.0d0, -5.0d0,  3.0d0,&
              1.0d0,-45.0d0, -7.0d0, -5.0d0,  1.0d0,  1.0d0,& 
!    13 25 26 DCPI(0o): Nafie
              4.0d0,  8.0d0,&
             45.0d0,  1.0d0,  5.0d0, 45.0d0,  1.0d0,  5.0d0, -1.0d0,&
             -1.0d0,-45.0d0, -1.0d0, -5.0d0, -1.0d0, +1.0d0,&
!    14 27 28 DCPI(90o): Nafie
              2.0d0,  2.0d0, &
             45.0d0, 13.0d0, 15.0d0, 45.0d0, 13.0d0, 15.0d0, -1.0d0,&
             +1.0d0,-45.0d0,-13.0d0,-15.0d0, -1.0d0, -1.0d0,&
!    15 29 30 DCPI(180o): Nafie
             24.0d0, 16.0d0,&
              0.0d0,  1.0d0,  0.0d0,  0.0d0,  3.0d0,  0.0d0,  1.0d0,&
              0.0d0,  0.0d0, -3.0d0,  0.0d0,  1.0d0,  0.0d0,&
!    16 31 32 DCPII(0o): Nafie
              24.0d0, 16.0d0,&
               0.0d0,  1.0d0,  0.0d0,  0.0d0,  3.0d0,  0.0d0,  1.0d0,&
               0.0d0,  0.0d0,  3.0d0,  0.0d0, -1.0d0,  0.0d0,&
!    17 33 34 DCPII(90o): 
               2.0d0,  2.0d0,&
              45.0d0, 13.0d0, 15.0d0, 45.0d0, 13.0d0, 15.0d0, -1.0d0,&
               1.0d0, 45.0d0, 13.0d0, 15.0d0,  1.0d0, 1.0d0,&
!    18 35 36 DCPII(180o): 
               4.0d0,  8.0d0,&
              45.0d0,  1.0d0,  5.0d0, 45.0d0,  1.0d0,  5.0d0, -1.0d0,&
              -1.0d0, 45.0d0,  1.0d0,  5.0d0,  1.0d0, -1.0d0,&
!    19 special
               1.0d0,  0.0d0,&
               1.0d0,   .0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,&
               0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0,  0.0d0/
!             prefix      a2  bs(a)2  ba(a)2      aG  bs(G)2 ba(G)2
!             bs(A)2  ba(A)2      aG  bs(G)2  ba(G)2  bs(A)2 ba(A)2
!     script:            ^          ^    ^       ^        ^     ^
!     invariants:
!     1 a2     =(1/9)Re (as_aa as*_bb)
!     2 bs(a)2 =(1/2)Re(3as_ab as*_ab-as_aa as*_bb)
!     3 ba(a)2 =(1/2)Re(3as_ab as*_ab-as_aa as*_bb)
!     4 aG     =(1/9)Im(as_aa Gs*_bb)
!     5 bs(G)2 =(1/2)Im(3as_ab Gs*_ab-as_aa Gs*_bb)
!     6 ba(G)2 =(3/2)Im(3aa_ab Ga*_ab)
!     7 bs(A)2 =(w/2)Im(i as_ab eas*_ab)
!                ea_ab=eps_adg A_d,gb
!     8 ba(A)2 =(w/2)Im(i aa_ab eaa*_ab)+i aa_ab eap_ab
!                eap_ab=eps_abg A_d,gd
!       script tensors:
!     9 aG     =(1/9)Im(as_aa Gs*_bb)
!    10 bs(G)2 =(1/2)Im(3as_ab Gs*_ab-as_aa Gs*_bb)
!    11 ba(G)2 =(3/2)Im(3aa_ab Ga*_ab)
!    12 bs(A)2 =(w/2)Im(i as_ab eas*_ab)
!    13 ba(A)2 =(w/2)Im(i aa_ab eaa*_ab)+i aa_ab eap_ab
      
      CM=219474.63d0
      AMU=1822.0d0
      BOHR=0.529177d0
      EXCA=EXCNM*10.0d0
      ECM=(ef-ei)*CM
      gpisvejc=(AMU*BOHR**5)*1.0d4*2.0d0*4.0d0*atan(1.0d0)/EXCA
      clight=137.03599d0
      e0=1d0/(4d0*pi)
      
      a_arr=0d0
      
      Kcon=((pi/e0)**2)*((wr/(2*pi*cc_AU))**4)/90
      ! Kcon=((pi/e0)**2)*((w-ei/(2*pi*cc_AU))**4)/90
      if(.not.spectrum_temp)Kcon=Kcon*bf
      
!     eap_ab = eps_abc A_d,cd
      call calcep(eapr ,atr )
      call calcep(eapi ,ati )
      call calcep(eapci,atci)
      call calcep(eapcr,atcr)
!      
!     ea_ab = eps_adc A_d,cb
      ear =0.0d0
      eai =0.0d0
      eacr=0.0d0
      eaci=0.0d0
      do ia=1,3
      id=ia+1
      if(id.gt.3)id=1
      e=id+1
      if(e.gt.3)e=1
      do b=1,3
         ear( ia,b)=atr (id,e,b)-atr (e,id,b)
         eai( ia,b)=ati (id,e,b)-ati (e,id,b)
         eacr(ia,b)=atcr(id,e,b)-atcr(e,id,b)
         eaci(ia,b)=atci(id,e,b)-atci(e,id,b)
      end do
      end do

!     symmetric and antisymmetric tensor combinations
      call dsa(alr , ali, alsr, alsi, alar, alai)
      call dsa(alr_fcht , ali_fcht, alsr_fcht, alsi_fcht, alar_fcht, alai_fcht)
      call dsa(gtr , gti, gtsr, gtsi, gtar, gtai)
      call dsa(gtcr,gtci,gtcsr,gtcsi,gtcar,gtcai)
      call dsa(ear  ,eai  ,easr  ,easi  ,eaar  ,eaai  )
      call dsa(eapr ,eapi ,eapsr ,eapsi ,eapar ,eapai )
      call dsa(eacr ,eaci ,eacsr ,eacsi ,eacar ,eacai )
      call dsa(eapcr,eapci,eapcsr,eapcsi,eapcar,eapcai)
      inv=0.0d0
!     a2: (1/9) Re (als_aa als*_bb)
      aaar=alsr(1,1)+alsr(2,2)+alsr(3,3)
      aaai=alsi(1,1)+alsi(2,2)+alsi(3,3)
      
      aaar_fcht=alsr_fcht(1,1)+alsr_fcht(2,2)+alsr_fcht(3,3)
      aaai_fcht=alsi_fcht(1,1)+alsi_fcht(2,2)+alsi_fcht(3,3)
      inv(1)=(aaar_fcht*aaar_fcht+aaai_fcht*aaai_fcht)/9.0d0
!     beta_s(alpha)2:(1/2)Re(3als_ab als*_ab-als_aa als*_bb)
      call abab(tr,ti,alsr_fcht,alsi_fcht,alsr_fcht,alsi_fcht)
      inv(2)=1.5d0*tr-4.5d0*inv(1)
!     beta_a(alpha)2:(3/2)Re(als_ab als*_ab)
      call abab(tr,ti,alar_fcht,alai_fcht,alar_fcht,alai_fcht)
      inv(3)=1.5d0*tr

      if(luseg)then
!      aG: (1/9) Im(als_aa gs*_bb)
       gtaar=gtsr(1,1)+gtsr(2,2)+gtsr(3,3)
       gtaai=gtsi(1,1)+gtsi(2,2)+gtsi(3,3)
       inv(4)=(aaai*gtaar-aaar*gtaai)/9.0d0
!      beta_s(G)2:(1/2)Im(3als_ab Gs*_ab-als_aa Gs*_bb)
       call abab(tr,ti,alsr,alsi,gtsr,gtsi)
       inv(5)=1.5d0*ti-4.5d0*inv(4)
!      beta_a(G)2: (3/2)Im(ala_ab Ga*_ab)
       call abab(tr,ti,alar,alai,gtar,gtai)
       inv(6)=1.5d0*ti
!      aGc:
       gtcaar=gtcsr(1,1)+gtcsr(2,2)+gtcsr(3,3)
       gtcaai=gtcsi(1,1)+gtcsi(2,2)+gtcsi(3,3)
       inv(9)=(aaai*gtcaar-aaar*gtcaai)/9.0d0
!      beta_s(Gc)2:
       call abab(tr,ti,alsr,alsi,gtcsr,gtcsi)
       inv(10)=1.5d0*ti-4.5d0*inv(9)
!      beta_a(Gc)2: (3/2)Im(ala_ab Gca*_ab)
       call abab(tr,ti,alar,alai,gtcar,gtcai)
       inv(11)=1.5d0*ti
      endif

      if(lusea)then
!      betasA2: (w/2)Im(i als_ab e_agd As*g,db)
       call abab(tr,ti,alsr,alsi,easr,easi)
       inv(7)=0.5d0*w*tr/clight
!      betaaA2: (w/2)Im(i[ala_ab[e_agd Aa*g,db+ e_abg Aa*d,gd)
       call abab(tr,ti,alar,alai,eaar,eaai)
       inv(8)=tr
       call abab(tr,ti,alar,alai,eapar,eapai)
       inv(8)=0.5d0*w*(inv(8)+tr)/clight
!      betasAc2: (w/2)Im(i als_ab e_agd Acs*g,db)
       call abab(tr,ti,alsr,alsi,eacsr,eacsi)
       inv(12)=0.5d0*wR*tr/clight
!      betaaAc2:
       call abab(tr,ti,alar,alai,eacar,eacai)
       inv(13)=tr
       call abab(tr,ti,alar,alai,eapcar,eapcai)
       inv(13)=0.5d0*wR*(inv(13)+tr)/clight
      endif
      
      
      YDX=0D0
      YDY=0d0
!     spectral intensities:
      do is=1,ns0
      if(ldo(is))then
!      Raman:
       a(1)= &
       Kcon*co(1,is)*(co(3,is)*inv(1)+co(4,is)*inv(2)+co(5,is)*inv(3))
!      ROA:
       a(2)=0.0d0
       do ii=4,ni0
         a(2)=a(2)+inv(ii)*co(2+ii,is)
       end do
       a(2)=co(2,is)*a(2)*kcon/cc_au
       a_arr(1,is)=a(1)
       a_arr(2,is)=a(2)
       WRITE(tabOutputs(is),3001)nt,ECM,YDX,YDY,a(1),a(2),bf,TR(str_v1),TR(str_v3)
3001   FORMAT(I7,f9.2,2f3.0,g12.4,' 0 0 0 0',2g12.4,A,' --> ',A)
       flush(tabOutputs(is))
      end if
      end do
      
      return
   end subroutine wrram3_noSpectra
   
   function MakeLameInvariants(Ap,G,Gc,A,Ac,w)result(inv)
      double complex Ap(3,3),G(3,3),Gc(3,3),A(3,3,3),Ac(3,3,3)
      double complex inv(15)
      double precision w

      
      !From Hecht 1991
      !a0,a1,a2
      inv(1)=TensorProduct_AABB(Ap,Ap)
      inv(2)=TensorProduct_ABAB(Ap,Ap)
      inv(3)=TensorProduct_BAAB(Ap,Ap)
      
      !G0,G1,G2
      inv(4)=TensorProduct_AABB(Ap,G)
      inv(5)=TensorProduct_ABAB(Ap,G)
      inv(6)=TensorProduct_BAAB(Ap,G)
      
      !Gc0,Gc1,Gc2
      inv(7)=TensorProduct_AABB(Ap,Gc)
      inv(8)=TensorProduct_ABAB(Ap,Gc)
      inv(9)=TensorProduct_BAAB(Ap,Gc)
      
      !A1,A2,A3
      inv(10)=TensorProduct_ABCDB(Ap,A,w)
      inv(11)=TensorProduct_BACDB(Ap,A,w)
      inv(12)=TensorProduct_ABDCD(Ap,A,w)
      
      !Ac1,Ac2,Ac3
      inv(13)=TensorProduct_ABCDB(Ap,Ac,w)
      inv(14)=TensorProduct_BACDB(Ap,Ac,w)
      inv(15)=TensorProduct_ABDCD(Ap,Ac,w)
   end function MakeLameInvariants
   
   function LameInvariants2ComplicatedInvariants(inv)result(inv2)
      double complex inv(15)
      double precision inv2(13)
      
      !TODO: Check the first three cool invariants cause the first 3 lame invariants ought to be purely real
      inv2(1)=realpart(inv(1))/9d0 !a2
      inv2(2)=(-2*realpart(inv(1))+3*realpart(inv(2))+3*realpart(inv(3)))/4d0 !Bs(a)2
      inv2(3)=3d0/4d0*(realpart(inv(2))-realpart(inv(3))) !Ba(a)2
      
      inv2(4)=imagpart(inv(4))/9d0 !aG
      inv2(5)=imagpart(-2*inv(4)+3*inv(5)+3*inv(6))/4d0 !Bs(G)2
      inv2(6)=imagpart(inv(5)-inv(6))*3d0/4d0 !Ba(G)2
      
      inv2(7)=imagpart(inv(7))/9d0
      inv2(8)=imagpart(-2*inv(7)+3*inv(8)+3*inv(9))/4d0
      inv2(9)=imagpart(inv(8)-inv(9))*3d0/4d0
      
      inv2(10)=imagpart(inv(10)+inv(11))*3d0/4d0
      inv2(11)=imagpart(inv(10)-inv(11)+2*inv(12))*3d0/4d0
      
      inv2(12)=imagpart(inv(13)+inv(14))*3d0/4d0
      inv2(13)=imagpart(inv(13)-inv(14)+2*inv(15))*3d0/4d0
   end function LameInvariants2ComplicatedInvariants
   
   function TensorProduct_AABB(A,B)result(res)
      double complex A(3,3),B(3,3),res
      res=(A(1,1)+A(2,2)+A(3,3))*conjg(B(1,1)+B(2,2)+B(3,3))
   end function TensorProduct_AABB
   
   function TensorProduct_ABAB(A,B)result(res)
      integer i,j
      double complex A(3,3),B(3,3),res
      res=0d0
      do i = 1,3
         do j = 1,3
            res=res+A(i,j)*conjg(B(i,j))
         end do
      end do
   end function TensorProduct_ABAB
   
   function TensorProduct_BAAB(A,B)result(res)
      integer i,j
      double complex A(3,3),B(3,3),res
      res=0d0
      do i = 1,3
         do j = 1,3
            res=res+A(j,i)*conjg(B(i,j))
         end do
      end do
   end function TensorProduct_BAAB
   
   function TensorProduct_ABCDB(At,Bt,w)result(res)
      integer a,b,c,d
      double complex At(3,3),Bt(3,3,3),res
      double precision w
      
      res=0d0
      do a=1,3
         do b=1,3
            do c=1,3
               do d=1,3
                  res=res+eps(a,c,d)*At(a,b)*conjg(Bt(c,d,b))
               end do
            end do
         end do
      end do
      res=1d0/3d0*iu*w*res
   end function TensorProduct_ABCDB
   
   function TensorProduct_BACDB(At,Bt,w)result(res)
      integer a,b,c,d
      double complex At(3,3),Bt(3,3,3),res
      double precision w
      
      res=0d0
      do a=1,3
         do b=1,3
            do c=1,3
               do d=1,3
                  res=res+eps(a,c,d)*At(b,a)*conjg(Bt(c,d,b))
               end do
            end do
         end do
      end do
      res=1d0/3d0*iu*w*res
   end function TensorProduct_BACDB
   
   function TensorProduct_ABDCD(At,Bt,w)result(res)
      integer a,b,c,d
      double complex At(3,3),Bt(3,3,3),res
      double precision w
      
      res=0d0
      do a=1,3
         do b=1,3
            do c=1,3
               do d=1,3
                  res=res+eps(a,b,c)*At(a,b)*conjg(Bt(d,c,d))
               end do
            end do
         end do
      end do
      res=1d0/3d0*iu*w*res
   end function TensorProduct_ABDCD
   
   subroutine dsa(tr,ti,tsr,tsi,tar,tai)
   implicit none
   real*8 tr(3,3),ti(3,3),tsr(3,3),tsi(3,3),tar(3,3),tai(3,3)
   integer*4 a,b
   do a=1,3
      do b=1,3
         tsr(a,b)=0.5d0*(tr(a,b)+tr(b,a))
         tsi(a,b)=0.5d0*(ti(a,b)+ti(b,a))
         tar(a,b)=0.5d0*(tr(a,b)-tr(b,a))
         tai(a,b)=0.5d0*(ti(a,b)-ti(b,a))
      end do
   end do
   return
   end subroutine dsa
   
   subroutine calcep(e,a)
   !  eap_ab = eps_abc A_d,cd
      implicit none
      real*8 e(3,3),a(3,3,3)
      e(1,1)=0.0d0
      e(2,2)=0.0d0
      e(3,3)=0.0d0
      e(1,2)= a(1,3,1)+a(2,3,2)+a(3,3,3)
      e(1,3)=-a(1,2,1)-a(2,2,2)-a(3,2,3)
      e(2,3)= a(1,1,1)+a(2,1,2)+a(3,1,3)
      e(2,1)=-e(1,2)
      e(3,1)=-e(1,3)
      e(3,2)=-e(2,3)
      return
   end subroutine calcep
   
   subroutine abab(tr,ti,ar,ai,br,bi)
   !  t = a_ab b*_ab
      implicit none
      integer*4 a
      real*8 tr,ti,ar(3,3),ai(3,3),br(3,3),bi(3,3)
      tr=0.0d0
      ti=0.0d0
      do a=1,3
      tr=tr+ar(a,1)*br(a,1)+ai(a,1)*bi(a,1) &
         +ar(a,2)*br(a,2)+ai(a,2)*bi(a,2) &
         +ar(a,3)*br(a,3)+ai(a,3)*bi(a,3) 
      ti=ti-ar(a,1)*bi(a,1)+ai(a,1)*br(a,1) &
         -ar(a,2)*bi(a,2)+ai(a,2)*br(a,2) &
         -ar(a,3)*bi(a,3)+ai(a,3)*br(a,3)
      END DO
      return
   end subroutine abab
      
   function DoPrescreen(v1,v1_pos,v1_n,td,wexc,nexc,w_from,w_ps,ignore_ps,check,uncoup_modes,ignore_imag_modes,n_thr)result(exc_m_v)
      type(ExcState) td
      integer v1_n,v1(v1_n),nexc,n_thr
      integer(int16) v1_pos(v1_n)
      double precision w_from,wexc(nexc)
      integer exc_m_v(td%nq,td%mc_ms_ps,nexc)
      logical w_ps,ignore_ps,check,uncoup_modes(td%nq),ignore_imag_modes
      
      double precision curThresh,wn,thr_arr(nexc)
      integer ii
      
      curThresh=DEM_1v_pickThr(v1,v1_n,td%thr_v2,td%thr_v2_0,td%thr_v2_over,td%thr_v2_comb,td%thr_v2_other)
      wn=w_from
      
      if(w_ps)then
         do ii = 1,nexc
            !thr_arr(ii)=curThresh*sqrt(()/(td%e_00-wn-wexc(ii))**2+td%gamma**2)
            thr_arr(ii)=curThresh
         end do
      else
         thr_arr(1)=curThresh
      end if
      if(ignore_ps)then
         do ii = 1,td%mc_ms_ps
            exc_m_v(:,ii,:)=td%max_v2s(ii)
         end do
      else
         exc_m_v=DEM_1v_rroa(v1,v1_pos,v1_n,td%max_v2s_ps,td%mc_ms_ps,td,thr_arr,check,uncoup_modes,w_ps,wn,wexc,nexc,n_thr,td%gamma,td%fc_sum_min_ps,td%fc_sum_min_ps_its,ignore_imag_modes)
      end if
   end function DoPrescreen
   
   
   subroutine ReadCustomFreqs(modes,n_m)
      type(Mode),allocatable :: modes(:)
      integer n_m
      integer :: v(3),v_pos(3),idx,idx2,v_n,i
      double precision :: w
      character(200) s200
      logical fex
      
      inquire(file='FILE.W.TR',exist=fex)
      if(.not.fex)then
         write(output_unit,*)'Reading custom frequencies was requested but FILE.W.TR was not found!'
         call exit(667)
      end if
      
      open(667,file='FILE.W.TR')
      read(667,*)n_m
      read(667,*)
      do i = 1,n_m
         v_pos=0
         v=0
         read(667,*)s200
         idx=index(s200,':')
         idx2=index(s200,',')
         read(s200(1:idx-1),*)v_pos
         read(s200(idx+1:idx2-1),*)v
         read(s200(idx2+1:200),*)w
         v_n=count(v/=0)
         modes(i)%v=v(1:v_n)
         modes(i)%v_pos=v_pos(1:v_n)
         modes(i)%w=w
      end do
      
      close(667)
   end subroutine ReadCustomFreqs
   
   
   subroutine OrderModesToW(fc_ss,v_modes)
      type(FC_storSys) fc_ss
      type(Mode) v_modes(:)
      
      integer i,v1_class
      integer,allocatable :: v1(:)
      integer(int16),allocatable :: v1_pos(:),v_pos_demapped(:)
      double precision :: wm
      double precision,allocatable :: ws(:)
      integer,allocatable :: order(:)
      type(Mode),allocatable :: v_modes_copy(:)
      
      allocate(order(fc_ss%fc_arr_size),ws(fc_ss%fc_arr_size))
      
      do i = 1,fc_ss%fc_arr_size
         v1=v_modes(i)%v
         v1_pos=v_modes(i)%v_pos
         v1_class=count(v1>0)
         v_pos_demapped=DemapVPos(v1_pos,v1_class,fc_ss%v_map,fc_ss%red_n)
         wm=getFreq_short(v1,v_pos_demapped,fc_ss%w,v1_class)*au_2_cm
         ws(i)=wm
         v_modes(i)%w=wm
      end do
      
      call linsort_D_big(ws,fc_ss%fc_arr_size,order)
      v_modes_copy=v_modes
      do i = 1,fc_ss%fc_arr_size
         v_modes(i)=v_modes_copy(order(i))
      end do
      
      deallocate(order,ws,v1,v1_pos,v_modes_copy,v_pos_demapped)
      ! str_v1=FC2Str_new(v1,v1_pos_mapped,v1_class,.true.,.false.)
   end subroutine OrderModesToW
   
   subroutine IgnoreImagModesInExcm(exc_m,mc,nq,imag_modes)
      integer mc,nq,exc_m(nq,mc)
      logical imag_modes(nq)
      
      integer i
      do i = 1,nq
         if(imag_modes(i))THEN
            exc_m(i,:)=0
         end if
      end do
   end subroutine IgnoreImagModesInExcm

   subroutine DoRROA(n_thr,runtype,tds,td_n,wgg,wgg_nq,wexc,nexc,kbt_lim,temp,wmin,wmax,npx,FWHM,rroa_spr_do,st,ht,ht2,only_ht,useA,useG, &
                     check,use_gauss,write_excm,write_ten,write_inv,write_fc_sys,write_fcarr,norm_sp,UnCoup_lim,uncoup_modes_maxval,w_ps,ignore_ps,sel_rules, &
                     sel_rules2,correct_gr_freqs,add_pol_coeff,output_moments,output_polars,ignore_imag_modes,output_polContrs,alphafc,wr_is_e00, &
                     spectrum_temp,doModes,trs,fundLeadCount,TD_approach,td_N_points,td_tmax,td_fs,td_sparse,td_alt,td_fixphase,J_tol,td_batch_n,trueground,write_corrf,write_fft_cf,num_integ, &
                     norm_fft,correctPhaseX,correctPhaseX_abs,interpolateFFT,wexc_adapt,w_ad_zero,evcd,cpl,contr)
      integer td_n,nexc
      type(ExcState),target :: tds(td_n)
      type(Transition_R),allocatable :: trs(:)
      type(Transition_R),allocatable :: trs_fund(:,:)
      type(Transition_R) :: cur_tr
      integer :: trs_c,trs_i,trs_fund_c,fundleadcount
      logical correct_gr_freqs,output_polContrs,fundLead,wexc_adapt
      logical write_corrf,write_fft_cf,num_integ,correctPhaseX,correctPhaseX_abs,interpolateFFT
      
      logical :: TD_approach,contr(9),td_alt,td_fixphase
      double precision,allocatable :: gamma_gr(:,:),v_td(:)
      integer, allocatable :: br_order(:),batch(:)
      double precision kgk,td_tmax,td_fs
      integer td_N_points,nz
      integer td_batch_n,batch_n,batches_c,bc
      type(Polar_exc),allocatable,target :: polars(:)
      
      logical td_sparse,trueground,norm_fft,w_ad_zero
      double precision J_tol
      type(list_arr_int32),allocatable :: J_nz(:)
      
      double precision temp,FWHM,UnCoup_lim,add_pol_coeff,wr,td_wmax,w_nat
   
      integer c_gs,c_ms,c_gfs,i,ii,iii,iv,j,idx,jj,k,kk,nq,jjj
      integer nq_r,v_start
      integer v,v_pos
      integer(int64) comb_n,i_64,v1_idx,v3_idx
      integer(int64),allocatable :: combs_n(:)
      
      type(modes_arr_t),allocatable :: modes_arr_v2(:),modes_arr_v1(:),modes_arr_v2_c1(:)
      type(Mode_col) :: vs_td(td_n)
      type(Polar_exc) :: curpolar,curPolar_anti
      type(ExcState),pointer :: td
      type(fc_storSys) fc_ss1,fc_ss2,fc_ss_gr,fc_ss2_c1
      type(fc_storSys_dusch) fc_ss_dusch,fc_ss_dusch_c1
      type(fc_storSys_derivs) fc_ss_derivs,fc_ss_derivs_c1
      type(Mode),allocatable :: v1_modes(:)
      type(Mode),allocatable :: corrected_modes(:)
      type(Mode) :: vcur
      type(Mode),allocatable :: m_ps_done(:,:)
      type(Mode),allocatable :: m_ps_done_gr(:)
      integer :: m_ps_done_c,m_ps_done_c_gr,m_ps_done_gr_c
      type(v_col),allocatable :: v2_arr(:)
      integer,parameter :: ns0 = 19,ni0=13
      
      integer :: gp_v1_c = 5,n_m
      integer, parameter :: gp_v3_c = 5000 !this will do for now
      type(Mode),allocatable :: ground_polars_v1(:)
      type(Mode),allocatable :: ground_polars_v3(:)
      type(Polar_exc),allocatable :: ground_polars(:,:)
      integer :: pol_idx_v1,pol_idx_v3
      integer :: pol_idx_v1_last,pol_idx_v3_last
      integer(int16),allocatable :: v1_pos_mapped(:),v3_pos_mapped(:)
      logical :: mode_found,w_ps,ignore_ps,sel_rules,sel_rules2,output_moments,output_polars
      logical :: ignore_imag_modes,alphafc,spectrum_temp,wr_is_e00
      logical :: doIn,doFin
      logical,allocatable :: mode_psed(:,:)
      integer ,allocatable :: mode_psed_from(:),doModes(:)
      
      integer :: wgg_nq,n_thr,maxVal1,maxVal2
      double precision :: wgg(wgg_nq),wm
      
      character(3) :: runtype
      character(2) :: TMExpand
      character(20) str_v1,str_v3
      character(20) :: wexc_nm_str(nexc)
      character(9),parameter :: rroa_exp(ns0) = &
        ['ICP_0    ','ICP_x_90 ','ICP_z_90 ','ICP_*_90 ','ICP_u_90 ','ICP_180  ',&
         'SCP_0    ','SCP_x_90 ','SCP_z_90 ','SCP_*_90 ','SCP_u_90 ','SCP_180  ',&
         'DCPI_0   ','DCPI_90  ','DCPI_180 ',&
         'DCPII_0  ','DCPII_90 ','DCPII_180',&
         'SPECIAL  ']
      character(80) c80,c80_2
      
      logical :: isOMP,check,st,ht,ht2,full_ht,gse,usea,useg,only_ht
      logical :: use_gauss,write_excm,write_ten,write_inv,write_fc_sys,norm_sp,write_fcarr
      integer :: td_i,nthr,n_gr,i_gr,i_gr2,tr_idx,npx,iexc,tr_count
      integer :: mode_i,mode_f,v13_n,v13_mc,fundLeadCount_current
      
      integer,allocatable :: exc_m2(:,:,:,:,:),exc_m1(:,:),exc_m13(:,:),exc_m_temp(:,:,:)
      integer,allocatable :: exc_m13_big(:,:),exc_m_c1(:,:)
      integer,allocatable :: exc_m_union(:,:) !or maybe overlap
      integer,allocatable :: v1(:),v3(:),v2(:)
      integer(int16),allocatable :: v1_pos(:),v2_pos(:),v3_pos(:),v13_pos(:)
      integer(int16),allocatable :: v1_pos_short(:),v3_pos_short(:)
      integer,allocatable :: v13(:)
      integer :: v1_class,v3_class,v13_class,wmin,wmax,v2_class,v1_n,grCount,i_iexc
      integer(int64) :: HighOVIdx(2),counts(16)
      double precision :: w1,w3,kBT,kbt_anti,kbt_lim,wexc(nexc),w,fc_sum1,fc_sum3,time1,time2,w_ve,wn
      double precision,allocatable :: fc_arr(:,:),fc_arr_c1(:,:)
      
      double precision,allocatable :: rroa_spr(:,:,:,:),rroa_spr_temp(:,:,:),rroa_inv(:)
      double precision,allocatable :: x(:),thr_arr(:),jgj(:,:)
      type(Polar_exc), allocatable :: polar_ders(:)
      double precision intgRAM,intgROA,curThresh,fc_sum1_tot,fc_sum3_tot
      double complex rroa_inv2(15)
      integer,allocatable :: rroa_tabOutputs(:,:)
      logical :: rroa_spr_do(ns0),evcd,cpl
      logical,allocatable :: canBeLow(:),canBeLow_gr(:),uncoup_modes(:),doGround(:),v_success(:)
      
      integer,allocatable :: uncoup_modes_excm(:),uncoup_modes_vmap(:),uncoup_modes_reddims(:)
      integer :: uncoup_modes_maxval,uncoup_modes_c
      double precision,allocatable :: uncoup_fcarr(:)
      double precision :: fcsum1_un,fcsum3_un
      
      logical :: success,success_gr
      double precision :: a_arr(2,ns0)
      
      allocate(rroa_spr(npx,2,ns0,nexc),rroa_inv(ni0),rroa_tabOutputs(ns0,nexc),x(npx))
      do iexc = 1,nexc
         do i = 1,ns0
            rroa_tabOutputs(i,iexc)=3000+i+(iexc-1)*ns0
         end do
      end do
      rroa_spr=0d0
      rroa_inv=0d0
      rroa_inv2=0d0
      !rroa_spr_do=.false.
      !rroa_spr_do(12)=.true.
      fundlead=fundLeadCount>0
      sel_rules=sel_rules .or. fundlead
      
      
      if(TD_approach)then
         br_order=BRMake(td_N_points)
         td_wmax=td_fs/2d0
      end if
      
      do i = 1,npx
         x(i)=dble(wmax-wmin)/(npx-1)*(i-1)+wmin
      end do
      
      do i = 1,nexc
         if(wexc(i)==0d0)then
            write(wexc_nm_str(i),'(F6.1,"nm")')0.0d0
         else
            write(wexc_nm_str(i),'(F6.1,"nm")')dble(NINT(10d0**7/(wexc(i))*cm_2_au*10))/10d0
         end if
      end do
      if(.false.)then
         do iexc = 1,nexc
            do i = 1,ns0
               if(rroa_spr_do(i))then
                  open(rroa_tabOutputs(i,iexc),file=TR(rroa_exp(i))//'_'//TR(wexc_nm_str(iexc))//'.TAB')
               end if
            end do
         end do
      end if
      allocate(rroa_spr_temp(npx,2,ns0))
      
      isOMP=.false.
      nthr=1
      !$ isOMP=.true.
      !$ call Polar_set_nexc(nexc)
      tr_idx=1
      
      gp_v1_c=maxval(tds(:)%n_gr)
      allocate(ground_polars(gp_v1_c,gp_v3_c))
      allocate(ground_polars_v1(gp_v1_c))
      allocate(ground_polars_v3(gp_v3_c))
      pol_idx_v1=0
      pol_idx_v3=0
      pol_idx_v1_last=0
      pol_idx_v3_last=0
      call AllocatePolar(curPolar,nexc)
      call AllocatePolar(curPolar_anti,nexc)
      maxval1=maxval(tds(:)%mc_ms_ps)
      maxval2=maxval(tds(:)%mc_ms)
      call Prepare_nck_arr(maxval(tds(:)%nq),max(maxval1,maxval2))
      
      ! if(add_pol_coeff/=0d0)then
         ! allocate(polar_ders(wgg_nq))
         ! call ReadPolarDers(polar_ders,wgg_nq)
      ! end if
      
      if(correct_gr_freqs)then
         call ReadCustomFreqs(corrected_modes,n_m)
      end if
      if(write_inv)then
         open(99,file='INV.TXT')
         open(66,file='INV2.TXT')
         write(99,'(A13,1X,A13,1X,A8,13(1X,A13))')'init','final','w','a^2','bs(a)^2','ba(a)^2','aG','bs(G)^2','ba(G)^2','bs(A)^2','ba(A)^2','aGc','bs(Gc)^2','ba(Gc)^2','bs(Ac)^2','ba(Ac)^2'
         write(66,'(A13,1X,A13,1X,A8,15(2X,A26))')'init','final','w','a_0','a_1','a_2','G_0','G_1','G_2','Gc_0','Gc_1','Gc_2','A_1','A_2','A_3','Ac_1','Ac_2','Ac_3'
      end if
      
      ! if(allocated(trs))then
         ! trs_c=size(trs,dim=1)
         ! trs_mc_gs=0
         ! trs_mc_fs=0
         ! do i = 1,trs_c
            ! trs_mc_gs=max(trs_mc_gs,size(trs(i)%vi,dim=1))
            ! trs_mc_fs=max(trs_mc_fs,size(trs(i)%vf,dim=1))
         ! end do
         ! trs_mc_gfs=max(trs_mc_gs,trs_mc_fs)
         ! tds(:)%mc_gfs=trs_mc_gfs
      ! end if
      
      if(td_approach)tr_count=tds(1)%nq
      if(write_ten)then
         do iexc=1,nexc
            open(1000+iexc,file='FILE.'//TR(wexc_nm_str(iexc))//'.POLARS')
            write(1000+iexc,*)tds(1)%nq,tds(1)%e_00*au_2_cm
         end do
      end if
      
      if(output_polContrs .and. .not.td_approach)then
         do iexc=1,nexc
            open(op_polcontrs_unit-1+iexc,file='FILE.'//TR(wexc_nm_str(iexc))//'.POLCON')
            write(op_polcontrs_unit-1+iexc,'(A11)')'Ram SCP 180'
         end do
      end if
      
      tr_count=0
      do td_i=1,td_n
         td=>tds(td_i)
         nq=td%nq
         
         w_nat=td%e_00
         do ii=1,nq
            w_nat=w_nat-0.5d0*td%K(ii)**2*(td%wg(ii))**2
         end do
         write(output_unit,'(A,I4)')' TD_ROOT=',td%root
         write(output_unit,'(A,I4)')'   TD_NQ=',td%nq
         write(output_unit,'(A,G10.3)')  '  <0|0> =',td%fc_00
         write(output_unit,'(A,F10.2,A)')'  TD_E00=',td%e_00*au_2_cm,' cm-1'
         write(output_unit,'(A,F10.2,A)')'  TD_E00=',1d7/(td%e_00*au_2_cm),' nm'
         write(output_unit,'(A,F10.2,A)')'   GAMMA=',td%gamma*au_2_cm,' cm-1'
         flush(output_unit)
         
         if(TD_approach)then
            write(output_unit,'(A,F10.2,A)')'   W_NAT=',w_nat*au_2_cm,' cm-1'
            write(output_unit,'(A,F10.2,A)')'TD_T_MAX=',td_tmax,' au'
            !write(output_unit,'(A,G10.3,A)')'   TD_dW=',2*pi/td_tmax*au_2_cm,' cm-1'
            write(output_unit,'(A,I12,A)')'X_POINTS=',(2_8**TD_N_POINTS)
            flush(output_unit)
            nz=0
            allocate(gamma_gr(nq,nq))
            gamma_gr=0d0
            allocate(J_nz(nq))
            
            do i = 1,nq
               nz=0
               gamma_gr(i,i)=td%wg(i)
               do ii = 1,nq
                  if(abs(td%J(ii,i))>=J_tol)nz=nz+1
               end do
               J_nz(i)%nz=nz
               allocate(J_nz(i)%arr(nz))
               iii=1
               do ii = 1,nq
                  if(abs(td%J(ii,i))>=J_tol)then
                     J_nz(i)%arr(iii)=ii
                     iii=iii+1
                  end if
               end do
            end do
            iii=sum(J_nz(:)%nz)
            write(output_unit,'(A,F6.2,A)')'Duschinsky matrix - significant elements: ',dble(iii)/dble(nq**2)*100d0,'%'
!            allocate(J_nonzero(2,nz))
            nz=1
            nz=nz-1
            kgk=0d0
            do i = 1,nq
               kgk=kgk+td%K(i)*td%wg(i)*td%K(i)
            end do
            jgj=matmul(transpose(td%J),matmul(gamma_gr,td%J))
            v_td=matmul(transpose(td%J),matmul(gamma_gr,td%K))
            td%J_i=transpose(td%J)
            td%K_i=-matmul(td%J_i,td%K)
            
            if(.not.td_alt)write(output_unit,*)'TD, batches'
            nq_r=0
            if(.not.td_alt)then
               do i = 1,nq
                  wr=td%wg(i)
                  if(wr*au_2_cm>td%max_freq_mode)exit
                  nq_r=nq_r+1
               end do
            end if
            v_start=0
            if(.not.td_alt)then
               do i = 1,nq
                  wr=td%wg(i)
                  if(wr*au_2_cm>=td%min_freq_mode)exit
                  v_start=v_start+1
               end do
            end if
            if(td_alt)then
               batches_c=1
            else
               batches_c=CEILING(dble(nq_r-v_start)/dble(td_batch_n)) !only fundamentals right now so this is sufficent
            end if
            bc=0
            batch_n=0
            v1=[0]
            v1_pos=[0]
            v1_class=0
            v3=[1]
            v3_class=1
            w1=0d0
            kbt=1
            str_v1='|0 >'
            if(num_integ)then
               write(output_unit,*)'Trapezoidal integration'
            else
               write(output_unit,*)'FFT integration'
               if(interpolateFFT)write(output_unit,*)'FFT interpolation'
            end if
            
            do i = 1,batches_c
               if(td_alt)then
                  batch_n=nq
               elseif(i==batches_c .and. mod(nq_r,td_batch_n)/=0)then!last batch, ie. batch_n <= td_batch_n
                  batch_n=mod(nq_r,td_batch_n)
               else
                  batch_n=td_batch_n
               end if
               allocate(batch(batch_n),polars(batch_n))
               
               batch=0
               do ii = 1,batch_n
                  iii=(i-1)*td_batch_n+ii
                  batch(ii)=iii+v_start
                  write(output_unit,'(1X,I3)',advance='no')batch(ii)
                  !call Polar_new_nexc(polars(ii),nexc)
               end do
               write(output_unit,*)
               flush(output_unit)
               polars=DoRROA_TD_k(0,nq,batch,batch_n,td%wg,td%we,td%u_gr,td%m_gr,td%q_gr,td%du_gr,td%dm_gr,td%dq_gr,td%u_ex,td%m_ex,td%q_ex,td%du_ex,td%dm_ex,td%dq_ex,td%du2_ex,td%dm2_ex,td%dq2_ex, &
               v_td,gamma_gr,kgk,jgj,td%J,td%K,td%J_i,td%K_i,td%gamma,td%theta,td%eps,td%e_00,td_N_points,td_tmax,td_fs,td_wmax,td_sparse,td_alt,J_nz,nz, &
               ht,ht2,td_fixphase,[td%TMExpand(1:1)=='G',td%TMExpand(2:2)=='G'],st,contr,n_thr,wexc,nexc,br_order,8,write_corrf,write_fft_cf,num_integ,norm_fft,correctPhaseX,correctPhaseX_abs,interpolateFFT,wexc_adapt,w_ad_zero)
               ! do ii = 1,batch_n
                  ! I dunno, TD spectra for Br2 molecule had higher intensity but same excitation wavelength distribution.
                  ! This magic value is the division of the two intensities.
                  ! The intensities still differ quite a lot off-resonance, but in resonance they are pretty well fit.
                  ! polars(ii)=polars(ii)*(1d0/2.25581969334231)
               ! end do
               
               a_arr=0d0
               do ii = 1,batch_n
                  v3_pos=[batch(ii)]
                  w3=td%wg(batch(ii))
                  str_v3=FC2Str_new(v3,v3_pos,v3_class,.true.,.false.)
                  curPolar=polars(ii)
                  curPolar%ap=conjg(curPolar%ap)
                  curPolar%G=conjg(curPolar%G)
                  curPolar%Gc=conjg(curPolar%Gc)
                  curPolar%A=conjg(curPolar%A)
                  curPolar%Ac=conjg(curPolar%Ac)
                  do iexc=1,nexc
                     if(write_ten)then
                        call WriteTenPretty(1000+iexc,wexc(iexc),PolarExc2Polar(curPolar,iexc),v1,v1_pos,v1_class,w1*au_2_cm,v3,v3_pos,v3_class,w3*au_2_cm)
                     end if
                     if(write_inv)then
                        rroa_inv2=MakeLameInvariants(curPolar%Ap(:,:,iexc),curPolar%G(:,:,iexc),curPolar%Gc(:,:,iexc),curPolar%A(:,:,:,iexc),curPolar%Ac(:,:,:,iexc),wexc(iexc))
                        write(99,'(A13,1X,A13,1X,F8.2,13(1X,G13.6),F7.2," nm")')TR(str_v1),TR(str_v3),(w3-w1)*au_2_cm,rroa_inv,10d0**7/wexc(iexc)*cm_2_au
                        write(66,'(A13,1X,A13,1X,F8.2,30(1X,G13.6),F7.2," nm")')TR(str_v1),TR(str_v3),(w3-w1)*au_2_cm,rroa_inv2,10d0**7/wexc(iexc)*cm_2_au
                     end if
                  end do
               end do
               call DeallocatePolar(curPolar)
               bc=bc+1
               do ii = 1,batch_n
                  call DeallocatePolar(polars(ii))
               end do
               deallocate(batch,polars)
            end do
            cycle
         end if
         
         if(fundLead)then
            trs_fund_c=0
            allocate(trs_fund(nq,nexc))
         end if
         
         
         allocate(canBeLow_gr(td%nq),canBeLow(td%nq))
         
         if(uncoup_lim<1d0 .and. .not. TD_approach)then
            uncoup_modes=DUM_EX(td%J,td%nq,.true.,UnCoup_lim)
            uncoup_modes_c=count(uncoup_modes)
            allocate(uncoup_modes_excm(uncoup_modes_c),uncoup_modes_vmap(uncoup_modes_c),uncoup_modes_reddims(uncoup_modes_c))
            uncoup_modes_reddims=0
            uncoup_modes_vmap=0
            ii=1
            do i = 1,nq
               if(uncoup_modes(i))then
                  uncoup_modes_vmap(ii)=i
                  ii=ii+1
               end if
            end do
         else
            allocate(uncoup_modes(td%nq))
            uncoup_modes_c=0
            uncoup_modes=.false.
         end if
         
         
         allocate(exc_m1(td%nq,td%mc_gfs))
         exc_m1=GroundExcM(td)
         if(ignore_imag_modes)call IgnoreImagModesInExcm(exc_m1,td%mc_gfs,nq,td%im_modes_gr)
         
         fc_ss_gr=FC_storSys_make(exc_m1,td%mc_gfs,td%nq,.true.,td%mc_gfs,td%wg)
         allocate(v1_modes(fc_ss_gr%fc_arr_size))
         allocate(modes_arr_v1(fc_ss_gr%mc_v))
         
         v1_modes(1)=MakeMode([0],[0_2],1) !Im just too lazy, so lets collect all included vibrational ground states in an array
         do c_gfs = 1,fc_ss_gr%mc_v
            allocate(v1(c_gfs))
            comb_n=nck(fc_ss_gr%red_n,c_gfs)
            i=1
            ii=1
            idx=1
            modes_arr_v1(c_gfs)%arr=LOMC_Wrap(fc_ss_gr%red_n,c_gfs,comb_n)
            do i_64 = 1,comb_n
               ii=1
               call Excm_v1_modes(ii,v1,modes_arr_v1(c_gfs)%arr(:,i_64),c_gfs,v1_modes,fc_ss_gr,td)
            end do
            deallocate(v1)
         end do
         if(.not.fundlead)call OrderModesToW(fc_ss_gr,v1_modes)
         write(output_unit,'(A)')'No. of possible vibrational states in ground:'
         write(output_unit,'(I8)')fc_ss_gr%fc_arr_size
         allocate(exc_m2(td%nq,td%mc_ms_ps,fc_ss_gr%fc_arr_size,nexc,td%n_gr))
         exc_m2=0
         !determine excited state excitation maxima for all ground state vibrational modes
         !write(output_unit,'(A)')'Prescreening'
         write(output_unit,'(A10,I3)')'MC_MS_PS=',td%mc_ms_ps
         write(output_unit,'(A17,G10.3)')'THR_0/FC_00 =',td%thr_v2_0/td%fc_00
         write(output_unit,'(A17,G10.3)')'THR_FUND/FC_00 =',td%thr_v2/td%fc_00
         write(output_unit,'(A17,G10.3)')'THR_COMB/FC_00 =',td%thr_v2_comb/td%fc_00
         write(output_unit,'(A17,G10.3)')'THR_OVER/FC_00 =',td%thr_v2_over/td%fc_00
         write(output_unit,'(A17,G10.3)')'THR_OTHER/FC_00 =',td%thr_v2_other/td%fc_00
         if(ignore_ps)write(output_unit,*)'Prescreening will be ignored, modes set to max_v2s'
         flush(output_unit)
         time1=getTime()
         
         if(w_ps)then
            allocate(thr_arr(nexc))
         else
            allocate(thr_arr(1))
         end if
                  
         time1=getTime()-time1
         write(output_unit,*)
         write(output_unit,'(3X,A)')TR(GetTimeStr(time1))
         write(output_unit,*)
         !do transitions and polarizabilities
         i_gr=0
         if(write_fc_sys)then
            open(10,file='FCSYS_TD'//TR(I2STR(td%root,3))//'.TXT')
         end if
         
         allocate(mode_psed(fc_ss_gr%fc_arr_size,td%n_gr),mode_psed_from(td%n_gr))
         mode_psed=.false.
         mode_psed_from=HUGE(1)
         trs_c=0
         do i = 1,fc_ss_gr%fc_arr_size !loop over initial states
            i_gr=i_gr+1
            if(i_gr>td%n_gr)exit
            doIn=.true.
            
            v1=v1_modes(i)%v
            v1_pos_mapped=v1_modes(i)%v_pos
            v1_class=count(v1/=0)
            if(i>1)then
               if(maxval(v1)>td%max_v1s(v1_class))cycle
               if(allocated(doModes))then
                  doIn=.false.
                  allocate(v_success(v1_class))
                  v_success=.false.
                  do ii = 1,v1_class
                     do iii=1,size(domodes,dim=1)
                        if(doModes(iii)==0)exit
                        if(doModes(iii)==v1_pos_mapped(ii))then
                           v_success(ii)=.true.
                        end if
                     end do
                  end do
                  if(all(v_success))then
                     doIn=.true.
                  else
                     i_gr=i_gr-1
                  end if
                  deallocate(v_success)
               end if

               v1_pos=DemapVPos(v1_pos_mapped,v1_class,fc_ss_gr%v_map,fc_ss_gr%red_n)
            else
               v1_pos=[0]
            end if
            if(.not.doIn)cycle
            !wm=getFreq_short(v1,v1_pos,fc_ss_gr%w,v1_class)*au_2_cm
            !v1_pos_short=GetSeq(v1_class)
            curThresh=DEM_1v_pickThr(v1,v1_class,td%thr_v2,td%thr_v2_0,td%thr_v2_over,td%thr_v2_comb,td%thr_v2_other)
            
            if(v1_class>td%mc_gs)cycle !perhaps exit
            if(i==1)then
               w1=0d0
               str_v1='|0 >'
            else
               w1=getFreq_short(v1,v1_pos_mapped,td%wg,v1_class)
               str_v1=FC2Str_new(v1,v1_pos_mapped,v1_class,.true.,.false.)
            end if
            ! if(w1*au_2_cm<td%min_freq_mode .or. w1*au_2_cm>td%max_freq_mode)then
               ! cycle
            ! end if
            kBT=exp(-h*cc*w1*au_2_cm*100d0/(kb*temp)) !check this
            if(kbt<kbt_lim)then
               i_gr=i_gr-1
               cycle
            end if
            canBeLow_gr=.false.
            do ii = 1,v1_class
               canBeLow_gr(v1_pos_mapped(ii))=.true.
            end do
            
            fundLeadCount_current=0
            do j = 2,fc_ss_gr%fc_arr_size !loop over final states
               if(i==j)cycle
               
               v3=v1_modes(j)%v
               v3_pos_mapped=v1_modes(j)%v_pos
               v3_class=count(v3/=0)
               doFin=.true.
               if(allocated(doModes))then
                  allocate(v_success(v3_class))
                  v_success=.false.
                  doFin=.false.
                  do ii = 1,v3_class
                     do iii=1,size(domodes,dim=1)
                        if(doModes(iii)==0)exit
                        if(doModes(iii)==v3_pos_mapped(ii))then
                           v_success(ii)=.true.
                        end if
                     end do
                  end do
                  if(all(v_success))THEN
                     doFin=.true.
                  end if
                  deallocate(v_success)
               end if
               if(.not.doFin)cycle
               
               v3_pos=DemapVPos(v3_pos_mapped,v3_class,fc_ss_gr%v_map,fc_ss_gr%red_n)
               !v3_pos_short=GetSeq(v3_class)
               if(maxval(v3)>td%max_v3s(v3_class))cycle
               if(v3_class>td%mc_fs)cycle
               w3=getFreq_short(v3,v3_pos_mapped,td%wg,v3_class)
               if(w3<w1)cycle
               ! if(w3*au_2_cm<td%min_freq_mode .or. w3*au_2_cm>td%max_freq_mode)then
                  ! cycle
               ! end if
               kbt_anti=exp(-h*cc*w3*au_2_cm*100d0/(kb*temp))
               str_v3=FC2Str_new(v3,v3_pos_mapped,v3_class,.true.,.false.)
               
               if((w3-w1)*au_2_cm<td%min_freq_mode .or. (w3-w1)*au_2_cm>td%max_freq_mode)cycle
               if(sel_rules)then
                  if(.not.IsAllowedTransition(v1,v1_pos,v1_class,v3,v3_pos,v3_class))cycle
               end if
               
               write(output_unit,*)'------------------------------'
               write(output_unit,*)
               write(output_unit,*)TR(str_v1),' --> ',TR(str_v3)
               write(output_unit,'(F9.2,A)')(w3-w1)*au_2_cm,' cm-1'
               write(output_unit,'(A,G10.2)')'Boltzmann lower: ',kbt
               write(output_unit,'(A,G10.2)')'Boltzmann upper: ',kbt_anti
               flush(output_unit)
               if(sel_rules .and. fundLead .and. v3_class==2 .and. v1_class==1)then
                  !if(findloc(v3_pos_mapped,v1_pos_mapped(1))==0)
                  if(fundLeadCount_current>=fundLeadCount)then
                     doFin=.false.
                     cycle
                  end if
                  idx=0
                  if(v3_pos_mapped(1)==v1_pos_mapped(1))then
                     idx=FindTransition_R_c1(trs_fund,nq,trs_fund_c,[0],[0_2],0,v3(2:2),v3_pos_mapped(2:2),1,1,nexc)
                  else if(v3_pos_mapped(2)==v1_pos_mapped(1))then
                     idx=FindTransition_R_c1(trs_fund,nq,trs_fund_c,[0],[0_2],0,v3(1:1),v3_pos_mapped(1:1),1,1,nexc)
                  end if
                  if(idx>fundLeadCount .or. idx==0)then
                     doFin=.false. 
                     write(output_unit,*)'Skipping this transition, fundamental is not intense enough.'
                     cycle
                  end if
                  fundLeadCount_current=fundLeadCount_current+1
               end if
               
               if(.not.mode_psed(i,i_gr))then !maybe?
                  exc_m2(:,:,i,:,i_gr)=DoPrescreen(v1,v1_pos_mapped,v1_class,td,wexc,nexc,w1,w_ps,ignore_ps,check,uncoup_modes,ignore_imag_modes,n_thr)
                  mode_psed_from(i_gr)=i
                  mode_psed(i,i_gr)=.true.
               end if
               
               if(.not.mode_psed(j,i_gr))then !maybe?
                  exc_m2(:,:,j,:,i_gr)=DoPrescreen(v3,v3_pos_mapped,v3_class,td,wexc,nexc,w1,w_ps,ignore_ps,check,uncoup_modes,ignore_imag_modes,n_thr)
                  mode_psed(j,i_gr)=.true.
               end if
               
               call v_union(v1,v1_pos_mapped,v1_class,v3,v3_pos_mapped,v3_class,v13,v13_pos,v13_n,exc_m13,v13_mc,td%nq)
               fc_ss1=FC_storSys_make(exc_m13,v13_mc,td%nq,.true.,v13_mc,td%wg)
               fc_ss1%v_map=v13_pos
               if(v1_class==0)then
                  v1_pos_short=[0]
               else
                  allocate(v1_pos_short(v1_class))
                  do ii=1,v1_class
                     v1_pos_short(ii)=findloc(fc_ss1%v_map,v1_pos_mapped(ii),dim=1)
                  end do
               end if
               if(v3_class==0)then
                  v3_pos_short=[0]
               else
                  allocate(v3_pos_short(v3_class))
                  do ii=1,v3_class
                     v3_pos_short(ii)=findloc(fc_ss1%v_map,v3_pos_mapped(ii),dim=1)
                  end do
               end if
               
               
               v1_idx=ReduceFC_idx_fcss(v1,v1_pos_short,v1_class,fc_ss1)
               v3_idx=ReduceFC_idx_fcss(v3,v3_pos_short,v3_class,fc_ss1)
               call Polar_new(curpolar)
               if(td%antiStokes)call Polar_new(curPolar_anti)
               trs_c=trs_c+1
               
               !unionize excited vibrational states which are large for initial or final state
               do iexc=1,nexc
                  if(output_polContrs)then
                     if(.not.w_ps)then
                        do i_iexc = 1,nexc
                           write(op_polcontrs_unit-1+i_iexc,'(1X,A24)',advance='no')'Init.state -> Fin.state'
                           if(v1_class==0)then
                              write(op_polcontrs_unit-1+i_iexc,'(1X,I3,"^",I2)',advance='no')0,0
                           else
                              do ii=1,v1_class
                                 write(op_polcontrs_unit-1+i_iexc,'(1X,I3,"^",I2)',advance='no')v1_pos(ii),v1(ii)
                              end do
                           end if
                           do ii=1,v3_class
                              write(op_polcontrs_unit-1+i_iexc,'(1X,I3,"^",I2,1X,F10.4," cm-1")',advance='no')v3_pos(ii),v3(ii),(w3-w1)*au_2_cm
                           end do
                           write(op_polcontrs_unit-1+i_iexc,*)
                        end do
                     else
                        write(op_polcontrs_unit-1+iexc,'(1X,A24)',advance='no')'Init.state -> Fin.state'
                        if(v1_class==0)then
                           write(op_polcontrs_unit-1+iexc,'(1X,I3,"^",I2)',advance='no')0,0
                        else
                           do ii=1,v1_class
                              write(op_polcontrs_unit-1+iexc,'(1X,I3,"^",I2)',advance='no')v1_pos(ii),v1(ii)
                           end do
                        end if
                        do ii=1,v3_class
                           write(op_polcontrs_unit-1+iexc,'(1X,I3,"^",I2,1X,F10.4," cm-1")',advance='no')v3_pos(ii),v3(ii),(w3-w1)*au_2_cm
                        end do
                        write(op_polcontrs_unit-1+iexc,*)
                     end if
                  end if
                  write(output_unit,*)
                  write(output_unit,*)
                  write(output_unit,'("w_exc == ",A)')TR(wexc_nm_str(iexc))
                  exc_m_union=MakeExcMaxUnion_1ExcSt(exc_m2(:,:,i,iexc,i_gr),td%mc_ms_ps,exc_m2(:,:,j,iexc,i_gr),td%mc_ms_ps,td%nq)
                  
                  do ii = 1,min(td%mc_ms_ps,td%mc_ms)
                     iii=0
                     jjj=0
                     do jj = 1,td%nq
                        exc_m_union(jj,ii)=max(exc_m_union(jj,ii),td%min_v2s_excm(ii))
                        if(uncoup_modes(jj))then
                           iii=iii+1
                           uncoup_modes_reddims(iii)=0
                           if(findloc(v1_pos_mapped,jj,dim=1) > 0 .or. findloc(v3_pos_mapped,jj,dim=1) > 0)cycle
                           if(uncoup_modes_maxval<0)then
                              uncoup_modes_reddims(iii)=maxval(exc_m_union(jj,:),dim=1)
                           else
                              uncoup_modes_reddims(iii)=uncoup_modes_maxval
                           end if
                           exc_m_union(jj,ii)=0
                        end if
                        
                     end do
                  end do
                  if(uncoup_modes_c>0)then
                     allocate(uncoup_fcarr(sum(uncoup_modes_reddims)))
                     iii=0
                     do ii=1,uncoup_modes_c
                        v_pos=uncoup_modes_vmap(ii)
                        do v = 1,uncoup_modes_reddims(ii)
                           iii=iii+1
                           uncoup_fcarr(iii)=FC_c1(v,td%fc_00,td%C(v_pos,v_pos),td%D(v_pos))
                        end do
                     end do
                  end if
                  write(output_unit,*)
                  write(output_unit,'(A)')'#FC system...'
                  flush(output_unit)
                  time1=getTime()
                  fc_ss2=FC_storSys_make_minmaxVs(exc_m_union,td%mc_ms_ps,td%nq,td%min_v2s,td%max_v2s,.true.,full_ht,td%we,td%mc_ms,td%we, &
                   td%max_osc_v2,td%max_osc_v2_algo,td%max_osc_v2_coeff,td%B)
                           
                  if(write_fc_sys)then
                     write(10,'(A,A,A)')TR(str_v1),' --> ',TR(str_v3)
                     write(10,'(F9.2,A)')(w3-w1)*au_2_cm,' cm-1'
                     
                     write(10,'(A)')'## GR. FC-SYS:'
                     call WriteFCSys(10,fc_ss1)
                     write(10,*)
                     write(10,'(A)')'## EX. FC-SYS:'
                     call WriteFCSys(10,fc_ss2)
                     write(10,'(A)')'------------------------------'
                     write(10,*)
                  end if
                  
                  time1=getTime()-time1
                  write(output_unit,'(3X,A)')GetTimeStr(time1)
                  write(output_unit,'(A,I4)')'RED_N2=',fc_ss2%red_n
                  write(output_unit,'(A)')'MAX_VS RED_DIMS:'
                  do iii = 1,fc_ss2%mc_v
                     write(output_unit,'(I4)',advance='no')maxval(fc_ss2%red_dims(:,iii))
                  end do
                  write(output_unit,*)
                  
                  write(output_unit,*)GetSizeStr(FC_storSys_size(fc_ss2))
                  write(output_unit,*)
                  
                  write(output_unit,'(A,I4,A,I12,A)')'(',fc_ss1%fc_arr_size,',',FC_ss2%fc_arr_size,' )'
                  write(output_unit,*)GetSizeStr(fc_ss1%fc_arr_size*fc_ss2%fc_arr_size*kind(fc_arr))
                  if(runtype=='SI1')then
                     call FC_storSys_dispose(fc_ss1)
                     call FC_storSys_dispose(fc_ss2)
                     cycle
                  end if
                  
                  time1=getTime()
                  fc_ss_dusch=FC_storSys_Dusch_make(fc_ss1,fc_ss2,td%A,td%B,td%C,td%D,td%E,td%nq)
                  !fc_ss_derivs=FC_storSys_derivs_make(td,fc_ss2)
                  write(output_unit,*)
                  write(output_unit,'(A)',advance='no')'#FC matrix...'
                  flush(output_unit)
                  fc_arr=MakeFCArr_recalculate(fc_ss1,fc_ss2,fc_ss_dusch,td%fc_00,modes_arr_v2,n_thr,check,.true.,write_fcarr)
                  time1=getTime()-time1
                  write(output_unit,*)
                  write(output_unit,'(3X,A)')GetTimeStr(time1)
                  write(output_unit,*)GetSizeStr(Modes_arr_t_size(modes_arr_v2))
                  
                  !Find highest OV integral
                  HighOVIdx(1) = maxloc(abs(fc_arr(v1_idx,:)),dim=1,kind=8) !slicing the wrong dim, I know
                  HighOVIdx(2) = maxloc(abs(fc_arr(v3_idx,:)),dim=1,kind=8)
                  call UnReduceFC_idx(HighOVIdx(1),v2,v2_pos,v2_class,fc_ss2%red_dims,fc_ss2%mc_v,fc_ss2%red_n,fc_ss2%cs,fc_ss2%fc_ci)
                  write(output_unit,'(A,A,A,G10.2)')' MAX OV INIT: ',TR(FC2Str_new(v2,v2_pos,v2_class,.false.,.false.)),TR(DFC(str_v1,'|')),fc_arr(v1_idx,HighOVIdx(1))
                  deallocate(v2,v2_pos)
                  call UnReduceFC_idx(HighOVIdx(2),v2,v2_pos,v2_class,fc_ss2%red_dims,fc_ss2%mc_v,fc_ss2%red_n,fc_ss2%cs,fc_ss2%fc_ci)
                  write(output_unit,'(A,A,A,G10.2)')'MAX OV FINAL: ',TR(FC2Str_new(v2,v2_pos,v2_class,.false.,.false.)),TR(DFC(str_v3,'|')),fc_arr(v3_idx,HighOVIdx(2))
                  deallocate(v2,v2_pos)
                  
                  fc_sum1_tot=fc_arr(v1_idx,1)**2
                  fc_sum3_tot=fc_arr(v3_idx,1)**2
                  write(output_unit,*)
                  write(output_unit,'(A)')'FC Sums:'
                  write(output_unit,'("Class ", I2," - Gr: ",G10.3," %   ", "Ex: ",G10.3," %")')0,fc_sum1_tot*100,fc_sum3_tot*100
                  do ii = 1,fc_ss2%mc_v
                     fc_sum1=0d0
                     fc_sum3=0d0
                     do iii = fc_ss2%cs(ii)+1,fc_ss2%cs(ii+1)
                        fc_sum1=fc_sum1+fc_arr(v1_idx,iii)**2
                        fc_sum3=fc_sum3+fc_arr(v3_idx,iii)**2
                     end do
                     fc_sum1_tot=fc_sum1_tot+fc_sum1
                     fc_sum3_tot=fc_sum3_tot+fc_sum3
                     write(output_unit,'("Class ", I2," - Gr: ",G10.3," %   ", "Ex: ",G10.3," %")')ii,fc_sum1*100,fc_sum3*100
                  end do
                  write(output_unit,'("Total: "," - Gr: ",G10.3," %   ", "Ex: ",G10.3," %")')fc_sum1_tot*100,fc_sum3_tot*100
                  write(output_unit,*)
                  
                  
                  counts=0
                  do ii = 1,FC_ss1%fc_arr_size
                     !$OMP PARALLEL DO DEFAULT(NONE) &
                     !$OMP PRIVATE(iii) &
                     !$OMP SHARED(fc_ss2,ii,curThresh,fc_arr) &
                     !$OMP REDUCTION(+:counts)
                     do iii = 1,FC_ss2%fc_arr_size
                        if(abs(fc_arr(ii,iii))>=curThresh*100)then
                           counts(1)=counts(1)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*10)then
                           counts(2)=counts(2)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh)then
                           counts(3)=counts(3)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*1d-2)then
                           counts(4)=counts(4)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*1d-4)then
                           counts(5)=counts(5)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*1d-6)then
                           counts(6)=counts(6)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*1d-8)then
                           counts(7)=counts(7)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*1d-10)then
                           counts(8)=counts(8)+1
                        else if(abs(fc_arr(ii,iii))>=curThresh*1d-12)then
                           counts(9)=counts(9)+1
                        end if
                     end do
                     !$OMP END PARALLEL DO
                  end do
                  
                  write(output_unit,'(A,G10.4)')'Threshold: ',curThresh
                  write(output_unit,'(A,A16)')'THR    ','Count:'
                  write(output_unit,'(A,I16)')' * 100  : ',counts(1)
                  write(output_unit,'(A,I16)')' * 10   : ',counts(2)
                  write(output_unit,'(A,I16)')' * 1    : ',counts(3)
                  write(output_unit,'(A,I16)')' * 1d-2 : ',counts(4)
                  write(output_unit,'(A,I16)')' * 1d-4 : ',counts(5)
                  write(output_unit,'(A,I16)')' * 1d-6 : ',counts(6)
                  write(output_unit,'(A,I16)')' * 1d-8 : ',counts(7)
                  write(output_unit,'(A,I16)')' * 1d-10: ',counts(8)
                  write(output_unit,'(A,I16)')' * 1d-12: ',counts(9)
                  
                  
                  write(output_unit,*)
                  if(runtype=='SI2')then
                     call Modes_arr_t_dispose(modes_arr_v2,td%mc_ms)
                     !call FC_storSys_dispose(fc_ss1)
                     call FC_storSys_dispose(fc_ss2)
                     call FC_storSys_dusch_dispose(fc_ss_dusch)
                     !call FC_storSys_derivs_dispose(fc_ss_derivs)
                     deallocate(fc_arr)
                     cycle
                  end if
                  write(output_unit,'(A)')'#Polarizabilities...'
                  flush(output_unit)
                  ii=1
                  c_ms=0
                  time2=getTime()
                  allocate(v2(c_ms),v2_pos(c_ms))
                  
                  !<0|0> transition
                  fcsum1_un=0
                  fcsum3_un=0
                  canBeLow=.false.
                  call ExcM_v2_Polar_wps(ii,v3,v3_pos_short,v3_class,v3_idx,w3,v2,v2_pos,c_ms,v1,v1_pos_short,v1_class,v1_idx,w1,canBeLow_gr,canBeLow, &
                  td,fc_ss1,fc_ss2,fc_arr,curPolar,curPolar_anti,td%gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments,output_polContrs,alphafc, &
                  uncoup_modes_vmap,uncoup_modes_reddims,uncoup_modes_c,uncoup_fcarr,fcsum3_un,fcsum1_un,contr)
                  deallocate(v2,v2_pos)
                  fc_sum1_tot=fc_sum1_tot+fcsum1_un
                  fc_sum3_tot=fc_sum3_tot+fcsum3_un
                  if(uncoup_modes_c>0)then
                     write(output_unit,*)fc_sum1_tot*100
                     write(output_unit,*)fc_sum3_tot*100
                  end if
                  !class 1 transitions
                  allocate(exc_m_c1(td%nq,1))
                  exc_m_c1(:,1)=td%max_v2s(1)
                  !exc_m_c1(:,2)=td%max_v2s(1)
                  fc_ss2_c1=FC_storSys_make(exc_m_c1,1,td%nq,.true.,1,td%we)
                  fc_ss_dusch_c1=FC_storSys_Dusch_make(fc_ss1,fc_ss2_c1,td%A,td%B,td%C,td%D,td%E,td%nq)
                  !fc_ss_derivs_c1=FC_storSys_derivs_make(td,fc_ss2_c1)
                  fc_arr_c1=MakeFCArr_recalculate(fc_ss1,fc_ss2_c1,fc_ss_dusch_c1,td%fc_00,modes_arr_v2_c1,n_thr,check,.false.,.false.)
                  time1=GetTime()
                  write(output_unit,'("Class ",1X,I2,2X,I16," transitions",2X)',advance='no')1,(fc_ss2_c1%cs(2)-fc_ss2_c1%cs(1)+1)
                  flush(output_unit)
                  comb_n=size(modes_arr_v2_c1(1)%arr,dim=2)
                  allocate(v2(1),v2_pos(1))
                  fcsum1_un=0
                  fcsum3_un=0
                  do jj = 1,comb_n
                     v2_pos=modes_arr_v2_c1(1)%arr(:,jj)
                     ii=1
                     canBeLow(v2_pos(1))=.true.
                     call ExcM_v2_Polar_wps_c1(ii,v3,v3_pos_short,v3_class,v3_idx,w3,v2,v2_pos,v1,v1_pos_short,v1_class,v1_idx,w1,canBeLow_gr,canBeLow, &
                     td,fc_ss1,fc_ss2_c1,fc_arr_c1,curPolar,curPolar_anti,td%gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments,output_polContrs,alphafc, &
                     uncoup_modes_vmap,uncoup_modes_reddims,uncoup_modes_c,uncoup_fcarr,fcsum3_un,fcsum1_un,contr)
                     canBeLow(v2_pos(1))=.false.
                  end do
                  fc_sum1_tot=fc_sum1_tot+fcsum1_un
                  fc_sum3_tot=fc_sum3_tot+fcsum3_un
                  
                  time1=GetTime()-time1
                  write(output_unit,'(3X,A)')TR(GetTimeStr(time1))
                  if(uncoup_modes_c>0)then
                     write(output_unit,*)fc_sum1_tot*100
                     write(output_unit,*)fc_sum3_tot*100
                  end if
                  
                  flush(output_unit)
                  call Modes_arr_t_dispose(modes_arr_v2_c1,1)
                  call FC_storSys_dusch_dispose(fc_ss_dusch_c1)
                  call FC_storSys_dispose(fc_ss2_c1)
                  deallocate(exc_m_c1,fc_arr_c1,v2,v2_pos)
                  
                  !class 2 and higher transitions
                  do c_ms=2,td%mc_ms
                     time1=getTime()
                     write(output_unit,'("Class ",1X,I2,2X,I16," transitions",2X)',advance='no')c_ms,(fc_ss2%cs(c_ms+1)-fc_ss2%cs(c_ms)+1)
                     flush(output_unit)
                     allocate(v2(c_ms),v2_pos(c_ms))
                     ii=1
                     comb_n=size(modes_arr_v2(c_ms)%arr,dim=2)
                     fcsum1_un=0
                     fcsum3_un=0
                     !comb_n=nck(fc_ss2%red_n,c_ms)
                     !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
                     !$OMP SHARED(comb_n,v1_pos_short,v1_class,v3_pos_short,v3_class,i,w1,modes_arr_v2,c_ms,j,w3,td,fc_ss1,fc_ss2,wexc,nexc,st) &
                     !$OMP SHARED(v1_idx,v3_idx,fc_arr,fc_ss_derivs,fc_ss_dusch,ht,only_ht,iexc,canBeLow_gr,output_moments,output_polContrs,alphafc) &
                     !$OMP SHARED(uncoup_modes_vmap,uncoup_modes_reddims,uncoup_fcarr,uncoup_modes_c,w_ps,contr) &
                     !$OMP FIRSTPRIVATE(v1,v3,canBeLow,v2,v2_pos) &
                     !$OMP PRIVATE(ii,jj,kk) &
                     !$OMP REDUCTION(+:curpolar,curpolar_anti,fcsum1_un,fcsum3_un)
                     outer: do jj=1,comb_n !goddamn, I learn something new every day
                        ii=1
                        !success=.true.
                        v2_pos=modes_arr_v2(c_ms)%arr(:,jj)
                        do kk=1,c_ms
                           canBeLow(fc_ss2%v_map(v2_pos(kk)))=.true.
                        end do
                        call ExcM_v2_Polar_wps(ii,v3,v3_pos_short,v3_class,v3_idx,w3,v2,v2_pos,c_ms,v1,v1_pos_short,v1_class,v1_idx,w1,canBeLow_gr,canBeLow, &
                        td,fc_ss1,fc_ss2,fc_arr,curPolar,curpolar_anti,td%gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments,output_polContrs,alphafc, &
                        uncoup_modes_vmap,uncoup_modes_reddims,uncoup_modes_c,uncoup_fcarr,fcsum3_un,fcsum1_un,contr)
                        do kk=1,c_ms !what will canBeLow look like after exiting this OMP directive?
                           canBeLow(fc_ss2%v_map(v2_pos(kk)))=.false.
                        end do
                     end do outer
                     !$OMP END PARALLEL DO
                     fc_sum1_tot=fc_sum1_tot+fcsum1_un
                     fc_sum3_tot=fc_sum3_tot+fcsum3_un
                     !call Polar_new(curpolar)
                     !polars(i,j)=curPolar
                     deallocate(v2,v2_pos)
                     time1=getTime()-time1
                     write(output_unit,'(3X,A)')TR(GetTimeStr(time1))
                     if(uncoup_modes_c>0)then
                        write(output_unit,*)fc_sum1_tot*100
                        write(output_unit,*)fc_sum3_tot*100
                     end if
                  end do
                  call SymmetrizePolarAAc(curPolar)
                  curPolar%ap=conjg(curPolar%ap)
                  curPolar%G=conjg(curPolar%G)
                  curPolar%Gc=conjg(curPolar%Gc)
                  curPolar%A=conjg(curPolar%A)
                  curPolar%Ac=conjg(curPolar%Ac)

                  !deallocate(canBeLow)
                  time2=getTime()-time2
                  write(output_unit,'(3X,A)')TR(GetTimeStr(time2))
                  flush(output_unit)
                  tr_idx=tr_idx+1
                  
                  
                  
                  
                  call Modes_arr_t_dispose(modes_arr_v2,td%mc_ms)
                  call FC_storSys_dispose(fc_ss2)
                  call FC_storSys_dusch_dispose(fc_ss_dusch)
                  ! call FC_storSys_derivs_dispose(fc_ss_derivs)
                  deallocate(fc_arr)
                  a_arr=0d0
                  if(IsFundamentalTransition(v1,v1_pos_mapped,v1_class,v3,v3_pos_mapped,v3_class) .and. fundLead)then
                     cur_tr%vi=v1
                     cur_tr%vi_pos=v1_pos_mapped
                     cur_tr%vf=v3
                     cur_tr%vf_pos=v3_pos_mapped
                     cur_tr%wif=(w3-w1)*au_2_cm
                     if(rroa_spr_do(12))then
                        cur_tr%ram=a_arr(1,12)
                        cur_tr%roa=a_arr(2,12)
                     else
                        cur_tr%ram=maxval(a_arr(1,:))
                        cur_tr%roa=maxval(abs(a_arr(2,:))) !todo, not abs
                     end if
                     call InsertTransition(trs_fund,nq,trs_fund_c,cur_tr,iexc,nexc)
                  end if
                  if(write_ten)then
                     if(.not.w_ps)then
                        do i_iexc=1,nexc
                           call WriteTenPretty(1000+i_iexc,wexc(i_iexc),PolarExc2Polar(curPolar,i_iexc),v1,v1_pos_mapped,v1_class,w1*au_2_cm,v3,v3_pos_mapped,v3_class,w3*au_2_cm)
                        end do
                     else
                        call WriteTenPretty(1000+iexc,wexc(iexc),PolarExc2Polar(curPolar,iexc),v1,v1_pos_mapped,v1_class,w1*au_2_cm,v3,v3_pos_mapped,v3_class,w3*au_2_cm)
                     end if
                  end if
                  if(output_polContrs)then
                     if(.not.w_ps)then
                        do i_iexc=1,nexc
                           write(op_polcontrs_unit-1+i_iexc,*)
                           flush(op_polcontrs_unit-1+i_iexc)
                        end do
                     else
                        write(op_polcontrs_unit-1+iexc,*)
                        flush(op_polcontrs_unit-1+iexc)
                     end if
                  end if
                  if(uncoup_modes_c>0)deallocate(uncoup_fcarr)
                  if(.not.w_ps)exit
               end do
               deallocate(v3_pos_short)
               deallocate(v1_pos_short)
               if(runtype=='DEF')call FC_storSys_dispose(fc_ss1)
               ! if(.not.doFin)then
                  ! trs_c=trs_c-1
                  ! deallocate(v3_pos_short)
                  ! deallocate(v1_pos_short)
                  ! cycle
               ! end if
               mode_found=.true.
               pol_idx_v1=FindMode(ground_polars_v1,gp_v1_c,pol_idx_v1_last,v1,v1_pos_mapped,v1_class)
               if(pol_idx_v1==0)then
                  pol_idx_v1_last=pol_idx_v1_last+1
                  pol_idx_v1=pol_idx_v1_last
                  ground_polars_v1(pol_idx_v1_last)=MakeMode(v1,v1_pos_mapped,v1_class)
                  mode_found=.false.
               end if
               pol_idx_v3=FindMode(ground_polars_v3,gp_v3_c,pol_idx_v3_last,v3,v3_pos_mapped,v3_class)
               if(pol_idx_v3==0)then
                  pol_idx_v3_last=pol_idx_v3_last+1
                  pol_idx_v3=pol_idx_v3_last
                  ground_polars_v3(pol_idx_v3_last)=MakeMode(v3,v3_pos_mapped,v3_class)
                  mode_found=.false.
               end if
               
               if(output_polars)then
                  do ii = 1,nexc
                     write(670,*)td%root,TR(str_v1),TR(str_v3),(w3-w1)*au_2_cm," cm-1, ",1d7/(wexc(ii)*au_2_cm)," nm"
                     write(670,*)'Ap ',curPolar%Ap(1,1,ii),curPolar%Ap(1,2,ii),curPolar%Ap(1,3,ii),curPolar%Ap(2,1,ii),curPolar%Ap(2,2,ii),curPolar%Ap(2,3,ii),curPolar%Ap(3,1,ii),curPolar%Ap(3,2,ii),curPolar%Ap(3,3,ii)
                     write(670,*)'G  ',curPolar%G(1,1,ii),curPolar%G(1,2,ii),curPolar%G(1,3,ii),curPolar%G(2,1,ii),curPolar%G(2,2,ii),curPolar%G(2,3,ii),curPolar%G(3,1,ii),curPolar%G(3,2,ii),curPolar%G(3,3,ii)
                     write(670,*)'Gc ',curPolar%Gc(1,1,ii),curPolar%Gc(1,2,ii),curPolar%Gc(1,3,ii),curPolar%Gc(2,1,ii),curPolar%Gc(2,2,ii),curPolar%Gc(2,3,ii),curPolar%Gc(3,1,ii),curPolar%Gc(3,2,ii),curPolar%Gc(3,3,ii)
                     do jj = 1,3
                        write(670,*)'A  ',curPolar%A(jj,1,1,ii),curPolar%A(jj,1,2,ii),curPolar%A(jj,1,3,ii),curPolar%A(jj,2,1,ii),curPolar%A(jj,2,2,ii),curPolar%A(jj,2,3,ii),curPolar%A(jj,3,1,ii),curPolar%A(jj,3,2,ii),curPolar%A(jj,3,3,ii)
                     end do
                     do jj = 1,3
                        write(670,*)'Ac ',curPolar%Ac(jj,1,1,ii),curPolar%Ac(jj,1,2,ii),curPolar%Ac(jj,1,3,ii),curPolar%Ac(jj,2,1,ii),curPolar%Ac(jj,2,2,ii),curPolar%Ac(jj,2,3,ii),curPolar%Ac(jj,3,1,ii),curPolar%Ac(jj,3,2,ii),curPolar%Ac(jj,3,3,ii)
                     end do
                  end do
               end if
               
               if(mode_found)then !lazy workaround
                  if(allocated(ground_polars(pol_idx_v1,pol_idx_v3)%ap))then
                     ground_polars(pol_idx_v1,pol_idx_v3)=ground_polars(pol_idx_v1,pol_idx_v3)+curPolar
                  else
                     call AllocatePolar(ground_polars(pol_idx_v1,pol_idx_v3),nexc)
                     ground_polars(pol_idx_v1,pol_idx_v3)=curPolar
                     tr_count=tr_count+1
                  end if
               else
                  call AllocatePolar(ground_polars(pol_idx_v1,pol_idx_v3),nexc)
                  ground_polars(pol_idx_v1,pol_idx_v3)=curPolar
                  tr_count=tr_count+1
               end if
               !deallocate(canBeLow_gr)
            end do
         end do
         deallocate(v1_modes)
         deallocate(mode_psed,mode_psed_from)
         if(write_fc_sys)then
            close(10)
         end if
         deallocate(exc_m1,modes_arr_v1,v1)
         deallocate(exc_m2,thr_arr,canBeLow_gr,canBeLow)
         write(output_unit,'(A)')'------------------------------'
         if(fundLead)then
            deallocate(trs_fund)
         end if
         
         if(uncoup_lim<1d0)then
            deallocate(uncoup_modes_excm,uncoup_modes_vmap,uncoup_modes_reddims,uncoup_modes)
         else
            deallocate(uncoup_modes)
         end if
      end do
      if(write_ten)then
         do iexc=1,nexc
            close(1000+iexc)
         end do
      end if
      
      call Dispose_NCK_arr()
      if(runtype == 'SI2')return
      if(.not.td_approach)then
      write(output_unit,'(A)')'Summing polarizabilities...'
      flush(output_unit)
      tr_idx=1
      do i = 1,pol_idx_v1_last
         v1=ground_polars_v1(i)%v
         v1_pos=ground_polars_v1(i)%v_pos
         !v1_pos_mapped=MapVPos(v1_pos,v1_class,fc_ss_gr%v_map,fc_ss_gr%red_n)
         v1_class=size(v1,dim=1)
         w1=getFreq_short(v1,v1_pos,wgg,v1_class)*cm_2_au
         if(v1_class==0)then
            str_v1='|0 >'
         else
            str_v1=FC2Str_new(v1,v1_pos,v1_class,.true.,.false.)
         end if
         kBT=exp(-h*cc*w1*au_2_cm*100d0/(kb*temp)) !check this
         do j = 1,pol_idx_v3_last
            v3=ground_polars_v3(j)%v
            v3_pos=ground_polars_v3(j)%v_pos
            !v3_pos_mapped=MapVPos(v3_pos,v3_class,fc_ss_gr%v_map,fc_ss_gr%red_n)
            v3_class=size(v3,dim=1)
            w3=getFreq_short(v3,v3_pos,wgg,v3_class)*cm_2_au
            kbt_anti=exp(-h*cc*w3*au_2_cm*100d0/(kb*temp))
            curPolar=ground_polars(i,j)
            if(.not. allocated(curPolar%ap))cycle
            
            
            str_v3=FC2Str_new(v3,v3_pos,v3_class,.true.,.false.)
            do iexc = 1,nexc
               if(v1_class==0 .and. v3_class==1 .and. add_pol_coeff /= 0d0)then
                  if(v3(1)==1)then
                     curPolar%ap(:,:,i)=curPolar%ap(:,:,i)+add_pol_coeff*polar_ders(v3_pos(1))%ap(:,:,iexc)
                     curPolar%G(:,:,i)=curPolar%G(:,:,i)+add_pol_coeff*polar_ders(v3_pos(1))%G(:,:,iexc)
                     curPolar%A(:,:,:,i)=curPolar%A(:,:,:,i)+add_pol_coeff*polar_ders(v3_pos(1))%A(:,:,:,iexc)
                  end if
               end if
               rroa_spr_temp=0d0
               if(wr_is_e00)then
                  wr=sum(tds(:)%e_00)/dble(td_n)
               else
                  wr=wexc(iexc)-w3
               end if
               if(alphaFC)then
                  call wrram3(wexc(iexc),wr,temp,spectrum_temp,w1,w3,TINY(1d0),tr_idx,dble(wmin),dble(wmax),npx,use_gauss,FWHM,.true.,rroa_inv,&
                               realpart(curPolar%Ap(:,:,iexc)),imagpart(curPolar%Ap(:,:,iexc)), &
                               realpart(curPolar%G(:,:,iexc)),imagpart(curPolar%G(:,:,iexc)), &
                               realpart(curPolar%Gc(:,:,iexc)),imagpart(curPolar%Gc(:,:,iexc)), &
                               realpart(curPolar%A(:,:,:,iexc)),imagpart(curPolar%A(:,:,:,iexc)), &
                               realpart(curPolar%Ac(:,:,:,iexc)),imagpart(curPolar%Ac(:,:,:,iexc)), &
                               rroa_spr_temp,kBT,useA,useG,rroa_spr_do,rroa_tabOutputs(:,iexc),str_v1,str_v3,&
                               realpart(curPolar%Ap_fcht(:,:,iexc)),imagpart(curPolar%Ap_fcht(:,:,iexc)))
               else
                  call wrram3(wexc(iexc),wr,temp,spectrum_temp,w1,w3,TINY(1d0),tr_idx,dble(wmin),dble(wmax),npx,use_gauss,FWHM,.true.,rroa_inv,&
                               realpart(curPolar%Ap(:,:,iexc)),imagpart(curPolar%Ap(:,:,iexc)), &
                               realpart(curPolar%G(:,:,iexc)),imagpart(curPolar%G(:,:,iexc)), &
                               realpart(curPolar%Gc(:,:,iexc)),imagpart(curPolar%Gc(:,:,iexc)), &
                               realpart(curPolar%A(:,:,:,iexc)),imagpart(curPolar%A(:,:,:,iexc)), &
                               realpart(curPolar%Ac(:,:,:,iexc)),imagpart(curPolar%Ac(:,:,:,iexc)), &
                               rroa_spr_temp,kBT,useA,useG,rroa_spr_do,rroa_tabOutputs(:,iexc),str_v1,str_v3,&
                               realpart(curPolar%Ap(:,:,iexc)),imagpart(curPolar%Ap(:,:,iexc)))
               end if
               if(write_inv)then
                  rroa_inv2=MakeLameInvariants(curPolar%Ap(:,:,iexc),curPolar%G(:,:,iexc),curPolar%Gc(:,:,iexc),curPolar%A(:,:,:,iexc),curPolar%Ac(:,:,:,iexc),wexc(iexc))
                  write(99,'(A13,1X,A13,1X,F8.2,13(1X,G13.6),F7.2," nm")')TR(str_v1),TR(str_v3),(w3-w1)*au_2_cm,rroa_inv,10d0**7/wexc(iexc)*cm_2_au
                  write(66,'(A13,1X,A13,1X,F8.2,30(1X,G13.6),F7.2," nm")')TR(str_v1),TR(str_v3),(w3-w1)*au_2_cm,rroa_inv2,10d0**7/wexc(iexc)*cm_2_au
               end if
               
               rroa_spr(:,:,:,iexc)=rroa_spr(:,:,:,iexc)+rroa_spr_temp
            end do
            tr_idx=tr_idx+1
         end do
      end do
      if(write_inv)then
         close(99)
      end if
      if(output_polContrs)then
         do iexc=1,nexc
            close(op_polcontrs_unit-1+iexc)
         end do
      end if
      ! do i = 1,ns0
         ! if(.not.rroa_spr_do(i))cycle
         ! do j = 1,nexc
            ! c80='RAM'//TR(rroa_exp(i))//'_'//TR(wexc_nm_str(j))//'.PRN'
            ! c80_2='ROA'//TR(rroa_exp(i))//'_'//TR(wexc_nm_str(j))//'.PRN'
            ! open(77,file=c80)
            ! open(78,file=c80_2)
            ! do k = 1,npx
               ! w=dble(wmax-wmin)/npx*(k-1)+wmin
               ! write(77,'(f14.4,g25.12)')w,rroa_spr(k,1,i,j)
               ! write(78,'(f14.4,g25.12)')w,rroa_spr(k,2,i,j)
            ! end do
            ! close(77)
            ! write(output_unit,*)'Written file: '//c80
            ! close(78)
            ! write(output_unit,*)'Written file: '//c80_2
         ! end do
      ! end do
      deallocate(v1,v1_pos)
      deallocate(v3,v3_pos)
      end if
      if(.not.TD_approach)then
         deallocate(exc_m_union)
         call deallocatePolar(curPolar)
      end if
      deallocate(rroa_spr,rroa_inv,rroa_tabOutputs)
      deallocate(rroa_spr_temp)
      deallocate(ground_polars)
      deallocate(ground_polars_v1)
      deallocate(ground_polars_v3)
      !summedPols=SumUpTDPols(tds,tdn)
   end subroutine DoRROA
   
   function FindTransition_R_c1(trs_arr,n,trs_arr_c,v1,v1_pos,v1_class,v3,v3_pos,v3_class,iexc,nexc)result(idx)
      integer n,nexc,v1(v1_class),v1_class,v3(v3_class),v3_class,trs_arr_c,iexc
      integer(int16) v1_pos(v1_class),v3_pos(v3_class)
      type(Transition_R) trs_arr(n,nexc),cur_tr
      
      integer i,idx,v1_class_r,v3_class_r
      
      idx=0
      ! v1_class_r=count(v1>0)
      ! v3_class_r=count(v3>0)
      do i = 1,trs_arr_c
         cur_tr=trs_arr(i,iexc)

         if(v1_class /= count(cur_tr%vi>0) .or. v3_class /= count(cur_tr%vf>0))cycle
         if(v1_class==0)then
            if(Equal_IArr(v3,cur_tr%vf))then
               if(Equal_I2Arr(v3_pos,cur_tr%vf_pos,v3_class))then
                  idx=i
                  return
               end if
            end if
         else
            if(Equal_IArr(v1,cur_tr%vi) .and. Equal_IArr(v3,cur_tr%vf))then
               if(Equal_I2Arr(v1_pos,cur_tr%vi_pos,v1_class) .and. Equal_I2Arr(v3_pos,cur_tr%vf_pos,v3_class))then
                  idx=i
                  return
               end if
            end if
         end if
      end do
   end function FindTransition_R_c1
   
   subroutine InsertTransition(trs_arr,n,trs_arr_c,tr,iexc,nexc)
      integer iexc,nexc,n,trs_arr_c
      type(Transition_R) trs_arr(n,nexc),tr,cur_tr,cur_tr2
      
      integer i,j
      
      if(trs_arr_c==0)then
         trs_arr_c=trs_arr_c+1
         trs_arr(trs_arr_c,iexc)=tr
      else
         do i = 1,trs_arr_c
            cur_tr=trs_arr(i,iexc)
            if(abs(tr%roa)>abs(cur_tr%roa) .or. tr%ram>cur_tr%ram)then
               do j = trs_arr_c+1,i+1,-1
                  trs_arr(j,iexc)=trs_arr(j-1,iexc)
               end do
               trs_arr(i,iexc)=tr
               trs_arr_c=trs_arr_c+1
               return
            end if
         end do
         trs_arr_c=trs_arr_c+1
         trs_arr(trs_arr_c,iexc)=tr
      end if
   end subroutine InsertTransition
   
   function FindMode(modes,modes_n,lastModeIdx,v,v_pos,v_n)result(idx)
      integer v_n,v(v_n)
      integer(int16) :: v_pos(v_n)
      integer modes_n,idx,lastModeIdx
      type(Mode) modes(modes_n)
      
      integer i
      
      idx=0
      do i = 1,lastModeIdx
         if(size(modes(i)%v_pos,dim=1)==v_n)then
            if(Equal_I2Arr(modes(i)%v_pos,v_pos,v_n).and.Equal_IArr(modes(i)%v,v))then
               idx=i
               return
            end if
         end if
      end do
   end function FindMode
   
   function IsFundamentalTransition(v1,v1_pos,v1_class,v3,v3_pos,v3_class)result(res)
      integer v1_class,v3_class
      integer v1(v1_class),v3(v3_class)
      integer(int16) v1_pos(v1_class),v3_pos(v3_class)
      logical res
      
      res=.false.
      if(v1_class > 0 .or. v3_class > 1)return
      if(v3(1)>1)return
      res=.true.
   end function IsFundamentalTransition
   
   function IsAllowedTransition(v1,v1_pos,v1_class,v3,v3_pos,v3_class)result(res)
      integer v1_class,v3_class
      integer v1(v1_class),v3(v3_class)
      integer(int16) v1_pos(v1_class),v3_pos(v3_class)
      logical res
      integer idx
      
      res=.false.
      if(abs(v3_class-v1_class)>1)return
      select case(v1_class)
         case(0)
            res=.true.
            if(v3_class>1)then ! |0> to combination mode
               res=.false.
               return
            end if
            if(v3(1)>1)then ! |0> to higher overtone (eg. to v^3, v^4...)
               res=.false.
               return
            end if
         case(1)
            res=.true.
            if(v3_class>2)then !I assume that class 1 -> class 3 transition is forbidden
               res=.false.
               return
            end if
            idx=findloc(v3_pos,v1_pos(1),dim=1) 
            if(idx>0)then 
               if(v1(1)-v3(idx)==1)return !|1^2> --> |1^1 2^1>
            end if
            if(.not.ModeIsSubsetOf(v1,v1_pos,v1_class,v3,v3_pos,v3_class))then !other combination transitions
               res=.false.
               return
            end if
         case(2) !me not have clue
            res=.false.
            
         case default
            res = .false. !Yeah, I have not gotten that far yet
      end select
      
   end function IsAllowedTransition
   
   function ModeIsSubsetOf(v1,v1_pos,v1_class,v3,v3_pos,v3_class)result(res)
      integer v1_class,v3_class
      integer v1(v1_class),v3(v3_class)
      integer(int16) v1_pos(v1_class),v3_pos(v3_class)
      integer i,idx
      logical res
      
      if(v1_class==v3_class)then
         res=.true.
         do i = 1,v1_class
            idx=findloc(v3_pos,v1_pos(i),dim=1)
            if(idx>0)then
               if(v1(i)<=v3(idx))then
                  cycle
               end if
            end if
            res=.false.
            return
         end do
      else if(v1_class<v3_class)then
         res=.true.
         do i = 1,v1_class
            idx=findloc(v3_pos,v1_pos(i),dim=1)
            if(idx>0)then
               if(v1(i)<=v3(idx))then
                  cycle
               end if
            end if
            res=.false.
            return
         end do
      else if(v1_class>v3_class)then
         res=.true.
         do i = 1,v3_class
            idx=findloc(v1_pos,v3_pos(i),dim=1)
            if(idx>0)then
               if(v3(i)<=v1(idx))then
                  cycle
               end if
            end if
            res=.false.
            return
         end do
      end if
   end function ModeIsSubsetOf
   
   subroutine WriteFCSys(unitt,fc_ss1)
      integer unitt,ii,jj
      type(FC_storSys) fc_ss1
      
      write(unitt,'("MAX_CLASS: ",I2,"  REDUCED_NQ: ",I5)')fc_ss1%mc_v,fc_ss1%red_n
      write(unitt,'("V_MAP:")')
      do jj = 1,fc_ss1%red_n
         write(unitt,'(I4,1X)',advance='no')fc_ss1%v_map(jj)
      end do
      write(unitt,*)
      write(unitt,*)
      
      do ii = 1,fc_ss1%mc_v
         write(unitt,'("CLASS ",I2)')ii
         do jj = 1,fc_ss1%red_n
            write(unitt,'(I3,1X)',advance='no')fc_ss1%red_dims(jj,ii)
         end do
         write(unitt,*)
         write(unitt,*)
      end do
      write(unitt,*)
      write(unitt,*)'CLASS STARTS: '
      do ii = 1,fc_ss1%mc_v+1
         write(unitt,'(I19,1X)',advance='no')fc_ss1%cs(ii)
      end do
      write(unitt,*)
   end subroutine WriteFCSys
   
   
   function GroundExcM(td)result(exc_m1)
      type(ExcState) td
      integer :: c_gfs,i
      double precision freq
      integer,allocatable :: exc_m1(:,:)
      
      allocate(exc_m1(td%nq,td%mc_gfs))
      exc_m1=0
      
      do c_gfs = 1,td%mc_gfs
         do i = 1,td%nq
            if(td%mc_gs<c_gfs)then
               exc_m1(i,c_gfs)=td%max_v3s(c_gfs)
               cycle
            else if(td%mc_fs<c_gfs)then
               exc_m1(i,c_gfs)=td%max_v1s(c_gfs)
               cycle
            end if
            exc_m1(i,c_gfs)=max(td%max_v1s(c_gfs),td%max_v3s(c_gfs))
            
            !exc_m1(i,c_gfs)=max(exc_m1(i,c_gfs),td%min_v1s(c_gfs)) !what if min>max? :(
         end do
      end do
      
      !get rid of low frequency and high frequency modes
      !do i = 1,td%nq
         !freq=td%wg(i)*au_2_cm
         !if(freq<min_freq_mode .or. freq>max_freq_mode)exc_m1(i,:)=0
      !end do
   end function GroundExcM
   
   recursive subroutine Excm_v1_modes(i,v1,v1_pos,v1_clas,v1_modes,fc_ss1,td)
      type(ExcState) td
      integer i,ii,j,v1_clas,v1(v1_clas),v1_pos_normal(v1_clas)
      integer(int16) v1_pos(v1_clas),v1_pos_mapped(v1_clas)
      type(FC_storSys) fc_ss1
      type(Mode) v1_modes(:)
      integer :: v1_idx
      
      
      !<v| is set up
      if(i>v1_clas)then !Set final state mode excitations
         v1_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_clas,fc_ss1)
         !v1_pos_normal=ShortPos2NormPos(v1_pos,v1_clas,fc_ss1%v_map,fc_ss1%red_n)
         v1_modes(v1_idx)=MakeMode(v1,MapVPos(v1_pos,v1_clas,fc_ss1%v_map,fc_ss1%red_n),v1_clas)
         return
      end if
      
      ii=v1_clas-i+1
      do j = 1,fc_ss1%red_dims(v1_pos(ii),v1_clas)
         v1(ii)=j
         call Excm_v1_modes(i+1,v1,v1_pos,v1_clas,v1_modes,fc_ss1,td)
      end do
   end subroutine Excm_v1_modes
   
   !Determine Uncoupled Modes, EXcited state
   function DUM_EX(J,NQ,crit_sq,thr)result(isUncoupled)
      integer NQ,c,i
      double precision J(NQ,NQ),thr
      logical crit_sq,isUncoupled(NQ)
      
      isUncoupled=.false.
      if(crit_sq)then !the metric is square of the J-matrix element
         do i = 1,NQ
            if(J(i,i)**2>thr)isUncoupled(i)=.true.
         end do
      else !or the absolute value of J-matrix element
         do i = 1,NQ
            if(abs(J(i,i))>thr)isUncoupled(i)=.true.
         end do
      end if
   end function DUM_EX
   
   
   function DEM_1v_pickThr(v,v_n,thr,thr_0,thr_ov,thr_comb,thr_other)result(thr_pick)
      integer v_n,v(v_n),v_class
      double precision thr,thr_0,thr_ov,thr_comb,thr_other,thr_pick
      
      v_class=count(v/=0)
      select case(v_class)
         case(0)
            thr_pick=thr_0
         case(1)
            if(v(1)>1)then
               thr_pick=thr_ov
            else
               thr_pick=thr
            end if
         case(2) !I dunno yet about combination overtones
            thr_pick=thr_comb
         case default
            thr_pick=thr_other
      end select
   end function DEM_1v_pickThr
   
   !Determine (excited state) Excitation Maxima for 1 vibrational ground state
   !ie. prescreening of excited state modes for one ground state mode v1
   !done numerically
   function DEM_1v_RROA(v1,v1_pos,v1_n,max_v2s,mc_ms_tmp,td,thr,check,uncoup_m,w_ps,wn,wexc,nexc,n_thr,gamma,fc_sum_min,fc_sum_min_its,ignore_imag_modes)result(exc_m_v)
      type(ExcState),intent(in) :: td
      logical,intent(in) :: check,uncoup_m(td%nq),w_ps
      integer,allocatable :: exc_m2(:,:)
      integer,intent(in) :: mc_ms_tmp,nexc
      integer,intent(in) :: max_v2s(mc_ms_tmp),fc_sum_min_its
      double precision,intent(in) :: fc_sum_min
      
      type(FC_storSys):: fc_ss1,fc_ss2
      type(FC_storSys_dusch) :: FC_ss_dusch
      
      integer(1),allocatable :: fc_arr_check(:,:)
      double precision,allocatable :: FC_arr(:,:)
      
      integer,intent(in) :: v1_n
      integer :: v1_class,v3_class,n_thr
      integer,intent(inout) :: v1(v1_n)
      integer(int16),intent(inout) :: v1_pos(v1_n)
      integer :: v1_nz(v1_n)
      integer(int16),allocatable :: v1_pos_new(:),v1_pos_mem(:)
      integer,allocatable :: exc_m1(:,:),exc_m_v(:,:,:)
      integer,allocatable :: v1_mem(:)
      integer :: red_n
      double precision thr(nexc),summ,wn,wm,fc,wexc(nexc),w1,gamma,fc_sum(nexc),fc_sum_tot
      logical eqna,success(nexc),ignore_imag_modes
      double precision time1,time2,tot_size
      
      type(modes_arr_t),allocatable :: modes_arr_v2(:)
      integer(int64),allocatable :: combs_n(:)
      type(v_col),allocatable :: v2_arr(:)
      integer,allocatable :: v2(:)
      integer(int16),allocatable :: v2_pos(:)
      
      integer i,ii,ic,j,its
      integer,parameter :: ic_thr=5
      integer(int64) comb_n,i_64
      integer c_ms,c_gs,mc_gfs
      
      
      time1=gettime()
      write(output_unit,*)
      !write(output_unit,*)

      !thr=td%thr_v2
      !call v_union(v1,v1_pos,v1_n,v1_class,v3,v3_pos,v3_n,v3_class,v13,v13_pos,v13_n,exc_m1)
      v1_class=count(v1/=0)
      
      mc_gfs=v1_class
      wm=0d0
      if(v1_class>0)then
         allocate(exc_m1(td%nq,v1_class))
         exc_m1=0
         do i = 1,v1_class
            exc_m1(v1_pos(i),:)=v1(i)
            wm=wm+td%wg(v1_pos(i))*v1(i)
         end do
         wm=wm*au_2_cm
         fc_ss1=FC_storSys_make(exc_m1,v1_class,td%nq,.true.,mc_gfs)
         v1_pos_new=GetSeq(v1_class)
      else
         fc_ss1=FC_storSys_make_0(td%nq)
      end if
      write(output_unit,'(A)')'---'
      write(output_unit,'("State: ",A,1X,F7.2," cm-1",1X,"From: ",F7.2," cm-1")')TR(FC2Str_new(v1,v1_pos,v1_n,.true.,.false.)),wm,wn*au_2_cm
      flush(output_unit)
      allocate(exc_m2(td%nq,mc_ms_tmp))
      !allocate(v1_modes(fc_ss1%fc_arr_size))
      exc_m2=0
      do i = 1,mc_ms_tmp
         do j = 1,td%nq
            !if(dabs(td%D(j))>1 .or. dabs(td%E(j,j))>1.79)then
               exc_m2(j,i)=max_v2s(i)
               ! if(ignore_imag_modes)then
                  ! if(td%TMExpand=='EE' .and. td%im_modes_ex_ex(j))exc_m2(j,i)=0
                  ! if(td%TMExpand=='GG' .and. td%im_modes_ex_gr(j))exc_m2(j,i)=0
                  ! if((td%TMExpand=='GE' .or. td%TMExpand=='EG') .and. &
                  ! (td%im_modes_ex_gr(j) .or. td%im_modes_ex_ex(j)))exc_m2(j,i)=0
               ! end if
            !else
            !   exc_m2(j,i)=min(3,max_v2s(i))
            !end if
         end do
      end do
      
      fc_ss2=FC_storSys_make_uncoup(exc_m2,mc_ms_tmp,td%nq,.true.,mc_ms_tmp,uncoup_m,td%we)
      FC_ss_dusch=FC_storSys_Dusch_make(fc_ss1,fc_ss2,td%A,td%B,td%C,td%D,td%E,td%nq)
      
      
      ! max_v_tot=maxval(max_v2s)
      ! if(v1_n==1)then
         ! allocate(fc_arr(FC_ss1%fc_arr_size,max_v_tot+1))
         ! do i = 1,nq
            ! if(.not.uncoup_m(i))cycle
            
            
         ! end do
         ! deallocate(fc_arr)
      ! else
      
      ! end if
      
      !For some reason the whole fc_arr goes into RES memory with increasing class, even though it does not need to, it should all be in virtual memory
      !I guess this is because arrays in gfortran are contiguous and there are no memory gaps between elements
      tot_size=fc_ss1%fc_arr_size*fc_ss2%fc_arr_size*kind(fc_arr)
      write(output_unit,*)TR(GetSizeStr_D(tot_size))
      flush(output_unit)
      allocate(fc_arr(fc_ss1%fc_arr_size,fc_ss2%fc_arr_size))
 
#ifdef DEBUG
      if(check)then
         allocate(fc_arr_check(fc_ss1%fc_arr_size,fc_ss2%fc_arr_size))
         fc_arr_check=0
         fc_arr_check(1,1)=1
         fc_arr=HUGE(1d0)
      end if
#endif

      fc_arr(1,1)=td%fc_00
      ! if(w_ps)then
         allocate(exc_m_v(td%nq,mc_ms_tmp,nexc))
      ! else
         ! allocate(exc_m_v(td%nq,mc_ms_tmp,1))
      ! end if
      exc_m_v=0

      allocate(combs_n(mc_ms_tmp))
      allocate(v2_arr(mc_ms_tmp))
      do c_ms=1,mc_ms_tmp
         allocate(v2_arr(c_ms)%arr(c_ms))
         comb_n=nck(fc_ss2%red_n,c_ms)
         combs_n(c_ms)=comb_n
      end do
      
      fc_arr=MakeFCArr_recalculate(fc_ss1,fc_ss2,fc_ss_dusch,td%fc_00,modes_arr_v2,n_thr,check,.false.,.false.)
      fc_sum_tot=0d0
      do ii = 1,fc_ss2%fc_arr_size
         fc_sum_tot=fc_sum_tot+fc_arr(FC_ss1%fc_arr_size,ii)**2
      end do
      write(output_unit,*)'ps'
      flush(output_unit)
      if(v1_class>0)then
         
         its=1
         success=.false.
         do while(.not.all(success)  .and. its<=fc_sum_min_its)
            ii=1
            fc_sum=0d0
            exc_m_v=0
            call Excm_v1_FC_inner_detExcM_RROA(ii,v1,v1_pos_new,v1_class,v2_arr,modes_arr_v2,combs_n,mc_ms_tmp,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,exc_m_v,thr,ic_thr,w_ps,wn,wexc,nexc,td%e_00,gamma,fc_sum)
            do ii=1,nexc
               success(ii)=fc_sum(ii)*100>=fc_sum_min
               if(.not.success(ii))thr(ii)=AdjustThreshold(thr(ii))
               if(.not.w_ps)exit
            end do
            its=its+1
            if(.not.all(success) .and. its<=fc_sum_min_its .and. all(abs(fc_sum-fc_sum_tot)>1d-6))then
               write(output_unit,*)'Iter ',its
#ifdef DEBUG
               if(check)then
                  fc_arr_check=0
                  fc_arr_check(1,1)=1
               end if
#endif
            end if
         end do
         if(its>fc_sum_min_its)write(Output_unit,*)'Reached maximum iterations: ',its-1
      else
         its=1
         success=.false.
         do while(.not.all(success) .and. its<=fc_sum_min_its)
            fc_sum=fc_arr(1,1)**2
            exc_m_v=0
            do c_ms=1,mc_ms_tmp
               comb_n=nck(fc_ss2%red_n,c_ms)
               allocate(v2(c_ms),v2_pos(c_ms))
               !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
               !$OMP PRIVATE(ic,ii,i_64) &
               !$OMP FIRSTPRIVATE(v2,v2_pos) &
               !$OMP SHARED(modes_arr_v2,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,thr) &
               !$OMP SHARED(nexc,wexc,wn,w_ps,comb_n,c_ms,td,gamma) &
               !$OMP REDUCTION(MAX:exc_m_v) &
               !$OMP REDUCTION(+:fc_sum)
               do i_64 = 1,comb_n
                  ii=1
                  ic=0
                  v2_pos=modes_arr_v2(c_ms)%arr(:,i_64)
                  call ExcM_v2_FC_0_inner_detExcM_RROA(ii,v2,v2_pos,c_ms,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check, &
                  FC_ss1%fc_arr_size,fc_ss2%fc_arr_size,check,exc_m_v,thr,ic,ic_thr,w_ps,wn,wexc,nexc,td%e_00,gamma,fc_sum)
               end do
               !$OMP END PARALLEL DO
               deallocate(v2,v2_pos)
            end do
            do ii=1,nexc
               success(ii)=fc_sum(ii)*100>=fc_sum_min
               if(.not.success(ii))thr(ii)=AdjustThreshold(thr(ii))
               if(.not.w_ps)exit
            end do
            its=its+1
            if(.not.all(success).and. its<=fc_sum_min_its .and. all(abs(fc_sum-fc_sum_tot)>1d-6))then
               write(output_unit,*)'Iter ',its
#ifdef DEBUG
               if(check)then
                  fc_arr_check=0
                  fc_arr_check(1,1)=1
               end if
#endif
            end if
         end do
         if(its>fc_sum_min_its)write(Output_unit,*)'Reached maximum iterations: ',its-1
      end if
      
      
      deallocate(fc_arr)
#ifdef DEBUG
      if(check)then
         if(sum(fc_arr_check(FC_ss1%fc_arr_size,:))/=fc_ss2%fc_arr_size)then
            i = FC_ss1%fc_arr_size
            do ii = 1,fc_ss2%fc_arr_size
               if(fc_arr_check(i,ii)/=1)then
                  i=FC_ss1%fc_arr_size
               end if
            end do
         end if
         deallocate(fc_arr_check)
      end if
#endif
      
      do i = 1,mc_ms_tmp
         do j = 1,td%nq
            if(ignore_imag_modes)then
               if(td%TMExpand=='EE' .and. td%im_modes_ex_ex(j))exc_m_v(j,i,:)=0
               if(td%TMExpand=='GG' .and. td%im_modes_ex_gr(j))exc_m_v(j,i,:)=0
               if((td%TMExpand=='GE' .or. td%TMExpand=='EG') .and. &
               (td%im_modes_ex_gr(j) .or. td%im_modes_ex_ex(j)))exc_m_v(j,i,:)=0
            end if
         end do
      end do
      
      
      !write(output_unit,*)'11111'
      do i = 1,nexc
         !!$OMP CRITICAL
         if(w_ps)write(output_unit,'(F7.1," nm")')10d0**7/(wexc(i)*au_2_cm)
         do c_ms=1,mc_ms_tmp
            write(output_unit,'("Cl. ",1X,I2,I5," sig. modes; ",I4," modes reached max_v2s_ps (",I3,")")')c_ms,count(exc_m_v(:,c_ms,i)>0),count(exc_m_v(:,c_ms,i)==max_v2s(c_ms)),max_v2s(c_ms)
         end do
         call RepairExcM(exc_m_v(:,:,i),fc_ss2%nq,fc_ss2%mc_v,ii,max_v2s)
         write(output_unit,'(A,I0,A)')'Repaired Exc_m array ',ii,' times'
         write(output_unit,'(A,I4,A)')'Found: ',count(exc_m_v(:,1,i)/=0),' significant modes'
         write(output_unit,'("FC Sum: ",G10.4," %")')fc_sum(i)*100
         write(output_unit,'("Total prescreening FC Sum: ",G10.4," %")')fc_sum_tot*100
         write(output_unit,*)
         if(.not.w_ps)exit
         !!$OMP END CRITICAL
      end do
      
      do c_ms = 1,mc_ms_tmp
         deallocate(v2_arr(c_ms)%arr)
         deallocate(modes_arr_v2(c_ms)%arr)
      end do
      deallocate(v2_arr,modes_arr_v2,combs_n)
      
      !maybe pass this out instead?
      ! do i = 1,fc_ss1%fc_arr_size
         ! deallocate(v1_modes(i)%v,v1_modes(i)%v_pos)
      ! end do
      deallocate(exc_m2)
      if(v1_class>0)then
         deallocate(v1_pos_new)
         deallocate(exc_m1)
         call FC_storSys_dispose(fc_ss1)
      end if
      call FC_storSys_dispose(fc_ss2)
      time2=getTime()
      write(output_unit,'(3X,A)')TR(GetTimeStr(time2-time1))
      flush(output_unit)
      
   end function DEM_1v_RROA
   
   function AdjustThreshold(thr)result(res)
      double precision thr,res,ten
      integer power
      
      power=CEILING(log10(thr))-1
      
      res=thr-10d0**power !I dont know how to adjust it if I want a different decrement other than 10**-power
      !I guess add an if(res<0)
   end function AdjustThreshold
   
   subroutine RepairExcM(exc_m,nq,mc_v,repairCounter,max_vs)
      integer :: nq,mc_v,exc_m(nq,mc_v),max_vs(mc_v)
      integer :: exc_m_1mode(mc_v),i,j,jj,repairCounter
      
      repairCounter=0
      if(mc_v==1)return
      do i = 1,nq
         exc_m_1mode=exc_m(i,:)
         
         do j = 2,mc_v !some class
            do jj = 1,j-1 !lower classes than j
               !if mode in lower class is not excited high enough, excite it
               !this is done due to the recursive calculation scheme
               if(exc_m_1mode(jj)<exc_m_1mode(j))then
                  exc_m_1mode(jj)=exc_m_1mode(j)
                  repairCounter=repairCounter+1
               end if
               
               !If a mode in higher class reaches the prescreening maximum, get the excitation from a lower class
               if(exc_m_1mode(j)<exc_m_1mode(jj) .and. exc_m_1mode(j)==max_vs(j) .and. j-jj==1)then
                 exc_m_1mode(j)=exc_m_1mode(jj)
                 repairCounter=repairCounter+1
               end if
            end do
         end do
         exc_m(i,:)=exc_m_1mode
      end do
   end subroutine RepairExcM
   
      
   
   recursive subroutine Excm_v1_FC_inner_detExcM_RROA(i,v1,v1_pos,v1_clas,v2_arr,modes_arr_v2,combs_n,mc_v2,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,exc_m,thr,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
      integer mc_v2,c_v2
      integer :: i,nexc
      integer ii,j,v1_clas,v1(v1_clas)
      integer,allocatable :: v2(:)
      integer(int16) :: v1_pos(v1_clas)
      integer(int16),allocatable :: v2_pos(:)
      integer(int64) :: combs_n(mc_v2),ii_64
      type(v_col) v2_arr(mc_v2) !just some memory to pass, so it does not need to be allocated (slow-down)
      type(modes_arr_t) modes_arr_v2(:)
      type(FC_storSys) fc_ss2,fc_ss1
      type(FC_storSys_dusch) fc_ss_dusch
      integer :: v1_idx,ic,ic_thr
      integer exc_m(fc_ss2%red_n,fc_ss2%mc_v,nexc)
      integer :: v1_help(2),v2_pos_tmp(1),v2_tmp(1)
      
      double precision,allocatable,intent(inout) :: fc_arr(:,:)
      double precision thr(nexc),res,res2,wn,wexc(nexc),e_00,gamma,fc_sum(nexc)
      integer order
      integer(1),allocatable,intent(inout) :: fc_arr_check(:,:)
      logical check,eqna,w_ps
      
      !<v| is set up
      if(i>v1_clas)then !Set final state mode excitations
         v1_idx=ReduceFC_idx_fcss(v1,v1_pos,v1_clas,fc_ss1)
         res=fc_arr(v1_idx,1)
         fc_sum=fc_sum+res**2
         ! if(v1_clas>=1)then
            ! v2_pos_tmp(1)=0
            ! v2_tmp(1)=0
            ! res = FCOV1_stored_fcss_backtrack(v1,v1_pos,v1_clas,1,fc_ss1,v2_tmp,v2_pos_tmp,1,0,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,FC_ss1%fc_arr_size,FC_ss2%fc_arr_size) !<v|0>
         ! end if
#ifdef DEBUG
         if(check)fc_arr_check(v1_idx,1)=fc_arr_check(v1_idx,1)+1
         if(check.and.v1_clas<=2)then !mb use #ifdef at some point TODO
            select case(v1_clas)
               case(1)
                  eqna=.true.
                  if(v1(1)==13)then
                     res2=FC_v1_0_v2_13(fc_arr(1,1),fc_ss_dusch%C(v1_pos(1),v1_pos(1)),fc_ss_dusch%D(v1_pos(1)))
                     eqna=.false.
                  else if(v1(1)==20)then
                     res2=FC_v1_0_v2_20(fc_arr(1,1),fc_ss_dusch%C(v1_pos(1),v1_pos(1)),fc_ss_dusch%D(v1_pos(1)))
                     eqna=.false.
                  end if
                  !if(eqna)write(output_unit)'V1(1)=',v1(1),' ANALYTICAL EQ. NOT FOUND'
               case(2)
                  eqna=.true.
            end select
            if(.not.eqna)then !Analytical relation for NQi=13 is wrong
               if(res==0d0 .and. res2==0d0)then
               
               else if(abs(res/res2-1)>=1d-8)then
                  !$OMP CRITICAL
                  write(output_unit,*)'ANALYTICAL RELATION NOT MATCHING RECURSIVE, V1'
                  print *,v1
                  !call exit(42)
                  !$OMP END CRITICAL
               end if
            end if
         end if
#endif
         !set up |v>
         do c_v2 = 1,mc_v2
            allocate(v2_pos(c_v2),v2(c_v2))
            !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
            !$OMP PRIVATE(ic,ii_64,j) &
            !$OMP FIRSTPRIVATE(v2,v2_pos) &
            !$OMP SHARED(combs_n,c_v2,v1_pos,v1_clas,v1_idx,modes_arr_v2,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,thr,mc_v2,ic_thr) &
            !$OMP SHARED(nexc,wexc,wn,w_ps,e_00,gamma) &
            !$OMP REDUCTION(+:fc_sum) &
            !$OMP REDUCTION(MAX:exc_m)
            do ii_64=1,combs_n(c_v2)
               ! if(.not.w_ps)then
                  ! fc_sum(2:nexc)=0d0
                  ! exc_m(:,:,2:nexc)=0d0
               ! end if
               j=1
               ic=0
               v2_pos=modes_arr_v2(c_v2)%arr(:,ii_64)
               call ExcM_v2_FC_inner_detExcM_RROA(j,v1_idx,v2,v2_pos,c_v2,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,exc_m,thr,ic,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
            end do
            !$OMP END PARALLEL DO
            deallocate(v2_pos,v2)
         end do
         return
      end if
      
      ii=v1_clas-i+1
      do j = 1,fc_ss1%red_dims(v1_pos(ii),v1_clas)
         v1(ii)=j
         call Excm_v1_FC_inner_detExcM_RROA(i+1,v1,v1_pos,v1_clas,v2_arr,modes_arr_v2,combs_n,mc_v2,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,exc_m,thr,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
      end do
   end subroutine Excm_v1_FC_inner_detExcM_RROA
   
   recursive subroutine ExcM_v2_FC_0_inner_detExcM_RROA(i,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,fcs1,fcs2,check,exc_m,thr,ic,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
      integer,value,intent(in) :: i
      integer, intent(in) :: v2_clas,ic_thr,nexc
      integer(int16),intent(in) :: v2_pos(v2_clas)
      integer(int64) :: fcs1,fcs2
      integer, intent(inout) :: v2(v2_clas),ic
      type(FC_storSys),intent(in) :: fc_ss2
      type(FC_storSys_dusch),intent(in) :: fc_ss_dusch
      integer j
      integer ii,i2,idx_v2
      
      double precision,intent(inout) :: fc_arr(fcs1,fcs2)
      double precision,intent(in) :: thr(nexc),wn,wexc(nexc)
      double precision :: res,res2,den,wj,e_00,gamma,fc_sum(nexc),facr,faci
      integer,intent(inout) :: exc_m(fc_ss2%nq,fc_ss2%mc_v,nexc)
      
      integer(1),allocatable,intent(inout) :: fc_arr_check(:,:)
      logical,intent(in) :: check
      logical w_ps,eqna,success
      
      if(i>v2_clas)then
         idx_v2=ReduceFC_idx_fcss(v2,v2_pos,v2_clas,fc_ss2)
         res=fc_arr(1,idx_v2)
#ifdef DEBUG
         if(check)fc_arr_check(1,idx_v2)=fc_arr_check(1,idx_v2)+1
         if(check.and.v2_clas<=2)then !mb use #ifdef at some point TODO
            select case(v2_clas)
               case(1)
                  eqna=.true.
                  if(v2(1)==13)then
                     res2=FC_v1_0_v2_13(fc_arr(1,1),fc_ss_dusch%C(v2_pos(1),v2_pos(1)),fc_ss_dusch%D(v2_pos(1)))
                     eqna=.false.
                  else if(v2(1)==20)then
                     res2=FC_v1_0_v2_20(fc_arr(1,1),fc_ss_dusch%C(v2_pos(1),v2_pos(1)),fc_ss_dusch%D(v2_pos(1)))
                     eqna=.false.
                  end if
                  !if(eqna)write(output_unit)'V2(1)=',v2(1),' ANALYTICAL EQ. NOT FOUND'
               case(2)
                  eqna=.true.
                  !if(eqna)write(output_unit,*)'V2(1)=',v2(1),'V2(2)=',v2(2),' ANALYTICAL EQ. NOT FOUND'
            end select
            if(v2_clas==2)then
               if(v2(1)*v2(2)==66 .and. (v2(1)==6 .or. v2(1)==11))eqna=.true.
               if(v2(1)*v2(2)==72 .and. (v2(1)==18 .or. v2(1)==4))eqna=.true.
            end if
            
            if(.not.eqna)then !Analytical relation for NQi=13 is wrong (20 also I think)
               !order=floor(log10(min(abs(res),abs(res2))))
               !res=res*10**(-order)
               !res2=res2*10**(-order)
               if(res2==0d0 .and. res==0d0)then !edge case when K is zero: every odd excitation is gives a zero integral
               
               else if(abs(res/res2-1)>=1d-6)then
                  !$OMP CRITICAL
                  write(output_unit,*)'ANALYTICAL RELATION NOT MATCHING RECURSIVE, V2'
                  print *,v2
                  !call exit(42)
                  !$OMP END CRITICAL
               end if
               !res=res*10**(order)
            end if
         end if
#endif
         if(.not.w_ps)then
            if(abs(res)>=thr(1))then
               do j = 1,v2_clas
                  exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,1)=max(exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,1),v2(j))
               end do
               fc_sum(1)=fc_sum(1)+res**2
            else
               ic=ic+1
            end if
         else
            wj=e_00
            do ii=1,v2_clas
               wj=wj+fc_ss2%w(v2_pos(ii))*v2(ii)
            end do
            success=.false.
            do ii = 1,nexc
               ! den=1d0/sqrt(abs(wj-wn-wexc(ii)))
               den=(wj-wn-wexc(ii))**2 + gamma**2
               facr=sqrt(abs(wj-wn-wexc(ii))/den)
               faci=sqrt(gamma/den)
               if(abs(res*facr)>=thr(ii) .or. abs(res*faci)>=thr(ii))then
                  fc_sum(ii)=fc_sum(ii)+res**2
                  success=.true.
                  do j = 1,v2_clas
                     exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,ii)=max(exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,ii),v2(j))
                  end do
               end if
            end do
            if(.not.success)then
               ic=ic+1
            end if
         end if
         return
      end if
      
      ii=v2_clas-i+1
      do j = 1,fc_ss2%red_dims(v2_pos(ii),v2_clas)
         v2(ii)=j
         call ExcM_v2_FC_0_inner_detExcM_RROA(i+1,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,fcs1,fcs2,check,exc_m,thr,ic,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
         !if(ic>ic_thr)return
      end do
   end subroutine ExcM_v2_FC_0_inner_detExcM_RROA
   
   
   recursive subroutine ExcM_v2_FC_0_inner(i,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
      integer i,j,v2_clas,v2(v2_clas)
      integer(int16) v2_pos(v2_clas)
      type(FC_storSys) fc_ss2
      type(FC_storSys_dusch) fc_ss_dusch
      integer ii
      
      double precision fc_arr(:,:),res
      integer,allocatable :: fc_arr_checK(:,:)
      logical check
      
      if(i>v2_clas)then
         res=FCOV2_stored_fcss(v2,v2_pos,v2_clas,1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
         return
      end if
      
      ii=v2_clas-i+1
      do j = 1,fc_ss2%red_dims(v2_pos(ii),v2_clas)
         v2(ii)=j
         call ExcM_v2_FC_0_inner(i+1,v2,v2_pos,v2_clas,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
      end do
   end subroutine ExcM_v2_FC_0_inner
   
   recursive subroutine ExcM_v2_FC_inner_detExcM_RROA(i,v1_idx,v2,v2_pos,v2_clas,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,exc_m,thr,ic,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
      integer i,j,v2_clas,v2(v2_clas)
      integer(int16) v2_pos(v2_clas)
      type(FC_storSys) fc_ss1,fc_ss2
      type(FC_storSys_dusch) fc_ss_dusch
      integer ii,v1_idx,ic,ic_thr,nexc,v2_idx
      
      double precision,allocatable,intent(inout) :: fc_arr(:,:)
      integer(1),allocatable,intent(inout) :: fc_arr_check(:,:)
      double precision res,thr(nexc),den,wj,wn,wexc(nexc),e_00,gamma,fc_sum(nexc),facr,faci
      integer exc_m(fc_ss2%nq,fc_ss2%mc_v,nexc)
      logical check,w_ps,success
      
      if(i>v2_clas)then
         v2_idx=ReduceFC_idx_fcss(v2,v2_pos,v2_clas,fc_ss2)
         res=fc_arr(v1_idx,v2_idx)
#ifdef DEBUG
         if(check)fc_arr_check(v1_idx,v2_idx)=fc_arr_check(v1_idx,v2_idx)+1
#endif
         ! res = FCOV1_stored_fcss_backtrack(v1,v1_pos,v1_clas,1,fc_ss1,v2,v2_pos,v2_clas,1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,FC_ss1%fc_arr_size,fc_ss2%fc_arr_size)
         if(v1_idx==FC_ss1%fc_arr_size)then
            if(.not.w_ps)then
               if(abs(res)>=thr(1))then
                  do j = 1,v2_clas
                     exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,1)=max(exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,1),v2(j))
                  end do
                  fc_sum(1)=fc_sum(1)+res**2
               else
                  ic=ic+1
               end if
            else
               wj=e_00
               do ii=1,v2_clas
                  wj=wj+fc_ss2%w(v2_pos(ii))*v2(ii)
               end do
               success=.false.
               do ii = 1,nexc
                  ! den=1d0/sqrt(abs(wj-wn-wexc(ii)))
                  den=(wj-wn-wexc(ii))**2 + gamma**2
                  facr=sqrt(abs(wj-wn-wexc(ii))/den)
                  faci=sqrt(gamma/den)
                  if(abs(res*facr)>=thr(ii).or.abs(res*faci)>=thr(ii))then
                     fc_sum(ii)=fc_sum(ii)+res**2
                     success=.true.
                     do j = 1,v2_clas
                        exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,ii)=max(exc_m(fc_ss2%v_map(v2_pos(j)),v2_clas,ii),v2(j))
                     end do
                  end if
               end do
               if(.not.success)THEN
                  ic=ic+1
               end if
            end if
         end if
         
         ! if(v1_idx == fc_ss1%fc_arr_size .and. abs(res)>=thr)then
            ! do j = 1,v2_clas
               ! exc_m(v2_pos(j),v2_clas)=max(exc_m(v2_pos(j),v2_clas),v2(j))
            ! end do
         ! else
            ! ic=ic+1
         ! end if
         return
      end if
      
      ii=v2_clas-i+1
      do j = 1,fc_ss2%red_dims(v2_pos(ii),v2_clas)
         v2(ii)=j
         call ExcM_v2_FC_inner_detExcM_RROA(i+1,v1_idx,v2,v2_pos,v2_clas,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check,exc_m,thr,ic,ic_thr,w_ps,wn,wexc,nexc,e_00,gamma,fc_sum)
         !if(ic>ic_thr)return
      end do
   end subroutine ExcM_v2_FC_inner_detExcM_RROA
   
   recursive subroutine ExcM_v2_FC_inner(i,v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
      integer i,j,v2_clas,v2(v2_clas),v1_clas,v1(v1_clas)
      integer(int16) v1_pos(v1_clas),v2_pos(v2_clas)
      type(FC_storSys) fc_ss1,fc_ss2
      type(FC_storSys_dusch) fc_ss_dusch
      integer ii
      
      double precision fc_arr(:,:),res
      integer,allocatable :: fc_arr_checK(:,:)
      logical check
      
      if(i>v2_clas)then
         res = FCOV1_stored_fcss(v1,v1_pos,v1_clas,1,fc_ss1,v2,v2_pos,v2_clas,1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
         return
      end if
      
      ii=v2_clas-i+1
      do j = 1,fc_ss2%red_dims(v2_pos(ii),v2_clas)
         v2(ii)=j
         call ExcM_v2_FC_inner(i+1,v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,fc_ss1,fc_ss2,fc_ss_dusch,fc_arr,fc_arr_check,check)
      end do
   end subroutine ExcM_v2_FC_inner
         
   
   subroutine MakeTMs_FC(fc,td,expandInExc,rightIsExc,dip,mag,quad)
      double precision,intent(in) :: fc
      double precision,intent(out) :: dip(3),mag(3),quad(3,3)
      logical,intent(in) :: rightIsExc
      logical,intent(in) :: expandInExc
      type(ExcState), intent(in) :: td
      
      if(rightIsExc)then !TM from excited to ground, first TDM in transition polarizability first term
         if(expandInExc)then
            dip=td%u_ex*fc
            mag=-td%m_ex*fc
            ! mag=td%m_ex*fc
            quad=td%q_ex*fc
         else !expand in ground
            dip=td%u_gr_tr*fc
            mag=-td%m_gr_tr*fc
            ! mag=td%m_gr_tr*fc
            quad=td%q_gr_tr*fc
         end if
      else !TM from ground to excited, second TDM in transition polarizability first term
         if(expandInExc)then
            dip=td%u_ex_tr*fc
            mag=td%m_ex_tr*fc
            ! mag=-td%m_ex_tr*fc
            quad=td%q_ex_tr*fc
         else !expand in ground
            dip=td%u_gr*fc
            mag=td%m_gr*fc
            ! mag=-td%m_gr*fc
            quad=td%q_gr*fc
         end if
      end if
   end subroutine MakeTMs_FC
   
   !Make Transitions MomentS for Herzberg-Teller term
   subroutine MakeTMs_HT(k,ht,td,expandInExc,rightIsExc,dip,mag,quad)
      integer,intent(in) :: k
      double precision,intent(in) :: ht
      double precision,intent(inout) :: dip(3),mag(3),quad(3,3)
      logical,intent(in) :: rightIsExc
      logical,intent(in) :: expandInExc
      type(ExcState), intent(in) :: td
      
      if(rightIsExc)then !TM from excited to ground, first TDM in transition polarizability first term
         if(expandInExc)then
            dip=dip+td%du_ex(:,k)*ht
            mag=mag-td%dm_ex(:,k)*ht
            ! mag=mag+td%dm_ex(:,k)*ht
            quad=quad+td%dq_ex(:,:,k)*ht
         else !expand in ground
            dip=dip+td%du_gr_tr(:,k)*ht
            mag=mag-td%dm_gr_tr(:,k)*ht
            ! mag=mag+td%dm_gr_tr(:,k)*ht
            quad=quad+td%dq_gr_tr(:,:,k)*ht
         end if
      else !TM from ground to excited, second TDM in transition polarizability first term
         if(expandInExc)then
            dip=dip+td%du_ex_tr(:,k)*ht
            mag=mag+td%dm_ex_tr(:,k)*ht
            quad=quad+td%dq_ex_tr(:,:,k)*ht
         else !expand in ground
            dip=dip+td%du_gr(:,k)*ht
            mag=mag+td%dm_gr(:,k)*ht
            quad=quad+td%dq_gr(:,:,k)*ht
         end if
      end if
   end subroutine MakeTMs_HT
   
   
   recursive subroutine ExcM_v2_Polar_wps(i,v1,v1_pos,v1_clas,v1_idx,wm,v2,v2_pos,v2_clas,v3,v3_pos,v3_clas,v3_idx,wn,&
      canBeLow_gr,canBeLow,td,fc_ss1,fc_ss2,fc_arr,curpolar,curpolar_anti,gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments, &
      output_polContributions,alphafc,uncoup_modes_vmap,uncoup_modes_reddims,uncoup_c,uncoup_fcarr,fcsum1_un,fcsum3_un,contr)
      integer i,j,jj,k,l,ii,v,v_pos,kk
      integer v3_clas,v3(v3_clas)
      integer v1_clas,v1(v1_clas)
      integer v2_clas,v2(v2_clas),v2_idx,v2_add(v2_clas+1)
      integer(int64) :: v2_idx_low,v3_idx_low
      integer(int64) :: v3_idx,v1_idx
      integer(int16) v1_pos(v1_clas),v2_pos(v2_clas),v3_pos(v3_clas),v2_pos_add(v2_clas+1),k_pos
      integer nexc,idx_v2_local,iexc
      type(Polar_exc) :: curpolar,curpolar_anti
      type(ExcState) td
      type(fc_storSys) fc_ss1,fc_ss2
      double precision wm,wn,gamma,wexc(nexc),fc_arr(fc_ss1%fc_arr_size,fc_ss2%fc_arr_size)
      logical :: output_moments,d2=.false.,output_polContributions,w_ps,contr(9)
      integer uncoup_c,uncoup_modes_vmap(:),uncoup_modes_reddims(:)
      double precision uncoup_fcarr(:),Ci,Di,fcsum1_un,fcsum3_un
      
      double precision fc_12,fc_23,fc_12_low,fc_23_low,fc_12_high,fc_23_high
      double precision fc_un,fc_12_un,fc_23_un,fc_12_low_un,fc_23_low_un,fc_12_high_un,fc_23_high_un
      double precision sqrt_w,sqrt_wg,sqrt_vkp1,sqrt_vk,sqrt_vkp,sqrt_vk_gr,sqrt_vkp_gr,HT_12,HT_23
      
      logical neg,st,ht,only_ht,canBeLow(td%nq),canBeLow_gr(td%nq)
      logical canBeLow_2(td%nq),canbeLow_2_gr(td%nq),alphafc
      integer(int64) idx_is,idx_v2,idx_v2_un,idx_fs,fc_idx
      integer a,b,c,idx,idx_gr
      double precision dip_12(3),quad_12(3,3)
      double precision mag_12(3)
      double precision dip_23(3),quad_23(3,3)
      double precision mag_23(3)
      
      double precision dip_12_ht(3),quad_12_ht(3,3)
      double precision mag_12_ht(3)
      double precision dip_23_ht(3),quad_23_ht(3,3)
      double precision mag_23_ht(3)
      
      double precision dip_12_fc(3),quad_12_fc(3,3)
      double precision mag_12_fc(3)
      double precision dip_23_fc(3),quad_23_fc(3,3)
      double precision mag_23_fc(3)
      
      double precision wj,wj_inc,wj_un
      
      !v2 is set up (m => j transition)
      if(i>v2_clas)then
         idx_v2=ReduceFC_idx_fcss(v2,v2_pos,v2_clas,fc_ss2)
         fc_12=FC_arr(v1_idx,idx_v2)
         fc_23=FC_arr(v3_idx,idx_v2)
         wj=getFreq_short(v2,v2_pos,fc_ss2%w,v2_clas)+td%e_00
         
         if(only_ht)then
             dip_12_fc=0
             mag_12_fc=0
            quad_12_fc=0
             dip_23_fc=0
             mag_23_fc=0
            quad_23_fc=0
         else
            call MakeTMs_FC(fc_12,td,td%ExpandInExc1,.true.,dip_12_fc,mag_12_fc,quad_12_fc)
            call MakeTMs_FC(fc_23,td,td%ExpandInExc2,.false.,dip_23_fc,mag_23_fc,quad_23_fc)
         end if
          dip_12_ht=0
          mag_12_ht=0
         quad_12_ht=0
          dip_23_ht=0
          mag_23_ht=0
         quad_23_ht=0
         if(ht)then
            idx=1
            idx_gr=1
            do k = 1,td%nq
               if(canBeLow_gr(k))then
                  call LowerMode_inplace(v3,idx_gr,1,v3_clas,neg)
                  v3_idx_low=ReduceFC_idx_fcss(v3,v3_pos,v3_clas,fc_ss1)
                  fc_23_low=fc_arr(v3_idx_low,idx_v2)
                  ! if(ht2 .and. v3(idx_gr)>=1)then
                     ! call LowerMode_inplace(v3,idx_gr,1,v3_clas,neg)
                     ! v3_idx_low=ReduceFC_idx_fcss(v3,v3_pos,v3_clas,fc_ss1)
                     ! fc_23_low_2=fc_arr(v3_idx_low,idx_v2)
                     ! call IncrementMode_inplace(v3,idx_gr,1,v3_clas)
                  ! end if
                  call IncrementMode_inplace(v3,idx_gr,1,v3_clas)
                  fc_23_high=FCOV1_stored_HT_high_generic(v3,v3_pos,v3_clas,v2,v2_pos,v2_clas,k,idx_gr,td,fc_ss1,fc_ss2,v3_idx,idx_v2,fc_arr)
                  sqrt_vk_gr=sqrt_arr(v3(idx_gr)) !check this
                  sqrt_vkp_gr=sqrt_arr(v3(idx_gr)+1)
                  idx_gr=idx_gr+1
               else
                  fc_23_low=0
                  fc_23_high=FCOV1_stored_HT_high_generic(v3,v3_pos,v3_clas,v2,v2_pos,v2_clas,k,0,td,fc_ss1,fc_ss2,v3_idx,idx_v2,fc_arr)
                  sqrt_vk_gr=0
                  sqrt_vkp_gr=1
               end if
               
               if(canBeLow(k))then
                  call LowerMode_inplace(v2,idx,1,v2_clas,neg)
                  v2_idx_low=ReduceFC_idx_fcss(v2,v2_pos,v2_clas,fc_ss2)
                  fc_12_low=fc_arr(v1_idx,v2_idx_low)
                  call IncrementMode_inplace(v2,idx,1,v2_clas)
                  fc_12_high=FCOV2_stored_HT_high_generic(v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,k,idx,td,fc_ss1,fc_ss2,v1_idx,idx_v2,fc_arr)
                  sqrt_vk=sqrt_arr(v2(idx)) !check this
                  sqrt_vkp=sqrt_arr(v2(idx)+1)
                  idx=idx+1
               else
                  fc_12_low=0
                  fc_12_high=FCOV2_stored_HT_high_generic(v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,k,0,td,fc_ss1,fc_ss2,v1_idx,idx_v2,fc_arr)
                  sqrt_vk=0
                  sqrt_vkp=1
               end if
               
               ht_12=(sqrt_vk*fc_12_low+sqrt_vkp*fc_12_high)
               ht_23=(sqrt_vk_gr*fc_23_low+sqrt_vkp_gr*fc_23_high)
               call MakeTMs_HT(k,ht_12,td,td%ExpandInExc1,.true.,dip_12_ht,mag_12_ht,quad_12_ht)
               call MakeTMs_HT(k,ht_23,td,td%ExpandInExc2,.false.,dip_23_ht,mag_23_ht,quad_23_ht)
               
            end do
            
         end if         
         
         call DoPol_new_iexc(curpolar,dip_12_fc,dip_12_ht,mag_12_fc,mag_12_ht,quad_12_fc,quad_12_ht,&
               dip_23_fc,dip_23_ht,mag_23_fc,mag_23_ht,quad_23_fc,quad_23_ht,&
               wj,wn,wm,gamma,wexc,nexc,st,iexc, &
               alphafc,w_ps,output_polContributions,v2,v2_pos,v2_clas,contr)
         return
      end if
      
      !setup the v2 (|v2''>) mode
      do j = 1,fc_ss2%red_dims(v2_pos(v2_clas-i+1),v2_clas)
         v2(v2_clas-i+1)=j
         call Excm_v2_Polar_wps(i+1,v1,v1_pos,v1_clas,v1_idx,wm,v2,v2_pos,v2_clas,v3,v3_pos,v3_clas,v3_idx,wn,canBeLow_gr,canBeLow, &
         td,fc_ss1,fc_ss2,fc_arr,curpolar,curpolar_anti,gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments,output_polContributions,alphafc, &
         uncoup_modes_vmap,uncoup_modes_reddims,uncoup_c,uncoup_fcarr,fcsum1_un,fcsum3_un,contr)
      end do
   end subroutine Excm_v2_Polar_wps
   
   
   recursive subroutine ExcM_v2_Polar_wps_c1(i,v1,v1_pos,v1_clas,v1_idx,wm,v2,v2_pos,v3,v3_pos,v3_clas,v3_idx,wn, &
      canBeLow_gr,canBeLow,td,fc_ss1,fc_ss2,fc_arr,curpolar,curpolar_anti,gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments, &
      output_polContributions,alphafc,uncoup_modes_vmap,uncoup_modes_reddims,uncoup_c,uncoup_fcarr,fcsum1_un,fcsum3_un,contr)
      integer i,j,jj,k,ii,kk
      integer v3_clas,v3(v3_clas)
      integer v1_clas,v1(v1_clas)
      integer(int64) :: v3_idx,v1_idx
      integer,parameter :: v2_clas = 1
      integer v2(v2_clas),v2_idx,v2_add(v2_clas+1)
      integer(int16) v1_pos(v1_clas),v2_pos(v2_clas),v3_pos(v3_clas),v2_pos_add(v2_clas+1),v_pos
      integer nexc,iexc
      type(Polar_exc) :: curpolar,curpolar_anti
      type(ExcState) td
      type(fc_storSys) fc_ss1,fc_ss2
      double precision wm,wn,gamma,wexc(nexc),fc_arr(fc_ss1%fc_arr_size,fc_ss2%fc_arr_size)
      logical output_moments,output_polContributions,alphafc,w_ps,contr(9)
      integer uncoup_c,uncoup_modes_vmap(:),uncoup_modes_reddims(:)
      double precision uncoup_fcarr(:)
      double precision fc_12_low_arr(td%nq),fc_12_high_arr(td%nq)
      double precision fc_23_low_arr(td%nq),fc_23_high_arr(td%nq)
      double precision Ci,Di,fcsum1_un,fcsum3_un
      
      double precision fc_12,fc_23,fc_12_low,fc_23_low,fc_12_high,fc_23_high
      double precision fc_un,fc_12_un,fc_23_un,fc_12_low_un,fc_23_low_un,fc_12_high_un,fc_23_high_un
      double precision sqrt_w,sqrt_wg,sqrt_vkp1,sqrt_vk,sqrt_vk_gr,sqrt_vkp_gr,sqrt_vkp,HT_12,HT_23
      
      logical neg,st,ht,only_ht,canBeLow(td%nq),canBeLow_gr(td%nq)
      integer(int64) idx_is,idx_v2,idx_v2_un,idx_fs,fc_idx,v2_idx_low,v3_idx_low
      integer a,b,c,idx,k_counter,idx_gr,v
      double precision dip_12_fc(3),quad_12_fc(3,3)
      double precision mag_12_fc(3)
      double precision dip_23_fc(3),quad_23_fc(3,3)
      double precision mag_23_fc(3)
      
      double precision dip_12_ht(3),quad_12_ht(3,3)
      double precision mag_12_ht(3)
      double precision dip_23_ht(3),quad_23_ht(3,3)
      double precision mag_23_ht(3)
      
      double precision wj,wj_un,wj_inc
      
      !v2 is set up
      if(i>v2_clas)then
         idx_v2=ReduceFC_idx_fcss(v2,v2_pos,v2_clas,fc_ss2)
         fc_12=FC_arr(v1_idx,idx_v2)
         fc_23=FC_arr(v3_idx,idx_v2)
         wj=getFreq_short(v2,v2_pos,fc_ss2%w,v2_clas)+td%e_00
         
         if(only_ht)then
             dip_12_fc =0
             mag_12_fc =0
            quad_12_fc=0
             dip_23_fc =0
             mag_23_fc =0
            quad_23_fc=0
         else
            call MakeTMs_FC(fc_12,td,td%ExpandInExc1,.true.,dip_12_fc,mag_12_fc,quad_12_fc)
             ! dip_12=td%u*fc_12
             ! mag_12=-td%m*fc_12
            ! quad_12=td%q*fc_12
            call MakeTMs_FC(fc_23,td,td%ExpandInExc2,.false.,dip_23_fc,mag_23_fc,quad_23_fc)
             ! dip_23=td%u*fc_23
             ! mag_23=td%m*fc_23
            ! quad_23=td%q*fc_23
         end if
          dip_12_ht=0
          mag_12_ht=0
         quad_12_ht=0
          dip_23_ht=0
          mag_23_ht=0
         quad_23_ht=0
         if(ht)then
            idx=1
            idx_gr=1
            do k = 1,td%nq
               ! sqrt_w=td%sqrt_w(k)
               ! sqrt_wg=td%sqrt_wg(k)
               if(canBeLow_gr(k))then
                  call LowerMode_inplace(v3,idx_gr,1,v3_clas,neg)
                  v3_idx_low=ReduceFC_idx_fcss(v3,v3_pos,v3_clas,fc_ss1)
                  fc_23_low=fc_arr(v3_idx_low,idx_v2)
                  call IncrementMode_inplace(v3,idx_gr,1,v3_clas)
                  fc_23_high=FCOV1_stored_HT_high_generic(v3,v3_pos,v3_clas,v2,v2_pos,v2_clas,k,idx_gr,td,fc_ss1,fc_ss2,v3_idx,idx_v2,fc_arr)
                  sqrt_vk_gr=sqrt_arr(v3(idx_gr)) !check this
                  sqrt_vkp_gr=sqrt_arr(v3(idx_gr)+1)
                  idx_gr=idx_gr+1
               else
                  fc_23_low=0
                  fc_23_high=FCOV1_stored_HT_high_generic(v3,v3_pos,v3_clas,v2,v2_pos,v2_clas,k,0,td,fc_ss1,fc_ss2,v3_idx,idx_v2,fc_arr)
                  sqrt_vk_gr=0
                  sqrt_vkp_gr=1
               end if
               if(uncoup_c>0)then   
                  fc_23_low_arr(k)=fc_23_low
                  fc_23_high_arr(k)=fc_23_high
               end if
               
               if(canBeLow(k))then
                  call LowerMode_inplace(v2,idx,1,v2_clas,neg)
                  v2_idx_low=ReduceFC_idx_fcss(v2,v2_pos,v2_clas,fc_ss2)
                  fc_12_low=fc_arr(v1_idx,v2_idx_low)
                  call IncrementMode_inplace(v2,idx,1,v2_clas)
                  fc_12_high=FCOV2_stored_HT_high_generic(v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,k,idx,td,fc_ss1,fc_ss2,v1_idx,idx_v2,fc_arr)
                  sqrt_vk=sqrt_arr(v2(idx)) !check this
                  sqrt_vkp=sqrt_arr(v2(idx)+1)
                  idx=idx+1
               else
                  fc_12_low=0
                  fc_12_high=FCOV2_stored_HT_high_generic(v1,v1_pos,v1_clas,v2,v2_pos,v2_clas,k,0,td,fc_ss1,fc_ss2,v1_idx,idx_v2,fc_arr)
                  sqrt_vk=0
                  sqrt_vkp=1
               end if
               if(uncoup_c>0)then
                  fc_12_low_arr(k)=fc_12_low
                  fc_12_high_arr(k)=fc_12_high
               end if
               
               ht_12=(sqrt_vk*fc_12_low+sqrt_vkp*fc_12_high)
               ht_23=(sqrt_vk_gr*fc_23_low+sqrt_vkp_gr*fc_23_high)
               call MakeTMs_HT(k,ht_12,td,td%ExpandInExc1,.true.,dip_12_ht,mag_12_ht,quad_12_ht)
               call MakeTMs_HT(k,ht_23,td,td%ExpandInExc2,.false.,dip_23_ht,mag_23_ht,quad_23_ht)
               
            end do
         end if
         
         call DoPol_new_iexc(curpolar,dip_12_fc,dip_12_ht,mag_12_fc,mag_12_ht,quad_12_fc,quad_12_ht,&
               dip_23_fc,dip_23_ht,mag_23_fc,mag_23_ht,quad_23_fc,quad_23_ht,&
               wj,wn,wm,gamma,wexc,nexc,st,iexc, &
               alphafc,w_ps,output_polContributions,v2,v2_pos,v2_clas,contr)
               
         return
      end if
      
      !setup the v2 (|v2''>) mode
      do j = 1,fc_ss2%red_dims(v2_pos(v2_clas-i+1),v2_clas)
         v2(v2_clas-i+1)=j
         call Excm_v2_Polar_wps_c1(i+1,v1,v1_pos,v1_clas,v1_idx,wm,v2,v2_pos,v3,v3_pos,v3_clas,v3_idx,wn,canBeLow_gr,canBeLow, &
         td,fc_ss1,fc_ss2,fc_arr,curpolar,curpolar_anti,gamma,wexc,nexc,ht,only_ht,st,w_ps,iexc,output_moments,output_polContributions,alphafc, &
         uncoup_modes_vmap,uncoup_modes_reddims,uncoup_c,uncoup_fcarr,fcsum1_un,fcsum3_un,contr)
      end do
   end subroutine Excm_v2_Polar_wps_c1
   

   subroutine DoPol_new_iexc(pol,dip_12_fc,dip_12_ht,mag_12_fc,mag_12_ht,quad_12_fc,quad_12_ht,dip_23_fc,dip_23_ht,mag_23_fc,mag_23_ht,quad_23_fc,quad_23_ht,&
      wj,wn,wm,gamma,wexc,nexc,st,iexc,alphafc,w_ps,output_polContributions,v2,v2_pos,v2_n,contr)
      integer,intent(in) :: nexc,v2_n,v2(v2_n)
      integer(int16),intent(in) :: v2_pos(v2_n)
      logical,intent(in) :: st,alphafc,w_ps,output_polContributions,contr(9)
      double precision,intent(in) :: wexc(nexc),wj,wm,wn,gamma
      double precision,intent(in) :: dip_12_fc(3),dip_23_fc(3)
      double precision,intent(in) :: quad_12_fc(3,3),quad_23_fc(3,3)
      double precision, intent(in) :: mag_12_fc(3),mag_23_fc(3)
      double precision,intent(in) :: dip_12_ht(3),dip_23_ht(3)
      double precision,intent(in) :: quad_12_ht(3,3),quad_23_ht(3,3)
      double precision, intent(in) :: mag_12_ht(3),mag_23_ht(3)
      double precision :: dip_12(3),mag_12(3),quad_12(3,3),dip_23(3),mag_23(3),quad_23(3,3)
      
      type(Polar_exc),intent(inout) :: pol
      
      integer aa,bb,cc,iexc,i,ii
      double precision :: wjm,wjn,buf
      double complex :: den1,den2
      double complex :: Ap_inc(3,3),ap_fcht_inc(3,3),G_inc(3,3),Gc_inc(3,3),A_inc(3,3,3),Ac_inc(3,3,3)
      
      wjn=wj-wn
      
      dip_12=dip_12_fc+dip_12_ht
      mag_12=mag_12_fc+mag_12_ht
      quad_12=quad_12_fc+quad_12_ht
      
      dip_23=dip_23_fc+dip_23_ht
      mag_23=mag_23_fc+mag_23_ht
      quad_23=quad_23_fc+quad_23_ht
   


      if(.not.w_ps)then
         do i = 1,nexc
         den1=(wjn-wexc(i)-iu*gamma)
         
         Ap_inc=0d0
         ap_fcht_inc=0d0
         G_inc=0d0
         Gc_inc=0d0
         A_inc=0d0
         Ac_inc=0d0
         
         
         
         
         if(st)then
            !wjm=wj-wm
            den2=(wjn+(wexc(i)-wm)+iu*gamma)
            if(alphafc)then
               do bb=1,3
                  do aa=1,3
                     Ap_inc(aa,bb)=Ap_inc(aa,bb) + dip_12_fc(bb)*dip_23_fc(aa)/den2
                     Ap_fcht_inc(aa,bb)=Ap_fcht_inc(aa,bb) + dip_12(bb)*dip_23(aa)/den2
                     G_inc(aa,bb)=G_inc(aa,bb) + mag_12(bb)*dip_23(aa)*iu/den2
                     Gc_inc(aa,bb)=Gc_inc(aa,bb) + dip_12(bb)*mag_23(aa)*iu/den2
                     do cc=1,bb
                        A_inc(aa,bb,cc)=A_inc(aa,bb,cc) + quad_12(bb,cc)*dip_23(aa)/den2
                        Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc) + dip_12(aa)*quad_23(bb,cc)/den2
                     end do
                  end do
               end do
            else
               do bb=1,3
                  do aa=1,3
                     Ap_inc(aa,bb)=Ap_inc(aa,bb) + dip_12(bb)*dip_23(aa)/den2
                     G_inc(aa,bb)=G_inc(aa,bb) + mag_12(bb)*dip_23(aa)*iu/den2
                     Gc_inc(aa,bb)=Gc_inc(aa,bb) + dip_12(bb)*mag_23(aa)*iu/den2
                     do cc=1,bb
                        A_inc(aa,bb,cc)=A_inc(aa,bb,cc) + quad_12(bb,cc)*dip_23(aa)/den2
                        Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc) + dip_12(aa)*quad_23(bb,cc)/den2
                     end do
                  end do
               end do
            end if
         end if
         if(alphafc)then
            do bb=1,3
               do aa=1,3
                  Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12_fc(aa)*dip_23_fc(bb)/den1
                  Ap_fcht_inc(aa,bb)=Ap_fcht_inc(aa,bb)+dip_12(aa)*dip_23(bb)/den1
                  G_inc(aa,bb)=G_inc(aa,bb)+dip_12(aa)*mag_23(bb)*iu/den1
                  Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12(aa)*iu*dip_23(bb)/den1
                  do cc=1,bb
                     A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12(aa)*quad_23(bb,cc)/den1
                     Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12(bb,cc)*dip_23(aa)/den1
                  end do
               end do
            end do
         else
            do bb=1,3
               do aa=1,3
                  if(contr(1))then
                     Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12_fc(aa)*dip_23_fc(bb)/den1
                     G_inc(aa,bb)=G_inc(aa,bb)+dip_12_fc(aa)*mag_23_fc(bb)*iu/den1
                     Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12_fc(aa)*iu*dip_23_fc(bb)/den1
                  end if
                  if(contr(2))then
                     Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12_fc(aa)*dip_23_ht(bb)/den1
                     G_inc(aa,bb)=G_inc(aa,bb)+dip_12_fc(aa)*mag_23_ht(bb)*iu/den1
                     Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12_fc(aa)*iu*dip_23_ht(bb)/den1
                  end if
                  if(contr(3))then
                     Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12_ht(aa)*dip_23_fc(bb)/den1
                     G_inc(aa,bb)=G_inc(aa,bb)+dip_12_ht(aa)*mag_23_fc(bb)*iu/den1
                     Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12_ht(aa)*iu*dip_23_fc(bb)/den1
                  end if
                  if(contr(4))then
                     Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12_ht(aa)*dip_23_ht(bb)/den1
                     G_inc(aa,bb)=G_inc(aa,bb)+dip_12_ht(aa)*mag_23_ht(bb)*iu/den1
                     Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12_ht(aa)*iu*dip_23_ht(bb)/den1
                  end if
                  ! Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12(aa)*dip_23(bb)/den1
                  ! G_inc(aa,bb)=G_inc(aa,bb)+dip_12(aa)*mag_23(bb)*iu/den1
                  ! Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12(aa)*iu*dip_23(bb)/den1
                  do cc=1,bb
                     if(contr(1))then
                        A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12_fc(aa)*quad_23_fc(bb,cc)/den1
                        Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12_fc(bb,cc)*dip_23_fc(aa)/den1
                     end if
                     if(contr(2))then
                        A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12_fc(aa)*quad_23_ht(bb,cc)/den1
                        Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12_fc(bb,cc)*dip_23_ht(aa)/den1
                     end if
                     if(contr(3))then
                        A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12_ht(aa)*quad_23_fc(bb,cc)/den1
                        Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12_ht(bb,cc)*dip_23_fc(aa)/den1
                     end if
                     if(contr(4))then
                        A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12_ht(aa)*quad_23_ht(bb,cc)/den1
                        Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12_ht(bb,cc)*dip_23_ht(aa)/den1
                     end if
                     ! A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12(aa)*quad_23(bb,cc)/den1
                     ! Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12(bb,cc)*dip_23(aa)/den1
                  end do
               end do
            end do
         end if
         ! end if
         
         if(output_polContributions)then
            buf=RamanPolcon(Ap_inc)
            if(buf>1d-7)then
               !$OMP CRITICAL
               write(op_polcontrs_unit-1+i,'(1X,2(G11.4,1X))',advance='no')buf
               if(v2_n==0)then
                  write(op_polcontrs_unit-1+i,'(1X,I3,1X,I2)',advance='no')0,0
               else
                  do ii = 1,v2_n
                     write(op_polcontrs_unit-1+i,'(1X,I3,"^",I2)',advance='no')v2_pos(ii),v2(ii)
                  end do
               end if
               write(op_polcontrs_unit-1+i,'(1X,F11.4," cm-1")')wjn*au_2_cm
               !$OMP END CRITICAL
            end if
         end if
         
         if(alphaFC)then
            pol%Ap_fcht(:,:,i)=pol%Ap_fcht(:,:,i)+ap_fcht_inc
         end if
         pol%Ap(:,:,i)=pol%Ap(:,:,i)+Ap_inc
         pol%G(:,:,i)=pol%G(:,:,i)+G_inc
         pol%Gc(:,:,i)=pol%Gc(:,:,i)+Gc_inc
         pol%A(:,:,:,i)=pol%A(:,:,:,i)+A_inc
         pol%Ac(:,:,:,i)=pol%Ac(:,:,:,i)+Ac_inc
         end do
         return
      end if
      
      den1=(wjn-wexc(iexc)-iu*gamma)
      
      Ap_inc=0d0
      ap_fcht_inc=0d0
      G_inc=0d0
      Gc_inc=0d0
      A_inc=0d0
      Ac_inc=0d0
      
      if(st)then
         !wjm=wj-wm
         den2=(wjn+(wexc(iexc)-wm)+iu*gamma)
         if(alphafc)then
            do bb=1,3
               do aa=1,3
                  Ap_inc(aa,bb)=Ap_inc(aa,bb) + dip_12_fc(bb)*dip_23_fc(aa)/den2
                  Ap_fcht_inc(aa,bb)=Ap_fcht_inc(aa,bb) + dip_12(bb)*dip_23(aa)/den2
                  G_inc(aa,bb)=G_inc(aa,bb) + mag_12(bb)*dip_23(aa)*iu/den2
                  Gc_inc(aa,bb)=Gc_inc(aa,bb) + dip_12(bb)*mag_23(aa)*iu/den2
                  do cc=1,bb
                     A_inc(aa,bb,cc)=A_inc(aa,bb,cc) + quad_12(bb,cc)*dip_23(aa)/den2
                     Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc) + dip_12(aa)*quad_23(bb,cc)/den2
                  end do
               end do
            end do
         else
            do bb=1,3
               do aa=1,3
                  Ap_inc(aa,bb)=Ap_inc(aa,bb) + dip_12(bb)*dip_23(aa)/den2
                  G_inc(aa,bb)=G_inc(aa,bb) + mag_12(bb)*dip_23(aa)*iu/den2
                  Gc_inc(aa,bb)=Gc_inc(aa,bb) + dip_12(bb)*mag_23(aa)*iu/den2
                  do cc=1,bb
                     A_inc(aa,bb,cc)=A_inc(aa,bb,cc) + quad_12(bb,cc)*dip_23(aa)/den2
                     Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc) + dip_12(aa)*quad_23(bb,cc)/den2
                  end do
               end do
            end do
         end if
      end if
      if(alphafc)then
         do bb=1,3
            do aa=1,3
               Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12_fc(aa)*dip_23_fc(bb)/den1
               Ap_fcht_inc(aa,bb)=Ap_fcht_inc(aa,bb)+dip_12(aa)*dip_23(bb)/den1
               G_inc(aa,bb)=G_inc(aa,bb)+dip_12(aa)*mag_23(bb)*iu/den1
               Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12(aa)*iu*dip_23(bb)/den1
               do cc=1,bb
                  A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12(aa)*quad_23(bb,cc)/den1
                  Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12(bb,cc)*dip_23(aa)/den1
               end do
            end do
         end do
      else
         do bb=1,3
            do aa=1,3
               Ap_inc(aa,bb)=Ap_inc(aa,bb)+dip_12(aa)*dip_23(bb)/den1
               G_inc(aa,bb)=G_inc(aa,bb)+dip_12(aa)*mag_23(bb)*iu/den1
               Gc_inc(aa,bb)=Gc_inc(aa,bb)+mag_12(aa)*iu*dip_23(bb)/den1
               do cc=1,bb
                  A_inc(aa,bb,cc)=A_inc(aa,bb,cc)+dip_12(aa)*quad_23(bb,cc)/den1
                  Ac_inc(aa,bb,cc)=Ac_inc(aa,bb,cc)+quad_12(bb,cc)*dip_23(aa)/den1
               end do
            end do
         end do
      end if
      ! end if
      
      if(output_polContributions)then
         buf=RamanPolcon(Ap_inc)
         if(buf>1d-7)then
            !$OMP CRITICAL
            write(op_polcontrs_unit-1+iexc,'(1X,(G11.4,1X))',advance='no')buf
            if(v2_n==0)then
               write(op_polcontrs_unit-1+iexc,'(1X,I3,1X,I2)',advance='no')0,0
            else
               do ii = 1,v2_n
                  write(op_polcontrs_unit-1+iexc,'(1X,I3,"^",I2)',advance='no')v2_pos(ii),v2(ii)
               end do
            end if
            write(op_polcontrs_unit-1+iexc,'(1X,F11.4," cm-1")')wjn*au_2_cm
            !$OMP END CRITICAL
         end if
      end if
      
      if(alphaFC)then
         pol%Ap_fcht(:,:,iexc)=pol%Ap_fcht(:,:,iexc)+ap_fcht_inc
      end if
      pol%Ap(:,:,iexc)=pol%Ap(:,:,iexc)+Ap_inc
      pol%G(:,:,iexc)=pol%G(:,:,iexc)+G_inc
      pol%Gc(:,:,iexc)=pol%Gc(:,:,iexc)+Gc_inc
      pol%A(:,:,:,iexc)=pol%A(:,:,:,iexc)+A_inc
      pol%Ac(:,:,:,iexc)=pol%Ac(:,:,:,iexc)+Ac_inc
   end subroutine DoPol_new_iexc

   function Calc_a2(Ap)result(res)
      double complex Ap(3,3),buf1,buf2
      double precision res
      
      buf1=Ap(1,1)+Ap(2,2)+Ap(3,3)
      buf2=conjg(Ap(1,1))+conjg(Ap(2,2))+conjg(Ap(3,3))
      res=1d0/9d0*realpart(buf1*buf2)
      
   end function Calc_a2
   
   function Calc_Bsa2(Ap)result(res)
      integer a,b
      double complex Ap(3,3),Ap_s(3,3),buf
      double precision res
      
      Ap_s=Ts(Ap)
      
      buf=(0d0,0d0)
      do a = 1,3
         do b = 1,3
            buf=buf+3*Ap_s(a,b)*conjg(Ap_s(a,b))-Ap_s(a,a)*conjg(Ap_s(b,b))
         end do
      end do
      res=0.5d0*realpart(buf)
   end function Calc_Bsa2
   
   function Calc_Baa2(Ap)result(res)
      integer a,b
      double complex Ap(3,3),Ap_a(3,3),buf
      double precision res
      
      Ap_a=Ta(Ap)
      buf=(0d0,0d0)
      do a = 1,3
         do b = 1,3
            buf=buf+Ap_a(a,b)*conjg(Ap_a(a,b))
         end do
      end do
      res=1.5d0*realpart(buf)
   end function Calc_Baa2
   
   function RamanPolcon(ap)result(I)
      double complex ap(3,3)
      double precision a2,Bsa2,Baa2,I
      
      a2=Calc_a2(ap)
      Bsa2=Calc_bsa2(ap)
      Baa2=Calc_baa2(ap)
      
      I=4d0/90d0*(45*a2+7*Bsa2+5*Baa2)
   end function RamanPolcon
   
   function Ts(T)result(res)
      double complex T(3,3),res(3,3)
      integer a,b
      
      do a = 1,3
         do b = 1,3
            res(a,b)=0.5d0*(T(a,b)+T(b,a))
         end do
      end do
   end function ts
   
   function Ta(T)result(res)
      double complex T(3,3),res(3,3)
      integer a,b
      
      do a = 1,3
         do b = 1,3
            res(a,b)=0.5d0*(T(a,b)-T(b,a))
         end do
      end do
   end function ta
   
end module RROA


program fcov
   use constants
   use ANALYTICAL
   use Fcov_primitive
   use iso_fortran_env
   use FCOV_storage
   use FCOV_FCFuns
   use OMP_functions
   use strings
   use RROA
   use stuff
   use util
   !$ use omp_lib
   
   implicit none
   
   !----------------------
   !C functions
   !----------------------
   !cpu.c - copy-paste from NASA https://www.nas.nasa.gov/hecc/support/kb/instrumenting-your-fortran-code-to-check-processthread-placement_309.html
   integer, external :: findmycpu
   
   !----------------------
   
   integer,parameter :: maxExcStates = 100 !completely arbitrary number
   
   double precision,allocatable :: ens(:),wg(:)
   double precision,allocatable :: vv(:),dd(:),aa(:),qq(:)
   double precision,allocatable :: smat(:,:)
   type(ExcState),allocatable,target :: tds(:)
   type(ExcState),pointer :: td
   integer,allocatable :: doModes(:)
   
   
   double precision time1,time2
   double precision :: thresholdOrig_gr=0.1
   double precision :: thresholdOrig_ex=0.1
   integer :: clas,argc,jj,istat,nat,nat2,td_i,nroot
   integer :: nq_gr,sz
   integer(int64) :: comb_n,max_storage,i_64
   integer :: line
   integer :: mc_ms=3,mc_gs=0,mc_fs=1 !Max Class - Excited, Ground, Final State
   integer :: mc_gfs
   
   
   integer,allocatable :: max_v1s(:),max_v3s(:),max_v2s(:)
   integer,allocatable :: max_v2s_ps(:),min_v2s_ps(:)
   integer mc_ms_ps
   integer,allocatable :: min_v1s(:),min_v2s(:)
   integer,allocatable :: min_v2s_excm(:)

   
   integer :: v1_n,v3_n,td_counter
   
   integer :: idx,i,ii,j,k
   integer :: c_ms,c_gs
   character(2) TMExpand_def
   character(80) c80,c80_2
   character(400) :: outFiles(maxExcStates)
   integer :: nexc,summ
   double precision :: bufr,bigTime
   double precision :: min_freq_mode,max_freq_mode,eau,uncoup_lim
   integer :: uncoup_modes_maxval=-1
   
   double precision :: fc_sum_min_ps
   integer :: fc_sum_min_ps_its = 10
   
   double precision :: mult_tm(5),mult_dtm(5)
   logical :: isOMP = .false.,fex, isRROA=.true. !for now
   logical :: HT = .true.,ht2=.false.
   logical :: vel = .false.
   double precision :: kbt_lim = 0.7,fwhm=7,temp=300,add_pol_coeff=0d0
   logical :: rroa_st = .false., fcarrcheck=.false.
   logical :: isEVCD = .false., stokes=.true.
   logical :: verbose = .true.,full_ht=.false.
   logical :: useA=.true.,useg=.true.,only_ht=.false.
   logical :: write_ten=.true.,write_excm=.false.,write_fc_sys=.false.,write_inv=.false.,stdoutToFile=.false.
   logical :: write_fcarr=.false. !very false, this has the potential to waste GBs of HDD space
   logical :: use_gauss=.false.
   logical :: vertical_en=.false.
   logical :: w_ps=.true. , ignore_ps=.false.
   logical :: derivatives_in_excGeom=.false.
   logical :: evcd=.false.,cpl=.false.
   logical :: antistokes=.false.
   logical :: sel_rules=.false.,sel_rules2=.false.
   logical :: correct_gr_freqs = .false.
   logical :: output_moments = .false., output_polars=.false.
   logical :: ignore_imag_modes = .false.,output_polContrs=.false.
   logical :: isAH = .false.,MomentsExToGr=.false.
   logical :: AlphaFC=.false.,avg_e00=.false.,wexc_eq_e00=.false.
   logical :: spectrum_Temp=.true.,AHAS=.false.
   logical :: wr_is_e00=.false.,vel_grad=.true.
   logical :: coordinates=.true.
   logical :: transform_ground_tm=.false., transform_excited_tm=.false.
   logical :: write_corrf=.false.,write_fft_cf=.false.,td_trueground=.false.
   logical :: gr_exp = .false.,num_integ=.true.,norm_fft=.true.,elpol_exc=.false.
   logical :: correctPhaseX=.true.,correctPhaseX_abs=.false.,interpolateFFT=.false.
   logical :: wexc_adapt=.false.,SE_is_SG=.false.,elpol_grad=.true.,elpol_only=.false.
   logical :: w_ad_zero=.false.,td_alt=.false.,td_fixphase=.true.
   
   !In rroa_td_num, to get the same ROA sign as experiment and Cheeseman for Co(III)EDDS complex, I need to switch the signs on the magnetic dipoles. 
   !It does not make any sense, and I am not certain about it but whatever.
   logical :: switchMagneticSign = .false. 
   
   integer :: red_n2_max=huge(1),transform_tm_type=0
   integer :: red_n2_max_algo=0
   integer :: n_gr=1,fundLeadCount=0
   
   double precision :: overwrite_e00=0d0
   
   double precision :: thr_v2_0=0.1, thr_v2_comb=0.1, thr_v2_other=0.1, thr_v2_ov=0.1
   logical :: thr_v2_0_found=.false., thr_v2_comb_found=.false., thr_v2_other_found=.false., thr_v2_ov_found=.false.
   
   
   character(3) runType
   
   double precision,allocatable :: m
   double precision,allocatable :: G(:,:),GT(:,:),GP(:,:),GPT(:,:),JT(:,:),F(:,:),FI(:,:),TP(:,:)
   double precision,allocatable :: T(:,:),ee(:,:)
   integer ierr,ix
   
   !Threads
   integer :: n_thr
   
   !Molecule
   integer,allocatable :: z(:)
   double precision,allocatable :: mas(:)
   integer :: NQ = 0
   logical :: linearMol=.false.
   
   !Vibronic systems
   character(20) str_v1,str_v3
   double precision kfac_,kfac
   character(120) outfile_tmp
   character(80) fille
   character(120) dir_gr,dir_ex,dir_cwd
   character(20) sp_type
   logical contr(9)
   
   !Spectrum
   integer npx,wmin,wmax
   double precision freq_i,freq_j,kt,w
   
   !ECD
   double precision thr_d,thr_R
   integer fwhm_which
   
   !ROA
   double precision :: gamma=100,gamma_el=850,gamma_,theta=0,td_eps=0
   logical :: RGSF = .false.
   integer,parameter :: ns0=19,ni0=13
   integer iexc,td_n
   integer,allocatable :: rroa_tabOutputs(:,:),IArr(:)
   double precision,allocatable :: rroa_spr(:,:,:,:),rroa_inv(:)
   double precision,allocatable :: wexc(:),wexc_nm(:)
   character(6),allocatable :: wexc_nm_str(:)
   logical :: rroa_spr_do(ns0)
   character(9),parameter :: rroa_exp(ns0) = ['ICP_0    ','ICP_x_90 ','ICP_z_90 ','ICP_*_90 ','ICP_u_90 ','ICP_180  ',&
   'SCP_0    ','SCP_x_90 ','SCP_z_90 ','SCP_*_90 ','SCP_u_90 ','SCP_180  ','DCPI_0   ','DCPI_90  ','DCPI_180 ','DCPII_0  ','DCPII_90 ','DCPII_180','SPECIAL  ']
   
   integer :: storage_idx
   double precision :: popsum,minpop
   
   logical :: td_approach = .true.
   integer :: td_n_points = 12,td_batch_n=10,td_j_tol_int=3
   double precision :: td_tmax=10d0,td_fs=0,td_j_tol
   logical :: td_sparse=.true.
   
   type(Transition_R),allocatable :: Trs(:)
   
   contr=.true.
   td_j_tol=10d0**(-td_j_tol_int)
   call GETCWD(dir_cwd)
   bigTime=getTime()
   !program defaults
   n_thr=1
   !!$ n_thr=OMP_get_max_threads()
   isOMP=.false.
   !$ isOMP=.true.
   minpop = 0.5d0
   gamma=500*cm_2_au
   gamma_el=0
   theta=0
   write(sp_type,*)'RROA'
   mult_tm=[1d0,1d0,0.5d0,1d0,1d0]
   mult_dtm=[1d0,1d0,0.5d0,1d0,1d0]
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc<1)then
      write(output_unit,*)'Program for calculation of vibronic spectra.'
      write(output_unit,*)'U S A G E: '
      write(output_unit,*)'FCOV.out <ground freq> <excited_1 freq gr.geom> <excited_1 freq ex.geom> <excited_2 freq gr.geom> ...'
      call exit(1)
   end if
   td_n=(argc-1)/2
   
   call ReadOpt()
   inquire(file='FILE.TR',exist=fex)
   if(fex)then
      call ReadFileTR(Trs)
   end if
   
   if(stdoutToFile)open(output_unit,file='FCOV.out',action='WRITE')
   
#ifdef DEBUG
   write(output_unit,*)'Debug'
#endif   
   if(verbose)then
      write(output_unit,'(A)')'-------------------------------------------------'
      write(output_unit,'(A)')'|          Franck Condon OVerlap program        |'
      write(output_unit,'(A)')'-------------------------------------------------'
      write(output_unit,*)
      
      if(isOMP)then
         write(output_unit,*)'OpenMP is active'
         !$OMP PARALLEL !!NUM_THREADS(n_thr)
         !$ if(omp_get_thread_num()==0)then
         !$    n_thr=omp_get_num_threads()
         !$    write(output_unit,*)'Will use ',omp_get_num_threads(),' threads'
         !$    write(output_unit,*)'OMP_get_max_threads: ',omp_get_max_threads()
         !$    write(output_unit,*)'OMP_get_num_procs  : ',OMP_get_num_procs()
         !$ end if
         !$OMP BARRIER
         !$OMP CRITICAL
         !$    write(output_unit,*)'Thread ',omp_get_thread_num(),'is on CPU ',findmycpu()
         !$OMP END CRITICAL
         !$OMP END PARALLEL
      end if
   else
      if(isOMP)then
         write(output_unit,*)'OpenMP'
      end if
   end if
   write(output_unit,'(8X,5(1X,A7))')'LEN.DIP','VEL.DIP','MAG.DIP','Q.DIAG','Q.nDIAG'
   write(output_unit,'(A8,5(1X,F7.2))')'MULT_TM ',mult_tm
   write(output_unit,'(A8,5(1X,F7.2))')'MULT_DTM',mult_dtm
   call MakeSqrtArr(1000)
   do i = 1,argc
      call get_command_argument(i,outFiles(i))
   end do
   allocate(ens(argc),tds((argc-1)/2))
   ens=0
   
   call moleculeDims(outFiles(1),nat)
   NQ=3*nat-6
   if(linearMol)NQ=3*nat-5
   dir_gr=outfiles(1)
   call ReplaceSuffix(dir_gr,' ','/')
   istat = CHDIR(TR(dir_gr)//'/ground')
   if(istat/=0)then
      write(output_unit,*)'FORTRAN ERROR AT:'
      write(output_unit,*)'CHDIR(',TR(dir_gr),')'
      call exit(19)
   end if
   
   open(77,file='ENERGY',action='READ')
   read(77,*)ens(1)
   close(77)

   
   
   ! allocate(wg(NQ))
   ! call readsi(3*nat,wg,nq,z,mas,'F.INP')
   !!remove imaginary frequencies from 'wg'
   ! imag_c=count(wg<0d0)
   ! wg_old=wg
   ! deallocate(wg)
   ! wg=NegFreqRemove(wg_old,nq,imag_c)
   ! nq_gr=nq-imag_c
   allocate(wexc_nm_str(nexc))
   do i = 1,nexc
      write(wexc_nm_str(i),'(I4,"nm")')NINT(wexc_nm(i))
   end do
   
   call chdir(dir_cwd)
   td_counter=1
   if(gamma_el==0d0)gamma_el=gamma
   do i = 2,argc,2
      td_i=td_counter
      td=>tds(td_i)
      td_counter=td_counter+1
      
      write(output_unit,*)'Ground geometry TDM derivatives'
      call ReadTDFreq(outfiles(i),td,nat,.true.,ens(i),ens(1),gr_exp)
      write(output_unit,*)'Excited geometry TDM derivatives and Duschinsky objects'
      call ReadTDFreq(outfiles(i+1),td,nat,.false.,ens(i+1),ens(1),gr_exp)
      !tds(td_i)%e_00=ens(i+1)-ens(1)
      flush(output_unit)
      if(vel)then
         write(Output_unit,*)'Using velocity dipole'
         open(77,file='VEL.LEN.DIFF')
         write(77,*)'Ground TDM, Len. - Vel.: '
         write(77,*)(td%u_gr+0.00001)/(td%v_gr+0.00001)
         write(77,*)'Excited TDM, Len. - Vel.: '
         write(77,*)(td%u_ex+0.00001)/(td%v_ex+0.00001)
         write(77,*)'Ground DTDM, Len. - Vel.: '
         write(77,*)(td%du_gr+0.0001)/(td%dv_gr+0.00001)
         write(77,*)'Excited DTDM, Len. - Vel.: '
         write(77,*)(td%du_ex+0.0001)/(td%dv_ex+0.00001)
         close(77)
         
         
         td%u_gr=td%v_gr
         td%u_ex=td%v_ex
         td%du_gr=td%dv_gr
         td%du_ex=td%dv_ex
         if(ht2)then
         td%du2_gr=td%dv2_gr
         td%du2_ex=td%dv2_ex
         end if
      end if
      
      if(MomentsExToGr)then
         write(output_unit,*)'Transition moments in ground are inherited from excited'
         td%u_gr=td%u_ex
         td%v_gr=td%v_ex
         td%m_gr=td%m_ex
         td%q_gr=td%q_ex
         
         td%du_gr=td%du_ex
         td%dv_gr=td%dv_ex
         td%dm_gr=td%dm_ex
         td%dq_gr=td%dq_ex
      end if
   end do
   wg=tds(1)%wg !TODO workaround
   nq_gr=tds(1)%nq
   istat=chdir(dir_cwd)
   !qpar=qpar_vib(wg,NQ_gr,temp) 
   !thresholdOrig_gr=1d-1
   !thresholdOrig_ex=1d-1
   if(elpol_only)then
      write(output_unit,*)'Only ELPOL was requested.'
      call UnmakeSqrtArr()
      bigTime=getTime()-bigTime
      write(output_unit,*)
      write(output_unit,*)'----------------------------------------'
      write(output_unit,*)'Walltime:    ',GetTimeStr_butBig(bigTime,.false.)
      write(output_unit,*)'Finished on: ',ctime(time8())
      if(verbose)write(output_unit,*)'Program terminated normally.'
      return
   end if
   if(stokes .and. mc_fs==0)then
      write(output_unit,*)'The maximum class for final state (MAX_CLASS_FS) cannot be 0 for Stokes transitions.'
      call exit(5)
   end if
   
   
   mc_gfs=max(mc_gs,mc_fs)
   do i = 1,(argc-1)/2
      td_i=i
      td=>tds(td_i)
      
      td%thr_v1=thresholdOrig_gr
      td%thr_v2=thresholdOrig_ex
      td%thr_v2_0=thr_v2_0
      td%thr_v2_over=thr_v2_ov
      td%thr_v2_comb=thr_v2_comb
      td%thr_v2_other=thr_v2_other
      td%gamma=gamma
      td%theta=theta
      td%eps=td_eps
      td%TMExpand=TMExpand_def
      td%mc_gs=mc_gs
      td%mc_fs=mc_fs
      td%mc_ms=mc_ms
      td%mc_gfs=max(td%mc_gs,td%mc_fs)
      td%min_v1s=min_v1s
      td%min_v2s=min_v2s
      td%min_v2s_excm=min_v2s_excm
      td%max_v1s=max_v1s
      td%max_v3s=max_v3s
      td%max_v2s=max_v2s
      td%max_v2s_ps=max_v2s_ps
      td%mc_ms_ps=mc_ms_ps
      td%min_freq_mode=min_freq_mode
      td%max_freq_mode=max_freq_mode
      td%antiStokes=antistokes
      td%red_n2_max=red_n2_max
      td%red_n2_max_algo=red_n2_max_algo
      td%n_gr=n_gr
      if(fc_sum_min_ps>0)td%fc_sum_min_ps=fc_sum_min_ps
      td%fc_sum_min_ps_its=fc_sum_min_ps_its
      
      td%we=td%we*cm_2_au
      td%wg=td%wg*cm_2_au
      td%e_00=td%e_00-sum(td%wg)/2d0+sum(td%we)/2d0 !ZPE correction
      
      
      if(transform_excited_tm)then !unlikely to help with RROA, will implement when somebody wants fluorescence from this program
         call TransformTM(td,.false.,transform_tm_type)
      end if
      !transform gr. transtion moment to excited state du/dQ''=J*du/dQ'
      !this comes from the transformation of coordinates in derivatives dQ'/dQ'' = J
      !or from the TI approach when: <v'|du/dQ' * Q'|v''> = sum(<v'|du/dQ' * Q'|v''>)
      if(transform_ground_tm)then 
         call TransformTM(td,.true.,transform_tm_type)
      end if
      
      if(.not.td_approach)then
         call MakeTransTM(td)
      else 
         if(td_fs==0d0)then
            td_fs=(2_8**TD_N_POINTS)*2*pi/td_tmax
         else
            td_tmax=(2_8**TD_N_POINTS)/(td_fs/(2*pi))
         end if
      end if
      
      
      td%sqrt_w=dsqrt(1d0/(2d0*td%we))
      td%sqrt_wg=dsqrt(1d0/(2d0*td%wg))
      
      td%thr_v1=td%thr_v1*td%fc_00
      td%thr_v2=td%thr_v2*td%fc_00
      td%thr_v2_0=td%thr_v2_0*td%fc_00
      td%thr_v2_over=td%thr_v2_over*td%fc_00
      td%thr_v2_comb=td%thr_v2_comb*td%fc_00
      td%thr_v2_other=td%thr_v2_other*td%fc_00
      td%mc_gfs=max(td%mc_gs,td%mc_fs)
      td%ExpandInExc1=(td%TMExpand(1:1)=='E')
      td%ExpandInExc2=(td%TMExpand(2:2)=='E')
      
      
      sz=size(td%min_v2s,dim=1)
      if(sz<td%mc_ms)then !TODO: for other min maxes too
         allocate(Iarr(td%mc_ms))
         Iarr=0
         Iarr(1:sz)=td%min_v2s
         deallocate(td%min_v2s)
         td%min_v2s=Iarr
         deallocate(Iarr)
      end if
      
      if(stokes .and. td%mc_fs==0)then
         write(output_unit,*)'The maximum class for final state (MAX_CLASS_FS) cannot be 0 for Stokes transitions.'
         call exit(5)
      end if
      
       
      td%mc_v1=td%mc_gfs
      td%mc_v2=td%mc_ms
   end do
   
   if(overwrite_e00/=0d0)then
      tds(:)%e_00=overwrite_e00
   else if(avg_e00)then
      tds(:)%e_00=sum(tds(:)%e_00)/td_n
   end if
   
   if(wexc_eq_e00)then
      deallocate(wexc)
      if(avg_e00)then
         allocate(wexc(1))
         wexc=tds(1)%e_00
         nexc=1
      else
         allocate(wexc(td_n))
         wexc=tds(:)%e_00
         nexc=td_n
      end if
   end if
   
   if(TR(sp_type)=='RROA')then
      if(verbose)then
         write(output_unit,*)
         write(output_unit,*)
         write(output_unit,*)
         write(output_unit,'(15X,A)')'------------------------------'
         write(output_unit,'(15X,A)')'|       Resonance ROA        |'
         write(output_unit,'(15X,A)')'------------------------------'
      else
         write(output_unit,*)'RROA' 
      end if
      flush(output_unit)
     ! uncoup_modes_maxval=-1
      call DoRROA(n_thr,runtype,tds,size(tds,dim=1),wg,nq_gr,wexc,nexc,kbt_lim,temp,wmin,wmax,npx,fwhm,rroa_spr_do,rroa_st,ht,ht2,only_ht,useA,useG, &
                  fcarrcheck,use_gauss,write_excm,write_ten,write_inv,write_fc_sys,write_fcarr,.false.,uncoup_lim,uncoup_modes_maxval,w_ps,ignore_ps,sel_rules,sel_rules2,correct_gr_freqs, &
                  add_pol_coeff,output_moments,output_polars,ignore_imag_modes,output_polContrs,AlphaFC,wr_is_e00,spectrum_temp,doModes,trs,fundLeadCount, &
                  td_approach,td_n_points,td_tmax,td_fs,td_sparse,td_alt,td_fixphase,td_J_tol,td_batch_n,td_trueground,write_corrf,write_fft_cf,num_integ,norm_fft,correctPhaseX,correctPhaseX_abs,interpolateFFT,wexc_adapt,w_ad_zero, &
                  evcd,cpl,contr)
      write(output_unit,*)'Finished RROA'
   end if
   
   
   if(runtype/='DEF')then
      if(verbose)then
         write(output_unit,*)'No spectra will be made.'
         if(runtype=='SI1')write(output_unit,*)'Only FC matrix sizes were desired.'
         if(runtype=='SI2')write(output_unit,*)'Only FC sums were desired.'
         !if(runtype=='SI3')write(output_unit,*)'DEPRECATED'
         write(output_unit,*)'Program terminated normally.'
      else
         write(output_unit,*)runtype
      end if
      call UnmakeSqrtArr()
      bigTime=getTime()-bigTime
      write(output_unit,*)
      write(output_unit,*)'----------------------------------------'
      write(output_unit,*)'Walltime:    ',GetTimeStr_butBig(bigTime,.false.)
      write(output_unit,*)'Finished on: ',ctime(time8())
      if(verbose)write(output_unit,*)'Program terminated normally.'
      if(stdoutToFile)close(output_unit)
      
      call exit(0)
   end if
      
      
   
   
   !Program clean up
   ! if(.not.only_dusch_and_der)then
      ! if(verbose)write(output_unit,*)'Removing scratch directories: "ground" ,"excited"'
      ! istat = System('rm -r ground')
      ! istat = System('rm -r excited')
   ! end if
   call UnmakeSqrtArr()
   bigTime=getTime()-bigTime
   write(output_unit,*)
   write(output_unit,*)'----------------------------------------'
   write(output_unit,*)'Walltime:    ',GetTimeStr_butBig(bigTime,.false.)
   write(output_unit,*)'Finished on: ',ctime(time8())
   if(verbose)write(output_unit,*)'Program terminated normally.'
   if(stdoutToFile)close(output_unit)
   call exit(0)
   
   contains
   
   subroutine TransformTM(td,isGround,trfType)
      type(ExcState) td 
      logical isGround
      integer nq,i,j,a,b
      integer trfType
      double precision,allocatable :: JtK(:)
      double precision bufr
      
      nq=td%nq
      allocate(Jtk(nq))
      Jtk=0d0
      
      if(isGround)THEN
         !Transition moments
         td%u_ex=td%u_gr
         td%m_ex=td%m_gr
         td%q_ex=td%q_gr
         if(trfType>=1)then
            do i = 1,nq
               td%u_ex=td%u_ex+td%du_gr(:,i)*td%K(i)
               td%m_ex=td%m_ex+td%dm_gr(:,i)*td%K(i)
               td%q_ex=td%q_ex+td%dq_gr(:,:,i)*td%K(i)
               if(.not.ht2 .or. .not.trfType>=2)cycle
               do j = 1,nq
                  td%u_ex=td%u_ex+td%du2_gr(:,j,i)*td%K(j)*td%K(i)
                  td%m_ex=td%m_ex+td%dm2_gr(:,j,i)*td%K(j)*td%K(i)
                  td%q_ex=td%q_ex+td%dq2_gr(:,:,j,i)*td%K(j)*td%K(i)
               end do
            end do
         end if
         !Transition moment first derivatives
         do a = 1,nq
            td%du_ex(:,a)=0
            td%dm_ex(:,a)=0
            td%dq_ex(:,:,a)=0
            do i = 1,nq
               td%du_ex(:,a)=td%du_ex(:,a)+td%J(i,a)*td%du_gr(:,i)
               td%dm_ex(:,a)=td%dm_ex(:,a)+td%J(i,a)*td%dm_gr(:,i)
               td%dq_ex(:,:,a)=td%dq_ex(:,:,a)+td%J(i,a)*td%dq_gr(:,:,i)
            end do
            if(.not.ht2 .or. .not.trfType>=2)cycle
            do i = 1,nq
               do j = 1,nq
                  bufr=2*td%J(j,a)*td%K(i)
                  td%du_ex(:,a)=td%du_ex(:,a)+td%du2_gr(:,j,i)*bufr
                  td%dm_ex(:,a)=td%dm_ex(:,a)+td%dm2_gr(:,j,i)*bufr
                  td%dq_ex(:,:,a)=td%dq_ex(:,:,a)+td%dq2_gr(:,:,j,i)*bufr
               end do 
            end do 
         end do
         
         if(.not.ht2)return
         !Transition moment second derivatives
         td%du2_ex=0d0
         td%dm2_ex=0d0
         td%dq2_ex=0d0
         !$OMP PARALLEL DO DEFAULT(NONE) &
         !$OMP PRIVATE(a,b,i,j,bufr) &
         !$OMP SHARED(nq,td)
         do a = 1,nq
            do b = 1,nq
               do i =1,nq
                  do j = 1,nq
                     bufr=td%J(j,a)*td%J(i,b)
                     td%du2_ex(:,b,a)=td%du2_ex(:,b,a)+td%du2_gr(:,i,j)*bufr
                     td%dm2_ex(:,b,a)=td%dm2_ex(:,b,a)+td%dm2_gr(:,i,j)*bufr
                     td%dq2_ex(:,:,b,a)=td%dq2_ex(:,:,b,a)+td%dq2_gr(:,:,i,j)*bufr
                  end do
               end do
            end do
         end do
         !$OMP END PARALLEL DO
      else
         stop 66
      end if
      deallocate(jtk)
   end subroutine TransformTM
   
   subroutine ReadFileTR(trs)
      type(Transition_R),allocatable :: trs(:)
      character(80) :: s80
      integer row_c
      
      open(77,file='FILE.TR',status='old')
      
      row_c=0
      do while(.true.)
         read(77,'(A80)',end=200)s80
         row_c=row_c+1
      end do
      
200   allocate(trs(row_c))
      rewind(77)
      do i = 1,row_c
         read(77,'(A80)')s80
         trs(i)=Str2Tr_R(s80)
      end do
      
300   close(77)
   end subroutine ReadFileTR
   
   function Str2Tr_R(str)result(tr)
      integer i,vi_class,vf_class
      integer,allocatable :: vi(:),vi_pos(:),vf(:),vf_pos(:)
      character(*) str
      type(Transition_r) tr
      character(80) vi_str,vf_str
      character(80),allocatable :: Str_arr(:),str_arr2(:)
      
      str_arr=splitString(str,80,'>')
      vi_str=TR(str_arr(1))
      vf_str=TR(str_arr(2))
      !read to v_initial
      str_arr=splitString(TR(vi_str),len(TR(vi_str)),' ')
      vi_class=size(str_arr,dim=1)
      allocate(vi(vi_class),vi_pos(vi_class))
      do i = 1,vi_class
         if(index(str_arr(i),'^')==0 .and. index(str_arr(i),'0')/=0)then
            vi_pos=[0]
            vi=[0]
            exit
         end if
         str_arr2=splitString(str_arr(i),80,'^')
         read(str_arr2(1),*)vi_pos(i)
         read(str_arr2(2),*)vi(i)
      end do
      tr%vi=vi
      tr%vi_pos=vi_pos
      
      str_arr=splitString(TR(vf_str),len(TR(vf_str)),' ')
      vf_class=size(str_arr,dim=1)
      allocate(vf(vf_class),vf_pos(vf_class))
      do i = 1,vf_class
         str_arr2=splitString(str_arr(i),80,'^')
         read(str_arr2(1),*)vf_pos(i)
         read(str_arr2(2),*)vf(i)
      end do
      tr%vf=vf
      tr%vf_pos=vf_pos
      
      deallocate(str_arr,str_arr2)
      deallocate(vi,vi_pos,vf,vf_pos)
   end function Str2Tr_R
   
   subroutine ReadTDFreq(filee,td,nat,isGroundGeom,en,en_gr,gr_exp)
      integer istat,j,nroot,line,nat2,nat,i,ii,nq,n3,idxv,idxq
      logical isGroundGeom,fex,gr_exp,found_tms,found_tms2
      double precision en,eau,en_gr
      double precision u_copy(3),v_copy(3),m_copy(3),q_copy(3,3)
      double precision,allocatable :: dd(:),aa(:),vv(:),qq(:),smat(:,:),e(:,:)
      double precision,allocatable :: du2_c(:,:,:),dv2_c(:,:,:),dm2_c(:,:,:),dq2_c(:,:,:,:)
      double precision,allocatable :: du_c(:,:),dv_c(:,:),dm_c(:,:),dq_c(:,:,:)
      type(ExcState) td
      character(*) filee
      character(300) dir_ex
      
      found_tms=.false.
      dir_ex=filee
      n3=3*nat
      call ReplaceSuffix(dir_ex,' ','/')
      write(output_unit,*)'Reading files from ',dir_ex

      call moleculeDims(TR(filee),nat2)
      if(nat/=nat2)then
         write(output_unit,*)'Error: different atom count in ground '//TR(filee)//' and excited '//TR(filee)
         call exit(1)
      end if
      istat = CHDIR(TR(dir_ex)//'/excited')
      open(77,file='ENERGY',action='READ')
      read(77,*)en
      close(77)
      istat = CHDIR('..')
      td%e_00=en-en_gr
      allocate(dd(9*nat),aa(9*nat),vv(9*nat),qq(18*nat))
      
      if(isGroundGeom)then
         call readDusch(td%nq_gr,3*nat,td%A_gr,td%B_gr,td%C_gr,td%D_gr,td%E_gr,td%J_gr,td%K_gr,td%SG_gr,td%SE_gr,td%wg,td%we_gr,'DUSCH.OUT')
         nq=td%nq_gr
         allocate(td%wg_gr_idx(nq))
         open(792,file='KEPT_MODES_GR')
         do i = 1,nq
            read(792,*)td%wg_gr_idx(i)
         end do
         close(792)
         
         
         if(ignore_imag_modes)THEN
            allocate(td%im_modes_gr(nq))
            open(792,file='IMAG_MODES_GR')
            td%im_modes_gr=.false.
            do i = 1,nq
               read(792,*,end=1000)ii
               td%im_modes_gr(ii)=.true.
            end do
1000        close(792)

            allocate(td%im_modes_ex_gr(nq))
            open(792,file='IMAG_MODES_EX')
            td%im_modes_ex_gr=.false.
            do i = 1,nq
               read(792,*,end=1010)ii
               td%im_modes_ex_gr(ii)=.true.
            end do
1010        close(792)


         end if
         allocate(e(nq,nq*2),td%J_gr_i(nq,nq))
         call Inv(td%J_gr,td%J_gr_i,nq,e,ierr)
         deallocate(e)
         td%K_gr_i=-matmul(td%J_gr_i,td%K_gr)
      else
         call readDusch(td%nq_ex,3*nat,td%A,td%B,td%C,td%D,td%E,td%J,td%K,td%SG,td%SE,td%wg,td%we,'DUSCH.OUT')
         nq=td%nq_ex
         allocate(td%wg_ex_idx(nq))
         open(792,file='KEPT_MODES_GR')
         do i = 1,nq
            read(792,*)td%wg_ex_idx(i)
         end do
         close(792)
         
         if(ignore_imag_modes)THEN
            allocate(td%im_modes_ex_ex(nq))
            open(792,file='IMAG_MODES_EX')
            td%im_modes_ex_ex=.false.
            do i = 1,nq
               read(792,*,end=1020)ii
               td%im_modes_ex_ex(ii)=.true.
            end do
1020        close(792)
         end if
         
         allocate(e(nq,nq*2),td%J_i(nq,nq))
         call Inv(td%J,td%J_i,nq,e,ierr)
         deallocate(e)
         td%K_i=-matmul(td%J_i,td%K)
      end if
      
      
      
      
      call CW('Reading transition derivatives')
      
      allocate(td%du_raw(9*nat),td%dm_raw(9*nat),td%dv_raw(9*nat),td%dq_raw(18*nat))
      line=0
      inquire(file='TM.U',exist=found_tms)
      inquire(file='DTM2.U',exist=found_tms2)
      if(isGroundGeom)then
         smat=td%SG_gr
         if(found_tms)then
            call ReadTMS(td%u_gr,td%v_gr,td%m_gr,td%q_gr)
            call ReadDTMS(du_c,dv_c,dm_c,dq_c,n3)
            td%du_gr=Car2NM(n3,nq,du_c,smat)
            td%dv_gr=Car2NM(n3,nq,dv_c,smat)
            td%dm_gr=Car2NM(n3,nq,dm_c,smat)
            td%dq_gr=Car2NM_Q(n3,nq,dq_c,smat)
            td%root=-1
            
         else
            call Rdd(ShaveFront(filee,'/'),dd,nat,td%u_gr,td%m_gr,aa,vv,td%v_gr,qq,td%q_raw,line,nroot,eau,mult_tm,mult_dtm)
            td%q_gr=TracelessQ(td%q_raw)
            td%root=nroot
         
            !According to some paper scribbling ive done: du_el/du_x' = du_el/du_x''
            ! if(derivatives_in_excGeom)then
               ! smat=td%SE_gr
            ! else
            !smat=td%SG_gr
            ! end if
            
            call trafN_new(3,dd,td%du_raw,nat,nq,smat)
            call trafN_new(3,aa,td%dm_raw,nat,nq,smat)
            call trafN_new(6,qq,td%dq_raw,nat,nq,smat)
            call trafN_new(3,vv,td%dv_raw,nat,nq,smat)
            
            allocate(du_c(3,n3),dv_c(3,n3),dm_c(3,n3),dq_c(3,3,n3))
            do i = 1,n3
               idxv=(i-1)*3+1
               idxq=(i-1)*6+1
               du_c(:,i)=dd(idxv:idxv+2)
               dm_c(:,i)=aa(idxv:idxv+2)
               dv_c(:,i)=vv(idxv:idxv+2)
               dq_c(:,:,i)=TracelessQ(qq(idxq:idxq+5))
            end do
            
            call TidyDerivatives(td%du_raw,9*nat,td%du_gr,3,NQ)
            call TidyDerivatives(td%dm_raw,9*nat,td%dm_gr,3,NQ)
            call TidyDerivativesQ(td%dq_raw,18*nat,td%dq_gr,6,NQ)
            call TidyDerivatives(td%dv_raw,9*nat,td%dv_gr,3,NQ)
               
         end if
         
         if(found_tms2.and.ht2)then
            call ReadDTMS2(du2_c,dv2_c,dm2_c,dq2_c,n3)
            td%du2_gr=0.5d0*Car2NM_2(n3,nq,du2_c,smat)
            td%dv2_gr=0.5d0*Car2NM_2(n3,nq,dv2_c,smat)
            td%dm2_gr=0.5d0*Car2NM_2(n3,nq,dm2_c,smat)
            td%dq2_gr=0.5d0*Car2NM_Q_2(n3,nq,dq2_c,smat)
            
            !force symmetry
            do i = 1,nq
               do j = i+1,nq
                  td%du2_gr(:,j,i)=td%du2_gr(:,i,j)
                  td%dm2_gr(:,j,i)=td%dm2_gr(:,i,j)
                  td%dv2_gr(:,j,i)=td%dv2_gr(:,i,j)
                  td%dq2_gr(:,:,j,i)=td%dq2_gr(:,:,i,j)
               end do
            end do
         elseif(ht2)then
            write(output_unit,*)'HT2 was desired but no 2nd derivatives found'
         end if
         
         open(77,file='ENERGY.VERTICAL',status='old')
         read(77,*)td%e00_vertical
         close(77)
         
         call VelToLen_dip(td%v_gr,td%dv_gr,td%e00_vertical,nq,.false.,.false.,vel_grad,coordinates,.false.)
         call VelToLen_quad(td%q_gr,td%dq_gr,td%e00_vertical,nq,.false.,.false.,vel_grad,coordinates,.false.)
         ! if(vel)then
            ! write(Output_unit,*)'Using velocity dipole in ground'
            ! td%u_gr=td%v_gr
            ! td%du_gr=td%dv_gr
         ! end if
         
         if(.not.elpol_exc)then
            u_copy=td%u_gr
            m_copy=td%m_gr
            v_copy=td%v_gr
            q_copy=td%q_gr
            call ElPol(td,nat,smat,nq,td%e00_vertical,wexc,nexc,gamma_el,u_copy,v_copy,m_copy,q_copy,du_c,dv_c,dm_c,dq_c,ht2,du2_c,dv2_c,dm2_c,dq2_c,vel,.true.,elpol_grad,found_tms)
         end if
      else
         smat=td%SE
         if(found_tms)then
            call ReadTMS(td%u_ex,td%v_ex,td%m_ex,td%q_ex)
            call ReadDTMS(du_c,dv_c,dm_c,dq_c,n3)
            td%du_ex=Car2NM(n3,nq,du_c,smat)
            td%dv_ex=Car2NM(n3,nq,dv_c,smat)
            td%dm_ex=Car2NM(n3,nq,dm_c,smat)
            td%dq_ex=Car2NM_Q(n3,nq,dq_c,smat)
            !deallocate(du_c,dv_c,dm_c,dq_c)
            td%root=-1
         else
            call Rdd(ShaveFront(filee,'/'),dd,nat,td%u_ex,td%m_ex,aa,vv,td%v_ex,qq,td%q_raw,line,nroot,eau,mult_tm,mult_dtm)
            td%q_ex=TracelessQ(td%q_raw)
            td%root=nroot
         
            !According to some paper scribbling ive done: du_el/du_x' = du_el/du_x''
            ! if(derivatives_in_excGeom)then
            if(SE_is_SG)smat=td%SG_gr
            ! else
               ! smat=td%SG
            ! end if
            call trafN_new(3,dd,td%du_raw,nat,nq,smat)
            call trafN_new(3,aa,td%dm_raw,nat,nq,smat)
            call trafN_new(6,qq,td%dq_raw,nat,nq,smat)
            call trafN_new(3,vv,td%dv_raw,nat,nq,smat)
            
            allocate(du_c(3,n3),dv_c(3,n3),dm_c(3,n3),dq_c(3,3,n3))
            do i = 1,n3
               idxv=(i-1)*3+1
               idxq=(i-1)*6+1
               du_c(:,i)=dd(idxv:idxv+2)
               dm_c(:,i)=aa(idxv:idxv+2)
               dv_c(:,i)=vv(idxv:idxv+2)
               dq_c(:,:,i)=TracelessQ(qq(idxq:idxq+5))
            end do
            
            call TidyDerivatives(td%du_raw,9*nat,td%du_ex,3,NQ)
            call TidyDerivatives(td%dm_raw,9*nat,td%dm_ex,3,NQ)
            call TidyDerivativesQ(td%dq_raw,18*nat,td%dq_ex,6,NQ)
            call TidyDerivatives(td%dv_raw,9*nat,td%dv_ex,3,NQ)
         end if
         
         if(found_tms2.and.ht2)then
            call ReadDTMS2(du2_c,dv2_c,dm2_c,dq2_c,n3)
            td%du2_ex=0.5d0*Car2NM_2(n3,nq,du2_c,smat)
            td%dv2_ex=0.5d0*Car2NM_2(n3,nq,dv2_c,smat)
            td%dm2_ex=0.5d0*Car2NM_2(n3,nq,dm2_c,smat)
            td%dq2_ex=0.5d0*Car2NM_Q_2(n3,nq,dq2_c,smat)
            !force symmetry
            do i = 1,nq
               do j = i+1,nq
                  td%du2_gr(:,j,i)=td%du2_gr(:,i,j)
                  td%dm2_gr(:,j,i)=td%dm2_gr(:,i,j)
                  td%dv2_gr(:,j,i)=td%dv2_gr(:,i,j)
                  td%dq2_gr(:,:,j,i)=td%dq2_gr(:,:,i,j)
               end do
            end do
         elseif(ht2)then
            write(output_unit,*)'HT2 was desired but no 2nd derivatives found'
         end if
         
         
         open(77,file='ENERGY.VERTICAL',status='old')
         read(77,*)td%e00_vertical_exc
         close(77)
         
         call VelToLen_dip(td%v_ex,td%dv_ex,td%e00_vertical_exc,nq,isAH,AHAS,vel_grad,coordinates,.false.)
         call VelToLen_quad(td%q_ex,td%dq_ex,td%e00_vertical_exc,nq,isAH,AHAS,vel_grad,coordinates,.false.)
         ! if(vel)then
            ! write(Output_unit,*)'Using velocity dipole in ground'
            ! td%u_ex=td%v_ex
            ! td%du_ex=td%dv_ex
         ! end if
         open (77,file='OVERLAP',action='READ')
         read(77,*)td%FC_00
         close(77)
         
         
         td%nq=td%nq_ex
         if(elpol_exc)then
            u_copy=td%u_ex
            m_copy=td%m_ex
            v_copy=td%v_ex
            q_copy=td%q_ex
            call ElPol(td,nat,smat,nq,td%e00_vertical_exc,wexc,nexc,gamma_el,u_copy,v_copy,m_copy,q_copy,du_c,dv_c,dm_c,dq_c,ht2,du2_c,dv2_c,dm2_c,dq2_c,vel,.false.,elpol_grad,found_tms)
         end if
      end if
      deallocate(du_c,dv_c,dm_c,dq_c)
      if(ht2)then
         deallocate(du2_c,dv2_c,dm2_c,dq2_c)
      end if
      deallocate(aa,dd,vv,qq)
      deallocate(td%du_raw,td%dm_raw,td%dv_raw,td%dq_raw)
      call chdir(dir_cwd)
   end subroutine ReadTDFreq
   
   subroutine ReadTMS(u,v,m,q)
      double precision u(3),v(3),m(3),q(3,3)
      
      call ReadVec(77,'TM.U',u)
      call ReadVec(77,'TM.V',v)
      call ReadVec(77,'TM.M',m)
      call ReadVec(77,'TM.Qtr',q)
   end subroutine ReadTMS
   
   subroutine ReadDTMS(du,dv,dm,dq,n3)
      double precision,allocatable :: du(:,:),dv(:,:),dm(:,:),dq(:,:,:)
      integer n3
      
      call ReadDerVec(77,'DTM.U',du,n3)
      call ReadDerVec(77,'DTM.V',dv,n3)
      call ReadDerVec(77,'DTM.M',dm,n3)
      call ReadDerQ(77,'DTM.Qtr',dq,n3)
   end subroutine ReadDTMS
   
   subroutine ReadDTMS2(du,dv,dm,dq,n3)
      double precision,allocatable :: du(:,:,:),dv(:,:,:),dm(:,:,:),dq(:,:,:,:)
      integer n3
      
      call TDD2_read(77,'DTM2.U',du,n3)
      call TDD2_read(77,'DTM2.V',dv,n3)
      call TDD2_read(77,'DTM2.M',dm,n3)
      call TDD2_readq(77,'DTM2.Qtr',dq,n3)
   end subroutine ReadDTMS2
   
   
   
   subroutine MakeTransTM(td)
      integer i,k
      type(ExcState) td
      
      td%u_gr_tr=td%u_gr
      td%m_gr_tr=td%m_gr
      td%q_gr_tr=td%q_gr
      
      td%u_ex_tr=td%u_ex
      td%m_ex_tr=td%m_ex
      td%q_ex_tr=td%q_ex
      
      allocate(td%du_gr_tr(3,td%nq),td%dm_gr_tr(3,td%nq),td%dq_gr_tr(3,3,td%nq))
      allocate(td%du_ex_tr(3,td%nq),td%dm_ex_tr(3,td%nq),td%dq_ex_tr(3,3,td%nq))
      td%du_gr_tr=0d0
      td%dm_gr_tr=0d0
      td%dq_gr_tr=0d0
      
      td%du_ex_tr=0d0
      td%dm_ex_tr=0d0
      td%dq_ex_tr=0d0
      
      
      do i = 1,td%nq
         td%u_gr_tr=td%u_gr_tr+td%du_gr(:,i)*td%K(i) !Do I use the ground state geometry Duschinsky matrices? Probably not
         td%m_gr_tr=td%m_gr_tr+td%dm_gr(:,i)*td%K(i)
         td%q_gr_tr=td%q_gr_tr+td%dq_gr(:,:,i)*td%K(i)
         
         td%u_ex_tr=td%u_ex_tr+td%du_ex(:,i)*td%K_i(i)
         td%m_ex_tr=td%m_ex_tr+td%dm_ex(:,i)*td%K_i(i)
         td%q_ex_tr=td%q_ex_tr+td%dq_ex(:,:,i)*td%K_i(i)
         
         
         do k = 1,td%nq
            td%du_gr_tr(:,i)=td%du_gr_tr(:,i)+td%du_gr(:,k)*td%J(k,i)
            td%dm_gr_tr(:,i)=td%dm_gr_tr(:,i)+td%dm_gr(:,k)*td%J(k,i)
            td%dq_gr_tr(:,:,i)=td%dq_gr_tr(:,:,i)+td%dq_gr(:,:,k)*td%J(k,i)
            
            td%du_ex_tr(:,i)=td%du_ex_tr(:,i)+td%du_ex(:,k)*td%J_i(k,i)
            td%dm_ex_tr(:,i)=td%dm_ex_tr(:,i)+td%dm_ex(:,k)*td%J_i(k,i)
            td%dq_ex_tr(:,:,i)=td%dq_ex_tr(:,:,i)+td%dq_ex(:,:,k)*td%J_i(k,i)
         end do
         td%du_gr_tr(:,i)=td%du_gr_tr(:,i)*(1d0/(sqrt(2d0*td%we(i))))
         td%dm_gr_tr(:,i)=td%dm_gr_tr(:,i)*(1d0/(sqrt(2d0*td%we(i))))
         td%dq_gr_tr(:,:,i)=td%dq_gr_tr(:,:,i)*(1d0/(sqrt(2d0*td%we(i))))
         
         td%du_ex_tr(:,i)=td%du_ex_tr(:,i)*(1d0/(sqrt(2d0*td%wg(i))))
         td%dm_ex_tr(:,i)=td%dm_ex_tr(:,i)*(1d0/(sqrt(2d0*td%wg(i))))
         td%dq_ex_tr(:,:,i)=td%dq_ex_tr(:,:,i)*(1d0/(sqrt(2d0*td%wg(i))))
      end do
      
      do i = 1,td%nq
         td%du_gr(:,i)=td%du_gr(:,i)*(1d0/(sqrt(2d0*td%wg(i))))
         td%dm_gr(:,i)=td%dm_gr(:,i)*(1d0/(sqrt(2d0*td%wg(i))))
         td%dq_gr(:,:,i)=td%dq_gr(:,:,i)*(1d0/(sqrt(2d0*td%wg(i))))
         
         td%du_ex(:,i)=td%du_ex(:,i)*(1d0/(sqrt(2d0*td%we(i))))
         td%dm_ex(:,i)=td%dm_ex(:,i)*(1d0/(sqrt(2d0*td%we(i))))
         td%dq_ex(:,:,i)=td%dq_ex(:,:,i)*(1d0/(sqrt(2d0*td%we(i))))
      end do
   end subroutine MakeTransTM
   
   
   
   subroutine ReadOpt()
      character(400) str,bufC
      double precision wexc_helper(200),wexc_l,wexc_r,dw
      integer nexc_help,buf(20),i,ii,iii,semicolon,comma,dash,nexc_help_new
      
      !init
      wexc_helper=0d0
      nexc_help=0
      buf=0
      
      !program defaults
      runtype='DEF'
      sp_type='RROA'
      TMExpand_def='EG'
      nexc=1
      wexc_nm=[532d0]
      wexc=10**7/wexc_nm*4.5564d-6
      wmin=0d0
      wmax=2000d0
      npx=2001
      mc_gs=0
      mc_fs=1
      mc_ms=3
      min_v1s=[0]
      max_v1s=[1,1]
      max_v3s=[1,1]
      min_v2s=[0]
      min_v2s_excm=[0]
      max_v2s=[30,30,5]
      verbose=.true.
      thresholdOrig_gr=1d-2
      thresholdOrig_ex=1d-2
      kfac=-1.0
      mc_ms_ps=2
      max_v2s_ps=[20,13]
      min_freq_mode=0
      max_freq_mode=huge(1d0)
      n_gr=1
      kBT_lim=1d0
      fwhm=10d0
      use_gauss=.false.
      vertical_en=.false.
      n_thr=1
      thr_D=1
      thr_r=1d-4
      fwhm_which=1
      uncoup_lim=1d0
      w_ps=.true.
      ignore_ps=.false.
      stdoutToFile=.false.
      sel_rules=.false.
      add_pol_coeff=0
      correct_gr_freqs=.false.
      !coordinates=.true.
      ! thresholdOrig=1d-2
      ! mc_gs=0
      ! mc_fs=1
      ! mc_ms=1
      ! gamma=150d0
      ! runtype='DEF'
      ! rroa_spr_do=.false.
      ! wmin=0d0
      ! wmin=4000d0
      ! npx=4001
      ! max_v1s=[1]
      ! max_v2s=[30]
      inquire(file='FCOV.OPT',exist=fex)
      if(.not.fex)then
         if(verbose)then
            write(output_unit,*)'FCOV.OPT not found, FCOV will use defaults.'
            return
         end if
      end if
      open(77,file='FCOV.OPT')
1000  read(77,'(A400)',end=1001)str
         buf=0
         select case(TR(str))
            ! case('SG_SWITCHPHASE')
               ! allocate(switchPhase(1000))
               ! switchPhase=0
               ! read(77,*,err=223,end=223)switchPhase
               ! 223 continue
               ! backspace(77)
               ! switchPhase_c=count(switchPhase/=0)
            case('TERMS')
               contr=.false.
               buf=0
               read(77,*,err=33,end=33)buf
33             continue               
               do i = 1,9
                  if(buf(i)==0)exit
                  contr(buf(i))=.true.
               end do
               buf=0
            case('SWITCH_MAG_SIGN')
               switchMagneticSign=.true.
            case('LINEAR')
               linearMol=.true.
            case('W_AD_ZERO')
               w_ad_zero=.true.
            case('ELPOL_ONLY')
               elpol_only=.true.
            case('NO_ELPOL_GRAD')
               elpol_grad=.false.
            case('SE_IS_SG')
               SE_is_SG=.true.
            case('WEXC_ADAPT')
               wexc_adapt=.true.
            case('TD')
               read(77,*)td_approach
            case('TD_ALT')
               read(77,*)td_alt
            case('TD_EPS')
               read(77,*)td_eps
            case('TD_NO_REPAIR_PHASE')
               td_fixphase=.false.
            case('TD_N_POINTS')
               read(77,*)td_N_points
            case('TD_MAXTIME')
               read(77,*)td_tmax
               !td_tmax=td_tmax/(2*pi)
            case('TD_SAMPLING_FREQ')
               read(77,*)td_fs
            case('TD_BATCH_N')
               read(77,*)td_batch_n
            case('TD_JTOL_N')
               read(77,*)td_j_tol_int
               td_j_tol=10d0**(-td_j_tol_int)
            case('TD_SPARSE')
               read(77,*)td_sparse
            case('TD_TRUEGROUND')
               td_trueground=.true.
            case('TD_NUMERICAL')
               num_integ=.true.
            case('TD_NORMALIZE_FFT')
               norm_fft=.true.
            case('TD_NO_CORRECT_PHASE_X')
               correctPhaseX=.false.
            case('TD_CORRECT_PHASE_X_ABS')
               correctPhaseX_abs=.true.
            case('TD_INTERPOLATE_FFT')
               interpolateFFT=.true.
            case('ELPOL_EXC')
               elpol_exc=.true.
            case('TM_GR_EXP')
               gr_exp=.true.
            case('WRITE_CORRF')
               write_corrf=.true.
            case('WRITE_FFT_CORRF')
               write_fft_cf=.true.
            case('PS_FC_SUM_MIN')
               read(77,*)fc_sum_min_ps
            case('PS_FC_SUM_MIN_ITS')
               read(77,*)fc_sum_min_ps_its
            case('TRANSFORM_GROUND_TM')
               transform_ground_tm=.true.
            case('TRANSFORM_EXCITED_TM')
               transform_excited_tm=.true.
            case('TRANSFORM_TM_TYPE')
               read(77,*)transform_tm_type
            case('MULT_TM')
               read(77,*)mult_tm
            case('MULT_DTM')
               read(77,*)mult_dtm
            case('COORDINATES')
               read(77,*)coordinates
            case('VEL_GRAD')
               read(77,*)vel_grad
            case('FUNDAMENTAL_LEAD_COUNT')
               read(77,*)fundLeadCount
            case('OVERWRITE_E00')
               read(77,*)overwrite_e00
               overwrite_e00=(1d7/overwrite_e00)*cm_2_au
            case('AHAS')
               read(77,*)AHAS
            case('DO_MODES')
               allocate(doModes(1000))
               doModes=0
               read(77,*,err=222,end=222)doModes
               222 continue
               backspace(77)
            case('SPECTRUM_TEMP')
               read(77,*)spectrum_temp
            case('WEXC_EQ_E00')
               read(77,*)wexc_eq_e00
            case('WR_EQ_E00')
               read(77,*)wr_is_e00
            case('AVERAGE_E00')
               read(77,*)avg_e00
            case('ALPHA_FC')
               read(77,*)AlphaFC
            case('MOMENTS_EX_TO_GR')
               read(77,*)MomentsExToGr
            case('ADIABATIC')
               read(77,*)isAH
            case('WRITE_POLAR_CONTRS')
               output_polContrs=.true.
            case('IGNORE_IMAG')
               read(77,*)ignore_imag_modes
            case('WRITE_MOMENTS')
               output_moments=.true.
            case('CORRECT_GROUND_FREQS')
               read(77,*)correct_gr_freqs
            case('ADD_POL_COEFF')
               read(77,*)add_pol_coeff
            case('SELECTION_RULES')
               read(77,*)sel_rules
            case('IGNORE_PS')
               read(77,*)ignore_ps
            case('OUTPUT_TO_FILE')
               read(77,*)stdoutToFile
            case('EVCD')
               evcd=.true.
            case('CPL')
               cpl=.true.
            ! case('OMP_N_THR')
               ! read(77,*)n_thr
            case('VERTICAL_EN')
               read(77,'(L)')vertical_en
            case('LINE_GAUSS')
               read(77,'(L)')use_gauss
            case('WRITE_EXC_MAXES')
               write_excm=.true.
            case('WRITE_TENSORS')
               write_ten=.true.
            case('WRITE_INVARIANTS')
               write_inv=.true.
            case('WRITE_FC_SYSTEM')
               write_fc_sys=.true.
            case('WRITE_FC_ARR')
               write_fcarr=.true.
            case('FWHM')
               read(77,*)fwhm
            case('USE_A')
               read(77,'(L)')useA
            case('USE_G')
               read(77,'(L)')useG
            case('ONLY_HT')
               read(77,'(L)')only_ht
            case('RUNTYPE')
               read(77,'(A3)')runType
            case('FC_MAT_CHECK')
               fcarrcheck=.true.
            case('HT')
               read(77,'(L)')ht
            case('HT2')
               read(77,'(L)')ht2
            case('UNCOUP_LIM')
               read(77,*)uncoup_lim
            case('UNCOUP_V')
               read(77,*)uncoup_modes_maxval
            case('WEXC_PS')
               read(77,'(L)')w_ps
            case('RROA_ST')
               read(77,'(L)')rroa_st
            case('VEL')
               vel=.true.
            case('THRESHOLD_00_GR')
               read(77,*)thresholdOrig_gr
            case('THRESHOLD_00_EX')
               read(77,*)thresholdOrig_ex
            case('THRESHOLD_00_EX_0')
               read(77,*)thr_v2_0
               thr_v2_0_found=.true.
            case('THRESHOLD_00_EX_OV')
               read(77,*)thr_v2_ov
               thr_v2_ov_found=.true.
            case('THRESHOLD_00_EX_COMB')
               read(77,*)thr_v2_comb
               thr_v2_comb_found=.true.
            case('THRESHOLD_00_EX_OTHER')
               read(77,*)thr_v2_other
               thr_v2_other_found=.true.
            case('WEXC_NM')
               read(77,'(A)')str
               dash=index(str,'-')
               if(dash>0)then
                  comma=index(str,',')
                  semicolon=index(str,';')
                  
                  read(str(1:dash-1),*)wexc_l
                  read(str(dash+1:comma-1),*)wexc_r
                  if(semicolon>0)then
                     read(str(comma+1:semicolon-1),*)nexc
                  else
                     read(str(comma+1:),*)nexc
                  end if
                  dw=(wexc_r-wexc_l)/(nexc-1)
                  do ii = 1,nexc
                     wexc_helper(ii)=(wexc_l+(ii-1)*dw)
                     !wexc(ii)=1d7/wexc_nm(ii)*cm_2_au
                  end do
                  wexc_nm=wexc_helper(1:nexc)
                  !wexc=1d7/wexc_nm*cm_2_au
                  deallocate(wexc)
                  allocate(wexc(nexc))
                  do i = 1,nexc
                     if(wexc_nm(i)==0d0)then
                        wexc(i)=0d0
                     else
                        wexc(i)=10**7/wexc_nm(i)*cm_2_au
                     end if
                  end do
                  
                  if(semicolon>0)then
                     wexc_helper=0
                     read(str(semicolon+1:),*,err=29,end=29)wexc_helper
29                   continue
                     nexc_help=count(wexc_helper/=0d0)
                     i=1
                     ii=0
                     nexc_help_new=nexc_help
                     outer: do while(ii<nexc)
                        ii=ii+1
                        do iii = 1,nexc_help
                           if(wexc_helper(iii)==wexc_nm(ii))cycle outer
                        end do
                        wexc_helper(nexc_help+i)=wexc_nm(ii)
                        i=i+1
                        nexc_help_new=nexc_help_new+1
                     end do outer
                     nexc_help=nexc_help_new
                     wexc_nm=wexc_helper(1:nexc_help)
                     nexc=nexc_help
                     call linsort_D_noord(wexc_nm,nexc,.false.)
                     !wexc=1d7/wexc_nm*cm_2_au
                     deallocate(wexc)
                     allocate(wexc(nexc))

                     do i = 1,nexc
                        if(wexc_nm(i)==0d0)then
                           wexc(i)=0d0
                        else
                           wexc(i)=10**7/wexc_nm(i)*cm_2_au
                        end if
                     end do
                  end if
               else
                  str=TR(str)
                  nexc=CountSubstring(TR(str),' ')+1
                  read(str,*,err=20,end=20)wexc_helper
20                continue
   !               BACKSPACE(77) !I hope this is portable
                  !I do thins because the list-directed format skips the next line for some reason
                  nexc=count(wexc_helper/=0d0)
                  wexc_nm=wexc_helper(1:nexc)
                  deallocate(wexc)
                  allocate(wexc(nexc))
                  do i = 1,nexc
                     if(wexc_nm(i)==0d0)then
                        wexc(i)=0d0
                     else
                        wexc(i)=10**7/wexc_nm(i)*cm_2_au
                     end if
                  end do
               end if
            case('TM_EXPANSION')
               read(77,*)TMExpand_def
            case('GAMMA')
               read(77,*)gamma
               gamma=gamma*cm_2_au
            case('GAMMA_EL')
               read(77,*)gamma_el
               gamma_el=gamma_el*cm_2_au
            case('THETA')
               read(77,*)theta
               theta=theta*cm_2_au
            case('MAX_CLASS_GS') !gr state
               read(77,*)mc_gs
            case('MAX_CLASS_FS') !final state
               read(77,*)mc_fs
            case('MAX_CLASS_MS') !intermediate state
               read(77,*)mc_ms
            case('MAX_RED_N2')
               read(77,*)red_n2_max
            case('MAX_RED_N2_ALGO')
               read(77,*)red_n2_max_algo
            case('RROA_SPR_TYPES')
               read(77,'(A200)')bufC 
               do i = 1,ns0
                  if(index(bufC,TR(rroa_exp(i)))>0)then
                     rroa_spr_do(i)=.true.
                  else
                     rroa_spr_do(i)=.false.
                  end if
               end do
            case('WMIN')
               read(77,*)wmin
            case('WMAX')
               read(77,*)wmax
            case('NPOINTS')
               read(77,*)npx
            case('MAX_V1S_CLASSES')
               if(allocated(max_v1s))deallocate(max_v1s)
               allocate(max_v1s(mc_gs))
               read(77,*)max_v1s
            case('MAX_V3S_CLASSES')
               if(allocated(max_v3s))deallocate(max_v3s)
               allocate(max_v3s(mc_fs))
               read(77,*)max_v3s
            case('MAX_V2S_CLASSES')
               if(allocated(max_v2s))deallocate(max_v2s)
               allocate(max_v2s(mc_ms))
               read(77,*)max_v2s
            case('MIN_V1S_CLASSES')
               if(allocated(min_v1s))deallocate(min_v1s)
               allocate(min_v1s(max(mc_gs,mc_fs)))
               read(77,*)min_v1s
            case('MIN_V2S_CLASSES')
               if(allocated(min_v2s))deallocate(min_v2s)
               allocate(min_v2s(mc_ms))
               read(77,*)min_v2s
            case('MIN_V2S_EXC_MAXES')
               if(allocated(min_v2s_excm))deallocate(min_v2s_excm)
               allocate(min_v2s_excm(mc_ms))
               read(77,*)min_v2s_excm
            case('MAX_V2S_PRESCREEN_CLASS')
               read(77,*)mc_ms_ps
            case('MAX_V2S_PRESCREEN')
               if(allocated(max_v2s_ps))deallocate(max_v2s_ps)
               allocate(max_v2s_ps(mc_ms_ps))
               read(77,*)max_v2s_ps
            case('MIN_V2S_PRESCREEN')
               if(allocated(min_v2s_ps))deallocate(min_v2s_ps)
               allocate(min_v2s_ps(mc_ms_ps))
               read(77,*)min_v2s_ps
            case('MIN_FREQ_MODE_V1')
               read(77,*)min_freq_mode
            case('MAX_FREQ_MODE_V1')
               read(77,*)max_freq_mode
            case('N_GROUND')
               read(77,*)n_gr
            case('TEMPERATURE')
               read(77,*)temp
            case('ANTI_STOKES')
               read(77,*)antistokes
            case('KBT_LIM')
               read(77,*)kbt_lim
            case('VERBOSE')
               read(77,*)verbose
            !case('V1_CEILING')
            !case('V2_CEILING')
         end select
      goto 1000
1001  close(77)
      if(.not.thr_v2_ov_found)thr_v2_ov=thresholdOrig_ex
      if(.not.thr_v2_comb_found)thr_v2_comb=thresholdOrig_ex
      if(.not.thr_v2_other_found)thr_v2_other=thresholdOrig_ex
      if(.not.thr_v2_0_found)thr_v2_0=thresholdOrig_ex
   end subroutine ReadOpt
      
   
   
   function TracelessQ6(q_6)result(q_tr)
      double precision q_6(6),q_tr(6)
      
      !Q
      !1  2  3  4  5  6
      !xx xy xz yy yz zz
      q_tr(1)=q_6(1)-0.5d0*(q_6(4)+q_6(6)) !xx
      q_tr(4)=q_6(4)-0.5d0*(q_6(1)+q_6(6)) !yy 
      q_tr(6)=q_6(6)-0.5d0*(q_6(1)+q_6(4)) !zz
      q_tr(2)=1.5d0*q_6(2) !xy,yx
      q_tr(3)=1.5d0*q_6(3) !xz,zx
      !q_tr(2,1)=q_tr(1,2)
      q_tr(5)=1.5d0*q_6(5) !yz,zy
      !q_tr(3,1)=q_tr(1,3)
      !q_tr(3,2)=q_tr(2,3)
   end function TracelessQ6
      
   !returns purely electronic contribution to polarizability
   !like polderdip maybe
   subroutine ElPol(td,nat,smat,nq,w_jn,wexc,nexc,gamma,u,v,m,q,du_c,dv_c,dm_c,dq_c,ht2,du2_c,dv2_c,dm2_c,dq2_c,vel,ground,grad,found_tms)
      type(excstate) td
      logical ground,vel,grad,found_tms,ht2
      integer nat,nexc,dl_idx,nq_orig,n3,iz,q_idx
      integer a,b,c,i,iexc,j,ii,jj,k,l
      integer,value :: nq
      integer,allocatable :: z_at(:)
      double precision w_jn,wexc(nexc),gamma,sqrt_w,wr
      double precision u(3),v(3),m(3),q(3,3),du_c(3,nat*3),dv_c(3,nat*3),dm_c(3,nat*3),dq_c(3,3,nat*3)
      double precision du2_c(:,:,:),dm2_c(:,:,:),dv2_c(:,:,:),dq2_c(:,:,:,:)
      double precision,allocatable :: grad_ex(:),grad_gr(:)
      double precision,allocatable :: ff_ex(:,:),ff_gr(:,:)
      double precision,allocatable :: ff_c_ex(:,:),ff_c_gr(:,:)
      double precision,allocatable :: smat(:,:),wg_orig(:)
      double complex :: f,f2
      double precision,allocatable :: du(:,:),dm(:,:),dv(:,:),dq(:,:,:),r(:)
      double precision,allocatable :: du2(:,:,:),dm2(:,:,:),dv2(:,:,:),dq2(:,:,:,:)
      type(Polar_exc),allocatable,target :: polars(:),polars2(:,:)
      type(Polar_exc),pointer :: pol
      double complex :: df,df2
            
      
      if(ground)then
         deallocate(smat)
         call readsi(n3,smat,wg_orig,nq,'ground/F.INP',z_at,r,.true.,iz)
         !smat=smat*0.0234280d0
         td%nq_gr=nq
         td%wg_gr=wg_orig
      ! else
         ! call readsi(n3,smat,,nq_orig,'excited/F.INP',.true.,iz)
      end if
      allocate(grad_ex(nq),grad_gr(nq))
      call readgradq(nq,grad_gr,'FILE.Q.GR.ground')
      call readgradq(nq,grad_ex,'FILE.Q.GR.excited.ground')
      if(ht2)then !fundamentals + combinations,overtones
         allocate(polars2(nq,nq))
         do i = 1,nq
            do j = i,nq
               call AllocatePolar(polars2(i,j),nexc)
               call Polar_new_nexc(polars2(i,j),nexc)
            end do
         end do
         allocate(ff_c_gr(n3,n3),ff_c_ex(n3,n3))
         call readff(nq,ff_c_gr,'ground/FILE.FC')
         call readff(nq,ff_c_ex,'excited/FILE.FC')
         ff_ex=Car2NM_FF(n3,nq,ff_c_ex,smat)
         ff_gr=Car2NM_FF(n3,nq,ff_c_gr,smat)
         deallocate(ff_c_gr,ff_c_ex)
      end if
      allocate(polars(nq))
      do i = 1,nq
         call AllocatePolar(polars(i),nexc)
         call Polar_new_nexc(polars(i),nexc)
      end do
      
      du=Car2NM(n3,nq,du_c,smat)
      dv=Car2NM(n3,nq,dv_c,smat)
      dm=Car2NM(n3,nq,dm_c,smat)
      dq=Car2NM_q(n3,nq,dq_c,smat)
      
      if(ht2)then
         du2=Car2NM_2(n3,nq,du2_c,smat)
         dv2=Car2NM_2(n3,nq,dv2_c,smat)
         dm2=Car2NM_2(n3,nq,dm2_c,smat)
         dq2=Car2NM_Q_2(n3,nq,dq2_c,smat)
         
         !force symmetry
         do i = 1,nq
            do j = i+1,nq
               du2(:,j,i)  =du2(:,i,j)
               dv2(:,j,i)  =dv2(:,i,j)
               dm2(:,j,i)  =dm2(:,i,j)
               dq2(:,:,j,i)=dq2(:,:,i,j)
            end do
         end do
      end if
      
      if(.not.found_tms)then
         call VelToLen_dip(v,dv,w_jn,nq,.false.,.false.,vel_grad,coordinates,.true.)
         call VelToLen_quad(q,dq,w_jn,nq,.false.,.false.,vel_grad,coordinates,.true.)
      end if
      if(vel)then
         u=v
         du=dv
         if(ht2)then
            du2=dv2
         end if
      end if
      
      if(ht2)then
         !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) &
         !$OMP SHARED(polars2,nq,nexc,wexc,wg_orig,u,m,q,du,dm,dq,du2,dm2,dq2,ff_gr,ff_ex,grad_gr,grad_ex,gamma,w_jn) &
         !$OMP PRIVATE(k,l,iexc,a,b,c,wr,pol,sqrt_w)
         do k = 1,nq
            do l = k,nq
               pol=>polars2(k,l)
               sqrt_w=1d0/sqrt(4d0*wg_orig(k)*wg_orig(l)*cm_2_au**2)*0.5d0 !0.5d0 is from Taylor expansion
               do iexc=1,nexc
                  wr=wexc(iexc)-(wg_orig(k)+wg_orig(l))*cm_2_au
                  do a = 1,3
                     do b = 1,3
                        pol%ap(a,b,iexc)=D2Pol(u(a),u(b),du(a,l),du(a,k),du(b,l),du(b,k),du2(a,l,k),du2(b,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wexc(iexc),gamma,.false.) &
                         + D2Pol(u(b),u(a),du(b,l),du(b,k),du(a,l),du(a,k),du2(b,l,k),du2(a,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wr,gamma,.true.)
                        pol%G(a,b,iexc)=D2Pol(u(a),m(b),du(a,l),du(a,k),dm(b,l),dm(b,k),du2(a,l,k),dm2(b,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wexc(iexc),gamma,.false.) &
                         + D2Pol(-m(b),u(a),-dm(b,l),-dm(b,k),du(a,l),du(a,k),-dm2(b,l,k),du2(a,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wr,gamma,.true.)
                        pol%Gc(a,b,iexc)=D2Pol(-m(a),u(b),-dm(a,l),-dm(a,k),du(b,l),du(b,k),-dm2(a,l,k),du2(b,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wexc(iexc),gamma,.false.) &
                         + D2Pol(u(b),m(a),du(b,l),du(b,k),dm(a,l),dm(a,k),du2(b,l,k),dm2(a,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wr,gamma,.true.)
                        do c = 1,3
                           pol%A(a,b,c,iexc)=D2Pol(u(a),q(b,c),du(a,l),du(a,k),dq(b,c,l),dq(b,c,k),du2(a,l,k),dq2(b,c,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wexc(iexc),gamma,.false.) &
                            + D2Pol(q(b,c),u(a),dq(b,c,l),dq(b,c,k),du(a,l),du(a,k),dq2(b,c,l,k),du2(a,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wr,gamma,.true.)
                           pol%Ac(a,b,c,iexc)=D2Pol(q(b,c),u(a),dq(b,c,l),dq(b,c,k),du(a,l),du(a,k),dq2(b,c,l,k),du2(a,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wexc(iexc),gamma,.false.) &
                            + D2Pol(u(a),q(b,c),du(a,l),du(a,k),dq(b,c,l),dq(b,c,k),du2(a,l,k),dq2(b,c,l,k),w_jn,grad_ex(l),grad_ex(k),grad_gr(l),grad_gr(k),ff_ex(l,k),ff_gr(l,k),wr,gamma,.true.)
                        end do
                     end do
                  end do
               end do
               
               pol%ap=conjg(pol%ap)*sqrt_w
               pol%G=pol%G*sqrt_w*iu
               pol%G=conjg(pol%G)
               pol%Gc=pol%Gc*sqrt_w*iu
               pol%Gc=conjg(pol%Gc)
               pol%A=conjg(pol%A*sqrt_w)
               pol%Ac=conjg(pol%Ac*sqrt_w)
            end do
         end do
         !$OMP END PARALLEL DO
      end if
      
      do k = 1,nq
         sqrt_w=1d0/sqrt(2d0*wg_orig(k)*cm_2_au)
         pol=>polars(k)
         do iexc=1,nexc
            wr=wexc(iexc)-wg_orig(k)*cm_2_au
            do a = 1,3
               do b = 1,3
                  pol%ap(a,b,iexc)=DPol(u(a),u(b),du(a,k),du(b,k),w_jn,grad_ex(k),grad_gr(k),wexc(iexc),gamma,.false.) &
                                  + DPol(u(b),u(a),du(b,k),du(a,k),w_jn,grad_ex(k),grad_gr(k),wr,gamma,.true.)
                  pol%G(a,b,iexc)=DPol(u(a),m(b),du(a,k),dm(b,k),w_jn,grad_ex(k),grad_gr(k),wexc(iexc),gamma,.false.) &
                                  + DPol(-m(b),u(a),-dm(b,k),du(a,k),w_jn,grad_ex(k),grad_gr(k),wr,gamma,.true.)
                  pol%Gc(a,b,iexc)=DPol(-m(a),u(b),-dm(a,k),du(b,k),w_jn,grad_ex(k),grad_gr(k),wexc(iexc),gamma,.false.) &
                                  + DPol(u(b),m(a),du(b,k),dm(a,k),w_jn,grad_ex(k),grad_gr(k),wr,gamma,.true.)
                  do c = 1,3
                     pol%A(a,b,c,iexc)=DPol(u(a),q(b,c),du(a,k),dq(b,c,k),w_jn,grad_ex(k),grad_gr(k),wexc(iexc),gamma,.false.) &
                                     + DPol(q(b,c),u(a),dq(b,c,k),du(a,k),w_jn,grad_ex(k),grad_gr(k),wr,gamma,.true.)
                     pol%Ac(a,b,c,iexc)=DPol(q(b,c),u(a),dq(b,c,k),du(a,k),w_jn,grad_ex(k),grad_gr(k),wexc(iexc),gamma,.false.) &
                                     + DPol(u(a),q(b,c),du(a,k),dq(b,c,k),w_jn,grad_ex(k),grad_gr(k),wr,gamma,.true.)
                  end do
               end do
            end do
         end do
         pol%ap=conjg(pol%ap)*sqrt_w
         pol%G=pol%G*sqrt_w*iu
         pol%G=conjg(pol%G)
         pol%Gc=pol%Gc*sqrt_w*iu
         pol%Gc=conjg(pol%Gc)
         pol%A=conjg(pol%A*sqrt_w)
         pol%Ac=conjg(pol%Ac*sqrt_w)
      end do
      
      deallocate(du,dv,dm,dq,wg_orig,grad_ex,grad_gr)
      do iexc = 1,nexc
         open(555+iexc,file=TR(dir_cwd)//'/FILE.'//TR(wexc_nm_str(iexc))//'.POLARS.EL')
         write(555+iexc,*)td%nq_gr,td%e_00*au_2_cm
         do ii = 1,td%nq_gr
            call WriteTenPretty(555+iexc,wexc(iexc),PolarExc2Polar(polars(ii),iexc),[0],[0_2],0,0d0,[1],[int(ii,kind=2)],1,td%wg_gr(ii))
         end do
         if(.not.ht2)then
            close(555+iexc)
            write(output_unit,*)'Written '//'FILE.'//TR(wexc_nm_str(iexc))//'.POLARS.EL'
         end if
      end do
      do i = 1,nq
         call deallocatePolar(polars(i))
      end do
         
      if(ht2)then
         deallocate(du2,dv2,dm2,dq2,ff_ex,ff_gr)
         do iexc = 1,nexc
            do ii = 1,td%nq_gr
               do jj = ii,td%nq_gr
                  if(ii==jj)then !overtone
                     call WriteTenPretty(555+iexc,wexc(iexc),PolarExc2Polar(polars2(ii,jj),iexc),[0],[0_2],0,0d0,[2],[int(ii,kind=2)],1,td%wg_gr(ii)+td%wg_gr(jj))
                  else !combination
                     call WriteTenPretty(555+iexc,wexc(iexc),PolarExc2Polar(polars2(ii,jj),iexc),[0],[0_2],0,0d0,[1,1],[int(ii,kind=2),int(jj,kind=2)],2,td%wg_gr(ii)+td%wg_gr(jj))
                  end if
               end do
            end do
            close(555+iexc)
            write(output_unit,*)'Written '//'FILE.'//TR(wexc_nm_str(iexc))//'.POLARS.EL'
         end do
         do i = 1,nq
            do j = i,nq
               call deallocatePolar(polars2(i,j))
            end do
         end do
         deallocate(polars2)
      end if
      
      deallocate(polars)
   end subroutine ElPol
   
   pure function DPol(a_a,b_b,da_a,db_b,wjn,dwj,dwn,w0,gamma,secondTerm)result(res)
      logical,intent(in) :: secondTerm
      double precision,intent(in) :: a_a,b_b,da_a,db_b,wjn,dwj,dwn,w0,gamma
      double complex f,df
      double complex res
      
      if(secondTerm)then !w0 is wR here
         f=wjn+w0+iu*gamma
      else
         f=wjn-w0-iu*gamma
      end if
      res=(da_a*b_b+a_a*db_b)/f + (a_a*b_b/f**2)*(dwj-dwn)
   end function DPol
   
   pure function D2Pol(a_a,b_b,da_a_l,da_a_k,db_b_l,db_b_k,da2_a,db2_b,wjn,dwj_l,dwj_k,dwn_l,dwn_k,d2wj,d2wn,w0,gamma,secondTerm)result(res)
      logical,intent(in) :: secondTerm
      double precision,intent(in) :: a_a,b_b,da_a_l,da_a_k,db_b_l,db_b_k,da2_a,db2_b,wjn,dwj_l,dwj_k,dwn_l,dwn_k,d2wj,d2wn,w0,gamma
      double complex f,df
      double complex res
      
      if(secondTerm)then !w0 is wR here
         f=wjn+w0+iu*gamma
      else
         f=wjn-w0-iu*gamma
      end if
      res=(da2_a*b_b+da_a_k*db_b_l+da_a_l*db_b_k+a_a*db2_b)/f + (da_a_k*b_b+a_a*db_b_k)/f**2*(dwj_l-dwn_l) &
      + (da_a_l*b_b+a_a*db_b_l)/f**2*(dwj_k-dwn_k)+a_a*b_b/f**3*(dwj_l-dwn_l)*(dwj_k-dwn_k)+a_a*b_b/f**2*(d2wj-d2wn)
   end function D2Pol
   
   
   
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
      double precision, allocatable :: dq_new(:,:,:)
      
      allocate(dq_new(3,3,nq))
      i_new=1
      do i = 1,nq
         !if(findloc(del_modes,i,dim=1)>0)cycle
         idx=(i_new-1)*a+1
         dq_new(:,:,i_new)=TracelessQ(dq(idx:idx+a-1))
         i_new=i_new+1
      end do
   end subroutine TidyDerivativesQ
   
   
   
   subroutine CW(str)
      character(*) str
      write(output_unit,*)str
      !write(output_unit,*)'------------------------------'
   end subroutine CW
   
   
   function MakeMap(del_modes,nq,nq_new)result(mapp)
      integer nq,nq_new,mapp(nq_new)
      integer del_modes(nq)
      integer i,idx
      
      idx=1
      do i = 1,nq
         if(findloc(del_modes,i,dim=1) > 0)cycle
         mapp(idx)=i
         idx=idx+1
      end do
   end Function MakeMap
   
     
   recursive subroutine makeCombinations(sequence,arr,n,k,i,ii)
      integer n,arr(n),k,j
      integer :: i
      integer :: sequence(k)
      integer :: ii
      
      if(i>k)then
         write(output_unit,*)sequence
         return
      end if
      
      do j = ii,n
         sequence(i)=arr(j)
         call makeCombinations(sequence,arr,n,k,i+1,j+1)
      end do
   end subroutine makeCombinations
   
   !taken from guvcde...I think
   subroutine moleculeDims(fo,nat)
      integer n,l,I,q,nmo,nmo2,ncmm,it,nt
      integer,intent(out) :: nat
      double precision axdum
      character*80 s80
      character*70 st
      character*(*) fo
      logical lzmat
      n=0
      nat=0
      ncmm=0
      it=0
      lzmat=.true. !no standard orientation for now
   
      open(2,file=fo)
1     read(2,2000,end=99,err=99)s80
2000  format(A80)

      IF((lzmat.and.(s80(19:39).EQ.'Z-Matrix orientation:'.OR. &
                     s80(26:46).EQ.'Z-Matrix orientation:'.OR. &
                     s80(20:37).EQ.'Input orientation:'.OR. &
                     s80(27:44).EQ.'Input orientation:')) &
          .OR. &
          ((.not.lzmat).and. &
                    (s80(20:40).EQ.'Standard orientation:'.OR. &
                     s80(26:46).EQ.'Standard orientation:')))THEN
       if(verbose)write(6,2000)s80
       DO I=1,4
          READ(2,*)
       enddo
       l=0
2005   READ(2,2000)s80
       IF(s80(2:4).NE.'---')THEN
        l=l+1
        BACKSPACE 2
        READ(2,*)q,q
        IF(q.EQ.-1)l=l-1
        GOTO 2005
       ENDIF
       nat=l
      ENDIF
      GOTO 1
99    close(2)
   end subroutine moleculeDims
   
   
      
   subroutine VelToLen_dip(v,dv,e00,nq,isAH,AHAS,vel_grad,coordinates,onlyDer)
      integer nq,i
      double precision v(3),dv(3,nq),e00,e00_copy
      double precision gr_q(nq),gr_q2(nq)
      logical isAH,AHAS,vel_grad,coordinates,onlyDer
      
      e00_copy=e00
      if(isAH)then !adiabatic
         ! open(77,file='ENERGY.VERTICAL')
         ! read(77,*)e00_copy
         ! close(77)
         if(AHAS)then
            call readgradQ(nq,gr_q,'FILE.Q.GR.excited')
            call readgradQ(nq,gr_q2,'FILE.Q.GR.ground.excited')
            if(.not.vel_grad)then
               gr_q=0
               gr_q2=0
            end if
            if(onlyDer)v=v*e00_copy
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
      if(onlyDer)v=v*e00_copy
      do i = 1,nq
         dv(:,i)=(v(:)/e00_copy**2)*gr_q(i)-dv(:,i)/e00_copy
      end do
      v=-v/e00_copy
   end subroutine VelToLen_dip
   
   subroutine VelToLen_quad(q,dq,e00,nq,isAH,AHAS,vel_grad,coordinates,onlyDer)
      integer nq,i
      double precision q(3,3),dq(3,3,nq),e00,e00_copy
      double precision gr_q(nq),gr_q2(nq)
      logical isAH,AHAS,vel_grad,coordinates,onlyDer
      
      !dq_help=dv
      e00_copy=e00
      if(isAH)then
         ! open(77,file='ENERGY.VERTICAL')
         ! read(77,*)e00_copy
         ! close(77)
         if(AHAS)then
            call readgradQ(nq,gr_q,'FILE.Q.GR.excited')
            call readgradQ(nq,gr_q2,'FILE.Q.GR.excited.ground')
            if(.not.vel_grad)then
               gr_q=0
               gr_q2=0
            end if
            if(onlyDer)q=q*e00_copy
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
      if(onlyDer)q=q*e00_copy
      do i = 1,nq
         dq(:,:,i)=(q(:,:)/e00_copy**2)*gr_q(i)-dq(:,:,i)/e00_copy
      end do
      q=-q/e00_copy
   end subroutine VelToLen_quad
   
   
   
   
   subroutine readDusch_Gauss(NQ,N3,A,B,C,D,E,J,K,filee)
      character(80) c80
      character(*) filee
      integer bufi,i,col,NQ,N3
      double precision,allocatable :: A(:,:),B(:),C(:,:),D(:),E(:,:),J(:,:),K(:),SG(:,:),SE(:,:),wg(:),we(:)
      integer,parameter :: unitt=77
      
      open(unitt,file=filee)
10    read(unitt,'(A80)')c80
      if(c80(1:18)==' Duschinsky matrix')THEN
         nq=n3-6
         goto 11
      end if
      if(c80(1:15)=='  Reduced system')THEN
         call Forward(unitt,2)
         nq=0
12       read(unitt,*)c80
         nq=nq+1
         if(index(c80,'=')>0)goto 12
         call Forward(unitt,2)
         goto 11
      end if
      goto 10
11    continue
      call Forward(unitt,3)
      deallocate(J,K)
      allocate(J(nq,nq),K(NQ))
      call readMatrix(unitt,J,NQ)
      call Forward(unitt,3)
      call readVector(unitt,K,NQ)
      
13    read(unitt,'(A80)')c80  
      if(c80(1:9)==' A Matrix')goto 14
      goto 13
14    continue
      deallocate(A,B,C,D,E)
      allocate(A(NQ,NQ),B(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ))
      call Forward(unitt,1)
      call readMatrix(unitt,A,NQ)
      
      call Forward(unitt,3)
      call readVector(unitt,B,NQ)
      
      call Forward(unitt,3)
      call readMatrix(unitt,C,NQ)
      
      call Forward(unitt,3)
      call readVector(unitt,D,NQ)
      
      call Forward(unitt,3)
      call readMatrix(unitt,E,NQ)
           
      ! SG=SG*0.0234280d0 !guvcde does this, i dunno why
      ! SE=SE*0.0234280d0
      close(unitt)
   end subroutine readDusch_Gauss
   
           
end program fcov
