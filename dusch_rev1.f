      module legacy
      use util
      implicit none
      
      contains
      
      subroutine reorder_J(N,NQ,MO,wg,we,sg,se,J_mat,K,lonly,swe)
         integer N,NQ,MO(NQ)
         double precision J_mat(NQ,NQ),K(NQ),sg(N,NQ),se(N,NQ),
     1    wg(NQ),we(NQ)
         logical swe,lonly
         
         integer i,jmin
         double precision dmin
         
       open(77,file='CORREL_MODES',action='WRITE')
        if(swe)then
         write(77,'(A2)')'EX'
        else
         write(77,'(A2)')'GR'
        end if
       write(6,*)'Correlating ground and excited state modes'
        i=NQ
c       loop starting from higher modes that migh be more reliable:
        do while(i>=1)
c       find best match for i within i...1:
        if(swe)then
           jmin=maxloc(J_mat(i,:)**2,dim=1)
           dmin=J_mat(i,jmin)**2
        else
           jmin=maxloc(J_mat(:,i)**2,dim=1)
           dmin=J_mat(jmin,i)**2
        end if
        
        if(i.ne.jmin)then
         if(.not.lonly)then
         if(swe)then
          call swi_w(we,i,jmin,nq)
          !call swi_w(ep_orig,i,jmin,nq)
          call swi_s(se,i,jmin,N)
          call swi_mat(J_mat,i,jmin,nq,.false.)
          !call swi_arrd(K,i,jmin,nq)
          !call swi_arrd(ge_q,i,jmin,nq)
         else
          call swi_w(wg,i,jmin,nq)
          call swi_s(sg,i,jmin,N)
          call swi_mat(J_mat,i,jmin,nq,.true.)
          call swi_arrd(K,i,jmin,nq)
          !call swi_arrd(ge_q,i,jmin,nq)
         endif
         if(swe)then
          write(6,619)i,jmin,dmin
          call swi(MO,i,jmin,NQ)
         else
          call swi(MO,i,jmin,NQ)
          write(6,620)i,jmin,dmin
         end if
         end if
         !ma(i)=.false.
         i=i-1
619      format(' Exc. mode',i4,' corresponds to',i4,' s =',f10.4,
     1   ' modes switched')
620      format(' Gr. mode',i4,' corresponds to',i4,' s =',f10.4,
     1   ' modes switched')
6091    format(' Exc. mode',i4,' corresponds to',i4,' s =',f10.4)
        ! write(77,*)i,jmin,' Y'
        else
         write(6,6091)i,jmin,dmin
         !ma(i)=.false.
         i=i-1
         end if
         end do
      end subroutine reorder_J
      
      subroutine ExtrapolateGradient(n3,H,re,rg,bohr,gr)
         implicit none
         integer n3
         double precision H(n3,n3),re(n3),rg(n3),gr(n3),bohr
         gr=matmul(H,(re-rg)/bohr)
      end subroutine ExtrapolateGradient
      
      function VerticalEnergyAH(n3,E_ex,E_gr,H,re,rg,bohr)result(en)
         implicit none
         integer n3
         double precision H(n3,n3),re(n3),rg(n3),gr(n3),E_ex,E_gr
         double precision en,mat_res,r_dif(n3),bohr
         
         r_dif=(re-rg)/bohr
         mat_res=0.5d0*dot_product(r_dif,matmul(H,r_dif))
         en=E_ex-E_gr+mat_res
         
      end function VerticalEnergyAH
            
      subroutine readDusch_Gauss(NQ,N3,A,B,C,D,E,J,K,filee)
         implicit none
         character(80) c80
         character(*) filee
         integer bufi,i,col,NQ,N3
         double precision :: A(NQ,NQ),B(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ),
     1 J(NQ,NQ),K(NQ)
         integer,parameter :: unitt=77
      
         open(unitt,file=filee)
10       read(unitt,'(A80)')c80
         if(c80(1:18)==' Duschinsky matrix')THEN
            nq=n3-6
            goto 11
         end if
         if(c80(1:15)=='  Reduced system')THEN
            call Forward(unitt,2)
            nq=0
12         read(unitt,*)c80
            nq=nq+1
            if(index(c80,'=')>0)goto 12
            call Forward(unitt,2)
            goto 11
         end if
         goto 10
11       continue
         call Forward(unitt,3)
         call readMatrix(unitt,J,NQ)
         call Forward(unitt,3)
         call readVector(unitt,K,NQ)
      
13       read(unitt,'(A80)')c80  
         if(c80(1:9)==' A Matrix')goto 14
         goto 13
14       continue
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
         close(unitt)
      end subroutine readDusch_Gauss
            
      
      subroutine CalcExpGeom(nat,N,nq,sg,M,K,r,fac)
         implicit none
         integer :: nat,n,nq
         double precision sg(N,N),M(nat),K(nq),r(N),fac
         double precision, parameter :: bohr = 0.529177d0
         integer :: i,ii,ix,j,l
         
         do i = 1,nat
            do ii = 1,3
               ix=3*(i-1)+ii
               do j = 1,NQ
                  r(ix)=r(ix)+sg(ix,j)/sqrt(1822.89d0)*
     1             K(j)*bohr*fac
               end do
            end do
         end do
      end SUBROUTINE CalcExpGeom
      
c     ==============================================================
      subroutine ortog(N,NQ,s,m)
c     orthogonalize s-vectors, to correct numerical inaccuracies:
      implicit none
      integer*4 N,NQ,i,j,ii,a
      real*8 s(N,NQ),m(N/3),sij,sn,dm
      real*8,allocatable::st(:,:)
      allocate(st(N,N))
      do 1 a=1,N/3
      dm=dsqrt(m(a))
      do 1 i=1,3
      ii=i+3*(a-1)
      do 1 j=1,NQ
 1    st(ii,j)=dm*s(ii,j)

      do 15 i=1,NQ
      do 16 j=1,i-1
      sij=0.0d0
      do 17 ii=1,N
17    sij=sij+st(ii,i)*st(ii,j)
      do 16 ii=1,N
16    st(ii,i)=st(ii,i)-sij*st(ii,j)
      sn=0.0d0
      do 18 ii=1,N
18    sn=sn+st(ii,i)**2
      sn=1.0d0/dsqrt(sn)
      do 15 ii=1,N
15    st(ii,i)=st(ii,i)*sn

      do 2 a=1,N/3
      dm=dsqrt(m(a))
      do 2 i=1,3
      ii=i+3*(a-1)
      do 2 j=1,NQ
 2    s(ii,j)=st(ii,j)/dm
      return
      end
c     ==============================================================
      subroutine cre(N,NQ,m,rg,re,s,y)
      implicit none
      integer*4 N,NQ,i,j,ia,ii
      real*8 rg(N),re(N),s(N,N),xn,t(3,3),m(N/3),dm
      real*8, allocatable::x(:),xv(:),v(:),st(:,:)
      character*1 y
      allocate(x(N),xv(NQ),v(N),st(N,N))
      if(y.eq.'e')then
       open(44,file='TMX.TXT')
       do 9 i=1,3
9      read(44,*)(t(j,i),j=1,3)
       close(44)
c      excited state S-matrix in ground-like orientation:
       do 14 ia=1,N/3
       dm=dsqrt(m(ia))
       do 14 j=1,NQ
       do 14 i=1,3
       ii=3*(ia-1)
14     st(i+ii,j)=(s(1+ii,j)*t(1,i)+s(2+ii,j)*t(2,i)+s(3+ii,j)*t(3,i))
     1 *dm
      else
       do 141 ia=1,N/3
       dm=dsqrt(m(ia))
       do 141 j=1,NQ
       do 141 i=1,3
141    st(i+3*(ia-1),j)=s(i+3*(ia-1),j)*dm
      endif

      x=re-rg
      xn=sp(x,x,N)
      write(6,600)xn
600   format(' Correcting excited geometry, old shift',f12.3)
      do 1 i=1,NQ
      do 2 j=1,N
2     v(j)=st(j,i)
1     xv(i)=sp(x,v,N)
      x=0.0d0
      do 3 i=1,NQ
      do 4 j=1,N
4     v(j)=st(j,i)
      do 3 j=1,N
3     x(j)=x(j)+xv(i)*v(j)
      xn=sp(x,x,N)
      write(6,601)xn
601   format('                              new shift',f12.3)
      re=rg+x
      return
      end
c     ==============================================================
      subroutine wdxx(s,n,z,r)
      implicit none
      character*(*) s
      integer*4 n,z(*),ia,ix
      real*8 r(*)
      open(9,file=s)
      write(9,*)'geometry'
      write(9,*)n
      do 1 ia=1,n
1     write(9,90)z(ia),(r(3*(ia-1)+ix),ix=1,3)
90    format(i3,3f12.6)
      close(9)
      return
      end
c     ==============================================================
      subroutine rdxx(s,ic,n,z,r)
      implicit none
      character*(*) s
      integer*4 ic,n,z(*),ia,ix
      real*8 r(*)
      open(9,file=s)
      read(9,*)
      read(9,*)n
      if(ic.ne.0)then
       do 1 ia=1,n
1      read(9,*)z(ia),(r(3*(ia-1)+ix),ix=1,3)
      endif
      close(9)
      return
      end
c     ==============================================================
      function sp(A,C,N)
c     scalar product of two vectors
      implicit none
      integer*4 N,k
      real*8 A(*),C(*),sp
      sp=0.0d0
      do 1 k=1,N
1     sp=sp+A(k)*C(k)
      return
      end
c     ==============================================================
      subroutine norm(v,n)
      implicit none
      integer*4 n
      real*8 v(n),vn
      vn=dsqrt(sp(v,v,n))
      v=v/vn
      return
      end
c     ==============================================================
      FUNCTION LOGDET(A,N,NQ,sgn)
c     logarithm of a determinant of a matrix
      IMPLICIT none
      integer*4 N,I,J,K,NQ,sgn
      logical DETEXISTS
      real*8 A(N,N),M,TEMP,LOGDET
      real*8, allocatable::ELEM(:,:)
      allocate(ELEM(NQ,NQ))
      sgn=1
      DO 1 I=1,NQ 
      DO 1 J=1,NQ 
1     ELEM(I,J)=A(I,J)
      DETEXISTS=.TRUE.
c     L=1.0d0
      DO 2 K=1,NQ-1
      IF(ABS(ELEM(K,K)).LE.1.0d-20)THEN
       DETEXISTS=.FALSE.
       DO 3 I=K+1,NQ
       IF(ELEM(I,K).NE.0.0d0)THEN
        DO 4 J=1,NQ
        TEMP=ELEM(I,J)
        ELEM(I,J)=ELEM(K,J)
4       ELEM(K,J)=TEMP
        DETEXISTS=.TRUE.
c       L=-L
        EXIT
       ENDIF
3      CONTINUE
       IF (DETEXISTS.EQV..FALSE.)THEN 
        LOGDET = 0.0d0
        RETURN
       ENDIF
      ENDIF
      DO 2 J=K+1,NQ
      M=ELEM(J,K)/ELEM(K,K)
      DO 2 I=K+1,NQ
2     ELEM(J,I)=ELEM(J,I)-M*ELEM(K,I)
      LOGDET =0.0d0
      DO 5 I=1,NQ
      if(ELEM(I,I).lt.0.0d0)then
c      L=-L
       LOGDET=LOGDET+log(-ELEM(I,I))
       sgn=sgn*-1
      else
       LOGDET=LOGDET+log( ELEM(I,I))
      endif
5     continue
      RETURN
      END
c     ==============================================================
      subroutine mv(A,B,C,N)
c     multiplication of matrix B by vector C, A = B . C
      implicit none
      integer*4 N,i,k
      real*8 A(*),B(N,N),C(*)
      do 1 i=1,N
      A(i)=0.0d0
      do 1 k=1,N
1     A(i)=A(i)+B(i,k)*C(k)
      return
      end
c     ==============================================================
      subroutine mm(A,B,C,N)
c     matrix multiplication A = B x C
      implicit none
      integer*4 N,i,j,k
      real*8 A(N,N),B(N,N),C(N,N)
      do 1 i=1,N
      do 1 j=1,N
      A(i,j)=0.0d0
      do 1 k=1,N
1     A(i,j)=A(i,j)+B(i,k)*C(k,j)
      return
      end
c     ==============================================================
      SUBROUTINE wdv(io,N,V,s)
      IMPLICIT none
      integer*4 N,io,LN
      real*8 V(*)
      character*(*) s
      call wlabel(io,s)
      write(io,501)1
501   format(i14)
      DO 130 LN=1,N
130   write(io,500)LN,V(LN)
500   format(i7,1x,D14.6)
      return
      end
c     ==============================================================
      subroutine ms(A,B,C,N)
c     matrix sum A = B + C
      implicit none
      integer*4 N,i,j
      real*8 A(N,N),B(N,N),C(N,N)
      do 1 i=1,N
      do 1 j=1,N
1     A(i,j)=B(i,j)+C(i,j)
      return
      end
c     ==============================================================
      subroutine mt(AT,A,N)
c     matrix transposition AT = A^T
      implicit none
      integer*4 N,i,j
      real*8 A(N,N),AT(N,N)
      do 1 i=1,N
      do 1 j=1,N
1     AT(i,j)=A(j,i)
      return
      end
c     ==============================================================
      subroutine vz(A,N)
      implicit none
      integer*4 n,i
      real*8 A(*)
      do 3 i=1,N
3     A(i)=0.0d0
      return
      end
c     ==============================================================
      subroutine mz(A,N)
      implicit none
      integer*4 n,i,j
      real*8 A(N,N)
      do 3 i=1,N
      do 3 j=1,N
3     A(i,j)=0.0d0
      return
      end
c     ==============================================================
      subroutine trse(N,s,NQ,t)
      implicit none
      integer*4 N,NQ,i,a,ia,ix
      real*8 s(N,NQ),t(3,3),v(3)
      do 1 i=1,NQ
      do 1 ia=1,N/3
      a=3*(ia-1)
      do 2 ix=1,3
2     v(ix)=s(ix+a,i)
      do 1 ix=1,3
1     s(ix+a,i)=t(1,ix)*v(1)+t(2,ix)*v(2)+t(3,ix)*v(3)
      return
      end
c     ==============================================================
      subroutine wlabel(io,s)
      implicit none
      character*(*) s
      integer*4 J,io
      write(io,*)
      write(io,101)
      do 1 J=1,len(s)
1     write(io,103)s(J:J)
      write(io,*)
      write(io,101)
101   format(' ',$)
      do 2 J=1,len(s)
2     write(io,103)'-'
103   format(a1,$)
      write(io,*)
      return
      end
c     ==============================================================
      SUBROUTINE wwm(io,N,e,s)
      IMPLICIT none
      integer*4 io,N,i
      real*8 e(*),CM
      character*(*) s
      CM=219474.630d0
      call wlabel(io,s)
      write(io,100)(e(i)*CM,i=1,N)
100   format(6f11.2)
      return
      end
c     ==============================================================
      SUBROUTINE wds(io,N0,N,M,s)
      IMPLICIT none
      integer*4 N,N1,N3,io,LN,J,N0
      real*8 M(N0,N0)
      character*(*) s
      call wlabel(io,s)
      N1=1
1     N3=min(N1+4,N)
      write(io,100)(J,J=N1,N3)
100   format(10x,5i14)
      DO 130 LN=1,N0
130   write(io,200)LN,(M(LN,J),J=N1,N3)
200   format(i7,1x,5e14.6)
      N1=N1+5
      IF(N3.LT.N)GOTO 1
      return
      end
c     ==============================================================
      SUBROUTINE wdm(io,N,M,s)
      IMPLICIT none
      integer*4 N,N1,N3,io,LN,J
      real*8 M(N,N)
      character*(*) s
      call wlabel(io,s)
      N1=1
1     N3=min(N1+4,N)
      write(io,100)(J,J=N1,N3)
100   format(10x,5i14)
      DO 130 LN=1,N
130   write(io,200)LN,(M(LN,J),J=N1,N3)
200   format(i7,1x,5e14.6)
      N1=N1+5
      IF(N3.LT.N)GOTO 1
      return
      end
c     ==============================================================
      subroutine report(s)
      character*(*) s
      write(6,*)s
      stop
      end
c     ==============================================================
C       subroutine seq(arr,N)
C          implicit none
C          integer N,arr(N),i
C          do i = 1,N
C             arr(i)=i
C          end do
C       end subroutine seq
      
      subroutine swi(arr,i,j,N)
         implicit none
         integer N,arr(N),i,j,buf
         
         buf=arr(i)
         arr(i)=arr(j)
         arr(j)=buf
      end subroutine swi
      
      subroutine swi_arrd(arr,i,j,N)
         implicit none
         integer N,i,j
         double precision arr(N),buf
         
         buf=arr(i)
         arr(i)=arr(j)
         arr(j)=buf
      end subroutine swi_arrd
      
      subroutine swi_mat(mat,i,j,n,rows)
         implicit none
         integer i,j,n,k
         double precision mat(n,n),buf
         logical rows
         
         if(rows)then
            do k = 1,N
               buf=mat(i,k)
               mat(i,k)=mat(j,k)
               mat(j,k)=buf
            end do
         else
            do k = 1,N
               buf=mat(k,i)
               mat(k,i)=mat(k,j)
               mat(k,j)=buf
            end do
         end if
      end subroutine swi_mat
      
      subroutine swi_s(S,i,j,N)
         implicit none
         integer N,i,j,k
         double precision S(N,N),buf
         
         do k = 1,N
            buf=S(k,i)
            S(k,i)=S(k,j)
            S(k,j)=buf
         end do
      end subroutine swi_s
      
      subroutine swi_w(w,i,j,N)
         implicit none
         integer N,i,j,k
         double precision w(N),buf
         
            buf=w(i)
            w(i)=w(j)
            w(j)=buf
      end subroutine swi_w
      
      !insert column i to column j and shift all elements to right
      subroutine ins(vec,n,i,j)
         integer i,j,n,ii
         double precision vec(n),buf,buf2
         
         buf=vec(j)
         vec(j)=vec(i)
         
         do ii = j+1,n-1
            buf2=vec(ii)
            vec(ii)=buf
            buf=buf2 !wait wtf, fix this TODO
         end do
      end subroutine ins
      
c     ==============================================================
      subroutine readm(m)
      IMPLICIT none
      real*8 m(*)
      logical lex
      character*3 key
      integer*4 i
      inquire(file='DUSCH.OPT',exist=lex)
      if(lex)then
       open(11,file='DUSCH.OPT')
1      read(11,900,end=11,err=11)key
 900   format(a3)
       if(key.eq.'MAS')read(11,*)i,m(i)
       goto 1
11     close(11)
      endif
      return
      end
c     ==============================================================
      subroutine readopt(lit,maxit,lwr,wmin,lfix,lmscan,NB,NP,TOL,
     1lspec,del,e00,wsmin,wsmax,nps,lshifte,lshifti,lshiftq,ncut,
     1LEXCL,isu,llin,ifx,wlim,lcor,lonly,lsfix,lcre,kfac,lcom,wfix,
     1ldz,sslim,swe,vert_h,j_one,we_eq_wg,se_eq_sg,lnoreor,vert_h_cor,
     1 vert_h_w,delk,fix6,expGeom,expGeom_fac,Gauss_dusch,
     1 Gauss_dusch_file,isplanar,J_ad,deltol,reg,corenVH)
      IMPLICIT none
      character*80 Gauss_dusch_file
      integer*4 maxit,NB,NP,nps,ncut,LEXCL,isu,ifx,nblo,delk,expGeom,
     1 Gauss_dusch,corenVH
      real*8 wmin,TOL,e00,wsmin,wsmax,del,wlim,kfac,wfix,sslim,
     1 expGeom_fac,deltol
      logical lit,lwr,lex,lfix,lmscan,lspec,lshifte,lshifti,lshiftq,
     1llin,lcor,lonly,lsfix,lcre,lcom,ldz,swe,vert_h,j_one,we_eq_wg,
     1 se_eq_sg,lnoreor,vert_h_cor,vert_h_w,fix6,isplanar,J_ad,reg
      character*13 key
c     delete zero modes:
      ldz=.true.
c     replace negative modes by this value if required (in cm-1):
      wfix=100.0d0
c     ensure that Q Q' and K are in the same space
      lcom=.false.
c     factor for K shift
      kfac=1.0d0
      delk=0
c     fix final geometry according to available modes:
      lcre=.false.
      deltol=0.1d0
      reg=.false. 
c
c     fix se matrix for deleted modes:
      lsfix=.false.
      lnoreor=.false.
c     
c     ilinear shift = excited modes same as ground
      llin=.false.
c
c     type of FC sum
      isu=2
c     maximal excitation for statesum1:
      LEXCL=3
c     if ncut>0, modify excited q definition according to the shift
      ncut=0
c     various paths to get the normal modeshift:
      lshifte=.false.
      lshifti=.false.
      lshiftq=.false.
c     make a spectrum:
      lspec=.false.
c     spectral parameters:e00,wmin,wmax (all cm-1) and number of points:
      e00=0.0d0
      wsmin=-1000.0d0
      wsmax= 1000.0d0
      nps=4001
c     spectral bandwidth:
      del=10.0d0
c     reorder excited modes according to ground:
      lcor=.false.
      swe=.false.
      j_one=.false.
      we_eq_wg=.false.
      se_eq_sg=.false.
      lnoreor=.false.
      vert_h=.false.
      vert_h_w=.false.
      vert_h_cor=.false.
      corenVH=2
c     to CORREL, only show best match, no action:
      lonly=.false.
c     use fixed geometries, do not alling:
      lfix=.false.
c     scan modes -NP...NP:
      NP=10
c     quantum number biffer size:
      NB=10
c     take quantum number cut off:
      TOL=0.1d0
c     tolerance for J-metrix element
      sslim=0.9d0
      lit=.false.
      wmin=1.0d0
      maxit=100
      lwr=.true.
      lmscan=.false.
      fix6=.false.
      isplanar=.false.
      J_ad=.false.
      expGeom=0
      Gauss_dusch=0
      Gauss_dusch_file='no'
c     what to do with negative modes:
c     ifx = 1 make them +100 cm-1
c           2 delete
c           3 borrow ground or excited state frequency
c           4 delete if w < wlim
      ifx=1
c     limit for mode deletion cm-1:
      wlim=600.0d0
      inquire(file='DUSCH.OPT',exist=lex)
      if(lex)then
       open(11,file='DUSCH.OPT')
1      read(11,900,end=11,err=11)key
 900   format(a13)
       if(key(1:3).eq.'BLO')then
        write(6,*)' Bloccked modes will be read again in fixmodes'
        read(11,*)nblo
        read(11,*)
        if(nblo.lt.0)read(11,*)
       endif
       if(key(1:6).eq.'COR_VH')then
         read(11,900)key
         read(key,'(L)')vert_h_cor
       end if
       if(key(1:3).eq.'REG')read(11,*)reg
       if(key(1:6).eq.'PLANAR')read(11,*)isplanar
       if(key(1:3).eq.'LDZ')read(11,*)ldz
       if(key(1:5).eq.'K_DEL')read(11,'(I8)')delk
       if(key(1:4).eq.'COR ')read(11,*)lcor
       if(key(1:3).eq.'SWE')read(11,*)swe
       if(key(1:3).eq.'DEL')read(11,*)del
       if(key(1:3).eq.'DLT')read(11,*)deltol
       if(key(1:3).eq.'E00')read(11,*)e00
       if(key(1:4).eq.'FIX ')read(11,'(I4)')ifx
       if(key(1:3).eq.'ISU')read(11,*)isu
       if(key(1:3).eq.'ITE')read(11,*)lit
       if(key(1:3).eq.'MAX')read(11,*)maxit
       if(key(1:2).eq.'NB')read(11,*)NB
       if(key(1:4).eq.'NCUT')read(11,*)ncut
       if(key(1:2).eq.'NP')read(11,*)NP
       if(key(1:3).eq.'NSP')read(11,*)nps
       if(key(1:5).eq.'LEXCL')read(11,*)LEXCL
       if(key(1:4).eq.'KFAC')read(11,*)kfac
       if(key(1:4).eq.'LCOM')read(11,*)lcom
       if(key(1:4).eq.'WFIX')read(11,*)wfix
       if(key(1:3).eq.'LFI')read(11,*)lfix
       if(key(1:5).eq.'LINEA')read(11,*)llin
       if(key(1:3).eq.'LWR')read(11,*)lwr
       if(key(1:3).eq.'ONL')read(11,*)lonly
       if(key(1:3).eq.'SCA')read(11,*)lmscan
       if(key(1:5).eq.'LSFIX')read(11,*)lsfix
       if(key(1:4).eq.'LCRE')read(11,*)lcre
       if(key(1:6).eq.'SHIFTQ')read(11,*)lshiftq
       if(key(1:6).eq.'SHIFTI')read(11,*)lshifti
       if(key(1:6).eq.'SHIFTE')read(11,*)lshifte
       if(key(1:5).eq.'SPECT')read(11,*)lspec
       if(key(1:3).eq.'TOL')read(11,*)tol
       if(key(1:3).eq.'WLI')read(11,*)wlim
       if(key(1:5).eq.'SSLIM')read(11,*)sslim
       if(key(1:3).eq.'WMI')read(11,*)wmin
       if(key(1:5).eq.'WSMIN')read(11,*)wsmin
       if(key(1:5).eq.'WSMAX')read(11,*)wsmax
       if(key(1:6).eq.'VERT_H')read(11,'(L)')vert_h
       if(key(1:4).eq.'J_AD')read(11,'(L)')J_ad
       if(key(1:5).eq.'J_ONE')read(11,'(L)')j_one
       if(key(1:8).eq.'WE_IS_WG')read(11,'(L)')we_eq_wg
       if(key(1:8).eq.'SE_IS_SG')read(11,'(L)')se_eq_sg
       if(key(1:8).eq.'NO_REORI')read(11,'(L)')lnoreor
       if(key(1:8).eq.'VH_COREN')read(11,'(I1)')corenVH
       if(key(1:4).eq.'W_VH')read(11,'(L)')vert_h_w
       if(key(1:4).eq.'FIX6')read(11,*)fix6
       if(key(1:9).eq.'EXP_GEOM ')read(11,*)expGeom
       if(key(1:12).eq.'EXP_GEOM_FAC')read(11,*)expGeom_fac
       if(key(1:8) .eq.'DUSCH_G ')read(11,*)Gauss_dusch
       if(key(1:12) .eq.'DUSCH_G_FILE')read(11,*)Gauss_dusch_file
       goto 1
11     close(11)
      endif
      return
      end
c     ==============================================================
      subroutine fixmodes(N,NQ,NE,g,sg,p,se,amas,z,ifx,wlim,lcor,lonly,
     1 sgf,sef,NF,wfix,m,nat,sslim,swe,del_modes_gr,del_modes_ex,MO,
     1 vert_h,j_one,we_eq_wg,fix6,rg,re,delK,J_mat,K,iz,J_adiabatic,
     1 ep_orig,deltol)
      IMPLICIT none
      integer*4 N,NE,NQ,I,iz,ifx,j,ix,no,jmin,nblo,z(*),NU,NF,ii,
     1 imx,nat,del_modes_gr(NQ),del_modes_ex(NQ),MO(NQ),delk,idx
      real*8 CM,wlim,dij,dmin,t,sgf(N,NQ),
     1 sef(N,NQ),wfix,m(*),
     1 ss,sslim,sg_m(N,NQ),ge_q(NQ),ffq(NQ,NQ),m_invsq(N,N),
     1 rg(N),re(N),J_mat(nq,nq),K(nq),ffq_i(NQ,NQ),ep_orig(NQ)
      real*8,allocatable :: g(:),p(:),sg(:,:),se(:,:)
      real*8 amas(*)
      real*8 tol,deltol
      real*8,parameter :: amu2au = 1822.89d0, bohr=0.529177d0
      real*8, allocatable :: ws(:)
      integer*4,allocatable::bl(:),be(:)
      character*3 key
      logical lex,lcor,lonly,swe,vert_h,j_one,we_eq_wg,fix6,J_adiabatic
      logical,allocatable :: ma(:)
      integer,allocatable :: coup1(:),coup2(:)
      integer :: idx1,idx2
      real*8,allocatable::JT(:,:)
      CM=219474.630d0
      tol=1d-10
      
      no=NQ
      
      write(6,6009)ifx,NQ
6009  format(/,' Fixmodes',/,' IFIX = ',I4,/,' NQ = ',I4)

      if(NQ.ne.NE)then
       NU=min(NE,NQ)
       write(6,60091)NU
60091  format(' NE<>NQ, take minimum ',i4,' modes')
       NE=NU
       NQ=NU
      endif
      
      if(vert_h)then
         call makeSM(nat,N,nq,sg,m,sg_m,iz)
         allocate(ws(NQ))
         !allocate(J_mat(NQ,NQ),K(NQ),ws(NQ))
         J_mat=0d0
         K=0d0
         call mkj_vert(tol,amu2au,nat,N,NQ,m,sg_m,m_invsq,ge_q,ffq,
     1    ffq_i,J_mat,K,ws,j_one,we_eq_wg,g)
         !if(J_one)call eye(J_mat,nq)
         if(J_adiabatic)then
            J_mat=0d0
            call mkj(nat,N,NQ,J_mat,sg,se,m)
            J_mat=transpose(J_mat)
         end if
      else
         !allocate(J_mat(nq,nq),K(nq))
         J_mat=0d0
         !if(J_one)then
         !   call eye(J_mat,nq)
         !else
         call mkj(nat,N,NQ,J_mat,sg,se,m)
         !end if
         J_mat=transpose(J_mat)
         call mkk(bohr,nat,N,Nq,m,rg,re,sg,K)
      end if
      
      
      if(lcor)then
        open(77,file='CORREL_MODES',action='WRITE')
        if(swe)then
         write(77,'(A2)')'EX'
        else
         write(77,'(A2)')'GR'
        end if
       write(6,*)'Correlating ground and excited state modes'
       if(lonly)then
        write(6,*)'Correlation only, no action'
        do 121 i=NQ,1,-1
        call findbest(i,NQ,jmin,dmin,dij,sg,se,N,amas,z)
121     write(6,609)i,jmin,dmin
609     format(' Ground mode',i4,' corresponds to',i4,' s =',f10.4)
        write(6,*)
        do 122 i=NQ,1,-1
        call findbest(i,NQ,jmin,dmin,dij,se,sg,N,amas,z)
122     write(6,6091)i,jmin,dmin
6091    format(' Exc. mode',i4,' corresponds to',i4,' s =',f10.4)
       else
        allocate(ma(nq))
        i=NQ
        ma=.true.
c       loop starting from higher modes that migh be more reliable:
        do while(i>=1)
c       find best match for i within i...1:
        if(swe)then
           jmin=maxloc(J_mat(i,:)**2,dim=1,mask=ma)
           dmin=J_mat(i,jmin)**2
        else
           jmin=maxloc(J_mat(:,i)**2,dim=1,mask=ma)
           dmin=J_mat(jmin,i)**2
        end if
        
        if(i.ne.jmin)then
         if(swe)then
          call swi_w(p,i,jmin,nq)
          call swi_w(ep_orig,i,jmin,nq)
          call swi_s(se,i,jmin,N)
          call swi_mat(J_mat,i,jmin,nq,.false.)
          !call swi_arrd(K,i,jmin,nq)
          !call swi_arrd(ge_q,i,jmin,nq)
         else
          call swi_w(g,i,jmin,nq)
          call swi_s(sg,i,jmin,N)
          call swi_mat(J_mat,i,jmin,nq,.true.)
          call swi_arrd(K,i,jmin,nq)
          call swi_arrd(ge_q,i,jmin,nq)
         endif
         if(swe)then
          write(6,619)i,jmin,dmin
          !call swi(MO,i,jmin,NQ)
         else
          call swi(MO,i,jmin,NQ)
          write(6,620)i,jmin,dmin
         end if
         !ma(i)=.false.
         i=i-1
619      format(' Exc. mode',i4,' corresponds to',i4,' s =',f10.4,
     1   ' modes switched')
620      format(' Gr. mode',i4,' corresponds to',i4,' s =',f10.4,
     1   ' modes switched')
        ! write(77,*)i,jmin,' Y'
         
        else
         write(6,6091)i,jmin,dmin
         !ma(i)=.false.
         i=i-1
        ! write(77,*)i,jmin,' N'
        endif
        end do
       endif
        deallocate(ma)
        do i = 1,NQ
         write(77,*)i,MO(i)
        end do
        close(77)
      endif
            
c     ifx = 0   do nothing
c           1   change negative frequencies to wfix
c          -1   change negative frequencies to positive
c           2   delete negatives
c           3   borrow exc. freqs from ground and vice versa
c           4   delete smaller than a limit
c           5   delete smaller than a limit in ground and then in excited state
c          >5   modes which have sum of squares of J matrix elements smaller than sslim are removed (ie. modes are incomplete after deletion)
c     how many modes deleted:
      NF=0
      
C       allocate(coup1(nq),coup2(nq))
C       do i = 1,nq
C          coup1=0
C          coup2=0
C          coup1(1)=i
C          idx1=2
C          idx2=1
C          call GetCoupledModes_Box(i,idx1,idx2,nq,.false.,J_mat,0.1d0,
C      1    .true.,coup1,coup2)
C       end do
C       deallocate(coup1,coup2)
      
      
      open(77,file='IMAG_MODES_EX')
      do i=1,nq
         if(ep_orig(i)<0d0)then
            write(77,*)i
         end if
      end do
      close(77)
      
      open(77,file='IMAG_MODES_GR')
      do i=1,nq
         if(g(i)<0d0)then
            write(77,*)i
         end if
      end do
      close(77)
      
      if(abs(ifx).eq.1)then
       iz=0
       do 2 I=1,NQ
       if(g(I).lt.0.0d0)then
        if(ifx==-1)then
         g(I)=-g(I)
        else
         g(I)=wfix/CM
        end if
        iz=iz+1
       endif
       if(p(I).lt.0.0d0)then
        if(ifx==-1)then
         p(I)=-p(I)
        else
         p(I)=wfix/CM
        end if
        iz=iz+1
       endif
 2     continue
       if(iz.gt.0)then
         if(ifx==-1)then
            write(6,45)iz
         else
            write(6,46)iz,wfix
         end if
       end if
45     format(i6,' neg. modes made positive')       
46     format(i6,' neg. modes made',f12.2,' cm-1 ')
      endif

      if(ifx.eq.2)then
       open(789,file='DELETED_MODES_GR')
       allocate(coup1(nq),coup2(nq))
       iz=0
66     do 3 I=1,NQ 
       if(g(I).lt.0.0d0)then
        coup1=0
        coup1(1)=i
        coup2=0
        idx1=2
        idx2=1
        if(deltol==1)then
         coup2(1)=i
         idx2=2
        else
        call GetCoupledModes_Box(i,idx1,idx2,nq,.true.,J_mat,deltol,
     1  .true.,coup1,coup2,nf)
        end if
        if(idx1/=idx2)stop 33 !I assume that the Duschinsky matrix has to be square by some vague and unexplained physical definitions
        nf=nf+idx1-1
        do ii=1,nq
         if(coup1(ii)>0)then
            call DelMode_one(nat,n,nq,coup1(ii),.true.,
     1       J_mat,K,sg,se,g,p,MO)
         end if
         if(coup2(ii)>0)then
            call DelMode_one(nat,n,nq,coup2(ii),.false.,
     1       J_mat,K,sg,se,g,p,MO)
         end if
        end do
        cycle
c       save deleted sg vector:
        NF=NF+1
        del_modes_gr(nf)=i
        do 41 ix=1,N
41      sgf(ix,NF)=sg(ix,i)
        do 4 j=i,NQ-1
        g(j)=g(j+1)
        J_mat(j,:)=J_mat(j+1,:)
        K(j)=K(j+1)
        do 4 ix=1,N
4       sg(ix,j)=sg(ix,j+1)
        jmin=maxloc(J_mat(i,:)**2,dim=1)
!        jmin=i
        del_modes_ex(nf)=jmin
        do ix=1,N
         sef(ix,NF)=se(ix,jmin)
        end do
        do 44 j=jmin,NQ-1
        p(j)=p(j+1)
        J_mat(:,j)=J_mat(:,j+1)
        do 44 ix=1,N
44      se(ix,j)=se(ix,j+1)
        iz=iz+1
        NQ=NQ-1
        goto 66 !is this cycle or does it repeat the loop? sigh...
       endif
 3     continue
       
666    do 333 I=1,NQ
       if(p(i).lt.0.0d0)then
        coup1=0
        coup1(1)=i
        coup2=0
        idx1=2
        idx2=1
        if(deltol==1)then
         coup2(1)=i
         idx2=2
        else
        call GetCoupledModes_Box(i,idx1,idx2,nq,.false.,J_mat,deltol,
     1  .true.,coup1,coup2,nf)
        end if
        if(idx1/=idx2)stop 33
        nf=nf+idx1-1
        do ii=1,nq
         if(coup1(ii)>0)then
            write(789,*)coup1(ii)
            call DelMode_one(nat,n,nq,coup1(ii),.true.,
     1       J_mat,K,sg,se,g,p,MO)
         end if
         if(coup2(ii)>0)then
            call DelMode_one(nat,n,nq,coup2(ii),.false.,
     1       J_mat,K,sg,se,g,p,MO)
         end if
        end do
        cycle
c       save deleted sg vector:
        NF=NF+1
        del_modes_ex(nf)=i
        do 411 ix=1,N
411      sef(ix,NF)=se(ix,i)
        do 444 j=i,NQ-1
        p(j)=p(j+1)
        J_mat(:,j)=J_mat(:,j+1)
        do 444 ix=1,N
444       se(ix,j)=se(ix,j+1)
        jmin=maxloc(J_mat(:,i)**2,dim=1)
        !jmin=i
        do ix=1,N
         sgf(ix,NF)=sg(ix,jmin)
        end do
        del_modes_gr(nf)=jmin
        do 4444 j=jmin,NQ-1
        g(j)=g(j+1)
        J_mat(j,:)=J_mat(j+1,:)
        K(j)=K(j+1)
        do 4444 ix=1,N
4444    sg(ix,j)=sg(ix,j+1)
        iz=iz+1
        NQ=NQ-1
        goto 666 !is this cycle or does it repeat the loop? sigh...
       endif
 333   continue
       if(nf.gt.0)write(6,*)nf,' coupled modes to negative mode deleted'
       deallocate(coup1)
       close(789)
       if(fix6)ifx=7 !!!
      endif

      if(abs(ifx).eq.3)then
       iz=0
       do 5 I=1,NQ
       if(p(i).lt.0d0)then
        !idx=i
        if(g(i).gt.0d0)then !ground state freq. is a-ok
         if(ifx==-3)then
            p(i)=0d0
            do ii = 1,nq
               p(i)=p(i)+J_mat(ii,i)**2*g(ii)
            end do
         else
            p(i)=g(i)
         end if
         iz=iz+1
        else if(g(i).lt.0d0)then !ground state frequency is negative (like in water clusters)
         g(i)=-g(i)
         p(i)=g(i)
         iz=iz+1
        else
         write(6,*)I
         call report(' mode cannot be fixed')
        endif
       endif
       if(g(i).lt.0.0d0)then
        idx=i
        if(p(idx).gt.0.0d0)then
         if(ifx==-3)then
            g(i)=0d0
            do ii = 1,nq
               g(i)=g(i)+J_mat(i,ii)**2*p(ii)
            end do
         else 
            g(i)=p(idx)
         end if
         iz=iz+1
        else
         write(6,*)I
         call report(' mode cannot be fixed')
        endif
       else if(g(i).lt.0d0)then !ground state frequency is negative (like in water clusters)
        g(i)=-g(i)
       endif
5      continue
       if(iz.gt.0)write(6,*)iz,' neg. modes borrowed'
      endif

      if(abs(ifx).eq.4)then
       iz=0
67     do 6 i=1,NQ
       if(g(i)*CM.lt.wlim.or.p(i)*CM.lt.wlim .or.
     1  abs(g(i)*CM).lt.wlim.or.abs(p(i)*CM).lt.wlim.and.ifx==-4)then
c       save deleted sg vector:
        NF=NF+1
        do 71 ix=1,N
        sef(ix,nf)=se(ix,i)
71      sgf(ix,NF)=sg(ix,i)
        call DelMode(nat,n,no,i,i,g(i)*CM.lt.wlim,J_mat,K,sg,se,g,p)
        !call ddm(N,i,NQ,g,sg)
        !call ddm(N,i,NQ,p,se)
        NQ=NQ-1
        iz=iz+1
        goto 67
       endif
 6     continue
       if(iz.gt.0)write(6,*)iz,' small modes skipped'
       if(fix6)ifx=6
      endif

      if(ifx.eq.5)then

c      delete ground small modes
       iz=0
87     do 8 i=1,NQ
       if(g(i)*CM.lt.wlim)then
c       make Duschinsky J-matrix
        allocate(JT(NQ,NQ))
        call mkj(nat,N,NQ,JT,sg,se,m)
c       find corresponding excited
        imx=fce(i,NQ,JT)
        write(6,601)'  ground',i,g(i)*CM,imx,p(imx)*CM,' excited'
        call ddm(N,i,NQ,g,sg)
        call ddm(N,imx,NQ,p,se)
        iz=iz+1
        NQ=NQ-1
        nf=nf+1
        deallocate(JT)
        goto 87
       endif
 8     continue

c      delete excited small modes
107    do 10 i=1,NQ
       if(p(i)*CM.lt.wlim)then
        allocate(JT(NQ,NQ))
        call mkj(nat,N,NQ,JT,sg,se,m)
        imx=fcg(i,NQ,JT)
        write(6,601)' excited',i,p(i)*CM,imx,g(imx)*CM,'  ground'
601     format(a8,i4,f7.1,' deleted, with #',i5,f7.1,a8)
        call ddm(N,i,NQ,p,se)
        call ddm(N,imx,NQ,g,sg)
        iz=iz+1
        NQ=NQ-1
        nf=nf+1 !however sgf not saved
        deallocate(JT)
        goto 107
       endif
 10    continue


       if(iz.gt.0)write(6,*)iz,' small modes skipped'
       if(NQ.lt.1)call report('Nothing left')

      endif
      
      if(ifx.gt.5)then
c       check that big overlap exists
222     allocate(JT(NQ,NQ))
        call mkj(nat,N,NQ,JT,sg,se,m)
        do 166 i=1,NQ
        ss=0.0d0
        do 13 j=1,NQ
13      ss=ss+JT(j,i)**2
        if(ss.lt.sslim)then
         imx=fce(i,NQ,JT)
         write(6,602)'  ground',i,g(i)*CM,ss,imx,p(imx)*CM
602      format(a8,i5,f7.1,f10.2,' eliminated with',i5,f7.1)
         call ddm(N,i,NQ,g,sg)
         call ddm(N,imx,NQ,p,se)
         iz=iz+1
         NQ=NQ-1
         nf=nf+1 !however sgf not saved
         deallocate(JT)
         if(NQ.gt.0)goto 222
        endif
166      continue

        do 14 i=1,NQ
        ss=0.0d0
        do 17 j=1,NQ
17      ss=ss+JT(i,j)**2
        if(ss.lt.sslim)then
         imx=fcg(i,NQ,JT)
         write(6,602)' excited',i,p(i)*CM,ss,imx,g(imx)*CM
         call ddm(N,i,NQ,p,se)
         call ddm(N,imx,NQ,g,sg)
         iz=iz+1
         NQ=NQ-1
         nf=nf+1 !however sgf not saved
         deallocate(JT)
         if(NQ.gt.0)goto 222
        endif
14      continue
      endif
      
      if(delk>0 .and. sum(abs(K))>0.001)then
         do i = 1,delK
            idx=maxloc(abs(K),dim=1)
C             ii=maxloc(J_mat(idx,:)**2,dim=1)
            ii=idx
            call DelMode(nat,n,nq,idx,ii,.true.,J_mat,K,sg,se,g,p)
         end do
         NQ=NQ-delK
         NF=NF+delK
      end if
      if(NQ.lt.1)call report('Nothing left')
      
      nblo=0
      inquire(file='DUSCH.OPT',exist=lex)
      if(lex)then
       open(11,file='DUSCH.OPT')
1      read(11,900,end=11,err=11)key
 900   format(a3)
       if(key.eq.'BLO')then
        read(11,*)nblo
        allocate(bl(abs(nblo)),be(abs(nblo)))
        read(11,*)bl
        if(nblo.lt.0)read(11,*)be
       endif
       goto 1
11     close(11)
      endif

      if(nblo.ne.0)then
       call dm(N,NQ,abs(nblo),bl,sg,g)
       if(nblo.gt.0)then
c       blocked same modes in excited:
        call dm(N,NQ,abs(nblo),bl,se,p)
       else
c       special list for excited:
        call dm(N,NQ,abs(nblo),be,se,p)
       endif
       NQ=NQ-abs(nblo)
       if(NQ.ne.no)write(6,*)NQ,' modes now'
      endif
      
C       if(nf>0 .and. vert_h)then
         
C          call makeSM(nat,N,nq,sg,m,sg_m,iz)
C          deallocate(ws)
C          allocate(ws(NQ))
C          !allocate(J_mat(NQ,NQ),K(NQ),ws(NQ))
C          J_mat=0d0
C          K=0d0
C          call mkj_vert(tol,amu2au,nat,N,NQ,m,sg_m,m_invsq,ge_q,ffq,
C      1    ffq_i,J_mat,K,ws,j_one,we_eq_wg,g)
C          if(J_one)call eye(J_mat,nq)
C       end if
C       K=0
      if(J_one)call eye(J_mat,nq)
      if(vert_h)then
         if(we_eq_wg)then
            ws=g
         else
            ws=p
         end if
         call mkk_vert(NQ,ws,J_mat,ge_q,K)
         !call mkk_vert_ffq(NQ,ffq,ge_q,K)
      end if
      J_mat=transpose(J_mat)
      open(792,file='VEC.K')
      do i = 1,nq
         write(792,'(E16.8E2)')K(i)
      end do
      close(792)
      open(792,file='KEPT_MODES_GR')
      do i = 1,nq-nf
         write(792,*)MO(i)
      end do
      close(792)
      return
      end
      

      
      !Removes modes which have low J element
      !TODO
      subroutine RemoveDiffuseModes(nat,n,nq,tol,J,K,sg,se,wg,we,
     1 del_gr,del_ex)
         integer nat,n,nq,del_gr(nq),del_ex(nq),mode_c,mode_l(NQ,NQ),
     1   mode_l_i ,i,ii
         double precision tol,J(nq,nq),K(nq),sg(N,N),se(N,N),
     1    wg(N),we(N)
         
         
         !i dont think it matters whether you start with excited or ground modes
         !so im just going to follow cache-friendly approach
         do i = 1,NQ
            mode_c=0
            mode_l(:,i)=0
            mode_l_i=1
            do ii = 1,NQ
               if(J(ii,i)**2>tol)then
                  mode_c=mode_c+1
                  mode_l(mode_l_i,i)=ii
                  mode_l_i=mode_l_i+1
               end if
            end do
         end do
      end subroutine RemoveDiffuseModes
      
      !invHes is the inverted Hessian of either excited state or ground state
      !gr is the excited state gradient
      subroutine ExtrapolateGeometry(r,N,invHes,gr,fac)
         implicit none
         integer n
         double precision r(N),invHes(N,N),gr(N),fac
         r=r-fac*matmul(invHes,gr)*0.529177249
      end subroutine ExtrapolateGeometry
      
      !tol should be given like 0.1
      subroutine GetCoupledModes(i,nq,isGr,J,tol,useabs,coup,coup_n)
         integer i,nq,coup(nq),ii,idx,coup_n
         double precision j(nq,nq),tol
         logical isGr,useabs
         
         
         coup=0
         idx=1
         if(isGr)then
            do ii = 1,NQ
               if(useabs .and. abs(J(i,ii))>tol)then
                  coup(idx)=ii
                  idx=idx+1
               else if((J(i,ii))**2>tol)then
                  coup(idx)=ii
                  idx=idx+1
               end if
            end do
         else
            do ii = 1,NQ
               if(useabs .and. abs(J(ii,i))>tol)then
                  coup(idx)=ii
                  idx=idx+1
               else if((J(i,ii))**2>tol)then
                  coup(idx)=ii
                  idx=idx+1
               end if
            end do
         end if
         coup_n=count(coup>0)
      end subroutine GetCoupledModes
      
      recursive subroutine GetCoupledModes_Box(i,idx1,idx2,
     1 nq,isGr,J,tol,useabs,coup1,coup2,deletedC)
         integer, value :: i
         integer ii,nq,coup1(nq),coup2(nq),idx1,idx2,deletedC
         double precision j(nq,nq),tol
         logical isGr,useabs,iscoupled,atleastone
         
         
         atleastone=.false.
         do ii = 1,nq-deletedC
            iscoupled=DeterAddToCoup(coup2,idx2,i,ii,nq,useabs,isGr,
     1          tol,J)
            atleastone=atleastone.or.iscoupled
            if(iscoupled)then
               call GetCoupledModes_Box(ii,idx2,idx1,nq,.not.isGr,J,tol,
     1          useabs,coup2,coup1,deletedC)
            end if
         end do
         if(.not.atleastone)return
      end subroutine GetCoupledModes_Box
      
      function DeterAddToCoup(coup,idx,fixI,ii,nq,useabs,isGr,tol,J)
     1 result(res)
         logical useabs,res,isGr
         integer nq,idx,fixI,ii,coup(nq)
         double precision J(nq,nq),tol,vall
         
         res=.false.
         if(findloc(coup,ii,dim=1)>0)return
         if(isGr)then
            vall=J(fixI,ii)
         else
            vall=J(ii,fixI)
         end if
         if(useabs .and. abs(vall)>tol)then
            coup(idx)=ii
            idx=idx+1
            res=.true.
         else if((vall)**2>tol)then
            coup(idx)=ii
            idx=idx+1
            res=.true.
         end if
      end function DeterAddToCoup
      
      !deletes 1 mode (i) and then 1 mode (i2). IisGround determines if i is an index of the ground state or the excited state
      subroutine DelMode(nat,N,NQ,i,i2,IisGround,J,K,sg,se,wg,we)
         integer nat,n,nq,i,ii,i2
         double precision J(Nq,nq),k(NQ),sg(N,N),se(N,N),wg(N),we(N)
         logical IisGround
         
         do ii = i,NQ-1
            if(IisGround)then
               wg(ii)=wg(ii+1)
               sg(:,ii)=sg(:,ii+1)
               J(ii,:)=J(ii+1,:)
               K(ii)=K(ii+1)
            else
               we(ii)=we(ii+1)
               se(:,ii)=se(:,ii+1)
               J(:,ii)=J(:,ii+1)
            end if
         end do
         
         do ii = i2,NQ-1
            if(IisGround)then
               we(ii)=we(ii+1)
               se(:,ii)=se(:,ii+1)
               J(:,ii)=J(:,ii+1)
            else
               wg(ii)=wg(ii+1)
               sg(:,ii)=sg(:,ii+1)
               J(ii,:)=J(ii+1,:)
               K(ii)=K(ii+1)
            end if
         end do
         
         
      end subroutine DelMode

      !deletes 1 mode (i), IisGround determines if i is an index of the ground state or the excited state
      subroutine DelMode_one(nat,N,NQ,i,IisGround,J,K,sg,se,wg,we,MO)
         integer nat,n,nq,i,ii,MO(nq)
         double precision J(Nq,nq),k(NQ),sg(N,N),se(N,N),wg(N),we(N)
         logical IisGround
         
         do ii = i,NQ-1
            if(IisGround)then
               wg(ii)=wg(ii+1)
               sg(:,ii)=sg(:,ii+1)
               J(ii,:)=J(ii+1,:)
               K(ii)=K(ii+1)
               MO(ii)=MO(ii+1)
            else
               we(ii)=we(ii+1)
               se(:,ii)=se(:,ii+1)
               J(:,ii)=J(:,ii+1)
            end if
         end do
         
      end subroutine DelMode_one


      subroutine dm(N,NQO,nblo,b,s,e)
      implicit none
      integer*4 nblo,N,NQ,i,j,ix,b(nblo),NQO
      real*8 s(N,N),e(N)
      NQ=NQO
      do 1 i=1,NQ
      do 1 j=1,nblo
1     if(i.eq.b(j))e(i)=-3.333d0
68    do 2 i=1,NQ
      if(e(i).eq.-3.333d0)then
       do 3 j=i,NQ-1
       e(j)=e(j+1)
       do 3 ix=1,N
3      s(ix,j)=s(ix,j+1)
       NQ=NQ-1
       goto 68
      endif
2     continue
      write(6,*)nblo,' modes blocked (deleted)'
      return
      end

      SUBROUTINE readsi_dusch(N3,S,E,NQ,fn,r,ldz,iz)
c     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,iz,fuck1,fuck2
      real*8 CM,r(:)
      real*8,allocatable :: S(:,:), S_help(:,:)
      real*8,allocatable :: E(:),E_help(:)
      logical ldz
      character*(*) fn
      CM=219474.630d0
C     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn)
      read(4,*)NQ,nat,nat
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
c
      write(6,*)NQ,' modes found'
c     delete zero modes if exist:
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
         write(6,*)iz,' zero modes: deleted'
         allocate(S_help(N3,NQ),E_help(NQ))
         S_help=S(:,1:NQ)
         E_help=E(1:NQ)
         deallocate(S,E)
         S=S_help
         E=E_help
         deallocate(S_help,E_help)
       endif 
      endif
      
      
      write(6,*)NQ,' vibrational modes considered'
      
      DO 3 I=1,NQ
3     E(I)=E(I)/CM
      
      RETURN
      end
c     ==============================================================
      SUBROUTINE writesi_dusch(NAT,N3,S,E,NQ,fn,x,z)
c     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,z(*)
      real*8 S(N3,N3),E(*),CM,x(*)
      character*(*) fn
      CM=219474.630d0
C     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn)
      write(4,*)NQ,nat,nat
      do 1 i=1,NAT
1     write(4,401)z(i),(x(ix+3*(i-1)),ix=1,3)
401   format(i6,3f12.6)
      write(4,*)'Atom mode x y z'
      DO 2 I=1,NAT
      DO 2 J=1,NQ
2     write(4,402)I,J,(s(3*(i-1)+ix,NQ-J+1),ix=1,3)
402   format(2i6,3f12.6)
      write(4,402)NQ
      write(4,4000)(E(NQ-J+1)*CM,J=1,NQ)
4000  FORMAT(6F11.3)
      close(4)
      RETURN
      end
c     ==============================================================
      subroutine xst(nat,rgo,re,m,z,t,isplanar)
c     standard orientation of re with respect to rg
c     J. Phys. Chem. A 2001, 105, 5326-5333
      implicit none
      logical isplanar
      integer*4 nat,i,j,ia,IERR,ix,iy,iz,ii,ib,ic,z(*)
      real*8 rgo(*),re(*),c(3,3),m(*),cc(3,3),r(3,3),l(3),ci(3,3),
     1TOL,e(3,6),lrc(3,3),bl(3),rc(3,3),err,t0(3,3),t(3,3),em,
     1cmg(3),cme(3),mass
      real*8,allocatable::rn(:),rg(:)
      TOL=1.0d-10
      allocate(rg(3*nat))
c     calculate mass centers:
      do 14 i=1,3
      cmg(i)=0.0d0
      cme(i)=0.0d0
      mass=0.0d0
      do 16 ia=1,nat
      mass=mass+m(ia)
      cmg(i)=cmg(i)+m(ia)*rgo(i+3*(ia-1))
16    cme(i)=cme(i)+m(ia)* re(i+3*(ia-1))
      cmg(i)=cmg(i)/mass
14    cme(i)=cme(i)/mass
      write(6,603)cmg,cme
603   format(' Mass centers: g: ',3f12.4,/,
     1       '               e: ',3f12.4)
c     mass-center coordinates:
      do 15 ia=1,nat
      do 15 i=1,3
      rg(i+3*(ia-1))=rgo(i+3*(ia-1))-cmg(i)
15    re(i+3*(ia-1))= re(i+3*(ia-1))-cme(i)
c     form C =X M X':
      do 1 i=1,3
      do 1 j=1,3
      c(i,j)=0.0d0
      do 1 ia=1,nat
1     c(i,j)=c(i,j)+rg(i+3*(ia-1))*m(ia)*re(j+3*(ia-1))
c     correction for planar molecules:
      do 2 i=1,3
2     if(isplanar)c(i,i)=1.0d0
c     C^TC:
      call mm3(c,c,cc)
      CALL TRED12(3,cc,r,l,2,IERR)
      IF(IERR.NE.0)call report(' cannot diagonalize C^TC')
      call INV(3,c,ci,3,TOL,e,IERR)
      IF(IERR.NE.0)call report(' cannot invert C')
c     R^TC^-1:
      call mm3(r,ci,rc)
c     temporary coordinate string:
      allocate(rn(3*nat))
c     exlore all eight axis switching possibilities:
      ib=0
      ic=0
      do 4 ix=-1,1,2
      do 4 iy=-1,1,2
      do 4 iz=-1,1,2
      ic=ic+1
      bl(1)=dble(ix)
      bl(2)=dble(iy)
      bl(3)=dble(iz)
c     L l^(1/2)R^TC^-1:
      do 6 i=1,3
      do 6 j=1,3
6     lrc(i,j)=bl(i)*dsqrt(l(i))*rc(i,j)
      call md(r,lrc,t0)
      call trx(nat,rn,re,t0)
      err=0.0d0
      do 8 ii=1,3*nat
8     err=err+(rn(ii)-rg(ii))**2
      if(ic.eq.1)then
       em=err
       do 9 i=1,3
       do 9 j=1,3
9      t(i,j)=t0(i,j)
       ib=ic
      endif
      if(det(t0).gt.0.0d0.and.err.lt.em)then
       em=err
       do 11 i=1,3
       do 11 j=1,3
11     t(i,j)=t0(i,j)
       ib=ic
      endif
4     write(6,600)ix,iy,iz,err,det(t0)
600   format(3i2,'err:',f10.3,' det:',f10.3)
      write(6,601)ib
601   format(' Best trial number ',i2)
      write(6,604)t
604   format(' Transformation matrix:',/,3(3f10.3,/))
      open(44,file='TMX.TXT')
      do 17 i=1,3
17    write(44,605)(t(j,i),j=1,3)
605   format(3f10.3)
      close(44)
      call trx(nat,rn,re,t)
      do 131 ia=1,nat
      do 131 i=1,3
131   re(i+3*(ia-1))=cmg(i)+rn(i+3*(ia-1))
      call wrx(nat,re,z,'exc.as.ground.x',
     1'Excited orientation to match the ground:')
      write(6,*)'exc.as.ground.x written'
      return
      end
c     ==============================================================
      subroutine wrx(nat,r,z,s,sp)
      implicit none
      integer*4 nat,z(nat),a,i
      real*8 r(3*nat)
      character*(*) s,sp
      open(44,file=s)
      write(44,'(a)')sp
      write(44,*)nat
      do 13 a=1,nat
13    write(44,602)z(a),(r(i+3*(a-1)),i=1,3)
602   format(i4,3f12.6)
      close(44)
      return
      end
c     ==============================================================
      subroutine trx(nat,n,e,t)
      implicit none
      integer*4 nat,ia,i,a
      real*8 n(*),e(*),t(3,3)
      do 7 ia=1,nat
      a=3*(ia-1)
      do 7 i=1,3
7     n(i+a)=t(1,i)*e(1+a)+t(2,i)*e(2+a)+t(3,i)*e(3+a)
      return
      end
c     ==============================================================
      function det(t)
      implicit none
      real*8 det,t(3,3)
      det=t(1,1)*(t(2,2)*t(3,3)-t(2,3)*t(3,2))
     1   -t(1,2)*(t(2,1)*t(3,3)-t(2,3)*t(3,1))
     1   +t(1,3)*(t(2,1)*t(3,2)-t(2,2)*t(3,1))
      return
      end
c     ==============================================================
      subroutine md(a,b,c)
      implicit none
c     c=a.b
      real*8 a(3,3),b(3,3),c(3,3)
      integer*4 i,j,ii
      do 3 i=1,3
      do 3 j=1,3
      c(i,j)=0.0d0
      do 3 ii=1,3
3     c(i,j)=c(i,j)+a(i,ii)*b(ii,j)
      return
      end
c     ==============================================================
      subroutine mm3(a,b,c)
      implicit none
c     c=aT.b
      real*8 a(3,3),b(3,3),c(3,3)
      integer*4 i,j,ii
      do 3 i=1,3
      do 3 j=1,3
      c(i,j)=0.0d0
      do 3 ii=1,3
3     c(i,j)=c(i,j)+a(ii,i)*b(ii,j)
      return
      end
c     ==============================================================
      SUBROUTINE TRED12(N,A,Z,D,IEIG,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N),Z(N,N),D(*)
      real*8,allocatable::E(:)
      real*8 h
      allocate(E(N))
      IF (IEIG .LT. 0) GO TO 110
      DO 100 I = 1, N
C
      DO 100 J = 1, I
      Z(I,J) = A(I,J)
  100 CONTINUE
C
  110 ICODE = IABS(IEIG)
      IF (N .EQ. 1) GO TO 320
C     :::::::::: FOR I=N STEP -1 UNTIL 2 DO -- ::::::::::
      DO 300 II = 2, N
      I = N + 2 - II
      L = I - 1
      H = 0.0D0
      SCALE = 0.0D0
      IF (L .LT. 2) GO TO 130
C     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
      DO 120 K = 1, L
  120 SCALE = SCALE + DABS(Z(I,K))
C
      IF (DABS(SCALE) .GT. 1.0D-10) GO TO 140
  130 E(I) = Z(I,L)
      GO TO 290
C
  140 DO 150 K = 1, L
      F = Z(I,K) / SCALE
      Z(I,K) = F
      H = H + F * F
  150 CONTINUE
C
      G = -DSIGN(DSQRT(H),F)
      E(I) = SCALE * G
      H = H - F * G
      Z(I,L) = F - G
      F = 0.0D0
C
      DO 240 J = 1, L
      IF (ICODE .EQ. 2) Z(J,I) = Z(I,J) / (SCALE * H)
      G = 0.0D0
C     :::::::::: FORM ELEMENT OF A*U ::::::::::
      DO 180 K = 1, J
  180 G = G + Z(J,K) * Z(I,K)
C
      JP1 = J + 1
      IF (L .LT. JP1) GO TO 220
C
      DO 200 K = JP1, L
  200 G = G + Z(K,J) * Z(I,K)
C     :::::::::: FORM ELEMENT OF P ::::::::::
  220 E(J) = G / H
      F = F + E(J) * Z(I,J)
  240 CONTINUE
C
      HH = F / (H + H)
C     :::::::::: FORM REDUCED A ::::::::::
      DO 260 J = 1, L
      F = Z(I,J)
      G = E(J) - HH * F
      E(J) = G
C
      DO 260 K = 1, J
      Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260 CONTINUE
C
      DO 280 K = 1, L
  280 Z(I,K) = SCALE * Z(I,K)
C
  290 D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0D0
      E(1) = 0.0D0
C     :::::::::: ACCUMULATION OF TRANSFORMATION MATRICES ::::::::::
      IF (ICODE .NE. 2) GO TO 600
      DO 500 I = 1, N
      L = I - 1
      IF (DABS(D(I)) .LT. 1.0D-10) GO TO 380
C
      DO 360 J = 1, L
      G = 0.0D0
C
      DO 340 K = 1, L
  340 G = G + Z(I,K) * Z(K,J)
C
      DO 360 K = 1, L
      Z(K,J) = Z(K,J) - G * Z(K,I)
  360 CONTINUE
C
  380 D(I) = Z(I,I)
      Z(I,I) = 1.0D0
      IF (L .LT. 1) GO TO 500
C
      DO 400 J = 1, L
      Z(I,J) = 0.0D0
      Z(J,I) = 0.0D0
  400 CONTINUE
C
  500 CONTINUE
C
  620 CALL TQL12 (N,Z,D,IERR,ICODE,E)
      RETURN
C     :::::::::: ALTERNATE FINAL LOOP FOR EIGENVALUES ONLY ::::::::::
  600 DO 610 I=1,N
  610 D(I) = Z(I,I)
      GO TO 620
      END
C
      SUBROUTINE TQL12(N,Z,D,IERR,ICODE,E)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION Z(N,N),D(*),E(*)
      REAL*8 MACHEP,H
      EPS = 1.0D0
   10 EPS = 0.50D0*EPS
      TOL1 = EPS + 1.0D0
      IF((TOL1.GT.1.0D0).AND.(TOL1-EPS.EQ.1.0D0)) GO TO 10
      IF(TOL1-EPS.EQ.1.0D0)  EPS = EPS + EPS
      MACHEP = EPS
C
C     MACHEP=16.0D0**(-13)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
      J = 0
      H = MACHEP * (DABS(D(L)) + DABS(E(L)))
      IF (B .LT. H) B = H
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
      DO 110 M = L, N
      IF (DABS(E(M)) .LE. B) GO TO 120
C     :::::::::: E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ::::::::::
  110 CONTINUE
C
  120 IF (M .EQ. L) GO TO 220
  130 IF (J .EQ. 30) GO TO 1000
      J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
      L1 = L + 1
      G = D(L)
      P = (D(L1) - G) / (2.0D0 * E(L))
      R = DSQRT(P*P+1.0D0)
      D(L) = E(L) / (P + DSIGN(R,P))
      H = G - D(L)
C
      DO 140 I = L1, N
  140 D(I) = D(I) - H
C
      F = F + H
C     :::::::::: QL TRANSFORMATION ::::::::::
      P = D(M)
      C = 1.0D0
      S = 0.0D0
      MML = M - L
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
      DO 200 II = 1, MML
      I = M - II
      G = C * E(I)
      H = C * P
      IF (DABS(P) .LT. DABS(E(I))) GO TO 150
      C = E(I) / P
      R = DSQRT(C*C+1.0D0)
      E(I+1) = S * P * R
      S = C / R
      C = 1.0D0 / R
      GO TO 160
  150 C = P / E(I)
      R = DSQRT(C*C+1.0D0)
      E(I+1) = S * E(I) * R
      S = 1.0D0 / R
      C = C * S
  160 P = C * D(I) - S * G
      D(I+1) = H + S * (C * G + S * D(I))
C     :::::::::: FORM VECTOR ::::::::::
      IF (ICODE .NE. 2) GO TO 200
      DO 180 K = 1, N
      H = Z(K,I+1)
      Z(K,I+1) = S * Z(K,I) + C * H
      Z(K,I) = C * Z(K,I) - S * H
  180 CONTINUE
C
  200 CONTINUE
C
      E(L) = S * P
      D(L) = C * P
      IF (DABS(E(L)) .GT. B) GO TO 130
  220 D(L) = D(L) + F
  240 CONTINUE
C     :::::::::: ORDER EIGENVALUES AND EIGENVECTORS ::::::::::
      DO 300 II = 2, N
      I = II - 1
      K = I
      P = D(I)
C
      DO 260 J = II, N
      IF (D(J) .LE. P) GO TO 260
      K = J
      P = D(J)
  260 CONTINUE
C
      IF (K .EQ. I) GO TO 300
      D(K) = D(I)
      D(I) = P
C
      IF (ICODE .NE. 2) GO TO 300
      DO 280 J = 1, N
      P = Z(J,I)
      Z(J,I) = Z(J,K)
      Z(J,K) = P
  280 CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
 1000 IERR = L
 1001 RETURN
      END
C
      SUBROUTINE INV(nr,a,ai,n,TOL,e,IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      dimension ai(nr,nr),a(nr,nr),e(nr,2*nr)
C
10000 IERR=0
      DO 1  ii=1,n
      DO 1  jj=1,n
      e(ii,jj)=a(ii,jj)
1     e(ii,jj+n)=0.0D0
      do 13 ii=1,n
13    e(ii,ii+n)=1.0D0
c
      DO 2 ii=1,n-1
      if (ABS(e(ii,ii)).LE.TOL) then
       DO 3 io=ii+1,n
3      if (ABS(e(io,ii)).GT.TOL) goto 11
       IERR=1
       write(6,*)ii
       write(6,*)'tol = ',tol
       tol=tol*0.50d0
       if(tol.gt.1.0d-20)goto 10000
       RETURN
c
11     CONTINUE
       DO 4 kk=1,2*n
       w=e(ii,kk)
       e(ii,kk)=e(io,kk)
4      e(io,kk)=w
      ENDIF
      eii=e(ii,ii)
      DO 5 jj=ii+1,n
      e1=e(jj,ii)/eii
      DO 6 kk=ii+1, 2*n
6     e(jj,kk)=e(jj,kk)-e(ii,kk)*e1
5     e(jj,ii)=0.0D0
2     CONTINUE
c
      DO 7  i2=n,2,-1
      eii=e(i2,i2)
      DO 7 j2=i2-1,1,-1
      e1=e(j2,i2)/eii
      DO 9 kk=1, n
9     e(j2,kk+n)=e(j2,kk+n)-e(i2,kk+n)*e1
7     e(j2,i2)=0.0d0
c
      DO 10 ii=1,n
      ei=1.0d0/e(ii,ii)
      DO 12  jj=1,n
12    ai(ii,jj)=e(ii,jj+n)*ei
10    CONTINUE
c
      RETURN
      END
c     ============================================================
      subroutine makeG(NQ,G,GP,GH,GHP,eg,ep)
      implicit none
      real*8 eg(*),ep(*),G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ)
      integer*4 ix,NQ
      call mz(G,NQ)
      call mz(GP,NQ)
      call mz(GH,NQ)
      call mz(GHP,NQ)
      do 10 ix=1,NQ
      G(ix,ix)=eg(ix)
      GP(ix,ix)=ep(ix)
      GH(ix,ix) =dsqrt(abs(eg(ix)))
10    GHP(ix,ix)=dsqrt(abs(ep(ix)))
      return
      end
      
      !zero modes last, imaginary modes before them, then the rest
      subroutine repairH(W,J,N)
         integer N,i
         real*8 J(N,N),W(N)
         
         do i = N,1,-1
            if(W(i)<1d-10)THEN
               
            end if
         end do
      end subroutine repairH
c     ============================================================
      function overlap(ic,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac,
     1 vert_h,j_one,ifx,wfix,cor,swe,MO,sg_full,se_full,we_eq_wg,
     1 vert_h_w,Gauss_dusch,Gauss_dusch_file)
c     calculates <0|0'>
c     ic  ... 0   only the overlap
c             1   also write DUSCH.OUT
      implicit none
      integer*4 ic,NQ,N,ig,ia,nat,ix,jx,ii,IERR,i,i2,imag_c,
     1 j2,jj,ifx,idx,MO(NQ),nqn,kk,dels,Gauss_dusch
      real*8 eg(NQ),ep(NQ),G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ),
     1 J(NQ,NQ),JT(NQ,NQ),sg(N,NQ),se(N,NQ),m(nat),k(NQ),rg(N),re(N),
     2 bohr,T(NQ,NQ),TP(NQ,NQ),F(NQ,NQ),FI(NQ,NQ),E2(NQ,2*NQ),TOL,
     3 A(NQ,NQ),TV(*),TU(*),B(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ),overlap,
     4 fac,kfac,toll,wfix,sg_full(N,N),se_full(N,N),d2,curfreq,CFN,
     5 J2_mean,J2_std,dw_mean,dw_std
      character(*) Gauss_dusch_file
      logical :: vert_h,j_one,cor,swe,we_eq_wg,vert_h_w
      real*8, allocatable ::
     5 grad_gr(:),grad_exc(:),grad_gr_q(:),grad_exc_q(:),
     6 ffc(:,:),ffq(:,:),ffq_i(:,:),W_VH2(:),sg_m(:,:),se_m(:,:),
     6 help(:,:),
     7 m_invsq(:,:),ffc_a(:,:),W_VH2_inv(:,:),K_VH(:),J_VH(:,:),
     8 test(:,:),ws(:),js(:,:),ks(:)
      real*8,parameter :: amu2au=1822.89d0
      real*8,parameter :: au2cm=1d0/4.556431403d-6
      
      call makeG(NQ,G,GP,GH,GHP,eg,ep)
C       allocate(sg_m(N,N),se_m(N,N))
C       sg_m=0
C       se_m=0
C       do i = 1,nat
C          ix=(i-1)*3+1 
C          do i2 = 0,2 !x y z
C             do jj = 1,NQ
C                !jx=(jj-1)*3+1
C                !do j2 = 0,2 !x y z
C                   sg_m(ix+i2,jj)=sg(ix+i2,jj)*sqrt(m(i))
C                   se_m(ix+i2,jj)=se(ix+i2,jj)*sqrt(m(i))
C                !end do
C             end do
C          end do
C       end do
      !https://doi.org/10.1039/C2CP41169E <-- good as overview, however the VH relations are not fulfilling sanity checks
      !https://doi.org/10.1002/9781118008720.ch8 <-- VH relations fulfill sanity checks
      !What I mean as sanity checks (if no mode is deleted):
      ! grad_q has dimension of 3N but only 3N-6 are non-zero
      ! J_VH has dimensions of 3Nx3N but only (3N-6)x(3N-6) are not purely one (the eigenvectors ie. columns do not have an element that equals exactly 1)
      ! K_VH has dimension of 3N but only 3N-6 are non-zero
      ! also N variable here is actually the 3x(no. of atoms)

c     K = s'^-1 sqrt(M) dX      
         do 4 ig=1,NQ
         k(ig)=0
         do 4 ia=1,nat
         do 4 ix=1,3
         ii=ix+3*(ia-1)
4        k(ig)=k(ig)+sg(ii,ig)*(re(ii)-rg(ii))/bohr
     1   *m(ia)*dsqrt(1822.89d0)

c     Duschinsky transformation:
c     Q' = JQ'' + K    (Q' initial, Q'' final)

c     J = s'^-1 s'', JT = J^t
         if(j_one)then
            call eye(jt,nq)
         else
            call mkj(nat,N,NQ,JT,sg,se,m)
         end if
c     scale the shift if kfac <> 1:

      if(Gauss_dusch>=2)then
         call readDusch_Gauss(NQ,N,A,B,C,D,E,J,K,Gauss_dusch_file)
      end if
      do 41 ig=1,NQ
41    k(ig)=k(ig)*kfac
      
      call mt(J,JT,NQ)
      call mm(T,G,J,NQ)
      call mm(TP,JT,T,NQ)
      call ms(F,TP,GP,NQ)
c     F = Jt G J + GP
      call INV(NQ,F,FI,NQ,TOL,E2,IERR)
      if(IERR.ne.0)call report('Inversion error')
      call mm(T,JT,GH,NQ)
      call mm(TP,FI,T,NQ)
      call mm(T,J,TP,NQ)
      call mm(TP,GH,T,NQ)
      
      if(.not.Gauss_dusch>2)then
      do 11 ix=1,NQ
      do 11 jx=1,NQ
11    A(ix,jx)=2.0d0*TP(ix,jx)
      do 12 ix=1,NQ
12    A(ix,ix)=A(ix,ix)-1.0d0
c     A=2 G^1/2 J Fi Jt G^1/2 - I
      end if
      
      call mm(T,JT,G,NQ)
      call mm(TP,FI,T,NQ)
      call mm(T,J,TP,NQ)
      do 13 ix=1,NQ
13    T(ix,ix)=T(ix,ix)-1.0d0
      call mv(TV,T,K,NQ)
      call mv(TU,GH,TV,NQ)
      
      if(.not.Gauss_dusch>2)then
      do 14 ix=1,NQ
14    B(ix)=-2.0d0*TU(ix)
c     B= - 2 G^1/2 ( J Fi Jt G - I) K
      end if
      
      call mm(T,FI,GHP,NQ)
      call mm(TP,GHP,T,NQ)
      
      if(.not.Gauss_dusch>2)then
      do 15 ix=1,NQ
      do 15 jx=1,NQ
15    C(ix,jx)=2.0d0*TP(ix,jx)
      do 16 ix=1,NQ
16    C(ix,ix)=C(ix,ix)-1.0d0
c     C= 2 GP^1/2 Fi GP^1/2 - I
      end if
      
      call mv(TV,G,K,NQ)
      call mv(TU,JT,TV,NQ)
      call mv(TV,FI,TU,NQ)
      call mv(TU,GHP,TV,NQ)
      if(.not.Gauss_dusch>2)then
      do 17 ix=1,NQ
17    D(ix)=-2.0d0*TU(ix)
c     D= - 2 GP^1/2 Fi Jt G K
      end if
      call mm(T,JT,GH,NQ)
      call mm(TP,FI,T,NQ)
      call mm(T,GHP,TP,NQ)
      if(.not.Gauss_dusch>2)then
      do 18 ix=1,NQ
      do 18 jx=1,NQ
18    E(ix,jx)=4.0d0*T(ix,jx)
c     E= 4 GP^1/2 Fi Jt G^1/2
      end if
      write(6,6555)'A ',1,1,NQ,NQ,A(1,1),A(NQ,NQ)
      write(6,6555)'B ',1,1,NQ,NQ,B(1),B(NQ)
      write(6,6555)'C ',1,1,NQ,NQ,C(1,1),C(NQ,NQ)
      write(6,6555)'D ',1,1,NQ,NQ,D(1),D(NQ)
      write(6,6555)'E ',1,1,NQ,NQ,E(1,1),E(NQ,NQ)
      write(6,6555)'F ',1,1,NQ,NQ,F(1,1),F(NQ,NQ)
      write(6,6555)'Fi',1,1,NQ,NQ,FI(1,1),FI(NQ,NQ)
      write(6,6555)'G ',1,1,NQ,NQ,G(1,1),G(NQ,NQ)
      write(6,6555)'GP',1,1,NQ,NQ,GP(1,1),GP(NQ,NQ)
      write(6,6555)'J ',1,1,NQ,NQ,J(1,1),J(NQ,NQ)
      write(6,6555)'JT',1,1,NQ,NQ,JT(1,1),JT(NQ,NQ)
      write(6,6555)'K ',1,1,NQ,NQ,K(1),K(NQ)
6555  format(1x,a2,2i3,1x,2i3,2f12.6)
      d2=0
      do i = 1,NQ
         d2=d2+K(i)**2
      end do
      
      J2_mean=0
      do i = 1,NQ
         J2_mean=J2_mean+maxval(J(:,i)**2,dim=1)
      end do
      J2_mean=J2_mean/nq
      J2_std=0d0
      do i = 1,NQ
         J2_std=J2_std+(maxval(J(:,i)**2,dim=1)-J2_mean)**2
      end do
      J2_std=sqrt(J2_std/nq)
      dw_mean=0
      do i = 1,nq
         dw_mean=dw_mean+abs(eg(i)-ep(i))
      end do
      dw_mean=dw_mean/nq
      dw_std=0
      do i = 1,NQ
         dw_std=dw_std+(abs(eg(i)-ep(i))-dw_mean)**2
      end do
      dw_std=sqrt(dw_std/nq)
      
      write(6,'(A,G15.9)')'D^2 = ',d2
      write(6,'(A,F12.6," +- ",F12.6)')'< (J_max)^2 > = ',J2_mean,
     1    J2_std
      write(6,'(A,F12.6," +- ",F12.6)')'< dW > = ',dw_mean*au2cm,
     1 dw_std*au2cm
      write(6,*)
      
      call mm(T,G,GP,NQ)
      if(Gauss_dusch==1)THEN
         call readDusch_Gauss(NQ,N,A,B,C,D,E,J,K,Gauss_dusch_file)
      end if
      call dofac(NQ,fac,J,JT,F,FI,G,GP,K)
      overlap=fac
      
      if(ic.eq.1) !but why write the J^T instead of J and also not mentioning it in DUSCH.OUT...omg
     1call wdusch(NQ,N,A,B,C,D,E,J,K,eg,ep,sg,se,rg,re,m)
      
      return
      end
      
      function overlap_ABCDE00(ic,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1 m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,Gauss_dusch,
     1 Gauss_dusch_file)
     1 result(ff)
         implicit none
         integer*4 ic,NQ,N,ig,ia,nat,ix,jx,ii,IERR,i,i2,imag_c,
     1   j2,jj,ifx,idx,MO(NQ),nqn,kk,dels,Gauss_dusch
         real*8 eg(NQ),ep(NQ),G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ),
     1   J(NQ,NQ),JT(NQ,NQ),sg(N,NQ),se(N,NQ),m(nat),k(NQ),rg(N),re(N),
     2   bohr,T(NQ,NQ),TP(NQ,NQ),F(NQ,NQ),FI(NQ,NQ),E2(NQ,2*NQ),TOL,
     3   A(NQ,NQ),TV(*),TU(*),B(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ),overlap,
     4   fac,kfac,toll,wfix,sg_full(N,N),se_full(N,N),d2,curfreq,CFN,ff,
     5   J2_mean,J2_std,dw_mean,dw_std
         character(*) Gauss_dusch_file
         logical :: vert_h,j_one,cor,swe,we_eq_wg,vert_h_w
         real*8, allocatable ::
     5    grad_gr(:),grad_exc(:),grad_gr_q(:),grad_exc_q(:),
     6   ffc(:,:),ffq(:,:),ffq_i(:,:),W_VH2(:),sg_m(:,:),se_m(:,:),
     6   help(:,:),
     7   m_invsq(:,:),ffc_a(:,:),W_VH2_inv(:,:),K_VH(:),J_VH(:,:),
     8   test(:,:),ws(:),js(:,:),ks(:)
         real*8,parameter :: amu2au=1822.89d0
         real*8,parameter :: au2cm=1d0/4.556431403d-6
      
         call makeG(NQ,G,GP,GH,GHP,eg,ep)
         
         if(Gauss_dusch>=2)then
            call readDusch_Gauss(NQ,N,A,B,C,D,E,J,K,Gauss_dusch_file)
         end if
         call mt(J,JT,NQ)
         call mm(T,G,J,NQ)
         call mm(TP,JT,T,NQ)
         call ms(F,TP,GP,NQ)
c        F = Jt G J + GP
         call INV(NQ,F,FI,NQ,TOL,E2,IERR)
         if(IERR.ne.0)call report('Inversion error')
         call mm(T,JT,GH,NQ)
         call mm(TP,FI,T,NQ)
         call mm(T,J,TP,NQ)
         call mm(TP,GH,T,NQ)
         
         if(.not.Gauss_dusch>2)then
         do 11 ix=1,NQ
         do 11 jx=1,NQ
11       A(ix,jx)=2.0d0*TP(ix,jx)
         do 12 ix=1,NQ
12       A(ix,ix)=A(ix,ix)-1.0d0
c        A=2 G^1/2 J Fi Jt G^1/2 - I
         end if
         
         call mm(T,JT,G,NQ)
         call mm(TP,FI,T,NQ)
         call mm(T,J,TP,NQ)
         do 13 ix=1,NQ
13       T(ix,ix)=T(ix,ix)-1.0d0
         call mv(TV,T,K,NQ)
         call mv(TU,GH,TV,NQ)
         
         if(.not.Gauss_dusch>2)then
         do 14 ix=1,NQ
14       B(ix)=-2.0d0*TU(ix)
c     B= - 2 G^1/2 ( J Fi Jt G - I) K
         end if
         
         call mm(T,FI,GHP,NQ)
         call mm(TP,GHP,T,NQ)
         
         if(.not.Gauss_dusch>2)then
         do 15 ix=1,NQ
         do 15 jx=1,NQ
15       C(ix,jx)=2.0d0*TP(ix,jx)
         do 16 ix=1,NQ
16       C(ix,ix)=C(ix,ix)-1.0d0
c     C= 2 GP^1/2 Fi GP^1/2 - I
         end if
         
         call mv(TV,G,K,NQ)
         call mv(TU,JT,TV,NQ)
         call mv(TV,FI,TU,NQ)
         call mv(TU,GHP,TV,NQ)
         
         if(.not.Gauss_dusch>2)then
         do 17 ix=1,NQ
17       D(ix)=-2.0d0*TU(ix)
c     D= - 2 GP^1/2 Fi Jt G K
         end if
         
         call mm(T,JT,GH,NQ)
         call mm(TP,FI,T,NQ)
         call mm(T,GHP,TP,NQ)
         
         if(.not.Gauss_dusch>2)then
         do 18 ix=1,NQ
         do 18 jx=1,NQ
18       E(ix,jx)=4.0d0*T(ix,jx)
c     E= 4 GP^1/2 Fi Jt G^1/2
         end if
         
         write(6,6555)'A ',1,1,NQ,NQ,A(1,1),A(NQ,NQ)
         write(6,6555)'B ',1,1,NQ,NQ,B(1),B(NQ)
         write(6,6555)'C ',1,1,NQ,NQ,C(1,1),C(NQ,NQ)
         write(6,6555)'D ',1,1,NQ,NQ,D(1),D(NQ)
         write(6,6555)'E ',1,1,NQ,NQ,E(1,1),E(NQ,NQ)
         write(6,6555)'F ',1,1,NQ,NQ,F(1,1),F(NQ,NQ)
         write(6,6555)'Fi',1,1,NQ,NQ,FI(1,1),FI(NQ,NQ)
         write(6,6555)'G ',1,1,NQ,NQ,G(1,1),G(NQ,NQ)
         write(6,6555)'GP',1,1,NQ,NQ,GP(1,1),GP(NQ,NQ)
         write(6,6555)'J ',1,1,NQ,NQ,J(1,1),J(NQ,NQ)
         write(6,6555)'JT',1,1,NQ,NQ,JT(1,1),JT(NQ,NQ)
         write(6,6555)'K ',1,1,NQ,NQ,K(1),K(NQ)
6555     format(1x,a2,2i3,1x,2i3,2f12.6)
         d2=0
         do i = 1,NQ
            d2=d2+K(i)**2
         end do
         write(6,'(A,G15.9)')'D^2 = ',d2
         
         J2_mean=0
         do i = 1,NQ
            J2_mean=J2_mean+maxval(J(:,i)**2,dim=1)
         end do
         J2_mean=J2_mean/nq
         J2_std=0d0
         do i = 1,NQ
            J2_std=J2_std+(maxval(J(:,i)**2,dim=1)-J2_mean)**2
         end do
         J2_std=sqrt(J2_std/nq)
         
         dw_mean=0
         do i = 1,nq
            dw_mean=dw_mean+abs(eg(i)-ep(i))
         end do
         dw_mean=dw_mean/nq
         dw_std=0
         do i = 1,NQ
            dw_std=dw_std+(abs(eg(i)-ep(i))-dw_mean)**2
         end do
         dw_std=sqrt(dw_std/nq)
         
         write(6,'(A,F12.6," +- ",F12.6)')'< (J_max)^2 > = ',J2_mean,
     1    J2_std
         write(6,'(A,F12.6," +- ",F12.6)')'< dW > = ',dw_mean*au2cm,
     1    dw_std*au2cm
         write(6,*)
         
         call mm(T,G,GP,NQ)

         call dofac(NQ,fac,J,JT,F,FI,G,GP,K)
         ff=fac
         
         if(Gauss_dusch==1)then
            call readDusch_Gauss(NQ,N,A,B,C,D,E,J,K,Gauss_dusch_file)
         end if

         if(ic.eq.1) !but why write the J^T instead of J and also not mentioning it in DUSCH.OUT...omg
     1   call wdusch(NQ,N,A,B,C,D,E,J,K,eg,ep,sg,se,rg,re,m)

      end function overlap_ABCDE00
      !Calculate Frequency Number (i)
      function CFN(ii,ffq_i,J_VH,N)result(freq)
         integer N,ii,jj,kk
         double precision ffq_i(N,N),J_VH(N,N),freq,help(N,N)
         
         freq=0d0
         do jj = 1,N
            help=0d0
            do kk = 1,N
               help(ii,jj)=help(ii,jj)+J_VH(kk,ii)*ffq_i(kk,jj)
            end do
            freq=freq+help(ii,jj)*J_VH(jj,ii)
         end do
         if(freq<0)then
            freq=-1d0/sqrt(abs(freq))
         else if(abs(freq)<1d-10)then
            freq=0d0
         else
            freq=1d0/sqrt(abs(freq))
         end if
      end function CFN
      
      function CKN(ii,ffq_i,grad_q,n)result(k_res)
         integer ii,n
         double precision ffq_i(n,n),grad_q(n),k_res
         
      end function CKN
      
      subroutine ReorderVec_R(arr,N,order)
         implicit none
         integer N,order(N),i
         real*8 arr(N),arr_help(N)
         
         arr_help=arr
         do i = 1,N
            arr(i)=arr_help(order(i))
         end do
      end subroutine ReorderVec_R
      
      !maybe???
      subroutine ReorderMat_R(mat,N,order)
         implicit none
         integer N,order(N),i,j
         real*8 mat(N,N),mat_help(N,N)
         
         mat_help=mat
         do i = 1,N
            do j = 1,N
               mat(i,j)=mat_help(order(i),order(j))
            end do
         end do
      end subroutine ReorderMat_R
      
c     ============================================================
      subroutine mkj(nat,N,NQ,JT,sg,se,m)
      implicit none
      integer*4 NQ,ig,ie,ia,ix,ii,nat,N
      real*8 JT(NQ,NQ),sg(N,N),se(N,N),m(*)
      do 3 ig=1,NQ
      do 3 ie=1,NQ
      JT(ie,ig)=0.0d0
      do 3 ia=1,nat
      do 3 ix=1,3
      ii=ix+3*(ia-1)
3     JT(ie,ig)=JT(ie,ig)+sg(ii,ig)*m(ia)*se(ii,ie)
      return
      end
      
      subroutine mkk(bohr,nat,N,NQ,m,rg,re,sg,k)
         implicit none
         integer nat,nq,N,ig,ia,ix,ii
         double precision bohr,rg(N),re(N),k(nq),m(nat),sg(N,N)
c        K = s'^-1 sqrt(M) dX      
         do 4 ig=1,NQ
         k(ig)=0
         do 4 ia=1,nat
         do 4 ix=1,3
         ii=ix+3*(ia-1)
4        k(ig)=k(ig)+sg(ii,ig)*(re(ii)-rg(ii))/bohr
     1   *m(ia)*dsqrt(1822.89d0)
         open(792,file='VEC.K')
         do ii = 1,NQ
            write(792,'(1X,D14.6)')K(ii)
         end do
         close(792)
      end subroutine mkk
      
      subroutine makeSM(nat,N,nq,s,m,s_m,iz)
         implicit none
         integer nat,N,i,ix,i2,jj,nq,iz
         double precision s(N,NQ),s_m(N,NQ),m(nat)
         
         s_m=0
         do i = 1,nat
            ix=(i-1)*3+1 
            do i2 = 0,2 !x y z
               do jj = 1,NQ
                  s_m(ix+i2,jj)=s(ix+i2,jj)*sqrt(m(i))
               end do
            end do
         end do
      end subroutine makeSM
      
      subroutine mk_ffqgr(nat,N,NQ,s,m,m_invsq,ffc,ffq,gr,gr_q,amu2au,
     1 iz,exc,s_m)
         implicit none
         integer nat,n,nq,iz
         double precision :: s(N,NQ),ffc(N,N),ffq(NQ,NQ),gr(N),gr_q(NQ),
     1    m(nat),m_invsq(N,N),amu2au,s_m(N,NQ),FFc_A(N,N)
         logical exc
         
         
         call makeSM(nat,n,nq,s,m,s_m,iz)
         if(exc)then
            call readgrad_dusch(nat,gr,'excited/FILE.GR') !read excited state gradient
            call READFF(N,ffc,'excited/FILE.FC') !read excited state force field
            call PROJECT(FFc_A,FFc,N,'excited/FILE.X') 
         else
            call readgrad_dusch(nat,gr,'ground/FILE.GR') !read excited state gradient
            call READFF(N,ffc,'ground/FILE.FC') !read excited state force field
            call PROJECT(FFc_A,FFc,N,'ground/FILE.X') 
         end if
         
         call grad2Q(nat,nq,gr,gr_q,m,m_invsq,s_m,amu2au) !convert gradient to normal mode coordinates
         call FFc2q(nat,nq,ffc,ffq,m_invsq,s_m) !convert force field to normal mode coordinates
      end subroutine mk_ffqgr
      
      subroutine mkj_vert(tol,amu2au,nat,N,NQ,m,sg_m,m_invsq,
     1 grad_exc_q,ffq,ffq_i,J,K,ws,J_one,we_eq_wg,wg)
         implicit none
         integer*4 NQ,N,nat,ierr,nqn,imag_c,NQo,i,imag_c_real
         double precision J(NQ,NQ),K(NQ),sg_m(N,NQ),m_invsq(N,N),m(nat),
     1   ffc(N,N),ffc_a(N,N),ffq(NQ,NQ),ffq_i(NQ,NQ),
     1   grad_exc(N),grad_exc_q(NQ),
     1   help(N,2*N),J_VH(NQ,NQ),K_VH(NQ),W_VH2(NQ),W_VH2_mat(NQ,NQ),
     1   ks(nq),js(nq,nq),ws(nq),amu2au,wg(NQ),buf,tol,toll,
     1   ffq_ch(nq,nq),J_VHi(NQ,NQ),W_VHi2(NQ)
         logical J_one,we_eq_wg
         
         call readgrad_dusch(nat,grad_exc,'excited/FILE.GR') !read excited state gradient
         call grad2Q(nat,NQ,grad_exc,grad_exc_q,m,m_invsq,sg_m,amu2au) !convert gradient to normal mode coordinates
         call READFF(N,ffc,'excited/FILE.FC') !read excited state force field
         call PROJECT(FFc_A,FFc,N,'excited/FILE.X') 
         
         call FFc2q(nat,nq,ffc,ffq,m_invsq,sg_m) !convert force field to normal mode coordinates
         !call WRITEFF()
         help=0d0
         ffq_i=0d0
         toll=tol
         call INV(NQ,ffq,ffq_i,NQ,toll,help,ierr)
         if(ierr/=0)then
            write(6,*)
     1      'Unable to invert Hessian in NM coordinates'
            call exit(101)
         end if
         
         call TRED12(NQ,ffq_i,J_VHi,W_VHi2,2,IERR)
         call TRED12(NQ,ffq,J_VH,W_VH2,2,IERR)
         if(ierr/=0)then
            write(6,*)
     1      'Unable to diagonalize FFq_i'
            call exit(102)
         end if
         
C          !check that J columns (or rows, it doesnt matter) are normalized and non zero
C          do i = 1,NQ
C             buf=0d0
C             buf=sum(J_VH(:,i)**2)
C             if(abs(buf-1d0)>1d-6)THEN
C                stop 3
C             end if
C          end do
         
         
         
         imag_c=count(W_VH2<0d0)
         !imag_c=findloc((W_VH2)<1d-6,.true.,dim=1,back=.true.)-nq+1
         nqn=NQ
         ws=0
         js=0
         ks=0
         
         
C          !check that J columns (or rows, it doesnt matter) are normalized and non zero
C          do i = 1,NQ
C             buf=0d0
C             buf=sum(J_VH(:,i)**2)
C             if(abs(buf-1d0)>1d-6)THEN
C                stop 3
C             end if
C          end do
         
         
         call ReverseArr(W_VH2,nq)
         W_VH2_mat=0d0
         do i = 1,NQ
           W_VH2_mat(i,i)=(W_VH2(i))
         end do
         call ReverseJ(J_VH,nq)
         
         ffq_ch=matmul(matmul(J_VH,W_VH2_mat),transpose(J_VH))
C          if(J_one)then
C             call Eye(J_VH,NQ)
C          end if
         if(we_eq_wg)then
            W_VH2=0d0
            W_VH2(1:NQ)=(wg(1:NQ))**2
            imag_c_real=count(wg<0d0)
         end if
         call mkk_vert(NQ,W_VH2,J_VH,grad_exc_q,K_VH)
C          W_VH2_mat=0d0
C          do i = 1,NQ
C            W_VH2_mat(i,i)=(1d0/W_VH2(i))
C          end do
         !W_VH2_mat=matmul(transpose(J_VH),matmul(ffq_i,(J_VH)))!correct?
C          K_VH=-1*matmul(matmul(J_VH,matmul(W_VH2_mat,transpose(J_VH)))
C      1    ,grad_exc_q)
C          K_VH=-1*matmul(ffq_i,grad_exc_q)
C          K_VH=-1*matmul(ffq_i,grad_exc_q)
         
         !K_VH=-1*matmul(ffq_i,grad_exc_q)
         !ws(imag_c+1:)=1d0/sqrt(W_VH2(1:nqn-imag_c))
C          ws(imag_c+1:)=(W_VH2(1:nqn-imag_c))
C          js(1:nqn,imag_c+1:nqn)=J_VH(1:nqn,1:nqn-imag_c)
C          if(imag_c>0)then
C             !ws(1:imag_c)=-1d0/sqrt(abs(W_VH2(N-imag_c+1:N)))
C             ws(1:imag_c)=abs(W_VH2(NQ-imag_c+1:NQ))
C             js(1:nqn,1:imag_c)=J_VH(1:nqn,NQ-imag_c+1:NQ)
C          end if
C          call eye(J_VH,NQ)
C          J_VH(1:NQ,1:NQ)=js
C          ks=K_VH(1:nqn)
         !if(imag_c>0)then
         !end if
         ws=sgn(W_VH2)*sqrt(abs(W_VH2))
         J=J_VH
         K=K_VH
      end subroutine mkj_vert
      
      elemental function sgn(num)
         double precision,intent(in) :: num
         integer sgn
         if(num>=0d0)then
            sgn=1
         else
            sgn=-1
         end if
      end function sgn
      
      subroutine ReverseArr(arr,n)
         integer n
         double precision arr(n),buf
         
         integer i
         
         do i = 1,n/2
            buf=arr(i)
            arr(i)=arr(n-i+1)
            arr(n-i+1)=buf
         end do
      end subroutine ReverseArr
      
      subroutine ReverseJ(J,n)
         integer n
         double precision J(n,n),buf(n)
         
         integer i
         
         do i = 1,n/2
            buf=J(:,i)
            J(:,i)=J(:,n-i+1)
            J(:,n-i+1)=buf
C             buf=J(i,:)
C             J(i,:)=J(n-i+1,:)
C             J(n-i+1,:)=buf
         end do
      end subroutine ReverseJ
      
      
      
      function sort(arr,order,n)result(arr_new)
         integer n,order(n)
         double precision arr(n),arr_new(n)
         
         integer i,j
         
         
      end function sort
      
      subroutine mkk_vert(NQ,wg,J,grad_exc_q,K)
         implicit none
         integer nq,ii,i
         double precision K(NQ),grad_exc_q(NQ),J(NQ,NQ),wg(nq)
         double precision W_VH2_mat(NQ,NQ)
         
         W_VH2_mat=0d0
         do i = 1,NQ
           W_VH2_mat(i,i)=(1d0/wg(i)**2)
         end do
         K=-1*matmul(matmul(J,matmul(W_VH2_mat,transpose(J)))
     1    ,grad_exc_q)
      end subroutine mkk_vert
      
      subroutine mkk_vert_ffq(NQ,ffq,grad_exc_q,K)
         implicit none
         integer nq,ii,i
         double precision K(NQ),grad_exc_q(NQ),ffq(NQ,NQ)
         
         K=-1*matmul(ffq,grad_exc_q)
      end subroutine mkk_vert_ffq
      
c     ============================================================
      subroutine iterate(maxit,NQ,lwr,wmin,
     1N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac,
     1 vert_h,j_one,ifx,wfix,cor,swe,MO,sg_full,se_full,WE_IS_WG,
     1 vert_h_w)
      implicit none
      integer*4 iter,maxit,iw,NQ,ichange,IERR,N,nat,ifx,MO(NQ)
      real*8 dd,wp,wm,o,om,op,w,CM,wmin,delw,oold,kfac,wfix
      real*8 eg(*),ep(*),G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ),
     1J(NQ,NQ),JT(NQ,NQ),sg(N,N),se(N,N),m(*),k(*),rg(*),re(*),
     2bohr,T(NQ,NQ),TP(NQ,NQ),F(NQ,NQ),FI(NQ,NQ),E2(NQ,2*NQ),TOL,
     3A(NQ,NQ),TV(*),TU(*),B(*),C(NQ,NQ),D(*),E(NQ,NQ),
     4sg_full(N,N),se_full(N,N)
      logical lwr,vert_h,j_one,cor,swe,WE_IS_WG,vert_h_w
      CM=219474.630d0
      dd=0.5d0
      oold=1.0d0
      do 1 iter=1,maxit
      ichange=0
      if(lwr)write(6,*)'Iteration ',iter,' of ',maxit
      if(lwr)write(6,601)
601   format(3x,'iw',11x,'w',10x,'wp',10x,'wm',11x,'o',10x,'op',
     110x,'om')
      do 2 iw=1,NQ
      w=ep(iw)
      delw=max(wmin/CM,dd*w)
      o=overlap(0,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac, 
     1 vert_h,j_one,ifx,wfix,cor,swe,MO,sg_full,se_full,WE_IS_WG,
     1 vert_h_w,0,'')
      wp=w+delw
      ep(iw)=wp
      op=overlap(0,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac,
     1 vert_h,j_one,ifx,wfix,cor,swe,MO,sg_full,se_full,WE_IS_WG,
     1 vert_h_w,0,'')
      wm=max(wmin/CM,w-delw)
      ep(iw)=wm
      om=overlap(0,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac, 
     1 vert_h,j_one,ifx,wfix,cor,swe,MO,sg_full,se_full,WE_IS_WG,
     1 vert_h_w,0,'')
      if(lwr)write(6,600)iw,w*CM,wp*CM,wm*CM,o,op,om
600   format(i5,3f12.3,3g12.3)
      if(om.gt.o.or.op.gt.o)then
       if(om.gt.op)then
        ep(iw)=wm
       else
        ep(iw)=wp
       endif
       ichange=ichange+1
      else
       ep(iw)=w
      endif
2     continue
      if(ichange.eq.0)dd=dd/2.0d0
      if(oold.eq.o)return
1     oold=o
      return
      end
c     ============================================================
      function spsv(sg,se,i,j,amas,z,N)
      implicit none
      integer*4 N,z(*),i,j,ix
      real*8 spsv,sg(N,N),se(N,N),dij
      real*8 amas(*)
      dij=0.0d0
      do 1 ix=1,N
1     dij=dij+sg(ix,i)*se(ix,j)*dble(amas(z((ix-1)/3+1)))
      spsv=dij
      return
      end
c     ============================================================
      subroutine findbest(i,jstart,jmin,dmin,dij,sg,se,N,amas,z)
      implicit none
      integer*4 i,jstart,jmin,N,z(*),j
      real*8 dmin,dij,sg(N,N),se(N,N)
      real*8 amas(*)
c     find best match for i within jstart...1
      jmin=0
      dmin=0.0d0
      dij=0.0d0
      do 14 j=jstart,1,-1
      dij=spsv(sg,se,i,j,amas,z,N)
      if(dabs(dij).gt.dmin.or.j.eq.jstart)then
       dmin=dabs(dij)
       jmin=j
      endif
14    continue
      return
      end
      
      subroutine mkj_col(ie,nat,N,NQ,J_col,sg,se,m) !ie. expansion of 1 excited state to NQ ground states
         implicit none
         integer*4 NQ,ig,ie,ia,ix,ii,nat,N
         real*8 J_col(NQ),sg(N,N),se(N,N),m(*)
         
         do 3 ig=1,NQ
         !do 3 ie=1,NQ
         J_col(ig)=0.0d0
         do 3 ia=1,nat
         do 3 ix=1,3
         ii=ix+3*(ia-1)
   3     J_col(ig)=J_col(ig)+sg(ii,ig)*m(ia)*se(ii,ie)
         return
      end subroutine mkj_col
      
      subroutine mkj_row(ig,nat,N,NQ,J_row,sg,se,m) !ie. expansion of 1 ground state to NQ excited states
         implicit none
         integer*4 NQ,ig,ie,ia,ix,ii,nat,N
         real*8 J_row(NQ),sg(N,N),se(N,N),m(*)
         
         !do 3 ig=1,NQ
         do 3 ie=1,NQ
         J_row(ie)=0.0d0
         do 3 ia=1,nat
         do 3 ix=1,3
         ii=ix+3*(ia-1)
   3     J_row(ie)=J_row(ie)+sg(ii,ig)*m(ia)*se(ii,ie)
         return
      end subroutine mkj_row
      
      
c     ============================================================
      subroutine gradq(nat,N,m,ldz)
      implicit none
      integer*4 nat,ia,ix,N,NQ,iq,i,j,iz
      real*8 SFAC,CM,wep,a,ap,x0,ef,o,m(*)
      logical lexg,lexe,ldz
      real*8,allocatable::g(:),q(:),u(:,:),qe(:),qt(:),r(:),
     1sg(:,:),se(:,:),eg(:),ep(:)
      allocate(sg(N,N),se(N,N),eg(N),ep(N))
      SFAC=0.0234280d0
      CM=219474.0d0
      inquire(file='ground/FILE.GR',exist=lexg)
      if(lexg)then
       write(6,*)
       write(6,*)' Ground state gradient found'
       write(6,*)' GG.TXT opened'
       open(66,file='GG.TXT')
       allocate(g(N),q(N),r(N))
       open(8,file='ground/FILE.GR')
       do 1 ia=1,nat
1      read(8,*)(g(ix+3*(ia-1)),ix=1,3)
       close(8)
       call readsi_dusch(N,sg,eg,NQ,'ground/F.INP',r,ldz,iz)
       do 2 iq=1,NQ
       q(iq)=0.0d0
       do 2 ia=1,nat
       do 2 ix=1,3
2      q(iq)=q(iq)+sg(ix+3*(ia-1),iq)*g(ix+3*(ia-1))*SFAC
       write(66,*)' ground state:'
       write(66,*)' mode        wg gradient (au)'
       do 3 iq=1,NQ
3      write(66,600)iq,eg(iq)*CM,q(iq)
600    format(i4,f12.2,f12.7)
       inquire(file='excited/FILE.GR',exist=lexe)
       if(lexe)then
        write(6,*)' Excited state gradient found, too'
        allocate(qe(N),qt(N))
        open(8,file='excited/FILE.GR')
        do 4 ia=1,nat
4       read(8,*)(g(ix+3*(ia-1)),ix=1,3)
        close(8)
c       reload because se might be transformed:
        call readsi_dusch(N,se,ep,NQ,'excited/F.INP',r,ldz,iz)

        do 5 iq=1,NQ
        qe(iq)=0.0d0
        do 5 ia=1,nat
        do 5 ix=1,3
5       qe(iq)=qe(iq)+se(ix+3*(ia-1),iq)*g(ix+3*(ia-1))*SFAC
        write(66,*)' excited state:'
        write(66,*)' mode        we gradient (au)'
        do 6 iq=1,NQ
6       write(66,600)iq,ep(iq)*CM,qe(iq)
        write(6,*)' .... gradient written'
        write(66,*)' Independent mode approximation:'
        write(66,*)' Excited state in ground state modes:'
c       u=se.s transformation matrix between ground and excited modes:
        allocate(u(NQ,NQ))
        do 7 i=1,NQ
        do 7 j=1,NQ
        u(i,j)=0.0d0
        do 7 ia=1,nat
        do 7 ix=1,3
7       u(i,j)=u(i,j)+sg(ix+3*(ia-1),i)*m(ia)*se(ix+3*(ia-1),j)
        write(66,602)
602     format(' mode         wg         wep gradient (au)'
     1  //' q0    <f00|F0e>')
        do 8 iq=1,NQ
        qt(iq)=0.0d0
c       diagonal part of Wij:
        wep=0.0d0
        do 9 i=1,NQ
        qt(iq)=qt(iq)+u(iq,i)*qe(i)
9       wep   =wep   +u(iq,i)*ep(i)*u(iq,i)
        a =eg(iq)/2.0d0
        ap=wep   /2.0d0
        x0=-qe(iq)/wep**2
        if(wep.lt.0.0d0)x0=999999999.99d0
        if(a.gt.0.0d0.and.ap.gt.0.0d0)then
         ef=x0**2*a*ap/(a+ap)
         if(ef.gt.20.0d0)then
          o=0.0d0
         else
          o=dsqrt(2.0d0*dsqrt(a*ap)/(a+a))*exp(-ef)
         endif
        else
         o=-999.99d0
        endif
8       write(66,601)iq,eg(iq)*CM,wep*CM,qe(iq),x0,o
601     format(i4,2f12.2,3g12.4)
        write(6,*)' .... transformed into ground modes'
       endif
       close(66)
       write(6,*)' GG.TXT closed'
       write(6,*)
      endif
      return
      end
c     ============================================================
      subroutine grade(nat,N,ldz)
      implicit none
      integer*4 nat,ia,ix,N,NQ,iq,ii,iz
      real*8 SFAC,CM,y,ani,w,bohr
      logical lex,ldz
      real*8,allocatable::g(:),q(:),rg(:),x(:),qt(:),ee(:),r(:),
     1sg(:,:),eg(:)
      integer*4 ,allocatable::z(:)
      allocate(sg(N,N),eg(N))
      SFAC=0.0234280d0
      CM=219474.0d0
      y=50.0d0/CM
      bohr=0.529177d0
      inquire(file='groundtd/FILE.GR',exist=lex)
      if(lex)then
       write(6,*)' Excited state gradient for ground state geometry'
       allocate(g(N),q(N),qt(N),ee(N),r(N))
c      cartesian gradient:
       open(8,file='groundtd/FILE.GR')
       do 1 ia=1,nat
1      read(8,*)(g(ix+3*(ia-1)),ix=1,3)
       close(8)
       call readsi_dusch(N,sg,eg,NQ,'groundtd/F.INP',r,ldz,iz)
c      normal mode gradient:
       do 2 iq=1,NQ
       q(iq)=0.0d0
       do 2 ii=1,N
2      q(iq)=q(iq)+sg(ii,iq)*g(ii)*SFAC
       open(9,file='PGRAD')
       write(6,*)' mode        wg          wgp gradient (au)  Qt   n'
       write(9,*)' mode        wg          wgp gradient (au)  Qt   n'
       do 3 iq=1,NQ
       if(eg(iq).lt.y)then
        w=(q(iq)**2/2.0d0)**(1.0d0/3.0d0)
       else
        w=eg(iq)
       endif
       ee(iq)=w
       qt(iq)=q(iq)/w**2
       ani=w*qt(iq)**2/2.0d0-0.50d0
       write(9,600)iq,eg(iq)*CM,w*cm,q(iq),qt(iq),nint(ani+0.50d0)
3      write(6,600)iq,eg(iq)*CM,w*cm,q(iq),qt(iq),nint(ani+0.50d0)
600    format(i4,2f12.2,2g12.4,i6)
       close(9)
       write(6,*)'PGRAD written'

       allocate(z(nat),rg(N),x(N))
       call rdxx('ground/FILE.X',1,  nat,z,rg)
       do 4 ii=1,N
       x(ii)=rg(ii)
       do 4 iq=1,NQ
4      x(ii)=x(ii)+qt(iq)*sg(ii,iq)*SFAC*bohr
       call system('mkdir virtual')
       call wdxx('virtual/FILE.X',nat,z,x)
       call writesi_dusch(nat,N,sg,ee,NQ,'virtual/F.INP',x,z)
       write(6,*)'virtual geometry made'
      endif
      return
      end
c     ============================================================
      subroutine shiftq(nat,N,NQ,se,ep,sg,eg)
      implicit none
      integer*4 nat,N,NQ,i,j,ii,IERR,nt,it,jc,ia,jp
      real*8 se(N,N),ep(*),CM,ani,d,bohr,SFACR,TOL,t(3,3),pt,
     1sg(N,N),eg(*),dc,pj,sgn,sn
      real*8,allocatable::sms(:,:),qt(:),xt(:),re(:),rg(:),x(:),
     1v(:),si(:,:),e(:,:),pc(:),st(:,:),u(:)
      integer*4,allocatable::z(:),iind(:)
      TOL=1.0d-10
      CM=219474.0d0
      bohr=0.529177d0
      SFACR=0.02342179d0
      allocate(sms(NQ,NQ),qt(NQ),xt(N),z(nat),re(N),rg(N),x(N),
     1v(NQ),si(NQ,NQ),e(NQ,2*NQ),pc(NQ),u(N),iind(NQ),st(N,NQ))

      call rdxx('ground/FILE.X',1,  nat,z,rg)
      call rdxx('exc.as.ground.x',1,nat,z,re)
      d=0.0d0
      do 5 i=1,N
      x(i)=(rg(i)-re(i))/bohr
5     d=d+x(i)**2
      write(6,607)d
607   format(' d = ',g12.4)

      do 1  i=1,NQ
      v(i)=0.0d0
      do 11 ii=1,N
11    v(i)=v(i)+x(ii)*se(ii,i)*SFACR
      do 1 j=1,NQ
      sms(i,j)=0.0d0
      do 1 ii=1,N
1     sms(i,j)=sms(i,j)+se(ii,i)*SFACR**2*se(ii,j)

      call INV(NQ,sms,si,NQ,TOL,e,IERR)
      if(IERR.ne.0)call report('inverson matrix cannot be found')

      do 2 i=1,NQ
      qt(i)=0.0d0
      do 2 j=1,NQ
2     qt(i)=qt(i)+si(i,j)*v(j) 

      d=0.0d0
      do 7 ii=1,N
      xt(ii)=re(ii)/bohr
      do 71 j=1,NQ
71    xt(ii)=xt(ii)+se(ii,j)*SFACR*qt(j)
      d=d+(rg(ii)/bohr-xt(ii))**2
7     xt(ii)=xt(ii)*bohr
      write(6,607)d
      open(44,file='test1.x')
      write(44,*)'Excited orientation to match the ground:'
      write(44,*)nat
      do 13 i=1,nat
13    write(44,608)z(i),(xt(ii+3*(i-1)),ii=1,3)
608   format(i4,3f12.6)
      close(44)
      
      nt=0
      open(9,file='QSTATE')
      write(9,600)
600   format(/,' mode       freq      Qt      n  nround:')
      do 3 i=1,NQ
      ani=ep(i)*qt(i)**2/2.0d0-0.50d0
      nt=nt+nint(ani+0.50d0)
3     write(9,601)i,ep(i)*CM,qt(i),ani,nint(ani+0.50d0)
 601  format(i4,f12.2,2G12.4,i6)
      close(9)
      write(6,*)nt,' - excited'
      write(6,*)'QSTATE written'

      write(6,610)'Excited state vibrations'
610   format(a24,'Projection of geometry difference to normal modes')
      open(44,file='TMX.TXT')
      do 9 i=1,3
9     read(44,*)(t(j,i),j=1,3)
      close(44)
c     excited state S-matrix in ground-like orientation:
      do 14 j=1,NQ
      do 14 ia=1,N/3
      do 14 i=1,3
      ii=3*(ia-1)
14    st(i+ii,j)=se(1+ii,j)*t(1,i)+se(2+ii,j)*t(2,i)+se(3+ii,j)*t(3,i)
c     as_ground_ai=tj,i exc_aj
      do 4 j=1,NQ
      iind(j)=j
c     unit vector along the mode:
      do 17 i=1,N
17    u(i)=st(i,j)
      call norm(u,N)
c     projection to the geometry change:
      pc(j)=dabs(sp(u,x,N))
c     in units of normal mode amplitude (~ 1/sqrt(w)):      
4     pc(j)=pc(j)*dsqrt(ep(j))
c     normalize:
      call norm(pc,NQ) !if the ground and excited geometries are the same, this here throws an arithmetic exception (DIV BY 0)
c     order
      do 10 j=1,NQ-1
      do 10 jp=j+1,NQ
      if(pc(jp).gt.pc(j))then
       pt=pc(jp)
       pc(jp)=pc(j)
       pc(j)=pt
       it=iind(jp)
       iind(jp)=iind(j)
       iind(j)=it
      endif
10    continue
      write(6,609)
609   format(' Mode exc   contrinution (%)  closest ground' )
      do 12 ii=1,NQ
      j=iind(ii)
c     most similar ground state mode:
      sgn=0.0d0
      sn=0.0d0
      jc=1
      dc=0.0d0
      do 15 i=1,N
      sgn=sgn+sg(i,1)**2
      sn =sn +st(i,j)**2
15    dc=dc+sg(i,1)*st(i,j)
      dc=dabs(dc)/dsqrt(sn*sgn)
      do 19 jp=2,NQ
      pj=0.0d0
      sgn=0.0d0
      sn=0.0d0
      do 16 i=1,N
      sgn=sgn+sg(i,jp)**2
      sn =sn +st(i,j)**2
16    pj=pj+sg(i,jp)*st(i,j)
      pj=dabs(pj)/dsqrt(sn*sgn)
      if(pj.gt.dc)then
       dc=pj
       jc=jp
      endif
19    continue
12    write(6,612)ii,j,ep(j)*CM,pc(ii)**2*100.0d0,jc,eg(jc)*CM
612   format(2i5,'(',f12.3,' cm-1)',f10.2,i6,'(',f12.2,' cm-1)')

      write(6,610)'Ground state vibrations '
      do 20 j=1,NQ
      iind(j)=j
      do 21 i=1,N
21    u(i)=sg(i,j)
      call norm(u,N)
      pc(j)=dabs(sp(u,x,N))
20    pc(j)=pc(j)*dsqrt(eg(j))
      call norm(pc,NQ)
      do 24 j=1,NQ-1
      do 24 jp=j+1,NQ
      if(pc(jp).gt.pc(j))then
       pt=pc(jp)
       pc(jp)=pc(j)
       pc(j)=pt
       it=iind(jp)
       iind(jp)=iind(j)
       iind(j)=it
      endif
24    continue
      write(6,619)
619   format(' Mode gr   contrinution (%)  closest exc' )
      do 25 ii=1,NQ
      j=iind(ii)
      sgn=0.0d0
      sn=0.0d0
      jc=1
      dc=0.0d0
      do 26 i=1,N
      sgn=sgn+st(i,1)**2
      sn =sn +sg(i,j)**2
26    dc=dc+st(i,1)*sg(i,j)
      dc=dabs(dc)/dsqrt(sn*sgn)
      do 27 jp=2,NQ
      pj=0.0d0
      sgn=0.0d0
      sn=0.0d0
      do 28 i=1,N
      sgn=sgn+st(i,jp)**2
      sn =sn +sg(i,j)**2
28    pj=pj+st(i,jp)*sg(i,j)
      pj=dabs(pj)/dsqrt(sn*sgn)
      if(pj.gt.dc)then
       dc=pj
       jc=jp
      endif
27    continue
25    write(6,612)ii,j,eg(j)*CM,pc(ii)**2*100.0d0,jc,ep(jc)*CM

      return
      end
c     ============================================================
      subroutine shifti(nat,N,NQ,se,ep)
c     iterative way
      implicit none
      integer*4 nat,ia,N,NQ,i,j,ii
      real*8 se(N,N),ep(*),CM,ani,sx,bohr,s2,SFACR,d
      real*8,allocatable::qt(:),rg(:),re(:),x(:),xt(:)
      integer*4,allocatable::z(:)
      CM=219474.0d0
      bohr=0.529177d0
      SFACR=0.02342179d0

      allocate(z(nat),rg(N),re(N),x(N),xt(N))
      call rdxx('ground/FILE.X',1  ,nat,z,rg)
      call rdxx('exc.as.ground.x',1,nat,z,re)
      d=0.0d0
      do 5 i=1,N
      x(i)=(rg(i)-re(i))/bohr
5     d=d+x(i)**2
      write(6,607)d
607   format(' d = ',g12.4)

c     delta(i,j)=se(ii,j)*se(ii,i)*m(ia)
c     m is in g/mol

      allocate(qt(NQ))
      call vz(qt,NQ)
      do 1 i=1,NQ

      sx=0.0d0
      do 2 ii=1,N
      sx=sx+x(ii)*se(ii,i)*SFACR
      do 2 j=1,i-1
2     sx=sx-se(ii,i)*SFACR**2*se(ii,j)*qt(j)

      s2=0.0d0
      do 6 ii=1,N
6     s2=s2+se(ii,i)**2*SFACR**2
1     qt(i)=sx/s2

      d=0.0d0
      do 7 ii=1,N
      xt(ii)=re(ii)/bohr
      do 71 j=1,NQ
71    xt(ii)=xt(ii)+se(ii,j)*SFACR*qt(j)
      d=d+(rg(ii)/bohr-xt(ii))**2
7     xt(ii)=xt(ii)*bohr
      write(6,607)d
      open(44,file='test.x')
      write(44,*)'Excited orientation to match the ground:'
      write(44,*)nat
      do 13 ia=1,nat
13    write(44,608)z(ia),(xt(ii+3*(ia-1)),ii=1,3)
608   format(i4,3f12.6)
      close(44)

      open(9,file='ISTATE')
      write(9,600)
600   format(/,' mode       freq      Qt      n  nround:')
      do 3 i=1,NQ
      ani=ep(i)*qt(i)**2/2.0d0-0.50d0
3     write(9,601)i,ep(i)*CM,qt(i),ani,nint(ani+0.50d0)
 601  format(i4,f12.2,2G12.4,i6)
      close(9)
      write(6,*)'ISTATE written'
      return
      end
c     ============================================================
      subroutine shifte(nat,N,NQ,sg,se,eg,ep,m,lmscan,fac,
     1A,B,C,D,E,F,FI,G,GP,GH,GHP,J,JT,K,NB,NP,TOL,
     2lspec,del,e00,wsmin,wsmax,nps,ncut,LEXCL,
     1isu)
      implicit none
      integer*4 nat,N,NQ,i,ii,nt,ia,ix,NB,NP,nps,ncut,LEXCL,isu
      real*8 se(N,N),ep(*),CM,bohr,m(*),fac,TOL,EH,
     1del,e00,wsmin,wsmax,sg(N,N),eg(*)
      real*8 A(NQ,NQ),B(*),C(NQ,NQ),D(*),E(NQ,NQ),F(NQ,NQ),FI(NQ,NQ),
     1G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ),J(NQ,NQ),JT(NQ,NQ),K(*)
      real*8,allocatable::qt(:),re(:),rg(:),x(:),KT(:),EK(:),ani(:)
      integer*4,allocatable::z(:),nopt(:),noptnew(:)
      logical lmscan,lspec
      CM=219474.0d0
      bohr=0.529177d0
      allocate(qt(NQ),z(nat),re(N),rg(N),x(N),EK(NQ))
      call vz(EK,NQ)

c     x = xe - xg:
      call rdxx('ground/FILE.X',1,  nat,z,rg)
      call rdxx('exc.as.ground.x',1,nat,z,re)
      do 5 ii=1,N
5     x(ii)=(re(ii)-rg(ii))/bohr

c     qt(e) = -se^(-1) x:
      do 2 i=1,NQ
      qt(i)=0.0d0
      do 2 ia=1,nat
      do 2 ix=1,3
      ii=ix+3*(ia-1)
2     qt(i)=qt(i)-se(ii,i)*x(ii)*m(ia)*dsqrt(1822.89d0)

c     Energy = 0.5 sum(i) wi^2 qti^2
      EH=0.d0
      do 4 i=1,NQ
4     EH=EH+0.50d0*qt(i)**2*ep(i)**2
      write(6,6007)EH,EH*219474.0d0
6007  format('EH = ',f12.6,'au = ',f10.2,' cm-1')

      nt=0
      allocate(nopt(NQ),noptnew(NQ),ani(NQ))
      do 3 i=1,NQ
      ani(i)=0.5d0*(ep(i)*qt(i)**2-1.0d0)
      nopt(i)=max(0,nint(ani(i)+0.50d0))
      noptnew(i)=nopt(i)
3     nt=nt+nopt(i)
      write(6,*)nt,' - excited'

      if(ncut.gt.0)then
       allocate(KT(NQ))
       do 12 i=1,NQ
12     KT(i)=K(i)
       do 1 i=1,NQ
       if(nopt(i).gt.ncut)then
c       redefinition of excited normal mode Qk -> Qk
        write(6,6008)i,nopt(i)
6008    format(i4,':',i4,' -> 0')
        EK(i)=dble(nopt(i))*ep(i)
        noptnew(i)=0
        do 11 ii=1,NQ
11      KT(ii)=KT(ii)+J(ii,i)*qt(i)
       endif
1      continue
       call dofac(NQ,fac,J,JT,F,FI,G,GP,KT)
       write(6,300)fac,fac**2*100.0d0
300    format('<0|0*p> = ',g12.4,' (',f6.2,'%)')
       call newbd(NQ,G,KT,J,JT,FI,GH,GHP,B,D)
       call wdusch(NQ,N,A,B,C,D,E,JT,KT,eg,ep,sg,se,rg,re,m)
      endif

      open(9,file='PSTATE')
      write(9,600)
600   format(/,' mode       freq      Qt      ne  nround nnew  EK:')
      do 31 i=1,NQ
31    write(9,601)i,ep(i)*CM,qt(i),ani(i),nopt(i),noptnew(i),EK(i)
 601  format(i4,f12.2,2G12.4,2i6,g14.6)
      write(9,601)NB
      close(9)
      write(6,*)'PSTATE written'

      if(lmscan)call mscan(noptnew,NQ,fac,C,D,NB,NP,TOL,
     1lspec,del,e00,wsmin,wsmax,nps,ep,EK,LEXCL,isu)

      return
      end
c     ============================================================
      subroutine mscan(so,NQ,fac,C,D,NB,NP,TOL,
     1lspec,del,e00,wsmin,wsmax,nps,ep,EK,LEXCL,isu)
      implicit none
      real*8 fac,TOL,r,ojmax,del,e00,wsmin,wsmax,ep(*),
     1C(NQ,NQ),D(NQ),of,EK(*)
      integer*4 Ni,Nj,so(*),NQ,i,np0,nstart,ii,j,
     1kk,nnmax,jj,nend,NB,NP,nps,Nm,LEXCL,isu
      logical lspec
      integer*4,allocatable::si(:),sj(:),st(:),nbi(:),nn(:,:)
      real*8,allocatable::ob(:),oj(:)

      np0=1000000
      write(6,*)'scan called'

      write(6,*)'mother state:'
      do 3 kk=1,NQ
3     if(so(kk).ne.0)write(6,607)kk,so(kk)
607   format(i4,':',i2,$)
      write(6,*)
      Nm=les(NQ,so)
      write(6,*)Nm,' excited'
      allocate(sj(Nm+1))
      call viz(sj,Nm+1)
      call puts(NQ,so,sj)
      of=FC(fac,sj,Nm,np0,C,D,NQ)
      deallocate(sj)
      write(6,6093)of,tol
6093  format(' < 0 | mo> = ',e12.4,' tol: ',g12.4)
      write(6,*)


      Ni=0
      allocate(si(Ni+1),st(NQ),nbi(NQ),nn(NQ,NB),ob(NB),oj(2*NP+1))
      call viz(si,Ni+1)

      open(9,file='PSTATE.TAB')
      open(91,file='PROFILES.TXT')
      do 1 i=1,NQ
c     reference state to st:
      nstart=max(0,so(i)-NP)
      do 11 ii=-NP,NP
c     <0|j>:
      oj(ii+NP+1)=0.0d0
      if(ii+so(i).ge.nstart)then
       do 10 j=1,NQ
10     st(j)=so(j)
       st(i)=so(i)+ii
c      put to short record:
       Nj=les(NQ,st)
       allocate(sj(Nj+1))
       call viz(sj,Nj+1)
       call puts(NQ,st,sj)
       oj(ii+NP+1)=FC(fac,sj,Nj,np0,C,D,NQ)
       deallocate(sj)
      endif
11    continue
c
c     select maximum N:
      nnmax=1
      ojmax=0.0d0
      do 12 ii=-NP,NP
      if(dabs(oj(ii+NP+1)).gt.ojmax)then
       ojmax=dabs(oj(ii+NP+1))
       nnmax=so(i)+ii
      endif
12    continue

c     buffer with quantum number giving the best overlaps:
      nbi(i)=0
      do 13 ii=-NP,NP
      if(ojmax.gt.0.0d0)then
       r=dabs(oj(ii+NP+1))/ojmax
       if(r.gt.tol)then
        if(nbi(i).eq.0)then
         nbi(i)=nbi(i)+1
         nn(i,nbi(i))=so(i)+ii
         ob(nbi(i))=r
        else
         do 14 jj=1,nbi(i)
14       if(r.gt.ob(jj))goto 141
         if(nbi(i).lt.NB)then
          nbi(i)=nbi(i)+1
          nn(i,nbi(i))=so(i)+ii
          ob(nbi(i))=r
         endif
         goto 151
141      nend=min(NB,nbi(i)+1)
         do 15 kk=nend,jj+1,-1
         ob(kk)=ob(kk-1)
15       nn(i,kk)=nn(i,kk-1)
         ob(jj)=r
         nn(i,jj)=so(i)+ii
         nbi(i)=nend   
151      continue
        endif
       endif
      endif
13    continue

      write(9 ,601)i,nnmax,nbi(i),(nn(i,kk),kk=1,nbi(i))
601   format(100i4)
1     write(91,691)(oj(ii+NP+1),ii=-NP,NP)
691   format(100E10.2)
      close(9)
      write(6,*)'PSTATE.TAB written'

      
      if(isu.eq.1)call statesum1(fac,NQ,C,D,
     1lspec,del,e00,wsmin,wsmax,nps,ep,EK,LEXCL)

      if(isu.eq.2)call statesum2(fac,NQ,np0,C,D,NB,nbi,nn,
     1lspec,del,e00,wsmin,wsmax,nps,ep,EK)
      return
      end

c     ==============================================================
      subroutine viz(v,n)
      integer*4 v(*),i,n
      do 1  i=1,n
1     v(i)=0
      return 
      end
c     ==============================================================
      subroutine  puts(NQ,sd,sj)
c     rewrite vib state from long to compact notation:
      integer*4 NQ,sd(*),sj(*),kk,i,ii
      kk=0
      do 21 i=1,NQ
      do 21 ii=1,sd(i)
      kk=kk+1
21    sj(kk)=i
      return
      end
c     ==============================================================
      function les(N,s)
      integer*4 les,N,s(*),u,i
      u=0
      do 1 i=1,N
1     u=u+s(i)
      les=u
      return
      end
c     ==============================================================
      function FC(fac,si,Nexc,np0,C,D,N)
c     FC = Franc Condon factor for <0|si>
c     fac = <0|0*>
c     np0 ... dimension of working buffer
      implicit none
      integer*4 si(*),Nexc,np,iex,ii,nu,jj,ip,N,ic,jold,
     1nuj,jc,kk,np0,iq,jq
      real*8 FC,fac,D(*),C(N,N),pini
      real*8,allocatable::p(:)
      integer*4,allocatable::NN(:),ie(:,:),it(:)
      np=1
      allocate(p(np0),NN(np0),ie(np0,Nexc),it(np0))
      p(np)=1.0d0
      NN(np)=Nexc
      do 102 iex=1,NN(np)
102   ie(np,iex)=si(iex)
c     expand reccurent formula into a sum and shrink back to <0|0>
777   do 101 ii=1,np
c     term ii - reduce excitation one by one:
      pini=p(ii)
      do 101 iex=1,NN(ii)
      iq=ie(ii,iex)
      if(iq.gt.0)then
c      <0| si_nu> reduce to <0| i-1>, <0|i-2>, <0|i-1,j-1>(j<>i)
       nu=0
       do 1032 jj=1,NN(ii)
1032   if(ie(ii,jj).eq.iq)nu=nu+1

c      write term <0|v'-1_iex> into  string it:
       ip=0
       ic=0
       do 1031 jj=1,NN(ii)
       if(iq.eq.ie(ii,jj).and.ic.lt.1)then
        ic=ic+1
       else
        ip=ip+1
        it(ip)=ie(ii,jj)
       endif
1031   continue
c       
       call digest(np,ie,np0,Nexc,p,NN,it,NN(ii)-1,
     1 pini*D(iq)/dsqrt(dble(2*nu)))

       if(nu.gt.1)then
c       write term <0|v'-2_iex> into string it:
        ic=0
        ip=0
        do 1033 jj=1,NN(ii)
        if(ie(ii,jj).eq.iq.and.ic.lt.2)then
         ic=ic+1
        else
         ip=ip+1
         it(ip)=ie(ii,jj)
        endif
1033    continue
        call digest(np,ie,np0,Nexc,p,NN,it,NN(ii)-2,
     1  pini*dsqrt(dble(nu-1)/dble(nu))*C(iq,iq))
       endif

       jold=0
       do 106 jj=1,NN(ii)
       jq=ie(ii,jj)
       if(jq.ne.iq.and.jq.ne.jold)then
        nuj=0
        do 1034 kk=1,NN(ii)
1034    if(jq.eq.ie(ii,kk))nuj=nuj+1
c       write term <0|v'-1_iex-1_j, j<>iex> into it:
        ic=0
        jc=0
        ip=0
        do 1035 kk=1,NN(ii)
        if(ie(ii,kk).eq.iq.and.ic.lt.1)then
         ic=ic+1
        else
         if(ie(ii,kk).eq.jq.and.jc.lt.1)then
          jc=jc+1
         else
          ip=ip+1
          it(ip)=ie(ii,kk)
         endif
        endif
1035    continue
        call digest(np,ie,np0,Nexc,p,NN,it,NN(ii)-2,
     1  pini*dsqrt(dble(nuj)/dble(nu))*C(iq,jq))
       endif
106    jold=jq

c      eliminate the old term and start over
       do 105 jj=ii,np-1
       p(jj)=p(jj+1)
       NN(jj)=NN(jj+1)
       do 105 kk=1,NN(jj)
105    ie(jj,kk)=ie(jj+1,kk)
       np=np-1

       goto 777
       
      endif
101   continue

      if(np.ne.1)call report('np <> 1')

      FC=p(1)*fac
      return
      end
c     ==============================================================
      subroutine digest(np,ie,np0,LEXCL,p,NN,it,iexc,T)
      implicit none
      integer*4 je,np,np0,LEXCL,NN(*),it(*),ie(np0,LEXCL),jj,iexc
      real*8 T,p(*)

c     does it exist already within 1 ... np?:
      je=jeje(np,iexc,it,NN,ie,np0,LEXCL)

      if(je.ne.0)then
c      the term already exist in the expansion as je^th - just add it
       p(je)=p(je)+T
      else
c      term not found - add as new term
       np=np+1

       if(np.gt.np0)call report('Too many terms')

       p(np)=T
       NN(np)=iexc
       do 1 jj=1,iexc
1      ie(np,jj)=it(jj)
      endif
      return
      end
c     ==============================================================
      function jeje(np,nx,it,NN,ee,np0,LEXCL)
      implicit none
      integer*4 np,nx,it(*),NN(*),np0,LEXCL,ee(np0,LEXCL),je,jeje,jj,ie
      je=0
      do 104 jj=1,np
      if(nx.ne.NN(jj))goto 104
      do 1041 ie=1,NN(jj)
1041  if(it(ie).ne.ee(jj,ie))goto 104
      je=jj
      goto 1042
104   continue
1042  jeje=je
      return
      end
c     ==============================================================
      subroutine statesum2(fac,NQ,np0,C,D,NB,nbi,nn,
     1lspec,del,e00,wsmin,wsmax,nps,ep,EK)
      implicit none
      integer*4 NQ,ii,NB,nbi(*),LEXP,nn(NQ,NB),np0,nps
      real*8 fac,proc,procold,del,e00,wsmin,wsmax,CM,dx,EJ,
     1C(NQ,NQ),D(NQ),y,ep(*),ECM,emin,emax,EK(*),estart,p8
      logical lspec
      integer*4,allocatable::sj(:),st(:),is(:)
      integer*8 nt,i8
      real*8,allocatable::s(:)
      allocate(st(NQ),is(NQ))
      call viz(is,NQ)

      estart=e00
      do 100 ii=1,NQ
100   estart=estart+EK(ii)

      proc=0.0d0
      CM=219474.630d0
      procold=0.0d0

      if(lspec)then
       emin=1.0d99
       emax=0.0d0
       allocate(s(nps))
       call vz(s,nps)
       dx=(wsmax-wsmin)/dble(nps-1)
      endif

      nt=1
      do 1 ii=1,NQ
1     nt=int(nt,8)*int(nbi(ii),8)
      write(6,*)nt,' terms'

      do 2 i8=1,nt
      if(proc.gt.procold+0.1d0)then
       p8=100.0d0*dble(i8)/dble(nt)
       write(6,6000)i8,p8,proc
6000   format(i20,f9.4,'%',f9.4,'%')
       procold=proc
      endif
      call getis(i8,is,st,NQ,nn,NB,nbi)
      
c     rewrite final state st into short form in sj:
      LEXP=les(NQ,st)
      allocate(sj(LEXP+1))
      call puts(NQ,st,sj)
c

c     <0|*>
      y=FC(fac,sj,LEXP,np0,C,D,NQ)**2
      proc=proc+100.0d0*y

      if(lspec)then
c      energy:
       EJ=estart
       DO 1002 ii=1,LEXP
1002   EJ=EJ+ep(sj(ii))
       ECM=EJ*CM
       if(ECM.lt.emin)emin=ECM
       if(ECM.gt.emax)emax=ECM
       call spread(nps,y,dx,wsmin,ECM,del,s)
      endif

2     deallocate(sj)

      if(lspec)then
       call ws(s,wsmin,dx,nps)
       write(6,604)emin,emax
604    format(' emin: ',f12.2,' emax: ',f12.2,' cm-1')
      endif


      return
      end
c     ==============================================================
      subroutine spread(nps,y,dx,wmin,E,del,s)
      implicit none
      real*8 y,dx,wmin,E,del,s(*),x
      integer*4 nps,i,nc,n3,nstart,nend
c     closest point, approx:
      nc=nint((E-wmin)/dx)
c     points over three delta:
      n3=nint(3.0d0*del/dx)
      nstart=max(nc-n3,1)
      nend  =min(nc+n3,nps)
      x=wmin+dx*dble(nstart-2)
      do 1 i=nstart,nend
      x=x+dx
1     s(i)=s(i)+y*exp(-((E-x)/del)**2)
      return
      end 
c     ==============================================================
      subroutine ws(s,wmin,dx,nps)
      implicit none
      real*8 s(*),wmin,dx,x,y
      integer*4 nps,i
      open(9,file='S.PRN')
      x=wmin-dx
      do 1 i=1,nps
      x=x+dx
      y=s(i)
      if(dabs(y).lt.1.0d-10)y=0.0d0
 1    write(9,90)x,y
 90   format(f10.2,e12.4)
      close(9)
      return
      end
c     ==============================================================
      subroutine dofac(NQ,fac,J,JT,F,FI,G,GP,K)
      implicit none
      integer*4 NQ,sgn,sgn2,info
      real*8 fac,G(NQ,NQ),GP(NQ,NQ),K(NQ),lT,lF,lQ,sp1,sp2,
     1 F(NQ,NQ),J(NQ,NQ),JT(NQ,NQ),FI(NQ,NQ),lJ,f1,gn,gpn,tn
      real*8,allocatable:: T(:,:),TV(:),TU(:)
      real*8 detT,detJ,detF,e1(1,1),e2(1,1),G_c(NQ,NQ),Gp_c(NQ,NQ),
     1 F_new(NQ,NQ),FI_new(NQ,NQ),help(NQ,2*NQ),tol,K_m(NQ,1)
      real*8 fac2
      
      allocate(T(NQ,NQ),TV(NQ),TU(NQ))
      
      detJ=DetFromLU(J,NQ)
      write(6,'(A,F8.4)')'|J| = ',detJ
      if(abs(detJ-1.0) > 0.05)write(6,*)
     1 'Warning: Determinant of J is different from one'
     
      goto 10
      
      !alternate route for <0|0>
      gn=NormalizeW(G,NQ)
      gpn=NormalizeW(Gp,NQ)
      G_c=G/gn
      Gp_c=Gp/gpn
      T=matmul(G_c,Gp_c)
      detT=(DiagProd(T,NQ))**(1d0/4d0)*gn**(NQ/4d0)*gpn**(NQ/4d0)
      F_new=matmul(JT,matmul(G,J))+Gp
      detF=DetFromLU(F_new,NQ)
      K_m(:,1)=K
      e1=-0.5d0*matmul(transpose(K_m),matmul(G,K_m))
      tol=1d-10
      !F_new=matmul(JT,matmul(G,J))+Gp_c
      call INV(NQ,F_new,FI_new,NQ,tol,help,info)
      e2=0.5d0*matmul(transpose(K_m),matmul(G,matmul(J,matmul(FI_new,
     1 matmul(JT,matmul(G,K_m))))))
      
      
      fac2=2d0**(dble(NQ)/2)*detT*sqrt(abs(detJ)/detF)*
     1 exp(e1(1,1)+e2(1,1))
      fac=fac2
      
10    sgn2=0
      call mm(T,G,Gp,NQ)
      lT=LOGDET(T,NQ,NQ,sgn)
      sgn2=sgn+sgn2
      lF=LOGDET(F,NQ,NQ,sgn)
      sgn2=sgn+sgn2
      lJ=LOGDET(J,NQ,NQ,sgn)
      sgn2=sgn+sgn2
      lQ=dble(NQ)*log(2.0d0)
      f1=exp((lQ+lt/2.0d0+lJ-lF)/2.0d0)
      call mv(TV,G,K,NQ)
      sp1=sp(K,TV,NQ)
      call mv(TU,G,K,NQ)
      call mv(TV,JT,TU,NQ)
      call mv(TU,FI,TV,NQ)
      call mv(TV,J,TU,NQ)
      call mv(TU,G,TV,NQ)
      sp2=sp(K,TU,NQ)
      fac=f1*exp(-0.50d0*(sp1-sp2))
      return
      end
      
      function DetFromLU(mat,N)result(det)
         implicit none
         integer*4 N,piv(N),info,detP
         real*8 mat(N,N),det
         real*8 mat_c(N,N)
         
         mat_c=mat
         
         call dgetrf(N,N,mat_c,N,Piv,info)
         if(info/=0)then
            stop 2
         end if
         
         call DetPermMat(Piv,N,detP)
         det=detP*DiagProd(mat_c,N)
      end function DetFromLU
      
      !diagonal product of a matrix
      !ie. the determinant of diagonal/upper tr./lower tr. matrix
      function DiagProd(mat,N)result(prod)
         implicit none
         integer*4 N,i
         double precision mat(N,N),prod
         
         prod=1
         do i = 1,N
            prod=prod*mat(i,i)
         end do
      end function DiagProd
      
      Function NormalizeW(W,NQ)result(prod)
         implicit none
         integer NQ,i
         double precision W(NQ,NQ),prod
         
         prod=1d0
         do i = 1,NQ
            prod=prod*W(i,i)**(1d0/dble(NQ))
         end do
      end function NormalizeW
      
      subroutine DetPermMat(Perm,N,det)
         implicit none
         integer*4 N,Perm(N),det
         integer*4 i,j,buf
         
         det=1
         do i = 1,N
            if(Perm(i)==i)then
               cycle
            end if
            
            do j = i+1,N
               if(Perm(j)==i)then
                  buf=Perm(i)
                  Perm(i)=Perm(j)
                  Perm(j)=buf
                  det=-det
                  exit
               end if
            end do
         end do
      end subroutine DetPermMat
      
c     ==============================================================
      subroutine newbd(NQ,G,K,J,JT,FI,GH,GHP,B,D)
      implicit none
      integer*4 NQ,ix
      real*8 G(NQ,NQ),K(*),JT(NQ,NQ),FI(NQ,NQ),GHP(NQ,NQ),D(*),
     1B(*),GH(NQ,NQ),J(NQ,NQ)
      real*8,allocatable::TU(:),TV(:),T(:,:),TP(:,:)
      allocate(TV(NQ),TU(NQ),T(NQ,NQ),TP(NQ,NQ))
      call mv(TV,G,K,NQ)
      call mv(TU,JT,TV,NQ)
      call mv(TV,FI,TU,NQ)
      call mv(TU,GHP,TV,NQ)
      do 17 ix=1,NQ
17    D(ix)=-2.0d0*TU(ix)
c     D= - 2 GP^1/2 Fi Jt G K

      call mm(T,JT,G,NQ)
      call mm(TP,FI,T,NQ)
      call mm(T,J,TP,NQ)
      do 13 ix=1,NQ
13    T(ix,ix)=T(ix,ix)-1.0d0
      call mv(TV,T,K,NQ)
      call mv(TU,GH,TV,NQ)
      do 14 ix=1,NQ
14    B(ix)=-2.0d0*TU(ix)
c     B= - 2 G^1/2 ( J Fi Jt G - I) K
      return
      end
c     ==============================================================
      subroutine statesum1(fac,NQ,C,D,
     1lspec,del,e00,wsmin,wsmax,nps,ep,EK,LEXCL)
      implicit none
      integer*4 NQ,ii,NB,nm1,nm2,i,ic,iex,iq,jq,Nexc,
     1nps,LEXCL,nne,nu,nuj
      real*8 fac,proc,procold,del,e00,wsmin,wsmax,CM,dx,EJ,
     1C(NQ,NQ),D(NQ),y,ep(*),ECM,emin,emax,EK(*),estart,fci
      logical lspec
      real*8,allocatable::s(:),fcs(:),fcm1(:),fcm2(:)
      integer*4,allocatable::si(:),sls(:,:),slm1(:,:),slm2(:,:),sl(:)

      estart=e00
      do 100 ii=1,NQ
100   estart=estart+EK(ii)

      proc=0.0d0
      CM=219474.630d0
      procold=0.0d0

      call statereport(NQ,LEXCL)

      if(lspec)then
       emin=1.0d99
       emax=0.0d0
       allocate(s(nps))
       call vz(s,nps)
       dx=(wsmax-wsmin)/dble(nps-1)
      endif

c     LEXCL: maximal number of excitations
      nb=NQ**LEXCL
      allocate(si(LEXCL+1),fcm1(nb),fcm2(nb),slm1(nb,NQ),
     1slm2(nb,NQ),fcs(nb),sls(nb,NQ),sl(NQ))

      si(1)=0
c     lc:current number of excitations
c     <0|0*>:
      y=fac**2
      proc=100.0d0*y

      if(lspec)then
c      energy:
       EJ=estart
       ECM=EJ*CM
       if(ECM.lt.emin)emin=ECM
       if(ECM.gt.emax)emax=ECM
       call spread(nps,y,dx,wsmin,ECM,del,s)
      endif
      nm2=0
      nm1=0
      nne=1
      do 77 ii=1,NQ
77    sls(1,ii)=0
      fcs(1)=fac
      fcm1(1)=0
      slm1(1,1)=0

      do 30000 Nexc=1,LEXCL
      write(6,6000)NExc,proc
6000  format(i2,f9.4,'%')
c     rewrite strings:
      do 71 ii=1,nm1
      fcm2(ii)=fcm1(ii)
      do 71 iq=1,NQ
71    slm2(ii,iq)=slm1(ii,iq)
      nm2=nm1
      do 81 ii=1,nne
      fcm1(ii)=fcs(ii)
      do 81 iq=1,NQ
81    slm1(ii,iq)=sls(ii,iq)
      nm1=nne
c     number of states Nexc excited
      nne=0
c     distribute Nexc excitations upon centers 1..NQ:
c     initial indices - all to the first center:
      do 41 iex=1,Nexc
41    si(iex)=1
50000 nne=nne+1
c     transcript to long:
      call trl(si,Nexc,sl,NQ)
      fci=0.0d0
c     reduce first excitation on center iq
      iq=si(1)
      nu=sl(iq)
c     from <0|i-1>, find it in previous list:
      do 44 ii=1,nm1
      do 441 jq=1,NQ
      if(jq.ne.iq)then
       if(sl(jq).ne.slm1(ii,jq))goto 44
      else
       if(nu.ne.slm1(ii,jq)+1)goto 44
      endif
441   continue
      fci=fci+D(iq)*fcm1(ii)/dsqrt(dble(2*nu))
      goto 442
44    continue
442   continue
c     from <0|i-2>, find it in pre-previous list:
      if(nu.gt.1)then
       do 54 ii=1,nm2
       do 541 jq=1,NQ
       if(jq.ne.iq)then
        if(sl(jq).ne.slm2(ii,jq))goto 54
       else
        if(nu.ne.slm2(ii,jq)+2)goto 54
       endif
541    continue
       fci=fci+C(iq,iq)*fcm2(ii)*dsqrt(dble(nu-1)/dble(nu))
       goto 542
54     continue
542    continue
      endif
c     from <0|i-1 j-1>
      do 64 ii=1,nm2
      nuj=0
      do 641 jq=1,NQ
      if(jq.ne.iq)then
       if(sl(jq).ne.slm2(ii,jq))then
        if(nuj.eq.0.and.sl(jq).eq.slm2(ii,jq)+1)then
         nuj=sl(jq)
         goto 641
        else
         goto 64
        endif
       else
        goto 64
       endif
      else
       if(nu.ne.slm2(ii,jq)+1)goto 64
      endif
641   continue
      fci=fci+C(iq,jq)*fcm2(ii)*dsqrt(dble(nuj)/dble(nu))
      goto 642
64    continue
642   continue
c     save for later:
      do 42 ii=1,NQ
42    sls(nne,ii)=sl(ii)
      fcs(nne)=fci
c
c     <0|*>
      y=fci**2
      proc=proc+100.0d0*y
      if(proc.gt.procold+1.0d0)then
       write(6,6000)Nexc,proc
       procold=proc
      endif
      if(lspec)then
c      energy:
       EJ=estart
       DO 1002 ii=1,Nexc
1002   EJ=EJ+ep(si(ii))
       ECM=EJ*CM
       if(ECM.lt.emin)emin=ECM
       if(ECM.gt.emax)emax=ECM
       call spread(nps,y,dx,wsmin,ECM,del,s)
      endif
c
c     find index to be changed
      do 80000 ic=Nexc,1,-1
80000 if(si(ic).lt.NQ)goto 90000
      goto 30000
90000 do 10000 i=ic+1,Nexc
10000 si(i)=si(ic)+1
      si(ic)=si(ic)+1
      goto 50000
30000 continue

      if(lspec)then
       call ws(s,wsmin,dx,nps)
       write(6,604)emin,emax
604    format(' emin: ',f12.2,' emax: ',f12.2,' cm-1')
      endif

      return
      end
c     ==============================================================
      subroutine trl(si,Nexc,sl,NQ)
      integer*4 si(*),Nexc,sl(*),NQ,i
      do 1 i=1,NQ
1     sl(i)=0
      do 2 i=1,Nexc
2     sl(si(i))=sl(si(i))+1
      return
      end
c     ==============================================================
      subroutine statereport(N,LEXCL)
      implicit none
      integer*4 N,LEXCL
      integer*8 i,sum
      integer*8 ,allocatable::NS(:)
      allocate(NS(LEXCL+1))
      NS(1)=1
      do 1 i=1,LEXCL
1     NS(i+1)=(NS(i)*(N+i-1))/i
      WRITE(6,4012)0,NS(1),NS(1)
4012  FORMAT(' Number of states excited',i3,'-times: ',2i20)
      sum=1
      do 2 i=1,LEXCL
      sum=sum+NS(i+1)
2     WRITE(6,3012)i,NS(i+1),sum
3012  FORMAT('                         ',i3,'-times: ',2i20)
      WRITE(6,5012)sum
5012  FORMAT('                           -------- ',/,
     1       '                             total: ',i20)
      return
      end
c     ==============================================================
      subroutine wdusch(NQ,N,A,B,C,D,E,JT,K,eg,ep,sg,se,rg,re,m)
      implicit none
      integer*4 NQ,N,iq,ix,ia
      real*8 A(NQ,NQ),B(*),C(NQ,NQ),D(*),E(NQ,NQ),JT(NQ,NQ),m(*),
     1eg(*),ep(*),K(*),sg(N,N),se(N,N),rg(*),re(*),SFAC,bohr
      real*8,allocatable::ci(:),cj(:)
      SFAC=0.0234280d0
      bohr=0.529177d0
      allocate(ci(NQ),cj(NQ))
      open(10,FILE='DUSCH.OUT')
      write(10,*)NQ
      call wdm(10,NQ,JT,'Duschinsky Matrix')
      call wdv(10,NQ,K,'Shift Vector')
      call wdm(10,NQ,A,'A Matrix')
      call wdv(10,NQ,B,'B Vector')
      call wdm(10,NQ,C,'C Matrix')
      call wdv(10,NQ,D,'D Vector')
      call wdm(10,NQ,E,'E Matrix')
      call wwm(10,NQ,eg,'Ground w')
      call wwm(10,NQ,ep,'Excited w')
      call wds(10,N,NQ,sg,'Ground S-Matrix')
      call wds(10,N,NQ,se,'Excited S-Matrix')
      do 20 iq=1,NQ
      ci(iq)=0.0d0
      cj(iq)=0.0d0
      do 20 ix=1,N
      ia=(ix+2)/3
      cj(iq)=cj(iq)+(rg(ix)-re(ix))/bohr
     1*se(ix,iq)*SFAC*m(ia)*1822.89d0
20    ci(iq)=ci(iq)+(rg(ix)-re(ix))/bohr
     1*sg(ix,iq)*SFAC*m(ia)*1822.89d0
      write(10,*)
      write(10,*)' mode  SG*dX   SE*dX  ng    ne'
      do 21 iq=1,NQ
21    write(10,621)iq,ci(iq),cj(iq),ci(iq)**2*eg(iq)/2.0d0,
     1cj(iq)**2*ep(iq)/2.0d0
621   format(i5,2f12.6,2f8.1)
      close(10)
      write(6,*)' DUSCH.OUT written'
      return
      end
c     ==============================================================
      subroutine getis(i8,is,st,NQ,nn,NB,nbi)
      implicit none
      integer*4 is(*),st(*),NQ,NB,nn(NQ,NB),nbi(*),ii,k,kk
      integer*8 i8,isum,ifa,ifa2
      do 21 ii=NQ,1,-1
      isum=i8-int(1,8)
      do 221 k=1,NQ-ii
      ifa=int(nbi(1),8)
      do 223 kk=2,NQ-k
223   ifa=ifa*int(nbi(kk),8)
221   isum=isum-ifa*int(is(NQ-k+1)-1,8)
      ifa2=1
      do 222 k=1,ii-1
222   ifa2=ifa2*int(nbi(k),8)
21    is(ii)=int(isum/ifa2+int(1,8),4)
      do 224 k=1,NQ
224   st(k)=nn(k,is(k))
      return
      end
c     ==============================================================
      subroutine seproj(N,NF,NQ,se,sgf,m)
      implicit none
      integer*4 N,NF,NQ,i,j,ix
      real*8 se(:,:),sgf(:,:),norm,m(*),sij
c     each mode in se:
      do 1 i=1,NQ
c     each deleted mode:
      do 2 j=1,NF
      sij=spm(sgf,se,j,i,N,NQ,m)
      do 2 ix=1,N
2     se(ix,i)=se(ix,i)-sij*sgf(ix,j)
c     previous modes
      do 3 j=1,i-1
      sij=spm(se,se,j,i,N,NQ,m)
      do 3 ix=1,N
3     se(ix,i)=se(ix,i)-sij*se(ix,j)
c     renormalize:
      norm=1.0d0/dsqrt(spm(se,se,i,i,N,NQ,m))
      do 1 ix=1,N
1     se(ix,i)=se(ix,i)*norm
      write(6,*)NF,' modes projected from se'
      return
      end
c     ==============================================================
      function spm(s1,s2,i,j,N,NQ,m)
      implicit none
      integer*4 i,j,N,ix,ia,nq
      real*8 s1(:,:),s2(:,:),a,m(*),spm
      a=0.0d0
      do 1 ia=1,N/3
      do 1 ix=1,3
1     a=a+s1(ix+3*(ia-1),i)*s2(ix+3*(ia-1),j)*m(ia)
      spm=a
      return
      end
c     ==============================================================
      subroutine zzz(NQ,JT,K,G,GP,GH,GHP,J,T,TP,F,FI,E2,A,B,TV,TU,C,D,
     1E,K0,vh)
      implicit none
      integer*4 NQ
      real*8 JT(NQ,NQ),K(NQ),G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ),
     1J(NQ,NQ),T(NQ,NQ),TP(NQ,NQ),F(NQ,NQ),FI(NQ,NQ),E2(NQ,2*NQ),
     1A(NQ,NQ),B(NQ),TV(NQ),TU(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ),K0(NQ)
      logical vh
      A=0.0d0
      B=0.0d0
      C=0.0d0
      D=0.0d0
      E=0.0d0
      E2=0.0d0
      F=0.0d0
      FI=0.0d0
      G=0.0d0
      GH=0.0d0
      GHP=0.0d0
      GP=0.0d0
      J=0.0d0
      if(.not.VH)JT=0.0d0
      if(.not.VH)K=0.0d0
      K0=0.0d0
      T=0.0d0
      TP=0.0d0
      TV=0.0d0
      TU=0.0d0
      return
      end
      
      subroutine rpl(N,NQ,nat,m,sg,se)
      implicit none
      integer*4  N,NQ,nat,ie,ig,ix,ii,ia
      real*8 m(nat),sg(N,N),se(N,N),cl,sp
      integer*4 ,allocatable::ige(:),it(:)
      allocate(ige(NQ),it(NQ))
      write(6,600)
600   format(' excited mode / closest ground product')
      it=0
      ie=1
111   ige(ie)=0    
      cl=0.0d0
      do 2 ig=1,NQ
      if(it(ig).eq.0)then
       sp=0.0d0
       do 21 ia=1,nat
       do 21 ix=1,3
       ii=ix+3*(ia-1)
21     sp=sp+sg(ii,ig)*m(ia)*se(ii,ie)
       if(dabs(sp).gt.cl)then
        ige(ie)=ig
        cl=dabs(sp)
       endif
      endif
2     continue
      if(ige(ie).eq.0)call report('mode cannot be assigned')
      it(ige(ie))=ie
      write(6,601)ie,ige(ie),cl
601   format(i5,' /',i5,f6.3,$)
      if(mod(ie,3).eq.0)write(6,*)
      do 1 ia=1,nat
      do 1 ix=1,3
      ii=ix+3*(ia-1)
1     se(ii,ie)=sg(ii,ige(ie))
      if(ie.lt.NQ)then
       ie=ie+1
       goto 111
      endif
      write(6,602)NQ
602   format(i6,' S-vectors replaced')
      return
      end

      subroutine ddm(N,i,NQ,e,s)
c     delete mode i from energies e and s-matrix s
      implicit none
      integer*4 N,i,NQ,j,ix,ii
      real*8,allocatable :: s(:,:),e(:)
C       real*8 :: s_help(N,NQ-1),e_help(NQ-1)
      
C       e_help=0
C       s_help=0
C       do ii = 1,i-1
C          e_help(ii)=e(ii)
C          s_help(:,ii)=s(:,ii)
C       end do
      do 10 j=i,NQ-1
      e(j)=e(j+1)
C       e_help(j)=e(j)
      do 10 ix=1,N
10    s(ix,j)=s(ix,j+1)
C 10    s_help(ix,j)=s(ix,j)
      return
      end

      function fcg(i,NQ,JT)
c     for excited i find corresponding ground index
c     in Duschinsky matrix
      implicit none
      integer*4 fcg,i,NQ,imx,j
      real*8 JT(NQ,NQ),smx
      imx=1
      smx=dabs(JT(i,1))
      do 1 j=2,NQ
      if(dabs(JT(i,j)).gt.smx)then
       smx=dabs(JT(i,j))
       imx=j
      endif
1     continue
      fcg=imx
      return
      end

      function fce(i,NQ,JT)
c     for ground i find corresponding excited index
c     in Duschinsky matrix
      implicit none
      integer*4 fce,i,NQ,imx,j
      real*8 JT(NQ,NQ),smx
      imx=1
      smx=dabs(JT(1,i))
      do 1 j=2,NQ
      if(dabs(JT(j,i)).gt.smx)then
       smx=dabs(JT(j,i))
       imx=j
      endif
1     continue
      fce=imx
      return
      end

      subroutine grad2Q(nat,NQ,gr,gr_q,M,M_invsq,L,amu2au)
         implicit none
         integer*4 nat,i,ix,NQ
         real*8 gr(3*nat),gr_q(NQ),M(nat),M_invsq(3*nat,3*nat),
     1   L(3*nat,NQ),amu2au!,hold(3*nat)
         
         M_invsq=0d0
         do i = 1,nat
            ix=(i-1)*3+1
            M_invsq(ix,ix)=1d0/sqrt(M(i)*amu2au)
            M_invsq(ix+1,ix+1)=M_invsq(ix,ix)
            M_invsq(ix+2,ix+2)=M_invsq(ix,ix)
         end do
         
         
         gr_q=matmul(transpose(L),matmul(M_invsq,gr))
         !call mv(hold,L,gr,3*nat)
         !call mv(gr_q,m_invsq,hold,3*nat)
      end subroutine grad2Q
      
      subroutine readgrad_dusch(nat,gr,filee)
         implicit none
         integer*4 nat,i,ix
         real*8 gr(3*nat)
         character(*) filee
         
         open(77,file=filee)
         do i = 1,nat
            ix=(i-1)*3+1
            read(77,*)gr(ix),gr(ix+1),gr(ix+2)
         end do
         close(77)
      end subroutine readgrad_dusch
      
     
      subroutine FFc2q(nat,nq,ffc,ffq,M_invsq,L)
         integer*4 nat,i,ix,nq
         real*8 ffc(3*nat,3*nat),ffq(nq,nq),
     1   M_invsq(3*nat,3*nat),L(3*nat,nq)
              
C          ffq=matmul(transpose(L),matmul(M_invsq,matmul(ffc,
C      1   matmul(L,M_invsq))))
         ffq=matmul(transpose(L),matmul(M_invsq,matmul(ffc,
     1   matmul(M_invsq,L))))
      end subroutine FFc2q
      
      !ripped from new4
      SUBROUTINE WRITEFF(N,FCAR)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FCAR(N,N)
      !CONST=4.359828d0/0.5291772d0**2
      CONST=1
      N1=1
1     N3=N1+4
      IF(N3.GT.N)N3=N
      DO 130 LN=N1,N
130   WRITE(20,17)LN,(FCAR(LN,J)/CONST,J=N1,MIN(LN,N3))
      N1=N1+5
      IF(N3.LT.N)GOTO 1
17    FORMAT(I4,5D14.6)
      RETURN
      END
      
      subroutine svgr(g,n,fille)
      implicit none
      character(*) fille
      real*8 g(*)
      integer*4 n,i,ix
      open(70,file=fille)
      do 1 i=1,n
1     write(70,700)(g(ix+3*(i-1)),ix=1,3)
700   format(3f15.9)
      close(70)
      write(6,*)'Gradient written to '//fille
      return
      end
      
      subroutine svgrq(g,nq,fille)
      implicit none
      character(*) fille
      real*8 g(*)
      integer*4 nq,i,ix
      open(70,file=fille)
      do 1 i=1,nq
1     write(70,700)g(i)
700   format(f15.9)
      close(70)
      write(6,*)'Gradient written to '//fille
      return
      end
      
      !ripped from new4
      ! SUBROUTINE READFF(N,FCAR,file_fc,file_x)
      ! IMPLICIT INTEGER*4 (I-N)
      ! IMPLICIT REAL*8 (A-H,O-Z)
      ! DIMENSION FCAR(N,N)
      ! CHARACTER(*) file_fc,file_x
      ! OPEN(20,FILE=file_fc,STATUS='OLD')
      ! N1=1
! 1     N3=N1+4
      ! IF(N3.GT.N)N3=N
      ! DO 130 LN=N1,N
! 130   READ(20,17)(FCAR(LN,J),J=N1,MIN(LN,N3))
      ! N1=N1+5
      ! IF(N3.LT.N)GOTO 1
! 17    FORMAT(4X,5D14.6)
! C
! c      CONST=4.359828d0/0.5291772d0/0.5291772d0 ! conversion to whatever
      ! CONST=1
      ! DO 3 I=1,N
      ! DO 3 J=I,N
! 3     FCAR(J,I)=FCAR(J,I)*CONST
      ! DO 31 I=1,N
      ! DO 31 J=I+1,N
   ! 31 FCAR(I,J)=FCAR(J,I)
      ! CLOSE(20)
      ! WRITE(*,*)' Cartesian FF read in ... '
      ! RETURN
      ! END
      
      SUBROUTINE TRERR(F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION F(N,N)
      err=0.0d0
      NAT=N/3
      do 1 I=1,N
      s1=F(I,1)
      s2=F(I,2)
      s3=F(I,3)
      do 2 ia=2,nat
      s1=s1+F(I,1+3*(ia-1))
      s2=s2+F(I,2+3*(ia-1))
2     s3=s3+F(I,3+3*(ia-1))
1     err=err+s1**2+s2**2+s3**2
      write(6,600)sqrt(err)
600   format(' Translational FF error ',f20.10)
      return
      end

      SUBROUTINE PROJECT(A,F,N,file_x)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N),F(N,N)
      real*8,allocatable::TEM(:,:)
      character(*) file_x
C
      call TRERR(F,N)
      allocate(TEM(N,N))
      WRITE(*,*)'Projecting transl/rotations from the force field...'
      CALL DOMA(A,N,file_x)
      WRITE(*,*)'Matrix done'
      DO 11 I=1,N
      DO 11 J=1,N
11    TEM(I,J)=0.0d0

      DO 1 I=1,N
      DO 1 II=1,N
      if(dabs(A(I,II)).gt.1.0d-10)then
       DO 12 J=1,N
12     TEM(I,J)=TEM(I,J)+A(I,II)*F(II,J)
      endif
1     continue

      DO 22 I=1,N
      DO 22 J=1,N
22    F(I,J)=0.0d0

      DO 2 J=1,N
      DO 2 II=1,N
      if(dabs(A(J,II)).gt.1.0d-10)then
       DO 23 I=1,N
23     F(I,J)=F(I,J)+TEM(I,II)*A(J,II)
      endif
2     continue
      WRITE(*,*)'  ... done.'
      call TRERR(F,N)
c      F=F/(4.359828d0/0.5291772d0/0.5291772d0)
      RETURN
      END

      SUBROUTINE ORT(A,N6,NAT)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(3*NAT,N6)
      AMACH=0.00000000000001d0
      NAT3=3*NAT
      C=0.0d0
      DO 1 I=1,NAT3
1     C=C+A(I,1)**2
      IF(C.GT.AMACH)C=1.0d0/SQRT(C)
      DO 2 I=1,NAT3
2     A(I,1)=A(I,1)*C
      DO 3 IC=2,N6
C     Column IC, orthonormalize to previous columns:
      DO 3 ICP=1,IC-1
      SPL=0.0d0
      DO 5 J=1,NAT3
5     SPL=SPL+A(J,IC)*A(J,ICP)
C     Subtract projection to ICP:
      DO 6 J=1,NAT3
6     A(J,IC)=A(J,IC)-SPL*A(J,ICP)
C     Normalize :
      SPL=0.0d0
      DO 7 J=1,NAT3
7     SPL=SPL+A(J,IC)**2
      IF(SPL.GT.AMACH)SPL=1.0d0/SQRT(SPL)
      DO 3 J=1,NAT3
3     A(J,IC)=A(J,IC)*SPL
      RETURN
      END

      SUBROUTINE DOMA(A,N,file_x)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N)
      real*8,allocatable::C(:,:),TEM(:,:)
      character(*) file_x
C
      N=N
      allocate(C(3,N/3),TEM(N,N))
      AMACH=0.00000000000001d0
      NAT=N/3
      OPEN(4,FILE=file_x,STATUS='OLD')
      READ(4,*)
      READ(4,*)NAT
      DO 11 I=1,NAT
11    READ(4,*)C(1,I),(C(II,I),II=1,3)
      CLOSE(4)
      DO 12 I=1,N
      DO 12 J=1,6
12    A(I,J)=0.0d0
      DO 1 I=1,NAT
      A(3*I-2,1)= 1.00d0
      A(3*I-1,2)= 1.00d0
      A(3*I  ,3)= 1.00d0
      A(3*I-2,4)=-C(2,I)
      A(3*I-1,4)= C(1,I)
      A(3*I-1,5)=-C(3,I)
      A(3*I  ,5)= C(2,I)
      A(3*I-2,6)= C(3,I)
1     A(3*I  ,6)=-C(1,I)
      N6=6
      CALL ORT(A,N6,NAT)
      ICOL=1
5     IF(ICOL.GT.N6)GOTO 8
      S=0.0d0
      DO 6 J=1,N
6     S=S+A(J,ICOL)**2
      IF(S.LT.1000.0d0*AMACH)THEN
       DO 7 J=1,N
7      A(J,ICOL)=A(J,N6)
       N6=N6-1
       GOTO 5
      ENDIF
      ICOL=ICOL+1
      GOTO 5
8     CONTINUE
      WRITE(*,*)N6,' coordinates projected'
      DO 2 I=1,N
      DO 2 J=1,N
      TEM(I,J)=0.0d0
      DO 2 II=1,N6
2     TEM(I,J)=TEM(I,J)+A(I,II)*A(J,II)
      DO 3 I=1,N
      DO 4 J=1,N
4     A(I,J)=-TEM(I,J)
3     A(I,I)=A(I,I)+1.0d0
      RETURN
      END
      
      
      subroutine eye(A,N)
         integer*4 N,i
         real*8 A(N,N)
         A=0d0
         do i = 1,N
            A(i,i)=1d0
         end do
      end subroutine eye
      
      end module legacy

      program dusch
c     Duschinsky transformation
      use legacy
      use util
      implicit none
      real*8 u(3,3),fac,TOL,facv,wlim,wfix,sslim,
     1 wmin,TOLN,del,e00,wsmin,wsmax,kfac,
     1 expGeom_fac,TOLL,e_vert,e_ex,e_gr,deltol,e_disp,e_disp2
      integer*4 ia,N,nat,NE,NQ,IERR,maxit,NB,NP,nps,ncut,LEXCL,isu,
     1 ifx,NF,arg_stat,nqn,kk,i,nqi_g,nqi_e,delk,ii,expGeom,ix,
     1 Gauss_dusch,iz,corenVH
      real*8,allocatable::rg(:),re(:),se(:,:),sg(:,:),JT(:,:),sgf(:,:),
     1eg(:),ep(:),K(:),m(:),F(:,:),G(:,:),GP(:,:),T(:,:),TP(:,:),
     1sef(:,:),FI(:,:),E2(:,:),GH(:,:),GHP(:,:),A(:,:),TV(:),TU(:),
     1B(:),J(:,:),C(:,:),D(:),E(:,:),ref(:),rgf(:),K0(:),
     1 sg_full(:,:),se_full(:,:),s_help(:,:),ffc(:,:),ffc_a(:,:),
     1 freqs(:,:),help_v(:),help_m(:,:),e_help(:,:),rexp(:),grad(:),
     1 help_se(:,:),help_sg(:,:),ffq_i(:,:),ge_q(:),w_new(:),
     1 eg_help(:),ep_help(:),r_disp(:),K_orig(:)
      double precision,allocatable :: sg_m(:,:),ffq(:,:),gr(:),gr_q(:),
     1 gr_exp(:),gr_exp_q(:),m_invsq(:,:),se_m(:,:),ep_orig(:)
      integer*4,allocatable::z(:),del_modes_gr(:),del_modes_ex(:),
     1 MO(:)
      real*8,parameter :: au2cm=1d0/4.556431403d-6
      real*8,parameter :: amu2au = 1822.89d0
      real*8,parameter :: bohr=0.529177d0

      character(10) argg
      character(4) hes_type
      character(80) Gauss_dusch_file
C       real amas(89)
C       data amas/1.007825,4.002603,
C      2  6.941, 9.012,   10.810,12.000,14.003074,15.9949,18.998,20.179,
C      3 22.990,24.305,   26.981,28.086,30.974,32.060,35.453,39.948,
C      4 39.098,40.080,44.956,47.900,50.941,51.996,54.938,55.847,
C      4 58.933,58.700,63.546,65.380,
C      4                  69.720,72.590,74.922,78.960,79.904,83.800,
C      5 85.468,87.620,88.906,91.220,92.906,95.940,98.906,101.070,
C      5 102.906,106.400,107.868,112.410,
C      5                 114.82,118.69,121.75,127.600,126.905,131.300,
C      6 132.905,137.330,138.906,
C      6                 140.120,140.908,144.240,145.000,150.400,
C      6 151.960,157.250,158.925,162.500,164.930,167.260,168.934,
C      6 173.040,174.970,
C      6 178.490,180.948,183.850,186.207,190.200,192.220,195.090,
C      6 196.967,207.590,204.370,207.200,208.980,210.000,210.001,
C      6 222.02,
C      7 223.000,226.025,227.028/     
      logical lit,lwr,lfix,lmscan,lspec,lshifte,lshifti,lshiftq,
     1llin,lcor,lonly,lsfix,lcre,lcom,ldz,swe,vert_h,j_one,we_eq_wg,
     1 se_eq_sg,lnoreor,vert_h_cor,vert_h_w,fix6,isplanar,J_ad,reg

      TOL=1.0d-10
      TOLL=TOL

      call readopt(lit,maxit,lwr,wmin,lfix,lmscan,NB,NP,TOLN,
     1lspec,del,e00,wsmin,wsmax,nps,lshifte,lshifti,lshiftq,
     1ncut,LEXCL,isu,llin,ifx,wlim,lcor,lonly,lsfix,lcre,kfac,
     1lcom,wfix,ldz,sslim,swe,vert_h,j_one,we_eq_wg,se_eq_sg,
     1 lnoreor,vert_h_cor,vert_h_w,delk,fix6,expGeom,expGeom_fac,
     1 Gauss_dusch,Gauss_dusch_file,isplanar,J_ad,deltol,reg,
     1 corenVH)
     
      call get_command_argument(1,value=argg,status=arg_stat)
      if(arg_stat==0)read(argg,*)kfac
      
      call get_command_argument(2,value=argg,status=arg_stat)
      if(arg_stat==0)then
         read(argg,*)hes_type
         select case(trim(adjustl(hes_type)))
            case('AH')
               vert_h=.false.
               we_eq_wg=.false.
               se_eq_sg=.false.
               J_one=.false.
               lnoreor=.false.
               write(6,*)'ADIABATIC HESSIAN'
            case('ASF')
               vert_h=.false.
               we_eq_wg=.false.
               se_eq_sg=.false.
               J_one=.true.
               lnoreor=.false.
               write(6,*)'ADIABATIC SHIFT with FREQUENCIES'
            case('AS')
               vert_h=.false.
               we_eq_wg=.true.
               se_eq_sg=.false.
               J_one=.true.
               lnoreor=.false.
               write(6,*)'ADIABATIC SHIFT'
            case('VH')
               vert_h=.true.
               we_eq_wg=.false.
               se_eq_sg=.false.
               J_one=.false.
               lnoreor=.true.
               write(6,*)'VERTICAL HESSIAN'
            case('VGF')
               vert_h=.true.
               we_eq_wg=.false.
               se_eq_sg=.false.
               J_one=.true.
               lnoreor=.true.
               write(6,*)'VERTICAL GRADIENT FREQUENCIES'
            case('VG')
               vert_h=.true.
               we_eq_wg=.true.
               se_eq_sg=.false.
               J_one=.true.
               lnoreor=.true.
               write(6,*)'VERTICAL GRADIENT'
            case default
               write(6,*)'READING P.E.S. SPECS FROM DUSCH.OPT'
         end select
      end if
c
c     Read Geometry of the ground and excited states:
      open(77,file='ground/ENERGY')
      read(77,*)e_gr
      close(77)
      open(77,file='excited/ENERGY')
      read(77,*)e_ex
      close(77)
      
      allocate(z(1),rg(1))
      if(lfix)then
       call rdxx('ground/fixed.x',0,nat,z,rg)
      else
       call rdxx('ground/FILE.X',0,nat,z,rg)
      endif
      deallocate(z,rg)
      N=3*nat
      allocate(rg(N),z(nat),re(N),sg(N,N),se(N,N),eg(N),ep(N),m(nat),
     1ref(N),rgf(N),sgf(N,N),sef(N,N))
      if(lfix)then
       call rdxx('ground/fixed.x',1,nat,z,rg)
C        if(vert_h)then
C          re=rg
C        else
         call rdxx('excited/fixed.x',1,nat,z,re)
C        end if
      else
       call rdxx('ground/FILE.X',1,nat,z,rg)
C        if(vert_h)then
C          re=rg
C        else
         call rdxx('excited/FILE.X',1,nat,z,re)
C        end if
      endif
c     assign masses, eventually do isotopic substitution:
      do 6 ia=1,nat
6     m(ia)=dble(amas(z(ia)))
      call readm(m)

c     read the ground (sg) and excited (se) S-matrices:
      write(6,*)'reading ground state S-matrix'
      call readsi_dusch(N,sg,eg,NQ,'ground/F.INP',rgf,ldz,iz)
      if(llin)then
       write(6,*)'excited state S-matrix same as ground'
       call readsi_dusch(N,se,ep,NE,'ground/F.INP',ref,ldz,iz)
      else
       write(6,*)'reading excited state S-matrix'
       call readsi_dusch(N,se,ep,NE,'excited/F.INP',ref,ldz,iz)
      endif
      
      nqi_g=count(eg<0d0)
      nqi_e=count(ep<0d0)
      
      if(.not.lnoreor)then
      call ortog(N,NE,se,m)
      call ortog(N,NQ,sg,m)
      write(6,*)'S-matrices corrected'
      end if
      ep_orig=ep
      if(we_eq_wg)then!
       !  ep_orig=ep
         ep=eg!
         se_eq_sg=.true.
      end if!
      if(se_eq_sg)then!
         se=sg!
      end if!
      
      
      
      allocate(del_modes_gr(NQ),del_modes_ex(NQ),MO(NQ))!
      del_modes_gr=0!
      del_modes_ex=0!
      sg_full=sg!
      se_full=se!
      MO=seq(NQ,.false.)
      nqn=NQ
      
      if(.not.lnoreor)then
      if(lfix)then
c      Find the transformation matrix from rgf to rg
       call xst(nat,rg,rgf,m,z,u,isplanar)
c      transform ground S-matrix to the rg system:
       call trse(N,sg,NQ,u)
       call trse(N,sg_full,NQ,u)!
c      Find the transformation matrix from ref to re
       call xst(nat,re,ref,m,z,u,isplanar)
c      transform excited S-matrix to the re system:
       call trse(N,se,NQ,u)
       call trse(N,se_full,NQ,u)!
      else
c      Find the transformation matrix, transform and shift re to rg:
       call xst(nat,rg,re,m,z,u,isplanar)
c      transform se to the ground system:
       call trse(N,se,NQ,u)
       call trse(N,se_full,NQ,u)!
      endif
      end if
      
      call writesi(N,nat,nq,sg,eg*au2cm,'F.INP.ground',z,rgf,.false.)
      call writesi(N,nat,nq,se,ep*au2cm,'F.INP.excited',z,ref,.false.)
      
      
      allocate(JT(NQ,NQ),K(NQ))
      call fixmodes(N,NQ,NE,eg,sg,ep,se,amas,z,ifx,wlim,lcor,lonly,
     1sgf,sef,NF,wfix,m,nat,sslim,swe,del_modes_gr,del_modes_ex,MO,
     1 vert_h,J_one,we_eq_wg,fix6,rg,re,delk,JT,K,iz,J_ad,ep_orig,
     1 deltol)
      if(nf>0 .and. vert_h)then
         allocate(help_m(N-6,N-6),help_v(N-6), 
     1    help_sg(N,NQ),help_se(N,NQ),
     1    eg_help(NQ),ep_help(NQ))
         nq=nq-nf
         help_M=JT
         help_V=K
         help_sg=sg
         help_se=se
         eg_help=eg
         ep_help=ep
         deallocate(JT,K,sg,se,eg,ep)
         allocate(JT(NQ,NQ),K(NQ),sg(N,NQ),se(N,NQ),eg(NQ),ep(nq))
         do i = 1,nq
            K(i)=help_V(i)
            eg(i)=eg_help(i)
            ep(i)=ep_help(i)
            do ii = 1,n
               sg(ii,i)=help_sg(ii,i)
               se(ii,i)=help_se(ii,i)
            end do
            do ii = 1,nq
               JT(i,ii)=help_m(i,ii)
            end do
         end do
         deallocate(help_m,help_v,help_se,help_sg,ep_help,eg_help)
         
         if(lsfix)call seproj(N,NF,NQ,se,sgf,m)
         if(lsfix)call seproj(N,NF,NQ,sg,sef,m)
         if(reg)then
            allocate(sg_m(N,NQ),m_invsq(N,N),ge_q(NQ),ffq(NQ,NQ))
            allocate(ffq_i(NQ,NQ))
            call makeSM(nat,N,nq,sg,m,sg_m,iz)
            allocate(w_new(NQ))
            JT=0d0
            K=0d0
            call mkj_vert(tol,amu2au,nat,N,NQ,m,sg_m,m_invsq,ge_q,
     1       ffq,ffq_i,JT,K,w_new,j_one,we_eq_wg,eg)
            deallocate(sg_m,m_invsq,ge_q,ffq,ffq_i)
            if(lcor)then
               call Reorder_J(N,NQ,MO,eg,ep,sg,se,JT,K,lonly,swe)
            end if
            JT=transpose(JT)
         end if
         !call mkk_vert(nat,N3,NQ,K,)
      else if(nf>0 .and. .not.vert_h)then !not implemented
         stop 2!deallocate(JT,K)
      end if
      if(lcor)then !I have no idea why this is here
         if(swe)then
            s_help=se_full
            do i = 1,nqn
               do kk = 1,N
                  se_full(kk,i)=s_help(kk,MO(i))
               end do
            end do
            deallocate(s_help)
         else
            s_help=sg_full
            do i = 1,nqn
               do kk = 1,N
                  sg_full(kk,i)=s_help(kk,MO(i))
               end do
            end do
            deallocate(s_help)
         end if
      end if
c     fix se for deleted modes:
      
      allocate(ffq(NQ,NQ),gr_q(NQ),ffc(N,N),gr(N),sg_m(N,NQ),se_m(N,NQ)
     1 ,m_invsq(N,N))
      ffq=0
      gr_q=0
      
      !Ground state gradient in normal mode coordinates of the ground states
      call mk_ffqgr(nat,N,NQ,sg,m,m_invsq,ffc,ffq,gr,gr_q,amu2au,iz,
     1 .false.,sg_m)
      !Gradient extrapolation when Adiabatic Hessian cause origin Independence, see RROA-Baiardi-2018
      if(.not.vert_h)then
         allocate(gr_exp(N),gr_exp_q(NQ))
         call ExtrapolateGradient(N,ffc,re,rg,bohr,gr_exp)
         call grad2Q(nat,nq,gr_exp,gr_exp_q,m,m_invsq,sg_m,amu2au)
         call svgrq(gr_exp_q,NQ,'FILE.Q.GR.EXTRAP')
         e_vert=VerticalEnergyAH(N,e_ex,e_gr,ffc,re,rg,bohr)
         open(77,file='ENERGY.VERTICAL')
         write(77,'(G19.12)')e_vert
         close(77)
         deallocate(gr_exp,gr_exp_q)
      else
         call readgrad_dusch(nat,gr,'excited/FILE.GR') !read excited state gradient
         call READFF(N,ffc,'excited/FILE.FC') !read excited state force field
         allocate(FFc_A(N,N))
         call PROJECT(FFc_A,FFc,N,'excited/FILE.X') 
         deallocate(FFc_A)
         K_orig=K
         K=K*kfac
         r_disp=matmul(M_invsq,matmul(sg_m,K))
         e_disp=dot_product(gr,r_disp)
         e_disp2=0.5d0*dot_product(r_disp,matmul(ffc,r_disp)) !For some reason, this term equals to 1/2 of negative gradient term in one case
         !ie. -dr^T.Hx.dr = dr^T.gx
         K=K_orig
         deallocate(r_disp,K_orig)
         
         !save the TDDFT "vertical approximation" energy
         open(77,file='ENERGY.VERTICAL')
         write(77,'(G19.12)')e_ex-e_gr
         close(77)
         
         if(corenVH==1)then
            e_ex=e_ex+e_disp
         elseif(corenVH==2)then
            e_ex=e_ex+e_disp+e_disp2
         else
            e_ex=e_ex
         end if
         
         open(77,file='excited/ENERGY')
         write(77,'(G19.12)')e_ex
         close(77)
         write(6,*)'Shifted excitation energy in VH appraoch.'
      end if
      call svgrq(gr_q,nQ,'FILE.Q.GR.ground')
      deallocate(ffq)
      deallocate(gr_q)
      
      !Excited state gradient in normal mode coordinates of the excited state
      allocate(ffq(NQ,NQ),gr_q(NQ))
      ffq=0
      gr_q=0
      call mk_ffqgr(nat,N,NQ,se,m,m_invsq,ffc,ffq,gr,gr_q,amu2au,iz,
     1 .true.,se_m)
      call svgrq(gr_q,nq,'FILE.Q.GR.excited')
      deallocate(ffq,gr_q)
      
      !Excited state gradient in normal mode coordinates of the ground state (which is basically used when using Vertical Hessian K)
      sg_m=0
      allocate(ffq(NQ,NQ),gr_q(NQ))
      ffq=0
      gr_q=0
      call mk_ffqgr(nat,N,NQ,sg,m,m_invsq,ffc,ffq,gr,gr_q,amu2au,iz,
     1 .true.,sg_m)
      call svgrq(gr_q,nq,'FILE.Q.GR.excited.ground')
      deallocate(ffq,gr_q)
      
      
c     experimental option-correct excited state so that it can be 
c     reached with available modes
      if(lcre)then
       call cre(N,NQ,m,rg,re,se,'e')
       call wrx(nat,re,z,'test2.x','test')
       call cre(N,NQ,m,rg,re,sg,'g')
       call wrx(nat,re,z,'test3.x','test')
      endif

c     if required, replace excited S-vectors by ground:
      if(lcom)call rpl(N,NQ,nat,m,sg,se)
c     Make Duschinsk's and derived matrices and vectors:
C       if(.not.vert_h)then
C          allocate(JT(NQ,NQ),K(NQ))
C       end if
      allocate(G(NQ,NQ),GP(NQ,NQ),GH(NQ,NQ),GHP(NQ,NQ),
     1J(NQ,NQ),T(NQ,NQ),TP(NQ,NQ),F(NQ,NQ),FI(NQ,NQ),E2(NQ,2*NQ),
     1A(NQ,NQ),B(NQ),TV(NQ),TU(NQ),C(NQ,NQ),D(NQ),E(NQ,NQ),K0(NQ))
      call zzz(NQ,JT,K,G,GP,GH,GHP,J,T,TP,F,FI,E2,A,B,TV,TU,C,D,E,K0,
     1 vert_h)

      IERR=0
      if(vert_h)then
      K=K*kfac
      fac=overlap_ABCDE00(1,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1 m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,
     1 Gauss_dusch,Gauss_dusch_file)
      else
      fac=overlap(1,NQ,N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1 m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac,
     1 vert_h,j_one,ifx,wfix,lcor,swe,MO,sg_full,se_full,WE_EQ_WG,
     1 vert_h_w,Gauss_dusch,Gauss_dusch_file)
      end if
      write(6,300)fac,fac*fac*100.0d0
300   format('<0|0*> = ',g12.4,' (',f6.2,'%)')
      open(77,file='OVERLAP')
      write(77,'(E22.15E2)')fac
      close(77)
      
      if(expGeom>0)then
         if(.not. vert_h)then
            write(6,*)
     1       'Extrapolated geometry skipped, not Vertical Hessian'
         else
            allocate(rexp(N))
            rexp=re
            call CalcExpGeom(nat,n,nq,sg,M,K,rexp,expGeom_fac)
            open(77,file='EXTRAP_GEOM')
            do i = 1,nat
               ix=3*(i-1)+1
               write(77,'(3F13.6))')rexp(ix),rexp(ix+1),rexp(ix+2)
            end do
            close(77)
            write(6,*)'EXTRAP_GEOM written'
            deallocate(rexp)
         end if
      end if

      
      write(6,*)' When K = 0:'
      call dofac(NQ,facv,J,JT,F,FI,G,GP,K0)
      write(6,300)facv,facv*facv*100.0d0
      open(77,file='OVERLAP_K0')
      write(77,'(E22.15E2)')facv
      close(77)


      if(lit)call iterate(maxit,NQ,lwr,wmin,
     1N,nat,G,GP,GH,GHP,eg,ep,JT,J,sg,se,
     1m,k,rg,re,bohr,T,TP,F,FI,E2,TOL,IERR,A,TV,TU,B,C,D,E,kfac,
     1 vert_h,j_one,ifx,wfix,lcor,swe,MO,sg_full,se_full,WE_EQ_WG,
     1 vert_h_w)

c     calculate normal mode gradient:
      call gradq(nat,N,m,ldz)
c     normal mode gradient in excited state but ground state geometry:
      call grade(nat,N,ldz)
c     calculate normal mode shift corresponding to geometry change:
      if(lshiftq)call shiftq(nat,N,NQ,se,ep,sg,eg)
c     different route:
      if(lshifti)call shifti(nat,N,NQ,se,ep)
c     different route, excited n's:
      if(lshifte)call shifte(nat,N,NQ,sg,se,eg,ep,m,lmscan,fac,
     1A,B,C,D,E,F,FI,G,GP,GH,GHP,J,JT,K,NB,NP,TOLN,
     2lspec,del,e00,wsmin,wsmax,nps,ncut,LEXCL,
     1isu)
      write(6,'(A,I5)')'    NQ = ',NQ
      write(6,'(A,I5)')'NQi_gr = ',nqi_g
      write(6,'(A,I5)')'NQi_ex = ',nqi_e
      write(6,'(A,I5)')'   DEL = ',NF
      write(6,'(A,F10.3,A)')'ZPE_E - ZPE_G = ',0.5*(sum(ep(1:nq))
     1 -sum(eg(1:nq)))*au2cm,' cm-1'
      
      
      end program dusch
