      PROGRAM GAR9
C     this program reads the archive of the Gaussian 94 output
C     version adapted for f77 on RS6000, from the garsgi.f that
C     ran in Calgary, on Silicon Graphics Irix; 2-12-1996; 
C     Petr Bour, UOCHB, Praha
      IMPLICIT none
      integer*4 i,ix,iy,iz,j,k,li,N,M,NAT,NB,ifr,ii1,ii2,n3,n4,iat,
     1N7,i1,j1,k1,i6,ia,iaa,id,ig,ii,in,iq,jaa,jq,MFR,MFRDET,NQ,kx,
     1MFRA
      real*8 amu,adum,CM,a,cqcm,dqcm,EPS,sum,APOL(9),GPOL(9),
     1ORT(9),FU,w,xx,yy,zz,xy,xz,yz,u(3,3),o(3),t27(27),
     2iargc,t81(81),AAPOL(27)
      CHARACTER*80 filename
      CHARACTER*1 BUFF(140),NUMBER(80),L,LS,LE
      CHARACTER*5 SSTR
      CHARACTER*2 LG,L2
      character*120 radek
      COMMON/LINE/I
      LOGICAL LVCD,LGEO,LARCH,auto,lex
      character*2 atsy(89)
      data atsy/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     3'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     4'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     4          'Ga','Ge','As','Se','Br','Kr',
     5'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     5          'In','Sn','Sb','Te','I ','Xe',
     6'Cs','Ba','La',
     6               'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
     6               'Er','Tm','Yb','Lu',
     6'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     6          'Tl','Pb','Bi','Po','At','Rn',
     7'Fr','Ra','Ac'/
      DIMENSION iaa(6),jaa(6)
      data iaa/1,2,3,1,1,2/
      data jaa/1,2,3,2,3,3/
      logical lsz,rdgeo,lnmr,lspin
      real*8,allocatable::ALPHA09(:,:),GTENS09(:,:),
     1ATENS09(:,:),F(:,:),ALPHA(:,:,:),GTENS(:,:,:),
     1ATENS(:,:,:,:),xg09(:,:),P(:,:,:),X(:,:),AAT(:,:,:),tm(:),
     1grad(:),t(:),s(:,:),e(:),ALPHAA(:,:,:),GGA(:,:,:),AAA(:,:,:),
     1C(:),D(:),ecd(:),rs(:),rz(:),APOLF(:,:),ww(:),AAPOLF(:,:),
     1ORTF(:,:),GPOLF(:,:),nmr(:),spin(:,:)
      CHARACTER*2,allocatable:: AT(:)
      logical,allocatable::lwr(:,:)
      integer*4,allocatable::iwr(:,:),qz(:)
c
      LG='%%'
      inquire(file='AUTO',exist=auto)
      if(auto)then
       FILENAME='FRE.OUT'
      else
       if(iargc().gt.0)then
        call getarg(1,filename)
       else
        WRITE(*,*)' Filename ? '
        READ(*,'(A)')filename
       endif
      endif

      NAT=0
c     allocate and define to avoid compiler messages:
      allocate(rs(1),rz(1),qz(1))
      qz(1)=0
      rs(1)=0.0d0
      rz(1)=0.0d0
      lsz=rdgeo(filename,nat,rs,rz,0,qz,atsy)
      deallocate(rs,rz,qz)
      if(lsz)then
       allocate(rs(3*nat),rz(3*nat),qz(nat))
       lsz=rdgeo(filename,nat,rs,rz,1,qz,atsy)
c      find the standard - Z-matrix transformation,by inertia:
       call xst(nat,rs,rz,u,o)
      else
       if(auto)then
        open(7,file='FILE.X')
        read(7,*) 
        read(7,*)NAT
        close(7)
       else
        WRITE(*,*)' How many atoms ?'
        READ(*,*)NAT
        WRITE(*,*)
     1  ' Geometry looked for after string ',LG,' (only for G94)'
       endif
       call vz(u,9)
       u(1,1)=1.0d0
       u(2,2)=1.0d0
       u(3,3)=1.0d0
      endif


c     determine and write the perturbation frequencies:
      MFR=mfrdet(FILENAME)
      call mfrdetw(FILENAME,MFR)
      N=3*NAT
c
c     if not frequencies,use one for allocation:
      MFRA=max(MFR,1)

      allocate (F(N,N),P(NAT,3,3),X(3,NAT),ALPHA(3,3,N),
     1AAT(NAT,3,3),GTENS(3,3,N),ATENS(3,3,3,N),AT(NAT),
     2xg09(6,3*N),grad(N),ALPHA09(MFRA,9*N),GTENS09(MFRA,9*N),
     1ATENS09(MFRA,27*N),lwr(100,MFRA),ecd(N),iwr(100,MFRA),
     1APOLF(MFRA,9),ww(MFRA),AAPOLF(MFRA,27),ORTF(MFRA,9),GPOLF(MFRA,9))

      do 24 i=1,100
      do 24 j=1,MFRA
      iwr(i,j)=0
24    lwr(i,j)=.false.
      amu=1822.88d0
      CM=219470.0d0
      call vz(GPOL,9)
      M=0
      n3=0
      n4=0
      N7=0
      call vz(AAT,3*3*NAT)
      call vz(P,3*3*NAT)
      call vz(ALPHA,3*3*N)
      call vz(GTENS,3*3*N)
      call vz(ATENS,3*3*3*N)
      NB=70
      LVCD=.FALSE.
      LGEO=.FALSE.
      lnmr=.false.
      lspin=.false.
      LS='\\'
      LE='='
      L2='\\\\'
      call vz(ALPHA09,MFRA*3*3*N)
      call vz(GTENS09,MFRA*3*3*N)
      call vz(ATENS09,MFRA*3*3*3*N)

      OPEN(7,FILE=FILENAME,STATUS='OLD',FORM='FORMATTED')
C
8000  read(7,8889,end=8999,err=8999)radek
8889  format(a120)

      if(radek(38:59).eq.'Forces (Hartrees/Bohr)')then
       if(.not.lwr(6,1))then
        lwr(6,1)=.true.
        read(7,*)
        read(7,*)
        do 105 i=1,NAT
105     read(7,*)(grad(ix+3*(i-1)),ix=1,2),(grad(ix+3*(i-1)),ix=1,3)
        write(6,*)' Gradient found'
        call svgr(grad,NAT)
       endif
      endif

      if(radek(1:37).eq.' Property number 1 -- Alpha(-w,w) fre')then
       read(radek(44:46),*)ifr
       read(radek(47:58),*)ww(ifr)
       iwr(1,ifr)=iwr(1,ifr)+1
       lwr(1,ifr)=.true.
       read(7,*)
       do 2 ix=1,3
2      read(7,*)ii1,(APOLF(ifr,ix+3*(iy-1)),iy=1,3)
      endif

      if(radek(1:52).eq.
     1 ' Property number 2 -- FD Optical Rotation Tensor fre')then
       read(radek(59:61),*)ifr
       iwr(2,ifr)=iwr(2,ifr)+1
       lwr(2,ifr)=.true.
       read(7,*)
       do 11 ix=1,3
11     read(7,*)ii1,(GPOL(ix+3*(iy-1)),iy=1,3)
       call rewr(0,APOL,9,APOLF,MFRA,ifr)
       call TG(APOL,GPOL,o)
       call smooth(GPOL,9)
       call rewr(1,GPOL,9,GPOLF,MFRA,ifr)
      endif

      if(radek(1:44).eq.
     1' Property number 4 -- D-Q polarizability fre')then
       read(radek(51:53),*)ifr
       iwr(4,ifr)=iwr(4,ifr)+1
       lwr(4,ifr)=.true.
       read(7,*)
       do 13 ii1=1,6
       ix=iaa(ii1)
       iy=jaa(ii1)
       read(7,*)ii2,(AAPOL(ix+3*(iy-1)+9*(iz-1)),iz=1,3)
       do 13 iz=1,3
13     AAPOL(iy+3*(ix-1)+9*(iz-1))=AAPOL(ix+3*(iy-1)+9*(iz-1))
       call rewr(0,APOL,9,APOLF,MFRA,ifr)
       call TA(APOL,AAPOL,o)
       call smooth(AAPOL,27)
       call rewr(1,AAPOL,27,AAPOLF,MFRA,ifr)

       if(lwr(1,ifr).and.lwr(2,ifr).and.(.not.lwr(3,ifr)))then
        lwr(3,ifr)=.true.
        iwr(3,ifr)=iwr(3,ifr)+1
c       optical rotation tensor
        call rewr(0,GPOL,9,GPOLF,MFRA,ifr)
        do 14 i=1,3
        do 14 j=1,3
        sum=GPOL(i+3*(j-1))+GPOL(j+3*(i-1))
        do 12 ig=1,3
        do 12 id=1,3
12      sum=sum-(EPS(i,ig,id)*AAPOL(id+3*(j-1)+9*(ig-1))
     1          +EPS(j,ig,id)*AAPOL(id+3*(i-1)+9*(ig-1)))/3.0d0
14      ORTF(ifr,i+3*(j-1))=sum/2.0d0
       endif
      endif

      if(radek(23:45).eq.'Alpha(-w,w) derivatives')then
       read(radek(57:59),*)ifr
       read(radek(62:71),*)ww(ifr)
       iwr(7,ifr)=iwr(7,ifr)+1
       lwr(7,ifr)=.true.
       ii1=0
8002   read(7,*)radek
       do 9 ix=1,3
9      read(7,8001)(xg09(ix,ii2),ii2=ii1+1,min(ii1+5,9*nat))
8001   format(8x,5e14.6)
       ii1=ii1+5
       if(ii1.lt.9*nat) goto 8002
       do 8 ix=1,3
       do 8 iy=1,3
       do 8 iat=1,n
8      ALPHA09(ifr,iy+3*(ix-1)+9*(iat-1))=xg09(ix,3*(iat-1)+iy)
      endif

      if(radek(23:60).eq.'FD Optical Rotation Tensor derivatives')then
       read(radek(72:74),*)ifr
       iwr(8,ifr)=iwr(8,ifr)+1
       lwr(8,ifr)=.true.
       ii1=0
8003   read(7,*)radek
       do 7 ix=1,3
7      read(7,8001)(xg09(ix,ii2),ii2=ii1+1,min(ii1+5,9*nat))
       ii1=ii1+5
       if(ii1.lt.9*nat) goto 8003

       do 6 ia=1,n/3

       do 61 kx=1,3
       iat=kx+3*(ia-1)
       do 61 ix=1,3
       do 61 iy=1,3
61     t27(ix+3*(iy-1)+9*(kx-1))=xg09(ix,iy+3*(iat-1))

       do 62 kx=1,3
       iat=kx+3*(ia-1)
       do 63 ix=1,3
       do 63 iy=1,3
       APOL(ix+3*(iy-1))=ALPHA09(ifr,ix+3*(iy-1)+9*(iat-1))
63     GPOL(ix+3*(iy-1))=t27(ix+3*(iy-1)+9*(kx-1))
       call TG(APOL,GPOL,o)
       do 62 ix=1,3
       do 62 iy=1,3
62     t27(ix+3*(iy-1)+9*(kx-1))=GPOL(ix+3*(iy-1))

       call smooth(t27,27)
       do 6 kx=1,3
       iat=kx+3*(ia-1)
       do 6 ix=1,3
       do 6 iy=1,3
6      GTENS09(ifr,iy+3*(ix-1)+9*(iat-1))=t27(ix+3*(iy-1)+9*(kx-1))

      endif

      if(radek(23:52).eq.'D-Q polarizability derivatives')then
       read(radek(64:66),*)ifr
       iwr(5,ifr)=iwr(5,ifr)+1
       lwr(5,ifr)=.true.
       ii1=0
8004   read(7,*)radek
       do 15 ix=1,6
15     read(7,8001)(xg09(ix,ii2),ii2=ii1+1,min(ii1+5,9*nat))
       ii1=ii1+5
       if(ii1.lt.9*nat) goto 8004

       do 16 ia=1,n/3

       do 161 kx=1,3
       iat=kx+3*(ia-1)
       do 161 iz=1,3
       ii2=iz+3*(iat-1)
       do 161 ii1=1,6
       ix=iaa(ii1)
       iy=jaa(ii1)
       t81(ix+3*(iy-1)+9*(iz-1)+27*(kx-1))=xg09(ii1,ii2)
 161   t81(iy+3*(ix-1)+9*(iz-1)+27*(kx-1))=xg09(ii1,ii2)

       do 162 kx=1,3
       iat=kx+3*(ia-1)
       do 163 ix=1,3
       do 163 iy=1,3
       APOL(ix+3*(iy-1))=ALPHA09(ifr,ix+3*(iy-1)+9*(iat-1))
       do 163 iz=1,3
 163   AAPOL(ix+3*(iy-1)+9*(iz-1))=t81(iy+3*(ix-1)+9*(iz-1)+27*(kx-1))
       call TA(APOL,AAPOL,o)
       do 162 ix=1,3
       do 162 iy=1,3
       do 162 iz=1,3
 162   t81(iy+3*(ix-1)+9*(iz-1)+27*(kx-1))=AAPOL(ix+3*(iy-1)+9*(iz-1))

       call smooth(t81,81)
       do 16 kx=1,3
       iat=kx+3*(ia-1)
       do 16 iz=1,3
       do 16 ix=1,3
       do 16 iy=1,3
16     ATENS09(ifr,iy+3*(ix-1)+9*(iz-1)+27*(iat-1))=
     1 t81(ix+3*(iy-1)+9*(iz-1)+27*(kx-1))

      endif

c     "Anharmonic" derivatives:
      if(radek(9:49).eq.'QUADRATIC FORCE CONSTANTS IN NORMAL MODES')then
       if(.not.lwr(9,1))then
        lwr(9,1)=.true.
        call s8
        M=0
3       read(7,*,end=88,err=88)adum,adum,adum
        M=M+1
        ecd(M)=adum/CM
        goto 3
88      write(6,*)M,' modes '
        allocate(C(M**3),D(M**3))
        call vz(C,M**3)
        call vz(D,M**3)
       endif
      endif

      if(radek(11:41).eq.'CUBIC FORCE CONSTANTS IN NORMAL'
     1 .and.lwr(9,1))then
       call s8
20     read(7,*,end=77,err=77)i,j,k,a
       c(i+M*(j-1)+M**2*(k-1))=a
       c(i+M*(k-1)+M**2*(j-1))=a
       c(j+M*(k-1)+M**2*(i-1))=a
       c(j+M*(i-1)+M**2*(k-1))=a
       c(k+M*(i-1)+M**2*(j-1))=a
       c(k+M*(j-1)+M**2*(i-1))=a
       n3=n3+1
       goto 20
77     if(n3.gt.0)then
        N7=N-M+1
        open(23,file='CQQ.SCR.TXT')
        do 64 i=N7,N
        i1=M-i+7
        do 64 j=i,N
        j1=M-j+7
        do 64 k=j,N
        k1=M-k+7
        cqcm=c(i1+M*(j1-1)+M**2*(k1-1))
        if(abs(cqcm).eq.0.0d0)goto 64
        write(23,2423)I,J,K,cqcm/CM,cqcm
2423    format(3I5,G16.8,' ',F9.1)
64      continue
        close(23)
        write(6,*)n3,' cubic constants in CQQ.SCR.TXT '
       endif
      endif

      if(radek(10:42).eq.'QUARTIC FORCE CONSTANTS IN NORMAL'
     1 .and.lwr(9,1))then
       call s8
21     read(7,*,end=66,err=66)i,j,k,li,a
       if(i.eq.j)then
        i1=i
        j1=k
        k1=li
       else
        if(i.eq.k)then
         i1=i
         j1=j
         k1=li
        else
         if(i.eq.li)then
          i1=i
          j1=j
          k1=k
         else
          if(j.eq.k)then
           i1=j
           j1=i
           k1=li
          else
           if(j.eq.li)then
            i1=j
            j1=i
            k1=k
           else
            if(k.eq.li)then
             i1=k
             j1=i
             k1=j
            else
             write(6,*)i,j,k,li,' general quartic not implemented'
             stop
            endif
           endif
          endif
         endif
        endif
       endif
       d(i1+M*(j1-1)+M**2*(k1-1))=a
       d(i1+M*(k1-1)+M**2*(j1-1))=a
       n4=n4+1
       goto 21
66     if(n4.gt.0)then
        N7=N-M+1
        open(23,file='DQQ.SCR.TXT')
        do 621 i=n7,n
        i1=M-i+7
        do 621 j=n7,n
        j1=M-j+7
        do 621 k=j,n
        k1=M-k+7
        dqcm=d(i1+M*(j1-1)+M**2*(k1-1))
        if(abs(dqcm).eq.0.0d0)goto 621
        write(23,5333)I,J,K,dqcm/CM,dqcm
5333    format(3I6,G16.8,' ',f9.1)
621     continue
        close(23)
        write(6,*)n4,' quartic constants in DQQ.SCR.TXT'
       endif
      endif
c
c     Mixed normal mode - Cartesian alpha derivatives 
      if(radek(23:56).eq.'DAlpha(-w,w) numerical derivatives')then
       read(radek(68:70),*)ifr
       read(radek(71:82),*)w
       if(.not.lwr(10,ifr))then
        write(6,*)' Mixed alpha second derivatives found ',ifr
        lwr(10,ifr)=.true.
        write(6,609)w
609     format(' frequency: ',f12.6,' au')
        N=3*nat
        M=3*nat-6
        allocate(t(9*N*M),s(N,N),e(N))
        call rgt(N,9,M,t)
        inquire(file='F.INP',exist=lex)
        if(lex)then
         call readsi(N,S,E,NQ,'F.INP')
         if(NQ.eq.M)then
          allocate(tm(9*M*M))
          call trt(t,tm,N,M,s,9)
          if(ifr.eq.1)then
           allocate(ALPHAA(N**2,9,MFRA),GGA(N**2,9,MFRA),
     1     AAA(N**2,27,MFRA))
c          use redundant space for the first index to be compatible with S4
           call vz(ALPHAA,9*MFRA*N**2)
           call vz(GGA,9*MFRA*N**2)
           call vz(AAA,27*MFRA*N**2)
          endif
          do 17 iq=1,M
          do 17 jq=1,M
          a=1.0d0/amu/dsqrt(dabs(e(iq)*e(jq)))
          do 17 i=1,9
17        ALPHAA(iq+N*(jq-1),i,ifr)=tm(iq+M*(jq-1)+M*M*(i-1))*a
          N7=N-M+1
          call wrttqaq(M,N,ALPHAA,GGA,AAA,MFRA,ifr,w*CM,N7)
          deallocate(tm)
         else
          write(6,*)'F.INP found, but number of modes does not match'
         endif
        else
         write(6,*)' F.INP not found, derivatives forgotten'
        endif
        deallocate(t,s,e)
       endif
      endif
c
c     Mixed normal mode - Cartesian G derivatives 
      if(radek(23:71).eq.
     1'DFD Optical Rotation Tensor numerical derivatives')then
       read(radek(83:85),*)ifr
       read(radek(86:97),*)w
       if(.not.lwr(11,ifr))then
        write(6,*)' Mixed G second derivatives found ',ifr
        lwr(11,ifr)=.true.
        write(6,609)w
        M=3*nat-6
        allocate(t(9*N*M),s(N,N),e(N))
        call rgt(N,9,M,t)
        inquire(file='F.INP',exist=lex)
        if(lex)then
         call readsi(N,S,E,NQ,'F.INP')
         if(NQ.eq.M)then
          allocate(tm(9*M*M))
          call trt(t,tm,N,M,s,9)
          do 18 iq=1,M
          do 18 jq=1,M
          a=1.0d0/amu/dsqrt(dabs(e(iq)*e(jq)))
          do 18 ix=1,3
          do 18 iy=1,3
          i=ix+3*(iy-1)
          j=iy+3*(ix-1)
18        GGA(iq+N*(jq-1),j,ifr)=tm(iq+M*(jq-1)+M*M*(i-1))*a
          N7=N-M+1
          call wrttqaq(M,N,ALPHAA,GGA,AAA,MFRA,ifr,w*CM,N7)
          deallocate(tm)
         else
          write(6,*)'F.INP found, but number of modes does not match'
         endif
        else
         write(6,*)' F.INP not found, derivatives forgotten'
        endif
        deallocate(t,s,e)
       endif
      endif
c
c     Mixed normal mode - Cartesian A derivatives 
      if(radek(23:62).eq.'DD-Q polarizability numerical derivative')then
       read(radek(75:77),*)ifr
       read(radek(78:89),*)w
       if(.not.lwr(12,ifr))then
        lwr(12,ifr)=.true.
        write(6,*)' Mixed A second derivatives found ',ifr
        M=3*nat-6
        allocate(t(18*N*M),s(N,N),e(N))
        call rgt(N,18,M,t)
        inquire(file='F.INP',exist=lex)
        if(lex)then
         call readsi(N,S,E,NQ,'F.INP')
         if(NQ.eq.M)then
          allocate(tm(18*M*M))
          call trt(t,tm,N,M,s,18)
          do 19 iq=1,M
          do 19 jq=1,M
          a=1.5d0/amu/dsqrt(dabs(e(iq)*e(jq)))
          ii=0
          do 19 i=1,3
          xx=tm(iq+M*(jq-1)+M*M*(ii+1-1))*a
          yy=tm(iq+M*(jq-1)+M*M*(ii+2-1))*a
          zz=tm(iq+M*(jq-1)+M*M*(ii+3-1))*a
          xy=tm(iq+M*(jq-1)+M*M*(ii+4-1))*a
          xz=tm(iq+M*(jq-1)+M*M*(ii+5-1))*a
          yz=tm(iq+M*(jq-1)+M*M*(ii+6-1))*a
          AAA(iq+N*(jq-1),1+3*(1-1)+9*(i-1),ifr)=xx
          AAA(iq+N*(jq-1),1+3*(2-1)+9*(i-1),ifr)=xy
          AAA(iq+N*(jq-1),1+3*(3-1)+9*(i-1),ifr)=xz
          AAA(iq+N*(jq-1),2+3*(1-1)+9*(i-1),ifr)=xy
          AAA(iq+N*(jq-1),2+3*(2-1)+9*(i-1),ifr)=yy
          AAA(iq+N*(jq-1),2+3*(3-1)+9*(i-1),ifr)=yz
          AAA(iq+N*(jq-1),3+3*(1-1)+9*(i-1),ifr)=xz
          AAA(iq+N*(jq-1),3+3*(2-1)+9*(i-1),ifr)=yz
          AAA(iq+N*(jq-1),3+3*(3-1)+9*(i-1),ifr)=zz
19        ii=ii+6
          N7=N-M+1
          call wrttqaq(M,N,ALPHAA,GGA,AAA,MFRA,ifr,w*CM,N7)
          deallocate(tm)
         else
          write(6,*)'F.INP found, but number of modes does not match'
         endif
        else
         write(6,*)' F.INP not found, derivatives forgotten'
        endif
        deallocate(t,s,e)
       endif
      endif

c     NMR parameters
      if((radek(11:21).eq.'Isotropic ='.or.
     1 radek(15:25).eq.'Isotropic =').and..not.lnmr.and.nat.gt.0)then
       lnmr=.true.
       write(6,*)'NMR shifts found'
       allocate(nmr(nat))
       backspace 7
       do 22 i=1,nat
       read(7,8889,end=8999,err=8999)radek
       if(radek(11:21).eq.'Isotropic =')then
        read(radek(22:32),*)nmr(i)
       else
        read(radek(26:36),*)nmr(i)
       endif
       do 22 j=1,4
22     read(7,*)
       call wrnmr(nat,nmr,qz)
      endif

c     NMR Spin-Spin (J) coupling
      if(radek(2:41).eq.'Total nuclear spin-spin coupling J (Hz):'.
     1 and..not.lspin)then
       lspin=.true.
       write(6,*)'NMR spin-spin couplings found'
       allocate(spin(nat,nat))
       call rspin(spin,nat)
       call wspin(spin,nat)
      endif

      goto 8000
8999  rewind(7)


c     Reading the archive:
c 
      LARCH=.FALSE.
5444  READ(7,7000,END=8888)(BUFF(I),I=1,NB)
      READ(7,7000,END=8888)(BUFF(I),I=NB+1,2*NB)
7000  FORMAT(1X,70A1)
      BACKSPACE 7
      DO 51 I=1,2*NB-5
      DO 52 J=1,5
52    SSTR(J:J)=BUFF(I+J-1) 
      IF(SSTR(1:2).EQ.L2)LARCH=.TRUE.
c     IF(SSTR.EQ.'leDer')THEN
      if(I.gt.3)then
       IF(SSTR.EQ.'leDer'.and.BUFF(I-3).eq.'i')THEN
        IF(.NOT.LARCH)GOTO 5444
        WRITE(*,*)' Atomic polar tensor found'
        GOTO 6555
       ENDIF
      endif
      IF(SSTR(3:5).EQ.'AAT')THEN
       IF(.NOT.LARCH)GOTO 5444
       WRITE(*,*)' Atomic axial tensor found'
       LVCD=.TRUE.
       GOTO 9555
      ENDIF
      IF(SSTR.EQ.'arDer')THEN
       IF(.NOT.LARCH)GOTO 5444
       WRITE(*,*)' Polarization Derivatives found'
       GOTO 8555
      ENDIF
      IF(SSTR.EQ.'olDer')THEN
       IF(.NOT.LARCH)GOTO 5444
       WRITE(*,*)' A-tensor Derivatives found'
       GOTO 11555
      ENDIF
      IF(SSTR.EQ.'otDer'.or.SSTR.EQ.'nsDer')THEN
       IF(.NOT.LARCH)GOTO 5444
       WRITE(*,*)' G-tensor Derivatives found'
       GOTO 10555
      ENDIF
C     NImag for g94, NIMAG for g92
      IF(SSTR.EQ.'NImag'.OR.SSTR.EQ.'NIMAG')THEN
       IF(.NOT.LARCH)GOTO 5444
       WRITE(*,*)' Second derivatives found'
       GOTO 5555
      ENDIF
      IF((.NOT.LGEO).AND.SSTR(1:2).EQ.LG)THEN
       IF(.NOT.LARCH)GOTO 5444
       WRITE(*,*)' Geometry found'
       LGEO=.TRUE.
       GOTO 7555
      ENDIF
51    CONTINUE
      GOTO 5444
C
7555  I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.NB)THEN
       DO 74 K=1,NB
74     BUFF(K)=BUFF(K+NB)
       I=I-NB
       READ(7,*)
      ENDIF
      DO 73 IA=1,NAT
      AT(IA)(1:1)=' '
      AT(IA)(2:2)=' '
C
C     Get element:
      I=I+1
      CALL GETLET(L,BUFF,NB)
      AT(IA)(1:1)=L
      I=I+1
      CALL GETLET(L,BUFF,NB)
      if(L.ne.',')then
       AT(IA)(2:2)=L
       I=I+1
       CALL GETLET(L,BUFF,NB)
      endif
C
      DO 773 IX=1,3
C     READ X(IX,IA):
      DO 75 J=1,20
75    NUMBER(J)=' '
      IN=0
7558  IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.LS)GOTO 7558
      NUMBER(IN)=' '
      X(IX,IA)=FU(NUMBER)
773   CONTINUE
73    CONTINUE
c
      OPEN(15,FILE='FILE.X')
      WRITE(15,*)' FILE.X from Gaussian output'
      WRITE(15,*)NAT
      DO 777 IA=1,NAT
      K=0
      do 1 ii=1,89
1     IF(AT(IA).EQ.atsy(ii))K=ii
      IF(K.EQ.0)call report(' Unknown atom '//AT(IA))
777   WRITE(15,1501)K,(X(IX,IA),IX=1,3)
1501  FORMAT(I4,3F15.8,' 0 0 0 0 0 0 0 0.0')
      CLOSE(15)
      WRITE(*,*)' X written into FILE.X'
c
      BACKSPACE 7
      GOTO 5444
C
C
C     G' derivatives
10555 I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.NB)THEN
       DO 786 K=1,NB
786    BUFF(K)=BUFF(K+NB)
       I=I-NB
       READ(7,*)
      ENDIF
85561 I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LE)GOTO 85561
      DO 801 IA=1,N
      DO 801 IX=1,3
      DO 801 IY=1,3
C
C     READ G'(IX,IY,IA)
      DO 791 J=1,20
791   NUMBER(J)=' '
      IN=0
85581 IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.'\\'.AND.L.NE.'\')GOTO 85581
      NUMBER(IN)=' '
801   GTENS(IX,IY,IA)=FU(NUMBER)
      CALL WRITETTT(NAT,ALPHA,GTENS,ATENS,N)
      WRITE(*,*)' GTENS written into FILE.TTT'
      BACKSPACE 7
      GOTO 5444
C
C
C     A (dipole-quadrupole) derivatives
11555 I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.NB)THEN
       DO 787 K=1,NB
787    BUFF(K)=BUFF(K+NB)
       I=I-NB
       READ(7,*)
      ENDIF
85562 I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LE)GOTO 85562
      DO 802 IA=1,N
      DO 802 IZ=1,3
      DO 802 I6=1,6
      IX=iaa(I6)
      IY=jaa(I6)
C
C     READ A(IX,IY,IZ,IA)
      DO 792 J=1,20
792   NUMBER(J)=' '
      IN=0
85582 IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.'\\'.AND.L.NE.'\')GOTO 85582
      NUMBER(IN)=' '
      ATENS(IY,IX,IZ,IA)=FU(NUMBER)
802   ATENS(IX,IY,IZ,IA)=FU(NUMBER)
      CALL WRITETTT(NAT,ALPHA,GTENS,ATENS,N)
      WRITE(*,*)' ATENS written into FILE.TTT'
      BACKSPACE 7
      GOTO 5444
C
C
C     Polarization derivatives
8555  I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.NB)THEN
       DO 78 K=1,NB
78     BUFF(K)=BUFF(K+NB)
       I=I-NB
       READ(7,*)
      ENDIF
8556  I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LE)GOTO 8556
      DO 80 IA=1,N
      DO 80 IX=1,3
      DO 80 IY=1,IX
C
C     READ ALPHA(IX,IY,IA)
      DO 79 J=1,20
79    NUMBER(J)=' '
      IN=0
8558  IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.'\\'.AND.L.NE.'\')GOTO 8558
      NUMBER(IN)=' '
      ALPHA(IY,IX,IA)=FU(NUMBER)
80    ALPHA(IX,IY,IA)=FU(NUMBER)
      CALL WRITETTT(NAT,ALPHA,GTENS,ATENS,N)
      WRITE(*,*)' ALPHA written into FILE.TTT'
      BACKSPACE 7
      GOTO 5444
C
C     Dipole derivatives     
6555  I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.NB)THEN
       DO 54 K=1,NB
54     BUFF(K)=BUFF(K+NB)
       I=I-NB
       READ(7,*)
      ENDIF
6556  I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LE)GOTO 6556
      DO 53 IA=1,NAT
      DO 53 IX=1,3
      DO 53 IY=1,3
C
C     READ P(IA,IX,IY):
      DO 55 J=1,20
55    NUMBER(J)=' '
      IN=0
6558  IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.'\\'.AND.L.NE.'\')GOTO 6558
      NUMBER(IN)=' '
53    P(IA,IX,IY)=FU(NUMBER)
      CALL WRITETEN(NAT,P,AAT,LVCD,NAT)
      BACKSPACE 7
      GOTO 5444
C
C
C     Magnetic dipole derivatives     
9555  I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.NB)THEN
       DO 56 K=1,NB
56     BUFF(K)=BUFF(K+NB)
       I=I-NB
       READ(7,*)
      ENDIF
9556  I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LE)GOTO 9556
      DO 58 IA=1,NAT
      DO 58 IX=1,3
      DO 58 IY=1,3
C
C     READ AAT(IA,IX,IY):
      DO 59 J=1,20
59    NUMBER(J)=' '
      IN=0
9558  IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.'\\'.AND.L.NE.'\')GOTO 9558
      NUMBER(IN)=' '
58    AAT(IA,IX,IY)=FU(NUMBER)
      CALL WRITETEN(NAT,P,AAT,LVCD,NAT)
      WRITE(*,*)' AAT written into FILE.TEN'
      BACKSPACE 7
      GOTO 5444
C
C     Force Field:
5555  I=I+1
C     I .. index of the letter to be gotten
      IF(I.GT.70)THEN
       DO 4 K=1,70
4      BUFF(K)=BUFF(K+70)
       I=I-70
       READ(7,*)
      ENDIF
C     NB .. line length
5556  I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LS)GOTO 5556
5557  I=I+1
      CALL GETLET(L,BUFF,NB)
      IF(L.NE.LS)GOTO 5557
c     write(6,*)'numbers',I

      DO 23 IX=1,N
      DO 23 IY=1,IX
C
C     READ F(IX,IY):
      DO 5 J=1,20
5     NUMBER(J)=' '
      IN=0
5558  IN=IN+1
      I=I+1
      CALL GETLET(L,BUFF,NB)
      NUMBER(IN)=L
      IF(L.NE.','.AND.L.NE.'\\'.and.L.NE.'\')GOTO 5558
      NUMBER(IN)=' '
      F(IX,IY)=FU(NUMBER)
23    F(IY,IX)=FU(NUMBER)

      OPEN(20,FILE='FILE.FC')
      CALL WRITEFF(N,N,F)
      CLOSE(20)
      WRITE(*,*)' FF written into FILE.FC'
      BACKSPACE 7
      GOTO 5444
C
C
8888  CLOSE(7)

      do 65 ifr=1,MFR
      write(6,6082)ifr
6082  format(/,' Frequency ',i3,':')
      if(lwr(1,ifr))write(6,6081)' Alpha            ',iwr(1,ifr)
      if(lwr(2,ifr))write(6,6081)' GP               ',iwr(2,ifr)
      if(lwr(3,ifr))write(6,6081)' ORT              ',iwr(3,ifr)
      if(lwr(4,ifr))write(6,6081)' A                ',iwr(4,ifr)
      if(lwr(5,ifr))write(6,6081)' A     derivatives',iwr(5,ifr)
      if(lwr(7,ifr))write(6,6081)' Alpha derivatives',iwr(7,ifr)
      if(lwr(8,ifr))write(6,6081)' GP    derivatives',iwr(8,ifr)
6081  format(A18,' read ',i4,' - times')

      call rewr(0,APOL , 9,APOLF ,MFRA,ifr)
      call rewr(0,GPOL , 9,GPOLF ,MFRA,ifr)
      call rewr(0,ORT  , 9,ORTF  ,MFRA,ifr)
      call rewr(0,AAPOL,27,AAPOLF,MFRA,ifr)
      if(lwr(1,ifr))call WRITEPOL(APOL ,ifr,'POL.TTT',ww(ifr))
      if(lwr(2,ifr))call WRITEPOL(GPOL ,ifr,'GP.TTT' ,ww(ifr))
      if(lwr(3,ifr))call WRITEPOL(ORT,  ifr,'ORT.TTT',ww(ifr))
      if(lwr(4,ifr))call WRITEPOL(AAPOL,ifr,'A.TTT'  ,ww(ifr))
      if(lwr(5,ifr).or.lwr(7,ifr).or.lwr(8,ifr))
     1CALL WRITETTT09(NAT,ALPHA09,GTENS09,ATENS09,N,MFRA,ifr,ww(ifr))
      
65    continue
      
      END
 
      SUBROUTINE WRITEFF(MX3,N,FCAR)
      IMPLICIT INTEGER*4 (I-N) 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FCAR(MX3,N)
C     CONST=4.359828/0.5291772**2
      CONST=1.0d0
      DO 6 I=1,N
      DO 6 J=1,N
6     FCAR(I,J)=FCAR(I,J)/CONST
      N1=1
1     N3=N1+4
      IF(N3.GT.N)N3=N
      DO 130 LN=N1,N
130   WRITE(20,17)LN,(FCAR(LN,J),J=N1,MIN(LN,N3))
      N1=N1+5
      IF(N3.LT.N)GOTO 1
17    FORMAT(I4,5D14.6)
      RETURN
      END

      SUBROUTINE GETLET(L,BUFF,NB)
      IMPLICIT INTEGER*4 (I-N) 
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 BUFF(140)
      CHARACTER*1 L
      COMMON/LINE/I
      IF(I.GT.NB)THEN
       IF(I.GT.NB+1)THEN
        WRITE(*,*)' Cannot read letter ',I
c       STOP
       ENDIF
       I=1
       READ(7,7000)(BUFF(J),J=1,70)
7000   FORMAT(1X,70A1)
      ENDIF
      L=BUFF(I)
      RETURN
      END
        
      SUBROUTINE WRITETEN(N0,P,A,LVCD,NAT)
      IMPLICIT INTEGER*4 (I-N) 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(N0,3,3),A(N0,3,3)
      LOGICAL LVCD
c     BOHR=0.52917705993d0
      Z0=0.0d0
      OPEN(15,FILE='FILE.TEN')
      WRITE(15,1500) NAT,NAT-6,0
1500  FORMAT(3I5)
      DO 10 L=1,NAT
      DO 10 J=1,3
10    WRITE(15,1501) (P(L,J,I),I=1,3),L
1501  FORMAT(3F14.8,I5)
      IF(.NOT.LVCD)THEN
       DO 220 L=1,NAT                                                 
       DO 220 J=1,3                                                    
c      write(6,*)L,J,nat
220    WRITE(15,1501) (Z0,I=1,3),L
      ELSE
       DO 221 L=1,NAT                                                 
       DO 221 J=1,3                                                    
c      write(6,*)L,J
221    WRITE(15,1501) (A(L,J,I),I=1,3),L
      ENDIF
      DO 230 L=1,NAT                                                 
      DO 230 J=1,3                                                    
230   WRITE(15,1501) (Z0,I=1,3),L
      DO 100 L=1,NAT
      DO 100 J=1,3
100   WRITE(15,1501) (P(L,J,I),I=1,3),L
      WRITE(*,*)' Dipole derivatives written into FILE.TEN'
      CLOSE(15)
      RETURN
      END
C                                                                      
      FUNCTION EPS(I,J,K)
      IMPLICIT INTEGER*4 (I-N) 
      REAL*8 EPS
      INTEGER*4 I,J,K
      EPS=0.0d0
      IF (I.EQ.1.AND.J.EQ.2.AND.K.EQ.3)EPS= 1.0d0
      IF (I.EQ.1.AND.J.EQ.3.AND.K.EQ.2)EPS=-1.0d0
      IF (I.EQ.2.AND.J.EQ.3.AND.K.EQ.1)EPS= 1.0d0
      IF (I.EQ.2.AND.J.EQ.1.AND.K.EQ.3)EPS=-1.0d0
      IF (I.EQ.3.AND.J.EQ.1.AND.K.EQ.2)EPS= 1.0d0
      IF (I.EQ.3.AND.J.EQ.2.AND.K.EQ.1)EPS=-1.0d0
      RETURN                                                          
      END                                                             
C 
c     SUBROUTINE VCDD0(VCD,VEL,DIP,C,NAT)
C     STOLEN FROM CADPAC, INDICES ARE DIFFERENT !!, VEL WITH NUCLEI
c     IMPLICIT REAL*8(A-H,O-Z)
c     IMPLICIT INTEGER*4 (I-N)
c     INTEGER*4 ALPHA,BETA,GAMMA,DELTA,LAMBDA                           
c     PARAMETER (MX=100)
c     DIMENSION VCD(MX,3,3),VEL(MX,3,3),DIP(MX,3,3),                  
c    1C(3,NAT)
c     DO 90 BETA =1,3                                                 
c     DO 90 GAMMA=1,3                                                 
c     DO 90 DELTA=1,3                                                 
c     SKEW=0.25D0*EPS(BETA,GAMMA,DELTA)                               
c     IF (DABS(SKEW).GT.1.0D-10)THEN
c     DO 80 ALPHA=1,3                                                 
c     DO 80 LAMBDA=1,NAT                                              
c80    VCD(LAMBDA,ALPHA,BETA)=VCD(LAMBDA,ALPHA,BETA)                   
c    1 -SKEW*C(GAMMA,LAMBDA)                                          
c    2 *(VEL(LAMBDA,ALPHA,DELTA)-DIP(LAMBDA,ALPHA,DELTA))             
c     ENDIF                                                           
c90    CONTINUE
c     RETURN
c     END
C     =======================================
      SUBROUTINE WRITETTT(NAT,ALPHA,GTENS,ATENS,MX3)
      IMPLICIT none
      INTEGER*4 NAT,MX3,I,L,IX,IIND,J,K,N,NQ,IA,IM,IQ,N7
      REAL*8 ALPHA(3,3,MX3),GTENS(3,3,MX3),ATENS(3,3,3,MX3),
     1A0(3,3),G0(3,3),AT0(3,3,3),a,amu
      logical lex
      real*8,allocatable::e(:),s(:,:),ALPHAQ(:,:,:),GTENSQ(:,:,:),
     1ATENSQ(:,:,:,:)
      integer*4,allocatable::nml(:)
      OPEN(2,FILE='FILE.TTT')
      WRITE(2,2000)NAT
2000  FORMAT(' ROA tensors, cartesian derivatives',/,I4,' atoms',/,
     1' The electric-dipolar electric-dipolar polarizability:',/,
     2' Atom/x    jx           jy           jz')
      DO 1 I=1,3
      WRITE(2,2002)I
2002  FORMAT(' Alpha(',I1,',J):')
      DO 1 L=1,NAT
      DO 1 IX=1,3
      IIND=3*(L-1)+IX
1     WRITE(2,2001)L,IX,(ALPHA(I,J,IIND),J=1,3)
2001  FORMAT(I5,1H ,I1,3g15.7)
      WRITE(2,2004)
2004  FORMAT(' The electric dipole magnetic dipole polarizability:',/,
     1       ' Atom/x    jx(Bx)       jy(By)       jz(Bz)')
      DO 2 I=1,3
      WRITE(2,2003)I
2003  FORMAT(' G(',I1,',J):')
      DO 2 L=1,NAT
      DO 2 IX=1,3
      IIND=3*(L-1)+IX
2     WRITE(2,2001)L,IX,(GTENS(I,J,IIND),J=1,3)
      WRITE(2,2005)
2005  FORMAT(' The electric dipole electric quadrupole polarizability:',
     2/,     ' Atom/x    kx           ky           kz')
      DO 3 I=1,3
      DO 3 J=1,3
      WRITE(2,2006)I,J
2006  FORMAT(' A(',I1,',',I1,',K):')
      DO 3 L=1,NAT
      DO 3 IX=1,3
      IIND=3*(L-1)+IX
3     WRITE(2,2007)L,IX,(ATENS(K,J,I,IIND)*3.d0/2.d0,K=1,3),L,IX,I,J
2007  FORMAT(I5,1H ,I1,3F15.7,' ',4i3)
      write(2,*)
      write(2,*)'dummy alpha v:'
      DO 4 I=1,3
      WRITE(2,2002)I
      DO 4 L=1,NAT
      DO 4 IX=1,3
      IIND=3*(L-1)+IX
4     WRITE(2,2001)L,IX,(ALPHA(I,J,IIND),J=1,3)
      CLOSE(2)

c     normal mode derivatives in S4 format:
      N=3*NAT
      inquire(file='F.INP',exist=lex)
      if(lex)then
       amu=1822.88d0
       allocate(s(N,N),e(N),nml(N),
     1 ALPHAQ(3,3,N),GTENSQ(3,3,N),ATENSQ(3,3,3,N))
       call readsi(N,S,E,NQ,'F.INP')
       N7=N-NQ+1
       do 101 I=1,3
       do 101 J=1,3
       AT0(I,J,1)=0.0d0
       AT0(I,J,2)=0.0d0
       AT0(I,J,3)=0.0d0
       A0( I,J  )=0.0d0
       G0( I,J  )=0.0d0
       do 101 IQ=1,NQ
       IM=N7+IQ-1
       nml(NQ-IQ+1)=IM
       ATENSQ(I,J,1,IM)=0.0d0
       ATENSQ(I,J,2,IM)=0.0d0
       ATENSQ(I,J,3,IM)=0.0d0
       ALPHAQ(I,J,  IM)=0.0d0
       GTENSQ(I,J,  IM)=0.0d0
       do 101 IA=1,N
       ATENSQ(I,J,1,IM)=ATENSQ(I,J,1,IM)+ATENS(I,J,1,IA)*s(IA,IQ)
       ATENSQ(I,J,2,IM)=ATENSQ(I,J,2,IM)+ATENS(I,J,2,IA)*s(IA,IQ)
       ATENSQ(I,J,3,IM)=ATENSQ(I,J,3,IM)+ATENS(I,J,3,IA)*s(IA,IQ)
       ALPHAQ(I,J,  IM)=ALPHAQ(I,J,  IM)+ALPHA(I,J,  IA)*s(IA,IQ)
101    GTENSQ(I,J,  IM)=GTENSQ(I,J,  IM)+GTENS(I,J,  IA)*s(IA,IQ)
       do 17 IQ=1,NQ
       IM=N7+IQ-1
       a=1.0d0/dsqrt(dabs(amu*e(IQ)))
       do 17 I=1,3
       do 17 J=1,3
       ATENSQ(I,J,1,IM)=ATENSQ(I,J,1,IM)*a*1.5d0
       ATENSQ(I,J,2,IM)=ATENSQ(I,J,2,IM)*a*1.5d0
       ATENSQ(I,J,3,IM)=ATENSQ(I,J,3,IM)*a*1.5d0
       ALPHAQ(I,J,IM)=ALPHAQ(I,J,IM)*a
17     GTENSQ(I,J,IM)=GTENSQ(I,J,IM)*a
       call wrttq(A0,NQ,nml,ALPHAQ,ALPHAQ,G0,AT0,
     1 GTENSQ,ATENSQ,N)
      endif

      RETURN
      END
c jk-start
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wrttq(polar0,ndif,nml,ALPHA0Q,alphan,gten0,aten0,
     1G0Q,A0Q,NT3)
      implicit none
      integer*4 NT3
      real*8 polar0(3,3),gten0(3,3),aten0(3,3,3),A0Q(3,3,3,NT3),
     1ALPHA0Q(3,3,NT3),G0Q(3,3,NT3),alphan(3,3,NT3)
      integer*4 i,j,k,nml(*),l,ndif
      open(9,file='TTTQ.TXT.SCR')
      write(9,9004)((polar0(i,j),i=1,3),j=1,3)
9004  format(3(3f15.7,' au',/),
     1' Normal mode polarizability derivatives num anal')
      write(9,*)ndif
      DO 19 I=1,ndif
      do 19 j=1,3
      do 19 k=1,3
19    write(9,9005)I,nml(I),j,k,alphan(j,k,nml(I)),ALPHA0Q(j,k,nml(I))
9005  format(4i5,2f12.6)
      write(9,9006)((gten0(i,j),i=1,3),j=1,3)
9006  format(3(3f15.7,' G au',/),
     1' Normal mode G-tensor derivatives:')
      DO 20 I=1,ndif
      do 20 j=1,3
      do 20 k=1,3
20    write(9,9005)I,nml(I),j,k,G0Q(j,k,nml(I))
      write(9,9008)(((aten0(i,j,k),i=1,3),j=1,3),k=1,3)
9008  format(9(3f15.7,' A au',/),
     1' Normal mode A-tensor derivatives:')
      DO 21 I=1,ndif
      do 21 j=1,3
      do 21 k=1,3
      do 21 l=1,3
21    write(9,9007)I,nml(I),j,k,l,A0Q(j,k,l,nml(I))
9007  format(5i5,2f13.6)
      close(9)
      write(6,*)'nm pol der in TTTQ.TXT.SCR'
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE WRITETTT09(NAT,ALPHA09,GTENS09,ATENS09,MX3,MFR,ifr,fr)
      IMPLICIT none
      INTEGER*4 MX3,ifr,MFR,I,J,K,L,IX,IIND,IM,IA,N,NQ,NAT,IQ,N7
      real*8 ALPHA09(MFR,9*MX3),GTENS09(MFR,9*MX3),ATENS09(MFR,27*MX3),
     1amu,a,A0(3,3),G0(3,3),AT0(3,3,3),fr
      character*80 filename
      character*10 fchar
      real*8,allocatable::e(:),s(:,:),ALPHAQ(:,:,:),GTENSQ(:,:,:),
     1ATENSQ(:,:,:,:),ALPHA(:,:,:),GTENS(:,:,:),ATENS(:,:,:,:)
      integer*4,allocatable::nml(:)
      logical lex

      write(fchar,'(i10)')ifr
      filename='FILE.TTT.f'//fchar(10:10)
      if(ifr.gt.9)filename='FILE.TTT.f'//fchar(9:10)
      if(ifr.gt.99)filename='FILE.TTT.f'//fchar(8:10)
      if(ifr.gt.999)filename='FILE.TTT.f'//fchar(7:10)
      write(*,*)
      write(*,*)filename
      OPEN(2,FILE=filename)
c     for the first frequency, write the file twice:
      if(ifr.eq.1)OPEN(21,FILE='FILE.TTT')
      WRITE(2,2000)NAT,ifr,fr
2000  FORMAT(' ROA tensors, cartesian derivatives',/,
     1I4,' atoms, freq. ',i2,f11.6/,
     1' The electric-dipolar electric-dipolar polarizability:',/,
     2' Atom/x    jx           jy           jz')
      if(ifr.eq.1)WRITE(21,2000)NAT,ifr,fr
      DO 1 I=1,3
      WRITE(2,2002)I
      if(ifr.eq.1)WRITE(21,2002)I
2002  FORMAT(' Alpha(',I1,',J):')
      DO 1 L=1,NAT
      DO 1 IX=1,3
      IIND=3*(L-1)+IX
      WRITE(2 ,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1)),J=1,3)
1     if(ifr.eq.1)
     1WRITE(21,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1)),J=1,3)
2001  FORMAT(I5,1H ,I1,3g15.7)
      WRITE(2,2004)
      if(ifr.eq.1)WRITE(21,2004)
2004  FORMAT(' The electric dipole magnetic dipole polarizability:',/,
     1       ' Atom/x    jx(Bx)       jy(By)       jz(Bz)')
      DO 2 I=1,3
      WRITE(2,2003)I
      if(ifr.eq.1)WRITE(21,2003)I
2003  FORMAT(' G(',I1,',J):')
      DO 2 L=1,NAT
      DO 2 IX=1,3
      IIND=3*(L-1)+IX
      WRITE(2 ,2001)L,IX,(GTENS09(ifr,I+3*(J-1)+9*(IIND-1)),J=1,3)
2     if(ifr.eq.1)
     1WRITE(21,2001)L,IX,(GTENS09(ifr,I+3*(J-1)+9*(IIND-1)),J=1,3)
      WRITE(2,2005)
      if(ifr.eq.1)WRITE(21,2005)
2005  FORMAT(' The electric dipole electric quadrupole polarizability:',
     2/,     ' Atom/x    kx           ky           kz')
      DO 3 I=1,3
      DO 3 J=1,3
      WRITE(2,2006)I,J
      if(ifr.eq.1)WRITE(21,2006)I,J
2006  FORMAT(' A(',I1,',',I1,',K):')
      DO 3 L=1,NAT
      DO 3 IX=1,3
      WRITE(2,2007)L,IX,
     1(ATENS09(ifr,K+3*(J-1)+9*(I-1)+27*(3*(L-1)+IX-1))
     2*3.d0/2.d0,K=1,3),L,IX,I,J
3     if(ifr.eq.1)WRITE(21,2007)L,IX,
     1(ATENS09(ifr,K+3*(J-1)+9*(I-1)+27*(3*(L-1)+IX-1))
     2*3.d0/2.d0,K=1,3),L,IX,I,J
2007  FORMAT(I5,1H ,I1,3g15.7,' ',4i3)
      write(2,*)
      if(ifr.eq.1)write(21,*)
      write(2,*)'dummy alpha v:'
      if(ifr.eq.1)write(21,*)'dummy alpha v:'
      DO 4 I=1,3
      WRITE(2,2002)I
      if(ifr.eq.1)WRITE(21,2002)I
      DO 4 L=1,NAT
      DO 4 IX=1,3
      IIND=3*(L-1)+IX
      WRITE(2 ,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1)),J=1,3)
4     if(ifr.eq.1)
     1WRITE(21,2001)L,IX,(ALPHA09(ifr,I+3*(J-1)+9*(IIND-1)),J=1,3)
      CLOSE(2)
      if(ifr.eq.1)CLOSE(21)

c     normal mode derivatives in S4 format:
      
      inquire(file='F.INP',exist=lex)
      if(lex.and.ifr.eq.1)then
       N=3*NAT
       allocate(ALPHA(3,3,N),GTENS(3,3,N),ATENS(3,3,3,N),
     1 ALPHAQ(3,3,N),GTENSQ(3,3,N),ATENSQ(3,3,3,N))
       do 201 I=1,3
       do 201 J=1,3
       do 201 IA=1,N
       ATENS(I,J,1,IA)=ATENS09(ifr,1+3*(J-1)+9*(I-1)+27*(IA-1))
       ATENS(I,J,2,IA)=ATENS09(ifr,2+3*(J-1)+9*(I-1)+27*(IA-1))
       ATENS(I,J,3,IA)=ATENS09(ifr,3+3*(J-1)+9*(I-1)+27*(IA-1))
       ALPHA(I,J,  IA)=ALPHA09(ifr,I+3*(J-1)+9*(IA-1))
201    GTENS(I,J,  IA)=GTENS09(ifr,I+3*(J-1)+9*(IA-1))
       amu=1822.88d0
       allocate(s(N,N),e(N),nml(N))
       call readsi(N,S,E,NQ,'F.INP')
       do 101 I=1,3
       do 101 J=1,3
       AT0(I,J,1)=0.0d0
       AT0(I,J,2)=0.0d0
       AT0(I,J,3)=0.0d0
       A0( I,J  )=0.0d0
       G0( I,J  )=0.0d0
       N7=N-NQ+1
       do 101 IQ=1,NQ
       IM=N7+IQ-1
       nml(NQ-IQ+1)=IM
       ATENSQ(I,J,1,IM)=0.0d0
       ATENSQ(I,J,2,IM)=0.0d0
       ATENSQ(I,J,3,IM)=0.0d0
       ALPHAQ(I,J,  IM)=0.0d0
       GTENSQ(I,J,  IM)=0.0d0
       do 101 IA=1,N
       ATENSQ(I,J,1,IM)=ATENSQ(I,J,1,IM)+ATENS(I,J,1,IA)*s(IA,IQ)
       ATENSQ(I,J,2,IM)=ATENSQ(I,J,2,IM)+ATENS(I,J,2,IA)*s(IA,IQ)
       ATENSQ(I,J,3,IM)=ATENSQ(I,J,3,IM)+ATENS(I,J,3,IA)*s(IA,IQ)
       ALPHAQ(I,J,  IM)=ALPHAQ(I,J,  IM)+ALPHA(I,J,  IA)*s(IA,IQ)
101    GTENSQ(I,J,  IM)=GTENSQ(I,J,  IM)+GTENS(I,J,  IA)*s(IA,IQ)
       do 17 IQ=1,NQ
       IM=N7+IQ-1
       a=1.0d0/dsqrt(dabs(amu*e(IQ)))
       do 17 I=1,3
       do 17 J=1,3
       ATENSQ(I,J,1,IM)=ATENSQ(I,J,1,IM)*a*1.5d0
       ATENSQ(I,J,2,IM)=ATENSQ(I,J,2,IM)*a*1.5d0
       ATENSQ(I,J,3,IM)=ATENSQ(I,J,3,IM)*a*1.5d0
       ALPHAQ(I,J,IM)=ALPHAQ(I,J,IM)*a
17     GTENSQ(I,J,IM)=GTENSQ(I,J,IM)*a
       call wrttq(A0,NQ,nml,ALPHAQ,ALPHAQ,G0,AT0,
     1 GTENSQ,ATENSQ,N)
      endif

      RETURN
      END

c     ===============================      
      SUBROUTINE WRITEPOL(APOL,ifr,ty,fr)
      IMPLICIT INTEGER*4 (I-N) 
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 APOL(*)
      character*80 filename
      character*1 fch1
      character*2 fch2
      character*3 fch3
      character*(*) ty

      if(ifr.lt.10)then
       write(fch1,'(i1)')ifr
       filename=ty//'.f'//fch1
      else
       if(ifr.lt.100)then
        write(fch2,'(i2)')ifr
        filename=ty//'.f'//fch2
       else
        if(ifr.lt.1000)then
         write(fch3,'(i3)')ifr
         filename=ty//'.f'//fch3
        else
         write(6,*)' too many frequencies, > 999 skipped !!'
        endif
       endif
      endif

      OPEN(90,FILE=filename)
      if(ty.eq.'A.TTT'.or.ty.eq.'A.s.TTT')then
       write(90,701)fr
701    format(' A, dipole-quadrupole polarizability, w = ',f11.6)
       write(90,902)
       do 1 IZ=1,3
1      write(90,900)((3.0d0/2.0d0*
     1 APOL(IY+(IX-1)*3+(IZ-1)*9),IY=1,IX),IX=1,3)
      else
       if(ty.eq.'POL.TTT'.or.ty.eq.'POL.s.TTT')then
        write(90,702)fr
702     format(' Polarizability, w = ',f11.6)
        write(90,902)
        write(90,900)((APOL(IY+(IX-1)*3),IY=1,IX),IX=1,3)
       else
        if(ty.eq.'ORT.TTT')then
         write(90,703)fr
703      format(' Optical rotation tensor, w = ',f11.6)
         write(90,902)
         write(90,900)((APOL(IY+(IX-1)*3),IY=1,IX),IX=1,3)
        else 
         if(ty.eq.'GP.TTT'.or.ty.eq.'GP.s.TTT')then
          write(90,704)fr
704       format(' G tensor, w = ',f11.6)
          write(90,901)
          write(90,900)((APOL(IY+(IX-1)*3),IY=1,3),IX=1,3)
         else
          write(90,705)fr
705       format(' Unspecified polarizability, w = ',f11.6)
          write(90,901)
          write(90,900)((APOL(IY+(IX-1)*3),IY=1,3),IX=1,3)
         endif
        endif
       endif
      endif
901   format(
     1'            XX            XY            XZ            YX',
     2'            YY            YZ            ZX            ZY',
     3'            ZZ')
902   format(
     1'            XX            XY            YY            XZ',
     2'            YZ            ZZ')
900   format(9G14.6)
      close(90)

      do 3 iend=len(filename),1,-1
3     if(filename(iend:iend).ne.' ')goto 4
4     WRITE(*,*)filename(1:iend)//' written'

c     for f1 write it twice:
      if(ifr.eq.1)then
       OPEN(90,FILE=ty)
       if(ty.eq.'A.TTT'.or.ty.eq.'A.s.TTT')then
        write(90,701)fr
        write(90,902)
        do 2 IZ=1,3
2       write(90,900)((3.0d0/2.0d0*
     1  APOL(IY+(IX-1)*3+(IZ-1)*9),IY=1,IX),IX=1,3)
       else
        if(ty.eq.'POL.TTT'.or.ty.eq.'POL.s.TTT')then
         write(90,702)fr
         write(90,902)
         write(90,900)((APOL(IY+(IX-1)*3),IY=1,IX),IX=1,3)
        else
         if(ty.eq.'ORT.TTT')then
          write(90,703)fr
          write(90,902)
          write(90,900)((APOL(IY+(IX-1)*3),IY=1,IX),IX=1,3)
         else 
          if(ty.eq.'GP.TTT'.or.ty.eq.'GP.s.TTT')then
           write(90,704)fr
           write(90,901)
           write(90,900)((APOL(IY+(IX-1)*3),IY=1,3),IX=1,3)
          else
           write(90,705)fr
           write(90,901)
           write(90,900)((APOL(IY+(IX-1)*3),IY=1,3),IX=1,3)
          endif
         endif
        endif
       endif
       close(90)
       WRITE(*,*)ty//' written'
      endif

      return

      end

      function FU(NUMBER)
      implicit none
      real*8 FU,x
      CHARACTER*1 NUMBER(*)
      CHARACTER*20 s20
      integer*4 i
      do 1 i=1,20
1     s20(i:i)=NUMBER(i)
      read(s20,*,err=99)x
      FU=x
      return
99    write(6,*)s20
      call report('format error')
      end
      
      function mfrdet(FILENAME)
      implicit none
      integer*4 mfrdet,i,n
      character*99 s99
      character*(*) FILENAME
      n=0
      open(7,FILE=FILENAME,STATUS='OLD',FORM='FORMATTED')
1     read(7,700,end=99,err=99)s99
700   format(a99)
      if(s99(2:32).eq.'Using perturbation frequencies:')then
       do 2 i=0,4
2      if(s99(38+i*12:38+i*12).eq.'.')n=n+1
      endif
      if(s99(2:4).eq.'   '.and.n.gt.0)goto 99
      goto 1
99    close(7)
      mfrdet=n
      write(6,*)n,' perturbation frequencies'
      return
      end
     
      subroutine mfrdetw(FILENAME,n0)
      implicit none
      integer*4 i,n,n0
      character*99 s99
      character*(*) FILENAME
      real*8,allocatable::e(:)
      allocate(e(n0))
      n=0
      open(7,FILE=FILENAME,STATUS='OLD',FORM='FORMATTED')
1     read(7,700,end=99,err=99)s99
700   format(a99)
      if(s99(2:32).eq.'Using perturbation frequencies:')then
       do 2 i=0,4
       if(s99(38+i*12:38+i*12).eq.'.')then
        n=n+1
        read(s99(38+i*12-5:38+i*12+6),*)e(n)
       endif
2      continue
      endif
      if(s99(2:4).eq.'   '.and.n.gt.0)goto 99
      goto 1
99    close(7)
      open(7,file='EP.LST')
      write(7,*)n,' perturbation frequencies'
      write(7,701)(e(i),i=1,n)
701   format(5f12.6)
      close(7)
      write(6,*)n,' perturbation frequencies in EP.LST'
      return
      end

      subroutine svgr(g,n)
      implicit none
      real*8 g(*)
      integer*4 n,i,ix
      open(70,file='FILE.GR')
      do 1 i=1,n
1     write(70,700)(g(ix+3*(i-1)),ix=1,3)
700   format(3f15.9)
      close(70)
      write(6,*)'Gradient written to FILE.GR'
      return
      end
c     ============================================================
      SUBROUTINE readsi(N3,S,E,NQ,fn)
c     switch order of modes to match Gaussian convention:
      IMPLICIT none
      integer*4 N3,NQ,nat,i,J,ix,iz
      real*8 S(N3,N3),E(*),CM
      character*(*) fn
      CM=219474.630d0
C     1 a.u.= 2.1947E5 cm^-1
      open(4,file=fn)
      read(4,*)NQ,nat,nat
      do 1 i=1,NAT
1     read(4,*)
      read(4,*)
      DO 2 I=1,NAT
      DO 2 J=1,NQ
2     read(4,*)(s(3*(i-1)+ix,NQ-J+1),ix=1,2),
     1(s(3*(i-1)+ix,NQ-J+1),ix=1,3)
      read(4,*)
      READ(4,4000)(E(NQ-J+1),J=1,NQ)
4000  FORMAT(6F11.6)
      close(4)
c
      write(6,*)NQ,' modes found'

      iz=0
c     delete zero modes if exist:
66    do 6 i=1,NQ
      if(dabs(E(i)).lt.0.1d0)then
       do 7 j=i,NQ-1
       E(j)=E(j+1)
       do 7 ix=1,N3
7      s(ix,j)=s(ix,j+1)
       iz=iz+1
       NQ=NQ-1
       goto 66
      endif
6     continue
      if(iz.gt.0)write(6,*)iz,' zero modes: deleted'

      write(6,*)NQ,' vibrational modes considered'

      DO 3 I=1,NQ
3     E(I)=E(I)/CM
      
      RETURN
      end
c     ============================================================
      subroutine  rgt(N,L,M,t)
c     read a N x M matrix from Gaussian output
      implicit none
      integer*4 N,M,i,j,L,ii,j1,j2
      real*8 t(*)
      j1=1
2     j2=min(j1+4,M)
      read(7,*)
      do 1 i=1,N
      do 1 ii=0,L-1
1     read(7,*)t(i+N*(j1-1)+N*M*ii),(t(i+N*(j-1)+N*M*ii),j=j1,j2)
c     order
c     alpha:
c     x x d/dx1
c     x y d/dx1
c     x z d/dx1
c     ...
c     z z d/dzNat
      if(j2.lt.M)then
       j1=j2+1
       goto 2
      endif
      return
      end
c     ============================================================
      subroutine trt(t,tm,N,M,s,L)
c     transform first index Cartesian -> normal modes
      implicit none
      integer*4 N,M,i,j,k,L,ii
      real*8 t(*),tm(*),s(N,N)
      do 1 ii=0,L-1
      do 1 i=1,M
      do 1 j=1,M
      tm(i+M*(j-1)+M*M*ii)=0.0d0
      do 1 k=1,N
1     tm(i+M*(j-1)+M*M*ii)=tm(i+M*(j-1)+M*M*ii)
     1+s(k,i)*t(k+N*(j-1)+N*M*ii)
      return
      end
c     ============================================================
      subroutine wrttqaq(ndif,N,ALPHAA,GGA,AAA,MFR,ifr,w,N7)
      implicit none
      integer*4 ndif,i,j,N,ix,k,l,ifr,MFR,N7
      real*8 ALPHAA(N**2,9,MFR),GGA(N**2,9,MFR),AAA(N**2,27,MFR),w
      character*80 filename
      character*10 fchar
      write(fchar,'(i10)')ifr
      if(ifr.eq.1)then
       filename='TTTAQ.TXT.SCR'
      else
       filename='TTTAQ.TXT.SCR.f'//fchar(10:10)
       if(ifr.gt.9)filename='TTTAQ.TXT.SCR.f'//fchar(9:10)
       if(ifr.gt.99)filename='TTTAQ.TXT.SCR.f'//fchar(8:10)
       if(ifr.gt.999)filename='TTTAQ.TXT.SCR.f'//fchar(7:10)
      endif
      write(6,*)filename
      open(21,file=filename)
      write(21,211)ifr,w
211   format(' NM pol sec derivatives iw:',i2,' w: ',f10.1,' cm-1')
      write(21,*)ndif
      do 23 i=1,ndif
      do 23 j=i,ndif
      k=ndif-i+1
      l=ndif-j+1
      write(21,2124)i,j,k+N7-1,l+N7-1
2124  format(4i5,' modes')
23    write(21,2121)(ALPHAA(k+N*(l-1),ix,ifr),ix=1,9)
2121  format(3f15.6)
      do 24 i=1,ndif
      do 24 j=i,ndif
      k=ndif-i+1
      l=ndif-j+1
      write(21,2124)i,j,k+N7-1,l+N7-1
24    write(21,2121)(GGA(k+N*(l-1),ix,ifr),ix=1,9)
      do 25 i=1,ndif
      do 25 j=i,ndif
      k=ndif-i+1
      l=ndif-j+1
      write(21,2124)i,j,k+N7-1,l+N7-1
25    write(21,2121)(AAA(k+N*(l-1),ix,ifr),ix=1,27)
      close(21)
      write(6,*)'TTTAQ.TXT.SCR written'
      end
c     ============================================================
      subroutine s8
      integer*4 i
      do 2 i=1,8
2     read(7,*)
      return
      end
c     ============================================================
      subroutine vz(v,n)
      real*8 v(*)
      integer*4 i,n
      do 1 i=1,n
1     v(i)=0.0d0
      return
      end
c     ============================================================
      function rdgeo(filename,NAT,rs,rz,ic,iz,atsy)
c     read geometry, standard or/and z-matrix
c     ic = 0  .. only fin dth enumbe rof atoms
c          1  .. also read the geometry
      implicit none
      character*2 atsy(89)
      logical rdgeo,ls,lz,lsnow,lznow,loff,lchk
      character*(*) filename
      character*50 FN
      character*80 s80
      integer*4 NAT,ic,NGS,NGZ,l,ig98,IX,KA,ica,i,ic1,ic2
      real*8 rs(*),rz(*),x,y,z
      integer*4 iz(*)
      NGS=0
      NGZ=0
      ls=.false.
      lz=.false.
      loff=.false.
      OPEN(2,FILE=filename)
1     READ(2,2000,END=1000)FN
2000  FORMAT(A50)
      if(FN(2:20).EQ.'Symmetry turned off')loff=.true.
      lsnow=FN(20:40).EQ.'Standard orientation:'.OR.
     2      FN(26:46).EQ.'Standard orientation:'
      lchk=FN(2:37).EQ.'Redundant internal coordinates found'
      lznow=FN(19:39).EQ.'Z-Matrix orientation:'.OR.
     1      FN(26:46).EQ.'Z-Matrix orientation:'.OR.
     1      FN(20:37).EQ.'Input orientation:'.OR.
     1      FN(27:44).EQ.'Input orientation:'.OR.lchk
      if(lsnow.or.lznow)then
       if(lsnow)NGS=NGS+1
       if(lznow)NGZ=NGZ+1
       ls=ls.or.lsnow
       lz=lz.or.lznow

       l=0
c      Z-matrix from checkpoint, listed in output
       if(lchk)then
283     read(2,280)s80
280     format(a80)
        if(s80(3:3).ne.','.and.s80(4:4).ne.',')goto 281
        l=l+1
        ica=0
        ic1=0
        ic2=0
        do 282 i=1,len(s80)
        if(s80(i:i).eq.',')then
         ica=ica+1
         s80(i:i)=' '
c        remember position of the first and second comma:
         if(ica.eq.1)ic1=i
         if(ica.eq.2)ic2=i
        endif
282     continue
        read(s80(ic2+1:len(s80)),*)x,y,z
        KA=0
        do 284 i=1,89
284     IF(s80(2:ic1-1).EQ.atsy(i))KA=i
        IF(KA.EQ.0)call report(' Unknown atom '//s80(2:ic1-1))
        if(ic.eq.1)then
         iz(l)=KA
         rz(3*(l-1)+1)=x
         rz(3*(l-1)+2)=y
         rz(3*(l-1)+3)=z
        endif
        goto 283
281     continue
       else
c       gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
c       Z-matrix further from output:
        ig98=0
        if(FN(26:46).EQ.'Z-Matrix orientation:'.OR.
     1     FN(27:44).EQ.'Input orientation:'.OR.
     1     FN(26:46).EQ.'Standard orientation:')ig98=1
        READ(2,*)
        READ(2,*)
        READ(2,*)
        READ(2,*)
5       READ(2,2000)FN
        IF(FN(2:4).NE.'---')THEN
         l=l+1
         BACKSPACE 2
         if(ig98.eq.0)then
          READ(2,*)KA,KA,x,y,z
         else
          READ(2,*)KA,KA,x,x,y,z
         endif
c        ignore dummy atoms:
         IF(KA.EQ.-1)then
          l=l-1
         else
          if(ic.eq.1)then
           iz(l)=KA
           if(lsnow)then
            rs(3*(l-1)+1)=x
            rs(3*(l-1)+2)=y
            rs(3*(l-1)+3)=z
           else
            rz(3*(l-1)+1)=x
            rz(3*(l-1)+2)=y
            rz(3*(l-1)+3)=z
           endif
          endif
         endif
         GOTO 5
        ENDIF
c      gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
       ENDIF
       nat=l

      ENDIF
      GOTO 1
1000  CLOSE(2)
      if(NGZ.eq.0)call report('z-matrix geometry not found')
      if(NGS.eq.0)then
       if(loff)then
        if(ic.eq.0)write(6,*)'standard orientation not used'
        if(ic.eq.1)then
         do 2 l=1,3*nat
2        rs(l)=rz(l)
        endif
       else
        call report('Standard geometry not found')
       endif
      endif

      if(ic.eq.0)then
       write(6,*)NGZ,' Z-matrix geometries'
       write(6,*)NGS,' standard geometries'
       write(*,*)nat,' atoms'
      endif

      rdgeo=NGZ.gt.0

      if(ic.eq.1)then
       OPEN(15,FILE='FILE.X')
       WRITE(15,*)'Last Z-matrix Geometry from Gaussian output'
       WRITE(15,*)nat
       DO 777 l=1,nat
777    WRITE(15,1501)iz(l),(rz(3*(l-1)+IX),IX=1,3)
1501   FORMAT(I4,3F15.8,' 0 0 0 0 0 0 0 0.0')
       CLOSE(15)
       WRITE(*,*)' X written into FILE.X'
      endif

      return
      END
c     ============================================================
      subroutine report(s)
      character*(*) s
      write(6,*)s
      stop
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
      subroutine xst(nat,r1,r2,t,o)
c     standard orientation of r2 with respect to r1
c     J. Phys. Chem. A 2001, 105, 5326-5333
      implicit none
      integer*4 nat,i,j,ia,IERR,ix,iy,iz,ii,ic
      real*8 c(3,3),cc(3,3),r(3,3),l(3),ci(3,3),r1(*),r2(*),
     1TOL,e(3,6),lrc(3,3),bl(3),rc(3,3),err,t0(3,3),det,t(3,3),em,
     1cmg(3),cme(3),o(3)
      real*8,allocatable::rn(:),rg(:),re(:)
      TOL=1.0d-10
      allocate(rg(3*nat),re(3*nat))
      do 3 i=1,3*nat
      rg(i)=r1(i)
3     re(i)=r2(i)
c     calculate geometrical centers:
      do 14 i=1,3
      cmg(i)=0.0d0
      cme(i)=0.0d0
      do 16 ia=1,nat
      cmg(i)=cmg(i)+rg(i+3*(ia-1))
16    cme(i)=cme(i)+re(i+3*(ia-1))
      cmg(i)=cmg(i)/dble(nat)
14    cme(i)=cme(i)/dble(nat)
      o(1)=cme(1)-cmg(1)
      o(2)=cme(2)-cmg(2)
      o(3)=cme(3)-cmg(3)
c     center coordinates:
      do 15 ia=1,nat
      do 15 i=1,3
      rg(i+3*(ia-1))=rg(i+3*(ia-1))-cmg(i)
15    re(i+3*(ia-1))=re(i+3*(ia-1))-cme(i)
c     form C =X X':
      do 1 i=1,3
      do 1 j=1,3
      c(i,j)=0.0d0
      do 1 ia=1,nat
1     c(i,j)=c(i,j)+rg(i+3*(ia-1))*re(j+3*(ia-1))
c     correction for planar molecules:
      do 2 i=1,3
2     if(dabs(c(i,i)).lt.1.0d-9)c(i,i)=1.0d0
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
c     transform re to rn:
      call trx(nat,rn,re,t0)
      err=0.0d0
      do 8 ii=1,3*nat
8     err=err+(rn(ii)-rg(ii))**2

c     write(6,*)'rg:'
c     do ia=1,nat
c     write(6,304)(rg(i+3*(ia-1)),i=1,3)
c04   format(3F12.6)
c     enddo
c     write(6,*)'re:'
c     do ia=1,nat
c     write(6,304)(re(i+3*(ia-1)),i=1,3)
c     enddo
c     write(6,*)'rn:'
c     do ia=1,nat
c     write(6,304)(rn(i+3*(ia-1)),i=1,3)
c     enddo

c     write(6,604)t,o,err

      if(ic.eq.1)then
       em=err
       t=t0
      endif
      if(det(t0).gt.0.0d0.and.err.lt.em)then
       em=err
       t=t0
      endif
4     continue
      write(6,604)t,o,em
604   format(' Transformation matrix:',/,3(3f10.3,/),
     1       ' Center distance      :',/,  3f10.3  ,/,
     1       ' Error                :'  ,   f10.3  )
      return
      end
c     ==============================================================
      SUBROUTINE TRED12(N,A,Z,D,IEIG,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N),Z(N,N),D(*)
      real*8,allocatable::E(:)
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
c     ==============================================================
C
      SUBROUTINE TQL12(N,Z,D,IERR,ICODE,E)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION Z(N,N),D(*),E(*)
      REAL*8 MACHEP
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
c     ==============================================================
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

c     ============================================================
      SUBROUTINE TG(AL,G,o)
      IMPLICIT none
      integer*4 IA
      real*8 t1,t2,t3, AL(9),G(9),o(3)
      t1=-o(1)/0.529177D0/2.0d0
      t2=-o(2)/0.529177D0/2.0d0
      t3=-o(3)/0.529177D0/2.0d0
      DO 2 IA=1,3
c     G(IA+3*(1-1))=G(IA+3*(1-1))+t2*AL(IA+6)-t3*AL(IA+3)
c     G(IA+3*(2-1))=G(IA+3*(2-1))+t3*AL(IA  )-t1*AL(IA+6)
c     G(IA+3*(3-1))=G(IA+3*(3-1))+t1*AL(IA+3)-t2*AL(IA  )
      G(1+3*(IA-1))=G(1+3*(IA-1))+t2*AL(IA+6)-t3*AL(IA+3)
      G(2+3*(IA-1))=G(2+3*(IA-1))+t3*AL(IA  )-t1*AL(IA+6)
2     G(3+3*(IA-1))=G(3+3*(IA-1))+t1*AL(IA+3)-t2*AL(IA  )
      return
      end
c     ============================================================
      SUBROUTINE TA(AL,A,o)
      IMPLICIT none
      integer*4 IA,IB,IC,ii
      real*8 t(3),AL(9),A(27),o(3),f
      f=2.0d0/3.0d0
      t(1)=-o(1)/0.529177D0*f
      t(2)=-o(2)/0.529177D0*f
      t(3)=-o(3)/0.529177D0*f
      DO 3 IA=1,3
      DO 3 IB=1,3
      DO 3 IC=1,3
c     IB IC ... quadrupole indices
      ii=IB+3*(IC-1)+9*(IA-1)
      if(IB.EQ.IC)A(ii)=A(ii)+t(1)*AL(IA)+t(2)*AL(IA+3)+t(3)*AL(IA+6)
3     A(ii)=A(ii)-1.5d0*(t(IB)*AL(IA+3*(IC-1))+t(IC)*AL(IA+3*(IB-1)))
      RETURN
      END
c     ============================================================
      subroutine smooth(t,n)
      implicit none
      integer*4 n,i
      real*8 t(*)
      do 1 i=1,n
1     if(dabs(t(i)).lt.1.0d-6)t(i)=0.0d0
      return
      end
c     ============================================================
      subroutine rewr(ic,A,n,AW,MFR,ifr)
      implicit none
      integer*4 ic,n,ifr,i,MFR
      real*8 A(*),AW(MFR,n)
      if(ic.eq.0)then
       do 1 i=1,n
1      A(i)=AW(ifr,i)
      else
       do 2 i=1,n
2      AW(ifr,i)=A(i)
      endif
      end
c     ============================================================
      subroutine wrnmr(nat,nmr,iz)
      implicit none
      integer*4 nat,i,iz(*)
      real*8 nmr(*)
      open(9,file='FILE.NMR')
      write(9,900)nat
900   format(i5,'   NMR isotropic shieldings',/,
     1       ' atom    sigma')
      write(9,901)
901   format(60(1h-))
      do 1 i=1,nat
1     write(9,902)i,iz(i),nmr(i)
902   format(i5,i3,f12.4)
      write(9,901)
      close(9)
      write(6,*)'FILE.NMR written'
      return
      end
c     ============================================================
      subroutine rspin(s,N)
      implicit none
      integer*4 N,N1,N3,LN,I,J
      real*8 s(N,N)
      N1=1
111   N3=N1+4
      IF(N3.GT.N)N3=N
      read(7,*)
      DO 130 LN=N1,N
130   READ(7,*)s(LN,N1),(s(LN,J),J=N1,MIN(LN,N3))
      N1=N1+5
      IF(N3.LT.N)GOTO 111
      DO 31 I=1,N
      DO 31 J=I,N
   31 s(I,J)=s(J,I)
      return
      end
c     ============================================================
      subroutine wspin(s,N)
      implicit none
      integer*4 N,N1,N3,LN,J
      real*8 s(N,N)
      open(88,file='FILE.SPI')
      write(88,880)
880   format(' NMR spin-spin coupling constants, Hz')
      N1=1
111   N3=N1+4
      IF(N3.GT.N)N3=N
      write(88,882)(J,J=N1,N3)
882   format(3x,5i14)
      DO 130 LN=N1,N
130   write(88,881)LN,(s(LN,J),J=N1,MIN(LN,N3))
881   format(i7,5d14.6)
      N1=N1+5
      IF(N3.LT.N)GOTO 111
      close(88)
      write(6,*)'FILE.SPI written'
      return
      end
c     ============================================================
c jk-end
C
C\Version=SGI-G94RevB.3\HF=-966.6745649\RMSD=0.000e+00\RMSF=1.196e-01\D
Cipole=0.,0.,0.\PG=C01 [X(C38H52N16O24P4)]\NImag=84\\0.76752455,0.00053
C432,0.03734306,0.23733561,0.00042079,0.13291623,-0.70920516,0.00065060
C,-0.21702343,1.33224819,0.00170226,-0.03103517,0.00051274,0.02526703,1
C.09059550,-0.21867837,0.00010553,-0.09804719,0.00890088,-0.02513206,1.
