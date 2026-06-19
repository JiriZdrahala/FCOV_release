      PROGRAM NEW1
C     THIS PROGRAM READS .XYZ FILE (CARTESIAN COORDINATES OF A MOLECUL
C     IN THE PMODEL FORMAT) AND MAKES THE UMATIN INPUT FILES FOR THE
C     VCD PROGRAMS
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      integer*4 ISTRE,IBEND,ITOR,N,IOP,MENDELEV,TBUF(6),narg
c     N7 . maximum number of bonds on one atom
      PARAMETER (MENDELEV=89,narg=11,N7=10)
      real*8 RAD(MENDELEV),u(3,3)
      COMMON/DUMMY/RAD0,RAD
      logical l(narg)
      CHARACTER*80 MOLNAM
      CHARACTER*120 s120
      CHARACTER*1 st
      CHARACTER*5 PGROUP
      CHARACTER*10 sr
      CHARACTER*2 atsy(MENDELEV),mm(6)
      data atsy/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     3'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     4' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     4          'Ga','Ge','As','Se','Br','Kr',
     5'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     5          'In','Sn','Sb','Te',' I','Xe',
     6'Cs','Ba','La',
     6               'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',
     6               'Er','Tm','Yb','Lu',
     6'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     6          'Tl','Pb','Bi','Po','At','Rn',
     7 'Fr','Ra','Ac'/
      data mm/'Tx','Ty','Tz','Rx','Ry','Rz'/
      dimension bonding(MENDELEV)
      data bonding/0.32,0.98,1.28,0.95,0.87,0.82,0.80,0.78,0.77,0.76,
     1             1.59,1.41,1.23,1.16,1.11,1.20,1.04,1.03,2.08,1.79,
     2             1.49,1.37,1.27,1.23,1.22,1.22,1.21,1.20,1.22,1.30,
     3             1.31,1.27,1.25,1.21,1.19,1.17,2.21,1.96,1.67,1.50,
     4             1.39,1.35,1.32,1.30,1.30,1.33,1.39,1.53,1.49,1.46,
     5             1.45,1.41,1.38,1.36,2.40,2.03,1.74,1.70,1.70,1.69,
     6             1.68,1.67,1.90,1.66,1.64,1.64,1.63,1.62,1.61,1.75,
     7             1.61,1.49,1.39,1.35,1.33,1.31,1.32,1.35,1.39,1.54,
     8             1.53,1.52,1.51,1.51,1.50,1.49,1.50,1.50,1.50/
c                     H   He   Li   Be    B    C    N    O    F   Ne
c                     Na  Mg   Al   Si    P    S   Cl   Ar    K   Ca   
c                     Sc  Ti    V   Cr   Mn   Fe   Co   Ni   Cu   Zn
c                     Ga  Ge   As   Se   Br   Kr   Rb   Sr    Y   Zr
c                     Nb  Mo   Tc   Ru   Rh   Pd   Ag   Cd   In   Sn 
c                     Sb  Te    I   Xe   Cs   Ba   La   Ce   Pr   Nd
c                     Pm  Sm   Eu   Gd   Tb   Dy   Ho   Er   Tm   Yb  
c                     Lu   Hf   Ta    W   Re  Os   Ir   Pt   Au   Hg 
c                     Tl   Pb   Bi   Po   At  Rn   Fr   Ra   Ac
c     helix coordinates:
      parameter (ih0=10)
      character*10 HC(ih0)
      data HC/' H-TX     ',' H-TY     ',' H-TZ     ',' H-RX     ',
     1        ' H-RY     ',' H-RZ     ',' H-stretch',' H-breath.',
     1        ' H-def.   ',' H-unwind.'/
      character*1 str(narg)
      data str/'b','a','A','B','H','M','r','p','h','t','s'/
      character*3,allocatable::stype(:)
      real*8,allocatable::di(:),X(:),Y(:),Z(:),Xg(:),Yg(:),Zg(:)
      integer*4,allocatable::nlist(:),mlist(:,:),ilist(:),blist(:,:),
     1hl(:),NBT(:),NBTg(:),BT(:,:),IBONDg(:),ICON(:,:),ICONg(:,:),
     1IST(:),JST(:),IBD(:),JBD(:),NITOR(:),NLTOR(:),JTOR(:),KTOR(:),
     3IBUF(:,:),JBUF(:,:),IOOP(:),JOOP(:),KOOP(:),LOOP(:),KBD(:),
     1IBOND(:),KATOM(:),BTg(:,:)
      CHARACTER*2,allocatable:: ATOMT(:)

      WRITE(*,*)
      WRITE(*,*)' VCD PROGRAM PACKAGE INPUT FILE GENERATION '
      WRITE(*,*)
      WRITE(*,*)'      PETR BOUR, UOCHB CSAV 1992-96'
      WRITE(*,*)
      WRITE(*,*)'       FILE X(YZ)  -->   FILE UMAT'

      l=.false.
      rl0=0.0d0
      I=1
2     if(I.le.iargc())then
       call getarg(I,st)
       do 21 j=1,narg
21     if(st.eq.str(j))l(j)=.true.
       if(st.eq.'r')then
        call getarg(I+1,sr)
        read(sr,*)rl0
        I=I+1
       endif 
       I=I+1
       goto 2
      endif

      if(l(3))l(2)=.true.
      if(l(8))then
       open(8,file='CELL.TXT',status='old')
       read(8,*)(u(1,i),i=1,3)
       read(8,*)(u(2,i),i=1,3)
       read(8,*)(u(3,i),i=1,3)
       close(8)
      endif
c                  a    A    M    b    B    H    r     p   h    s
      WRITE(6,6001)l(2),l(3),l(6),l(1),l(4),l(5),l(7),l(8),l(9),l(11),
     1l(10)
c      t
6001  format('Inline options: a ... generate all possible coords ',l2,/,
     1       '                A ... as a with intermolecular     ',l2,/,
     1       '                M ... monomolecular coords         ',l2,/,
     2       '                b ... use bond table in FILE.X     ',l2,/,
     3       '                B ... use bond table for intermol  ',l2,/,
     4       '                H ... helix coords for HA.LST      ',l2,/,
     5       '                      (with a only)                ',   /,
     6       '                r <rlim> ... redef. intermol rlim  ',l2,/,
     6       '                p ... periodicity, use CELL.TXT    ',l2,/,
     6       '                h ... include H-bonds              ',l2,/,
     6       '                s ... sort                         ',l2,/,
     6       '                t ... full torsion                 ',l2,/,
     6       '                n ... generate non red. (=default)',/)
      if(l(7))write(6,6003)rl0
6003  format(' Limit for intermolecular interactions ',f10.2,' A')

      IMM=0
      IHH=0
      nha=0
      IERR=0
      M1=0

      OPEN(55,FILE='FILE.X',STATUS='OLD')
      READ(55,888)MOLNAM
888   format(a80)
      read(55,*)N
c     maximal number of stretches ~ 4x Nat:
      if(l(8))then
       ns0=4*N*27
      else
       ns0=4*N
      endif
      allocate(IST(ns0),jst(ns0),BT(N,N7),NBT(N),IBD(ns0),JBD(ns0),
     1NITOR(2*ns0),NLTOR(2*ns0),JTOR(2*ns0),KTOR(2*ns0),IBUF(2*ns0,6),
     1JBUF(2*ns0,6),
     4IOOP(ns0),JOOP(ns0),KOOP(ns0),LOOP(ns0),ICON(N+2,8),
     1KBD(ns0),IBOND(N),ATOMT(27*N),KATOM(27*N),X(N),Y(N),Z(N))
      allocate(stype(ns0),di(ns0))
      stype=' - '
      BT=0
      im=0
      do 1111 I=1,N
c     for the first atom see if bond table or not:
      if(I.eq.1)then
       READ(55,120)s120
120    format(a120)
       ip=0
       do 1112 j=1,len(s120)
       if(s120(j:j).eq.'.')ip=ip+1
       if(ip.eq.3)then
        do 11113 k=j,len(s120)-1
11113   if(s120(k:k).eq.' '.and.s120(k+1:k+1).ne.' ')im=im+1
        goto 1114
       endif
1112   continue
1114   backspace 55
      endif
      if(im.ge.7)then
       READ(55,*)KATOM(I),X(I),Y(I),Z(I),(BT(I,J),J=1,7)
      else
       READ(55,*)KATOM(I),X(I),Y(I),Z(I)
      endif
1111  ATOMT(I)=ATSY(KATOM(I))
      CLOSE(55)
c
      RAD0=0.25d0
      DO 1113 I=1,89
1113  RAD(I)=bonding(I)
      RAD(1)=0.4d0
      RAD(6)=0.83d0
      RAD(9)=0.75d0
      RAD(12)=1.5d0
      RAD(14)=1.5d0
      RAD(15)=1.5d0
      RAD(16)=1.2d0
      RAD(17)=1.5d0
      RAD(35)=1.5d0
      RAD(80)=2.2d0
      
C     ************************* BOND GENERATION: *********************
      nb=0
      nbg=0
      IBOND=0
      ICON=0
      NBT=0

      if(l(1))then
       write(6,*)'Bond generated from bondtable'
       do 1117 I=1,N
       do 1117 J=1,N7
 1117  if(BT(I,J).ne.0)NBT(I)=NBT(I)+1
       DO 12 I=1,N
       DO 12 k=1,NBT(I)
       J=BT(I,k)
       if(J.gt.I)then
        nb=nb+1
        IBOND(I)=IBOND(I)+1
        IBOND(J)=IBOND(J)+1
        ICON(I,IBOND(I))=J
        ICON(J,IBOND(J))=I
       ENDIF
12     CONTINUE
      else
       call genb(N,X,Y,Z,KATOM,RAD0,RAD,NBT,BT,nb,IBOND,ICON,l(9),N7)
       write(6,*)nb,' bonds generated from distances'

       if(l(8))then
c       3x3x3 geometry
        Ng=27*N
        allocate(Xg(Ng),Yg(Ng),Zg(Ng),NBTg(Ng),IBONDg(Ng),
     1  ICONg(Ng+2,8),BTg(Ng,N7))
        do 13 I=1,N
        Xg(I)=X(I)
        Yg(I)=Y(I)
13      Zg(I)=Z(I)
        ia=N
        do 11 i1=-1,1
        do 11 i2=-1,1
        do 11 i3=-1,1
        if(i1.ne.0.or.i2.ne.0.or.i3.ne.0)then
         do 111 I=1,N
         ia=ia+1
         Xg(ia)=X(I)+dble(i1)*u(1,1)+dble(i2)*u(2,1)+dble(i3)*u(3,1)
         Yg(ia)=Y(I)+dble(i1)*u(1,2)+dble(i2)*u(2,2)+dble(i3)*u(3,2)
         Zg(ia)=Z(I)+dble(i1)*u(1,3)+dble(i2)*u(2,3)+dble(i3)*u(3,3)
         ATOMT(ia)=ATSY(KATOM(I))
111      KATOM(ia)=KATOM(I)
        endif
11      continue
        nbg=0
        IBONDg=0
        ICONg=0
        NBTg=0
        call genb(Ng,Xg,Yg,Zg,KATOM,RAD0,RAD,NBTg,BTg,nbg,IBONDg,ICONg,
     1  l(9),N7)
        write(6,*)nbg,' bonds in supercube'
        deallocate(X,Y,Z)
        allocate(X(Ng),Y(Ng),Z(Ng))
        X=Xg
        Y=Yg
        Z=Zg
       endif

      endif

      write(6,*)nb,' bonds'

      if(l(2).or.l(3))then
       write(6,*)'all possible coordinates generated'

       if(l(8))then
        call genstr(ISTRE,ns0,NBTg,Ng,N,BTg,IST,JST,N7)
        call genbe( IBEND,ns0,NBTg,Ng,N,BTg,IBD,JBD,KBD,N7)
        call gento(  ITOR,ns0,NBTg,Ng,N,BTg,NITOR,IBUF,JTOR,KTOR,
     1  NLTOR,JBUF,l(10),N7)
        call giop( IOP   ,ns0,NBTg,Ng,N,BTg,IOOP,JOOP,KOOP,LOOP,
     1  Xg,Yg,Zg,N7)
       else
        call genstr(ISTRE,ns0,NBT,N,N,BT,IST,JST,N7)
        call genbe(IBEND ,ns0,NBT,N,N,BT,IBD,JBD,KBD,N7)
        call gento(ITOR  ,ns0,NBT,N,N,BT,NITOR,IBUF,JTOR,KTOR,
     1  NLTOR,JBUF,l(10),N7)
        call giop( IOP   ,ns0,NBT,N,N,BT,IOOP,JOOP,KOOP,LOOP,
     1  X,Y,Z,N7)
       endif

        


c      divide system into molecules:
       if(l(3).or.l(6))then
        allocate(mlist(N,N),nlist(N))
        call gmol(N,NMOL,X,Y,Z,KATOM,nlist,mlist,l,BT,N7)
       endif

       if(l(3))then
c       for each molecule, find neighboring ones
        if(NMOL.gt.0)then
         allocate(blist(NMOL,NMOL),ilist(NMOL))
         call bmol(NMOL,X,Y,Z,N,nlist,mlist,blist,ilist,KATOM,mi,rl0)
        endif
        IMM=6*mi
       endif

       if(l(6))M1=6*NMOL

       if(l(5))then
c       helical/cigar coordinates:
c       determine number of helix atoms:
        call hmol1(nha,N)
        allocate(hl(nha))
c       determine helix atoms:
        call hmol2(hl,nha,N)
        IHH=ih0
       endif

       NTOT=ISTRE+IBEND+ITOR+IOP+IMM+IHH+M1
       call overview(ISTRE,IBEND,ITOR,IOP,M1,IMM,IHH,NTOT,N)

       goto 1001
      endif
      
C     ENDIF
C     ************************** COORDINATE DEFINITION: **************
      CALL INCO(IERR,N,ns0,ICON,ISTRE,IBEND,ITOR,IST,JST,IBD,
     1JBD,KBD,NITOR,NLTOR,JTOR,KTOR,IBUF,JBUF,IBOND,IOP)
C     ************************************************************
      IF(IERR.EQ.1)THEN
       WRITE(*,*)' Internal coordinates cannot be used'
       NTOT=3*N
       GOTO 1001
      ENDIF
      NTOT=ISTRE+IBEND+ITOR+IOP
      call overview(ISTRE,IBEND,ITOR,IOP,M1,IMM,IHH,NTOT,N)
      IF(NTOT.NE.3*N-6)THEN
       WRITE(*,*)' ** NUMBER OF FOUND COORDINATES IS NOT 3N-6  !! **'
       WRITE(*,*)' **     Manual redefinition required **'
      ENDIF
      WRITE(*,*)
C     **************************** OUTPUT: ***********************
1001  WRITE(*,*)'FILE.UMA - IS THE OUPTUT FILE'
      IPNCH=-1
      PGROUP='C1   '
      OPEN(7,FILE='FILE.UMA',STATUS='UNKNOWN')
      WRITE(7,888)MOLNAM
      WRITE(7,888)MOLNAM
      NFCS=((3*N-6)*((3*N-6)+1))/2
      WRITE(7,222)N,NTOT,NFCS,IPNCH,PGROUP
222   FORMAT(4I12,A5)
      DO 100 I=1,N
100   WRITE(7,3333)X(I),Y(I),Z(I),I,KATOM(I),ATOMT(I)
3333  FORMAT(3F12.6,2I4,A2)
      IF(IERR.EQ.1)GOTO 1002

      if(l(2).or.l(1).or.l(11))then
       write(6,*)'Sorting coordinates'

       DO 1041 I=1,ISTRE
       ia=KATOM(abs(IST(I)))
       ja=KATOM(abs(JST(I)))
1041   if(ja.lt.ia)call iswp(JST(I),IST(I))

c      stretch, first light atoms
       DO 104 I=1,ISTRE
       do 104 J=I+1,ISTRE
       ia=KATOM(abs(IST(I)))
       ja=KATOM(abs(IST(J)))
       if(ja.lt.ia)then
        call iswp(IST(J),IST(I))
        call iswp(JST(J),JST(I))
       ENDIF
104    continue

c      stretch, if firts same, then first heavy atoms
       DO 105 I=1,ISTRE
       do 105 J=I+1,ISTRE
       ib=KATOM(abs(IST(I)))
       jb=KATOM(abs(IST(J)))
       ia=KATOM(abs(JST(I)))
       ja=KATOM(abs(JST(J)))
       if(ia.lt.ja.and.ib.eq.jb)then
        call iswp(IST(J),IST(I))
        call iswp(JST(J),JST(I))
       ENDIF
105    continue

c      if same type, shorter bonds first
       DO 113 I=1,ISTRE
       do 113 J=I+1,ISTRE
       ia=abs(IST(I))
       ja=abs(IST(J))
       ib=abs(JST(I))
       jb=abs(JST(J))
       d1=(X(ia)-X(ib))**2+(Y(ia)-Y(ib))**2+(Z(ia)-Z(ib))**2
       d2=(X(ja)-X(jb))**2+(Y(ja)-Y(jb))**2+(Z(ja)-Z(jb))**2
       ibk=KATOM(ia)
       jbk=KATOM(ja)
       iak=KATOM(ib)
       jak=KATOM(jb)
       if(iak.eq.jak.and.ibk.eq.jbk.and.d1.gt.d2)then
        call iswp(IST(J),IST(I))
        call iswp(JST(J),JST(I))
       ENDIF
113    continue

       DO 115 I=1,ISTRE
       ia=abs(IST(I))
       ja=abs(JST(I))
       d=dsqrt((X(ia)-X(ja))**2+(Y(ia)-Y(ja))**2+(Z(ia)-Z(ja))**2)
       di(i)=d
       ia=KATOM(ia)
       ib=KATOM(ja)
c      C=O bond
       if(ia.eq.6.and.ib.eq.8.and.d.lt.1.300d0)stype(I)=' = '
c      aromatic C=C bond
       if(ia.eq.6.and.ib.eq.6.and.d.lt.1.454d0)stype(I)=' A '
c      C=C bond
       if(ia.eq.6.and.ib.eq.6.and.d.lt.1.358d0)stype(I)=' = '
c      tripple C\\\C bond
       if(ia.eq.6.and.ib.eq.6.and.d.lt.1.268d0)stype(I)=' T '
c      aromatic CN bond
       if(ia.eq.6.and.ib.eq.7.and.d.lt.1.397d0)stype(I)=' A '
c      C=N bond
       if(ia.eq.6.and.ib.eq.7.and.d.lt.1.309d0)stype(I)=' = '
c      triple C\\\N bond
       if(ia.eq.6.and.ib.eq.7.and.d.lt.1.214d0)stype(I)=' T '
c      N=N bond
       if(ia.eq.6.and.ib.eq.7.and.d.lt.1.335d0)stype(I)=' = '
c      triple N\\\N bond
       if(ia.eq.6.and.ib.eq.7.and.d.lt.1.176d0)stype(I)=' = '
115    continue

       DO 1061 I=1,IBEND
       ia=abs(IBD(I))
       ka=abs(KBD(I))
1061   if(KATOM(ka).lt.KATOM(ia))call iswp(KBD(I),IBD(I))

c      bend, first
       DO 106 I=1,IBEND
       do 106 J=I+1,IBEND
       ii=abs(IBD(I))
       ij=abs(IBD(J))
       if(KATOM(ij).lt.KATOM(ii))then
        call iswp(IBD(J),IBD(I))
        call iswp(JBD(J),JBD(I))
        call iswp(KBD(J),KBD(I))
       ENDIF
106    continue

c      bend, second
       DO 107 I=1,IBEND
       do 107 J=I+1,IBEND
       ii=abs(IBD(I))
       ij=abs(IBD(J))
       ji=abs(JBD(I))
       jj=abs(JBD(J))
       if(KATOM(jj).lt.KATOM(ji).and.KATOM(ii).eq.KATOM(ij))then
        call iswp(IBD(J),IBD(I))
        call iswp(JBD(J),JBD(I))
        call iswp(KBD(J),KBD(I))
       ENDIF
107    continue

c      bend, third
       DO 118 I=1,IBEND
       do 118 J=I+1,IBEND
       ii=abs(IBD(I))
       ij=abs(IBD(J))
       ji=abs(JBD(I))
       jj=abs(JBD(J))
       ki=abs(KBD(I))
       kj=abs(KBD(J))
       ia=KATOM(ii)
       ja=KATOM(ij)
       ib=KATOM(ji)
       jb=KATOM(jj)
       ic=KATOM(ki)
       jc=KATOM(kj)
       if(jc.lt.ic.and.ib.eq.jb.and.ia.eq.ja)then
        call iswp(IBD(J),IBD(I))
        call iswp(JBD(J),JBD(I))
        call iswp(KBD(J),KBD(I))
       ENDIF
118    continue

       DO 1181 I=1,ITOR
       ji=abs(JTOR(I))
       ki=abs(KTOR(I))
       if(KATOM(ki).lt.KATOM(ji))then
        call iswp(KTOR(I),JTOR(I))
        do 2081 ii=1,NITOR(I)
2081    TBUF(ii)=IBUF(I,ii)
        do 2082 ii=1,NLTOR(I)
2082    IBUF(I,ii)=JBUF(I,ii)
        do 2083 ii=1,NITOR(I)
2083    JBUF(I,ii)=TBUF(ii)
        call iswp(NITOR(I),NLTOR(I))
       endif
1181   continue

c      torsion, second
       DO 108 I=1,ITOR
       do 108 J=I+1,ITOR
       ji=abs(JTOR(I))
       jj=abs(JTOR(J))
       ia=KATOM(ji)
       ja=KATOM(jj)
       if(ja.lt.ia)then
        call iswp(JTOR(I),JTOR(J))
        call iswp(KTOR(I),KTOR(J))
        call iswp(NITOR(I),NITOR(J))
        call iswp(NLTOR(I),NLTOR(J))
        do 1081 ii=1,NITOR(J)
1081    TBUF(ii)=IBUF(I,ii)
        do 1082 ii=1,NITOR(I)
1082    IBUF(I,ii)=IBUF(J,ii)
        do 1083 ii=1,NITOR(J)
1083    IBUF(J,ii)=TBUF(ii)
        do 1084 ii=1,NLTOR(J)
1084    TBUF(ii)=JBUF(I,ii)
        do 1085 ii=1,NLTOR(I)
1085    JBUF(I,ii)=JBUF(J,ii)
        do 1086 ii=1,NLTOR(J)
1086    JBUF(J,ii)=TBUF(ii)
       ENDIF
108    continue

c      torsion, third
       DO 109 I=1,ITOR
       do 109 J=I+1,ITOR
       ji=abs(JTOR(I))
       jj=abs(JTOR(J))
       ki=abs(KTOR(I))
       kj=abs(KTOR(J))
       ia=KATOM(ji)
       ja=KATOM(jj)
       ib=KATOM(ki)
       jb=KATOM(kj)
       if(ia.eq.ja.and.jb.lt.ib)then
        call iswp(JTOR(I),JTOR(J))
        call iswp(KTOR(I),KTOR(J))
        call iswp(NITOR(I),NITOR(J))
        call iswp(NLTOR(I),NLTOR(J))
        do 1091 ii=1,NITOR(J)
1091    TBUF(ii)=IBUF(I,ii)
        do 1092 ii=1,NITOR(I)
1092    IBUF(I,ii)=IBUF(J,ii)
        do 1093 ii=1,NITOR(J)
1093    IBUF(J,ii)=TBUF(ii)
        do 1094 ii=1,NLTOR(J)
1094    TBUF(ii)=JBUF(I,ii)
        do 1095 ii=1,NLTOR(I)
1095    JBUF(I,ii)=JBUF(J,ii)
        do 1096 ii=1,NLTOR(J)
1096    JBUF(J,ii)=TBUF(ii)
       ENDIF
109    continue

c      oop first
       DO 110 I=1,IOP
       do 110 J=I+1,IOP
       ia=KATOM(abs(IOOP(I)))
       ja=KATOM(abs(IOOP(J)))
       if(ia.gt.ja)then
        call iswp(IOOP(I),IOOP(J))
        call iswp(JOOP(I),JOOP(J))
        call iswp(KOOP(I),KOOP(J))
        call iswp(LOOP(I),LOOP(J))
       ENDIF
110    continue
c
c      oop second
       DO 14 I=1,IOP
       do 14 J=I+1,IOP
       ia=KATOM(abs(IOOP(I)))
       ja=KATOM(abs(IOOP(J)))
       ib=KATOM(abs(JOOP(I)))
       jb=KATOM(abs(JOOP(J)))
       if(ia.eq.ja.and.jb.lt.ib)then
        call iswp(IOOP(I),IOOP(J))
        call iswp(JOOP(I),JOOP(J))
        call iswp(KOOP(I),KOOP(J))
        call iswp(LOOP(I),LOOP(J))
       ENDIF
14     continue

      endif

      IIT=0
      A1=1.0d0
      I0=0
      DO 101 I=1,ISTRE
      IIT=IIT+1
      ii=abs(IST(I))
      ji=abs(JST(I))
      iw=(ii-N*((ii-1)/N))
      jw=(ji-N*((ji-1)/N))
101   WRITE(7,444)1,IST(I),JST(I),0,0,0,0,1.0d0,
     1ATOMT(ii),iw,stype(I),ATOMT(ji),           jw, IIT, di(I)
444   FORMAT(7I6,F12.6,' STRETCH  ',A2,I4,A3,A2, I4, I6,  F8.4)
      DO 102 I=1,IBEND
      IIT=IIT+1
      ii=abs(IBD(I))
      ji=abs(JBD(I))
      ki=abs(KBD(I))
      iw=(ii-N*((ii-1)/N))
      jw=(ji-N*((ji-1)/N))
      kw=(ki-N*((ki-1)/N))
102   WRITE(7,555)2,IBD(I),JBD(I),KBD(I),I0,I0,I0,A1,
     1ATOMT(ii),iw,ATOMT(ji),jw,ATOMT(ki),kw,IIT
555   FORMAT(7I6,F12.6,' BEND   ',A2,I4,
     1' - ',A2,I4,' - ',A2,I4,I6)
      DO 103 I=1,ITOR
      IIT=IIT+1
      ji=abs(JTOR(I))
      ki=abs(KTOR(I))
      jw=(ji-N*((ji-1)/N))
      kw=(ki-N*((ki-1)/N))
      WRITE(7,666)4,NITOR(I),JTOR(I),KTOR(I),NLTOR(I),
     1I0,I0,A1,
     2ATOMT(ji),jw,ATOMT(ki),kw,IIT
      WRITE(7,666)(IBUF(I,J),J=1,NITOR(I))
103   WRITE(7,666)(JBUF(I,J),J=1,NLTOR(I))
666   FORMAT(7I6,F12.6,' TORSION   ',A2,I4,'  -  ',A2,I4,I6)
      DO 112 I=1,IOP
      IIT=IIT+1
      ii=abs(IOOP(I))
      ji=abs(JOOP(I))
      ki=abs(KOOP(I))
      li=abs(LOOP(I))
      iw=(ii-N*((ii-1)/N))
      jw=(ji-N*((ji-1)/N))
      kw=(ki-N*((ki-1)/N))
      lw=(li-N*((li-1)/N))
112   WRITE(7,333)3,IOOP(I),JOOP(I),KOOP(I),LOOP(I),I0,I0,A1,
     1ATOMT(ii),iw,ATOMT(ji),jw,ATOMT(ki),kw,ATOMT(li),lw,IIT
333   FORMAT(7I6,F12.6,' OPLA   ',A2,I4,
     1' - ',A2,I4,' - ',A2,I4,' - ',A2,I4,I6)
c
c     monomolecular:
      if(l(6).and.NMOL.gt.0)then
c      translations,rotations: 6 types 24 25 26  27 28 29
       do 116 ic=1,6
       do 116 i=1,NMOL
       IIT=IIT+1
       WRITE(7,7773)23+ic,nlist(i),i,0,0,0,0,1.0d0,i,mm(ic),IIT
7773   FORMAT(7I6,F12.6,' MOLECULE ',I4,' ',A2,', coord',I6)
       do 117 ia=1,nlist(i)
117    write(7,7771)mlist(i,ia)
116    write(7,*)
      endif

c     inter-molecular:
      if(IMM.ne.0)then
c
c      stretches a1 a2 b1 b2 torsion: 6 types 8 9 10 11 12 13
       do 114 ic=1,6
       do 114 i=1,NMOL
       do 114 ii=1,ilist(i)
       j=blist(i,ii)
       IIT=IIT+1
       WRITE(7,7777)7+ic,nlist(i),nlist(j),i,j,0,0,1.0d0,i,j,IIT
7777   FORMAT(7I6,F12.6,' INTERMOL.',I4,'  -',I4,', coord',I6)
       do 1141 ia=1,nlist(i)
1141   write(7,7771)mlist(i,ia)
7771   format(i6,$)
       write(7,*)
       do 1142 ia=1,nlist(j)
1142   write(7,7771)mlist(j,ia)
114    write(7,*)
      endif

c     helical/cigar
c     10 types 14 15 16 17 18  19 20 21 22 23
      do 1143 ic=1,IHH
      IIT=IIT+1
      WRITE(7,7778)13+ic,nha,nha,HC(ic),IIT
7778  FORMAT(3I6,4(5x,'0'),4x,'1.',6('0'),A10,3x,'0 ',3x,'0, coord',I6)
c     write the atom list only for the first coordinate:
1143  if(ic.eq.1)write(7,7772)(hl(ia),ia=1,nha)
7772  format(12i6)

      WRITE(7,*)' 0'
      A=1.0d0
      DO 771 I=1,3*N-6
771   WRITE(7,777)I,A
777   FORMAT(' 1',/,I5,F5.1)
1002  CLOSE(7)
      WRITE(*,*)' PROGRAM TERMINATED'
      END

      SUBROUTINE INCO(IERR,N,ns0,IBOND,ISTRE,IBEND,ITOR,IST,JST,IBD,
     1JBD,KBD,NITOR,NLTOR,JTOR,KTOR,IBUF,JBUF,IE,IOP)
      implicit none
      integer*4 IERR,N,ns0,IBOND(N+2,8),ISTRE,IBEND,ITOR,IST(*),JST(*),
     1IBD(*),JBD(*),KBD(*),NITOR(*),NLTOR(*),JTOR(*),KTOR(*),
     3IBUF(2*ns0,6),JBUF(2*ns0,6),IE(*),IOP,I0,IT,IR,
     1IB0,I,J,ICO,L2,L,K,K3,IA,KK,IB,icontr,IC,iw
      integer*4,allocatable::IAA(:,:),ICON(:,:),IDT(:),IDB(:),IBACK(:)
      allocate (IAA(8,N),ICON(N+2,8),IDT(N),IDB(N),IBACK(N+2))
      IOP=0
      IERR=0
      DO 1 I=1,N
      DO 1 J=1,8
1     ICON(I,J)=IBOND(I,J)
      DO 2 I=N+1,N+2
      DO 2 J=1,8
      IBOND(I,J)=0
2     ICON(I,J)=0
      IBACK=0
      IBACK(1)=N+1
24    ICO=0
      DO 1111 L2=1,N
      IF (IBACK(L2).EQ.0)GOTO 1111
      DO 2222 K3=1,4
      L=L2
      K=K3
      IA=IBOND(L,K)
      IF (IA.NE.0)THEN
       ICO=1
       IBOND(L,K)=0
5555   CONTINUE
       IF (IBACK(IA).EQ.0)IBACK(IA)=L
       DO 33 KK=1,8
       IF (IBOND(IA,KK).EQ.L)IBOND(IA,KK)=0
33     CONTINUE
       DO 44 KK=1,8
       IF (IBOND(IA,KK).NE.0)THEN
        L=IA
        IA=IBOND(IA,KK)
        IBOND(L,KK)=0
        GOTO 5555
       ENDIF
44     CONTINUE
      ENDIF
2222  CONTINUE
1111  CONTINUE
      IF (ICO.EQ.1)GOTO 24
      IBACK(N+1)=N+2
      IB=0
      DO 5 I=2,N
      IB=IB+1
      IAA(4,IB)=I
      IAA(3,IB)=IBACK(max(IAA(4,IB),1))
      IAA(2,IB)=IBACK(max(IAA(3,IB),1))
5     IAA(1,IB)=IBACK(max(IAA(2,IB),1))
c     if(nat.gt.100)then
c      write(*,8000)nat
c8000  format(I5,' atoms seems a lot -',/,
c     1 ' Do you still want the internal coordinates (0/1) ?')
c       read(5,*)iok
c       if(iok.ne.1)then
c        ierr=1
c        return
c       endif
c      endif
      WRITE(*,999)
999   FORMAT(1X,' BOND 1 - 2 - 3 -> (4) (IB) WAS DEFINED')
C
C     IF AN ATOM IA1 OF THE BOND IB IS NOT YET DEFINED, THE BOND IS PU
C     TO THE BACK OF THE LIST OF BONDS:
      icontr=0
47    IC=0
      IB0=IB
      DO 45 I=2,IB-1
      DO 46 J=1,I-1
46    IF (IAA(3,I).EQ.IAA(4,J).OR.IAA(3,I).EQ.1)GOTO 45
      IC=1
      DO 51 J=1,4
      I0=IAA(J,I)
      IAA(J,I)=IAA(J,I+1)
51    IAA(J,I+1)=I0
45    CONTINUE
      icontr=icontr+1
      if (icontr.gt.5000) then
       WRITE(*,*)'INFINITE LOOP'
       IERR=1
       RETURN
      endif
      IF (IC.EQ.1)GOTO 47
      WRITE(*,*)' BONDS ARRANGED'
      ISTRE=0
      IBEND=0
      ITOR=0
      DO 3 I=1,N
      IDB(I)=0
3     IDT(I)=0
      DO 54 IB=1,IB0
         ISTRE=ISTRE+1
         IST(ISTRE)=IAA(4,IB)
         JST(ISTRE)=IAA(3,IB)
          IBEND=IBEND+1
          IBD(IBEND)=IAA(4,IB)
          JBD(IBEND)=IAA(3,IB)
          KBD(IBEND)=IAA(2,IB)
          IDB(IBD(IBEND))=KBD(IBEND)
          IF (KBD(IBEND).GT.N) THEN
           IW=0
           DO 55 I=1,8
           IR=ICON(JBD(IBEND),I)
           IF (IR.EQ.IBD(IBEND))GOTO 55
           IF (IR.NE.0) IW=IR
55         CONTINUE
           IF (IW.EQ.0)IBEND=IBEND-1
           IF (IW.NE.0)THEN
            KBD(IBEND)=IW
            IDB(IBD(IBEND))=IW
            IF (IDB(IW).EQ.IBD(IBEND))IBEND=IBEND-1
           ENDIF
          ENDIF
          ITOR=ITOR+1
          if(ITOR.gt.2*ns0)call report('too many torsions')
          NITOR(ITOR)=1
          NLTOR(ITOR)=1
          JTOR(ITOR)=IAA(2,IB)
          KTOR(ITOR)=IAA(3,IB)
          IBUF(ITOR,1)=IAA(1,IB)
          JBUF(ITOR,1)=IAA(4,IB)
          IDT(JBUF(ITOR,1))=IBUF(ITOR,1)
          IF (IBUF(ITOR,1).EQ.N+1)THEN
           IW=0
           DO 56 I=1,8
           IR=ICON(JTOR(ITOR),I)
           IF (IR.EQ.KTOR(ITOR))GOTO 56
           IF (IR.NE.0) IW=IR
56         CONTINUE
           IF (IW.EQ.0)ITOR=ITOR-1
           IF (IW.NE.0)THEN
            IDT(JBUF(ITOR,1))=IW
            IBUF(ITOR,1)=IW
           ENDIF
           if(KTOR(ITOR+1).lt.1)goto 111
           IF ((IW.EQ.0).AND.(IE(KTOR(ITOR+1)).GT.2)) THEN
            DO 59 I=1,IE(KTOR(ITOR+1))
            IT=ICON(KTOR(ITOR+1),I)
            IF (IT.EQ.JTOR(ITOR+1).OR.IT.EQ.JBUF(ITOR+1,1))GOTO 59
            IF (IDB(IT).NE.JBUF(ITOR+1,1))THEN
             IBEND=IBEND+1
             IBD(IBEND)=IT
             JBD(IBEND)=KTOR(ITOR+1)
             KBD(IBEND)=JBUF(ITOR+1,1)
             IDB(JBUF(ITOR+1,1))=IT
            ENDIF
            GOTO 60
59          CONTINUE
           ENDIF
111        continue
60         CONTINUE
          ENDIF
          IF (IBUF(ITOR,1).EQ.N+2)THEN
           K=KTOR(ITOR)
           IW=0
           DO 57 I=1,8
           IR=ICON(KTOR(ITOR),I)
           IF (IR.EQ.JBUF(ITOR,1))GOTO 57
           IF (IR.NE.0) IW=IR
57         CONTINUE
           IF (IW.EQ.0)ITOR=ITOR-1
           IF (IW.EQ.0)GOTO 54
           IC=0
           DO 58 I=1,8
           IR=ICON(JBUF(ITOR,1),I)
           IF (IR.EQ.KTOR(ITOR)) GOTO 58
           IF (IR.NE.0)IC=IR
58         CONTINUE
           IF ((IDT(IW).LE.N).AND.(IDT(IC).GE.1))IC=0
           IF (IC.EQ.0)ITOR=ITOR-1
           IF (IC.EQ.0)GOTO 54
           JTOR(ITOR)=KTOR(ITOR)
           KTOR(ITOR)=JBUF(ITOR,1)
           JBUF(ITOR,1)=IC
           IBUF(ITOR,1)=IW
           IDT(IW)=IC
          ENDIF
54    CONTINUE
      RETURN
      END
c     ===========================================================     
      subroutine iswp(i,j)
      integer*4 ii,i,j
      ii=i
      i=j
      j=ii
      return
      end
c     ===========================================================     
      function sp(a,b)
      IMPLICIT none
      real*8 a(*),b(*),sp
      sp=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
c     ===========================================================     
      subroutine norm(v)
      real*8 v(3),n,sp
      n=dsqrt(sp(v,v))
      v(1)=v(1)/n
      v(2)=v(2)/n
      v(3)=v(3)/n
      return
      end
c     ===========================================================     
      subroutine vp(A,X,Y)
      implicit none
      real*8 A(*),X(*),Y(*)
      A(1)=X(2)*Y(3)-X(3)*Y(2)
      A(2)=X(3)*Y(1)-X(1)*Y(3)
      A(3)=X(1)*Y(2)-X(2)*Y(1)
      return
      end
c     ===========================================================     
      subroutine MAKEBONDS(nat,X,Y,Z,iz,nb,nbt,N7)
      IMPLICIT none
      real*8 bonding(88),X(*),Y(*),Z(*),rb,x1,y1,z1,rt,x2,y2,z2
      integer*4 i,nat,j,N7,nbt(nat,N7),iz(*),nb(*),
     1ti,tj
      real*4 bondin4(88)
c       Du    H   He   Li   Be    B    C    N    O    F
c       Ne   Na   Mg   Al   Si    P    S   Cl   Ar    K   Ca   Sc   Ti
c        V   Cr   Mn   Fe   Co   Ni   Cu   Zn   Ga   Ge   As 
      data bondin4/
     10.50,0.32,0.98,1.28,0.95,0.87,0.82,0.80,0.78,0.77,
     10.76,1.59,1.41,1.23,1.16,1.11,1.20,1.04,1.03,2.08,1.79,1.49,1.37,
     11.27,1.23,1.22,1.22,1.21,1.20,1.22,1.30,1.31,1.27,1.25,1.21,1.19,
     11.17,2.21,1.96,1.67,1.50,1.39,1.35,1.32,1.30,1.30,1.33,1.39,1.53,
     11.49,1.46,1.45,1.41,1.38,1.36,2.40,2.03,1.74,1.70,1.70,1.69,1.68,
     11.67,1.90,1.66,1.64,1.64,1.63,1.62,1.61,1.75,1.61,1.49,1.39,1.35,
     11.33,1.31,1.32,1.35,1.39,1.54,1.53,1.52,1.51,1.51,1.50,1.49,0.10/
     
      do 1 i=1,88
1     bonding(i)=dble(bondin4(i))
      do 4 i=1,nat
      nb(i)=0
      do 4 j=1,N7
4     nbt(i,j)=0
c
      do 3 i=1,nat
      ti=iz(i)+1
      if(ti.lt.1.or.ti.gt.89)call report('type out of range')
      rb=bonding(ti)
      x1=X(i)
      y1=Y(i)
      z1=Z(i)
      do 3 j=i+1,nat
      tj=iz(j)+1
      if(ti.ne.2.or.tj.ne.2)then
       if(tj.lt.1.or.tj.gt.88)call report('type out of range')
       rt=(rb+bonding(tj))**2
       x2=(x1-X(j))**2
       if(x2.lt.rt)then 
        y2=(y1-Y(j))**2+x2
        if(y2.lt.rt)then
         z2=(z1-Z(j))**2+y2
         if(z2.lt.rt)then
          nb(i)=nb(i)+1
          nb(j)=nb(j)+1
          nbt(i,nb(i))=j
          nbt(j,nb(j))=i
         endif
        endif
       endif
      endif
3     continue
c
      return
      end
c     ===========================================================     
      subroutine bmol(NMOL,X,Y,Z,nat,nlist,mlist,blist,ilist,iz,mi,rl0)
      implicit none
      integer*4 NMOL,nat,nlist(*),mlist(nat,nat),ia,ib,ii,
     1ilist(*),blist(NMOL,NMOL),im,mi,ip,jj,iz(*)
      real*8 x(*),Y(*),Z(*),rlim,rab,r1,wdr(118),rl0
c     van der waals radii in pm, unknown: 169 pm:
      integer*4 vdw(118)
      data vdw/120,140,182,169,169,170,155,152,147,154,
     1         227,173,169,210,180,180,175,188,275,169,
     2         169,169,169,169,169,169,169,163,140,139,
     3         187,169,185,190,185,202,169,169,169,169,
     4         169,169,169,169,169,163,172,158,193,217,
     5         169,206,198,216,169,169,169,169,169,169,
     6         169,169,169,169,169,169,169,169,169,169,
     7         169,169,169,169,169,169,169,175,166,155,
     8         196,202,169,169,169,169,169,169,169,169,
     9         169,186,169,169,169,169,169,169,169,169,
     1         169,169,169,169,169,169,169,169,169,169,
     2         169,169,169,169,169,169,169,169/
      
      do 10 ia=1,118
10    wdr(ia)=dble(vdw(ia))/100.0d0

      do 1 im=1,NMOL 
1     ilist(im)=0
      mi=0

      do 3 im=1,NMOL 
      do 3 ip=im+1,NMOL 

      do 4 ii=1,nlist(im)
      ia=mlist(im,ii)
      r1=wdr(iz(ia))
      do 4 jj=1,nlist(ip)
      ib=mlist(ip,jj)
      if(rl0.lt.1.0d-3)then
       rlim=r1+wdr(iz(ib))+0.4d0
      else
       rlim=rl0
      endif
      rab=dsqrt((X(ia)-X(ib))**2+(Y(ia)-Y(ib))**2+(Z(ia)-Z(ib))**2)
      if(rab.lt.rlim)then
       ilist(im)=ilist(im)+1
       blist(im,ilist(im))=ip
       mi=mi+1
       goto 3
      endif
4     continue

3     continue

      write(6,6003)mi
6003  format(i4,' molecular interactions found:')
      do 2 im=1,NMOL
      if(ilist(im).gt.0)then
       write(6,6001)im
6001   format('Molecule',i4,' with:')
       write(6,6002)(blist(im,ii),ii=1,ilist(im))
6002   format(8x,20i4)
      endif
2     continue
      return
      end
c     ===========================================================     
      subroutine gmol(N,NMOL,X,Y,Z,iz,nlist,mlist,l,BT,N7)
c     divide a system to individual molecules      
      implicit none
      integer*4 N,NMOL,iz(*),nlist(*),mlist(N,N),ia,ib,ii,im,
     1nleft,N7,BT(N,N7)
      real*8 x(*),Y(*),Z(*)
      integer*4,allocatable::nb(:),nbt(:,:),ind(:)
      logical l(*)
      NMOL=0
      if(N.eq.0)return
      allocate(nb(N),nbt(N,N7))
      if(l(4))then
c      take given bond table:
       do 8 ia=1,N
       nb(ia)=0
       do 81 ii=1,N7
       if(BT(ia,ii).ne.0)nb(ia)=nb(ia)+1
81     nbt(ia,ii)=BT(ia,ii)
8      continue
      else
c      make bondtable:
       call MAKEBONDS(N,X,Y,Z,iz,nb,nbt,N7)
      endif
c     find molecules:
      allocate(ind(N))
      do 1 ia=1,N
1     ind(ia)=0
      NMOL=1
      nlist(NMOL)=1
      mlist(NMOL,1)=1
      ind(1)=NMOL
      nleft=N-1
c     find all atoms in the current molecule:
99    do 3 ia=1,N

      do 4 im=1,nlist(NMOL)
c     atom ia is already in:
4     if(mlist(NMOL,im).eq.ia)goto 3

      do 5 ii=1,nb(ia)
      ib=nbt(ia,ii)
      do 5 im=1,nlist(NMOL)
c     atom ia is not in, but is bound to ib which is already in:
      if(mlist(NMOL,im).eq.ib)then
c      include atom ia in molecule NMOL, and start over:
       nlist(NMOL)=nlist(NMOL)+1
       mlist(NMOL,nlist(NMOL))=ia
       ind(ia)=NMOL
       nleft=nleft-1
       goto 99
      endif
5     continue

3     continue

      if(nleft.gt.0)then
c      some unassigned atoms left, start a new molecule
       do 6 ia=1,N
       if(ind(ia).eq.0)then
        NMOL=NMOL+1
        nlist(NMOL)=1
        mlist(NMOL,1)=ia
        nleft=nleft-1
        ind(ia)=NMOL
        goto 99
       endif
6      continue
      endif

      write(6,*)NMOL,' molecules'
      do 7 ii=1,NMOL
      write(6,6001)ii
6001  format(i4,':')
7     write(6,6002)(mlist(ii,ia),ia=1,nlist(ii))
6002  format(20i4)


      return
      end
c     ===========================================================     
      subroutine report(s)
      character*(*) s
      write(6,*)s
      stop
      end
c     ===========================================================     
      subroutine  hmol1(nha,N)
      implicit none
      integer*4 nha,N
      logical lex
c     determine number of helix atoms:
      inquire(file='HA.LST',exist=lex)
      if(lex)then
       open(9,file='HA.LST')
       read(9,*)nha
       close(9)
      else
       nha=N
      endif
      write(6,*)nha,' helix atoms'
      if(nha.gt.N.or.nha.lt.1)call report('Invalid number')
      return
      end
c     ===========================================================     
      subroutine  hmol2(hl,nha,N)
c     determine helix atoms:
      implicit none
      integer*4 nha,N,i,hl(*)
      logical lex
      inquire(file='HA.LST',exist=lex)
      if(lex)then
       open(9,file='HA.LST')
       read(9,*)hl(1),(hl(i),i=1,nha)
       close(9)
      else
       do 1 i=1,N
1      hl(i)=i
      endif
      return
      end
c     ===========================================================     
      subroutine overview(ISTRE,IBEND,ITOR,IOP,M1,IMM,IHH,NTOT,N)
      integer*4 ISTRE,IBEND,ITOR,IOP,IMM,IHH,NTOT,N,M1
      write(*,*)
      if(ISTRE.ne.0)WRITE(*,*) ISTRE,' OF STRETCHES'
      if(IBEND.ne.0)WRITE(*,*) IBEND,' OF VALENCE BENDS'
      if(ITOR .ne.0)WRITE(*,*)  ITOR,' OF TORSIONS'
      if(IOP  .ne.0)WRITE(*,*)   IOP,' OUT OF PLANE WAGGING'
      if( M1  .ne.0)WRITE(*,*)    M1,' MOLECULAR'
      if(IMM  .ne.0)WRITE(*,*)   IMM,' INTERMOLECULAR'
      if(IHH  .ne.0)WRITE(*,*)   IHH,' HELIX DEF'
      WRITE(*,*)' ------------------------------------'
      WRITE(*,*)NTOT,' TOTAL;  3N-6 = ',3*N-6
      return
      end
c     ===========================================================     
      subroutine addb(N,nb,NBT,I,J,BT,IBOND,ICON,N7)
      integer*4 N,N7,nb,NBT(*),I,J,BT(N,N7),IBOND(*),ICON(N+2,8)
      NBT(I)=NBT(I)+1
      NBT(J)=NBT(J)+1
      BT(I,NBT(I))=J
      BT(J,NBT(J))=I
      nb=nb+1
      IBOND(I)=IBOND(I)+1
      IBOND(J)=IBOND(J)+1
      ICON(I,IBOND(I))=J
      ICON(J,IBOND(J))=I
      return
      end
c     ===========================================================     
      function acidic(kk)
      integer*4 kk
      logical acidic
      acidic=kk.eq.7.or.kk.eq.8.or.kk.eq.16.or.kk.eq.9.or.kk.eq.17.
     1or.kk.eq.35.or.kk.eq.53
      return
      end
c     ===========================================================     
      subroutine genb(N,X,Y,Z,KATOM,RAD0,RAD,NBT,BT,nb,IBOND,ICON,lh,
     1N7)
      implicit none
      integer*4 N7,KATOM(*),N,nb,NBT(*),I,J,K,BT(N,N7),IBOND(*),
     1ICON(N+2,8),ki,kj,kk,nh
      real*8 X(*),Y(*),Z(*),RAD0,RAD(*),R9,rh1,rh2,xi,yi,zi,rik,
     1xh,yh,zh
      logical lh,acidic,ac1
      rh1=1.6d0**2
      rh2=2.2d0**2
      nh=0
      DO 1 I=1,N
      ki=KATOM(I)
      xi=X(I)
      yi=Y(I)
      zi=Z(I)
      DO 1 J=I+1,N
      kj=KATOM(J)
      R9=(X(J)-xi)**2+(Y(J)-yi)**2+(Z(J)-zi)**2
      IF (R9.LE.( RAD(ki)+ RAD(kj) )**2) then
       call addb(N,nb,NBT,I,J,BT,IBOND,ICON,N7)
      else
c      Hydrogen bonds:
       if(lh.and.(ki.eq.1.or.kj.eq.1))then
        if(ki.eq.1)then
         kk=kj
         xh=xi
         yh=yi
         zh=zi
        else
         kk=ki
         xh=X(J)
         yh=Y(J)
         zh=Z(J)
        endif
        ac1=.false.
        do 2 K=1,N
        if(K.ne.I.and.K.ne.J)then
         rik=(X(K)-xh)**2+(Y(K)-yh)**2+(Z(K)-zh)**2
         if(rik.lt.1.44d0.and.acidic(KATOM(K)))then
          ac1=.true.
          goto 99
         endif
        endif
 2      continue

99      if(ac1.and.acidic(kk))then
          if(R9.gt.rh1.and.R9.lt.rh2)then
           call addb(N,nb,NBT,I,J,BT,IBOND,ICON,N7)
           nh=nh+1
          endif
        endif
       endif     
      ENDIF
1     CONTINUE
      if(lh.and.nh.gt.0)write(6,*)nh,' H-bonds'
      return
      end
c     ===========================================================     
      subroutine genstr(ISTRE,ns0,NBT,N,Nto,BT,IST,JST,N7)
      implicit none
      integer*4 N,Nto,ISTRE,ns0,NBT(*),N7,BT(N,N7),IST(*),JST(*),I,JJ,J
      ISTRE=0
      DO 13 I=1,Nto
      DO 13 JJ=1,NBT(I)
      J=BT(I,JJ)
      if(J.gt.I)then
       ISTRE=ISTRE+1
       if(ISTRE.gt.ns0)call report('too many stretches')
       if(J.gt.Nto)J=-J
       IST(ISTRE)=I
       JST(ISTRE)=J
      ENDIF
13    CONTINUE
      return
      end
c     ===========================================================     
      subroutine genbe(IBEND,ns0,NBT,N,Nto,BT,IBD,JBD,KBD,N7)
      implicit none
      integer*4 IBEND,ns0,NBT(*),N,N7,BT(N,N7),IBD(*),JBD(*),I,JJ,J,
     1KBD(*),Nto,ix,KK,K
      IBEND=0
      DO 133 I=1,Nto
      DO 133 JJ=1,NBT(I)
      J=BT(I,JJ)
      DO 133 KK=1,NBT(J)
      K=BT(J,KK)
      if(K.gt.I)then
       IBEND=IBEND+1
       if(IBEND.gt.ns0)then
        do 1 ix=1,IBEND-1
1       write(6,600)ix,IBD(ix),JBD(ix),KBD(ix)
600     format(4i6)
        call report('too many bends')
       endif
       IBD(IBEND)=I
       if(J.gt.Nto)then
        JBD(IBEND)=-J
       else
        JBD(IBEND)=J
       endif
       if(K.gt.Nto)K=-K
       KBD(IBEND)=K
      ENDIF
133   CONTINUE
      return
      end
c     ===========================================================     
      subroutine gento(ITOR,ns0,NBT,N,Nto,BT,NITOR,IBUF,JTOR,KTOR,
     1NLTOR,JBUF,lt,N7)
      implicit none
      integer*4 ITOR,N,ns0,NBT(*),N7,BT(N,N7),I,JJ,J,KK,K,LL,L,NITOR(*),
     1IBUF(2*ns0,6),JTOR(*),KTOR(*),NLTOR(*),
     1JBUF(2*ns0,6),Nto,ih,jh,
     1II,ni,nj
      logical lt
      ITOR=0
      if(lt)then
c      full torsion (bond rotations)
c            I --- J
c         ni         nj   ... numbers of other hanging atoms on i,j
       DO  1 I=1,Nto
       DO  1 KK=1,NBT(I)
       J=BT(I,KK)
       if(J.gt.I)then
        ni=0
        nj=0
        DO 2 II=1,NBT(I)
2       if(BT(I,II).ne.J)ni=ni+1
        DO 3 JJ=1,NBT(J)
3       if(BT(J,JJ).ne.I)nj=nj+1
        if(ni.gt.0.and.nj.gt.0)then
         ITOR=ITOR+1
         if(ITOR.gt.2*ns0)call report('too many torsions')
         NITOR(ITOR)=ni
         NLTOR(ITOR)=nj
         ni=0
         DO 4 II=1,NBT(I)
         ih=BT(I,II)
         if(ih.ne.J)then
          ni=ni+1
          IBUF(ITOR,ni)=ih
          if(ih.gt.Nto)IBUF(ITOR,ni)= -ih
         endif
4        continue
         nj=0
         DO 5 JJ=1,NBT(J)
         jh=BT(J,JJ)
         if(jh.ne.I)then
          nj=nj+1
          JBUF(ITOR,nj)=jh
          if(jh.gt.Nto)JBUF(ITOR,nj)= -jh
         endif
5        continue
         JTOR(ITOR)=I
         KTOR(ITOR)=J
         if(J.gt.Nto)KTOR(ITOR) = -J
        endif
       endif
1      continue
      else
c      I --- J --- K --- L
       DO 135 I=1,Nto
       DO 135 JJ=1,NBT(I)
       J=BT(I,JJ)
       DO 135 KK=1,NBT(J)
       K=BT(J,KK)
       if(K.ne.I)then
        DO 137 LL=1,NBT(K)
        L=BT(K,LL)
        if(L.gt.I.and.L.ne.J)then
         ITOR=ITOR+1
         if(ITOR.gt.2*ns0)call report('too many torsions')
         NITOR(ITOR)=1
         IBUF(ITOR,1)=I
         if(J.gt.Nto)then
          JTOR(ITOR)=-J
         else
          JTOR(ITOR)=J
         endif
         if(K.gt.Nto)then
          KTOR(ITOR)=-K
         else
          KTOR(ITOR)=K
         endif
         NLTOR(ITOR)=1
         if(L.gt.Nto)L=-L
         JBUF(ITOR,1)=L
        endif
137     continue
       endif
135    CONTINUE
      endif
      return
      end
c     ===========================================================     
      subroutine giop(IOP,ns0,NBT,N,Nto,BT,IOOP,JOOP,KOOP,LOOP,
     1X,Y,Z,N7)
      implicit none
      integer*4 IOP,ns0,NBT(*),N,N7,BT(N,N7),I,J,KK,K,LL,L,Nto,
     1IOOP(ns0),JOOP(ns0),KOOP(ns0),LOOP(ns0)
      real*8 X(*),Y(*),Z(*),sp,KL(3),KI(3),KJ(3),av(3),oopl
      oopl=0.1d0
      IOP=0
      DO 138 I=1,Nto
      if(NBT(I).EQ.1)then
       J=BT(I,1)
       DO 139 KK=1,NBT(J)
       K=BT(J,KK)
       if(K.ne.I)then
        DO 140 LL=1,NBT(J)
        L=BT(J,LL)
        if(L.ne.I.and.L.gt.K)then
c                K
c               /
c         I - J           , I is the wagging atom
c               \
c                L
         KL(1)=X(K)-X(L)
         KL(2)=Y(K)-Y(L)
         KL(3)=Z(K)-Z(L)
         KI(1)=X(K)-X(I)
         KI(2)=Y(K)-Y(I)
         KI(3)=Z(K)-Z(I)
         KJ(1)=X(K)-X(J)
         KJ(2)=Y(K)-Y(J)
         KJ(3)=Z(K)-Z(J)
c        vector perpendicular to plane JKL:
         call vp(av,KL,KJ)
         call norm(av)
c        deviation of I from the plane < oopl:
         if(dabs(sp(av,KI)).lt.oopl)then
          IOP=IOP+1
          if(IOP.gt.ns0)call report('too many oops')
          IOOP(IOP)=I
          if(J.gt.Nto)J=-J
          JOOP(IOP)=J
          if(K.gt.Nto)K=-K
          KOOP(IOP)=K
          if(L.gt.Nto)L=-L
          LOOP(IOP)=L
         endif
        endif
140     continue
       endif
139    continue   
      endif
138   CONTINUE
      return
      end
