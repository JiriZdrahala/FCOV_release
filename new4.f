      PROGRAM NEW4
C     Harmonic force field diagonalization,
c     works in cartesian coordinates, for one molecule only
c     TRED12 diagonalization replaced by JACOBI_EIGENVALUE
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*1 OK
      CHARACTER*4 JOB
      LOGICAL LBMAT
      integer*4,allocatable::NC(:,:),NM(:)
      real*8,allocatable::FCAR(:,:),FINT(:,:),FNEW(:,:),SMAT(:,:),
     1A(:,:),AMULT(:),SCFAC(:),ZM(:),Q(:),ZNU(:),AS(:,:),ZIM(:),
     2CHR(:)
C
      wtol=0.01d0
      c1=1302.828d0
      it0=mclock()

      nlinepar=iargc()
      ili=0

      OPEN(16,FILE='FTRY.OUT')
      WRITE(*,60000)
      WRITE(16,60000)
60000 FORMAT('                    FORCE FIELD REFINEMENT PROGRAM',/,
     1       '                PC version - in cartesian coordinates',/,
     2       '                         Scaling of masses!',/,/,
     3       '                       Petr Bour, Prague 2-95',/,/)
C
      OPEN(17,FILE='FILE.X')
      READ(17,*)
      READ(17,*)NOAT
      CLOSE(17)
      NA=3*NOAT
      allocate(FCAR(NA,NA),FINT(NA,NA),FNEW(NA,NA),SMAT(NA,NA),
     1A(NA,NA),AMULT(NA),SCFAC(NOAT),ZM(NA),Q(NA),ZNU(NA),
     1NC(NOAT,NOAT),NM(NOAT),AS(3,NA),ZIM(NOAT),CHR(NOAT))
      AMULT=0.0D0
      Q=0.0d0
      FINT=0.0d0
      FNEW=0.0d0
      SMAT=0.0d0
      AS=0.0d0
C
C     Read-in cartesian FF from FILE.FC:
C     FCAR is converted into au/au^2
      CALL READFF(NA,FCAR)
C
C     Read-in scalefactors etc from FTRY.INP:
      call readinp(JOB,ISYM,NOSC,SCFAC,NM,NC,CHR,NOAT,NA,NNU,ZIM,ZM,ZNU)

      write(6,6009)mclock()-it0
6009  format(' Time: ',i10,' msec')
C
      WRITE(*,*)'Project out zero vibrations (y/n)?'
      OK='Y'
      if(nlinepar.gt.ili)then
       ili=ili+1
       call getarg(ili,OK)
       write(6,*)OK
      else
       IF(JOB.NE.'AUTO')READ(*,'(A)')OK
      endif

      IF(OK.NE.'n'.and.OK.NE.'N')CALL PROJECT(A,FCAR,NA)
      write(6,6009)mclock()-it0
C
C     Generate the scale constants at the start of iteration loop
      ITRY=0
111   ITRY=ITRY+1
      WRITE(16,*)' TRIAL ',ITRY
      WRITE(16,*)NOSC, ' THE SCALE FACTORS'
      WRITE(16,*)' SF       ATOMS'
      DO 200 I=1,NOSC
200   WRITE(16,10400) SCFAC(I),(NC(I,J),J=1,NM(I))
10400 FORMAT(F9.6,100I3)

C     FCAR is now in (au/au^2)/(amu^1/2*amu^1/2)
      CALL FAST(NA,FCAR,FNEW,ZM)

      ie=0
      if(ISYM.eq.1)then
       call SYMROUTE(NOAT,NA,FNEW,SMAT,AS)
      else
       ADUM=SUMk(NOSC,1,ie,NA,FNEW,SMAT,ZNU,AMULT,Q,AS,SCFAC,NM,NC)
      endif

      write(6,6009)mclock()-it0

      if(ie.gt.0)then
       WRITE(6,6667)
6667   FORMAT('  NEW SCALEFACTORS (Y/N/E) ?')
       OK='N'
       if(nlinepar.gt.ili)then
        ili=ili+1
        call getarg(ili,OK)
        write(6,*)OK
       else
        IF(JOB.NE.'AUTO') READ(*,'(A)')OK
       endif

       IF(OK.EQ.'y'.OR.OK.EQ.'Y')THEN
        DO 601 I=1,NOSC
        WRITE(*,6668)I,SCFAC(I)
6668    FORMAT(1X,'NEW VALUE FOR THE ',I3,'. SCALE FACTOR (OLD ',
     1  f5.3,' ) :')
601     READ(*,*)SCFAC(I)
        GOTO 111
       ENDIF

       IF(OK.EQ.'E'.or.ok.eq.'e')then
        write(6,*)'Direct FF scaling to experiment, interpolate (Y/N)?'
        if(nlinepar.gt.ili)then
         ili=ili+1
         call getarg(ili,OK)
         write(6,*)OK
        else
         read(5,'(A)')OK
        endif
       if(OK.eq.'y'.or.OK.eq.'Y')then
         do 2 iq=1,NA
         if(AS(2,iq).lt.0.01d0)then
          i1=0
          do 3 is=iq-1,1,-1
          if(AS(2,is).gt.wtol)then
           i1=is
           goto 33
          endif 
3         continue
33        i2=0
          do 4 ie=iq+1,NA
          if(AS(2,ie).gt.wtol)then
           i2=ie
           goto 44
          endif 
4         continue
44        if(i1.gt.0.and.i2.gt.0)AS(2,iq)=AS(2,is)
     1    +(AS(1,iq)-AS(1,is))/(AS(1,ie)-AS(1,is))*(AS(2,ie)-AS(2,is))
          if(i1.gt.0.and.i2.eq.0)AS(2,iq)=AS(2,is)+AS(1,iq)-AS(1,is)
          if(i1.eq.0.and.i2.gt.0)AS(2,iq)=AS(2,ie)+AS(1,iq)-AS(1,ie)
          if(AS(1,iq).lt.wtol)AS(2,iq)=0.0d0
          write(6,6001)iq,AS(2,iq)
6001      format(i6,' set to ',f8.2)
         endif
2        continue
        endif
        do 6 i=1,NA
        do 6 j=1,NA
        FCAR(i,j)=0.0d0
        do 6 jp=1,NA
6       FCAR(i,j)=FCAR(i,j)
     1  +SMAT(i,jp)*(AS(2,jp)/c1)**2*SMAT(j,jp)/ZM(i)/ZM(j)
        OPEN(20,FILE='NEW.FC')
        CALL WRITEFF(NA,FCAR)
        CLOSE(20)
        write(6,*)'NEW.FC written'
        stop
       endif
       
       WRITE(*,*)' AUTOMATIC REFINEMENT (Y/N) ?'
       OK='N'
       if(nlinepar.gt.ili)then
        ili=ili+1
        call getarg(ili,OK)
        write(6,*)OK
       else
        IF(JOB.NE.'AUTO') READ(*,'(A)')OK
       endif
       IF (OK.EQ.'y'.OR.OK.EQ.'Y')THEN
        DEL=0.05d0
333     WRITE(*,*)' MAXIMUM NUMBER OF ITERATIONS:'
        READ(*,*)NIT
        IF (NIT.LT.0)THEN
         NIT=-NIT
         WRITE(*,*)' STEP:'
         READ(*,*)DEL
        ENDIF
        WRITE(*,*)' ITERATION IN PROCESS ...'
        CALL IMPROVE(NIT,S0,NOSC,DEL,FNEW,SMAT,NA,
     1  SCFAC,AMULT,Q,AS,NM,NC,ZNU)
        WRITE(*,*)' AGAIN (Y/N) ?'
        READ(*,'(A)')OK
        IF (OK.EQ.'Y'.OR.OK.EQ.'y')GOTO 333
        GOTO 111
       ENDIF
      endif
C
C     Calculate the true mass-weighted S-matrix:
      DO 210 IA=1,NA
C     ZM is in 1/amu^1/2
      T=ZM(IA)
      DO 210 IM=1,NA
C     SMAT now is in 1/amu^1/2      
210   SMAT(IA,IM)=SMAT(IA,IM)*T
      CALL SMATOUT(SMAT,AS,NA,JOB,ili,nlinepar,i1,i2)
      WRITE(*,*)' S-matrix written out ...'
      CALL FPC(Q,CHR,SMAT,NOAT)
      WRITE(*,*)' FPC.TAB written ...'
      CALL TERMO(Q,NA,i1,i2)
C
C     Read-in B-matrix:
      INQUIRE(FILE='BBI.MAT',EXIST=LBMAT)
      IF(LBMAT)THEN
       OPEN(15,FILE='BBI.MAT',FORM='UNFORMATTED')
C      SMAT used for B-matrix:
C      FNEW used for B^-1:
       READ(15) NA,NA,((SMAT(I,J),J=1,NA),I=1,NA),
     1 ((FNEW(I,J),J=1,NA),I=1,NA)
       CLOSE(15)
       WRITE(*,*)' B and B^-1 read in ...'
       CALL CIFF(NA,FCAR,FNEW,FINT,3,AMULT)
       CALL CIFF(NA,FINT,SMAT,FNEW,4,ZM)
      ENDIF
C
      write(6,6009)mclock()-it0
      END
C
      SUBROUTINE SMATOUT(S,AS,NAT3,JOB,ili,nlinepar,IRANGE,JRANGE)
C     dX  =  S dQ   S(nat3xnint)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      CHARACTER*4 JOB
      CHARACTER*6 num
      DIMENSION S(NAT3,NAT3),AS(3,NAT3)
      WRITE(*,*)
     1'RANGE OF MODES IN F.INP (Use the bracket#,negative delete zero):'
      IRANGE=1
      JRANGE=NAT3
      if(nlinepar.gt.ili+1)then
       ili=ili+1
       call getarg(ili,num)
       write(6,*)num
       read(num,*)IRANGE
       ili=ili+1
       call getarg(ili,num)
       write(6,*)num
       read(num,*)JRANGE
      else
       IF(JOB.NE.'AUTO')READ(*,*) IRANGE,JRANGE
      endif

      if(IRANGE.lt.0)then
       ize=1
      else
       ize=0
      endif
      IRANGE=ABS(IRANGE)
      JRANGE=ABS(JRANGE)

      if(ize.eq.1)then
       NI=0
       do 2 J=IRANGE,JRANGE
2      if(dabs(AS(1,J)).gt.0.1d0)NI=NI+1
      else
       NI=JRANGE-IRANGE+1
      endif
c
c     look at the phase and make the biggest element always positive:
      DO 101 J=IRANGE,JRANGE
      amx=S(1,J)
      DO 102 I=2,NAT3
102   if(abs(S(I,J)).gt.abs(amx))amx=S(I,J)
      if(amx.lt.0)then
       do 103 I=1,NAT3
103    S(I,J)=-S(I,J)
      endif
101   continue
c
      OPEN(34,FILE='F.INP')
      WRITE(34,10)NI,NAT3,NAT3/3
10    FORMAT(3I7)
      OPEN(35,FILE='FILE.X',STATUS='OLD')
      READ(35,*)
      READ(35,*)NAT
      DO 3 I=1,NAT
      READ(35,*)IAT,X,Y,Z
3     WRITE(34,11)IAT,X,Y,Z
11    FORMAT(I7,3F12.6)
      CLOSE(35)
      WRITE(34,14)
14    FORMAT(' Atom Mode    X-disp.    Y-disp.      Z-disp.')
      DO 1 I=1,NAT3/3
      DO 1 J=IRANGE,JRANGE
      if(dabs(AS(1,J)).gt.0.1d0.or.ize.eq.0)then
       IU=3*(I-1)
       K=NAT3-J+1
       WRITE(34,15)I,K,S(IU+1,J),S(IU+2,J),S(IU+3,J)
15     FORMAT(2I7,3F11.6)
      endif
1     continue
      WRITE(34,11)NI
      II=0
      do 4 I=IRANGE,JRANGE
      if(dabs(AS(1,I)).gt.0.1d0.or.ize.eq.0)then
       II=II+1
       WRITE(34,13)AS(1,I)
13     format(F11.3,$)
       if(mod(II,6).eq.0)write(34,*)
      endif
4     continue
      write(34,*)
      CLOSE(34)
      RETURN
      END

      SUBROUTINE IMPROVE(NIT,S0,NOSC,DEL,FNEW,SMAT,NA,
     1SCFAC,AMULT,Q,AS,NM,NC,ZNU)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      dimension FNEW(NA,NA),SMAT(NA,NA),
     1Q(*),AS(3,NA),AMULT(*),SCFAC(*),NM(*),
     2NC(NA/3,NA/3),ZNU(*)
      IIT=0
 1    IIT=IIT+1
      IF (IIT.GT.NIT)RETURN
      IF (MOD(IIT,10).EQ.0)THEN
       WRITE(16,*)' ITERATION ',IIT
       WRITE(16,*)' THE SCALE FACTORS:'
       WRITE(16,1616)(SCFAC(I),I=1,NOSC)
1616   FORMAT(6G12.6)
      ENDIF
      ICHANGE=0
      ie=0
      DO 2 I=1,NOSC
      X0=SCFAC(I)
      S0=SUMk(NOSC,0,ie,NA,FNEW,SMAT,ZNU,
     1AMULT,Q,AS,SCFAC,NM,NC)
      SCFAC(I)=X0-DEL
      SM=SUMk(NOSC,0,ie,NA,FNEW,SMAT,ZNU,
     1AMULT,Q,AS,SCFAC,NM,NC)
      SCFAC(I)=X0+DEL
      SP=SUMk(NOSC,0,ie,NA,FNEW,SMAT,ZNU,
     1AMULT,Q,AS,SCFAC,NM,NC)
      SCFAC(I)=X0
      WRITE(16,1600)X0,X0-DEL,X0+DEL
1600  format(3g15.5)
      WRITE(16,1600)S0,SM,SP
      IF (SP.LT.S0)THEN
       S0=SP
       SCFAC(I)=X0+DEL
       ICHANGE=1
      ENDIF
      IF (SM.LT.S0)THEN
       S0=SM
       SCFAC(I)=X0-DEL
       ICHANGE=1
      ENDIF
2     CONTINUE
      IF (ICHANGE.EQ.0)DEL=DEL/2.0D0
      WRITE(*,6000)IIT,DEL,S0
6000  FORMAT(I4,' del = ',F15.9,',  error = ',F15.4)
      IF (DEL.GT.0.0001d0)GOTO 1
      RETURN
      END

      FUNCTION SUMk(NOSC,iwr,iee,NA,FNEW,SMAT,ZNU,
     1AMULT,Q,AS,SCFAC,NM,NC)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension SMAT(NA,NA),NM(*),
     1AMULT(*),Q(*),AS(3,NA),FNEW(NA,NA),SCFAC(*),
     2NC(NA/3,NA/3),ZNU(*)
      real*8,allocatable::FINT(:,:)
c
      wtol=0.01d0
      c1=1302.828d0
      IF(NOSC.GT.0)THEN
       DO 201 I=1,NA
201    AMULT(I)=1.0d0
       DO 202 I=1,NOSC
       SF=1.0d0/SQRT(SCFAC(I))
       DO 202 J=1,NM(I)
       DO 202 K=1,3
202    AMULT(3*(NC(I,J)-1)+K)=SF
       allocate(FINT(NA,NA))
       CALL FAST(NA,FNEW,FINT,AMULT)
       do i=1,NA
        do j=1,NA
         FNEW(i,j)=FINT(i,j)
        enddo
       enddo
      ENDIF
      if(IWR.gt.0)WRITE(6,*)' New cartesian force field done '

      CALL TRED12(NA,FNEW,SMAT,Q,2,IERR)
      IF(IERR.NE.0)call report('diagonalization not ok')
c      
      RMSW=0.0D0
      ie=0
      DO 310 I=1,NA
      T=Q(I)
      SIGN=1.0d0
C     Make imaginary frequencies negative:
      IF(T.LT.0.0d0)SIGN=-1.0d0
      T=SIGN*SQRT(ABS(T))*c1
      Q(I)=T
      TP=ZNU(I)
      DEL=T-TP
      if(ABS(TP).le.wtol)then
       DEL=0.0d0
      else
       ie=ie+1
      endif
      RMSW=RMSW+DEL**2
      WRITE(16,1230) I, T,TP,DEL
1230  FORMAT(I6,3F15.2)
      AS(1,I)=T
      AS(2,I)=TP
310   AS(3,I)=T-TP
      iee=ie
      if(ie.gt.0)then
       RMSW=DSQRT(RMSW/DFLOAT(ie))
      else
       if(NA.gt.0)RMSW=DSQRT(RMSW/DFLOAT(NA))
      endif
      SO1=0.0D0
      SO2=0.0D0
      DO 5151 J=1,NA
      SO1=SO1+AS(1,J)*AS(2,J)
5151  SO2=SO2+AS(1,J)*AS(1,J)
      AK=SO1/SO2
C
      if(iwr.gt.0)WRITE(*,5555)
5555  FORMAT(1X,'************  NEW FREQUENCIES ******************')
      DO 311 J=NA,1,-1
      a3=0.0d0
      if(AS(2,NA-J+1).gt.wtol)a3=AS(3,NA-J+1)
311   if(iwr.gt.0)WRITE(*,6666)J,NA-J+1,(AS(I,NA-J+1),I=1,2),a3
6666  FORMAT(1X,
     1 I5,' (',I4,') ...... ',F8.2,'   -',F8.2,'/',F8.2)
      if(iwr.gt.0)WRITE(*,5551)
5551  FORMAT(1X,'************************************************')
      if(iwr.gt.0)WRITE(6,6667)ie,RMSW,AK
6667  FORMAT(i6,'experimental frequencies found',/,
     1 ' >>>>><< >> ',F8.2,' <<<<<<  ',F8.6)
c
      SUMk=RMSW
      RETURN
      END

      SUBROUTINE FPC(V,CHR,S,NOAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      real*8 V(*),S(3*NOAT,3*NOAT),X(3,NOAT),CHR(*),DS(3),DR(3)
      DCOR = 388891.4D0
      RCOR = 12221.7D0
      OPEN(35,FILE='FILE.X',STATUS='OLD')
      READ(35,*)
      READ(35,*)
      DO 51 J=1,NOAT
51    READ(35,*)X(1,J),(X(I,J),I=1,3)
      CLOSE(35)
      OPEN(13,FILE='FPC.TAB')
      WRITE(13,50)
  50  FORMAT(3X,'FREQ',8X,'DIP. STRN.*E39',6X,'ROT. STRN.*E45'/
     1 3X,'CM-1',2X,2(6X,' (FR  CM)**2  '),/,'---------------')
      DO 42 I=1,3*NOAT
      DO 784 IKY=1,3
      DS(IKY)=0.0D0
      DR(IKY)=0.0D0
      DO 784 JR=1,NOAT
      DS(IKY)=DS(IKY)+CHR(JR)*S(3*(JR-1)+IKY,I)
      DO 784 IB=1,3
      DO 784 IG=1,3
784   DR(IKY)=DR(IKY)+EPS(IKY,IB,IG)*X(IB,JR)*S(3*(JR-1)+IG,I)
     1                       *CHR(JR)
      DT =DS(1)*DS(1)+DS(2)*DS(2)+DS(3)*DS(3)
      ROT=DS(1)*DR(1)+DS(2)*DR(2)+DS(3)*DR(3)
      IF(V(I).GT.1.0d0)then
       WRITE(13,85) I,V(I),DT*DCOR/ABS(V(I)),ROT*RCOR
85     FORMAT(I4,F8.2,' ',F13.6,' ',f13.6)
      endif
42    continue
      WRITE(13,86)
86    FORMAT('---------------------------------------------')
      CLOSE(13)
      RETURN
      END

      FUNCTION EPS(I,J,K)
      integer*4 I,J,K
      REAL*8 EPS
      EPS=0.0D0
      IF (I.EQ.1.AND.J.EQ.2.AND.K.EQ.3)EPS= 1.0D0
      IF (I.EQ.1.AND.J.EQ.3.AND.K.EQ.2)EPS=-1.0D0
      IF (I.EQ.2.AND.J.EQ.3.AND.K.EQ.1)EPS= 1.0D0
      IF (I.EQ.2.AND.J.EQ.1.AND.K.EQ.3)EPS=-1.0D0
      IF (I.EQ.3.AND.J.EQ.1.AND.K.EQ.2)EPS= 1.0D0
      IF (I.EQ.3.AND.J.EQ.2.AND.K.EQ.1)EPS=-1.0D0
      RETURN
      END

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
C
      SUBROUTINE CIFF(NA,A,B,C,IC,U)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(NA,NA),B(NA,NA),C(NA,NA),U(*)
      real*8,allocatable::TEM(:,:),UM(:)
      allocate(TEM(NA,NA),UM(NA))
C
C     IC: 1 ... calculate intrinsic FF      C   =     U.Bt  .  A .B  . U
C                                          FINT =    SF.B-1t.FCAR.B-1.SF
C         2 ... calculate mass-weighted cartesian FF
C                                          FCAR=sqrt(M)-1.Bt  .FINT.B.sqrt(M)-1
C         3 ... calculate Intrinsic FF and write it
C         4 ... calculate cartesian FF and write it in right units
C                                          FCAR=        Bt  .FINT.B
C
      DO 1 I=1,NA
1     UM(I)=U(I)
      IF(IC.EQ.4)THEN
       DO 2 I=1,NA
2      UM(I)=1.0d0
      ENDIF
      DO 212 J=1,NA
      UJ=UM(J)
      DO 212 I=1,NA
      TEM(I,J)=0.0d0
      DO 212 II=1,NA
212   TEM(I,J)=TEM(I,J)+A(I,II)*B(II,J)*UJ
      DO 214 I=1,NA
      UI=UM(I)
      DO 214 J=1,NA
      C(I,J)=0.0d0
      DO 214 II=1,NA
214   C(I,J)=C(I,J)+B(II,I)*TEM(II,J)*UI
      IF(IC.EQ.3)THEN
       OPEN(7,FILE='INTY.FC')
       WRITE(7,3337)(NA*(NA+1))/2,NA,NA
3337   FORMAT(/,I6,' INTRINSIC FORCE FIELD',I5,' x (',I5,' + 1) / 2')
       M=6
       MI=-5
2071   MI=MI+6
       MJ=M-5
       WRITE(7,28)(I,I=MJ,MIN(NA,M))
28     FORMAT(8X,6(I7,5X))
       DO 226 I=MI,NA
226    WRITE(7,7000)I,(C(J,I),J=MJ,MIN(I,M))
7000   FORMAT(I5,3X,6F12.6)
       IF (M.LT.NA)THEN
        M=M+6
        GOTO 2071
       ENDIF
       CLOSE(7)
       WRITE(*,*)' INTY.FC with intrinsic coordinates written ...'
      ENDIF
      IF(IC.EQ.4)THEN
       OPEN(20,FILE='NEW.FC')
       CALL WRITEFF(NA,C)
       CLOSE(20)
       WRITE(*,*)' NEW.FC written out ...'
      ENDIF
      RETURN
      END
C
      SUBROUTINE READFF(N,FCAR)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FCAR(N,N)
      OPEN(20,FILE='FILE.FC',STATUS='OLD')
      N1=1
1     N3=N1+4
      IF(N3.GT.N)N3=N
      DO 130 LN=N1,N
130   READ(20,17)(FCAR(LN,J),J=N1,MIN(LN,N3))
      N1=N1+5
      IF(N3.LT.N)GOTO 1
17    FORMAT(4X,5D14.6)
C
      CONST=4.359828d0/0.5291772d0/0.5291772d0 !conversion (d^2E/dx^2) attoJoule/Angstrom^2 => au(Hartree)/au(Bohr)^2 (I think)
      DO 3 I=1,N
      DO 3 J=I,N
3     FCAR(J,I)=FCAR(J,I)*CONST
      DO 31 I=1,N
      DO 31 J=I+1,N
   31 FCAR(I,J)=FCAR(J,I)
      CLOSE(20)
      WRITE(*,*)' Cartesian FF read in ... '
      RETURN
      END

      SUBROUTINE WRITEFF(N,FCAR)
      IMPLICIT INTEGER*4 (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FCAR(N,N)
      CONST=4.359828d0/0.5291772d0**2
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

      SUBROUTINE FAST(NA,FCAR,FNEW,ZM)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4(I-N)
      DIMENSION FCAR(NA,NA),FNEW(NA,NA),ZM(*)
      DO 214 I=1,NA
      UI=ZM(I)
      DO 214 J=I,NA
      FNEW(I,J)=FCAR(I,J)*ZM(J)*UI
214   FNEW(J,I)=FNEW(I,J)
      RETURN
      END
C
      SUBROUTINE TERMO(W,NA,i1,i2)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION W(NA)
      character*2 key
      logical lex
C     Thermodynamics according to Mopac notation
      Cl=1.4388d0
      w0=0.0d0
      T=298.15d0
C     1cm-1=1.4388 K
c 
      inquire(file='THERMO.OPT',exist=lex)
      if(lex)then
       open(99,file='THERMO.OPT')
3      read(99,990,end=99,err=99)key
990    format(a2)
       if(key.eq.'w0')read(99,*)w0
       if(key.eq.'TE')read(99,*)T
       goto 3
99     close(99)
      endif
C
      Ecm=0.0d0
      Cl=1.4388d0/T
      Qv=0.0d0
      Ucm=0.0d0
      Scm=0.0d0
      C=0.0d0
c
c     avoiding divergence of entropy of the low-frequency vibration
c     modes, see Grimme, S. Chem. Eur. J. 2012, 18, 9955-9964, Supramolecular 
c     Binding Thermodynamics by Dispersion-Corrected Density Functional Theory
      DO 2 I=i1,i2
      o=dabs(max(W(I),0.0d0))
      ot=o+1.0d-9
      EWJ=0.0d0
      EWJt=0.0d0
      IF(o*Cl.LT.200.0d0)EWJ=EXP(-o*Cl)
      IF(o*Cl.LT.200.0d0)EWJt=EXP(-ot*Cl*(1.0d0+w0/ot))
      Ecm=Ecm+o* 0.5d0
      if(dabs(1.0d0-EWJ).gt.1.0d-10)then
       Qv=Qv+1.0d0/(1.0d0-EWJ)
       Ucm=Ucm+o*(0.5d0+EWJ/(1.0d0-EWJ))
       Scm=Scm+ot*EWJt/(1.0d0-EWJt)/T-0.69503476d0*Log(1.0d0-EWJt)
       C=C+o**2*EWJ/(1.0d0-EWJ)**2
      endif
2     continue
c
      Gcm=Ucm-T*Scm
      Eau=Ecm/219470.0d0
      Uau=Ucm/219470.0d0
      Sau=Scm/219470.0d0
      Gau=Gcm/219470.0d0
      Ukcal=Uau*627.5d0
      Skcal=Sau*627.5d0
      Gkcal=Gau*627.5d0
      Ekcal=Eau*627.5d0
C     cal/mol/K
      C=C*Cl**2*1.987216d0
C     cal/mol/K
      WRITE(16,16001)T,w0
16001 FORMAT(/,' Thermochemistry parameters',/,
     1       ' Temperature :',f12.3,' K',/,
     1       ' w0          :',f12.3,' cm-1',/,
     224x,'    kcal/mol        cm-1     hartree',/,80(1H-))
      WRITE(16,16000)'Internal energy 0K U0:',Ekcal,Ecm,Eau
      WRITE(16,16000)'Internal energy TK U :',Ukcal,Ucm,Uau
      WRITE(16,16000)'Entropy            S :',Skcal,Scm,Sau
      WRITE(16,16000)'Gibbs              G :',Gkcal,Gcm,Gau
      WRITE(16,16000)'Heat capacity      C :',C/1000.0d0
16000 FORMAT(a22,F14.5,F12.1,F12.6)
      WRITE(16,16002)
16002 format(80(1H-))
      RETURN
      END
C
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

      SUBROUTINE PROJECT(A,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N),F(N,N)
      real*8,allocatable::TEM(:,:)
C
      call TRERR(F,N)
      allocate(TEM(N,N))
      WRITE(*,*)'Projecting transl/rotations from the force field...'
      CALL DOMA(A,N)
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
      RETURN
      END

      SUBROUTINE DOMA(A,N)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION A(N,N)
      real*8,allocatable::C(:,:),TEM(:,:)
C
      N=N
      allocate(C(3,N/3),TEM(N,N))
      AMACH=0.00000000000001d0
      NAT=N/3
      OPEN(4,FILE='FILE.X',STATUS='OLD')
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
C     ============================================================ 
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
      SP=0.0d0
      DO 5 J=1,NAT3
5     SP=SP+A(J,IC)*A(J,ICP)
C     Subtract projection to ICP:
      DO 6 J=1,NAT3
6     A(J,IC)=A(J,IC)-SP*A(J,ICP)
C     Normalize :
      SP=0.0d0
      DO 7 J=1,NAT3
7     SP=SP+A(J,IC)**2
      IF(SP.GT.AMACH)SP=1.0d0/SQRT(SP)
      DO 3 J=1,NAT3
3     A(J,IC)=A(J,IC)*SP
      RETURN
      END

      function mclock()
      integer*4 mclock
      mclock=0
      return
      end
C     ============================================================ 
      subroutine readinp(JOB,ISYM,NOSC,SCFAC,NM,NC,CHR,NOAT,NA,NNU,ZIM,
     1ZM,ZNU)
      implicit none
      character*4 JOB
      integer*4 ISYM,NOSC,NM(*),NOAT,I,J,NC(NOAT,NOAT),NMOL,NNU,NA
      real*8 SCFAC(*),CHR(*),CHATOT,ZIM(*),BM,ZM(*),ZNU(*)
      OPEN(15,FILE='FTRY.INP')
      READ(15,1500)JOB,ISYM
1500  FORMAT(A4,I4)
      READ(15,*)NOSC
      DO 222 I=1,NOSC
222   READ(15,*) SCFAC(I),NM(I),(NC(I,J),J=1,NM(I))
1040  FORMAT(6G12.6)
      WRITE(*,*)' Scale factors read in...'
C
C     Read-in atomic charges for FPC:
c     READ(15,1040) (CHR(I),I=1,NOAT)
      READ(15,*) (CHR(I),I=1,NOAT)
      CHATOT=0.0d0
      DO 77 I=1,NOAT
77    CHATOT=CHATOT+CHR(I)
      WRITE(*,777)CHATOT,NOAT
777   FORMAT(' TOTAL CHARGE',F20.6,/I6,' atoms')
C
C     Read-in atomic masses:
      READ(15,*)NMOL
      IF(NMOL.NE.1)THEN
       CLOSE(16)
       CLOSE(15)
       call report(' Only one molecule allowed here !')
      ENDIF
      READ(15,1040)(ZIM(I),I=1,NOAT)
      WRITE(16,1040)(ZIM(I),I=1,NOAT)
      BM=0.0d0
      DO 102 I=1,NOAT
102   BM=BM+ZIM(I)
      WRITE(*,778)BM
778   FORMAT(' TOTAL MASS  ',F20.6)
      DO 220 I=1,NOAT
      DO 220 J=1,3
220   ZM(3*I-3+J)=1/SQRT(ZIM(I))
C
C     Read-in experimental frequencies:
      READ(15,*)NNU
      IF(NNU.EQ.0)NNU=NA-6
      READ(15,1040)(ZNU(I),I=1,NNU)
      DO 211 I=NNU+1,NA
211   ZNU(I)=0.0d0
      CLOSE(15)
      WRITE(*,*)' FTRY.INP READ IN'
      return
      end
C     ============================================================ 
      subroutine SYMROUTE(nat,n,f,s,AS)
      implicit none
      real*8 f(n,n),s(n,n),AS(3,n),cm(3)
      real*8,allocatable::r(:),line(:),u(:,:),ft(:,:),fs(:,:),
     1fsmall(:,:),ssmall(:,:),esmall(:),ar(:,:),ac(:,:)
      integer*4,allocatable::q(:),table(:,:),nr(:),nc(:),ir(:,:),ic(:,:)
      integer*4 ns0,i,nat,n,j,ix,iy,iz,mx,n6,ia,iap,ii,ee,ib,ie,ifa,
     1ig,im,ixa,nb,nt,coord,ixp,iip,jj,ndim,k,nm,i1,j1,imax
      parameter (ns0=10,mx=12,n6=6)
      real*8 o(ns0,3,3),x,y,z,xs,ys,zs,xp,yp,zp,tol,d,ln,sn,t
      character*3 sst(ns0),SSP(ns0,mx),so(ns0),s3
      character*3,allocatable::label(:)
      character*1 xyz(3)
      data xyz/'x','y','z'/
      integer*4 CHT(ns0,mx,mx),KG(ns0),se(ns0,mx),nstart(mx),nend(mx)
      logical debug

      write(6,598)
 598  format(/,' Symmetry determination',/)
      allocate(r(n),q(nat),fs(n,n),ft(n,n),label(n))
      allocate(nr(n),nc(n),ar(n,mx),ac(n,mx),ic(n,mx),ir(n,mx))
c     distance limit = 10-2A, tol=(10-2)^2
      tol=1.0d-4
      debug=.false.

      cm=0.0d0
      open(9,file='FILE.X')
      read(9,*)
      read(9,*)nat
      do 1 i=1,nat
      read(9,*)q(i),(r(j+3*(i-1)),j=1,3)
      do 1 j=1,3
  1   cm(j)=cm(j)+r(j+3*(i-1))
      do 101 j=1,3
101   cm(j)=cm(j)/dble(nat)
      do 102 i=1,nat
      do 102 j=1,3
102   r(j+3*(i-1))=r(j+3*(i-1))-cm(j)
      close(9)

      call setel(ns0,o,so)
      call SETSYM(CHT,sst,KG,SSP,se)

c     loop over symmetry point groups:
c     find one with maximum symmetry elements:
      ifa=0
      im=0
      do 21 ig=1,n6

      if(debug)write(6,6001)ig,sst(ig)
6001  format(i2,' group ',A3) 

c     try all axis combinations:
      do 21 ix=1,3
      iy=ix+1
      if(iy.eq.4)iy=1
      iz=iy+1
      if(iz.eq.4)iz=1

      if(debug)write(6,6002)ix,iy,iz
6002  format('xyz:',3i2)

c     loop over symmetry elements:
      do 2 ie=1,KG(ig)
      ee=se(ig,ie)
      if(debug)write(6,6003)ie,ee
6003  format(i2,' element,',i2,' operation')
c     loop over atoms:
      do 31 ia=1,nat
      x=r(ix+3*(ia-1))
      y=r(iy+3*(ia-1))
      z=r(iz+3*(ia-1))
c     symmetry-related coordinate
      xs=o(ee,1,1)*x+o(ee,1,2)*y+o(ee,1,3)*z
      ys=o(ee,2,1)*x+o(ee,2,2)*y+o(ee,2,3)*z
      zs=o(ee,3,1)*x+o(ee,3,2)*y+o(ee,3,3)*z
      if(debug)write(6,6004)ia,x,y,z,xs,ys,zs
6004  format(i2,'atom',6f10.3)

c     is in this position of an atom?:
      do 3 iap=1,nat
      xp=r(ix+3*(iap-1))
      yp=r(iy+3*(iap-1))
      zp=r(iz+3*(iap-1))
      d=(xp-xs)**2+(yp-ys)**2+(zp-zs)**2
c     it is, go to another atom: 
3     if(d.le.tol.and.q(ia).eq.q(iap))goto 31
      if(debug)write(6,*)'not found'
c     it is not, try another symmetry group:
      goto 21
31    continue
c     ia
2     continue
c     ie

c     now for all symmetry elements atoms could be assigned
      write(6,601)sst(ig),KG(ig)
601   format(' ',a3,' group found, Nel = ',i2)
      if(KG(ig).gt.im)then
       im=KG(ig)
       ifa=ig
       ixa=ix
      endif
21    continue
c     ig,ix


      ix=ixa
      iy=ix+1
      if(iy.eq.4)iy=1
      iz=iy+1
      if(iz.eq.4)iz=1
      write(6,602)sst(ifa)
602   format(/,' ',a3,' is the best')
c     where atoms go under symmetry operations:
      allocate(table(n,im))
      if(debug)write(6,6002)ix,iy,iz
      table=0
      do 4 ie=1,KG(ifa)
      ee=se(ifa,ie)
      if(debug)write(6,6003)ie,ee
c     loop over atoms:
      do 4 ia=1,nat
      ii=3*(ia-1)
      x=r(ix+ii)
      y=r(iy+ii)
      z=r(iz+ii)
c     symmetry-related coordinate
      xs=o(ee,1,1)*x+o(ee,1,2)*y+o(ee,1,3)*z
      ys=o(ee,2,1)*x+o(ee,2,2)*y+o(ee,2,3)*z
      zs=o(ee,3,1)*x+o(ee,3,2)*y+o(ee,3,3)*z
      if(debug)write(6,6004)ia,x,y,z,xs,ys,zs
      do 4 iap=1,nat
      xp=r(ix+3*(iap-1))
      yp=r(iy+3*(iap-1))
      zp=r(iz+3*(iap-1))
      d=(xs-xp)**2+(ys-yp)**2+(zs-zp)**2
      if(d.le.tol.and.q(ia).eq.q(iap))then
c      look at symmetry-related coordinate changes
       if(o(ee,1,1).gt.0.0d0)then
        table(ix+ii,ie)= ix+3*(iap-1)
       else
        table(ix+ii,ie)=-ix-3*(iap-1)
       endif
       if(o(ee,2,2).gt.0.0d0)then
        table(iy+ii,ie)= iy+3*(iap-1)
       else
        table(iy+ii,ie)=-iy-3*(iap-1)
       endif
       if(o(ee,3,3).gt.0.0d0)then
        table(iz+ii,ie)= iz+3*(iap-1)
       else
        table(iz+ii,ie)=-iz-3*(iap-1)
       endif
      endif
4     continue

      write(6,603)(so(se(ifa,i)),i=1,KG(ifa))
603   format(/,'  coord   atom xyz  ',12(3x,a3))
      do 5 ia=1,nat
      do 5 ix=1,3
      ii=ix+3*(ia-1)
5     write(6,604)ii,ia,xyz(ix),(table(ii,i),i=1,KG(ifa))
604   format(2i7,3x,a1,2x,12i6)

c     transformation matrix:
      allocate(u(n,n))
      u=0.0d0
c     symmetry blocks:
      allocate(line(n))
      nt=0
      do 6 ib=1,KG(ifa)
      write(6,605)SSP(ifa,ib)
605   format(' Block ',a3)
      nb=0
      nstart(ib)=0
      do 71 ia=1,nat
      do 71 ix=1,3
      ii=ix+3*(ia-1)
      line=0.0d0
      do 7 i=1,KG(ifa)
      coord=iabs(table(ii,i))
      if(table(ii,i).lt.0)then
       sn=-1.0d0
      else
       sn= 1.0d0
      endif
7     line(coord)=line(coord)+sn*dble(CHT(ifa,ib,i))
      ln=0.0d0
      do 8 ii=1,n
 8    ln=ln+line(ii)**2
      ln=dsqrt(ln)
      if(ln.gt.tol)then
       do 9 ii=1,n
 9     line(ii)=line(ii)/ln


       nb=nb+1
       nt=nt+1
       write(6,606)nb
606    format(i7,$)
       do 10 iap=1,nat
       do 10 ixp=1,3
       ii=ixp+3*(iap-1)
 10    if(dabs(line(ii)).gt.tol)write(6,607)line(ii),iap,xyz(ixp)
607    format(f10.3,i6,a1,$)
c
c      record big elements:
       nr(1)=0
       do 112 jj=1,n
       if(dabs(line(jj)).gt.tol)then
        nr(1)=nr(1)+1
        if(nr(1).gt.mx)call report('nr > mx')
        ar(1,nr(1))=line(jj)
        ir(1,nr(1))=jj
       endif
112    continue

c      check if already defined
       do 11 iip=1,nt-1
       sn=0.0d0
       do 111 j1=1,nr(1)
111    sn=sn+ar(1,j1)*u(iip,ir(1,j1))
       if(dabs(sn).gt.tol)then
        write(6,610)iip
610     format(' same as',i7,', deleted')
        nb=nb-1
        nt=nt-1
        goto 71
       endif
11     continue

       if(nt.le.n)then
        if(nstart(ib).eq.0)nstart(ib)=nt
        nend(ib)=nt
        do 12 jj=1,n
12      u(nt,jj)=line(jj)
        write(6,*)' added'
       else
        write(6,*)'Error: too many symmetrized coordinates'
        stop
       endif
      endif
71    continue
c     ia,ix
6     write(6,608)nb,nstart(ib),nend(ib)
608   format(i8,' elements, from',i7,' to',i8)
c     ib

      write(6,*)
      write(6,609)nt
609   format(i8,' in total')

      call sparse(n,mx,nr,nc,ar,ac,tol,ir,ic,u)
      do 23 i=1,n
23    write(6,800)i,nr(i),(ir(i,j),j=1,nr(i))
800   format(20i7)

      write(6,*)'Symmetry transformation of force field'
      do 13 i=1,n
      do 13 jj=1,n
      ft(i,jj)=0.0d0
      do 13 i1=1,nr(i)
13    ft(i,jj)=ft(i,jj)+ar(i,i1)*f(ir(i,i1),jj)
      do 14 i=1,n
      do 14 j=1,n
      fs(i,j)=0.0d0
      do 14 j1=1,nr(j)
14    fs(i,j)=fs(i,j)+ar(j,j1)*ft(i,ir(j,j1))
      write(6,*)'Done'
c     OPEN(20,FILE='SYM.FC')
c     CALL WRITEFF(n,fs)
c     CLOSE(20)
c     write(6,*)'SYM.FC written'

      s=0.0d0
      nm=0
      do 15 ib=1,KG(ifa)
      write(6,605)SSP(ifa,ib)
      ndim=nend(ib)-nstart(ib)+1
      allocate(fsmall(ndim,ndim),ssmall(ndim,ndim),esmall(ndim))

      do 16 i=1,ndim
      do 16 j=1,ndim
 16   fsmall(i,j)=fs(nstart(ib)+i-1,nstart(ib)+j-1)

      call fdiag(ndim,fsmall,ssmall,esmall)

      do 17 i=1,ndim
      nm=nm+1
      label(nm)=SSP(ifa,ib)
      AS(1,nm)=esmall(i)
      do 17 j=1,ndim
      do 17 i1=1,nr(nstart(ib)+j-1)
      ii=ir(nstart(ib)+j-1,i1)
17    s(ii,nm)=s(ii,nm)+ar(nstart(ib)+j-1,i1)*ssmall(j,i)
c7    s(ii,nm)=s(ii,nm)+u(nstart(ib)+j-1,ii)*ssmall(j,i)

15    deallocate(fsmall,ssmall,esmall)

c     order frequencies:
      do 19 i=1,n
      imax=i
      do 191 j=i+1,n
191   if(AS(1,j).gt.AS(1,imax))imax=j
      if(imax.ne.i)then
       t=AS(1,imax)
       AS(1,imax)=AS(1,i)
       AS(1,i)=t
       s3=label(imax)
       label(imax)=label(i)
       label(i)=s3
       do 20 k=1,n
       t=s(k,imax)
       s(k,imax)=s(k,i)
20     s(k,i)=t
      endif
 19   continue

      WRITE(*,5555)
5555  FORMAT(/,' ************  NEW FREQUENCIES ******************')
      do 22 i=1,n
 22   write(6,6666)n-i+1,i,AS(1,i),label(i)
6666  FORMAT(1X,I5,' (',I4,') ...... ',F8.2,3x,a3)
      WRITE(*,5551)
5551  FORMAT(' ************************************************')

      return

      end

      subroutine fdiag(n,f,s,e)
      IMPLICIT none
      REAL*8 f(n,n),s(n,n),e(n),T,SIGN,c1
      integer*4 n,i,IERR,J
      c1=1302.828d0
      CALL TRED12(n,f,s,e,2,IERR)
      IF(IERR.NE.0)call report('diagonalization error')
      DO 310 I=1,n
      T=e(I)
      SIGN=1.0d0
      IF(T.LT.0.0d0)SIGN=-1.0d0
      T=SIGN*SQRT(ABS(T))*c1
310   e(I)=T
      DO 311 J=1,n
311   WRITE(*,6666)J,e(J)
6666  FORMAT(i10,F8.2)
      RETURN
      END

      SUBROUTINE SETSYM(CHT,sst,KG,SSP,se)
      integer*4 ns0,mx
      parameter (ns0=10,mx=12)
      integer*4 CHT(ns0,mx,mx),KG(ns0),se(ns0,mx)
      character*3 sst(ns0),SSP(ns0,mx)
      character*1 PRIME
      PRIME=CHAR(39)                               
      CHT=1
      sst='xxx'
      se=0
c     sst .. point group
c     KG .. number of symmetry elements
c     SSP .. symmetry

      sst(1)='C1 '
      KG(1)=1                 
      SSP(1,1)='A  '
      se(1,1)=1

      sst(2)='C2V'
      KG(2)=4                 
      CHT(2,2,3)=-1                                  
      CHT(2,2,4)=-1                                  
      CHT(2,3,2)=-1                                  
      CHT(2,3,4)=-1                                  
      CHT(2,4,2)=-1                                  
      CHT(2,4,3)=-1                                  
C       1 1 1 1 A1                                 
C       1 1-1-1 A2                                 
C       1-1 1-1 B1                                 
C       1-1-1 1 B2                                 
      SSP(2,1)='A1 '  
      SSP(2,2)='A2 ' 
      SSP(2,3)='B1 ' 
      SSP(2,4)='B2 ' 
c     E C2Z sigmaxz sigmazy
      se(2,1)=1
      se(2,2)=2
      se(2,3)=3
      se(2,4)=4

      sst(3)='C2H'
      KG(3)=4                 
      CHT(3,2,2)=-1                                  
      CHT(3,2,4)=-1                                  
      CHT(3,3,3)=-1                                  
      CHT(3,3,4)=-1                                  
      CHT(3,4,2)=-1                                  
      CHT(3,4,3)=-1                                  
C       1 1 1 1 Ag                                 
C       1-1 1-1 Bg                                 
C       1 1-1-1 Au                                 
C       1-1-1 1 Bu                                 
      SSP(3,1)='Ag '                               
      SSP(3,2)='Bg '                               
      SSP(3,3)='Au '                               
      SSP(3,4)='Bu '                               
c     E C2Z i sigmah
      se(3,1)=1
      se(3,2)=2
      se(3,3)=5
      se(3,4)=6

      sst(4)='C2 '
      KG(4)=2                                        
      CHT(4,2,2)=-1                                 
C       1 1 A                                      
C       1-1 B                                      
      SSP(4,1)='A  '                                
      SSP(4,2)='B  '                                
c     E C2
      se(4,1)=1
      se(4,2)=2

      sst(5)='CS '
      KG(5)=2                                        
      CHT(5,2,2)=-1                                 
C       1 1 A'                                     
C       1-1 A''                                    
      SSP(5,1)='AP '                                
      SSP(5,1)(2:2)=PRIME                           
      SSP(5,2)='APP'                                
      SSP(5,2)(2:2)=PRIME                           
      SSP(5,2)(3:3)=PRIME                           
c     E sigmah
      se(5,1)=1
      se(5,2)=6
     
      sst(6)='D2H'
      KG(6)=8
      CHT(6,2,3)=-1                                 
      CHT(6,2,4)=-1                                 
      CHT(6,2,7)=-1                                 
      CHT(6,2,8)=-1                                 
      CHT(6,3,2)=-1                                 
      CHT(6,3,4)=-1                                 
      CHT(6,3,6)=-1                                 
      CHT(6,3,8)=-1                                 
      CHT(6,4,2)=-1                                 
      CHT(6,4,3)=-1                                 
      CHT(6,4,6)=-1                                 
      CHT(6,4,7)=-1                                 
      CHT(6,5,5)=-1                                 
      CHT(6,5,6)=-1                                 
      CHT(6,5,7)=-1                                 
      CHT(6,5,8)=-1                                 
      CHT(6,6,3)=-1                                 
      CHT(6,6,4)=-1                                 
      CHT(6,6,5)=-1                                 
      CHT(6,6,6)=-1                                 
      CHT(6,7,2)=-1                                 
      CHT(6,7,4)=-1                                 
      CHT(6,7,5)=-1                                 
      CHT(6,7,7)=-1                                 
      CHT(6,8,2)=-1                                 
      CHT(6,8,3)=-1                                 
      CHT(6,8,5)=-1                                 
      CHT(6,8,8)=-1                                 
C      1 1 1 1 1 1 1 1 Ag   1                      
C      1 1-1-1 1 1-1-1 B1g  2                      
C      1-1 1-1 1-1 1-1 B2g  3                      
C      1-1-1 1 1-1-1 1 B3g  4                      
C      1 1 1 1-1-1-1-1 Au   5                      
C      1 1-1-1-1-1 1 1 B1u  6                      
C      1-1 1-1-1 1-1 1 B2u  7                      
C      1-1-1 1-1 1 1-1 B3u  8                      
      SSP(6,1)='Ag '                                
      SSP(6,2)='B1g'                                
      SSP(6,3)='B2g'                                
      SSP(6,4)='B3g'                                
      SSP(6,5)='Au '                                
      SSP(6,6)='B1u'                                
      SSP(6,7)='B2u'                                
      SSP(6,8)='B3u'                                
c     E c2z c2y c2x i xy xz yz
      se(6,1)=1
      se(6,2)=2
      se(6,3)=7
      se(6,4)=8
      se(6,5)=5
      se(6,6)=6
      se(6,7)=3
      se(6,8)=4
      RETURN
      END  
      
      subroutine setel(ns0,o,so)
      implicit none
      integer*4 ns0,j,i
      real*8 o(ns0,3,3)
      character*3 so(ns0)
      o=0.0d0
c     1 identity, E:
      do 1 i=1,ns0
      do 1 j=1,3
1     o(i,j,j) =  1.0d0
      so(1)='E  '
c     2 C2z:
      o(2,1,1) = -1.0d0
      o(2,2,2) = -1.0d0
      so(2)='C2z'
c     3 sigmazx:
      o(3,2,2) = -1.0d0
      so(3)='szx'
c     4 sigmazy:
      o(4,1,1) = -1.0d0
      so(4)='szy'
c     5 inversion
      o(5,1,1) = -1.0d0
      o(5,2,2) = -1.0d0
      o(5,3,3) = -1.0d0
      so(5)='i  '
c     6 sigmaxy:
      o(6,3,3) = -1.0d0
      so(6)='sxy'
c     7 C2y:
      o(7,1,1) = -1.0d0
      o(7,3,3) = -1.0d0
      so(7)='C2y'
c     8 C2x:
      o(8,2,2) = -1.0d0
      o(8,3,3) = -1.0d0
      so(8)='C2x'
      return
      end

      subroutine report(s)
      character*(*) s
      write(6,*)s
      stop
      end

      subroutine sparse(n,mx,nr,nc,ar,ac,tol,ir,ic,u)
      implicit none
      integer*4 n,nr(*),nc(*),mx,i,j,ir(n,mx),ic(n,mx)
      real*8 ar(n,mx),ac(n,mx),u(n,n),tol
      do 1 i=1,n
c     number of non-zero elements in row i
      nr(i)=0
      do 1 j=1,n
      if(dabs(u(i,j)).gt.tol)then
       nr(i)=nr(i)+1
       if(nr(i).gt.mx)call report('mx overflow in sparse')
       ar(i,nr(i))=u(i,j)
       ir(i,nr(i))=j
      endif
  1   continue
      do 2 i=1,n
c     number of non-zero elements in column i
      nc(i)=0
      do 2 j=1,n
      if(dabs(u(j,i)).gt.tol)then
       nc(i)=nc(i)+1
       if(nc(i).gt.mx)call report('mx overflow in sparse')
       ac(i,nc(i))=u(j,i)
       ic(i,nc(i))=j
      endif
  2   continue
      write(6,*)'sparse part selected'
      return
      end
