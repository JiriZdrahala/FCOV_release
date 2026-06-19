#define TR(arg) trim(adjustl(arg))
      program pmz
      implicit none
      integer*4 MENDELEV
      parameter (MENDELEV=89)
      character*80 s80,ts,sdalbas
      character*20 s20
      character*1 sgn,axis
      character*20 i2str,coord
      logical lex,lcp2k,lcastep,qdim,ldalton,lvstep,lschk,only_double,
     1 no_zero,nochk
      real*8,allocatable:: r(:,:),s(:,:),w(:),steps(:)
      real*8 step,x,y,z,dx,dy,step1,step2,
     1dz,drmax,wcm,ucm,dd1,dd2,xf,yf,zf,o11i,o12i,o13i,o22i,
     2o23i,o33i,dist,dd,dqs
      integer*4,allocatable::ig(:),iqlist(:),nl(:)
      integer*4 i,j,nat,idiff,ic,ix,ist,io,
     1ii,im,ndiff,nm,nmd,jm,jx,jst,nparallel,npart,iparallel,is,ie,
     1ip1,ia1,ia2,nn,ico1,ico2,ip2,iy,istart,jb1,jb2,
     1ndiff1,ndiff2,ib,jb,jo,fistart,n,io_start,at
      double precision, parameter :: cm_2_au=4.556431403d-6
      double precision, parameter :: amu_2_au=1822.89d0
      double precision, parameter :: amu_2_kg=1.6605402E-27
      double precision, parameter :: ang_2_au=1.8897259886d0
      double precision, parameter :: c_au=1.37036d2
      double precision, parameter :: hbar=1.054571628d-34
      double precision, parameter :: c=299792458
      double precision, parameter :: pi=4.0d0*atan(1.0d0)

      CHARACTER*2 sy(MENDELEV)
      data sy/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
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
     7'Fr','Ra','Ac'/
      character*7 round(2)
      character*1 r3(3)
      character*10 s10
      character*1 xyz(3)
      character*4 key,k2,tsn
      data round/ 'forward','   back'/
      data r3/'-',' ','+'/
      data xyz/'x','y','z'/
      character*80 tso

      write(6,6000)
6000  format('  Preparation of inputs using coordinate differentiations'
     1 ,/,/, '      Input: FILE.X     - geometry',/,
     2       '                G.TXT   - gaussian input head',/,
     3       '               CM.TXT   - charge, multiplicity',/,
     4       '            AFTER.TXT   - after coordinate options',/,
     5       '             CP2K.TXT   - CP2K input head',/,/,
     4       '        CP2KAFTER.TXT   - after coordinate options',/,
     6       '     Output: FILE.INP   - Gaussian Input',/,
     7       '                         (CP2K directory for CP2K)',/)
c
      open(7,file='FILE.X')
      read(7,80)s80
80    format(a80)
      read(7,*)nat
      n=3*nat
      allocate(r(3,nat),s(n,n),w(n),steps(n),ig(nat),iqlist(n),nl(n))
      do 1 i=1,nat
1     read(7,*)ig(i),(r(j,i),j=1,3)
      close(7)

      write(6,*)nat,' atoms'

      inquire(file='CP2K.TXT',exist=lcp2k)
      if(lcp2k)CALL system('mkdir CP2K')

      inquire(file='CASTEP.TXT',exist=lcastep)
      if(lcastep)then
       CALL system('mkdir CASTEP')
       call setcage(o11i,o12i,o13i,o22i,o23i,o33i)
      endif

      inquire(file='PMZ.PAR',exist=lex)
      lvstep=.false.
c     separate (differently-named) checkpoint files:
      lschk=.false.
      nochk=.false.
      step=0.05d0
      qdim=.false.
      dist=0.0d0
      idiff=2
      ic=0
      nparallel=0
      nmd=0
      ldalton=.false.
      only_double=.false.
      sdalbas='cc-pVDZ'
      if(lex)then
       open(7,file='PMZ.PAR')
72     read(7,700)key
700    format(a4)
       if(key.eq.'NO00')read(7,*)no_zero
       if(key.eq.'STEP')read(7,*)step
       if(key.eq.'QDIM')read(7,*)qdim
       if(key.eq.'PARA')read(7,*)nparallel
       if(key.eq.'DONL')read(7,*)only_double !Do only double steps d/dq_i^2, no cross terms like d/dq_i.dq_j
       if(key.eq.'DIST')read(7,*)dist
       if(key.eq.'IDIF')read(7,*)idiff
       if(key.eq.'DALB')read(7,*)sdalbas
       if(key.eq.'DALT')read(7,*)ldalton
       if(key.eq.'VSTE')read(7,*)lvstep
       if(key.eq.'LSCH')read(7,*)lschk
       if(key.eq.'NOCH')read(7,*)nochk
       if(key(1:2).eq.'IC')read(7,*)ic
       if(key(1:3).eq.'NMD')then
        read(7,*)nmd
        do 4 i=1,nmd
4       read(7,*)iqlist(i)
       endif
       if(key(1:3).eq.'END')goto 71
       backspace 7
       read(7,700)k2
       if(k2.eq.key)then
        write(6,700)k2
        call report('unknown option')
       endif
       goto 72
71     close(7)
      else
c      'Step (A) and degree (1,2,3) of the differentiation:
       step=0.05d0
       idiff=2
c      'cartesian or normal mode (0,1):'
       ic=0 
      endif
      
      if(no_zero)then
         io_start=1
      else
         io_start=0
      end if
      
c     idiff = 1  one-step diff
c             2  two-step diff (forward and back)
c             3  one-step in two coordinates
c             4  two-steps in two coordinates
c             5  <not used>
c             6  each atom two-step in x,y and z
c             8  two-step, ordered by coordinate (x1+,x1-,y1+,y1-, etc)
c                in relaxed normal modes
c            88  as 8, for two coordinates
c            22  two-step, ordered by coordinate (x1+,x1-,y1+,y1-, etc)
c            44  two-step, two coordinate, ordered by coordinate
c                (x+y+,x+y-,x-y+,x-y-,...), DIST keyword applied
c           100  each atom coordinate/mode differentiated N-times,
c                N is defined in file N.LST, single differentiation
c      	    200  each atom coordinate/mode differentiated N-times, 
c                N is defined in file N.LST, double differentiation
c
c     QDIM - use the dimensionless normal mode coordinates
c
c     VSTEP-variable step for each mode read as third column in N.LST,
c           for idif ? 100.
c

      if(ldalton)CALL system('mkdir DALTON')

      if(ic.eq.1)then
       write(6,*)'Normal mode differentiation'
       call loads(s,nat,w,nm)
       if(nmd.eq.0)then
        nmd=nm
        do 5 i=1,nm
5       iqlist(i)=i
       endif
       write(6,*)nmd,' modes to differentiate'
       do 6 i=1,nmd
6      write(6,600)i,iqlist(i),w(iqlist(i))
600    format(i4,i4,f8.1,' cm-1')
      else
       write(6,*)'Cartesian mode differentiation'
      endif

c     0000000000000000000000000000000000000000000000000000000
      if(idiff.eq.0)then
       io=0
       open(8,file='FILE.INP')
       call wgtxt(lschk,io)

       if(lcp2k)then
        write(ts,500)io
        istart=fistart(ts)
        CALL system('mkdir CP2K/'//ts(istart:len(ts)))
        open(88,file='CP2K/'//ts(istart:len(ts))//'/c.inp')
        open(89,file='CP2K/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt2
       endif

       write(ts,500)io
       istart=fistart(ts)

       if(lcastep)then
        CALL system('mkdir CASTEP/'//ts(istart:len(ts)))
        open(98,file='CASTEP/'//ts(istart:len(ts))//'/c.cell')
        open(99,file='CASTEP/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt3
       endif

       if(ldalton)then
        open(79,file='DALTON/'//ts(istart:len(ts))//'.INP')
        write(79,791)
791     format('BASIS')
        write(79,7911)(sdalbas(ii:ii),ii=fistart(sdalbas),len(sdalbas))
7911    format(80a1)
       endif

       write(8,*)
       write(8,80)s80
       write(6,3003)
       write(8,3003)
       if(lcp2k)write(89,3003)
       if(lcastep)write(99,3003)
       if(ldalton)then
        write(79,3003)
        write(79,3003)
        write(tsn,'(i4)')nat
        tso='AtomTypes='
     1//tsn(fistart(tsn):len(tsn))//' Angstroms Nosymmetry'
        write(79,7911)(tso(ii:ii),ii=fistart(tso),len(tso))
       endif
       write(8,*)
       call wcmtxt

       do 13333 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       
       if(lcp2k)write(88,8011)sy(ig(ii)),x,y,z
       
       if(lcastep)then
        xf=x*o11i+y*o12i+z*o13i
        yf=       y*o22i+z*o23i
        zf=              z*o33i      
        write(98,8012)sy(ig(ii)),xf,yf,zf
       endif

       if(ldalton)then
        if(ig(ii).lt.10)then
         write(79,3019)real(ig(ii)),ig(ii),x,y,z
3019     format('Charge=',f3.1,' Atoms=1',/,i1,3f12.6)
        else
         write(79,3020)real(ig(ii)),ig(ii),x,y,z
3020     format('Charge=',f4.1,' Atoms=1',/,i2,3f12.6)
        endif
       endif

13333  write(8,801)ig(ii),x,y,z

       call wafter
       write(8,*)
       if(lcastep)then
        call wafter3
        close(98)
        close(99)
       endif
       if(lcp2k)then
        call wafter2
        close(88)
        close(89)
       endif
       if(ldalton)close(79)
       close(8)
       write(6,*)1,' structure'
       stop
      endif
c     00000000000000000000000000000000000000000000000000000000000
c     1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 
      if(idiff.eq.1.or.idiff.eq.2)then
       ist=1
       i=1
       ix=0
       im=0

       if(ic.eq.1)then
        ndiff=nmd*idiff
       else
        ndiff=n*idiff
       endif

       do 2 io=0,ndiff
       if(io==0)then
         s20='equilibrium geometry'
       else
         if(io>ndiff/2)then
            sgn='-'
         else
            sgn='+'
         end if
         if(ic==0)then
            if(io>ndiff/2)then
               at=(io-ndiff/2-1)/3+1
            else
               at=(io-1)/3+1
            end if
            s20=i2str(at)
            if(mod(io,3)==1)then
               axis='x'
            else if(mod(io,3)==2)then
               axis='y'
            else
               axis='z'
            end if
            write(coord,'(a1,a,a1,a1)')'r',TR(s20),axis,sgn
         else
            write(coord,'(a1,a,a1)')'q',TR(s20),sgn
         end if
         s20=coord
       end if
       if(no_zero.and.io==0)cycle
       if(io.gt.0)then
        if(ic.eq.0)then
         ix=ix+1
         if(ix.gt.3)then
          ix=1
          i=i+1
          if(i.gt.nat)then
           i=1
           ist=ist+1
          endif
         endif
        else
         im=im+1
         if(im.gt.nmd)then
          ist=ist+1
          im=1
         endif
        endif
       endif

       open(8,file='FILE.INP')
       if(io.gt.io_start)write(8,802)
802    format('--link1--')
       call wgtxt(lschk,io)

       write(ts,500)io
500    format(i80)
       istart=fistart(ts)

       if(lcp2k)then
        CALL system('mkdir CP2K/'//ts(istart:len(ts)))
        open(88,file='CP2K/'//ts(istart:len(ts))//'/c.inp')
        open(89,file='CP2K/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt2
       endif

       if(lcastep)then
        CALL system('mkdir CASTEP/'//ts(istart:len(ts)))
        open(98,file='CASTEP/'//ts(istart:len(ts))//'/c.cell')
        open(99,file='CASTEP/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt3
       endif

       if(ldalton)then
        open(79,file='DALTON/'//ts(istart:len(ts))//'.INP')
        write(79,791)
        write(79,7911)(sdalbas(ii:ii),ii=fistart(sdalbas),len(sdalbas))
       endif

       if(io.gt.0 .and. .not. lschk .and. .not. nochk)write(8,803)
803    format('guess=checkpoint')
       write(8,*)
       if(io>0)write(8,'(a20)')s20

       if(io.eq.0)then
        write(6,3003)
3003    format('equilibrium geometry')
        write(8,3003)
        if(lcp2k)write(89,3003)
        if(lcastep)write(99,3003)
        if(ldalton)then
         write(79,3003)
         write(79,3003)
         write(tsn,'(i4)')nat
         tso='AtomTypes='
     1//tsn(fistart(tsn):len(tsn))//' Angstroms Nosymmetry'
         write(79,7911)(tso(ii:ii),ii=fistart(tso),len(tso))
        endif
       else
        if(ic.eq.0)then
         write(6,3004)round(ist),i,ix
3004     format(a7,', atom ',i4,', coord ',i2)
         if(lcp2k)write(89,3004)round(ist),i,ix
         if(lcastep)write(99,3004)round(ist),i,ix
         if(ldalton)then
          write(79,3004)round(ist),i,ix
          write(79,3004)round(ist),i,ix
          write(tsn,'(i4)')nat
          tso='AtomTypes='
     1//tsn(fistart(tsn):len(tsn))//' Angstroms Nosymmetry'
          write(79,7911)(tso(ii:ii),ii=fistart(tso),len(tso))
         endif
        else
         wcm=w(iqlist(im))
         write(6,3005)round(ist),im,iqlist(im),wcm
3005     format(a7,', mode ',i4,', (',i3,')',f8.1,' cm-1')
         write(8,3005)round(ist),im,iqlist(im),wcm
         if(lcp2k)write(89,3005)round(ist),im,iqlist(im),wcm
         if(lcastep)write(99,3005)round(ist),im,iqlist(im),wcm
         if(ldalton)then
          write(79,3005)round(ist),im,iqlist(im),wcm
          write(79,3005)round(ist),im,iqlist(im),wcm
          write(tsn,'(i4)')nat
          tso='AtomTypes='
     1//tsn(fistart(tsn):len(tsn))//' Angstroms Nosymmetry'
          write(79,7911)(tso(ii:ii),ii=fistart(tso),len(tso))
         endif
        endif
       endif

       write(8,*)
       call wcmtxt
       drmax=0.0d0

       do 3 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(io.gt.0)then
        if(ic.eq.0)then
         if(ii.eq.i)then
          if(ix.eq.1)x=x+(3.0d0-2.0d0*dble(ist))*step
          if(ix.eq.2)y=y+(3.0d0-2.0d0*dble(ist))*step
          if(ix.eq.3)z=z+(3.0d0-2.0d0*dble(ist))*step
         endif
        else
         if(qdim)then
C         dqs=step*dsqrt(1.0d3/wcm)
c         Conversion of s-matrix from amu^-1/2 to au^-1/2, dimensionless normal modes definition in atomic units and then conversion to angstrom
c         https://doi.org/10.1021/acs.jctc.3c01223
C           dqs=(1d0/sqrt(amu_2_au))*1d0/sqrt(2d0*pi*c_au*wcm*cm_2_au)
C      1     *step/ang_2_au
C           dqs=1d0/sqrt(2d0*pi*c_au*wcm*cm_2_au)
C      1     *step/ang_2_au
          dqs=sqrt(hbar/(2d0*pi*c*wcm*100))
     1     *step/sqrt(amu_2_kg)*1d10
         else
          dqs=step
         endif
         dx=s((ii-1)*3+1,iqlist(im))*dqs*(3.0d0-2.0d0*dble(ist))
         dy=s((ii-1)*3+2,iqlist(im))*dqs*(3.0d0-2.0d0*dble(ist))
         dz=s((ii-1)*3+3,iqlist(im))*dqs*(3.0d0-2.0d0*dble(ist))
         x=x+dx
         y=y+dy
         z=z+dz
         if(dabs(dx).gt.drmax)drmax=dabs(dx)
         if(dabs(dy).gt.drmax)drmax=dabs(dy)
         if(dabs(dz).gt.drmax)drmax=dabs(dz)
        endif
       endif
       if(lcp2k)then
        write(88,8011)sy(ig(ii)),x,y,z
8011    format(1x,a2,1x,3f20.10)
       endif
       if(lcastep)then
        xf=x*o11i+y*o12i+z*o13i
        yf=       y*o22i+z*o23i
        zf=              z*o33i      
        write(98,8012)sy(ig(ii)),xf,yf,zf
8012    format(1x,a2,3f21.16)
       endif
       if(ldalton)write(79,3019)real(ig(ii)),ig(ii),x,y,z
3      write(8,801)ig(ii),x,y,z

       if(ic.eq.1)write(6,6003)drmax
6003   format('drmax = ',f12.6,' A')
801    format(i3,3f12.6)
       call wafter
       write(8,*)
       if(lcastep)then
        call wafter3
        close(98)
        close(99)
       endif
       if(lcp2k)then
        call wafter2
        close(88)
        close(89)
       endif
       if(ldalton)close(79)
2      continue
       close(8)
       write(6,*)io,' structures'
       stop
      endif
c     1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 

c     3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 
      if(idiff.eq.3.or.idiff.eq.4)then
       write(6,*)' True double differentiation'
       if(ic.eq.1)then
        ndiff=((idiff-2)  *nmd)**2
       else
        ndiff=((idiff-2)*n)**2
       endif
       ist=1
       jst=1
       i=1
       j=1
       ix=0
       jx=1
       im=0
       jm=1
       do 7 io=0,ndiff
       if(no_zero.and.io==0)cycle
       if(io.gt.0)then
        if(ic.eq.0)then
         ix=ix+1
         if(ix.gt.3)then
          ix=1
          i=i+1
          if(i.gt.nat)then
           i=1
           jx=jx+1
           if(jx.gt.3)then
            jx=1
            j=j+1
            if(j.gt.nat)then
             j=1
             ist=ist+1
             if(ist.gt.2)then
              ist=1
              jst=jst+1
             endif
            endif
           endif
          endif
         endif
        else
         im=im+1
         if(im.gt.nmd)then
          im=1
          jm=jm+1
          if(jm.gt.nmd)then
           jm=1
           ist=ist+1
           if(ist.gt.2)then
            ist=1
            jst=jst+1
           endif
          endif
         endif
        endif
       endif
       if(only_double.and.im/=jm.and.io.gt.0)cycle
       if(only_double .and. ist/=jst)cycle
       open(8,file='FILE.INP')
       if(io.gt.io_start)write(8,802)
       call wgtxt(lschk,io)
       if(io.gt.0 .and. .not. nochk)write(8,803)
       write(8,*)
       write(8,80)s80
       if(io.eq.0)then
        write(6,3003)
        write(8,3003)
       else
        if(ic.eq.0)then
         write(6,3414)round(ist),i,ix,round(jst),j,jx
         write(8,3414)round(ist),i,ix,round(jst),j,jx
3414     format(2(a7,', atom',i5,', coord ',i2,1x))
        else
         wcm=w(iqlist(im))
         ucm=w(iqlist(jm))
         write(6,3415)round(ist),im,iqlist(im),wcm,round(jst),
     1   jm,iqlist(jm),ucm
         write(8,3415)round(ist),im,iqlist(im),wcm,round(jst),
     1   jm,iqlist(jm),ucm
3415     format(2(a7,', mode ',i4,', (',i3,')',f8.1,' cm-1 '))
        endif
       endif
       write(8,*)
       call wcmtxt
       drmax=0.0d0
       do 8 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(io.gt.0)then
        if(ic.eq.0)then
         if(ii.eq.i)then
          if(ix.eq.1)x=x+(3.0d0-2.0d0*dble(ist))*step
          if(ix.eq.2)y=y+(3.0d0-2.0d0*dble(ist))*step
          if(ix.eq.3)z=z+(3.0d0-2.0d0*dble(ist))*step
         endif
         if(ii.eq.j)then
          if(jx.eq.1)x=x+(3.0d0-2.0d0*dble(jst))*step
          if(jx.eq.2)y=y+(3.0d0-2.0d0*dble(jst))*step
          if(jx.eq.3)z=z+(3.0d0-2.0d0*dble(jst))*step
         endif
        else
         if(qdim)then
C         dqs=step*dsqrt(1.0d3/wcm)
c         Conversion of s-matrix from amu^-1/2 to au^-1/2, dimensionless normal modes definition in atomic units and then conversion to angstrom
c         https://doi.org/10.1021/acs.jctc.3c01223
C           dqs=(1d0/sqrt(amu_2_au))*1d0/sqrt(2d0*pi*c_au*wcm*cm_2_au)
C      1     *step/ang_2_au
          dqs=1d0/sqrt(2d0*pi*c_au*wcm*cm_2_au)
     1     *step/ang_2_au
         else
          dqs=step
         endif
         dx=s((ii-1)*3+1,iqlist(im))*dqs*(3.0d0-2.0d0*dble(ist))
     1     +s((ii-1)*3+1,iqlist(jm))*dqs*(3.0d0-2.0d0*dble(jst))
         dy=s((ii-1)*3+2,iqlist(im))*dqs*(3.0d0-2.0d0*dble(ist))
     1     +s((ii-1)*3+2,iqlist(jm))*dqs*(3.0d0-2.0d0*dble(jst))
         dz=s((ii-1)*3+3,iqlist(im))*dqs*(3.0d0-2.0d0*dble(ist))
     1     +s((ii-1)*3+3,iqlist(jm))*dqs*(3.0d0-2.0d0*dble(jst))
         x=x+dx
         y=y+dy
         z=z+dz
         if(dabs(dx).gt.drmax)drmax=dabs(dx)
         if(dabs(dy).gt.drmax)drmax=dabs(dy)
         if(dabs(dz).gt.drmax)drmax=dabs(dz)
        endif
       endif
8      write(8,801)ig(ii),x,y,z
       if(ic.eq.1)write(6,6003)drmax
       call wafter
       write(8,*)
7      continue
       close(8)
       write(6,*)io,' structures'
       stop
      endif
c     3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 3 4 

c     6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 
      if(idiff.eq.6)then
       write(6,*)' Two-atom quartic force field'
       if(ic.eq.1)call report('Cartesian coordinates allowed only')
       ndiff=36*nat
       if(nparallel.gt.1)then
        npart=ndiff/nparallel
       else
        npart=ndiff
       endif
       iparallel=1
       if(nparallel.gt.1)then
        write(s10,1010)iparallel
1010    format(i10)
        do 101 is=1,10
101     if(s10(is:is).ne.' ')goto 102
102     do 103 ie=10,1,-1
103     if(s10(ie:ie).ne.' ')goto 104
104     write(6,*)'FILE.36.'//s10(is:ie)//'.INP opened'
        open(8,file='FILE.36.'//s10(is:ie)//'.INP')
       else
        open(8,file='FILE.36.INP')
       endif
c      number of atoms
       i=1
c      coordinate
       ix=1
       jx=1
c      direction:
       ist=0
       jst=1
       do 9 io=0,ndiff
       if(io.gt.0)then
        if(io.gt.0)then
         ist=ist+1
         if(ist.gt.2)then
          ist=1
          jst=jst+1
          if(jst.gt.2)then
           jst=1
           ix=ix+1
           if(ix.gt.3)then
            ix=1
            jx=jx+1
            if(jx.gt.3)then
             jx=1
             i=i+1
            endif
           endif
          endif
         endif
        endif
       endif
       if(nparallel.gt.1)then
        if(mod(io,npart).eq.0.and.
     1   iparallel.lt.nparallel.and.io.gt.0)then
         write(6,*)'FILE.36.'//s10(is:ie)//'.INP closed'
         close(8)
         iparallel=iparallel+1
         write(s10,1010)iparallel
         do 201 is=1,10
201      if(s10(is:is).ne.' ')goto 202
202      do 203 ie=10,1,-1
203      if(s10(ie:ie).ne.' ')goto 204
204      open(8,file='FILE.36.'//s10(is:ie)//'.INP')
         write(6,*)'FILE.36.'//s10(is:ie)//'.INP opened'
         call wgtxt(lschk,io)
        else
         if(io.gt.0)write(8,802)
         call wgtxt(lschk,io)
         if(io.gt.0 .and. .not. nochk)write(8,803)
        endif
       else
        if(io.gt.0)write(8,802)
        call wgtxt(lschk,io)
        if(io.gt.0 .and. .not. nochk)write(8,803)
       endif
       write(8,*)
       write(8,80)s80
       if(io.eq.0)then
        write(6,3003)
        write(8,3003)
       else
        write(6,3016)io,i,xyz(ix),round(ist),xyz(jx),round(jst)
        write(8,3016)io,i,xyz(ix),round(ist),xyz(jx),round(jst)
3016    format(i5,1h:,' atom',i5,2(', ',a1,' ',a7))
       endif
       write(8,*)
       call wcmtxt
       do 10 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(io.gt.0)then
        if(ii.eq.i)then
         if(ix.eq.1)x=x+(3.0d0-2.0d0*dble(ist))*step
         if(ix.eq.2)y=y+(3.0d0-2.0d0*dble(ist))*step
         if(ix.eq.3)z=z+(3.0d0-2.0d0*dble(ist))*step
         if(jx.eq.1)x=x+(3.0d0-2.0d0*dble(jst))*step
         if(jx.eq.2)y=y+(3.0d0-2.0d0*dble(jst))*step
         if(jx.eq.3)z=z+(3.0d0-2.0d0*dble(jst))*step
        endif
       endif
10     write(8,801)ig(ii),x,y,z
       call wafter
       write(8,*)
9      continue
       close(8)
       write(6,*)io,' structures'
       stop
      endif
c     6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 6 
c     8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 
      if(idiff.eq.8.or.idiff.eq.88)then
c      relaxed normal mode coordinates
       open(91,file='g_opt')
       CALL system('mkdir diff')
       nn=0
c      zero geometry:
       write(ts,500)nn
       istart=fistart(ts)
       CALL system('mkdir diff/'//ts(istart:len(ts)))
       open(8,file='diff/'//ts(istart:len(ts))//'/FREE.INP')
       call wgtxtfres(8,'FREE.TXT')
       close(8)
       open(81,file='diff/'//ts(istart:len(ts))//'/G98.INP')
       call wgtxtfres(81,'G98.TXT')
       write(81,*)
       write(81,80)s80
       write(81,3003)
       write(6,3003)
       write(81,*)
       call wcmtxtn(81)

       if(lcp2k)then
        open(88,file='diff/'//ts(istart:len(ts))//'/c.inp')
        call wgtxt2
       endif

       if(lcastep)then
        open(88,file='diff/'//ts(istart:len(ts))//'/c.cell')
        call wgtxt3
       endif

       do 1004 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(lcp2k)then
        write(88,8011)sy(ig(ii)),x,y,z
       endif
       if(lcastep)then
        xf=x*o11i+y*o12i+z*o13i
        yf=       y*o22i+z*o23i
        zf=              z*o33i      
        write(98,8012)sy(ig(ii)),xf,yf,zf
       endif
1004   write(81,801)ig(ii),x,y,z

       call waftern(81)
       write(81,*)
       close(81)

       open(89,file='diff/'//ts(istart:len(ts))//'/t.inp')
       write(89,3003)
       close(89)

       if(lcastep)then
        call wafter3
        close(98)
        close(99)
       endif

       if(lcp2k)then
        call wafter2
        close(88)
        close(89)
       endif

       write(91,*)'cd diff/'//ts(istart:len(ts))
       write(91,*)'g09 G98.INP G98.OUT'
       write(91,*)'g09 FREE.INP FREE.OUT'
       write(91,*)'cd ../..'
       write(91,*)'cat diff/'//ts(istart:len(ts))//'/FREE.OUT > G.OUT'

c      coordinate cycles:
       if(ic.eq.1)then
        ndiff=nmd
       else
        ndiff=n
       endif
c      first coordinate
       do  3006 io=1,ndiff
       if(ic.eq.0)then
        i=(io+2)/3
        ix=io-3*(i-1)
       endif
c      first coordinate back,0,forward
       do  3006 ib=-1,1
       ist=ib+2

       if(idiff.eq.8)then
        ndiff1=1
        ndiff2=1
        jb1=0
        jb2=0
       else
        ndiff1=io
        ndiff2=ndiff
        jb1=-1
        jb2=1
       endif
c      second coordinate
       do  3006 jo=ndiff1,ndiff2
       dd=0.0d0
       if(ic.eq.0)then
        j=(jo+2)/3
        dd=
     1  dsqrt((r(1,i)-r(1,j))**2+(r(2,i)-r(2,j))**2+(r(3,i)-r(3,j))**2)
        jx=jo-3*(j-1)
       endif
c      second coordinate back,0,forward
       do  3006 jb=jb1,jb2
c      skip for too distant atoms:
       if(ic.eq.0.and.dist.ne.0.0d0.and.dd.gt.dd)goto 33333
c      skip equilibrium:
       if(ib.eq.0.and.idiff.eq.8)goto 33333
c      eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
       jst=jb+2
       nn=nn+1
       write(ts,500)nn
       istart=fistart(ts)
       CALL system('mkdir diff/'//ts(istart:len(ts)))
       open(8,file='diff/'//ts(istart:len(ts))//'/FREE.INP')
       call wgtxtfres(8,'FREE.TXT')
       close(8)
       open(81,file='diff/'//ts(istart:len(ts))//'/G98.INP')
       call wgtxtfres(81,'G98.TXT')
       open(89,file='diff/'//ts(istart:len(ts))//'/t.inp')
       if(lcp2k)then
        open(88,file='diff/'//ts(istart:len(ts))//'/c.inp')
        call wgtxt2
       endif

       if(lcastep)then
        open(98,file='diff/'//ts(istart:len(ts))//'/c.cell')
        call wgtxt3
       endif

       write(81,*)
       write(81,80)s80

       if(idiff.eq.8)then
        if(ic.eq.0)then
         write(6,3304)i,xyz(ix),r3(ist)
         write(81,3304)i,xyz(ix),r3(ist)
         write(89,3304)i,xyz(ix),r3(ist)
        else
         wcm=w(iqlist(io))
         write(6,3305)io,r3(ist),iqlist(io),wcm
         write(81,3305)io,r3(ist),iqlist(io),wcm
         write(89,3305)io,r3(ist),iqlist(io),wcm
        endif
       else
        if(ic.eq.0)then
         write(6,3314)i,xyz(ix),r3(ist),j,xyz(jx),r3(jst)
         write(81,3314)i,xyz(ix),r3(ist),j,xyz(jx),r3(jst)
         write(89,3314)i,xyz(ix),r3(ist),j,xyz(jx),
     1   r3(jst)
        else
         wcm=w(iqlist(io))
         ucm=w(iqlist(jo))
         write(6,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
         write(81,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
         write(89,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
        endif
       endif

       write(81,*)
       call wcmtxtn(81)
       drmax=0.0d0
       do 3009 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(ic.eq.0)then
        if(ii.eq.i)then
         if(ix.eq.1)x=x+dble(ib)*step
         if(ix.eq.2)y=y+dble(ib)*step
         if(ix.eq.3)z=z+dble(ib)*step
        endif
        if(ii.eq.j)then
         if(jx.eq.1)x=x+dble(jb)*step
         if(jx.eq.2)y=y+dble(jb)*step
         if(jx.eq.3)z=z+dble(jb)*step
        endif
       else
        dx=s((ii-1)*3+1,iqlist(io))*step*dble(ib)
     1    +s((ii-1)*3+1,iqlist(jo))*step*dble(jb)
        dy=s((ii-1)*3+2,iqlist(io))*step*dble(ib)
     1    +s((ii-1)*3+2,iqlist(jo))*step*dble(jb)
        dz=s((ii-1)*3+3,iqlist(io))*step*dble(ib)
     1    +s((ii-1)*3+3,iqlist(jo))*step*dble(jb)
        x=x+dx
        y=y+dy
        z=z+dz
        if(dabs(dx).gt.drmax)drmax=dabs(dx)
        if(dabs(dy).gt.drmax)drmax=dabs(dy)
        if(dabs(dz).gt.drmax)drmax=dabs(dz)
       endif
       if(lcp2k)write(88,8011)sy(ig(ii)),x,y,z
       if(lcastep)then
        xf=x*o11i+y*o12i+z*o13i
        yf=       y*o22i+z*o23i
        zf=              z*o33i      
        write(98,8012)sy(ig(ii)),xf,yf,zf
       endif
3009   write(81,801)ig(ii),x,y,z
       if(ic.eq.1)write(6,6003)drmax
       call waftern(81)
       write(81,*)
       if(lcastep)then
        call wafter3
        close(98)
       endif
       if(lcp2k)then
        call wafter2
        close(88)
       endif
c      eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
       close(81)
       close(89)

       open(89,file='diff/'//ts(istart:len(ts))//'/Q.OPT')
       call wgtxtfres(89,'Q.TXT')
       write(89,8902)iqlist(io)
8902   format('modes fixed',/,'1',/,i4,/,
     1 'hessian update',/,'f',/,
     1 's-matrix',/,'F.INP',/,'cartesian update',/,'f',/,'end')
       close(89)

       write(91,*)'cd diff/'//ts(istart:len(ts))
       write(91,*)'cp ../../F.INP .'
       write(91,*)'cp ../../FTRY.INP .'
       write(91,*)'runopt'
       write(91,*)'g09 FREE.INP FREE.OUT'
       write(91,*)'cd ../..'
       write(91,*)'cat diff/'//ts(istart:len(ts))//'/FREE.OUT >> G.OUT'

33333  continue
3006   continue

       close(91)

       write(6,*)nn+1,' structures'
       stop
      endif
c     8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 
c     22 44 22 44 22 44 22 44 22 44 22 44 22 44 22 44 22 44 22 44 
      if(idiff.eq.22.or.idiff.eq.44)then
       nn=0
       open(8,file='FILE.INP')
c      zero geometry:
       nn=nn+1
       call wgtxt(lschk,0)
       write(8,*)
       write(8,80)s80
       write(8,3003)
       write(6,3003)
       write(8,*)
       call wcmtxt

       if(lcp2k)then
        write(ts,500)io
        istart=fistart(ts)
        CALL system('mkdir CP2K/'//ts(istart:len(ts)))
        open(88,file='CP2K/'//ts(istart:len(ts))//'/c.inp')
        open(89,file='CP2K/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt2
        write(89,3003)
       endif

       if(lcastep)then
        write(ts,500)io
        istart=fistart(ts)
        CALL system('mkdir CASTEP/'//ts(istart:len(ts)))
        open(98,file='CASTEP/'//ts(istart:len(ts))//'/c.cell')
        open(99,file='CASTEP/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt3
        write(99,3003)
       endif

       do 1003 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(lcp2k)then
        write(88,8011)sy(ig(ii)),x,y,z
       endif
       if(lcastep)then
        xf=x*o11i+y*o12i+z*o13i
        yf=       y*o22i+z*o23i
        zf=              z*o33i      
        write(98,8012)sy(ig(ii)),xf,yf,zf
       endif
1003   write(8,801)ig(ii),x,y,z

       call wafter
       write(8,*)

       if(lcastep)then
        call wafter3
        close(98)
        close(99)
       endif

       if(lcp2k)then
        call wafter2
        close(88)
        close(89)
       endif

c      coordinate cycles:
       if(ic.eq.1)then
        ndiff=nmd
       else
        ndiff=n
       endif
c      first coordinate
       do  3002 io=1,ndiff
       if(ic.eq.0)then
        i=(io+2)/3
        ix=io-3*(i-1)
       endif
c      first coordinate back,0,forward
       do  3002 ib=-1,1
       ist=ib+2

       if(idiff.eq.22)then
        ndiff1=1
        ndiff2=1
        jb1=0
        jb2=0
       else
        ndiff1=io
        ndiff2=ndiff
        jb1=-1
        jb2=1
       endif
c      second coordinate
       do  3002 jo=ndiff1,ndiff2
       dd=0.0d0
       if(ic.eq.0)then
        j=(jo+2)/3
        dd=
     1  dsqrt((r(1,i)-r(1,j))**2+(r(2,i)-r(2,j))**2+(r(3,i)-r(3,j))**2)
        jx=jo-3*(j-1)
       endif
c      second coordinate back,0,forward
       do  3002 jb=jb1,jb2
c      skip for too distant atoms:
       if(ic.eq.0.and.dist.ne.0.0d0.and.dd.gt.dd)goto 3333
c      skip equilibrium:
       if(ib.eq.0.and.idiff.eq.22)goto 3333
c      eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
       jst=jb+2
       nn=nn+1
       write(8,802)
       write(ts,500)nn
       istart=fistart(ts)
       call wgtxt(lschk,jo)
       if(lcp2k)then
        CALL system('mkdir CP2K/'//ts(istart:len(ts)))
        open(88,file='CP2K/'//ts(istart:len(ts))//'/c.inp')
        open(89,file='CP2K/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt2
       endif

       if(lcastep)then
        CALL system('mkdir CASTEP/'//ts(istart:len(ts)))
        open(98,file='CASTEP/'//ts(istart:len(ts))//'/c.cell')
        open(99,file='CASTEP/'//ts(istart:len(ts))//'/t.inp')
        call wgtxt3
       endif

       if(.not.nochk)write(8,803)
       write(8,*)
       write(8,80)s80

       if(idiff.eq.22)then
        if(ic.eq.0)then
         write(6,3304)i,xyz(ix),r3(ist)
         write(8,3304)i,xyz(ix),r3(ist)
         if(lcp2k)write(89,3304)i,xyz(ix),r3(ist)
         if(lcastep)write(99,3304)i,xyz(ix),r3(ist)
3304     format(i4,2a1)
        else
         wcm=w(iqlist(io))
         write(6,3305)io,r3(ist),iqlist(io),wcm
         write(8,3305)io,r3(ist),iqlist(io),wcm
         if(lcp2k)write(89,3305)io,r3(ist),iqlist(io),wcm
         if(lcastep)write(99,3305)io,r3(ist),iqlist(io),wcm
3305     format(i4,a1,' (',i3,')',f8.1,' cm-1')
        endif
       else
        if(ic.eq.0)then
         write(6,3314)i,xyz(ix),r3(ist),j,xyz(jx),r3(jst)
         write(8,3314)i,xyz(ix),r3(ist),j,xyz(jx),r3(jst)
         if(lcp2k)write(89,3314)i,xyz(ix),r3(ist),j,xyz(jx),
     1   r3(jst)
         if(lcastep)write(99,3314)i,xyz(ix),r3(ist),j,xyz(jx),
     1   r3(jst)
3314     format(i5,2a1,' ',i5,2a1)
        else
         wcm=w(iqlist(io))
         ucm=w(iqlist(jo))
         write(6,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
         write(8,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
         if(lcp2k)write(89,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
         if(lcastep)write(99,3315)io,r3(ist),jo,r3(jst),
     1   iqlist(io),wcm,iqlist(jo),ucm
3315     format(2(i4,a1),' (',2(i3,',',f8.1,' cm-1'),')')
        endif
       endif

       write(8,*)
       call wcmtxt
       drmax=0.0d0
       do 3008 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(ic.eq.0)then
        if(ii.eq.i)then
         if(ix.eq.1)x=x+dble(ib)*step
         if(ix.eq.2)y=y+dble(ib)*step
         if(ix.eq.3)z=z+dble(ib)*step
        endif
        if(ii.eq.j)then
         if(jx.eq.1)x=x+dble(jb)*step
         if(jx.eq.2)y=y+dble(jb)*step
         if(jx.eq.3)z=z+dble(jb)*step
        endif
       else
        dx=s((ii-1)*3+1,iqlist(io))*step*dble(ib)
     1    +s((ii-1)*3+1,iqlist(jo))*step*dble(jb)
        dy=s((ii-1)*3+2,iqlist(io))*step*dble(ib)
     1    +s((ii-1)*3+2,iqlist(jo))*step*dble(jb)
        dz=s((ii-1)*3+3,iqlist(io))*step*dble(ib)
     1    +s((ii-1)*3+3,iqlist(jo))*step*dble(jb)
        x=x+dx
        y=y+dy
        z=z+dz
        if(dabs(dx).gt.drmax)drmax=dabs(dx)
        if(dabs(dy).gt.drmax)drmax=dabs(dy)
        if(dabs(dz).gt.drmax)drmax=dabs(dz)
       endif
       if(lcp2k)write(88,8011)sy(ig(ii)),x,y,z
       if(lcastep)then
        xf=x*o11i+y*o12i+z*o13i
        yf=       y*o22i+z*o23i
        zf=              z*o33i      
        write(98,8012)sy(ig(ii)),xf,yf,zf
       endif
3008   write(8,801)ig(ii),x,y,z
       if(ic.eq.1)write(6,6003)drmax
       call wafter
       write(8,*)
       if(lcastep)then
        call wafter3
        close(98)
        close(99)
       endif
       if(lcp2k)then
        call wafter2
        close(88)
        close(89)
       endif
c      eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
3333   continue
3002   continue

       close(8)
       write(6,*)nn,' structures'
       stop
      endif
c     22 44 22 44 22 44 22 44 22 44 22 44 22 44 22 44 22 44 22 44 

c     100 100 100 100 100 100 100 100 100 100 100 100 100 100 
      if(idiff.ge.100)then
       if(lvstep)
     1 write(6,*)'Differentiation steps predefined in N.LST'
       ndiff=0
       do 11 i=1,n
       ndiff=ndiff+2
       steps(i)=step
11     nl(i)=3
       open(33,file='N.LST')
       nn=0
12     if(lvstep)then
        read(33,*,end=333,err=333)i,j,steps(i)
       else
        read(33,*,end=333,err=333)i,j
       endif
       if(i.gt.n)call report('coordinate too high')
       if(mod(j+1,2).ne.0)call report('odd values allowed only')
       nl(i)=j
       ndiff=ndiff-3+nl(i)
       nn=nn+1
       goto 12
333    close(33)
       write(6,*)nn,' individual differentiations read from N.LST'
       if(qdim.and.ic.eq.1)write(6,*)'variable steps introduced'
c
c      zero point
       open(8,file='FILE.INP')
       call wgtxt(lschk,0)
       write(8,*)
       write(8,80)s80
       write(8,3003)
       write(8,*)
       call wcmtxt
       do 13 ii=1,nat
13     write(8,801)ig(ii),(r(ix,ii),ix=1,3)
       call wafter
       write(8,*)
c
c      Loop over coordinates:
       io=0
       step1=step
       do 14 ico1=1,n
c      skip modes that are not required in PMZ.PAR:
       if(ic.eq.1)then
        do 15 im=1,nmd
15      if(iqlist(im).eq.ico1)goto 151
        goto 14
151     continue
        step1=steps(ico1)
        if(qdim)step1=step1*dsqrt(1.0d3/w(ico1))
       endif
c      Loop over diff. points:
       do 141 ip1=1,nl(ico1)
       dd1=step1*dble((ip1-1-(nl(ico1)-1)/2))
c
c      Second loop over coordinates:
       do 214 ico2=ico1,n
       if(ico2.eq.ico1.and.idiff.eq.200)goto 214
       if(ico2.ne.ico1.and.idiff.eq.100)goto 214
c      skip modes that are not required:
       if(ic.eq.1)then
        do 215 im=1,nmd
215     if(iqlist(im).eq.ico2)goto 251
        goto 214
251     continue
        step2=steps(ico2)
        if(qdim)step2=step2*dsqrt(1.0d3/w(ico2))
       endif
c      Loop over diff. points:
       do 241 ip2=1,nl(ico2)
c      skip if only one-differentiation required:
       if(idiff.eq.100.and.ip2.ne.nl(ico2)/2+1)goto 241
       dd2=step2*dble((ip2-1-(nl(ico2)-1)/2))
c      skip if this produces zero geometry:
       if(ip2.eq.nl(ico2)/2+1.and.ip1.eq.nl(ico1)/2+1)goto 241

       io=io+1
       write(8,802)
       call wgtxt(lschk,io)
       if(.not.nochk)write(8,803)
       write(8,*)
       write(8,80)s80
c
       if(ic.eq.0)then
c       Cartesian:
        ia1=(ico1+2)/3
        ix=ico1-3*ia1+3
        write(6,3017)io,ia1,xyz(ix),ip1,dd1
        write(8,3017)io,ia1,xyz(ix),ip1,dd1
3017    format(i4,', atom',i4,1x,a1,', point',i3,f10.4)
        if(idiff.eq.200)then
         ia2=(ico2+2)/3
         iy=ico2-3*ia2+3
         write(6,3017)io,ia2,xyz(iy),ip2,dd2
         write(8,3017)io,ia2,xyz(iy),ip2,dd2
         write(8,30171)dd1,dd2
30171    format(' dd1 dd2',2f10.4)
        endif
       else
c       Normal Mode:
        wcm=w(ico1)
        write(6,3018)io,ico1,wcm,ip1,dd1
        write(8,3018)io,ico1,wcm,ip1,dd1
3018    format(i4,', mode',i4,1x,f8.1,' cm-1, point',i3,f10.4)
        if(idiff.eq.200)then
         write(6,3018)io,ico2,w(ico2),ip2,dd2
         write(8,3018)io,ico2,w(ico2),ip2,dd2
         write(8,30171)dd1,dd2
        endif
       endif
       write(8,*)
       call wcmtxt
       drmax=0.0d0
       do 16 ii=1,nat
       x=r(1,ii)
       y=r(2,ii)
       z=r(3,ii)
       if(ic.eq.0)then
c       Cartesian:
        if(ii.eq.ia1)then
         if(ix.eq.1)x=x+dd1
         if(ix.eq.2)y=y+dd1
         if(ix.eq.3)z=z+dd1
        endif
        if(idiff.eq.200)then
         if(ii.eq.ia2)then
          if(iy.eq.1)x=x+dd2
          if(iy.eq.2)y=y+dd2
          if(iy.eq.3)z=z+dd2
         endif
        endif
       else
c       Normal Mode:
        dx=s((ii-1)*3+1,ico1)*dd1+s((ii-1)*3+1,ico2)*dd2
        dy=s((ii-1)*3+2,ico1)*dd1+s((ii-1)*3+2,ico2)*dd2
        dz=s((ii-1)*3+3,ico1)*dd1+s((ii-1)*3+3,ico2)*dd2
        x=x+dx
        y=y+dy
        z=z+dz
        if(dabs(dx).gt.drmax)drmax=dabs(dx)
        if(dabs(dy).gt.drmax)drmax=dabs(dy)
        if(dabs(dz).gt.drmax)drmax=dabs(dz)
       endif
16     write(8,801)ig(ii),x,y,z
       if(ic.eq.1)write(6,6003)drmax
       call wafter
       write(8,*)
241    continue
214    continue

141    continue
14     continue
       close(8)
       write(6,*)io,' structures'
       stop
      endif
c     100 100 100 100 100 100 100 100 100 100 100 100 100 100 

      write(6,*)idiff
      write(6,*)'Unknown differentiation type requested'
      stop

      end
c     ============================================================
      subroutine loads(s,nat,w,nm)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION s(3*nat,3*nat),w(*)
      character*80 filename
      filename='F.INP'
      open(10,file=filename)
      read(10,*,end=9999,err=9999)nm
      do 2 i=1,nat+1
2     read(10,*)
      do 3 ia=1,nat
      do 3 im=1,nm
3     read(10,*)(s((ia-1)*3+ix,im),ix=1,2),(s((ia-1)*3+ix,im),ix=1,3)
      read(10,*)nm
      read(10,*)(w(im),im=1,nm)
      close(10)
      call msg('S matrix loaded from '//filename//'.')
      write(6,*)nm,' modes'
      return
9999  call bye(filename//'not found')
      end
c     ============================================================
      subroutine msg(ts)
      character*(*) ts
c
      do 1 i=min(80,len(ts)),1,-1
1     if(ts(i:i).ne.' ')goto 2
2     write(6,*)ts(1:i)
      return
      end
c ======================================================================
      subroutine bye(ts)
      character*(*) ts
      write(6,6001)
      write(6,'(a)')ts
      write(6,6001)
6001  format(80(1h*))
      write(6,6000)
6000  format('program stopped')
      close(6)
      stop
      end
c ======================================================================
      subroutine wgtxt2
      character*80 t80
      open(71,file='CP2K.TXT')
221   read(71,80,end=388,err=388)t80
      write(88,80)t80
80    format(a80)
      goto 221
388   close(71)
      return
      end
c ======================================================================
      subroutine wgtxt3
      character*80 t80
      open(71,file='CASTEP.TXT')
221   read(71,80,end=388,err=388)t80
      write(98,80)t80
80    format(a80)
      goto 221
388   close(71)
      return
      end
c ======================================================================
      subroutine wgtxtfres(io,s)
      character*(*) s
      character*80 t80
      integer*4 io
      open(71,file=s)
221   read(71,80,end=388,err=388)t80
      write(io,80)t80
80    format(a80)
      goto 221
388   close(71)
      return
      end
c ======================================================================
      subroutine wgtxt(lschk,io)
      implicit none
      integer*4 is,io,i
      character*80 t80
      character*10 tn
      logical lex,lschk
      write(tn,40)io+1
40    format(i10)
      do 1 is=1,len(tn)
1     if(tn(is:is).ne.' ')goto 2
2     inquire(file='G.TXT',exist=lex)
      if(lex)then
       open(71,file='G.TXT')
221    read(71,80,end=388,err=388)t80
       if(lschk.and.t80(1:4).eq.'%chk')then
        do 3 i=1,4
3       write(8,844)t80(i:i)
844     format(a1,$)
        write(8,844)'='
        do 4 i=is,len(tn)
4       write(8,844)tn(i:)
        write(8,844)'.'
        write(8,844)'c'
        write(8,844)'h'
        write(8,844)'k'
        write(8,*)
       else
        write(8,80)t80
       endif
80     format(a80)
       goto 221
388    close(71)
      else
       if(lschk)then
        write(8,*)'%chk='//tn(is:len(tn))//'.chk'
       else
        write(8,*)'%chk=xx.chk'
       endif
       write(8,8001)
8001   format('%mem=24000000',/,
     1 '#b3LYP/6-31++G** freq nosymm')
      endif
      return
      end
c ======================================================================
      subroutine wcmtxtn(io)
      character*80 t80
      logical lex
      integer*4 io
      inquire(file='CM.TXT',exist=lex)
      if(lex)then
       open(71,file='CM.TXT')
       read(71,80)t80
       close(71)
       write(io,80)t80
80     format(a80)
      else
       write(io,*)'0 1'
      endif
      return
      end
c ======================================================================
      subroutine wcmtxt
      character*80 t80
      logical lex
      inquire(file='CM.TXT',exist=lex)
      if(lex)then
       open(71,file='CM.TXT')
       read(71,80)t80
       close(71)
       write(8,80)t80
80     format(a80)
      else
       write(8,*)'0 1'
      endif
      return
      end
c ======================================================================
      subroutine wafter2
      character*80 t80
      open(71,file='CP2KAFTER.TXT')
8911  read(71,80,end=892,err=892)t80
80    format(a80)
      write(88,80)t80
      goto 8911
892   close(71)
      return
      end
c ======================================================================
      subroutine wafter3
      character*150 t80
      open(71,file='CASTEPAFTER.TXT')
8911  read(71,80,end=892,err=892)t80
80    format(a150)
      write(98,80)t80
      goto 8911
892   close(71)
      return
      end
c ======================================================================
      subroutine waftern(io)
      integer*4 io
      character*200 t200
      logical lex
      inquire(file='AFTER.TXT',exist=lex)
      if(lex)then
       write(io,*)
       open(71,file='AFTER.TXT')
8911   read(71,80,end=892,err=892)t200
80     format(a200)
       write(io,80)t200
       goto 8911
892    close(71)
      endif
      return
      end
c ======================================================================
      subroutine wafter
      character*200 t80
      logical lex
      inquire(file='AFTER.TXT',exist=lex)
      if(lex)then
       write(8,*)
       open(71,file='AFTER.TXT')
8911   read(71,80,end=892,err=892)t80
80     format(a200)
       write(8,80)t80
       goto 8911
892    close(71)
      endif
      return
      end
c ======================================================================
      subroutine report(s)
      character*(*) s
      write(6,'(a)')s
      stop
      end
c ======================================================================
      subroutine setcage(o11i,o12i,o13i,o22i,o23i,o33i)
      implicit none
      integer*4 i,j
      real*8 o11i,o12i,o13i,o22i,o23i,o33i,uv(3,3),a,b,c,ab,ac,bc,
     1alpha,beta,gamma,pi,ca,cb,cg,sg,v
c    ,o11,o12,o13,o22,o23,o33,sb
      character*80 s80
      
      pi=4.0d0*atan(1.0d0)
      open(88,file='CASTEP.TXT')
1     read(88,80,end=99,err=99)s80
80    format(a80)
      if(s80(1:19).eq.'%BLOCK LATTICE_CART')then
       do 201 i=1,3
201    read(88,*)(uv(i,j),j=1,3)
       write(6,*)'unit cell read'
       a=dsqrt(uv(1,1)**2+uv(1,2)**2+uv(1,3)**2)
       b=dsqrt(uv(2,1)**2+uv(2,2)**2+uv(2,3)**2)
       c=dsqrt(uv(3,1)**2+uv(3,2)**2+uv(3,3)**2)
       ab=uv(1,1)*uv(2,1)+uv(1,2)*uv(2,2)+uv(1,3)*uv(2,3)
       ac=uv(1,1)*uv(3,1)+uv(1,2)*uv(3,2)+uv(1,3)*uv(3,3)
       bc=uv(2,1)*uv(3,1)+uv(2,2)*uv(3,2)+uv(2,3)*uv(3,3)
       alpha=180.0d0/pi*acos(bc/b/c)
       beta=180.0d0/pi*acos(ac/a/c)
       gamma=180.0d0/pi*acos(ab/b/a)
       write(6,600)alpha,a,beta,b,gamma,c
600    format(/,' Crystall cell parameters:',/,/,
     1          '   alpha: ',f10.2,'   a: ',f10.4,/,
     1          '    beta: ',f10.2,'   b: ',f10.4,/,
     1          '   gamma: ',f10.2,'   c: ',f10.4,/,/)
       ca=cos(pi*alpha/180.0d0)  
c      sa=sin(pi*alpha/180.0d0)
       cb=cos(pi*beta /180.0d0)  
c      sb=sin(pi*beta /180.0d0)
       cg=cos(pi*gamma/180.0d0)  
       sg=sin(pi*gamma/180.0d0)
c      forward trafo:
c      o11=a
c      o12=b*cg
c      o13=c*cb
c      o22=b*sg
c      o23=-c*sb*(cb*cg-ca)/(sb*sg)
c      o33=c*sb*dsqrt(1.0d0-((cb*cg-ca)/(sb*sg))**2)
c      inverse trafo:
       o11i=1.0d0/a
       o12i=-cg/a/sg
       v=dsqrt(1.0d0-ca**2-cb**2-cg**2+2.0d0*ca*cb*cg)
       o13i=(ca*cg-cb)/a/v/sg
       o22i=1.0d0/b/sg
       o23i=(cb*cg-ca)/b/v/sg
       o33i=sg/c/v
       close(88)
       return
      endif
      goto 1
99    call report('cell not found')
      end

      function fistart(ts)
      implicit none
      integer*4 fistart,istart
      character*(*) ts
      do 5031 istart=1,len(ts)
5031  if(ts(istart:istart).ne.' ')goto 5041
5041  fistart=istart
      return
      end
      
      function i2str(i)result(str)
         implicit none
         integer i
         character(20) str
         write(str,*)i
         str=adjustl(str)
      end function i2str