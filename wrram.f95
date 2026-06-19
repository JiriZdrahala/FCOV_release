#define TR(arg) trim(adjustl(arg))
module wrram
   use constants
   implicit none
   
   integer,parameter :: ni0=13,ns0=19
   character(9),parameter :: rroa_exp(ns0) = &
     ['ICP_0    ','ICP_x_90 ','ICP_z_90 ','ICP_*_90 ','ICP_u_90 ','ICP_180  ',&
      'SCP_0    ','SCP_x_90 ','SCP_z_90 ','SCP_*_90 ','SCP_u_90 ','SCP_180  ',&
      'DCPI_0   ','DCPI_90  ','DCPI_180 ',&
      'DCPII_0  ','DCPII_90 ','DCPII_180',&
      'SPECIAL  ']
      
   public :: ni0,ns0,rroa_exp
      
   contains
   
   
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
   subroutine wrram3(w,temp,spectrum_temp,ei,ef,THRC,nt,EXCNM,wrin,wrax,npx,lglg,fwhh,ltab, &
            inv,alr,ali,gtr,gti,gtcr,gtci,atr,ati,atcr,atci,sr,bf,lusea,luseg,ldo,&
            tabOutputs,str_v1,str_v3,alr_fcht,ali_fcht)
!     EXCNM ... excitation frequency in nm
!     bf .. Boltzmann factor
      implicit none
      integer*4 ni,nf,nt,npx, &
      ia,b,id,e,ii,is
      integer tabOutputs(ns0)
!     ns0 : number of experimental setups
!     ni0 : number of invariants
      real*8 ef,ei,fwhh,THRC,CM,AMU,BOHR,ECM,EXCNM,wrax,wrin,w,clight, &
      YDY,YDX,sr(npx,2,ns0),gpisvejc,EXCA,roa1,ram1,tr,ti,a(2),bf, &
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
      alai_fcht(3,3),aaar_fcht,aaai_fcht,buf(3,3)
      double precision wR,temp,Kcon,e0,pi
      character(*) str_v1,str_v3
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
      gpisvejc=(AMU*BOHR**5)*1.0d4*2.0d0*4.0d0*atan(1.0d0)!/EXCA
      clight=137.03599d0
      pi=4*atan(1d0)
      e0=1d0/(4d0*pi)
      
      ! buf=Gtr
      ! Gtr=Gti
      ! Gti=buf
      
      ! buf=Gtcr
      ! Gtcr=Gtci
      ! Gtci=buf
      
      !Kcon=((pi/e0)**2)*((wr/(2*pi*cc_AU))**4)/90
      wR=w-ef
      !Kcon=((pi/e0)**2)*((wR/(2*pi*cc_AU))**4)/90
      Kcon=wR**4*gpisvejc/90
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

      !Okay, why is clight here? I cannot find it anywhere in literature
      !c usually appears in magnetic properties
      if(lusea)then
!      betasA2: (w/2)Im(i als_ab e_agd As*g,db)
       call abab(tr,ti,alsr,alsi,easr,easi)
       inv(7)=0.5d0*w*tr!/clight
!      betaaA2: (w/2)Im(i[ala_ab[e_agd Aa*g,db+ e_abg Aa*d,gd)
       call abab(tr,ti,alar,alai,eaar,eaai)
       inv(8)=tr
       call abab(tr,ti,alar,alai,eapar,eapai)
       inv(8)=0.5d0*w*(inv(8)+tr)!/clight
!      betasAc2: (w/2)Im(i als_ab e_agd Acs*g,db)
       call abab(tr,ti,alsr,alsi,eacsr,eacsi)
       inv(12)=0.5d0*wR*tr!/clight
!      betaaAc2:
       call abab(tr,ti,alar,alai,eacar,eacai)
       inv(13)=tr
       call abab(tr,ti,alar,alai,eapcar,eapcai)
       inv(13)=0.5d0*wR*(inv(13)+tr)!/clight
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
        if(ltab)then
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
   
   subroutine ap3(is,ns0,ECM,a,s,wrin,wrax,npx,fwhh,lglg,temp,temp_spr)
   integer*4 ns0,npx,i,is
   real*8 ECM,wrin,wrax,fwhh,w,dw,a(2),s(npx,2,ns0),t,temp,fac
   logical lglg,temp_spr
   dw=(wrax-wrin)/(npx-1)
   w=wrin-dw
   if(temp_spr)then
      ! fac=1d0/((ECM*cm_2_au)*(1d0-exp(-ECM*100*h*cc/(kb*temp))))
      fac=1d0/((1d0-exp(-ECM*100*h*cc/(kb*temp))))
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
   
end module wrram
