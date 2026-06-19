program elvib_harmonic
   use iso_fortran_env
   use util
   use constants
   implicit none
   integer nat,n3,nq,iz,argc,ifx,i
   double precision,allocatable :: we(:),wg(:)
   double precision,allocatable :: A(:,:),B(:),C(:,:),D(:),E(:,:),J(:,:),K(:),SG(:,:),SE(:,:)
   character(80) s80
   double precision e00,e00_pseudoharm,we_avg,eg,ee
   
   
   call readdusch(nq,n3,A,B,C,D,E,J,K,SG,Se,wg,we,'DUSCH.OUT')
   deallocate(a,b,c,d,e,j,k,sg,se)
   open(77,file='ground/ENERGY')
   read(77,*)eg
   close(77)
   open(77,file='excited/ENERGY')
   read(77,*)ee
   close(77)
   e00=ee-eg
   
   we_avg=0d0
   we=we*cm_2_au
   do i = 1,nq
      we_avg=we_avg+we(i)
   end do
   we_avg=we_avg/nq
   
   e00_pseudoharm=e00+0.5d0*sum(we)+we_avg
   e00_pseudoharm=1d7/(e00_pseudoharm*au_2_cm)
   write(output_unit,*)e00_pseudoharm
   deallocate(wg,we)
   contains
   
   
   
end program elvib_harmonic