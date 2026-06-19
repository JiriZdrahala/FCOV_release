program geten
   use iso_fortran_env
   use constants
   implicit none
   
   character(80) s80
   character(10) enUnit
   integer eq_idx,argc,i,nd,bs
   double precision :: En_SCF,En_TD,ev,en_tot_zpe
   logical fex
   
   En_SCF=HUGE(1d0)
   En_TD=0
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc==0)then
      write(output_unit,*)'Writes the molecular energy to file "ENERGY" from "SCF Done/EUMPx"'
      write(output_unit,*)'If there is a TD calculation as the last job, it will print the exc. state energy (E_SCF + E_TD) instead, from "This state for optimization..."'
      ! write(output_unit,*)'If there is a Freq (or TD-freq) calculation as the last job, it will print the energy with ZPC (E_SCF + E_ZPE (+ E_TD)) &
      ! instead, from "Sum of electronic and zero-point Energies"'
      write(output_unit,*)'USAGE: '
      write(output_unit,*)'geten <.out Gaussian file> [<cm-1/eV/au/nm/THz/kJmol/kcalmol>]'
      call exit(2)
   end if
   
   call get_command_argument(1,s80)
   
   enUnit='au'
   if(argc==2)then
      call get_command_argument(2,enUnit)
      call To_lower(enUnit)
   end if
   
   inquire(file=s80,exist=fex)
   if(.not.fex)then
      write(output_unit,*)'ERROR: File not found: '//s80
      call exit(1)
   end if
   open(77,file=s80,action='READ')
   
1  read(77,'(A80)',end=2)s80
      if(s80(2:9)=='SCF Done')then
         eq_idx = index(s80,'=')
         read(s80(eq_idx+2:eq_idx+2+19),*)En_SCF !is this losing precision, maybe not.
         !but it also checks if it's a valid number
         write(output_unit,*)'Read SCF Energy'
      end if
      if(s80(28:32)=='EUMP2')then
         eq_idx = index(s80,'=',.true.)
         read(s80(eq_idx+2:eq_idx+2+24),*)En_SCF
         write(output_unit,*)'Read MP2 Energy'
      end if
      if(s80(35:39)=='EUMP3')then
         eq_idx = index(s80,'=',.true.)
         read(s80(eq_idx+2:eq_idx+2+24),*)En_SCF
         write(output_unit,*)'Read MP3 Energy'
      end if
      if(s80(35:44)=='EUMP4(SDTQ)')then
         eq_idx = index(s80,'=',.true.)
         read(s80(eq_idx+1:eq_idx+1+19),*)En_SCF
         write(output_unit,*)'Read MP4(SDTQ) Energy, (I really do not know what that means btw)'
      end if
      if(s80(2:15)=='This state for')then
         !backspace 77
         !backspace 77
         bs=0
         do while(s80(2:15)/='Excited State ')
            backspace 77
            backspace 77
            bs=bs+1
            read(77,'(a80)')s80
         end do
         read(s80(51:58),*)ev
         ev=(1d7/ev)*cm_2_au
         ! do i=1,len(s80)-1
            ! if(s80(i:i)  .eq.':') read(s80(15:i-1),*)nd
            ! if(s80(i:i+1).eq.'eV')read(s80(i-10:i-2),*)ev
         ! end do
         ! En_TD=ev/(convertEnergy(1d0,'ev'))
         en_TD=ev
         write(output_unit,*)'Read TD energy of state: ',nd
         do i = 1,bs
            read(77,*)
         end do
      end if
      !I commented this because I deal with imaginary modes all the time
      !When there are no imag. modes in the ground state, but a large one in the excited state...
      !it makes the energy difference larger than it should be.
      ! if(s80(2:42)=='Sum of electronic and zero-point Energies')then
         ! eq_idx = index(s80,'=',.true.)
         ! read(s80(eq_idx+1:80),*)en_tot_zpe
         ! en_tot_zpe=convertEnergy(en_tot_zpe,enUnit)
         ! write(output_unit,*)'Read Freq energy'
         ! open(78,file='ENERGY')
         ! write(78,'(G19.12)')en_tot_zpe
         ! close(78)
         ! write(output_unit,*)'Written file ENERGY'
         ! return !terminate program
      ! end if
      goto 1
2  close(77)
   
   En_SCF=convertEnergy(en_scf+en_TD,enUnit)
   open(78,file='ENERGY')
   write(78,'(G19.12)')En_SCF
   close(78)
   write(output_unit,*)'Written file ENERGY'
   
   call exit(0)
   contains
   
   !to au
   function convertEnergy(ener,toUnit)result(res)
      !TAKEN FROM https://www.weizmann.ac.il/oc/martin/tools/hartree.html
      double precision,parameter :: cmm1=219474.63,evv=27.21138386,nm=45.5640,&
                                    thz=6579.683879634054
      double precision,parameter :: kjmol=2625.5002,kcalmol=627.5096080305927
      double precision ener,res
      character(*) toUnit
      
      select case(toUnit)
         case('cm-1')
            res=ener*cmm1
         case('ev')
            res=ener*evv
         case('au')
            res=ener
         case('nm')
            res=nm/ener
         case('thz')
            res=ener*thz
         case('kcalmol')
            res=ener*kcalmol
         case('kjmol')
            res=ener*kjmol
         case default
            write(output_unit,*)'ERROR: Unknown unit: '//toUnit
            call exit(3)
      end select
   end function convertEnergy
   
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
   
end program geten