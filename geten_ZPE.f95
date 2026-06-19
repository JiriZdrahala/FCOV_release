program geten
   use iso_fortran_env
   implicit none
   
   character(80) s80,s80_2
   character(10) enUnit
   integer eq_idx,argc,i,nd,bs
   double precision :: En_SCF,En_ZPE,ev,en_tot_zpe
   logical fex
   
   En_SCF=HUGE(1d0)
   En_ZPE=HUGE(1d0)
   
   argc=COMMAND_ARGUMENT_COUNT()
   if(argc==0)then
      write(output_unit,*)'Writes energy from Gaussian output file (1st arg.) to file "ENERGY"'
      write(output_unit,*)'Also corrects it with Zero-point correction from a frequency job file (2nd arg.)'
      write(output_unit,*)'USAGE: '
      write(output_unit,*)'geten_ZPE <.out Gaussian file> <.out Gaussian freq file> [<cm-1/eV/au/nm/THz/kJmol/kcalmol>]'
      call exit(2)
   end if
   
   call get_command_argument(1,s80)
   call get_command_argument(2,s80_2)
   
   enUnit='au'
   if(argc==3)then
      call get_command_argument(3,enUnit)
      call To_lower(enUnit)
   end if
   
   inquire(file=s80,exist=fex)
   if(.not.fex)then
      write(output_unit,*)'ERROR: File not found: '//s80
      call exit(1)
   end if
   open(77,file=s80,action='READ')
   
   inquire(file=s80_2,exist=fex)
   if(.not.fex)then
      write(output_unit,*)'ERROR: File not found: '//s80_2
      call exit(1)
   end if
   open(76,file=s80_2,action='READ')
   
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
      if(s80(2:30)=='Total Energy, E(TD-HF/TD-DFT)')then
         read(s80(33:48),*)En_SCF
         write(output_unit,*)'Read TD "Total Energy"'
      end if
      goto 1
2  close(77)
   
   open(76,file=s80_2)
11 read(76,'(A80)',end=22)s80
      if(s80(2:23)=='Zero-point correction=')then
         read(s80(45:58),*)En_ZPE
         goto 22
      end if
      goto 11
22 close(76)  
   
   En_SCF=convertEnergy(en_scf+En_ZPE,enUnit)
   open(78,file='ENERGY')
   write(78,'(G19.12)')En_SCF
   close(78)
   write(output_unit,*)'Written file ENERGY'
   
   call exit(0)
   contains
   
   !to au
   function convertEnergy(ener,toUnit)result(res)
      !TAKEN FROM https://www.weizmann.ac.il/oc/martin/tools/hartree.html
      double precision,parameter :: cmm1=219474.63,evv=27.211399,nm=45.5640,&
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