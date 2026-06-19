module orca_parser
    implicit none

   contains

   subroutine extract_last_geometry(orca_filename, x_filename,r)
        ! ---------------------------------------------------------------------
        ! Extracts the final geometry from an ORCA output file and writes it
        ! into the custom FILE.X format with atomic numbers and trailing zeroes.
        ! ---------------------------------------------------------------------
        character(len=*), intent(in) :: orca_filename
        character(len=*), intent(in) :: x_filename
        double precision,allocatable :: r(:,:)
        
        integer :: in_unit, out_unit, ierr, ierr2
        character(len=512) :: line
        
        ! Buffer sizing (10,000 atoms is more than sufficient for standard ORCA runs)
        integer, parameter :: MAX_ATOMS = 10000
        character(len=2) :: temp_sym(MAX_ATOMS), final_sym(MAX_ATOMS)
        real(8) :: temp_coords(3, MAX_ATOMS), final_coords(3, MAX_ATOMS)
        
        integer :: n_atoms_temp, n_atoms
        integer :: i, z_num
        
        n_atoms = 0
        
        open(newunit=in_unit, file=orca_filename, status='old', action='read')
        
        ! 1. Parse the entire file, overwriting the geometry buffers upon each match
        do
            read(in_unit, '(A)',end=20) line
            
            if (index(line, 'CARTESIAN COORDINATES (ANGSTROEM)') > 0) then
                ! Skip the dashes line
                read(in_unit, '(A)') line
                
                n_atoms_temp = 0
                do
                    read(in_unit, '(A)') line
                    if (len_trim(line) == 0) exit
                    
                    n_atoms_temp = n_atoms_temp + 1
                    if (n_atoms_temp > MAX_ATOMS) then
                        print *, "Error: Exceeded MAX_ATOMS limit of ", MAX_ATOMS
                        close(in_unit)
                        return
                    end if
                    
                    ! ORCA format: Symbol X Y Z
                    read(line, *) temp_sym(n_atoms_temp), &
                                                temp_coords(1, n_atoms_temp), &
                                                temp_coords(2, n_atoms_temp), &
                                                temp_coords(3, n_atoms_temp)
                    
                end do
                
                ! Update the final geometry buffer if read was successful
                if (n_atoms_temp > 0) then
                    n_atoms = n_atoms_temp
                    final_sym(1:n_atoms) = temp_sym(1:n_atoms)
                    final_coords(:, 1:n_atoms) = temp_coords(:, 1:n_atoms)
                end if
            end if
        end do
        
 20       close(in_unit)
        
        if (n_atoms == 0) then
            print *, "Error: No cartesian coordinates found in the ORCA file."
            return
        end if
        
        ! 2. Write the extracted geometry to the output file
        open(newunit=out_unit, file=x_filename, status='replace', action='write')
        
        write(out_unit, '(A)') ' Last Geometry Extracted from ORCA output'
        write(out_unit, '(I12)') n_atoms
        
        allocate(r(3,n_atoms))
        do i = 1, n_atoms
            z_num = get_atomic_number(final_sym(i))
            ! Format: atomic_number (I4), X Y Z (F15.8), followed by static padding
            write(out_unit, '(I4, 3F15.8, A)') z_num, &
                final_coords(1, i), final_coords(2, i), final_coords(3, i), &
                ' 0 0 0 0 0 0 0 0.0'
            r(:,i)=final_coords(:,i)
        end do
        
        close(out_unit)
        
    end subroutine extract_last_geometry


    function get_atomic_number(symbol) result(z)
        ! ---------------------------------------------------------------------
        ! Converts an atomic symbol string to its corresponding atomic number (Z).
        ! Supports the complete periodic table (Z = 1 to 118).
        ! ---------------------------------------------------------------------
        character(len=*), intent(in) :: symbol
        integer :: z
        character(len=2) :: u_sym
        integer :: i
        
        character(len=2), parameter :: elements(118) = [ &
            'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE', &
            'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA', &
            'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', &
            'GA', 'GE', 'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', &
            'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', &
            'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', &
            'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', &
            'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', &
            'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', &
            'PA', 'U ', 'NP', 'PU', 'AM', 'CM', 'BK', 'CF', 'ES', 'FM', &
            'MD', 'NO', 'LR', 'RF', 'DB', 'SG', 'BH', 'HS', 'MT', 'DS', &
            'RG', 'CN', 'NH', 'FL', 'MC', 'LV', 'TS', 'OG' ]
            
        ! Preprocess the string: convert to uppercase and left-align
        u_sym = '  '
        if (len_trim(symbol) == 1) then
            u_sym(1:1) = to_upper(symbol(1:1))
        else if (len_trim(symbol) >= 2) then
            u_sym(1:1) = to_upper(symbol(1:1))
            u_sym(2:2) = to_upper(symbol(2:2))
        end if
        
        z = 0
        do i = 1, 118
            if (u_sym == elements(i)) then
                z = i
                return
            end if
        end do
        
        ! Fallback if not found (will print as 0 in output, raising an obvious red flag)
        print *, "Warning: Unknown atomic symbol '", trim(symbol), "'"
    end function get_atomic_number


    function to_upper(c) result(uc)
        ! ---------------------------------------------------------------------
        ! Helper function to convert a single character to uppercase.
        ! ---------------------------------------------------------------------
        character(len=1), intent(in) :: c
        character(len=1) :: uc
        integer :: ic
        
        ic = iachar(c)
        if (ic >= iachar('a') .and. ic <= iachar('z')) then
            uc = achar(ic - iachar('a') + iachar('A'))
        else
            uc = c
        end if
    end function to_upper

end module orca_parser

program gar_orca
   use orca_parser
   implicit none
   character(500) s80
   double precision,allocatable :: r(:,:)
   
   call GET_COMMAND_ARGUMENT(1,s80)
   call extract_last_geometry(s80,'FILE.X',r)
   
contains


end program gar_orca