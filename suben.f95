program suben
   use iso_fortran_env
   implicit none
   
   character(500) file1,file2
   double precision :: en1,en2
   
   call GET_COMMAND_ARGUMENT(1,file1)
   call GET_COMMAND_ARGUMENT(2,file2)
   
   open(77,file=file1)
   read(77,*)en1
   close(77)
   
   open(77,file=file2)
   read(77,*)en2
   close(77)
   
   write(output_unit,*)en1-en2   
end program suben