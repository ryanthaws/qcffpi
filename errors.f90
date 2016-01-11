!  errors.f90
!  28th November 2005
!  S. Reed 
!
! A simple subroutine to deal with errors in opening files.
! Gives an error message and stops the program

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! apapad:: Parallelization
! NOTE: These routines are only called by the master
! so no need to wrap the write statements in an if block
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

subroutine fileerror(ierr,filename)
  USE qcfio , ONLY: gen_io 
  !use fileunits
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! apapad:: Parallelization
  !USE paralmod, ONLY: comm, mpierr
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! Get the error number and filename from the calling subroutine
  integer, intent(in)      :: ierr
  character(*), intent(in) :: filename

  ! Give a specific message if we know what the error number means
  ! if not, give a message including the error number and stop.


  ! File already exists but was required to be new
  if (ierr == 10) then
     write(gen_io,*) "Fatal error: The file ",trim(filename)," file already exists."

    ! File doesn't exist but should
  else if (ierr == 29) then
     write(gen_io,*) "Fatal error: The file ",trim(filename)," does not exist so I can't open it!" 

     ! Another error whose meaning we don't know
  else 
     write(gen_io,*) "Fatal error, ",ierr,": Problem opening the file ",trim(filename)
  end if
  
  !call MPI_Abort(comm, 1, mpierr)
  stop

end subroutine fileerror





subroutine filereaderror(ierr,errlocation,message)
  !use fileunits
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! apapad:: Parallelization
  !USE paralmod, ONLY: comm, mpierr
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  USE qcfio, ONLY: gen_io
  
  implicit none
  integer, intent(in) :: ierr
  character(*),intent(in) :: errlocation
  character(*), intent(in) :: message
  if (ierr == -1) then
     write(gen_io,*) "Fatal error: Reached the end of file whilst reading from",trim(errlocation),"."
     !call MPI_Abort(comm, 1, mpierr)
     stop
  else     
     write(gen_io,*) "Fatal error number ",ierr," in ",trim(errlocation)
  end if

  write(gen_io,*) "  ",trim(message)
  !call MPI_Abort(comm, 1, mpierr)
  stop
end subroutine filereaderror
