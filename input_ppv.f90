SUBROUTINE input
  
  USE qcfio, ONLY: gen_io
  USE ppv
  
  IMPLICIT NONE
  
  INTEGER :: ierr,i,j,k
  CHARACTER(len=50) :: line !APW for reading in from a file
  
  open(11,file='junction.dat',status='old',iostat=ierr)
  
  if (ierr .ne. 0) then
     call fileerror(ierr,"junction.dat")
     stop
  end if
  
  read(11,*,iostat=ierr) ntor_junc
  ALLOCATE(junction(6,ntor_junc))
  
  do i=1,ntor_junc
     read(11,*,iostat=ierr) ntor_junc(:,i)  
     if(ierr .ne. 0) call filereaderror(ierr,"param.inpt","while reading junction")
  end do
  
  close(11)

END SUBROUTINE input
