SUBROUTINE input
  
  USE commondata
  USE classical, ONLY: dfield
  USE quantum, ONLY: ntully,iexstate,nci,neig,osc,osc_vec
  USE qcfio, ONLY: gen_io
  USE ran_state, ONLY: ran_seed
  
  IMPLICIT NONE
  
  INTEGER :: ierr,i,j,k,nseeds
  INTEGER, ALLOCATABLE, DIMENSION(:) :: seeds
  CHARACTER(len=50) :: line !APW for reading in from a file
  INTEGER, DIMENSION(6) :: nt
  
  open(12,file='start.inpt',status='old',iostat=ierr)
  
  if (ierr .ne. 0) then
     call fileerror(ierr,"start.inpt")
     stop
  end if
  
  read(12,*,iostat=ierr)nsteps
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading nsteps")
  read(12,*,iostat=ierr)nnuc
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading nnuc")
  read(12,*,iostat=ierr)npi
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading npi")
  read(12,*,iostat=ierr)norb
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading norb")
  read(12,*,iostat=ierr)tmb
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading tmb")
  read(12,*,iostat=ierr) log_cicalc
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_cicalc")
  read(12,*,iostat=ierr)nci
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading nci")
  neig=nci
  !read(12,*,iostat=ierr)neig
  !if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading neig")
  read(12,*,iostat=ierr) log_exdyn
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_exdyn")
  read(12,*,iostat=ierr)ntully
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading ntully")
  read(12,*,iostat=ierr)iexstate
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading iexstate")
  read(12,*,iostat=ierr) log_restart
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_restart")  
  read(12,*,iostat=ierr)irestart
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading irestart")
  read(12,*,iostat=ierr) log_rescale_i
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_rescale_i")  
  read(12,*,iostat=ierr) log_scalevel
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_scalevel")
  read(12,*,iostat=ierr)iscale
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading iscale")
  read(12,*,iostat=ierr)ttol
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading ttol")  
  read(12,*,iostat=ierr)istat
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading istat")
  !read(12,*,iostat=ierr)iofreq
  !if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading iofreq")
  read(12,*,iostat=ierr)d_screen
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading d_screen")
  read(12,*,iostat=ierr)iseed
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading iseed")
  read(12,*,iostat=ierr) dfield
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading dfield")
  read(12,*,iostat=ierr)deltat
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading deltat")
  read(12,*,iostat=ierr) log_ppv
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_ppv")
  read(12,*,iostat=ierr) log_C60shift
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_C60shift")
  
  read(12,*,iostat=ierr)log_solv
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_solv")
  read(12,*,iostat=ierr)nsolve
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading nsolve")
  
  read(12,*,iostat=ierr)log_delre
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_delre")
  read(12,*,iostat=ierr) log_dumpvel
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_dumpvel")
  read(12,*,iostat=ierr) log_confine
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_confine")
  read(12,*,iostat=ierr)log_dderiv
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading log_dderiv")
  read(12,*,iostat=ierr)nexstate
  if(ierr .ne. 0) call filereaderror(ierr,"start.inpt","while reading nexstate")
    
  close(12)

  ! AY added
  ! The option of functions of alpha and gamma.
  ! If this option (log_addcos) is .T., alpha and gamma functions 
  ! depend on torsion angle. 
  ! (= original forms of PI/QCFF by Warshel & Karplus. JACS, 94, 5612 (1972))
  ! If it's .F., not depend on torsion angle.
  ! (= method of Sterpone & Roosky, JPC B 112,4983,(2008))
  log_addcos = .true.
!  log_addcos = .false.

  
  write(gen_io,*) "******************start.inpt**********************"
  write(gen_io,*) "nsteps = ",nsteps
  write(gen_io,*) "nnuc = ",nnuc
  write(gen_io,*) "npi = ",npi
  write(gen_io,*) "norb = ",norb
  write(gen_io,*) "tmb = ",tmb
  write(gen_io,*) "log_cicalc = ",log_cicalc
  write(gen_io,*) "nci = ",nci
  write(gen_io,*) "neig = ",neig
  write(gen_io,*) "log_exdyn = ",log_exdyn
  write(gen_io,*) "ntully = ",ntully  
  write(gen_io,*) "iexstate = ",iexstate
  write(gen_io,*) "log_restart = ",log_restart
  write(gen_io,*) "irestart = ",irestart
  write(gen_io,*) "log_rescale_i = ",log_rescale_i
  write(gen_io,*) "log_scalevel = ",log_scalevel
  write(gen_io,*) "iscale = ",iscale
  write(gen_io,*) "ttol = ",ttol
  write(gen_io,*) "istat = ",istat
  !write(gen_io,*) "iofreq = ",iofreq
  write(gen_io,*) "d_screen = ",d_screen
  write(gen_io,*) "iseed = ",iseed
  write(gen_io,*) "dfield = ",dfield
  write(gen_io,*) "deltat = ",deltat
  write(gen_io,*) "log_ppv = ",log_ppv
  
  write(gen_io,*) "log_solv = ",log_solv
  write(gen_io,*) "nsolve = ",nsolve

  write(gen_io,*) "log_delre = ", log_delre
  write(gen_io,*) "log_dumpvel = ",log_dumpvel
  write(gen_io,*) "confining = ",log_confine
  write(gen_io,*) "log_dderiv = ",log_dderiv
  write(gen_io,*) "nexstate = ", nexstate

  

  ! AY added
  write(gen_io,*) "log_addcos = ",log_addcos, "# (see input.f90)"
  if(log_addcos) then
     if(abs(d_screen-1.0d0).gt.1.0d-4) write(gen_io,*)"Warning: No support that d_screen .ne. 1.0"
     if(.not.log_ppv) write(gen_io,*)"Warning: log_addcos is not ON"
  endif
  
  write(gen_io,*) "************************************************"
  
  call ran_seed(sequence=iseed)

  !*  +---------------------------------------------------------------+
  !*  |  Sanity checks                                                |
  !*  +---------------------------------------------------------------+
  
  if ((.not.log_cicalc).and.(log_exdyn)) then
     write(*,*) 'You must calculate the CI states'
     write(*,*) 'if you are going to do excited state'
     write(*,*) 'dynamics.'
     STOP
  end if
  
  call npar(0)  !read in qcff potential parameters
   
END SUBROUTINE input
