!------------------------------------------------------------
!                                                           !
!  Main Program                                             !
!  NonAdiabatic MD Simulation for PPV/FULLERENE             !
!  Code by Lobaugh/Sterpone/Bedard-Hearn                    !
!-----------------------------------------------------------!
!                                                           !
!   Author Adam P. Willard 2010                             !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!



!!$* =====================================================================
!!$* []  DECK FULL_FULLERENE_MAIN --  PROGRAM FOR NONADIABATIC MD
!!$*     SIMULATION OF FULLERENE/PPV MIXTURES IN THE GAS PHASE
!!$*     2009 June - MBH added Girifalco Potential for fullerene spheres
!!$*                 that will (hopefully) prevent the OPV from wrapping
!!$*                 entirely around the C60, as happens with 1x(7)OPV+C60.
!!$*                 The Girifalco potential treates the excess C60s as
!!$*                 large single-atom spheres that interact with the
!!$*                 atomistic particles in an averaged way.
!!$* =====================================================================


PROGRAM main

  USE commondata, ONLY: nsteps,navg,log_exdyn,istep,istat,irestart,&
       eqcff,eppp,eki,ecore,ebetag,d_screen,nnuc
  USE classical, ONLY: fra,vel,pos,deriv
  USE quantum, ONLY: fexa,fexa_l
  USE qcfio, ONLY: gen_io

  IMPLICIT NONE

  LOGICAL :: twist
  INTEGER :: ierr,i
  INTEGER :: count1,count2,clock_rate,clock_max,count3,count4

  twist=.true.
  
  open (gen_io,file='run.out',status='new',iostat=ierr)
  if(ierr .ne. 0) then
     print *,"The file run.out already exists!"
     print *,"Please remove it and all other output files before trying again"
     !stop
  end if

  !-------------------------------------------------------------------

  write(gen_io,*) '--------------------------------------------------'
  write(gen_io,*) '      CIS/PPP/QCFF Dynamics     '
  write(gen_io,*) '--------------------------------------------------'
  write(gen_io,*) '---> Reading Input'
  call input
  write(gen_io,*) '---> Finished Reading Input'
  write(gen_io,*)"Dielectric screening: d_screen = ",d_screen

  !--- Set Initial Variables -----------------------------------------
  call setup_vars
  
  !-------------------------------------------------------------------
  ! Read QCFF PARAMETER. We set First step=.true.
  ! Inside QCFF first is set false for the dynamics
  ! Pre-loop calculation of classical PE and forces
  !-------------------------------------------------------------------
  istep=0
  
  call molecu(eqcff)


  !*  +---------------------------------------------------------------+
  !*  |  Pre-loop Quantum calculations here                           |
  !*  +---------------------------------------------------------------+
  
  if(twist) call twist_ppv(eppp,ebetag,ecore)
  write(gen_io,*) 'done with twist'
  
  !*  +---------------------------------------------------------------+
  !*  |  get forces before starting MD loop                           |
  !*  +---------------------------------------------------------------+

  
  call add_force_qcff
  
  !*
  !*  =============  START THE MD LOOP  ===============
  !* 
  call system_clock ( count1, clock_rate, clock_max )
  do istep=0,nsteps
     !*  +---------------------------------------------------------------+
     !*  |  do velocity verlet: position and half velocity calculation   |
     !*  +---------------------------------------------------------------+
     write(gen_io,*) '-----------------------------------------------'
     write(gen_io,*) ' MD Loop - step = ', istep
     write(gen_io,*) '-----------------------------------------------'
     call verlet_s1(eki)
     write(gen_io,*) '-----------------------------------------------'
     write(gen_io,*) 'Verlet first step done'
     write(gen_io,*) '-----------------------------------------------'

     if(mod(istep,istat).eq.0) then
        call stats
     end if
     
     if(mod(istep,istat).eq.0) then
        call output_struct
        call make_pdb
     end if

     !APW I don't know what prot is or where outcyc is?
     !if (prot) call outcyc 

     !APW now above
     !call new_step
     navg=navg+1
     fexa_l = fexa
     fexa = 0.0d0
     fra=0.0d0
     
     !c------------------------------------------------------------------+
     !c Internal forces (classical potential and forces for main MD loop)|
     !c------------------------------------------------------------------+

     !call qcffsol(eqcff,first)
     !call system_clock ( count3, clock_rate, clock_max )
     call molecu(eqcff)
     !call system_clock ( count4, clock_rate, clock_max )
     !print *,'molecu time',istep,(count4-count3)/real(clock_rate)
     
     write(gen_io,*) '---> Internal energy solute e = ',eqcff,'Kj/mol'
     
     !C       MBH SWAP ORDER OF TWIST AND ADD_FORCE CALLS March 2009
     !C       MAKES NO NUMERICAL DIFFERENCE FOR SHORT TEST RUNS (ground
     !C       and excited state runs of 1x(5)PPV)
     !*  +---------------------------------------------------------------+
     !*  |  Quantum calculations here for main MD loop                   |
     !*  +---------------------------------------------------------------+
     
     !call system_clock ( count3, clock_rate, clock_max )
     if(twist) call twist_ppv(eppp,ebetag,ecore)
     !call system_clock ( count4, clock_rate, clock_max )
     !print *,'twist time',istep,(count4-count3)/real(clock_rate)
     !*  +---------------------------------------------------------------+
     !*  |  add quantum (hellman-feynman) and classical (nuclear) forces |
     !*  +---------------------------------------------------------------+

     call add_force_qcff

     !*  +---------------------------------------------------------------+
     !*  |  do 2nd velocity verlet: velocity update                      |
     !*  +---------------------------------------------------------------+
     call verlet_s2
     write(gen_io,*) '-----------------------------------------------'
     write(gen_io,*) 'Verlet second step done'
     write(gen_io,*) '-----------------------------------------------'

     !*  +---------------------------------------------------------------+
     !*  |  do runge-kutta integration of the expansion coefficients     |  
     !*  |  and check for a transition                                   |  
     !*  +---------------------------------------------------------------+
     
     if(log_exdyn) then
        call rgkt
        call transition
        write(gen_io,*) 'transition analysis (Tully FSSH) done'
     end if
        
     if(istep.ne.0.and.mod(istep,irestart).eq.0) then
        write(gen_io,*) '**** WRITING RESTART FILE ****',istep
        call write_restart
        write(gen_io,*) '**** done writing restart file ****'
     end if
     
     call system_clock ( count2, clock_rate, clock_max )
     print *,'step',istep,(count2-count1)/real(clock_rate)
     count1=count2
     
  end do
  call cleanup

END PROGRAM main


