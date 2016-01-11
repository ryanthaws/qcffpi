!c------------------------------------------------------------
!c                                                           c
!c  subroutine output strcutral data                         c
!c  Torsion  : cosine of torsion around junction             c
!c  Bond     : bond distance of PI atoms                     c
!c-----------------------------------------------------------c
!c                                                           c
!c   New Subroutine for PPV code                             c
!c   Author Fabio Sterpone 2005                              c
!c   sterpone@cm.utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c

SUBROUTINE output_struct
  !subroutine anal_ppv(istep,pibond,nap_ppv)!,vx_pi,vy_pi,vz_pi)
  
  USE commondata, ONLY: istep,deltat
  USE qcfio, ONLY: tor_io,bond_io
  USE PPV, ONLY: torsion,pibond,ntor_junc
  USE quantum, ONLY: ibeta_pairs
  
  IMPLICIT NONE

  INTEGER :: i
  
  !c--- Dump Torsional angles -------------------------------------
  !c      if(istep.eq.1) write(79,'(a50)')
  !c     $     '    ------- Torsion angles around junction -------'
  
  !APW HERE need to make torsion a global variables
  write(tor_io,'(200(f12.8))') istep*deltat,(torsion(i),i=1,ntor_junc)
  !write(79,7901) istep*timestep,(torsion(i),i=1,ntor_jun)
  
  !c--- Dump Bond Distance ----------------------------------------
  !c      if(istep.eq.1) write(80,'(a50)')
  !c     $     '     ------- Bond Distance for PI atoms  -------  '
  write(bond_io,'(200(f12.8))') istep*deltat,(pibond(i),i=1,ibeta_pairs)
  
  !c--- Dump Velocity PI Atoms ------------------------------------
  !c     if(dump_vel) write(83,*) istep*timestep,(vx_pi(i),vy_pi(i),vz_pi(i
  !c    $     ),i=1,nap_ppv)
  
END SUBROUTINE output_struct

