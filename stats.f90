!c------------------------------------------------------------
!c                                                           c
!c  subroutine output statis data                            c
!c  Etot_ppv : Total Energy = Epot+Eki                       c
!c  Epot_ppv : Potential Energy = e_ppv+e_ppp                c
!c  E_ppv    : Potential Mechanical Energy (classical PE)    c
!c             = e_bond+e_bend+e_tors+e_LJ                   c
!c  E_ppp    : Potential electronic energy  (quantum PE)     c
!c             = e_core+epppkjmol
!c  epppkjmol: SCF energy + CI energy
!c-----------------------------------------------------------c
!c                                                           c
!c   New subroutien for PPV code                             c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm.utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c

SUBROUTINE stats
  !subroutine stat_ppv(istep,timestep,eki,epppkjmol,e_ppv,ebetag
  !   $     ,ecore,fmoinv_ppv,fkj_ppv,ene_gap,ene_ground,osc_ppv,u_conf
  !   $     ,osc_vec,vxa,vya,vza,dump_vel)
  
  USE commondata, ONLY: DP,istep,deltat,eki,eppp,ebetag,eqcff,eground,ecore,&
       a0,log_dumpvel,nnuc,log_cicalc
  USE classical, ONLY: vel
  USE qcfio, ONLY: stat_io,gap_io,vel_io,osc_io
  USE quantum, ONLY: osc,osc_vec,ciener,nci
  
  IMPLICIT NONE

  !local variables
  INTEGER :: i
  REAL(DP) :: etot,epot,efactor,eppp_in
  
  efactor=2625.5d0 !APW Hartree to kJ/mol

  etot    = eki  + eppp + ebetag + eqcff
  epot    =        eppp + ebetag + eqcff
  eppp_in =        eppp + ebetag
  
  ! AY comment 2011.01.03
  ! etot : total energy = kinetic energy (=eki) + potential energy (=epot)
  ! eki  : kinetic energy
  ! epot : potential energy = quantum part (=eppp_in) + classical part (=eqcff)
  ! eppp_in : potential energy of quantum part
  !           = SCF + CI energy (=eppp) + nuclear core-core repulsion (=ebetag)
  ! eppp, ebetag : see comments in modules.f90

  ! output to statis.dat. the unit of energies are kJ/mol
  write(stat_io,'(f10.6,2x,8g25.15)') istep*deltat,etot,eki,epot,eqcff,&
       eppp_in,eppp,ebetag,ecore
  
  if(log_cicalc) then
     do i=1,nci
        write(gap_io,'(f10.6,2x,g20.10,i3,g20.10,g20.10)') istep*deltat,&
             eground,i,ciener(i)*efactor,osc(i)*efactor/a0**2
        write(osc_io,'(f10.6,2x,i3,3(f20.6))') istep*deltat,i,&
             osc_vec(1,i)*efactor/a0,osc_vec(2,i)*efactor/a0,osc_vec(3,i)*efactor/a0
     end do
  end if
  
  if (log_dumpvel) then
     WRITE(vel_io,*) '#VELOCITY CONFIGS FOR STEP=',istep
     do i=1,nnuc
        WRITE(vel_io,'(3(f16.8,1x))') vel(1,i),vel(2,i),vel(3,i)
     end do
  end if
  
END SUBROUTINE stats
