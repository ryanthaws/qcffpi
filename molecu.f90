!-----------------------------------------------------------------!
!                                                                 !
! this subroutine calls all the (classical) potential energy      !
! functions cleaned up by mbedard.hearn on 24.oct.07 & 17.apr.08  !
! cleaned up again by adam willard June2010                       !
!    input variables:                                             !
!                                                                 ! 
! ihf: 1 : compute f matrix,bond orders and no c.i. calculation   !
! .... 0 : the same & solve the scf equation and print the input  !
! ....-1 : the same & solve the scf equation and print the output !
!                                                                 !
! isf: 0 : no scf                                                 !
! ...  1 : see above                                              !
!                                                                 !
! pos     : coordinates  -  units of Angstrom                      !
! deriv     : derivatives  -  units of kcal/(mol-Angstrom)           !
!          (grad U, not -grad U... see interface_qcffpi.f)        !
! eng_tot      : energy  -  units of kcal/mol                           !
!-----------------------------------------------------------------!
!                                                           !
!   Revised by Adam P. Willard 2010                         !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!

!subroutine molecu(ihf,isf,x_,d_,e,iblo_ppv,wb_ppv,wnb_ppv,wt_ppv
!  $     ,wp_ppv,ilink_ppv,wlin_ppv,iac2_ppv,nc_ppv,icg_ppv,istep)
SUBROUTINE MOLECU(eng_tot)

  USE commondata, ONLY: nnuc,DP,idsx,istep
  USE classical, ONLY:deriv,dderiv,b_sw,nb_sw,tp_sw,pos
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE

  !input variables
  REAL(DP), INTENT(OUT) :: eng_tot

  INTEGER :: i,j,k,l,i2,j2
  DOUBLE PRECISION :: eng_bond,eng_nonb,eng_theta,eng_phi,eng_dip,efactor
  INTEGER :: count1,count2,clock_rate,clock_max

  efactor  = 4.1868 !APW kcal/mol to kJ/mol
  idsx     =1
  eng_tot  =0.d0
  eng_bond =0.d0
  eng_nonb =0.d0
  eng_theta=0.d0
  eng_phi  =0.d0
  eng_dip  =0.d0
  
  deriv    =0.0d0
  dderiv   =0.0d0
  !call system_clock ( count1, clock_rate, clock_max )
  if( b_sw.eq.1) call het_bondp(eng_bond)
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'time1',(count2-count1)/real(clock_rate)
  !count1=count2
  if(nb_sw.eq.1) call het_nonbon(eng_nonb)
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'time2',(count2-count1)/real(clock_rate)
  !count1=count2
  call het_thetap(eng_theta)
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'time3',(count2-count1)/real(clock_rate)
  !count1=count2
  if(tp_sw.eq.1) call het_phip(eng_phi)
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'time4',(count2-count1)/real(clock_rate)


  call het_dipolex(eng_dip) !JJBL
  
  

  
  eng_tot = eng_bond + eng_nonb + eng_theta + eng_phi + eng_dip
  
  !print *,'shit1',eng_tot,eng_bond,eng_nonb,eng_theta,eng_phi,eng_dip

  write(gen_io,*) 'classical energies (kcal/mol) '
  write(gen_io,*) '----------------------------- '
  write(gen_io,*) 'bonding    :',eng_bond
  write(gen_io,*) 'non-bonding:',eng_nonb
  write(gen_io,*) 'theta      :',eng_theta
  write(gen_io,*) 'phi        :',eng_phi
  write(gen_io,*) 'dipole(dipx)',eng_dip,'or (in Hartree Units)',eng_dip/627.509391d0
  write(gen_io,*) 'total      :',eng_tot
  
  eng_tot = eng_tot*efactor
  
END SUBROUTINE MOLECU
