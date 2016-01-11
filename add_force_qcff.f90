!c------------------------------------------------------------
!c                                                           c
!c  subroutine adding quantum and classical forces           c
!c  and putting them in the correct units                    c
!c                                                           c
!c-----------------------------------------------------------c
!c                                                           c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm.utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c


SUBROUTINE add_force_qcff
!subroutine ADD_FORCE_PPV(fqcffx,fqcffy,fqcffz,fxa,fya,fza,fkj_ppv
  !     $     ,unitl_ppv,molcs)
  
  USE commondata, ONLY:DP,nnuc,navg
  USE classical, ONLY: deriv,fra
  
  IMPLICIT NONE

  !local variables
  REAL(DP) :: forcefac
  
  INTEGER :: i
  
  !c--- Executive Stantements -------------------------------------
  
  !* unitl has units A/internal length
  !* fkj_ppv has units kJ/internal energy
  
  !* multiply by 4.1868 J / cal, so converts:  kcal -> kJ
  !* divide by fkj_ppv converts:               kJ   -> int. energy
  !* multiply by unitl_ppv converts:           1/A  -> 1/int. leng.
  
  !* all this means that fqcffx(i) comes in with units of kcal / (mol-A)
  !* and fxa(i) enters and leaves with interal units of force (int. E/ int L)
  
  !APW the units of deriv come in as kcal/mol/angstrom
  forcefac=4.1868
  do i=1,nnuc
     
     fra(:,i)=fra(:,i) - deriv(:,i)*forcefac

  end do
  
    
END SUBROUTINE add_force_qcff
