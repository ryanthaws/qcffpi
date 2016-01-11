!c------------------------------------------------------------
!c  subroutine writeing PDB output frame                     c
!c-----------------------------------------------------------c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm/utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c-----------------------------------------------------------c

!subroutine pdb_ppv(xna,yna,zna,unitl_ppv,istep)
SUBROUTINE make_pdb

  USE commondata, ONLY: istep,nnuc
  USE classical, ONLY: species,atype,pos
  USE qcfio, ONLY: mov_io
  
  IMPLICIT NONE
  
  INTEGER :: i,itemp
  CHARACTER*6 :: str
!!$  
!!$  !c---  Executive Stantements ----------------------------------
!!$  str='ATOM  '
!!$
!!$  write(mov_io,*) '#COMMENT CONFIGURATION T=',istep
!!$  do i=1,nnuc
!!$     write(mov_io,1000) str,i,atm_pdb(i),res_pdb(i),ires_pdb(i),xna(i)
!!$     $        *unitl_ppv,yna(i)*unitl_ppv,zna(i)*unitl_ppv
!!$      enddo
!!$
!!$1000  format(a6,i5,2x,a2,2x,a3,i6,4x,3f8.3) 
!!$      
!!$      WRITE(60,*) 'END'

  write(mov_io,*) nnuc
  write(mov_io,*) " " 
  do i=1,nnuc
     itemp=atype(i)
     write(mov_io,*) species(itemp),pos(:,i)
  end do
  
END SUBROUTINE make_pdb
