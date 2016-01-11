!------------------------------------------------------------
!                                                           !
!   create a restart file                                   !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
!   Written by Adam P. Willard 2010                         !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!
SUBROUTINE write_restart
  
  USE commondata, ONLY: nnuc,log_exdyn,npi,istep
  USE classical, ONLY: pos,vel,nbonds,nthetas,nphis,atype,&
       bonds,thetas,phis,nonbs,nphis_in
  USE quantum, ONLY: fexa_l,ntully,bomat,iexstate,ctully,alpha_mu_i,gamma_screen_i

  IMPLICIT NONE
  
  INTEGER :: i,j,k,ierr,nnonb
  CHARACTER*17 :: tname
  
  if(istep.lt.10) write(tname,'(A16,I1)') 'tmprestart.00000',istep
  if(istep.ge.10 .and. istep.lt.100) write(tname,'(A15,I2)') 'tmprestart.0000',istep
  if(istep.ge.100 .and. istep.lt.1000) write(tname,'(A14,I3)') 'tmprestart.000',istep
  if(istep.ge.1000 .and. istep.lt.10000) write(tname,'(A13,I4)') 'tmprestart.00',istep
  if(istep.ge.10000 .and. istep.lt.100000) write(tname,'(A12,I5)') 'tmprestart.0',istep
  if(istep.ge.100000 .and. istep.lt.1000000) write(tname,'(A11,I6)') 'tmprestart.',istep
  
  
  open(55,file=tname,status='unknown',iostat=ierr)
  
  write(55,'(A5,2x,A3)') "atype","pos"
  do i=1,nnuc
     write(55,'(I1,2x,3(F22.16,2x))') atype(i),pos(:,i)
  end do

  write(55,'(A3)') "vel"
  do i=1,nnuc
     write(55,'(4(F22.16,2x))') vel(:,i)
  end do
  
  !APW gamma
  write(55,'(A5,2x,A5)') "gamma","alpha"
  do i=1,nnuc
     write(55,'(2(F22.16,2x))') gamma_screen_i(i),alpha_mu_i(i)
  end do
  
  write(55,'(A5)') "bonds"
  write(55,'(I5)') nbonds
  do i=1,nbonds
     write(55,'(2(I5,2x))') bonds(1,i),bonds(2,i)
  end do
  write(55,'(I5)') nthetas
  do i=1,nthetas
     write(55,'(3(I5,2x))') thetas(1,i),thetas(2,i),thetas(3,i)
  end do
  write(55,*) nphis,nphis_in
  do i=1,nphis
     write(55,'(4(I5,2x))') phis(1,i),phis(2,i),phis(3,i),phis(4,i)
  end do
  nnonb=0
  do i=1,nnuc
     do j=i,nnuc
        if(nonbs(i,j).eq.0) nnonb=nnonb+1
     end do
  end do
  write(55,*) nnonb
  do i=1,nnuc
     do j=i,nnuc
        if(nonbs(i,j).eq.0) write(55,*)i,j
     end do
  end do
  
  do i=1,npi
     write(55,*) (bomat(i,j),j=1,npi)
  end do

  write(55,*) "state",iexstate
  
  if(log_exdyn) then
     write(55,*) "ctull",ntully
     do i=0,ntully
        write(55,*) ctully(i)
     end do
  else
     write(55,*) "done",0
  end if
  
  close(55)
  
 END SUBROUTINE write_restart
