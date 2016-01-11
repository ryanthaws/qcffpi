!------------------------------------------------------------
!                                                           !
!                                                           !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
!   Revised by Adam P. Willard 2010                         !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!


SUBROUTINE molein()
  
  USE commondata, ONLY:nnuc,npi
  USE classical, ONLY:npispec,species,pispec,nbonds,pi_iden,atype,locate_atom
  USE qcfio, ONLY:gen_io
  USE quantum, ONLY: inv_pi,nnopi,ind_nopi
  
  IMPLICIT NONE

  !c:::  local vars
  LOGICAL :: ltemp
  INTEGER :: ichem3d,ierr
  INTEGER :: i,j,k
  INTEGER, DIMENSION(2) :: itemp
  INTEGER, DIMENSION(nnuc) :: tind_nopi
  
  call read3d()
    
  write(gen_io,*)'natom =',nnuc
  write(gen_io,*)'nbond =',nbonds

  call connect() !determines molecular bonding
  call allpak() !packs bonding info into arrays
  
  !c:::  checking pi atoms
  
  itemp(1)=0
  nnopi=0
  do i=1,nnuc
     ltemp=.false.
     do k=1,npispec
        if(pispec(k).eq.species(atype(i))) ltemp=.true.
     end do
     if(ltemp) then
        itemp(1) = itemp(1)+1
        pi_iden(itemp(1))=i !APW nuclear coordinate of pi-atom itemp(1)
     else
        !if(atype(i).ne.locate_atom('A')) then
        nnopi=nnopi+1
        tind_nopi(nnopi)=i
     end if
  end do
  
  ALLOCATE(ind_nopi(nnopi))
  do i=1,nnopi
     ind_nopi(i) = tind_nopi(i)
  end do
  
  do i=1,npi
     inv_pi(pi_iden(i))=i
  end do
  
  if(itemp(1).ne.npi) then
     write(gen_io,*) '*** WARNING: input number of pi atoms (npi) does not equal actual number of pi atoms ***'
     stop
  end if

  !APW what is the point of the code below
  !itemp(1)=0
  !do i=1,nnuc
  !   if(atype(i).eq.locate_atom('N')) itemp(1)+1
  !   if(atype(i).eq.locate_atom('X')) itemp(1)+1
  !   if(atype(i).eq.locate_atom('Z')) itemp(1)+1
  !end do
  
END SUBROUTINE molein
