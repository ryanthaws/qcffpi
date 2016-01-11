!------------------------------------------------------------
!                                                           !
!   this subroutine reads the chem3d file                   !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
!   Revised by Adam P. Willard 2010                         !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!
SUBROUTINE read3d()
  
  USE commondata, ONLY:nnuc
  USE qcfio, ONLY:gen_io,mol_io
  USE classical
  USE quantum, ONLY: gamma_screen,alpha_mu,gamma_screen_i,alpha_mu_i
  
  IMPLICIT NONE
  
  CHARACTER*3, DIMENSION(nnuc) :: acode_ppv
  LOGICAL :: ltemp
  INTEGER :: ierr
  INTEGER :: i,j,k
  INTEGER, DIMENSION(3) :: itemp
  INTEGER, DIMENSION(2,4000) :: bond_temp
  INTEGER, DIMENSION(nnuc) :: iachem_ppv
  INTEGER, DIMENSION(4,nnuc) :: ineigh_ppv
  DOUBLE PRECISION, DIMENSION(3*nnuc) :: x_ppv
  
  open(mol_io,file='initconfig.inpt',status='old',iostat=ierr)
  if (ierr .ne. 0) then
     call fileerror(ierr,"initconfig.inpt")
     stop
  end if
  
  do i=1,nnuc
     iachem_ppv(i)=0
     do j=1,4
        ineigh_ppv(j,i)=0
     end do
  end do
  
  read(mol_io,*,iostat=ierr) itemp(1)
  if(ierr .ne. 0) call filereaderror(ierr,"initconfig.inpt","while reading itemp(1)")
  if(itemp(1).ne.nnuc) then
     write(gen_io,*) "***error natom_ppv!=nnuc***"
     stop
  end if
  write(gen_io,*) 'natom_ppv=',nnuc,': nbond_ppv=',nbonds
  
  do i=1,nnuc
     read(mol_io,*,iostat=ierr) acode_ppv(i),itemp(1),(pos(k,i),k=1,3),iachem_ppv(i),(ineigh_ppv(k,i),k=1,4)
  end do
  
  do i=1,nnuc
     atype(i)=locate_atom(acode_ppv(i))
     !APW gamma
     gamma_screen_i(i)=gamma_screen(atype(i))
     alpha_mu_i(i)=alpha_mu(atype(i))
     if(acode_ppv(i).eq.'Z') then
        !set gamma parameter above but change everything else to type A parameters
        atype(i)=locate_atom('A')
     end if
     
     do j=1,3
        x_ppv(i*3-3+j)=pos(j,i)
     end do
     do j=1,4
        if(ineigh_ppv(j,i).ne.0) then
           ltemp=.true.
           do k=1,nbonds
              if(bond_temp(1,k).eq.i .and. bond_temp(2,k).eq.ineigh_ppv(j,i)) ltemp=.false.
              if(bond_temp(2,k).eq.i .and. bond_temp(1,k).eq.ineigh_ppv(j,i)) ltemp=.false.
           end do
           if(ltemp) then
              nbonds = nbonds+1
              bond_temp(1,nbonds)=i
              bond_temp(2,nbonds)=ineigh_ppv(j,i)
           end if
        end if
     end do
  end do
  
  !APW Allocate arrays that depend on the number of bonds
  ALLOCATE(bonds(2,nbonds))
  do i=1,nbonds
     bonds(1,i)=0
     bonds(2,i)=0
  end do
  
  do i=1,nbonds
     bonds(1,i)=bond_temp(1,i)
     bonds(2,i)=bond_temp(2,i)
  end do
  !c--- PPV VARIABLES -------------------------------------------
  
!!$
!!$      do i=1,natom
!!$         iachem_ppv(i)=iachem_ppv(i)
!!$         do j=1,4
!!$            ineigh_ppv(i,j)=ineigh_ppv(i,j)
!!$         enddo
!!$         acode_ppv(i)=acode(i)
!!$      enddo
!!$
!!$c       write(*,*)'now nbond = ',nbond
!!$c       stop
!!$
!!$      return
  
  close(mol_io)

END SUBROUTINE read3d
