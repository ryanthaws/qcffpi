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
SUBROUTINE read_restart
  
  USE commondata, ONLY: nnuc,log_exdyn,npi
  USE classical, ONLY: pos,vel,nbonds,nthetas,nphis,atype,&
       bonds,thetas,phis,species,pispec,npispec,pi_iden,&
       nonbs,nphis_in
  USE quantum, ONLY: fexa_l,ntully,nnopi,inv_pi,ind_nopi,iexstate,bomat,ctully,alpha_mu_i,gamma_screen_i
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  LOGICAL :: ltemp
  CHARACTER*3 :: atemp3
  CHARACTER*5 :: atemp5,atemp5b
  INTEGER :: i,j,k,ierr,itemp,ntbond,nnonb,ind1,ind2
  INTEGER, DIMENSION(nnuc) :: tind_nopi
  INTEGER, DIMENSION(2,1000) :: tbonds
  
  open(56,file='restart.inpt',status='old',iostat=ierr)
  if (ierr .ne. 0) then
     call fileerror(ierr,"restart.inpt")
     stop
  end if
  
  read(56,'(A5,2x,A3)') atemp5,atemp3
  if(atemp5.ne."atype") then
     write(gen_io,*) '**WARNING: read_restart, expected atype got',atemp5
     stop
  end if
  if(atemp3.ne."pos") then
     write(gen_io,*) '**WARNING: read_restart, expected pos got',atemp3
     stop
  end if
  
  do i=1,nnuc
     !read(56,'(I1,2x,3(F15.9,2x))') atype(i),pos(:,i)
     read(56,*) atype(i),pos(:,i)
  end do
  
  read(56,'(A3)') atemp3
  if(atemp3.ne."vel") then
     write(gen_io,*) '**WARNING: read_restart, expected pos got',atemp3
     stop
  end if
  do i=1,nnuc
     !read(56,'(4(F15.9,2x))') vel(:,i)
     read(56,*) vel(:,i)
  end do

  !APW gamma
  read(56,'(A5,2x,A5)') atemp5,atemp5b
  if(atemp5.ne."gamma") then
     write(gen_io,*) '**WARNING: read_restart, expected gamma got',atemp5
     stop
  end if
  if(atemp5b.ne."alpha") then
     write(gen_io,*) '**WARNING: read_restart, expected alpha got',atemp5b
     stop
  end if
  do i=1,nnuc
     !read(56,'(4(F15.9,2x))') vel(:,i)
     read(56,*) gamma_screen_i(i),alpha_mu_i(i)
  end do
     
  read(56,'(A5)') atemp5
  if(atemp5.ne."bonds") then
     write(gen_io,*) '**WARNING: read_restart, expected bonds got',atemp5
     stop
  end if
  
  nonbs=1  !APW nonbonded interaction matrix
  ntbond=0
  
  read(56,'(I5)') nbonds
  ALLOCATE(bonds(2,nbonds))
  do i=1,nbonds
     !read(56,'(2(I5,2x))') bonds(1,i),bonds(2,i)
     read(56,*) bonds(1,i),bonds(2,i)
     ntbond=ntbond+1
  end do
  read(56,'(I5)') nthetas
  ALLOCATE(thetas(3,nthetas))
  do i=1,nthetas
     !read(56,'(3(I5,2x))') thetas(1,i),thetas(2,i),thetas(3,i)
     read(56,*) thetas(1,i),thetas(2,i),thetas(3,i)
     ntbond=ntbond+1
  end do
  read(56,*) nphis,nphis_in
  ALLOCATE(phis(4,nphis))
  do i=1,nphis
     !read(56,'(4(I5,2x))') phis(1,i),phis(2,i),phis(3,i),phis(4,i)
     read(56,*) phis(1,i),phis(2,i),phis(3,i),phis(4,i)
     ntbond=ntbond+1
  end do
  read(56,*) nnonb
  do i=1,nnonb
     read(56,*) ind1,ind2
     nonbs(ind1,ind2)=0
     nonbs(ind2,ind1)=0
  end do

  do i=1,npi
     read(56,*,iostat=ierr) (bomat(i,j),j=1,npi)
     if(ierr .ne. 0) call filereaderror(ierr,"restart.inpt","while reading bomat")
  end do

  read(56,*) atemp5,itemp
  if(atemp5.ne."state") then
     write(gen_io,*) '**WARNING: read_restart, expected state and got',atemp5
     stop
  end if
  iexstate=itemp
  !c:::  checking pi atoms
  
  !APW read in tully coefs
  read(56,*) atemp5,itemp
  if(log_exdyn) then
     ALLOCATE(ctully(0:ntully))
     if(atemp5.eq."ctull") then
        if(itemp.ne.ntully) then
           write(gen_io,*) '**WARNING: read_restart, ntully does not match restart.inpt'
           write(gen_io,*) 'expected',ntully,'got',itemp
           stop
        end if
        do i=0,ntully
           read(56,*) ctully(i)
        end do
     else
        !APW somtimes you kick off an exicited state run from a non-excited state configuration
        write(gen_io,*) '**WARNING: read_restart, not reading in tully coefficients', atemp5
        !APW moved from setup_vars.f90
        do i=0,ntully
           ctully(i) = (0.0d0,0.0d0)
        end do
        ctully(iexstate) = (1.0d0,0.0d0)
     end if
  end if
        
  
  itemp=0
  nnopi=0
  !APW new arrays pi2nuc nopi2nuc
  do i=1,nnuc
     ltemp=.false.
     do k=1,npispec
        if(pispec(k).eq.species(atype(i))) ltemp=.true.
     end do
     if(ltemp) then
        itemp = itemp+1
        pi_iden(itemp)=i
        !pi2nuc(itemp)=i
     else
        nnopi=nnopi+1
        tind_nopi(nnopi)=i
     end if
  end do
  
  ALLOCATE(ind_nopi(nnopi))
  do i=1,nnopi
     ind_nopi(i) = tind_nopi(i)
     !nopi2nuc(i)=tind_nopi(i)
  end do
  
  do i=1,npi
     inv_pi(pi_iden(i))=i
  end do
  
  if(itemp.ne.npi) then
     write(gen_io,*) '*** WARNING: input number of pi atoms (npi) does not equal actual number of pi atoms ***'
     stop
  end if

END SUBROUTINE read_restart
