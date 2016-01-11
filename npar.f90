!------------------------------------------------------------
!                                                           !
!  this subroutine reads in and updates the parameters      !
!  of the classical potential functions.                    !
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

SUBROUTINE npar(ip)

  USE commondata, ONLY: log_delre,log_addcos,nnuc
  USE qcfio, ONLY: parm_io,gen_io
  USE classical
  USE quantum, ONLY: gamma_screen,alpha_mu,gamma_screen_i,alpha_mu_i
  
  IMPLICIT NONE
  
  !INTEGER :: locate_atom               !function
  !c:::  input vars
  INTEGER :: ip
  
  !c:::  local vars
  LOGICAL :: ltemp
  CHARACTER*1, DIMENSION(120) :: ctemp
  CHARACTER*1, ALLOCATABLE, DIMENSION(:,:) :: bond_spec,theta_spec,phi_spec,&
       beta_spec,nonb_spec
  INTEGER :: i,j,k,ii,jj,kk,lp,irow,icol
  INTEGER :: ierr,nlines
  INTEGER, DIMENSION(2) :: itemp
  DOUBLE PRECISION, DIMENSION(5) :: ftemp
  
  open(parm_io,file='qcff_parm.inpt',status='old',iostat=ierr)
  if (ierr .ne. 0) then
     call fileerror(ierr,"qcff_parm.inpt")
     stop
  end if

  read(parm_io,*,iostat=ierr) nspec
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nspec")
  ALLOCATE(species(nspec),theta_equiv(nspec),icvn(nspec),qatom(nspec),new_do(nspec))
  ALLOCATE(pair_code(nspec,nspec))
  
  do i=1,nspec
     read(parm_io,'(1A)',iostat=ierr) species(i)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading species")
  end do
  do i=1,nspec
     read(parm_io,*,iostat=ierr) theta_equiv(i),icvn(i)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading icvn")
  end do

  !APW assign a numerical value to each pair of atoms
  k=0
  do i=1,nspec
     do j=1,i
        k=k+1
        pair_code(j,i)=k           
        !with codes i,j respectively
        pair_code(i,j)=k
        write(gen_io,*) k,species(j),species(i)
     enddo
  enddo
  
  !APW substituting all ?-D (and ?-H) pair_codes with ?-H (and H-?)  
  if(locate_atom('D').ne.0) then
     do i=1,nspec 
        pair_code(i,locate_atom('D'))= pair_code(locate_atom('H'),i)
        pair_code(locate_atom('D'),i)= pair_code(locate_atom('H'),i)
     end do
  end if
  
  !APW this loop accomplishes the same as the commented call to stpar() in the older version of the code
  read(parm_io,*,iostat=ierr) nbond_type
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nbond")
  ALLOCATE(bond_pms(3,nbond_type),bond_code(nbond_type),bond_spec(2,nbond_type))
  
  do i=1,nbond_type
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,2),(ftemp(ii),ii=1,3)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading bond_pms")
     bond_spec(1,i)=ctemp(1)
     bond_spec(2,i)=ctemp(2)
     bond_code(i)=pair_code(locate_atom(ctemp(1)),locate_atom(ctemp(2)))
     do k=1,3
        bond_pms(k,i)=ftemp(k)
     end do
  end do

  
  !c:::  initialize theta param and update them
  read(parm_io,*,iostat=ierr) ntheta_type
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading ntheta")
  ALLOCATE(theta_pms(2,ntheta_type),theta_cub(3,ntheta_type),theta_spec(3,ntheta_type),theta_code(ntheta_type))
  
  do i=1,ntheta_type
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,3),(ftemp(ii),ii=1,5)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading theta pms")
     theta_spec(1,i)=ctemp(1)
     theta_spec(2,i)=ctemp(2)
     theta_spec(3,i)=ctemp(3)
     itemp(1) = pair_code(locate_atom(ctemp(1)),locate_atom(ctemp(3)))
     theta_code(i) = itemp(1)*100+locate_atom(ctemp(2))
     do k=1,2
        theta_pms(k,i)=ftemp(k)
     end do
     do k=1,3
        theta_cub(k,i)=ftemp(k+2)
     end do
     do k=1,3
        theta_equiv(locate_atom(ctemp(k)))=locate_atom(ctemp(k))
     end do
  end do

!!$ 70   ntp1=ntp+1
!!$      read (parm_io,1001) at1,at2,at3
!!$      if(at1.eq.bln(1)) go to 100
!!$      read(parm_io,*) (ct(i,ntp1),i=1,2),(cub(i,ntp1),i=1,3)
!!$ 
!!$      call stpar(at1,at3,at2,ct,itcod,2,ntp,lp)
!!$      if(lp.lt.ntp1) go to 80
!!$ 
!!$      icvt(iat(at1))=iat(at1)
!!$      icvt(iat(at2))=iat(at2)
!!$      icvt(iat(at3))=iat(at3)
!!$      go to 70
!!$ 
!!$   80 do 90 j=1, 3
!!$   90 cub(j,lp)=cub(j,ntp1)
!!$ 
!!$      go to 70
!!$ 
!!$  100 continue

  !c     initialize phi   param and update them
  read(parm_io,*,iostat=ierr) nphi_type
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nphi")
  ALLOCATE(phi_pms(3,nphi_type),phi_code(nphi_type),phi_spec(4,nphi_type))
  ALLOCATE(iatg(nphi_type),ncos(nphi_type),phi_in(nphi_type))
  
  do i=1,nphi_type
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,4),(ftemp(ii),ii=1,3),itemp(1),itemp(2)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading phi pms")
     do j=1,4
        phi_spec(j,i)=ctemp(j)
     end do
     phi_code(i)=pair_code(locate_atom(phi_spec(2,i)),locate_atom(phi_spec(3,i)))
     if(itemp(1) .eq. 0) itemp(1)=3
     if(itemp(2) .eq. 0) itemp(2)=1
     
     do k=1,3
        phi_pms(k,i)=ftemp(k)
     end do
     ncos(i)=itemp(1)
     iatg(i)=itemp(2)
  end do

  !APW for uniquely defined torsional pairs
  read(parm_io,*,iostat=ierr) nphib
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nphib")
  ALLOCATE(phib_pms(2,nphib),phib_pairs(2,nphib))
  
  do i=1,nphib
     read(parm_io,*,iostat=ierr) itemp(1),itemp(2),ftemp(1),ftemp(2)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading phib pms")
     
     phib_pairs(1,i)=itemp(1)
     phib_pairs(2,i)=itemp(2)
     
     phib_pms(1,i)=ftemp(1)
     phib_pms(2,i)=ftemp(2)
     
  end do

  itemp(1)=0
  do i=1,nphi_type
     if(iatg(i).ge.2) then
        itemp(1) = itemp(1)+1
        phi_in(itemp(1))=phi_code(i)
     end if
  end do
  nphi_in=itemp(1)
  do i=itemp(1)+1,nphi_type
     phi_in(i)=-1
  end do
  
  !c     initialize and update nonbond parameters
  read(parm_io,*,iostat=ierr) nnonb_type
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nnonb")
  ALLOCATE(nonb_pms(3,nnonb_type),nonb_code(nnonb_type),nonb_spec(2,nnonb_type))
  
  do i=1,nnonb_type
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,2),(ftemp(ii),ii=1,3)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nonb_pms")
     nonb_spec(1,i)=ctemp(1)
     nonb_spec(2,i)=ctemp(2)
     nonb_code(i)=pair_code(locate_atom(ctemp(1)),locate_atom(ctemp(2)))
     
     do k=1,3
        nonb_pms(k,i)=ftemp(k)
     end do
     icvn(locate_atom(ctemp(1)))=locate_atom(ctemp(1))
     icvn(locate_atom(ctemp(2)))=locate_atom(ctemp(2))
     
  end do
!!$  150 nbat1=nbat+1
!!$      read (parm_io,1000) at1,at2
!!$      if(at1.eq.bln(1)) go to 160
!!$      read(parm_io,*) (cnb(i,nbat1),i=1,3)
!!$c
!!$      call stpar(at1,at2,bln(1),cnb,nbcod,3,nbat,lp)
!!$
!!$      if(lp.lt.nbat1) go to 150
!!$c
!!$      icvn(iat(at1))=iat(at1)
!!$      icvn(iat(at2))=iat(at2)
!!$      go to 150
!!$c
!!$  160 continue


  !c     charges
  do i=1,nspec
     qatom(i)=0.d0
  end do
  !c  initialization and update of del re parameters
  
  read(parm_io,*,iostat=ierr) nlines
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nlines")
  do i=1,nlines
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,2),ftemp(1)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading gadrs")
     irow=locate_atom(ctemp(1))
     icol=locate_atom(ctemp(2))
     !APW delre
     !if(use_delre.eq.1) then
     !   gadr(irow,icol)=ftemp(1)
     !end if
  end do
  
  read(parm_io,*,iostat=ierr) nlines
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nlines")
  do i=1,nlines
     read(parm_io,*,iostat=ierr) ctemp(1),ftemp(1)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading new_do")
     itemp(1)=locate_atom(ctemp(1))
     new_do(itemp(1))=ftemp(1)
  end do

  read(parm_io,*,iostat=ierr) nlines
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nlines")
  do i=1,nlines
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,2),ftemp(1)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading eadrs")
     irow=locate_atom(ctemp(1))
     icol=locate_atom(ctemp(2))
  end do

  ALLOCATE( gamma_screen(nspec),alpha_mu(nspec))
  gamma_screen=0.0d0
  alpha_mu=0.0d0
  !APW gamma for different molecules with carbon
  ALLOCATE( gamma_screen_i(nnuc),alpha_mu_i(nnuc) )
  
  read(parm_io,*,iostat=ierr) itemp(2)
  
  if(ierr .ne.0) call filereaderror(ierr,"qcff_parm.inpt","while reading gamma_screens")
  do i=1,itemp(2)
     read(parm_io,*,iostat=ierr) ctemp(1),ftemp(1),ftemp(2)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading gamma_screen")
     itemp(1)=locate_atom(ctemp(1))
     gamma_screen(itemp(1))=ftemp(1)
     alpha_mu(itemp(1))=ftemp(2)
  end do
     
  !c
  !c     initialize alpha param and update them
  !c
  read(parm_io,*,iostat=ierr) npispec
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading npispec")
  ALLOCATE(pispec(npispec),alpha(3,npispec),pispec_code(npispec))
  
  do i=1,npispec
     read(parm_io,*,iostat=ierr) ctemp(1),(ftemp(ii),ii=1,3)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading alpha")
     pispec(i)=ctemp(1)
     itemp(1)=locate_atom(ctemp(1))
     pispec_code(i)=itemp(1)
     do j=1,3
        alpha(j,i)=ftemp(j)
     end do
  end do
     
  !c     initialize beta and gamma parameters
  !c
  !c
  !c     update beta parameters
  !c
  read(parm_io,*,iostat=ierr) nbeta_type
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nbeta_type")
  ALLOCATE(beta_spec(2,nbeta_type),beta(5,nbeta_type),beta_code(nbeta_type),gamma(4,nbeta_type))
  
  do i=1,nbeta_type
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,2),(ftemp(ii),ii=1,5)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading beta")
     beta_spec(1,i)=ctemp(1)
     beta_spec(2,i)=ctemp(2)
     beta_code(i)=pair_code(locate_atom(ctemp(1)),locate_atom(ctemp(2)))
     do j=1,5
        beta(j,i)=ftemp(j)
     end do
     
     !APW check to make sure beta atom is a pi atom
     ltemp=.true.
     do j=1,npispec
        if(ctemp(1).eq.pispec(j)) ltemp=.false.
     end do
     if(ltemp) then
        write(gen_io,*)'**** ',ctemp(1),' of beta(',ctemp(1),ctemp(2),' is not a pi atom (undefined by alpha) ***'
     end if
     ltemp=.true.
     do j=1,npispec
        if(ctemp(2).eq.pispec(j)) ltemp=.false.
     end do
     if(ltemp) then
        write(gen_io,*)'**** ',ctemp(2),' of beta(',ctemp(1),ctemp(2),' is not a pi atom (undefined by alpha) ***'
     end if
        
           
  end do
  
  !c
  !c     update gamma parameters
  !c
  read(parm_io,*,iostat=ierr) nlines
  if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading nlines")
  
  do i=1,nbeta_type
     read(parm_io,*,iostat=ierr) (ctemp(ii),ii=1,2),(ftemp(ii),ii=1,4)
     if(ierr .ne. 0) call filereaderror(ierr,"qcff_parm.inpt","while reading gamma")
     
     itemp(1)=pair_code(locate_atom(ctemp(1)),locate_atom(ctemp(2)))
     do j=1,nbeta_type
        if(itemp(1).eq.beta_code(j)) then
           do k=1,4
              gamma(k,j)=ftemp(k)
           end do
        end if
     end do
  end do
  close(parm_io)

  !APW ip determines if things should be printed
  if(ip.eq.0) return
  
  !c-----------------------------------------------------------------------
  !c
  !c     print the set of parameters that define the force field used
  !c
  !c
  !c-----------------------------------------------------------------------

  write(gen_io,*)'     bond streching energy constants for'
  write(gen_io,*)'     f(b)=kb1(dexp(-2kb3(b-kb2))-2dexp(-kb3(b-kb2)))'
  write(gen_io,*)'                          or for'
  write(gen_io,*)'     f(b)=kb1(b-kb2)**2-kb3'
  write(gen_io,*)'bond        kb1       kb2       kb3'
  do i=1,nbond_type
     write(gen_io,'(8a,3f10.3)') bond_spec(1,i),bond_spec(2,i),(" ",j=1,6),(bond_pms(j,i),j=1,3)
  end do
  
  write(gen_io,*)'     angle bending energy constants,for'
  write(gen_io,*)'     f(th)=kt1(th-kt2)**2+kq2(x1-p1)+kq1(q-kq3)**2'
  write(gen_io,*)'x1=th for ccc-angle,x1=q for all  ohers;p1=kt2 for ccc,p1=kq3 for others'
  write(gen_io,*)'angle       kt1        kt2       kq1       kq2      kq3'
  do i=1,ntheta_type
     write(gen_io,'(9a,5f10.3)')theta_spec(1,i),theta_spec(2,i),theta_spec(3,i),(" ",j=1,6),&
          (theta_pms(j,i),j=1,2),(theta_cub(j,i),j=1,3)
  end do

  write(gen_io,*)'     equivalent atoms in theta functions'
  write(gen_io,*)'     e.g. hch is equivalent to dcd,nnn equiv. to aaa'

  do i=1,nspec
     itemp(1)=0
     do j=1,nspec
        if(theta_equiv(j).eq.i) then 
           itemp(1)=itemp(1)+1
           ctemp(itemp(1))=species(j)
        end if
     end do
     if(itemp(1).gt.0) write(gen_io,'(15(1x,1a))') (ctemp(j),j=1,itemp(1))
  end do
  
  write(gen_io,*) '     torsional energy constants for'
  write(gen_io,*) '     f(phi)=kp1*dcos(n*phi)+kp3*dcos(phi)'
  write(gen_io,*) '            +p1(th1-th0)(th2-th0)*dcos(phi)'
  write(gen_io,*) '     p1=kp3 for hcch,p1=-dsqrt(kp2*kp3) for hccc,p1=kp2 for others'
  write(gen_io,*) '           kp1 is -kp1 for n=2'
  write(gen_io,*) 'phi-angle    kp1       kp2       kp3      n   iatg'
  do i=1,nphi_type
     write(gen_io,'(9a,3f10.3,2i5)') (" ",j=1,2),(phi_spec(j,i),j=1,4),(" ",j=1,3),&
          (phi_pms(j,i),j=1,3),ncos(i),iatg(i)
  end do

  write(gen_io,*) '     unique torsional energy constants'
  write(gen_io,*) '     f(phi)=g1*dcos(2*phi)+g2*dcos(4*phi)'
  write(gen_io,*) '       i    j     g1       g2 '
  do i=1,nphib
     write(gen_io,'(4a,2i5,2f10.3)') (" ",j=1,4),(phib_pairs(j,i),j=1,2),(phib_pms(j,i),j=1,2)
  end do

  write(gen_io,*) '     non-bond energy constants for'
  write(gen_io,*) '     f(r)=knb2*dexp(-knb3*r)-knb1/r**6'
  write(gen_io,*) 'atoms       knb1           knb2           knb3'

  do i=1,nnonb_type
     write(gen_io,'(10a,3e15.6)') (nonb_spec(j,i),j=1,2),(" ",j=1,8),(nonb_pms(j,i),j=1,3)
  end do

  write(gen_io,*) '     equivalent atoms in nonbond functions'
  do i=1,nspec
     itemp(1)=0
     do j=1,nspec
        if(icvn(j).eq.i) then
           itemp(1)=itemp(1)+1
           ctemp(itemp(1))=species(j)
        end if
     end do
     if(itemp(1).gt.0) write(gen_io,'(15(1x,1a))') (ctemp(j),j=1,itemp(1))
  end do
        
  write(gen_io,*) '     Atomic Charges'
  do i=1,6
     write(gen_io,'(4a,f7.3)') species(i),(" ",j=1,3),qatom(i)
  end do

  ! AY improved  2011.01  
  write(gen_io,*) '     Constants for Alpha Integrals'
!OLD  write(gen_io,*) '     alpha=al1-al2*dexp(-2bt3(b-bt2))dcos(phi)**2'
  write(gen_io,*) '     alpha_ii = al1 - sum_j[Z_i*g(i,j)]            '
  write(gen_io,*) '                - al2*exp(-2bt3(b-bt2))*cos(phi)**2'
  write(gen_io,*) 'atom     al1        al2       al3'
  do i=1,npispec
     write(gen_io,'(8a,6f10.5)') pispec(i),(" ",j=1,7),(alpha(j,i),j=1,3)
  end do

  ! AY improved  2011.01
  write(gen_io,*) '     constants for beta integrals'
  write(gen_io,*) '     beta=bt1*dexp(-bt3(b-bt2))(1+bt4(b-bt2))*cos(phi)'
!OLD  write(gen_io,*) '     beta=bt1*dexp(-bt3(b-bt2))(1+bt4(b-bt2)) dcos(phi)'
!OLD  write(gen_io,*) '          *(1+bt5*dcos(phi))/(1+bt5) '
  write(gen_io,*) 'atoms      bt1       bt2       bt3       bt4       bt5'
  do i=1, nbeta_type
     write(gen_io,'(8a,5f10.4)') (beta_spec(j,i),j=1,2),(" ",j=1,6),(beta(j,i),j=1,5)
  end do

  ! AY improved  2011.01
  write(gen_io,*) '     constants for gamma functions'
!OLD  write(gen_io,*) '     g(i,i)=al3+g4*dexp(-2bt3(b-bt2))dcos(phi)**2'
!OLD  write(gen_io,*) '                        and'
!OLD  write(gen_io,*) '     g(i,j)=(g1-g2)*dexp(-g3*r)+14.397/(r+14.397/g2)'
!OLD  write(gen_io,*) '             -g4/2*dexp(-2bt3(b-bt2))dcos(phi)**2'
!OLD  write(gen_io,*) '     ,but for atoms whose beta const. are not defined'
!OLD  write(gen_io,*) '     g(i,j)=14.397/(r+14.397*2/(al3(i)+al3(j)))'
  if(log_addcos) then
  write(gen_io,*) '     g(i,i)= al3 + g4*dexp(-2bt3(b-bt2))*cos(phi)**2  '
  write(gen_io,*) '     for j=i-1,i+1                                    '
  write(gen_io,*) '     g(i,j)= (g1-g2)*exp(-g3*r) + 14.397/(r+14.397/g2)'
  write(gen_io,*) '             -g4*exp(-2bt3(b-bt2))*cos(phi)**2        '
  write(gen_io,*) '     and for j is not equal to i,i-1,i+1              '
  write(gen_io,*) '     g(i,j)= (g1-g2)*exp(-g3*r) + 14.397/(r+14.397/g2)'
  else 
  write(gen_io,*) '     (Mataga-Nishimoto relation)           '
  write(gen_io,*) '     g(i,i) = al3                          '
  write(gen_io,*) '         and                               '
  write(gen_io,*) '     g(i,j) = 14.397/(r_ij + a_ij)         '
  write(gen_io,*) '     a_ij   = 2*14.397/(g(i,i)+g(j,j))     '
  endif
  write(gen_io,*)  'atoms         g1        g2        g3        g4'
  do i=1, nbeta_type
     write(gen_io,'(8a,4f10.4)') (beta_spec(j,i),j=1,2),(" ",j=1,6),(gamma(j,i),j=1,4)
  end do

  write(gen_io,*) '     constants in del re-sigma charges computation procedure'
  !APW delre
  !if(use_delre.eq.1) then
  !write(gen_io,*) 'atom   do   gadr '
  !write(gen_io,'(14x,15(a,5x))') (species(i),i=1,nspec)
  !do i=1,nspec
  !    write(gen_io,'(2x,a,2x,16f6.2)'),species(i),new_do(i),(gadr(i,j),j=1,15)
  !end do
  !write(gen_io,*) 'atom   eadr '
  !do i=1,nspec
  !   write(gen_io,'(2x,a,2x,15f6.2)'),species(i),(eadr(i,j),j=1,15)
  !end do
  !end if
  
  close(parm_io)
  
  DEALLOCATE(bond_spec,theta_spec,phi_spec,beta_spec,nonb_spec)
  
END SUBROUTINE npar

!c-----------------------------------------------------------------------
!c
!c     locate the atom in the list of given atoms.
!c
!c-----------------------------------------------------------------------

!APW moved into module classical.mod
!!$INTEGER FUNCTION locate_atom(atom)
!!$  
!!$  USE QCFIO, ONLY:gen_io
!!$  USE classical, ONLY:species,nspec
!!$  
!!$  IMPLICIT NONE
!!$  
!!$  CHARACTER*1 atom   ! character type
!!$  INTEGER i
!!$  
!!$  do i=1,nspec
!!$     if(atom.eq.species(i)) then
!!$        locate_atom=i
!!$        return
!!$     endif
!!$  enddo
!!$  
!!$  write(gen_io,*) "*** atom ",atom," is not one of the coded atoms ***"
!!$  stop 
!!$
!!$  return
!!$  
!!$END FUNCTION locate_atom
