!------------------------------------------------------------
!                                                           !
!  subroutine setting the intial values                     !
!  the global variables                                     !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
!   Author Adam P. Willard 2010                             !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!

SUBROUTINE setup_vars
  
  USE commondata, ONLY: nnuc,npi,DP,boltz,tmb,log_cicalc,log_delre,&
       log_dumpvel,log_exdyn,log_ppv,log_restart,log_rescale_i,norb,&
       log_addcos
  USE classical
  USE quantum
  USE qcfio
  USE PPV
  
  IMPLICIT NONE

  !local variables
  CHARACTER*1 :: c1,c2
  LOGICAL :: hit1,hit2,ltemp
  INTEGER :: i,j,k,ii,jj,kk,ib1,ib2,itemp,ind1,ind2,ierr,itheta
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ibeta_temp
  INTEGER, DIMENSION(nnuc) :: tind_nopi
  REAL(DP) :: gausvar,eng_KE,velscale,deltaabCC,deltaabCN,deltaabCS,deltaabCO,&
       deltaabCC13,deltaabCN13,deltaabCS13,deltaabCO13,efactor,tfactor,tempi
  REAL(DP), DIMENSION(3) :: ftemp,vtemp
!  REAL(DP), DIMENSION(5) :: ztemp
  REAL(DP), DIMENSION(9) :: ztemp    ! AY changed
  REAL :: ftemp2
  
  efactor=0.01 !(amu*(angstrom/ps)^2/(kJ/mol))
  tfactor=2.0d0*120.27/(3.d0*dfloat(nnuc)-3.d0) !to Kelvin
  
  !first allocate arrays
  !ALLOCATE (atoms(nnuc))
  ALLOCATE (mass(nspec))
  ALLOCATE (deriv(3,nnuc),dderiv(3,3,nnuc,nnuc))
  ALLOCATE (atype(nnuc),amass(nnuc),inv_pi(nnuc))
  !APW are pos_m and pos_l used as global variables
  ALLOCATE (pos(3,nnuc),vel(3,nnuc),vel_l(3,nnuc),pos_m(3,nnuc),pos_l(3,nnuc))
  ALLOCATE (fra(3,nnuc),nonbs(nnuc,nnuc))
  fra=0.0d0
    
  ALLOCATE (bomat(npi,npi),mocoef(npi,npi),pi_iden(npi),mocoef_old(npi,npi))
  
  if(log_cicalc) then
     ALLOCATE (ciener(nci),osc(nci),ciener_l(nci))
     ALLOCATE (osc_vec(3,nci),cieig(nci,nci),cieig_old(nci,nci),ciind(nci,2),ciind_old(nci,2))
  end if
  if(log_exdyn) then
     ALLOCATE(qr(3,npi,npi,npi),utmp(norb+1:npi,norb,npi,npi))
     ALLOCATE(aeainv(npi*(npi-1)/2,npi*(npi-1)/2))
  end if
  
  if(log_ppv) then
     open(13,file='junction.inpt',status='old',iostat=ierr)
     if (ierr .ne. 0) then
        call fileerror(ierr,"junction.inpt")
        stop
     end if
     read(13,*,iostat=ierr) ntor_junc
     if(ierr .ne. 0) call filereaderror(ierr,"junction.inpt","while reading ntor_junc")
     !APW inv_dih is a waste of memory
     ALLOCATE(junction(6,ntor_junc),inv_dih(nnuc,nnuc),junc_flg(nnuc),inv_dih_nopi(nnuc),torsion(ntor_junc),dcosdr(3,6,ntor_junc))
     junc_flg=.false.
     inv_dih=0
     do i=1,ntor_junc
        read(13,*,iostat=ierr) (junction(j,i),j=1,6)
        if(ierr .ne. 0) call filereaderror(ierr,"junction.inpt","while reading junction")
        inv_dih(junction(1,i),junction(2,i)) = i
        junc_flg(junction(1,i)) = .true.
        junc_flg(junction(2,i)) = .true.
     end do
     close(13)
  end if
  
  !APW read in positions and connectivity
  if(log_restart) then
     !read in restart file
     write(gen_io,*) 'started reading restart'
     call read_restart
     write(gen_io,*) 'finished reading restart'
     !open(13,file='bo.inpt',status='old',iostat=ierr)
     !if (ierr .ne. 0) then
     !   call fileerror(ierr,"bo.inpt")
     !   stop
     !end if
     
     !do i=1,npi
     !   read(13,*,iostat=ierr) (bomat(i,j),j=1,npi)
     !   if(ierr .ne. 0) call filereaderror(ierr,"bo.inpt","while reading bomat")
     !   end do
   !  close(13)
  else
     !read in chem3d file
     call molein
     !initialize bond-order matrix for PPP calculations
     bomat = 0.0d0
     do i= 1,npi
        bomat(i,i) = 1.0d0
     end do
     !APW we now read tully coefs from restart file
     if(log_exdyn) then
        ALLOCATE(ctully(0:ntully))
        do i=0,ntully
           ctully(i) = (0.0d0,0.0d0)
        end do
        ctully(iexstate) = (1.0d0,0.0d0)
     end if
  end if
  
  ALLOCATE(ibeta_temp(3,nbonds))
  
  !c  ibeta matrix for PI atoms.                               c
  !c  Gives connection between PI atoms. Ibeta(i,3) is a code. c
  !c  In PPV PI atoms are just CARBON. Code = 1                c
  itemp=0
  do i=1,nbonds
     !APW ibeta_pairs is said to be the number of carbon-carbon bonds
     !APW the code would suggest a more precise definition is the pi-pi bonds
     hit1 = .false.
     hit2 = .false.
     ind1 = bonds(1,i)
     ind2 = bonds(2,i)
     do j=1,npi
        if(pi_iden(j).eq.ind1) hit1 = .true.
        if(pi_iden(j).eq.ind2) hit2 = .true.
     end do
     !if both members of bond i are pi species they belong in beta
     if(hit1.and.hit2) then
        itemp = itemp + 1
        ibeta_temp(1,itemp) = ind1
        ibeta_temp(2,itemp) = ind2
        ibeta_temp(3,itemp) = 1
     end if
  end do
  ibeta_pairs = itemp

  !allocate actual beta array
  !APW changed
  ALLOCATE(ibeta(3,ibeta_pairs),pibond(ibeta_pairs),st_betalin(ibeta_pairs,2),st_fbetalin(ibeta_pairs,3),inv_ibeta(nnuc,nnuc))
  !AY added
  if(log_addcos) ALLOCATE(st_exp2mu(ibeta_pairs,2))
  write(gen_io,*) 'total pairs involved in beta integrals = ',ibeta_pairs
  write(gen_io,*) '------------------------------------------------------'
  write(gen_io,*) 'pairs involved are...'
  !fill beta array
  inv_ibeta = 0
  do i= 1,ibeta_pairs
     ibeta(:,i) = ibeta_temp(:,i)
     !APW changed
     inv_ibeta(ibeta(1,i),ibeta(2,i)) = i
     inv_ibeta(ibeta(2,i),ibeta(1,i)) = i
     write(gen_io,*) ibeta(1,i),'--',ibeta(2,i),'   ',ibeta(3,i)
  end do
  
  !c------------------------------------------------------------
  !c                                                           c
  !c  subroutine setting variables need for B-cos(theta) calc. c
  !c                                                           c
  !c  junction : atoms forming junction bond. Read in INPUT    c
  !c  ntor(i)  : for the ith junction return the number of     c
  !c             torsion involved around this bond.            c
  !c  ntor_jun : total number of bond involved in junction     c
  !c  at_dih(i,j,4) : for a given junction and a given torsion c
  !c                  return the 4 atoms involved in the died  c
  !c  inv_dih(i1,i2) : for a couple of atoms give index of     c
  !c                   junction bond if it exist, or zero      c
  !c  dih_noPI : # of noPI atoms involved in torsional junc.   c
  !c  inv_dih_noPI(j): for the Jth atoms in all list numbering c
  !c                   return the index in restrained list     c
  !c                   list                                    c
  !c-----------------------------------------------------------c
  if (log_ppv) then
     u_conf = 0.d0
     write(gen_io,*) 'setting up tortional degrees of freedom'
     write(gen_io,*) '---------------------------------------'
     write(gen_io,*) 'total diedral used in resonance modulation=',ntor_junc

     ! get no pi atoms index involving junction parts  AY changed
     call  get_nopi_junc
     if(log_exdyn) ALLOCATE(qr_nopi(3,npi,npi,ndih_nopi))

  end if  ! log_ppv  AY added

  !APW setup masses in units of unitm
  mass=0.0d0

  if(locate_atom('H').ne.0) mass(locate_atom('H'))=1.00794
  if(locate_atom('C').ne.0) mass(locate_atom('C'))=12.0107
  if(locate_atom('A').ne.0) mass(locate_atom('A'))=12.0107
  if(locate_atom('B').ne.0) mass(locate_atom('B'))=12.0107
  if(locate_atom('O').ne.0) mass(locate_atom('O'))=15.9994
  if(locate_atom('Q').ne.0) mass(locate_atom('Q'))=15.9994 ! AY added
  if(locate_atom('N').ne.0) mass(locate_atom('N'))=14.0067
  if(locate_atom('M').ne.0) mass(locate_atom('M'))=14.0067
  if(locate_atom('S').ne.0) mass(locate_atom('S'))=32.0650
  if(locate_atom('T').ne.0) mass(locate_atom('T'))=32.0650
  !if(locate_atom('G').ne.0) mass(locate_atom('G'))=12.0107*60.0d0
  if(locate_atom('G').ne.0) mass(locate_atom('G'))=12.0107*60.0d0
  if(locate_atom('E').ne.0) mass(locate_atom('E'))=12.0107
  
  !compute total mass
  mass_tot=0.d0
  do i=1,nnuc
     mass_tot = mass_tot + mass(atype(i))
     amass(i) = mass(atype(i))
     !check to make sure mass is defined
     if(mass(atype(i)).eq. 0.0) then
        write(gen_io,*) '***Warning: mass not defined for type',atype(i),species(atype(i)),'***'
        stop
     end if
  end do
  
  !APW old types were defined differntly in the ppv routine?
  !*  +---------------------------------------------------------------+
  !*  |  RECIPROCAL solute MASS:                                      |
  !*  +---------------------------------------------------------------+
  !C   Atom Type atype        ITYPE itype
  !c    A (carbon PI)    5  = 1
  !c    O (oxygen no PI) 2  = 2
  !c    C (carbon no PI) 4  = 3
  !c    H (hydrogen no PI)1 = 4
  !c    M (girifalco atom)8 = 5

  write(gen_io,*) '--------------------------------------------------'
  write(gen_io,*) 'setup maxwell distribution of velocities'
  if(log_rescale_i) then
     vel=0.0d0
     !sets the width of the velocity distribution
     gausvar  = dsqrt(boltz*Tmb)
     velscale = 2.454058E+11 !(sqrt(J/amu)/(Angstrom/picosecond))
     tempi=0
     do while((tempi-tmb)**2.gt.9.0)
        vtemp  = 0.0d0
        eng_KE = 0.0d0
        do i=1,nnuc
           ftemp(2) = dsqrt(amass(i))
           do k=1,3
              !gaussian random number
              call gasdev_s(ftemp2)
              vel(k,i) = velscale*ftemp2*gausvar/ftemp(2)
              vtemp(k) = vtemp(k) + vel(k,i)*ftemp(2)**2
           end do
           !for center of mass momentum (three dimensional array)
           !vtemp=vtemp+vel(:,i)*ftemp(2)**2
        end do
        
        
        do i=1,nnuc
           !subtract center of mass velocity
           vel(:,i) = vel(:,i) - vtemp/mass_tot
           !compute kinetic energy (amu*(angstrom/ps)^2)
           eng_KE =eng_KE+0.5*amass(i)*(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)
        end do
        write(gen_io,*) 'Kinetic energy(amu*(ang/ps)^2,kJ/mol) =',eng_KE,eng_KE*efactor
        tempi = eng_KE*efactor*tfactor
     end do
  end if
  
  ALLOCATE(fexa(3,nnuc,0:ntully,0:ntully),fexa_l(3,nnuc,0:ntully,0:ntully))
  ALLOCATE(fdv(0:ntully,0:ntully),fdv_l(0:ntully,0:ntully),fdv_m(ntully,ntully))
  ALLOCATE(dar(3,nnuc))
  
  fexa   = 0.0d0;
  fexa_l = 0.0d0;

  !APW set z_atom (formally facion), the nuclear charge on each atom
  ALLOCATE(z_atom(nnuc))
  !set atom type Z-values
  !APW this 'z_atom' vector is needlessly long, with object orientation it should be better though
  ztemp = 0.d0
  !APW this is the core charge on the atoms and should be equal to the number of atoms contributed to the pi-system.
  if(locate_atom('A').ne.0) then
     ztemp(locate_atom('A')) = 1.0d0
  end if
  if(locate_atom('S').ne.0) then
     ztemp(locate_atom('S')) = 2.0d0
  end if
  if(locate_atom('M').ne.0) then
     ztemp(locate_atom('M')) = 1.0d0
  end if
  if(locate_atom('N').ne.0) then
     ztemp(locate_atom('N')) = 2.0d0
  end if
  do i=1,nnuc
     if(log_delre) then
        z_atom(i) = sigma(i) + ztemp(atype(i))
     else
        z_atom(i) = ztemp(atype(i))
     end if
  end do
  
  !APW formally setup_delta_ppv
  ALLOCATE(deltaab(npi,npi))
  
  !APW adjust this for sulfers (or generalize)
  !c--- Set Variables --------------------------------------
  deltaabCC   =  0.21889
  deltaabCC13 =  3.8689e-02
  deltaabCN   =  0.19536
  deltaabCN13 =  2.6384E-02
  !APW uning N parameters for S for the time being
  deltaabCS   =  0.19536
  deltaabCS13 =  2.6384E-02
  deltaabCO   =  0.22207
  deltaabCO13 =  2.5509E-02
  
  deltaab     =  0.d0
  
  !c--- Interaction 1-2 in PI system ------------------------
  do i= 1,ibeta_pairs
     ltemp = .true.
     ib1=inv_pi(ibeta(1,i))
     ib2=inv_pi(ibeta(2,i))
     c1 = species(atype(ibeta(1,i)))
     c2 = species(atype(ibeta(2,i)))
     if(c1.eq.'A'.and.c2.eq.'A') then 
        deltaab(ib1,ib2) = deltaabCC
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'N').or.(c1.eq.'N'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCN
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'M').or.(c1.eq.'M'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCN
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'S').or.(c1.eq.'S'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCS
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'O').or.(c1.eq.'O'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCO
        ltemp = .false.
     end if
     !deltaab(ib1,ib2)=deltaabCC
     if(ltemp) then
        write(gen_io,*) '***warning: atom pair ',c1,c2,' not found in deltaab12***'
        call cleanup
        stop
     end if
  enddo
  
  !c--- Interaction 1-3 in PI system ------------------------
  do itheta= 1,nthetas
     ltemp = .true.
     i = thetas(1,itheta)
     j = thetas(2,itheta)
     k = thetas(3,itheta)
     ii = inv_pi(i)
     kk = inv_pi(k)
     c1 = species(atype(i))
     c2 = species(atype(k))
     if(c1.eq.'A'.and.c2.eq.'A') then
        deltaab(ii,kk) = deltaabCC13
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'N').or.(c1.eq.'N'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCN13
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'S').or.(c1.eq.'S'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCS13
        ltemp = .false.
     end if
     if((c1.eq.'A'.and.c2.eq.'O').or.(c1.eq.'O'.and.c2.eq.'A')) then 
        deltaab(ib1,ib2) = deltaabCO13
        ltemp = .false.
     end if
     !APW theta includes contributions from non-pi atoms
     !if(ltemp) then
     !   write(gen_io,*) '***warning: atom pair ',c1,c2,' not found in deltaab13***'
     !   call cleanup
     !   stop
     !end if
  end do
  
  !c--- Copy Upper Triangular in the Lower Part -------------
  do i= 1,npi
  do j= 1,i-1
     deltaab(i,j) = deltaab(j,i)
  end do
  end do
  
  !*  +---------------------------------------------------------------+
  !*  |  setup lookup table for Z vector method                       |
  !*  +---------------------------------------------------------------+
  !APW this may be a waste
  ALLOCATE(idxzvecq((npi-1)*npi/2),idxzvecp((npi-1)*npi/2))
  do i = 1,npi
  do j = 1,i-1
     idxzvecq((i-2)*(i-1)/2 + j) = i
     idxzvecp((i-2)*(i-1)/2 + j) = j
  enddo
  enddo
  
  !APW used
  
  !APW delre
  !if(use_delre.eq.1) then
  !   ALLOCATE (links(4,nnuc))
  !   links=0
  !end if

  !--- In/Out Files --------------------------------------------------
  
  !OPEN(unit=59,file='bo.out',status='unknown')
  open (mov_io,file='movie.xyz',     status='unknown',iostat=ierr)
  open(mo_io,  file='MO_coefs.dat',  status='unknown',iostat=ierr)
  open(stat_io,file='statis.dat',    status='unknown',iostat=ierr)
  open(tor_io, file='torsion.dat',   status='unknown',iostat=ierr)
  open(bond_io,file='bond.dat',      status='unknown',iostat=ierr)
  open(gap_io, file='gap.dat',       status='unknown',iostat=ierr)
  open(osc_io, file='osc_vec.dat',   status='unknown',iostat=ierr)
  open(temp_io,file='TEMP_K.dat',    status='unknown',iostat=ierr)
  open(orb_io, file='orb_energy.dat',status='unknown',iostat=ierr)
  open(bo_io, file='bo.dat',status='unknown',iostat=ierr)
  open(swap_io, file='swap.dat',status='unknown',iostat=ierr)
  
  if (log_exdyn) then
     open(tul_io,  file='tully.dat',     status='unknown',iostat=ierr)
     open(fdotv_io,file='FdotV.dat',     status='unknown',iostat=ierr)
     open(fssh_io, file='fssh_probs.dat',status='unknown',iostat=ierr)
     open(hop_io, file='hopping.dat',status='unknown',iostat=ierr)
  endif
  
  if (log_cicalc) open(cicoef_io,file='CI_coefs.dat',status='unknown',iostat=ierr)
  
  if(log_dumpvel) open(vel_io,file='velocity.dat',status='unknown',iostat=ierr)
  
  
end SUBROUTINE setup_vars

!AY
!-----------------------------------------------------------------
!  Get index of no pi atoms involving junction parts:
!     get variables of:
!         ndih_nopi
!         ind_dih_nopi(ndih_nopi)
!         inv_dih_nopi(nnuc)
!-----------------------------------------------------------------
SUBROUTINE get_nopi_junc

  USE commondata, ONLY: nnuc,npi,DP
  USE classical
  USE quantum
  USE qcfio
  USE PPV
  
  IMPLICIT NONE

  !local variables
  INTEGER :: i,j,k, i1,i2,i3,i4,i5,i6
  INTEGER, DIMENSION(nnuc) :: tind_nopi
!  REAL(DP) :: aaa

  ndih_nopi       = 0
  inv_dih_nopi(:) = 0

  do i= 1,ntor_junc

     i3 = junction(3,i)
     i4 = junction(4,i)
     i5 = junction(5,i)
     i6 = junction(6,i)

     write(gen_io,*) 'atoms involved in torsion',i,'are', i3,i4,i5,i6

     if( atype(i3).ne.locate_atom('A') .and.   &
         atype(i3).ne.locate_atom('O') .and.   &
         atype(i3).ne.locate_atom('S') )then

        if(inv_dih_nopi(i3)==0) then
           ndih_nopi             = ndih_nopi + 1
           inv_dih_nopi(i3)      = ndih_nopi
           tind_nopi(ndih_nopi)  = i3
        end if
     endif

     if( atype(i4).ne.locate_atom('A') .and.   &
         atype(i4).ne.locate_atom('O') .and.   &
         atype(i4).ne.locate_atom('S') )then

        if(inv_dih_nopi(i4)==0) then
           ndih_nopi            = ndih_nopi + 1
           inv_dih_nopi(i4)     = ndih_nopi
           tind_nopi(ndih_nopi) = i4
        end if
     endif
        
     if( atype(i5).ne.locate_atom('A') .and.   &
         atype(i5).ne.locate_atom('O') .and.   &
         atype(i5).ne.locate_atom('S') )then

        if(inv_dih_nopi(i5)==0) then
           ndih_nopi            = ndih_nopi + 1
           inv_dih_nopi(i5)     = ndih_nopi
           tind_nopi(ndih_nopi) = i5
        end if
     endif
        
     if( atype(i6).ne.locate_atom('A') .and.   &
         atype(i6).ne.locate_atom('O') .and.   &
         atype(i6).ne.locate_atom('S') )then

        if(inv_dih_nopi(i6)==0) then
           ndih_nopi            = ndih_nopi + 1
           inv_dih_nopi(i6)     = ndih_nopi
           tind_nopi(ndih_nopi) = i6
        end if
     end if

  end do

  ALLOCATE( ind_dih_nopi(ndih_nopi) )
  do i= 1,ndih_nopi
     ind_dih_nopi(i) = tind_nopi(i)
  end do

end SUBROUTINE get_nopi_junc

