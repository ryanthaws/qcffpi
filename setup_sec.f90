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

SUBROUTINE setup_sec
  
  USE commondata
  USE classical, ONLY: nbonds,bonds,pi_iden,mass,vel,atype,thetas,nthetas,&
       amass,mass_tot
  USE ppv
  USE quantum
  USE qcfio, ONLY: gen_io
  
  
  IMPLICIT NONE
  
  !local variables
  LOGICAL :: hit1,hit2
  INTEGER :: i,j,k,ii,jj,kk,ib1,ib2,itemp,ind1,ind2,ierr,itheta
  INTEGER, DIMENSION(3,nbonds) :: ibeta_temp
  INTEGER, DIMENSION(nnuc) :: tind_nopi
  REAL(DP) :: gausvar,ekjmol,eng_KE,velscale,deltaabcc,deltaabcc13,efactor,tfactor,tempi
  REAL(DP), DIMENSION(3) :: ftemp,vtemp
  REAL(DP), DIMENSION(5) :: ztemp
  REAL :: ftemp2
  
  efactor=0.01
  tfactor=2.0d0*120.27/(3.d0*dfloat(nnuc)-3.d0) !to Kelvin
  if(log_ppv) then
     open(13,file='junction.dat',status='old',iostat=ierr)
     
     if (ierr .ne. 0) then
        call fileerror(ierr,"junction.dat")
        stop
     end if
     
     read(13,*,iostat=ierr) ntor_junc
     if(ierr .ne. 0) call filereaderror(ierr,"param.inpt","while reading ntor_junc")
     !APW inv_dih is a waste of memory
     ALLOCATE(junction(6,ntor_junc),inv_dih(nnuc,nnuc),junc_flg(nnuc),inv_dih_nopi(nnuc))
     ALLOCATE(torsion(ntor_junc))
     do i=1,ntor_junc
        read(13,*,iostat=ierr) (junction(j,i),j=1,6)
        if(ierr .ne. 0) call filereaderror(ierr,"param.inpt","while reading junction")
        inv_dih(junction(1,i),junction(2,i))=i
        junc_flg(junction(1,i))=.true.
        junc_flg(junction(2,i))=.true.
     end do
     close(13)
  end if
  
  !unit mass
  unitm=1.007825     
  unitl=150.0        !angstrom  APW why this value?
  !unit energy
  unite=unitm*amu*(unitl*1.0E-10/(deltat*1.0E-12))**2
  !convert from dimensionless energy to kJ/mol
  ekjmol=.5*avogad*unite*1.0E-3
  
  !c  ibeta matrix for PI atoms.                               c
  !c  Gives connection between PI atoms. Ibeta(i,3) is a code. c
  !c  In PPV PI atoms are just CARBON. Code = 1                c
  itemp=0
  do i=1,nbonds
     !APW ibeta_pairs is said to be the number of carbon-carbon bonds
     !APW the code would suggest a more precise definition is the pi-pi bonds
     hit1=.false.
     hit2=.false.
     ind1=bonds(1,i)
     ind2=bonds(2,i)
     do j=1,npi
        if(pi_iden(j).eq.ind1) hit1=.true.
        if(pi_iden(j).eq.ind2) hit2=.true.
     end do
     !if both members of bond i are pi species they belong in beta
     if(hit1.and.hit2) then
        itemp=itemp+1
        ibeta_temp(1,itemp)=ind1
        ibeta_temp(2,itemp)=ind2
        ibeta_temp(3,itemp)=1
     end if
  end do
  ibeta_pairs=itemp
  !allocate actual beta array
  ALLOCATE(ibeta(3,ibeta_pairs),pibond(ibeta_pairs))
  write(gen_io,*) 'total pairs involved in beta integrals = ',ibeta_pairs
  write(gen_io,*) '------------------------------------------------------'
  write(gen_io,*) 'pairs involved are...'
  !fill beta array
  do i=1,ibeta_pairs
     ibeta(:,i)=ibeta_temp(:,i)
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
     u_conf=0.d0
     write(gen_io,*) 'setting up tortional degrees of freedom'
     write(gen_io,*) '---------------------------------------'
     write(gen_io,*) 'total diedral used in resonance modulation=',ntor_junc
     
     ndih_nopi=0
     inv_dih_nopi=0
     do i=1,ntor_junc
        write(gen_io,*) 'atoms involved in torsion',i,'are',junction(3,i),junction(4,i),junction(5,i),junction(6,i)
        if(atype(junction(3,i)).ne.5) then
           if(inv_dih_nopi(junction(3,i)).eq.0) then
              ndih_nopi=ndih_nopi+1
              inv_dih_nopi(junction(3,i))=ndih_nopi
              tind_nopi(ndih_nopi)=junction(3,i)
           end if
           
        else if(atype(junction(4,i)).ne.5) then
           if(inv_dih_nopi(junction(4,i)).eq.0) then
              ndih_nopi=ndih_nopi+1
              inv_dih_nopi(junction(4,i))=ndih_nopi
              tind_nopi(ndih_nopi)=junction(4,i)
           end if
           
        else if(atype(junction(5,i)).ne.5) then
           if(inv_dih_nopi(junction(5,i)).eq.0) then
              ndih_nopi=ndih_nopi+1
              inv_dih_nopi(junction(5,i))=ndih_nopi
              tind_nopi(ndih_nopi)=junction(5,i)
           end if
           
        else if(atype(junction(6,i)).ne.5) then
           if(inv_dih_nopi(junction(6,i)).eq.0) then
              ndih_nopi=ndih_nopi+1
              inv_dih_nopi(junction(6,i))=ndih_nopi
              tind_nopi(ndih_nopi)=junction(6,i)
           end if
        end if
     end do
  end if
  ALLOCATE(ind_dih_nopi(ndih_nopi))
  do i=1,ndih_nopi
     ind_dih_nopi(i)=tind_nopi(i)
  end do
  
  !APW setup masses in units of unitm
  mass=0.0d0
  
  mass(1)=unitm
  mass(2)=15.9994 !oxygen
  mass(3)=14.0067 !nitrogen
  mass(4)=12.0107 !carbon (C)
  mass(5)=12.0107 !carbon (A)
  mass(6)=12.0107 !carbon (B)
  !mass(2)=15.9994/unitm !oxygen
  !mass(3)=14.0067/unitm !nitrogen
  !mass(4)=12.0107/unitm !carbon (C)
  !mass(5)=12.0107/unitm !carbon (A)
  !mass(6)=12.0107/unitm !carbon (B)
  
  !compute total mass
  mass_tot=0.d0
  do i=1,nnuc
     mass_tot = mass_tot + mass(atype(i))
     amass(i) = mass(atype(i))
     !check to make sure mass is defined
     if(mass(atype(i)).eq. 0.0) then
        write(gen_io,*) '***Warning: mass not defined for type',atype(i),'***'
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
  
  !sets the width of the velocity distribution
  gausvar=dsqrt(boltz*Tmb)
  velscale=2.454058E+11 !(sqrt(J/amu)/(Angstrom/picosecond))
  ekjmol=.01 !(amu*(angstrom/ps)^2/(kJ/mol))
  tempi=0
  do while((tempi-tmb)**2.gt.9.0)
     vtemp=0.0d0
     eng_KE=0.0d0
     do i=1,nnuc
        ftemp(2)=dsqrt(amass(i))
        do k=1,3
           !gaussian random number
           call gasdev_s(ftemp2)
           vel(k,i)=velscale*ftemp2*gausvar/ftemp(2)
           vtemp=vtemp+vel(k,i)*ftemp(2)**2
        end do
        !for center of mass momentum (three dimensional array)
        !vtemp=vtemp+vel(:,i)*ftemp(2)**2
     end do
  
  
     do i=1,nnuc
        !subtract center of mass velocity
        vel(:,i)=vel(:,i)-vtemp/mass_tot
        !compute kinetic energy (amu*(angstrom/ps)^2)
        eng_KE=eng_KE+0.5*amass(i)*(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)
     end do
     write(gen_io,*) 'Kinetic energy(amu*(ang/ps)^2,kJ/mol) =',eng_KE,eng_KE*ekjmol
     tempi=eng_KE*efactor*tfactor
  end do
  
  
  !read in bond-order matrix for PPP calculations
  open(13,file='bo.dat',status='old',iostat=ierr)
  if (ierr .ne. 0) then
     call fileerror(ierr,"bo.dat")
     stop
  end if

  do i=1,npi
     read(13,*,iostat=ierr) (bomat(i,j),j=1,npi)
     if(ierr .ne. 0) call filereaderror(ierr,"bo.dat","while reading bomat")
  end do
  close(13)
  
  ALLOCATE(ctully(0:ntully),fexa(3,nnuc,0:ntully,0:ntully),fexa_l(3,nnuc,0:ntully,0:ntully))
  ALLOCATE(fdv(0:ntully,0:ntully),fdv_l(0:ntully,0:ntully))
  ALLOCATE(dar(3,nnuc))
  
  ctully=0.d0 !zero out tully coefficients
  ctully(iexstate)=(1.d0,0.d0)

  !APW set z_atom (formally facion), the nuclear charge on each atom
  ALLOCATE(z_atom(nnuc))
  !set atom type Z-values
  !APW this 'z_atom' vector is needlessly long, with object orientation it should be better though
  ztemp=0.d0
  ztemp(5)=1.d0
  do i=1,nnuc
     if(log_delre) then
        z_atom(i)=sigma(i)+ztemp(atype(i))
     else
        z_atom(i)=ztemp(atype(i))
     end if
  end do
  
  !APW formally setup_delta_ppv
  ALLOCATE(deltaab(npi,npi))
  
  !c--- Set Variables --------------------------------------
  deltaabCC =    0.21889
  deltaabCC13 =  3.8689e-02
  
  deltaab = 0.d0
  
  !c--- Interaction 1-2 in PI system ------------------------
  
  do i=1,ibeta_pairs
     ib1=inv_pi(ibeta(1,i))
     ib2=inv_pi(ibeta(2,i))
     deltaab(ib1,ib2)=deltaabCC
  enddo
  
  !c--- Interaction 1-3 in PI system ------------------------
  
  do itheta=1,nthetas
     i=thetas(1,itheta)
     j=thetas(2,itheta)
     k=thetas(3,itheta)
     
     !HERE
     if(atype(i).eq.5.AND.atype(k).eq.5) then
        ii=inv_pi(i)
        kk=inv_pi(k)
        deltaab(ii,kk)=deltaABCC13
     end if
  end do
  
  !c--- Copy Upper Triangular in the Lower Part -------------
  do i=1,npi
     do j=1,i-1
        deltaab(i,j)=deltaab(j,i)
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
  
end SUBROUTINE setup_sec
