!c------------------------------------------------------------
!c  subroutine computing the electronic matrix of PPP theory c
!c  Fock Matrix F                                            c
!c  Bond Order Matrix B                                      c
!c  One e core Matrix H                                      c
!c-----------------------------------------------------------c
!c                                                           c
!c   Modified from original subroutine of JL                 c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm/utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c   gets called from inside TWIST_PPV.F                     c
!c							     c
!c                                                           c
!c   added the DAMP SUCCESSIVE BO MATRIX part to control     c
!c   convergence issues -MBH 2007-nov-09                     c
!c                                                           c
!c   moved the MAKE SURE MOS HAVE THE SAME SIGN loops        c
!c   out of the IF (CICALC) condition so that MOs            c
!c   are assured to have the same phase between time         c
!c   steps - MBH 2007-nov-09                                 c
!*  +---------------------------------------------------------------+
!*  |  make sure the MO coefs have the same sign as before          |
!*  |   will have to fix this for restart configurations, too ??    |
!*  +---------------------------------------------------------------+
!c                                                           c
!c-----------------------------------------------------------c

!* =====================================================================
!* []  DECK PPP -- DIAGONALIZE THE PPP FOCK MATRIX
!* =====================================================================

!*  +---------------------------------------------------------------+
!*  |  fock(i,j) is the fock matrix                                 |  
!*  |  hmat(i,j) is the one electron matrix                         |  
!*  +---------------------------------------------------------------+

!c-------------------------------------------------------------------
!c IDIM is the numebr of PI atoms
!c is defined in paramPPP.h and enter in definition of nci
!c-------------------------------------------------------------------

SUBROUTINE ppp_ppv
  
  USE commondata, ONLY:pi,npi,DP,nnuc,navg,d_screen,norb,a0,log_cicalc,avogad,&
       log_exdyn,istep,eppp,ebetag,eground,ecore,log_ppv,navg,log_addcos,istat,log_C60shift
  USE quantum
  USE PPV, ONLY: ntor_junc,junction,torsion
  USE classical, ONLY: nspec,pi_iden,atype,pos,dfield,length,fra,locate_atom,species
  USE qcfio, ONLY: gen_io,mo_io,cicoef_io,orb_io,bo_io

  IMPLICIT NONE

  !local variables
  LOGICAL :: converged
  INTEGER :: i,j,k,ii,jj,kk,i1,j1,ii1,ii2,i2,j2,itemp,itemp2,imo,numiter,&
       iter,u,v,idh,indx,iatom,n,m
  REAL(DP), DIMENSION(npi,npi) :: gamma,fock,hmat,beta,bomat_old,&
       dotm1,gamma_pi,mo_old
  !REAL(DP), DIMENSION(3,npi,npi,npi) :: qr
  !REAL(DP), DIMENSION(npi,npi,npi,npi) :: utmp
  REAL(DP), DIMENSION(nnuc,nnuc) :: gamma_nuc, gamma_nuc_cos
  REAL(DP), DIMENSION(npi) :: ener,vtemp
  !REAL(DP), DIMENSION(ntor_junc) :: facppp
  !REAL(DP), DIMENSION(nspec) :: gamma_ppv
  REAL(DP), DIMENSION(npi*3) :: vtest
  REAL(DP), DIMENSION(3) :: dr,ftemp,ftmp,for_temp
  REAL(DP) :: rij,rij2,aij,rmse,told,enersum,enersumold,totener,&
       dte,dte1,dtbo1,dtbo,damp,e_qm_dip,trans,angle,chgion,&
       ebetaion,enb,fac,forcefac,scdfield,enerex,suma,sumb,sumn,&
       efactor
  !ci variables
  INTEGER, DIMENSION(norb*(npi-norb)) :: cib,cit
  INTEGER :: idum,jdum,lwork,iee,ierr
  REAL(DP), DIMENSION(nci) :: vtemp2,ene_gap
  REAL(DP), DIMENSION(nci,nci) :: cimat
  REAL(DP), DIMENSION(nci,nci) :: ci_old
  REAL(DP), DIMENSION(neig,neig) :: dotm2
  REAL(DP), DIMENSION(norb*(npi-norb)) :: cimatdiag
  !REAL(DP), DIMENSION(npi*(npi-1)/2,npi*(npi-1)/2) :: aeainv
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: work,iwork,ifail
  double precision DLAMCH
  REAL(DP) :: b_i1,b_i2, mu_b, r0, r12, gs, tmp,tmp11,tmp22,tmp12,timei,timef,sum
  CHARACTER(1) :: spec1, spec2
  REAL(DP), DIMENSION(npi,npi) :: cin,cio
  external DLAMCH
  
  !APW timer
  integer clock_max,clock_rate,counti,count1,count2,count3
  call system_clock ( counti, clock_rate, clock_max )
  !*  +---------------------------------------------------------------+
  !*  |  save the MO coefs from the previous time step                |
  !*  |   will need to be able to read this in for restart configs    |
  !*  |   will it come from bo_restart.dat?                           |
  !*  +---------------------------------------------------------------+
     
  if (navg.ne.0) then
     call dcopy(npi*npi,mocoef,1,mocoef_old,1)
     call dcopy(npi*npi,mocoef,1,mo_old,1)
  end if

  efactor=2625.5 !APW Hartree to kJ/mol
  
  !*  +---------------------------------------------------------------+
  !*  |  set up the gamma params which are independant of bond orders |  
  !*  |  note that the offdiag. gammas depend on distances between    |  
  !*  |  atoms                                                        |  
  !*  +---------------------------------------------------------------+
  !call gammaset_ppv(gamma,nap_ppv,nc_ppv)
  !* =====================================================================
  !* []  DECK GAMMASET -- ZEROTH ORDER APPROXIMATION TO FOCK MATRIX 
  !*                      ELEMENTS                                 
  !* =====================================================================
  !*  +---------------------------------------------------------------+  
  !*  |  Mataga-Nishimoto relationship for 2-e repulsion integral     |  
  !*  |                                                               |  
  !*  |  Y_uv = e^2 / (R_uv + a_uv)                                   |  
  !*  |  a_uv = 2*e^2 / (Y_uu + Y_vv)                                 |  
  !*  |                                                               |  
  !*  |  JL parameter is       gamma=0.409011 Hartree                 |
  !*  |  Wharshel parameter I-A is 0.360476 Hartree                   |
  !*  |  Wharshel parameter G0  is G0 0.188933 Hartree                |
  !*  |  Wharshel first terms (I-A)+G0 = 0.548053 Hartree             |
  !*  |                                                               |
  !*  | PPV needs two gammas if using delRe charge sum since          |
  !*  | H(n,n)=sum_(n.neq.r) Z_r gamma_nr is extended to no-PI site   |
  !*  | the extras are the gamma_ppvs, naturally, which have terms    |
  !*  | for all atoms in the solute, including non-pi-C and H atoms   |
  !*  |                                                               |
  !*  |  (If log_addcos=.T.)  AY added  2011.01                       |  
  !*  |  The terms which depend on cos(torsion angle) are added       |
  !*  |  to Gamma, according to the original PI/QCFF by Warshel.      |
  !*  +---------------------------------------------------------------+
  !SUBROUTINE GAMMASET_PPV(gamma,npi,nc_ppv)  

  gamma_nuc  = 0.d0    !matrix
  gamma_pi = 0.d0    !matrix
  
  !*  +---------------------------------------------------------------+
  !*  |  do diagonal gammas first                                     |
  !*  +---------------------------------------------------------------+
  
  !APW gamma
  do i=1,npi              !just over the pi-C atoms
     ii=pi_iden(i)
     !gamma_pi(i,i)=gamma_screen(atype(pi_iden(i)))/d_screen !dielectric screening MJBH Oct.2009
     gamma_pi(i,i)=gamma_screen_i(ii)/d_screen !dielectric screening MJBH Oct.2009
  end do

  do i=1,nnuc          !over all atoms in polymer, C and H
     !gamma_nuc(i,i)=gamma_screen(atype(i))/d_screen !dielectric screening MJBH Oct.2009
     gamma_nuc(i,i)=gamma_screen_i(i)/d_screen !dielectric screening MJBH Oct.2009
  end do
  
  !*  +---------------------------------------------------------------+
  !*  |  do the rest of the gammas (off-diagonals)                    |
  !*  +---------------------------------------------------------------+
  !c-------------------------------------------------------------------
  !c use nc_ppv array to convert PI index to molecular index
  !c-------------------------------------------------------------------
  
  !APW as this is written it applies to unitls of length of bohr radii
  do i=1,npi             !just over the pi-carbon atoms...
     do j = 1,npi        
        if (i.ne.j) then 
           ii=pi_iden(i)
           jj=pi_iden(j)
           aij = 2.0d0/(gamma_pi(i,i)+gamma_pi(j,j))
           dr=(pos(:,ii)-pos(:,jj))/a0 !APW to convert to bohr radii
           rij=dsqrt(length(dr))
           gamma_pi(i,j) = 1.0d0/(d_screen*(rij + aij))   !dielectric screening MJBH Oct.2009
        end if
     end do
  end do
  
  do i=1,nnuc         !over all atoms in the polymer...
     do j = 1,nnuc    !only used when delRe charges
        if (i.ne.j) then !  are enabled?
           aij = 2.0d0/(gamma_nuc(i,i)+gamma_nuc(j,j))
           dr=(pos(:,i)-pos(:,j))/a0 !APW to convert to bohr radii
           rij=dsqrt(length(dr))
           gamma_nuc(i,j) = 1.0d0/(d_screen*(rij + aij)) !dielectric screening MJBH Oct.2009
        end if
     end do
  end do
    

  !APW call betalin here so I don't have to compute it many times per timestep
  do i= 1,ibeta_pairs
     i1 = ibeta(1,i)
     j1 = ibeta(2,i)
     
     angle = 1.0d0
     if(log_ppv) then
        if(junc_flg(i1).and.junc_flg(j1)) then
           idh=inv_dih(i1,j1)
           if(idh.eq.0) then
              write(gen_io,*) '***Warning: error in subroutine ppp_ppv.f90, idh=0***'
              stop
           end if
           angle = torsion(idh)
        end if
     end if
     dr       = (pos(:,i1)-pos(:,j1))/a0
     ftemp(1) = dsqrt(length(dr))
     ftemp(2) = betalin(ftemp(1),species(atype(i1)),species(atype(j1)))
     ! (for beta)
     st_betalin(i,1) = ftemp(2)
     st_betalin(i,2) = angle
     st_fbetalin(i,:)= fbetalin(dr,species(atype(i1)),species(atype(j1)),ftemp(2),angle)
     ! (for alpha and gamma)   AY added
     if(log_ppv .and. log_addcos) then
     st_exp2mu(i,1) = exp2mu(ftemp(1),species(atype(i1)),species(atype(j1)))
     st_exp2mu(i,2) = angle  ! cos(phi)
!     st_exp2mu(i,2) = 1.0d0   ! test
     endif
  end do

  !AY  Terms including cos(phi) are added  
  gamma_nuc_cos(:,:)=0.0d0    
  if(log_addcos.and.log_ppv) call add_cos_term_to_gamma(gamma_nuc,gamma_pi,gamma_nuc_cos)

  call setup_ppp(fock,hmat,beta,gamma_nuc)
  
  !*  +---------------------------------------------------------------+
  !*  |  diagnolize the matrix, MO eigenvalues in array ener          |
  !*  +---------------------------------------------------------------+
  
  !APW compute selected eigenvalues of matrix
  call dcopy(npi*npi,fock,1,mocoef,1)
  !mocoef=fock
  
  itemp = 3*npi

  call dsyev('V','U',npi,mocoef,npi,ener,vtest,itemp,itemp2)
  
  if(itemp2.ne.0) then
     write(gen_io,*) '***warning: error in dsyev: itemp2=',itemp2
     stop
  end if

  !*  +---------------------------------------------------------------+
  !*  |  recalculate the bond-order matrix                            |  
  !*  |  the mos matrix is subdivided as                              |  
  !*  |                                                               |  
  !*  |  1st MO is given by mos(i = 1,...,N,1)                        |  
  !*  |  2nd MO is given by mos(i = 1,...,N,2)                        |  
  !*  |  etc. for N-electron system (where each electron comes from   |
  !*  |  a single, specific pi-atom                                   |
  !*  |                                                               |  
  !*  |  said another way, the array is mos(atom #,state #)           |  
  !*  |                                                               | 
  !*  |      bo(u,v) = P_uv = 2 * sum(i,N/2) C_u(i) C_v(i)^*          | 
  !*  |                                                               | 
  !*  |  said another way,                                            |  
  !*  |                                                               | 
  !*  |      bo(u,v) = P_uv = 2*sum(i)^N/2 mos(u,i)*mos(v,i)          | 
  !*  |                                                               | 
  !*  +---------------------------------------------------------------+

  call dgemm('N','T',npi,npi,norb,2.0d0,mocoef,npi,mocoef,npi,0.0d0,bomat,npi)
 
  !*  +---------------------------------------------------------------+
  !*  |  setup the fock and one e matrices with the Bond Order matrix |  
  !*  |  we use the bond order matrix from the previous time step     | 
  !*  |  if it's the first time step then it's from the config file   |  
  !*  +---------------------------------------------------------------+

  call setup_ppp(fock,hmat,beta,gamma_nuc)
  
  write(gen_io,*) '-------------------------------------------------'
  write(gen_io,*) 'WE START ITERATION FOR MOLECULAR ORBITALS'
  write(gen_io,*) '-------------------------------------------------'

  numiter = 400    !max number of iterations
  enersum = 0.0

  converged=.false.
  iter=1

  do while(.not.converged)
     !*  +---------------------------------------------------------------+
     !*  |  perform level shifting as described in                       |
     !*  |    K. Mehrotra, Theoret. Chim. Acta (Berl.) 46, 325-9 (1977)  |
     !*  +---------------------------------------------------------------+
     !c       lvl_shft=0.d
     !c       lvl_shft = 1.d-4
     !c       if (iter .gt. 180) lvl_shft = 2.d-2
     !c
     !c       if (iter .gt. 1) call dcopy(idim*idim,mos,1,tmpmos,1) 
     !c       if (iter .eq. 1) then
     !c         write(*,*)'using level shift ',lvl_shft
     !c         do i = 1,idim
     !c          do j = 1,idim
     !c            tmpmos(i,j) = 0.d0
     !c         enddo
     !c         enddo
     !c       endif
     !c
     !c       do i = 1,nap_ppv
     !c        do j = 1,nap_ppv
     !c          delta(i,j) = 0.d0
     !c        enddo
     !c       enddo
     !c
     !ccc     shift virtual orbitals up
     !c       do i = nap_ppv/2+1,nap_ppv
     !c          delta(i,i) = lvl_shft
     !c       enddo
     !c
     !ccc     shift occupied orbitals down
     !c        do i = 1,nap_ppv/2
     !c           delta(i,i) = -lvl_shft
     !c        enddo
     !c
     !c        call dgemm('n','n',nap_ppv,nap_ppv,nap_ppv,1.d0,delta,
     !c    x             nap_ppv,tmpmos,nap_ppv,0.d0,delta1,nap_ppv)
     !c        call dgemm('t','n',nap_ppv,nap_ppv,nap_ppv,1.d0,tmpmos,
     !c    x             nap_ppv,delta1,nap_ppv,0.d0,delta,nap_ppv)
     !c        do i = 1,nap_ppv
     !c         fock(i,i) = fock(i,i)+delta(i,i)
     !c        enddo

     !*  +---------------------------------------------------------------+
     !*  |  diagnolize the matrix                                        |
     !*  +---------------------------------------------------------------+

     

     call dcopy(npi*npi,fock,1,mocoef,1) 
     !mocoef=fock
!!$
!!$ 1122  format(400(f6.3,2x))
!!$
     itemp=3*npi
     call dsyev('V','U',npi,mocoef,npi,ener,vtest,itemp,itemp2)
     
     if (itemp2 .ne. 0) then
        write(gen_io,*) 'ERROR in ppp_ppv:dsyev; itemp2 = ', itemp2
     end if

     !*  +---------------------------------------------------------------+
     !*  |  recalc the BO matrix                                         |
     !*  +---------------------------------------------------------------+

     !bomat_old=bomat
     call dcopy(npi*npi,bomat,1,bomat_old,1)
     
     call dgemm('N','T',npi,npi,norb,2.0d0,mocoef,npi,mocoef,npi,0.0d0,bomat,npi)
     
     !*  +---------------------------------------------------------------+
     !*  |  Dampen successive BO matricies to kill oscillations          |
     !*  |  before re-calculating fock matrix                            |
     !*  +---------------------------------------------------------------+

!!$c         DAMP = 1.0d0
!!$c         DAMP = 0.8d0
!!$c         DAMP = 0.5d0
!!$c         DAMP = 0.05d0
!!$c         damp = (.9d0/dble(iter) + .1d0)
!!$c         damp = (.49d0/dble(iter) + .05d0)
!!$c         damp = 1.d0/(1.d0+dble(iter))
     !damp = -.7d0/dble(iter) + .9d0
     !APW change back
     damp = (numiter-iter)/dble(numiter)
     rmse = 1.0d-12
     told = 2.0d-7   

!!$c         dfactor = dble(iter-100)/60.d0
!!$c         dfactor = dble(iter)/30.d1
!!$c         damp = .99d0 - .95d0*exp(-dfactor**4)
!!$c         damp = .60d0 - .55d0*exp(-dfactor**4)
!!$c         if (istep .ge. 1) damp = -.7d0/dble(iter) + .9d0
!!$c         damp = (dble(iter)/(numiter+10))**4+.01
!!$c         damp = .55d0*exp(-(dble((iter-1))/100)**4)+.05d0
!!$c         if (istep .ge. 1) damp = 1.d0

     do i = 1,npi
        do j = 1,npi
           bomat(i,j) = damp*bomat(i,j) + (1.d0-damp)*bomat_old(i,j)
        end do
     end do
     if (iter .eq. 1) then
        write(gen_io,*)' using damping coef ',real(damp), real(rmsE),real(tolD)
     end if

     !*  +---------------------------------------------------------------+
     !*  |  setup the fock and one e matrices with the Bond Order matrix |
     !*  +---------------------------------------------------------------+
     call setup_ppp(fock,hmat,beta,gamma_nuc)
          
     enersumold = enersum

     !*  +---------------------------------------------------------------+
     !*  |  calculate the total energy                                   |  
     !*  |                                                               |  
     !*  |  E_0 = sum(u) sum(v) 0.5 * P_uv * (H_uv + F_uv)               | 
     !*  |      = Tr{ P (H + F)^*}                                       |  
     !*  +---------------------------------------------------------------+

     totener = 0.0
     do i = 1,npi
        do j = 1,npi
           totener = totener + bomat(i,j)*(hmat(i,j)+fock(i,j))
        end do
     end do

     totener = totener*0.5
     enersum = totener

!AY     dtE1 = dtE 
     dtE = enersumold-enersum 

     !*  +---------------------------------------------------------------+
     !*  |  Find largest difference between bond-order matrices          |
     !*  +---------------------------------------------------------------+
!AY     dtBO1 = dtBO
     dtBO = 0.0d0

     do j = 1, npi
        do i = 1, npi 
           ftemp(1) = dabs( bomat(i,j) - bomat_old(i,j) )
           dtBO = max(dtBO,ftemp(1))
        end do
     end do

     !*  +---------------------------------------------------------------+
     !*  |  Find deviation between successive density matricies          |
     !*  +---------------------------------------------------------------+
     dtBO1 = 0.d0
     do j = 1,npi
        do i = 1,npi
           ftemp(1) = bomat(i,j) - bomat_old(i,j)
           dtBO1 = dtBO1 + ftemp(1)**2
        end do
     end do
     dtBO1 = (dtBO1/(dble(npi**2)))**.5

!!$*  +---------------------------------------------------------------+
!!$*  | Here is where we check if the shit is hitting the fan or not  |
!!$*  +---------------------------------------------------------------+

     !if (mod(iter,10) .eq. 0) then
     !   write(*,6363)istep, iter, dtE, dtBO, enersum,damp!,lvl_shft, enersum
     !endif
     ! 6363 format(2(i5,1x),13(f16.8,1x))
     
     if(dabs(dtE) .lt. rmsE .and. dtBO .lt. tolD) converged=.true. 

     !*  +---------------------------------------------------------------+
     !*  |  Here is where the shit is hitting the fan...                 |
     !*  +---------------------------------------------------------------+

     if (iter.gt.numiter) then
        write(gen_io,*) '***warning: possible trouble***'
        write(gen_io,*) iter,' iteractions exceeded'
        write(gen_io,*) "enersum:enersumold=",enersum,":",enersumold
        write(gen_io,*) 'max difference in bo is dtbo',dtbo
        write(gen_io,*) dabs(dtE),rmsE,dtBO,tolD
        !if (enersum.gt.enersumold) converged=.true.
        stop
     end if
     iter=iter+1
  end do

  write(gen_io,*) 'number of SCF iterations =',iter
  write(gen_io,*) 'enersum and enersumold is ',enersum,enersumold
  WRITE(gen_io,*) '------------------------------------------------'
  write(gen_io,*) dtE, dtBO   ! AY changed
!  write(gen_io,*) dtE, dtE1, dtBO
  write(gen_io,*) 'converged = ', converged
  WRITE(gen_io,*) '------------------------------------------------'

  !*  +---------------------------------------------------------------+
  !*  |  Make sure the final MO coefs have the same sign as before    |
  !*  +---------------------------------------------------------------+

  !c         *****************************************************
  !c         ** goal of this next section will be to ensure the **
  !c         ** consistency of the MOs between time stpes such  **
  !c         ** that MO_i(t)*MO_i(t-dt) ~ 1  ... problems may   **
  !c         ** arrise during trajectories as an MO may switch  **
  !c         ** identities with (become lower in energy than) a **
  !c         ** neighboring MO, j.  Also, the MOs phases tend to**
  !c         ** oscillate, and we might see that the max overlap**
  !c         ** between timesteps is MO_i(t)*MO_j(t-dt) ~ -1    **
  !c         ** this stuff fixes all that...                    **
  !c         ** this is important if you plan on taking finite  **
  !c         ** differences (or doing other math) with the MOs  **
  !c         *****************************************************


  !c         ** save prev. mo moefs on 1st step, else calc dot  **
  !c         ** prod. for each MO with each of the previous MOs **
  if(navg.eq.0) then 
     call dcopy(npi*npi,mocoef,1,mocoef_old,1)   
  else
     !call dgemm('T','N',npi,npi,npi,1.0d0,mocoef,npi,mocoef_old,npi,0.0d0,dotm1,npi)
     call dgemm('T','N',npi,npi,npi,1.0d0,mocoef,npi,mo_old,npi,0.0d0,dotm1,npi)

     !c ** find the largest dotm1 to look for switched identity **
     do i =1,npi
        ftemp(1) = 0.0d0
        do k = 1,npi
           if (dabs(dotm1(i,k)) .gt. ftemp(1)) then 
              itemp = k
              ftemp(1) = dabs(dotm1(i,k))
           end if
        end do
        
        !c               ** rename the oldmos array for max over- **
        !c               ** lap if an identity has been switched  **
        if (itemp.ne.i) then
           
           vtemp=mo_old(:,itemp)
           call dcopy(npi,mo_old(1,itemp),1,vtemp,1) 
           call dcopy(npi,mo_old(1,i),1,mo_old(1,itemp),1)  
           call dcopy(npi,vtemp,1,mo_old(1,i),1)  
           
           vtemp=mocoef_old(:,itemp)
           call dcopy(npi,mocoef_old(1,itemp),1,vtemp,1) 
           call dcopy(npi,mocoef_old(1,i),1,mocoef_old(1,itemp),1)  
           call dcopy(npi,vtemp,1,mocoef_old(1,i),1)  

           !APW must reorder the CI indices as well
           if(log_cicalc) then
              do jdum=1,nci
                 if(ciind_old(jdum,1).eq.i) then
                    ciind_old(jdum,1)=itemp
                 else
                    if(ciind_old(jdum,1).eq.itemp) then
                       ciind_old(jdum,1)=i
                    end if
                 end if
                 if(ciind_old(jdum,2).eq.i) then
                    ciind_old(jdum,2)=itemp
                 else
                    if(ciind_old(jdum,2).eq.itemp) then
                       ciind_old(jdum,2)=i
                    end if
                 end if
              end do
           end if
           !c** recalc the dotm1 array **
           !call dgemm('T','N',npi,npi,npi,1.0d0,mocoef,npi,mocoef_old,npi,0.0d0,dotm1,npi)
           call dgemm('T','N',npi,npi,npi,1.0d0,mocoef,npi,mo_old,npi,0.0d0,dotm1,npi)
        end if
     end do
     !c** recalc the dotm1 array after all identity switches **
    !call dgemm('T','N',npi,npi,npi,1.0d0,mocoef,npi,mocoef_old,npi,0.0d0,dotm1,npi)
     call dgemm('T','N',npi,npi,npi,1.0d0,mocoef,npi,mo_old,npi,0.0d0,dotm1,npi)
     

     !c            ** ensure the phase of the MOs are the same for  **
     !c            ** each time step so that MO_i(t)*MO_i(t-dt) ~ 1 **
     do i =1,npi
        if (dotm1(i,i) .lt. 0.0d0) then
           !mocoef(:,i)=-mocoef(:,i)
           call dscal(npi,-1.0d0,mocoef(1,i),1) 
        end if
        !c         ** and save the MOs for next time step in oldmos **

        !call dcopy(npi,mocoef(1,i),1,mocoef_old(1,i),1)  
        call dcopy(npi,mocoef(1,i),1,mo_old(1,i),1)  
     end do
     
     !if(mod(istep+1,irestart).eq.0 .and. istep .gt. 0) THEN
!!$              write(93,*)' '
!!$              write(93,*)"MO COEFFICIENTS istep =",istep
!!$              do i = 1,nap_ppv
!!$                write(93,'(i5,1x,i5,500(1x,f16.8))')i,nc_ppv(i)
!!$     +                    ,(mos(i,j),j=1,nap_ppv)
!!$              enddo
!!$             ENDIF
  end if
  
  !*  +---------------------------------------------------------------+
  !*  |  recalc the BO matrix with the last iteration                 |
  !*  +---------------------------------------------------------------+
  call dgemm('N','T',npi,npi,norb,2.0d0,mocoef,npi,mocoef,npi,0.0d0,bomat,npi)

  !*  +---------------------------------------------------------------+
  !*  |  setup the fock matrix elements and one electron energies     |
  !*  +---------------------------------------------------------------+
  call setup_ppp(fock,hmat,beta,gamma_nuc)
  !call setupPPP_PPV(fock,beta,hmat,bo,nap_ppv,ibeta,nc_ppv,navg,mos,999,istep)
  !*  +---------------------------------------------------------------+
  !*  | calculate \langle \Psi_0 |(-dfield*\hat{x})|\Psi_0\rangle     |
  !*  +---------------------------------------------------------------+

  !do i= 1,nnuc   !JJBL
  !rposx(i) = xna(i)*unitl/0.529177   !converts pos into Bohr
  !rposy(i) = yna(i)*unitl/0.529177   !from internal length
  !rposz(i) = zna(i)*unitl/0.529177
  !enddo

  e_qm_dip = 0.d0 !JJBL
  do i = 1,npi
     ii = pi_iden(i)
     e_qm_dip = e_qm_dip - dfield * bomat(i,i) * pos(1,ii)/a0
  end do
  write(gen_io,*)'e_qm_dip = ',e_qm_dip,' Hartree'

  !*  +---------------------------------------------------------------+
  !*  |  calculate the total energy                                   |
  !*  |                                                               |  
  !*  |  E_0 = sum(u) sum(v) 0.5 * P_uv * (H_uv + F_uv)               | 
  !*  |      = Tr{ P (H + F)^*}                                       |  
  !*  +---------------------------------------------------------------+

  totener = 0.0
  do i = 1,npi
     do j = 1,npi
        totener = totener + bomat(i,j)*(hmat(i,j)+fock(i,j))
     end do
  end do

  totener = totener*0.5

  !*  +---------------------------------------------------------------+
  !*  |  the result is in Hartress so we convert to program units     |  
  !*  | fkjmol converts FROM DIMENSIONLESS ENERGY to kJ/mol           |
  !*  |                                                               |
  !*  |  epppkjmol = SCF ground state energy + CI eigenvalue of the   |  
  !*  |  currently occupied state                                     |  
  !*  +---------------------------------------------------------------+
  !check molcs che cosa e' se e' numero di solvente va cambiato -wtf?
  !APW translated check molcs to see if number of solvent is changed

  !do i = 1,npi
  !enerPPPkjmol(i) = ener(i)*4.3598d-18*1.0d-3/dfloat(molcs)*avogad !Hartrees to kJ/mol
  !enerPPPprog(i) = enerPPPkjmol(i)/fkjmol  !kj/mol -> program
  !end do
  
  eppp    = totener * efactor
  eground = eppp
  !ePPPkjmol = totener*4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !ePPPprog = ePPPkjmol/fkjmol
  !ene_ground= ePPPkjmol 

  !* =====================================================================
  !*     CALCULATE SINGLET TRANSITION ENERGIES                          
  !*                                  
  !*     don't use all the CI states since there are too many  
  !*     1.  Calculate the diagonal elements of the CI matrix
  !*     2.  Sort the diagonal elements according to their values  
  !*     3.  Calculate the off-diagonal elements of the CI matrix
  !*     4.  Diagnolize that matrix                        
  !*     5.  Keep the lowest NCI states 
  !* =====================================================================
  !*  +---------------------------------------------------------------+
  !*  | CICALC is defined in IN.PPV                                   |
  !*  +---------------------------------------------------------------+
  !*  |  save CI energies from the previous time step in order to     |  
  !*  |  integrate them; note that on the 0th time step cienerold( )  |  
  !*  |  will contain garbage since we haven't calculated the CI      |  
  !*  |  eigenvalues yet but since they are not used on the 0th time  |  
  !*  |  step it doesn't matter                                       |  
  !*  |                                                               |  
  !*  |                  dont forget that loop = nap_ppv              |
  !*  +---------------------------------------------------------------+
  !call system_clock ( count1, clock_rate, clock_max )
  !print *,'ppp_ppv.f90: MO',(count1-counti)/real(clock_rate)
  
  if(log_C60shift) then
     itemp=0
     do j=1,npi
        ftemp(1)=0.0d0
        do i=(npi-60),npi
           ftemp(1)=ftemp(1)+mocoef(i,j)*mocoef(i,j)
        end do
        if(ftemp(1) .gt. 0.9) then
           itemp=itemp+1
           if((itemp .gt. 33).and.(itemp .le. 36)) then
              ener(j)=ener(j)+0.022
           end if
        end if
     
        !ftemp(1)=0.0d0
        !do i=(npi-120),npi-60
        !   ftemp(1)=ftemp(1)+mocoef(i,j)*mocoef(i,j)
        !end do
        !if(ftemp(1) .gt. 0.9) then
        !   itemp=itemp+1
        !   if((itemp .gt. 33).and.(itemp .le. 36)) then
        !      ener(j)=ener(j)+0.022
        !   end if
        !end if
     end do
  end if
  
  if(log_cicalc) then
     WRITE(gen_io,*) '-------------------------------------------------'
     WRITE(gen_io,*) 'WE START ITERATION FOR C.I. ORBITALS'
     WRITE(gen_io,*) '-------------------------------------------------'

     lwork = 8*nci
     ALLOCATE(work(lwork),iwork(5*nci),ifail(norb*nci))

     !ciener_l=ciener
     call dcopy(nci,ciener,1,ciener_l,1)  
     
     
     
     !*  +---------------------------------------------------------------+
     !*  |  (1) set up diagonal elements of CIS matrix                   |
     !*  |                                                               |
     !*  |  < Psi(i->k) | H_el - V_pi | Psi(j->l) >                      |
     !*  |              = (E_k - E_i) delta(i,j) delta(k,l)              |
     !*  |              + sum(u;v) C_u^i C_v^k ( 2 C_u^k C_v^i           |
     !*  |                                       - C_u^i C_v^k) Y_uv     |
     !*  +---------------------------------------------------------------+

     jdum = 0 
     !APW cib and cit will always be the same, should they be global?
     do i=1,norb
        do k=norb+1,npi
           jdum=jdum+1
           cib(jdum)=i
           cit(jdum)=k
           cimatdiag(jdum) = ener(k)-ener(i)
           do u=1,npi
              do v=1,npi
                 cimatdiag(jdum)=cimatdiag(jdum)&
                      +(2.0*mocoef(u,k)*mocoef(u,i)*mocoef(v,i)*mocoef(v,k)&
                      -mocoef(u,k)*mocoef(u,k)*mocoef(v,i)*mocoef(v,i))&
                      *gamma_pi(u,v)
              end do
           end do
        end do
     end do
     
     !*  +---------------------------------------------------------------+
     !*  |(2) sort the array of diagonal CIS elements in ascending order |
     !*  +---------------------------------------------------------------+

     call sort2(jdum,cimatdiag,cit,cib)

     !c       do i = 1,jdum
     !c       if (ib(i) .eq. 57 .and. it(i) .eq. 61) 
     !c    x     write(*,*)"JJBL cimat(57,61)",i,cimatdiag(i),
     !c    x                       (ener(61)-ener(57))
     !c       enddo
     !c       stop

     !*  +---------------------------------------------------------------+
     !*  |  (3) calculate the off-diagonal elements:                     |  
     !*  |                                                               |
     !*  |  < Psi(i->k) | H_el - V_pi | Psi(j->l) >                      |
     !*  |              = sum(u;v) C_u^j C_v^k ( 2 C_u^l C_v^i           |
     !*  |              - C_u^i C_v^l) Y_uv                              |
     !*  +---------------------------------------------------------------+

     do jdum=1,nci
        cimat(jdum,jdum)=cimatdiag(jdum)
        ciind(jdum,1)=cib(jdum)
        ciind(jdum,2)=cit(jdum)
        do idum = jdum+1,nci
           cimat(idum,jdum)=0.0d0
           do u=1,npi
              do v=1,npi
                 cimat(idum,jdum)=cimat(idum,jdum)&
                      +(2.0*mocoef(u,cit(jdum))*mocoef(u,cib(jdum))&
                      *mocoef(v,cib(idum))*mocoef(v,cit(idum))&
                      -mocoef(u,cit(jdum))*mocoef(u,cit(idum))&
                      *mocoef(v,cib(jdum))*mocoef(v,cib(idum)))&
                      *gamma_pi(u,v)
              end do
           end do

           cimat(jdum,idum)=cimat(idum,jdum)
        end do
     end do
     
     !do i = 1,nci
     !write(*,1020)(cimat(j,i),j=1,nci)
     !enddo
     !1020 format(400(e12.5,2x))

     !*  +---------------------------------------------------------------+
     !*  | (4) Diagonalize CIS matrix                                    |
     !*  +---------------------------------------------------------------+
     !c         call dcopy(nci,cieig(1,ii),1,cieigold(1,ii),1)  !JJBL
     
     call dsyevx('V','I','U',nci,cimat,nci,1.0,1.0,1,neig,2*DLAMCH('S'),iee,ciener,cieig,nci,work,lwork,iwork,ifail,ierr)
     
     !*  +---------------------------------------------------------------+
     !*  | (5) Make sure the final CI coefs have the same sign as before |
     !*  |     also, make sure they have the maximum possible overlap    |
     !*  |     with a CI eigvec from prev. step in case they switch      |
     !*  |     exactly the same protocol is taken for the MO coefs       |
     !*  +---------------------------------------------------------------+

     if (navg.eq.0) then
        !cieig_old=cieig
        !call dcopy(nci*nci,cieig,1,cieig_old,1)  
        !call dcopy(nci*nci,cieig,1,ci_old,1)  
     else
        call dcopy(nci*nci,cieig_old,1,ci_old,1)
        call dgemm('T','N',neig,neig,nci,1.0d0,cieig,nci,ci_old,nci,0.0d0,dotm2,neig)
        
        do ii = 1,neig
           itemp=ii
           ftemp(1) = 0.0d0
           !APW make this more efficient
           cin=0.0d0
           do n=1,nci
              cin(ciind(n,1),ciind(n,2))=cieig(n,ii)
           end do
           do kk = 1,neig
              cio=0.0d0
              do n=1,nci
                 cio(ciind_old(n,1),ciind_old(n,2))=ci_old(n,kk)
              end do
              
              sum=0.0d0
              do n=1,norb
                 do m=norb+1,npi
                    sum=sum+cin(n,m)*cio(n,m)
                 end do
              end do
              
              !if (dabs(dotm2(ii,kk)) .gt. ftemp(1)) then  
              if (dabs(sum) .gt. ftemp(1)) then  
                 itemp = kk
                 ftemp(1) = dabs(dotm2(ii,kk))
              end if
           end do
           if (itemp.ne.ii) then
              !call dcopy(nci,cieig_old(1,itemp),1,vtemp2,1)  
              !call dcopy(nci,cieig_old(1,ii),1,cieig_old(1,itemp),1)  
              !call dcopy(nci,vtemp2,1,cieig_old(1,ii),1) 
              call dcopy(nci,ci_old(1,itemp),1,vtemp2,1)  
              call dcopy(nci,ci_old(1,ii),1,ci_old(1,itemp),1)  
              call dcopy(nci,vtemp2,1,ci_old(1,ii),1) 
              
              call dgemm('T','N',neig,neig,nci,1.0d0,cieig,nci,ci_old,nci,0.0d0,dotm2,neig)
           end if
        end do

        call dgemm('T','N',neig,neig,nci,1.0d0,cieig,nci,ci_old,nci,0.0d0,dotm2,neig)

        do ii =1,neig
           cin=0.0d0
           cio=0.0d0
           do n=1,nci
              cin(ciind(n,1),ciind(n,2))=cieig(n,ii)
              cio(ciind_old(n,1),ciind_old(n,2))=ci_old(n,ii)
           end do
           
           sum=0.0d0
           do n=1,norb
              do m=norb+1,npi
                 sum=sum+cin(n,m)*cio(n,m)
              end do
           end do
           
           !if (dotm2(ii,ii) .lt. 0.0d0) then
           if (sum .lt. 0.0d0) then
              call dscal(nci,-1.0d0,cieig(1,ii),1)
           end if
           call dcopy(nci,cieig(1,ii),1,ci_old(1,ii),1)  
        end do

        call dgemm('T','N',neig,neig,nci,1.0d0,cieig,nci,ci_old,nci,0.0d0,dotm2,neig)
     end if
     
     !*  +---------------------------------------------------------------+
     !*  |  write out 1st transition in nm                               |
     !*  +---------------------------------------------------------------+

     trans = ciener(1)
     write(99,*) 'ciener = ', trans 

     !*  +---------------------------------------------------------------+
     !*  |  convert hartrees to cm-1 to nm, write to fort.99             |
     !*  +---------------------------------------------------------------+
     trans = trans*219474.6
     write(99,*) 'Trans in cm-1 = ', trans    
     trans = 1.0 / trans * 1.0e+7
     write(99,*) 'S_0 -> S_1 transition = ', trans, ' nm'

     !c-------------------------------------------------------------------
     !c    Compute energy gap from ground state to each ex. states
     !c-------------------------------------------------------------------
     !APW commented out for now to determine why molcs=255
     !do i = 1,nci
     !   ene_gap(i)=ciener(i)*4.3598d-18*1.0d-3/dfloat(molcs)*avogad
     !end do

     !*  +---------------------------------------------------------------+
     !*  |  calculate the transition dipole                              |
     !*  +---------------------------------------------------------------+
  
     call transdoutput(cit,cib,cimatdiag)
     DEALLOCATE(work,iwork,ifail)
     
     if(navg.eq.0) then
        do i=1,nci
           ciind_old(i,1)=ciind(i,1)
           ciind_old(i,2)=ciind(i,2)
           do j=1,nci
              cieig_old(j,i)=cieig(j,i)
           end do
        end do
     end if
     
  end if  ! log_cicalc
  
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'ppp_ppv.f90: CI',(count2-count1)/real(clock_rate)
  !* =====================================================================
  !*     END OF CI CALCULATION                 
  !* =====================================================================
  !*     GROUND STATE FORCES  
  !* =====================================================================

  !*  +---------------------------------------------------------------+
  !*  |  electronic potential:                                        |
  !*  |                                                               |  
  !*  |      V_pi = 0.5 sum(u;v) P_uv (H_uv + F_uv)                   |  
  !*  |                                                               |
  !*  |               %%%%%%%%%%%%%%%%%%%%%%%%%%%                     |  
  !*  |                                                               | 
  !*  |      F_uv = H_uv - 0.5*P_uv*Y_uv      (u .ne.v)               |
  !*  |      F_uu = H_uu + 0.5*P_uu*Y_uu + sum(k .ne. u) P_kk *Y_uk   |
  !*  |      H_uu = alpha_u - sum(k .ne. u) Z_k*Y_uk                  |  
  !*  |                                                               |  
  !*  |               %%%%%%%%%%%%%%%%%%%%%%%%%%%                     |  
  !*  |                                                               |  
  !*  |      V_pi(u.ne.v) = sum(u .ne. v) P_uv H_uv                   |  
  !*  |                   - 0.5*0.5 sum(u .ne. v) P_uv P_uv Y_uv      |
  !*  |      V_pi(u.eq.v) = sum(u) P_uu H_uu                          |  
  !*  |                   - 0.5*0.5 sum(u) P_uu P_uu Y_uu             |  
  !*  |                   + 0.5 sum(u .ne. k) P_uu P_kk Y_uk          |  
  !*  |                                                               |
  !*  +---------------------------------------------------------------+
  !*  |  calculate force factors for TWIST                            |
  !*  +---------------------------------------------------------------+


  !*  +---------------------------------------------------------------+
  !*  |  Mataga-Nishimoto force:                                      | 
  !*  |  do it with my set of numberings first and translate back     |  
  !*  +---------------------------------------------------------------+

  forcefac=2625.5/a0 !from Hartree/bohr radius to kJ/(mol*angstrom)
  !forcefac = 3.2819590572086988E-003
  
  ebetaion = 0.0d0
  ebetag = 0.0d0
  
  !c----------------------------------------------------------------------
  !c
  !c  In case of MEH-PPV we separate loop over PI atom 
  !c  for matrix contribution to force 
  !c  to core-core electrostatic interaction. In this latter
  !c  we used loop over all atoms, because del re charge are on all atoms
  !c  First we loop over PI atoms. Then We loop Over NoPI atoms
  !c
  !c----------------------------------------------------------------------

  ! gradient terms of SCF energy (=0.5*P(H+F)) 
  ! except for the term involving beta
  do i =1,npi
     ii=pi_iden(i)
     do j = 1+i,npi

        fac=0.d0
        jj=pi_iden(j)

        dr=(pos(:,ii)-pos(:,jj))/a0
        rij=dsqrt(length(dr))

        fac = bomat(j,j)*bomat(i,i) - 0.5*bomat(i,j)**2

        fac = fac - z_atom(jj)*bomat(i,i)
        fac = fac - z_atom(ii)*bomat(j,j)


        !ftmp(:) =  fac*gamma_nuc(ii,jj)**2*dr(:)/rij      
        ftmp(:) =  fac*(gamma_nuc(ii,jj)-gamma_nuc_cos(ii,jj))**2*dr(:)/rij  !AY changed
        ftmp    =  ftmp * forcefac

        fra(:,ii) = fra(:,ii) + ftmp 
        fra(:,jj) = fra(:,jj) - ftmp  
     end do

     scdfield = dfield
     !scdfield = dfield * (1.d0-exp(-(dble(istep)/30.d0)**2))

     !force from both classical and qm interactions with E-field
     !fra(1,ii) = fra(1,ii) + forcefac*scdfield*(bomat(i,i)-1.d0)
     fra(1,ii) = fra(1,ii) + forcefac*scdfield*(bomat(i,i)-z_atom(ii))
  end do

  !c-------------------------------------------------------------------
  !c     Here I compute Zgamma taking into account internal loop over 
  !c     no PI atom use gamma_nuc for non-pi carbon atoms
  !c-------------------------------------------------------------------

  do i = 1,npi
     do j = 1,nnopi

        ii=pi_iden(i)
        jj=ind_nopi(j)

        dr=(pos(:,ii)-pos(:,jj))/a0
        rij=dsqrt(length(dr))

        fac = -z_atom(jj)*bomat(i,i)

        !-----------------------------------------------
        ! I use gamma_nuc because loop involve NO_PI atoms.
        !------------------------------------------------

        ftmp(:) = fac*gamma_nuc(ii,jj)**2*dr(:)/rij
        ftmp    = ftmp * forcefac

        fra(:,ii) = fra(:,ii) + ftmp
        fra(:,jj) = fra(:,jj) - ftmp  
     end do
  end do
  
  !*  +---------------------------------------------------------------+
  !*  |  force from core-core repulsion part of core-core             |  
  !*  |  potential:                                                   |  
  !*  |                                                               |  
  !*  |  V_c-c = sum(u;v) Y_uv Z_u Z_v + V_nb                         |  
  !*  +---------------------------------------------------------------+

  !*  +---------------------------------------------------------------+
  !*  |  get core-core repulsion energy and force                     |
  !*  +---------------------------------------------------------------+

  !write(62,*) 'repulsion force run over all atoms'
  do i=1,nnuc
     do j=i+1,nnuc
        dr=(pos(:,i)-pos(:,j))/a0
        rij=dsqrt(length(dr))

        chgion = z_atom(i)*z_atom(j)

        !APW initial value?
        ebetaion = ebetaion + chgion/rij
        ebetag = ebetag + chgion*gamma_nuc(i,j)

        !ftmp(:) =  chgion*gamma_nuc(i,j)**2*dr(:)/rij  !AY
        ftmp(:) =  chgion*(gamma_nuc(i,j)-gamma_nuc_cos(i,j))**2*dr(:)/rij  !AY changed
        ftmp    =  ftmp * forcefac

        fra(:,i) = fra(:,i) + ftmp
        fra(:,j) = fra(:,j) - ftmp
     end do
  end do

  !*  +---------------------------------------------------------------+
  !*  |  convert core-core energy to program units                    |  
  !*  +---------------------------------------------------------------+

  !APW this molcs again?
  !ebetaion = ebetaion*4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !ebetaion = ebetaion/fkjmol
  !ebetag = ebetag*4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !ebetag = ebetag/fkjmol
  ebetag = ebetag   * efactor
  ecore  = ebetaion * efactor

  !*  +---------------------------------------------------------------+
  !*  |  force from Lindberg approximation to the resonance integral  |
  !*  +---------------------------------------------------------------+
  !*  |  This is a delicate point. Torsion is value of tosion for     |
  !*  |  bond i1-i2 But in our case case, more torsions exist for the |
  !*  |  same bond. We defined teta as average over these possible    |
  !*  |  teta_i We use inv_dih: for bond index i1 and i2 give us the  |
  !*  |  index in the torsion list                                    |
  !*  +---------------------------------------------------------------+
  !AY  add force term including gradient betalin(:,1) in beta, 
  !    but does not add terms includiing gradient of cos(torsion angle) 
  !    in beta, here

  do i = 1, ibeta_pairs
     i1=ibeta(1,i)
     i2=ibeta(2,i)
     
     angle=1.0d0
     if(log_ppv) then
        if(junc_flg(i1).and.junc_flg(i2)) then
           idh=inv_dih(i1,i2)
           if(idh.eq.0) then
              write(gen_io,*)'***ERROR in subroutine ppp_ppv. idh=0***'
              stop
           end if
           angle=torsion(idh)
        end if
     end if

     dr       = (pos(:,i1)-pos(:,i2))/a0
     ftemp(1) = dsqrt(length(dr))
     ftemp(2) = betalin(ftemp(1),species(atype(i1)),species(atype(i2)))
     for_temp = fbetalin(dr,species(atype(i1)),species(atype(i2)),ftemp(2),angle)

     kk = inv_pi(ibeta(1,i))
     jj = inv_pi(ibeta(2,i))

     ftmp = 2.0d0 * bomat(kk,jj) * for_temp * forcefac

     fra(:,ibeta(1,i)) = fra(:,ibeta(1,i)) + ftmp
     fra(:,ibeta(2,i)) = fra(:,ibeta(2,i)) - ftmp
     
  end do


  
  !*  AY added
  !*  +------------------------------------------------------------------------+
  !*  |  force from the terms including cos(torsion angle) in alpha and gamma  |
  !*  |   (buta the force terms including dcos/dr are added later)             |
  !*  +------------------------------------------------------------------------+
  if(log_ppv .and. log_addcos) then

     do i= 1,ibeta_pairs
        !(atom number in all atoms)
        i1 = ibeta(1,i)
        i2 = ibeta(2,i)
        if( junc_flg(i1) .and. junc_flg(i2) ) then
           if(inv_dih(i1,i2)==0) then
              write(gen_io,*) '***Warning: error in ppp_ppv.f90 ***'
              stop
           end if

           !(atom number in pi atoms)
           ii1   = inv_pi(i1)
           ii2   = inv_pi(i2)

           spec1 = species(atype(i1))
           spec2 = species(atype(i2))
           call  get_mu_b_r0(spec1,spec2,mu_b,r0)

           dr(:)   = pos(:,i1)-pos(:,i2)  
           r12     = dsqrt(length(dr))
           tmp     = st_exp2mu(i,1) * st_exp2mu(i,2)*st_exp2mu(i,2) 
           tmp     = 2.0d0 * mu_b * tmp
           ftmp(:) = tmp * dr(:)/r12 * forcefac

           !=== Add force from the terms in Alpha in SCF energy gradient===
           call  get_beta_dash(spec1,b_i1)
           call  get_beta_dash(spec2,b_i2)
           tmp11   = b_i1 * bomat(ii1,ii1)
           tmp22   = b_i2 * bomat(ii2,ii2)

           ! d(term between i1 and i2 in alpha_i1)/dr_i1
           fra(:,i1) = fra(:,i1) + tmp11*ftmp(:)
           ! d(term between i1 and i2 in alpha_i1)/dr_i2
           fra(:,i2) = fra(:,i2) - tmp11*ftmp(:)
           ! d(term between i1 and i2 in alpha_i2)/dr_i1
           fra(:,i1) = fra(:,i1) + tmp22*ftmp(:)
           ! d(term between i1 and i2 in alpha_i2)/dr_i2
           fra(:,i2) = fra(:,i2) - tmp22*ftmp(:)


           !=== Add force from the terms including Gamma 
           !    in gradient from SCF energy and nuclear core-core replusion ===
           call  get_gs(spec1,spec2,gs)
           tmp11   = 0.25d0 * gs * bomat(ii1,ii1)*bomat(ii1,ii1)
           tmp22   = 0.25d0 * gs * bomat(ii2,ii2)*bomat(ii2,ii2)

           ! d(term between i1 and i2 in gamma(i1,i1))/dr_i1
           fra(:,i1) = fra(:,i1) + tmp11*ftmp(:)
           ! d(term between i1 and i2 in gamma(i1,i1))/dr_i2
           fra(:,i2) = fra(:,i2) - tmp11*ftmp(:)
           ! d(term between i1 and i2 in gamma(i2,i2))/dr_i1
           fra(:,i1) = fra(:,i1) + tmp22*ftmp(:)
           ! d(term between i1 and i2 in gamma(i2,i2))/dr_i2
           fra(:,i2) = fra(:,i2) - tmp22*ftmp(:)


           tmp11   = - bomat(ii1,ii1) * z_atom(i2)                &
                     + 0.5d0 * bomat(ii2,ii2) * bomat(ii1,ii1)
           tmp22   = - bomat(ii2,ii2) * z_atom(i1)                &
                     + 0.5d0 * bomat(ii1,ii1) * bomat(ii2,ii2)
           tmp12   = - 0.5d0 * bomat(ii1,ii2) * bomat(ii1,ii2)    &
                     + z_atom(i1) * z_atom(i2)

           tmp11   = - gs * tmp11
           tmp22   = - gs * tmp22
           tmp12   = - gs * tmp12

           ! d(term in gamma(i1,i2))/dr_i1 + d(term in gamma(i2,i1))/dr_i1
           fra(:,i1) = fra(:,i1) + ( tmp11+tmp22+tmp12 ) * ftmp(:)
           ! d(term in gamma(i1,i2))/dr_i2 + d(term in gamma(i2,i1))/dr_i2
           fra(:,i2) = fra(:,i2) - ( tmp11+tmp22+tmp12 ) * ftmp(:)

        end if
     enddo    ! i=1,ibeta_pairs

  endif

  
  !*  +---------------------------------------------------------------+
  !*  |  force from nonbonded interaction part of core-core           |
  !*  |  potential:                                                   |
  !*  |                                                               |
  !*  |  V_c-c = sum(u;v) Y_uv Z_u Z_v + V_nb                         |
  !*  +---------------------------------------------------------------+

  enb = 0.0
  indx = 0

  !* ====================================================================
  !*     EXCITED STATE FORCES                                          
  !* ====================================================================
  !*  +---------------------------------------------------------------+
  !*  |  at this point we have already added in the ground state      |  
  !*  |  forces; we just need to add in the excited state forces      |  
  !*  |  relative to the ground state                                 |  
  !*  +---------------------------------------------------------------+
  call system_clock ( count3, clock_rate, clock_max )
  !print *,'ppp_ppv.f90: AY',(count3-count2)/real(clock_rate)
  if (log_exdyn) then 

     write(gen_io,*) '-------------------------------------------------'
     write(gen_io,*) 'Starting Computation Forces for exited  dynamics '
     write(gen_io,*) '-------------------------------------------------'
     
     !print *,'in solvuary'
     call solvuary(gamma_nuc,gamma_pi,ener)
     write(gen_io,*) '---> solvuary'
     !call system_clock ( count1, clock_rate, clock_max )
     !print *,'solvuary.f90',(count1-count3)/real(clock_rate)
     
     !print *,'in add_force_ex'
     call add_force_ex(cit,cib,gamma_pi)
     write(gen_io,*) '---> excited force on solute'
     !call system_clock ( count2, clock_rate, clock_max )
     !print *,'add_force_ex.f90',(count2-count1)/real(clock_rate)
     enerex = 0.0d0

     !*  +---------------------------------------------------------------+
     !*  |  use the CI eigenvalue to see if this improves energy         |  
     !*  |  conservation; iexstat tells us which state we are in         |  
     !*  +---------------------------------------------------------------+
     if (iexstate.ne.0) then
        
        enerex = ciener(iexstate)
        
        do iatom = 1,nnuc
           fra(:,iatom) = fra(:,iatom) + fexa(:,iatom,iexstate,iexstate)
        end do
        
        eppp = eppp+enerex*efactor
        
     end if
  end if !if(exdyn)
  
  !**  write positions in angstrom of pi-carbon atoms   **
  !**  and ground-state MO coefs to MO_coefs.dat for    **
  !**  later use with plot-2p.f, ncgen, and the Chimera **
  !**  plotting software to see the MOs & chg. density  **
  !**                                                   **
  !**  output in CI_coefs.dat is:                       **
  !**      col 1: ex. state #                           **
  !**      col 2: CIS excitation #                      **_
  !**      col 3: remove electron from this GS MO       ** \ describes the
  !**      col 4: place electron in this virtual MO     **_/ virt. excitation
  !**      col 5: CI coeff: state (col 1),CIS ex (col 2)**_
  !**      col 6: says either 'large' or 'small'        ** \ tells you which
  !**      col 7: says 'component'                      **_/ CI coefs are big
  !**  at the end of each CI coefs (for each state)     **
  !**  are the values 'major', 'minor', and 'total'     **
  !**  which tells you how much the large and small     **
  !** components contribute to the total CI wave funct  **

  if (mod(istep,istat) .eq. 0) then
     
     do i = 1, npi                !write out orbital energies
        write(orb_io,'(i5,f16.8)') i, ener(i)       !to fort.99 in Hartrees
     enddo
     
     do i = 1,npi
        write(mo_io,1012)i,pi_iden(i),(pos(j,pi_iden(i)),j=1,3),(mocoef(i,j),j=1,npi)
     end do
     
     do i=1,npi
        write(bo_io,*) istep,i,bomat(i,i),z_atom(pi_iden(i))
        !do j=1,npi
        !   write(bo_io,1017) istep,i,j,bomat(i,j)
        !end do
     end do
     
     if (log_cicalc) then
        do j = 1,nci
           !write(cicoef_io,1016)istep,j
           sumb = 0.d0
           sumn = 0.d0
           do i = 1,nci
              write(cicoef_io,1013)istep,j,i,cib(i),cit(i),cieig(i,j)
             !    sumb = sumb + cieig(i,j)**2
             ! else
             !    write(cicoef_io,1014),j,i,cib(i),cit(i),cieig(i,j)
             !    sumn = sumn + cieig(i,j)**2
             ! end if
           end do
           !write(cicoef_io,1015)sumb,sumn,sumb+sumn
        end do
     end if
  end if
     
     
1013          format(5(i4,1x),f20.15)
1017 format(3(i4,1x),f15.10)
1012 format(i5,1x,i5,500(1x,f12.5))
!1013 format(4(i5,1x),f12.5,'  large component')
1014 format(4(i5,1x),f12.5,'  small component')
1015 format('major : ',f8.4,' minor : ',f8.4,' total : ',f8.4)
1016 format('istep = ',i6,1x,' CI coefs for state : ',i2)
  
END SUBROUTINE ppp_ppv


! AY added
!=====================================================================
! Add the terms including torsion angle to Gamma :
! This definition of the terms is given in the original PI/QCFF paper:
! (Eq.(31) in Warshel and Karplus, JACS, 94, 5612 (1972))
!
!---------------------------------------------------------------------
!  cos(phi) terms are added only to junction parts
!  (for Gamma_uu)
!  add  + Gs[ exp(-2mu(r_{u,u+1}-r^eq))*cos(phi_{u,u+1})^2
!            +exp(-2mu(r_{u,u-1}-r^eq))*cos(phi_{u,u-1})^2  ]
!  (for Gamma_uv, v=u-1,u+1)
!  add  - Gs*exp(-2mu(r_uv-r^eq))*cos(phi_uv)^2
!
!=====================================================================

SUBROUTINE  add_cos_term_to_gamma(gamma_nuc,gamma_pi,gamma_nuc_cos)

  USE commondata, ONLY: DP,npi,nnuc
  USE quantum
  USE PPV, ONLY: ntor_junc,junction,torsion
  USE classical, ONLY: atype,species

  IMPLICIT NONE
  INTEGER :: i,j,i1,ii1,i2,ii2,ijunc
  REAL(DP), DIMENSION(npi,npi)   :: gamma_pi
  REAL(DP), DIMENSION(nnuc,nnuc) :: gamma_nuc,gamma_nuc_cos
  REAL(DP) :: gs, mu_b, r0, r12, tmp
  CHARACTER(1) :: spec1, spec2


  do ijunc = 1,ntor_junc         

     !(atom number in all)
     i1  = junction(1,ijunc)
     i2  = junction(2,ijunc)
     !(atom number in pi atoms)
     ii1 = inv_pi(i1)
     ii2 = inv_pi(i2)

     do i= 1,ibeta_pairs
        if( (i1==ibeta(1,i) .and. i2==ibeta(2,i)) .or. &
            (i2==ibeta(1,i) .and. i1==ibeta(2,i)) ) then
            j=i
            exit
        endif
     enddo

     spec1 = species(atype(i1))
     spec2 = species(atype(i2))
     call  get_mu_b_r0(spec1,spec2,mu_b,r0)
     call  get_gs(spec1,spec2,gs)

     tmp     =  gs * st_exp2mu(j,1) * st_exp2mu(j,2)*st_exp2mu(j,2) 

     ! add to Gamma_uu 
     gamma_pi(ii1,ii1)  = gamma_pi(ii1,ii1)  + tmp
     gamma_pi(ii2,ii2)  = gamma_pi(ii2,ii2)  + tmp
     gamma_nuc(i1,i1)     = gamma_nuc(i1,i1)     + tmp
     gamma_nuc(i2,i2)     = gamma_nuc(i2,i2 )    + tmp
     gamma_nuc_cos(i1,i1) = gamma_nuc_cos(i1,i1) + tmp
     gamma_nuc_cos(i2,i2) = gamma_nuc_cos(i2,i2) + tmp


     ! add to Gamma_uv (v=u-1,u+1) 
     gamma_pi(ii1,ii2)  = gamma_pi(ii1,ii2)  - tmp
     gamma_pi(ii2,ii1)  = gamma_pi(ii2,ii1)  - tmp
     gamma_nuc(i1,i2)     = gamma_nuc(i1,i2)     - tmp
     gamma_nuc(i2,i1)     = gamma_nuc(i2,i1)     - tmp
     gamma_nuc_cos(i1,i2) = gamma_nuc_cos(i1,i2) - tmp        
     gamma_nuc_cos(i2,i1) = gamma_nuc_cos(i2,i1) - tmp 

  end do


  
  END SUBROUTINE add_cos_term_to_gamma

