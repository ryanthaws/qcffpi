!* =====================================================================
!* []  DECK SETUPPPP -- ZEROTH ORDER APPROXIMATION TO FOCK MATRIX
!*                      ELEMENTS
!* =====================================================================

SUBROUTINE setup_ppp(fock,hmat,beta,gamma_nuc)
  !SUBROUTINE SETUPPPP_PPV(FOCK,BETA,HMAT,BO,nap,ibeta,nc_ppv,navg,
  !x mos,iter,istep)

  USE commondata, ONLY:npi,dp,nnuc,a0,log_delre,log_ppv,log_addcos
  USE classical, ONLY: pi_iden,length,pos,dfield,species,atype
  USE ppv, ONLY: ntor_junc,torsion
  USE quantum
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  !input variables
  REAL(DP), DIMENSION(nnuc,nnuc), INTENT(IN) :: gamma_nuc
  REAL(DP), DIMENSION(npi,npi), INTENT(OUT) :: fock,hmat,beta
  REAL(DP), DIMENSION(2) :: ftemp
  REAL(DP), DIMENSION(3) :: dr
  !local variables
  INTEGER :: i,j,k,ii,jj,kk,idh,i1,i2, ii1, ii2
  REAL(DP) :: angle,scdfield,comm,dipx,e_qm_dip2,e_cl_dip2, tmp
  REAL(DP) :: b_i1, b_i2
  REAL(DP), DIMENSION(npi) :: alpha_cos_term


  !*  +---------------------------------------------------------------+
  !*  |  set all elements to zero then fill in the nonzero elements   |  
  !*  +---------------------------------------------------------------+

  fock=0.d0
  beta=0.d0
  hmat=0.d0
           
  !*  +---------------------------------------------------------------+
  !*  |  do diagonal elements first of the fock matrix and            |  
  !*  |  1-electron matrix                                            |  
  !*  |                                                               |  
  !*  |   [Simplified model (log_addcos=.F.)]                         |  
  !*  |      F_uu = H_uu + 0.5*P_uu*Y_uu + sum(k .ne. u) P_kk *Y_uk   |
  !*  |                                                               | 
  !*  |      H_uu = alpha_u - sum(k .ne. u) Z_k*Y_uk                  |  
  !*  |                                                               | 
  !*  |      P_uv = 2 * sum(i) c_u(i) c_v(i)                          |  
  !*  |                                                               |  
  !*  |      Y_uv = e^2 / (R_uv + a_uv)     Mataga-Nishimoto          |  
  !*  |      a_uv = 2*e^2 / (Y_uu + Y_vv)   calcd in gammaset_ppv     |  
  !*  |                                                               |
  !*  | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
  !*  |   [Model fully depending on torsion angle (log_addcos=.T.)]   |  
  !*  |      F_uu = H_uu + 0.5*P_uu*Y_uu + sum(k .ne. u) P_kk *Y_uk   |
  !*  |                                                               | 
  !*  |      H_uu = alpha_u - sum(k .ne. u) Z_k*Y_uk                  |  
  !*  |         + beta'[exp(-2myu_B(R_{u,u+1}-R0))*cos^2(phi_{u,u+1}) |
  !*  |               + exp(-2myu_B(R_{u,u-1}-R0))*cos^2(phi_{u,u-1})]|
  !*  |                                                               | 
  !*  |      P_uv = 2 * sum(i) c_u(i) c_v(i)                          |  
  !*  |                                                               |  
  !*  |      Y_uu = (I-A)                                             |
  !*  |           + Gs*[exp(-2myu_B(R_{u,u+1}-R0))*cos^2(phi_{u,u+1}) |
  !*  |               + exp(-2myu_B(R_{u,u-1}-R0))*cos^2(phi_{u,u-1})]|
  !*  |      (for v=u-1,u+1)                                          |
  !*  |      Y_uv = G'*exp(-myu_Y R_uv) + e^2 / (R_uv + D)            |
  !*  |           - Gs*exp(-2myu_B(R_{u,u+1}-R0))*cos^2(phi_{u,u+1})  |
  !*  |      (for v .ne. u,u-1,u+1)                                   |
  !*  |      Y_uv = G'*exp(-myu_Y R_uv) + e^2 / (R_uv + D)            |
  !*  |           G' = (I-A)-G0                                       |
  !*  |           D  = e^2 / G0                                       |
  !*  +---------------------------------------------------------------+
  
  !c-------------------------------------------------------------------
  !c In PPV, PI atom are just carbon...
  !c alpha_u = UC (defined in paramPPP.h)
  !c vs(ii) = solvent contribution to hamiltonian
  !c gamma_nuc(ii,jj) = two-electron repulsion integral (gammaset_ppv.f)
  !c facion(jj) = nuclear charge (Z_p = 1) equal to # of 2p electrons
  !c              that atom jj donates to the pi-system (diag.
  !c              elements of hamiltonian only) ... set in set_ppv.f
  !c bo(i,j) = bond-order matrix of expansion coefficients
  !c nap = # of pi atoms in systme
  !c nc_ppv(i) tells you the absolute index number of the ith pi-atom
  !c-------------------------------------------------------------------
  
  do i = 1,npi   ! Loop over the AOs
     ii=pi_iden(i)
     ftemp(1)=0.d0
     ftemp(2)=0.d0
     do j = 1,npi
        jj=pi_iden(j)    
        if (j.ne.i) then
           ftemp(1)=ftemp(1)+(bomat(j,j)*gamma_nuc(ii,jj))
           ftemp(2)=ftemp(2)+z_atom(jj)*gamma_nuc(ii,jj)
        end if
     end do
     
     !APW allow alpha_mu the valance state potential to be defined for each pi-atom type
     !hmat(i,i) = alpha_mu - ftemp(2)!+ vs(ii)  !SOLVENT CONTRIBUTIONS
     !APW gamma
     !hmat(i,i) = alpha_mu(atype(ii)) - ftemp(2)
     hmat(i,i) = alpha_mu_i(ii) - ftemp(2)
     fock(i,i) = hmat(i,i) + 0.5e0*bomat(i,i)*gamma_nuc(ii,ii)+ftemp(1) !+ vs(ii)
     
     !c---  This extended loop take into account del Re charge on NoPi Atoms
     !c---  I'm not taken into accoutn sigma charge now from del re
     if(log_delre) then
        do j=1,nnopi
           jj=ind_nopi(j)
           hmat(i,i) = hmat(i,i) - z_atom(jj)*gamma_nuc(ii,jj)
           fock(i,i) = fock(i,i) - z_atom(jj)*gamma_nuc(ii,jj)
        end do
     end if
  end do

  !AY added
  ! This loop take into account the terms depending on torsion angle
  if(log_addcos .and. log_ppv) then
     alpha_cos_term(:) = 0.0d0
     do i= 1,ibeta_pairs
        i1 = ibeta(1,i)
        i2 = ibeta(2,i)
        if( junc_flg(i1) .and. junc_flg(i2) ) then
           if(inv_dih(i1,i2)==0) then
              write(gen_io,*) '***Warning: error in setup_ppp.f90 ***'
              stop
           end if

           call  get_beta_dash(species(atype(i1)),b_i1)
           call  get_beta_dash(species(atype(i2)),b_i2)

           ii1 = inv_pi(i1)
           ii2 = inv_pi(i2)

           tmp = st_exp2mu(i,1) * st_exp2mu(i,2)* st_exp2mu(i,2) 
           alpha_cos_term(ii1) = alpha_cos_term(ii1) + b_i1 * tmp
           alpha_cos_term(ii2) = alpha_cos_term(ii2) + b_i2 * tmp

        end if
     enddo    ! i=1,ibeta_pairs

     do j = 1,npi
        hmat(j,j) = hmat(j,j) + alpha_cos_term(j)
        fock(j,j) = fock(j,j) + alpha_cos_term(j)
     enddo

  endif


  


  !*  +---------------------------------------------------------------+
  !*  |  now for da off diagonal elements - should I do him now boss? |
  !*  |  let me do him for ya boss? Please. He'll sleep wid da fishes |
  !*  |  good tonight  <-[italian gangster talk by FS?]               |
  !*  |                                                               |  
  !*  |      F_uv = H_uv - 0.5*P_uv*Y_uv      (u .ne.v)               |
  !*  |                                                               | 
  !*  |      H_uv = beta_uv        if u and v are covalent            |  
  !*  |      H_uv = 0              otherwise                          |  
  !*  |                                                               |  
  !*  +---------------------------------------------------------------+
  do i = 1,npi
     ii=pi_iden(i)
     do j = 1,npi
        jj=pi_iden(j)
        if (i.ne.j) fock(i,j) = -.5e0*bomat(i,j)*gamma_nuc(ii,jj)
     end do
  end do

  !*  +---------------------------------------------------------------+
  !*  |  now for off-diagonal resonance integrals!!!!                 |
  !*  |  this is legacy code here. At one point we used the variable  |  
  !*  |  beta method to the resonance terms where beta_ij = A_0       |  
  !*  |  + bo(i,j)*A_1;  This form of the PPP hamiltonian does not    |  
  !*  |  conserve energy and requires evaluation of gradients of      |  
  !*  |  the MO coefficients;  We use the lindeberg formula for       |  
  !*  |  beta_ij;  Under the Lindeberg formula beta_ijs are constant; |  
  !*  |  Hence the A1 values are all zero here                        |  
  !*  |                                                               |  
  !*  |  we use a new beta function with distance dependance here;    | 
  !*  |  this one has been found by a fit to the lindberg beta        |  
  !*  |                                                               |  
  !*  |      beta_uv = (h^2/m_e) * (1/R_uv) * (dS_uv/dR_uv)           |  
  !*  |                                                               |  
  !*  +---------------------------------------------------------------+
  
  !c-------------------------------------------------------------------
  !c For MEH-PPV system for index in betalin is 1: C-C Bond
  !c-------------------------------------------------------------------
  !c For PPV, the ring-twist is an important parameter, we use
  !c JUNCTION logical to check if atoms I and J are involved in 
  !c junction torsion. Torsion angle are sequentially oredered as bond
  !c for PI atoms. 
  !c-------------------------------------------------------------------
  !do i=1,nnuc             !recalls all atomic positions
  !   rposx(i) = xna(i)*unitl/0.529177 !from internal units to bohr
  !   rposy(i) = yna(i)*unitl/0.529177
  !   rposz(i) = zna(i)*unitl/0.529177
  !end do
  
  do i=1,ibeta_pairs     !loops over c-c bonds (pi-c only)
     
     i1=ibeta(1,i)          !tells you molecular index
     i2=ibeta(2,i)
     
     kk=inv_pi(ibeta(1,i)) !inverts from molecular index to
     jj=inv_pi(ibeta(2,i)) !the pi index; kk is the
                           !ibeta(i,1)th carbon atom.
                           !that is, nc_ppv(kk) = ibeta(i,1)
     
!!$     angle=1.0d0
!!$     if(log_ppv) then
!!$        if (junc_flg(i1).and.junc_flg(i2)) then
!!$           idh=inv_dih(i1,i2)
!!$           angle=torsion(idh)
!!$        end if
!!$     end if
!!$     
!!$     dr=(pos(:,i1)-pos(:,i2))/a0
!!$     ftemp(1)=dsqrt(length(dr))
!!$     
!!$     ftemp(2)=betalin(ftemp(1),species(atype(i1)),species(atype(i2)))
!!$     
!!$     beta(kk,jj)=ftemp(2)*angle
     beta(kk,jj)=st_betalin(i,1)*st_betalin(i,2)
  end do  !do i=1,ibetapairs_ppv
  !*******************************************************************
  !*                                                                 *
  !*  now for the presence of an external x-polarized electric field *
  !*                                                                 *
  !*  start by forming the x-component of the molecular dipole       *
  !*  moment (a one-electron operator)                               *
  !*                                                                 *
  !* Field strength, E = -1.73445x10-5 Hartree/bohr, from            *
  !*                                 Bittner, PRB 62, 11473 (2000)   *
  !*                                                                 *
  !* From Zhao, Chem Phys Lett, 351, 481 (2002), small and large     *
  !* field strengths run from 0.00024 -> 0.0058 Hartree/bohr         *
  !*                                                                 *
  !* V_e gets added directly to the fock matrix                      *
  !*                                                                 *
  !* 1 V/m   = 1.94465d-12 Hartree / (e*bohr)                        *
  !* 1 MV/cm = 1.94465d-04 Hartree/ (e*bohr)                         *
  !* 1 Hartree / (e*bohr) = 5.1423d+11 V/m                           *
  !* 1 Hartree / (e*bohr) = 5.1423d+3 MV/cm                          *
  !*                                                                 *
  !*******************************************************************
  !JJBL
  
  scdfield = dfield
  !c      scdfield = dfield * (1.d0-exp(-(dble(istep)/30.d0)**2))

  
  comm = 0.d0
  dipx = 0.d0
  e_qm_dip2 = 0.d0
  e_cl_dip2 = 0.d0
  do i = 1,npi
     ii=pi_iden(i)
     comm = comm + pos(1,ii)/a0
     dipx=dipx+pos(1,ii)*(1.d0-bomat(i,i))/a0
     e_qm_dip2 = e_qm_dip2-scdfield*pos(1,ii)*bomat(i,i)/a0
     e_cl_dip2=e_cl_dip2+scdfield*pos(1,ii)/a0
  end do
  comm = comm / dble(npi)
  !JJBL
  do i = 1,npi    !loop over the AOs
     ii = pi_iden(i)
     !c         comm = 0.d0
     fock(i,i) = fock(i,i) - scdfield*bomat(i,i) * (pos(1,ii)/a0-comm)
     !*(-1.d0)**i
     hmat(i,i) = hmat(i,i) - scdfield*bomat(i,i) * (pos(1,ii)/a0-comm)
     !* (-1.d0)**i
  end do

  !*  +---------------------------------------------------------------+
  !*  |  copy the upper triangular half of the matrix into the        |  
  !*  |  lower half                                                   |  
  !*  +---------------------------------------------------------------+
  
  do i=1,npi
     do j=1,i-1
        beta(i,j) = beta(j,i)
     end do
  end do

  do i=1,npi
     do j = 1,npi
        if (i.ne.j) then
           fock(i,j) = fock(i,j) + beta(i,j)
           hmat(i,j) = beta(i,j) 
        end if
     end do
  end do
END SUBROUTINE setup_ppp


!----------------------------------------------------------
!  Get Beta dash parameter from a species of atom
!----------------------------------------------------------
SUBROUTINE get_beta_dash(spec1,b)

  USE commondata, ONLY:dp
  USE quantum
  IMPLICIT NONE
  LOGICAL :: hit
  REAL(DP) :: b
  CHARACTER*1, INTENT(IN) :: spec1

  hit=.false.
  if(     spec1=='A') then
     b = beta_dash_c
  else if(spec1=='S') then
     b = beta_dash_s
  else if(spec1=='N') then
     !b = beta_dash_np 
     hit = .true.
  else if(spec1=='M') then
     !b = beta_dash_n 
     hit = .true.
  else if(spec1=='O') then
     b = beta_dash_o
  else 
     hit = .true.
  endif
  if(hit) then
     print *,'***WARNING: combo not defined in beta_dash ',spec1,' ***'
     stop
  end if
  
  RETURN
END SUBROUTINE get_beta_dash


!AY added
!-------------------------------------------------------------
!  Get mu_b and r^eq parameter from a pair of species of atom
!-------------------------------------------------------------
SUBROUTINE get_mu_b_r0(spec1,spec2,mu_b,r0)

  USE commondata, ONLY:dp
  USE quantum
  IMPLICIT NONE
  LOGICAL :: hit
  REAL(DP) :: mu_b, r0
  CHARACTER*1, INTENT(IN) :: spec1,spec2

  hit=.false.
  if(      spec1=='A' .and. spec2=='A') then
     mu_b = umubcc
     r0   = rcc
  else if((spec1=='A' .and. spec2=='S') .or.   &
          (spec1=='S' .and. spec2=='A')) then
     mu_b = umubcs
       r0   = rcs
  else if((spec1=='A' .and. spec2=='N') .or.   &
          (spec1=='N' .and. spec2=='A')) then
       mu_b = umubcn
       r0   = rcn
  else if((spec1=='A' .and. spec2=='M') .or.   &
          (spec1=='M' .and. spec2=='A')) then
       mu_b = umubcn
       r0   = rcn
  else if((spec1=='A' .and. spec2=='O') .or.   &
          (spec1=='O' .and. spec2=='A')) then
       mu_b = umubco
       r0   = rco
  else
     hit=.true.
  end if
    
  if(hit) then
    print *,'***WARNING: combo not defined in get_mu_b_r0 ',spec1,spec2,' ***'
    stop
  end if

RETURN
END SUBROUTINE get_mu_b_r0

!AY added
!----------------------------------------------------------
!  Get Gs parameter from a pair of species of atom
!----------------------------------------------------------
SUBROUTINE get_gs(spec1,spec2,gs)

  USE commondata, ONLY:dp,d_screen
  USE quantum
  IMPLICIT NONE
  LOGICAL :: hit
  REAL(DP) :: gs
  CHARACTER*1, INTENT(IN) :: spec1,spec2

  hit=.false.

  if(      spec1=='A' .and. spec2=='A') then
     gs   = gs_cc
  else if((spec1=='A' .and. spec2=='S') .or.   &
       (spec1=='S' .and. spec2=='A')) then
     gs=gs_cs
  else if((spec1=='A' .and. spec2=='N') .or.   &
       (spec1=='N' .and. spec2=='A')) then
     hit=.true.
     !       gs   = gs_cn ?
     !       gs   = gs_cnp ? which one ?
  else if((spec1=='A' .and. spec2=='M') .or.   &
       (spec1=='M' .and. spec2=='A')) then
     hit=.true.
     !       gs   = gs_cn ??  what is M?
  else if((spec1=='A' .and. spec2=='O') .or.   &
       (spec1=='O' .and. spec2=='A')) then
     gs   = gs_co
  else if( spec1=='O' .and. spec2=='O') then
     gs   = gs_oo
  else if( spec1=='N' .and. spec2=='N') then
       hit=.true.
       gs   = gs_nn  ! ??
    else
       hit=.true.
    end if
    
  if(hit) then
    print *,'***WARNING: combo not defined in get_gs ',spec1,spec2,' ***'
    stop
  end if
  
  gs=gs/d_screen
  RETURN
END SUBROUTINE get_gs


