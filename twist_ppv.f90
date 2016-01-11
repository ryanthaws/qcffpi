!c------------------------------------------------------------
!c  subroutine computing the dihedral contribution to forces c
!c                                                           c
!c  B*cos(teta)-> B*d(cos(teta))dx                           c
!c-----------------------------------------------------------c
!c                                                           c
!c   Modified from original subroutine of JL                 c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm/utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c


!* =====================================================================
!* []  DECK TWIST --                                                  
!* =====================================================================

!SUBROUTINE TWIST_PPV(nap_ppv,nc_ppv,pibond,ene_gap,ene_ground  
!!$     ,osc_ppv,osc_vec,istep)
SUBROUTINE twist_ppv

  USE commondata, ONLY:DP,a0,log_ppv,log_addcos
  USE classical, ONLY: length,pos,fra,atype,species
  USE qcfio, ONLY: gen_io
  USE quantum, ONLY: inv_pi,ibeta_pairs,ibeta,umubcc,linb2cc,&
       a0cc,rcc,bomat,st_betalin,inv_ibeta,st_exp2mu,z_atom
  USE PPV

  IMPLICIT NONE

  !local variables
  INTEGER :: ijunc,i,i1,i2,ii1,ii2,j
  INTEGER, DIMENSION(6) :: ind 
  REAL(DP), DIMENSION(3,6,ntor_junc):: for_temp
  REAL(DP), DIMENSION(ntor_junc) :: s
  REAL(DP), DIMENSION(2) :: ftemp
  REAL(DP), DIMENSION(3) :: dr
  REAL(DP) :: fcos,factor,forcefac
  REAL(DP) :: b_i1, b_i2, mu_b, r0, gs, tmp, tmp11,tmp22,tmp12
  CHARACTER(1) :: spec1, spec2

  !c--------------------------------------------------------------------
  !c This because we consider the beta function modulated by torsion angle
  !c--------------------------------------------------------------------
  !c  
  !c In MEH-PPV around a given junction bond we can have more than 
  !c one out of plane died. So we use Warshel definition:
  !c TETA = sum_i \phi_i
  !c 
  !c    B
  !c     \                 
  !c      B(i1)-i2         B
  !c     /        \       /
  !c    B          i3-B(i4)  
  !c                      \
  !c                       B 
  !c
  !c---------------------------------------------------------------------
  !c
  !c Loop ir run over all the bonds involved in torsion junction
  !c     - ntor gives how many torsional are involved for a specific 
  !c       junction bond
  !c     - atm_dih for each bond, and each diedral, gives atoms index in 
  !c       all sequences order
  !c     - we consider teta as average of n dihedral of each bond
  !c---------------------------------------------------------------------
  write(gen_io,*) '-------------------------------------------------'
  write(gen_io,*) 'WE START PPP COMPUTATION                         '
  write(gen_io,*) '-------------------------------------------------'
  
  !APW fix this
  forcefac=4962.34 !from Hartree/bohr radius to kJ/(mol*angstrom)
  
  !exd=.false.

  !AY calculate cosine of dihedral angle
  for_temp=0.0d0
  if(log_ppv) then
     do ijunc=1,ntor_junc
        ind=junction(:,ijunc)
        !i1=atm_d(ir,1)
        !i2=atm_d(ir,2)
        !i3=atm_d(ir,3)
        !i4=atm_d(ir,4)
        !ib1=junction(ir,1)
        !ib2=junction(ir,2)

        !AY (calculate cosine of dihedral angle and its gradient)
        call died_ppv(ind,for_temp(:,:,ijunc),fcos)
        
        if(fcos.lt.0) s(ijunc)=-1.d0               
        if(fcos.gt.0) s(ijunc)=1.d0
        
        torsion(ijunc)=abs(fcos)
     end do
  end if
  
  !AY compute the electronic matrix and energy and force of PPP theory
  call ppp_ppv
  
   
  !AY  
  !=======================================================
  !  Forces involving the gradient of cosine are added
  !=======================================================

  !AY  the forces including dcos(phi)/dr in alpha and gamma
  if(log_ppv .and. log_addcos) then
     do ijunc=1,ntor_junc         
        
        i1  = junction(1,ijunc)
        i2  = junction(2,ijunc)
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
        call  get_beta_dash(spec1,b_i1)
        call  get_beta_dash(spec2,b_i2)
        call  get_mu_b_r0(spec1,spec2,mu_b,r0)
        call  get_gs(spec1,spec2,gs)

        ind(:)  = junction(:,ijunc)
        tmp     = -2.0d0* s(ijunc)* st_exp2mu(j,1)* st_exp2mu(j,2)* forcefac

        !==== force from the terms in Alpha ===
        tmp11   = b_i1 * bomat(ii1,ii1)
        tmp22   = b_i2 * bomat(ii2,ii2)
        factor  = (tmp11 + tmp22 ) * tmp

        !==== force from the terms in Gamma ===
        !     in gradient from SCF energy and nuclear core-core replusion ===
        !(contribution from terms including gamma(u,u))
        tmp11   = 0.25d0 * gs * bomat(ii1,ii1)*bomat(ii1,ii1)
        tmp22   = 0.25d0 * gs * bomat(ii2,ii2)*bomat(ii2,ii2)
        factor  = factor + (tmp11 + tmp22) * tmp

        !(contribution from terms including gamma(u,v))
        tmp11   = - bomat(ii1,ii1) * z_atom(i2)                &
                  + 0.5d0 * bomat(ii2,ii2) * bomat(ii1,ii1)
        tmp22   = - bomat(ii2,ii2) * z_atom(i1)                &
                  + 0.5d0 * bomat(ii1,ii1) * bomat(ii2,ii2)
        tmp12   = - 0.5d0 * bomat(ii1,ii2) * bomat(ii1,ii2)    &
                  + z_atom(i1) * z_atom(i2)
        tmp11   = - gs * tmp11
        tmp22   = - gs * tmp22
        tmp12   = - gs * tmp12
        factor  = factor + (tmp11 + tmp22 + tmp12) * tmp


        !=== Add the forces ===
        do i=1,6
           fra(:,ind(i)) = fra(:,ind(i)) + factor * for_temp(:,i,ijunc)
        end do

        
     end do

  endif

  !AY  the forces including dcos(phi)/dr in beta
  if(log_ppv) then
     do ijunc=1,ntor_junc         
        
        !factor=-s(ir)*facppp(ir)
        !APW from ppp_ppv.f90, facppp(ir)=2.0*bomat(ii1,ii2)*f
        
        i1  = junction(1,ijunc)
        i2  = junction(2,ijunc)
        ii1 = inv_pi(i1)
        ii2 = inv_pi(i2)
        
        !dr  = (pos(:,i1)-pos(:,i2))/a0
        !ftemp(1) = dsqrt(length(dr))
        !APW changed
        !ftemp(2)=betalin(ftemp(1),species(atype(i1)),species(atype(i2)))
        ftemp(2) = st_betalin(inv_ibeta(i1,i2),1)
        
        factor = -s(ijunc) * 2.0 * bomat(ii1,ii2) * ftemp(2) * forcefac
        
        ind = junction(:,ijunc)
!!$     i1=atm_d(ir,1)
!!$     i2=atm_d(ir,2)
!!$     i3=atm_d(ir,3)
!!$     i4=atm_d(ir,4)
!!$
!!$     ib1=junction(ir,1)
!!$     ib2=junction(ir,2)

        do i=1,6
           fra(:,ind(i)) = fra(:,ind(i)) + factor * for_temp(:,i,ijunc)
        end do
        
     end do
  end if



  ! Compute Bond Distance between PI atoms
  do i=1,ibeta_pairs     

     i1=ibeta(1,i)
     i2=ibeta(2,i)

     dr=(pos(:,i2)-pos(:,i1))
     ftemp(1)=dsqrt(length(dr))
     !dr=(xna(i2)-xna(i1))**2
     !dr=dr+(yna(i2)-yna(i1))**2
     !dr=dr+(zna(i2)-zna(i1))**2
     
     !APW length units again
     pibond(i)=ftemp(1)
     !pibond(i)=sqrt(dr)*unitl

  end do
  
END SUBROUTINE twist_ppv



