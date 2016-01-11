!c------------------------------------------------------------
!c                                                           c
!c  subroutine check the occurring of transition between     c
!c             electronic states                             c
!c  TULLY FSSH PROCEDURE                                     c
!c                                                           c
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
!* []  DECK TRANSITION --   
!* =====================================================================

SUBROUTINE TRANSITION
  !SUBROUTINE TRANSITION_PPV

  USE commondata, ONLY: DP,nnuc,log_exdyn,eppp,istep
  USE classical, ONLY: vel,amass,fra
  USE quantum, ONLY: ifinish,iexstate,iexstate_l,ciener,dar,deltae,fexa
  USE qcfio, ONLY: gen_io,hop_io

  IMPLICIT NONE

  !local variables
  INTEGER :: iatom,k
  REAL(DP) :: efactor,aji,bji,gammaji,evib,enerex,enerex_l

  !*  +---------------------------------------------------------------+
  !*  |  converts energy in au (i.e. hartress) to program units       |
  !*  +---------------------------------------------------------------+
  !enerfac = 4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !enerfac = enerfac/fkjmol
  efactor = 2625.5

  !*  +---------------------------------------------------------------+
  !*  |  check to see if we are making a transition                   |
  !*  +---------------------------------------------------------------+
  aji=0.0d0
  bji=0.0d0
  gammaji=0.0d0


  if ((log_exdyn).and.(ifinish)) then
     !*  +---------------------------------------------------------------+
     !*  |  now get the gamma                                            |
     !*  +---------------------------------------------------------------+

     !c-------------------------------------------------------------------
     !c Added frame for Meh-PPV
     !c Note that in the case of Solute (PPV) the dax-day-daz are not 
     !c modified by the constrains procedure. IN the case of solvent are
     !c calling fixsolvent in rgkt_ppv
     !c-------------------------------------------------------------------

     do iatom = 1,nnuc
        do k=1,3
           aji = aji + dar(k,iatom)**2/(amass(iatom)*2.0d0)
           bji = bji + vel(k,iatom)*dar(k,iatom)
        end do
     end do

     !*  +---------------------------------------------------------------+
     !*  |  now for the solvent                                          |
     !*  +---------------------------------------------------------------+
     !c        if(solv) then
     !c           do i = 1,molcs
     !c              aji = aji + dsx(i)*dsx(i)/(2.0d0*wac)
     !c              aji = aji + dsy(i)*dsy(i)/(2.0d0*wac)
     !c              aji = aji + dsz(i)*dsz(i)/(2.0d0*wac)
     !c              
     !c              aji = aji + dsx(i+molcs)*dsx(i+molcs)/(2.0d0*wame)
     !c              aji = aji + dsy(i+molcs)*dsy(i+molcs)/(2.0d0*wame)
     !c              aji = aji + dsz(i+molcs)*dsz(i+molcs)/(2.0d0*wame)
     !c              
     !c              aji = aji + dsx(i+mol2)*dsx(i+mol2)/(2.0d0*waN)
     !c              aji = aji + dsy(i+mol2)*dsy(i+mol2)/(2.0d0*waN)
     !c              aji = aji + dsz(i+mol2)*dsz(i+mol2)/(2.0d0*waN)
     !c              
     !c              bji = bji + vx(i,np)*dsx(i)
     !c              bji = bji + vy(i,np)*dsy(i)
     !c              bji = bji + vz(i,np)*dsz(i)
     !c              
     !c              bji = bji + vx(i+molcs,np)*dsx(i+molcs)
     !c              bji = bji + vy(i+molcs,np)*dsy(i+molcs)
     !c              bji = bji + vz(i+molcs,np)*dsz(i+molcs)
     !c              
     !c              bji = bji + vx(i+mol2,np)*dsx(i+mol2)
     !c              bji = bji + vy(i+mol2,np)*dsy(i+mol2)
     !c              bji = bji + vz(i+mol2,np)*dsz(i+mol2)
     !c           enddo
     !c        endif

     !*  +---------------------------------------------------------------+
     !*  |  notice that we use -deltaE = final state - initial state     |
     !*  |  if we don't have sufficient energy to make the transition    |  
     !*  |  then switch back and only add to the velocities  bji/aji     |  
     !*  +---------------------------------------------------------------+
     write(hop_io,*) istep,iexstate_l,iexstate
     write(hop_io,*) 'bji=',bji
     write(hop_io,*) 'aji=',aji
     write(hop_io,*) 'deltaE=',-deltae*efactor
     if ((bji**2-4.0d0*aji*(-deltae*efactor*100.0d0)).le.0.0d0) then
        write(gen_io,*) ' '
        write(gen_io,*) 'unable to make transition'
        write(gen_io,*) ' '
        write(gen_io,*) ' bji,aji,-deltaE'
        write(gen_io,*) bji,aji,-deltae
        write(gen_io,*) '-deltaE is ',-deltae*27.21138,' eV'
        write(gen_io,*) ' '
        write(gen_io,*) 'rejected hop to state ',iexstate
        write(gen_io,*) 'remaining on state ',iexstate_l
        write(gen_io,*) 'NOT scaling velocities in any direction...'
        write(gen_io,*) ' '
        
        write(hop_io,*) '!rejecting hop!'
        write(hop_io,*) ''
        
        !c           pause
        !c           gammaji = bji/aji
        gammaji = 0.d0      !see comment above
        iexstate = iexstate_l
        ifinish = .false.

     else
        
        !*  +---------------------------------------------------------------+
        !*  |  we accept the move since we have enough energy to do it      |
        !*  |  so adjust the forces and energies for the change of state    |  
        !*  +---------------------------------------------------------------+
        !* use the CI eigenvalue to see if this improves energy conservation.

        if (iexstate_l .eq. 0) then
           enerex_l = 0.0d0
        else
           enerex_l = ciener(iexstate_l)
        end if
        
        if (iexstate .eq. 0) then
           enerex = 0.0d0
        else 
           enerex = ciener(iexstate)
        end if
        
        evib = 0.0d0*efactor
        
        !*  +---------------------------------------------------------------+
        !*  |  the forces for f(atom,0,0) should all be zero                |
        !*  +---------------------------------------------------------------+
        do iatom = 1,nnuc 
           
           fra(:,iatom) = fra(:,iatom) + fexa(:,iatom,iexstate,iexstate) - fexa(:,iatom,iexstate_l,iexstate_l)
           
        end do
        
        !c           if(solv) then
        !c              do iatom=1,matoms
        !c                 fx(iatom,1) = fx(iatom,1) 
        !c    x                 + fexx(iatom,iexstat,iexstat)
        !c    x                 - fexx(iatom,iexstatold,iexstatold) 
        !c                 fy(iatom,1) = fy(iatom,1) 
        !c    x                 + fexy(iatom,iexstat,iexstat)
        !c    x                 - fexy(iatom,iexstatold,iexstatold) 
        !c                 fz(iatom,1) = fz(iatom,1) 
        !c    x                 + fexz(iatom,iexstat,iexstat)
        !c    x                 - fexz(iatom,iexstatold,iexstatold) 
        !c              enddo
        !c           endif
        
        !*  +---------------------------------------------------------------+
        !*  |  adjust for 1 quantum of vibration for the transition         |
        !*  |  this will only work for this one time step -- you will get   |  
        !*  |  a lack of energy conservation in subsequent step --- the     |  
        !*  |  jump will be off by one quantum of vibration                 |  
        !*  +---------------------------------------------------------------+
        eppp = eppp + enerex*efactor - enerex_l*efactor + evib
        
        if (bji.gt.0.0d0) then
           !APW factor of 100 for proper units (gamma has units of ps^-1)
           gammaji = bji - dsqrt(bji**2-4.0d0*aji*(-deltae*efactor*100.0d0))
           gammaji = gammaji / (2.0d0*aji)
           
        else
           
           gammaji = bji + dsqrt(bji**2-4.0d0*aji*(-deltae*efactor*100.0d0))
           gammaji = gammaji / (2.0d0*aji)
           
        end if
        
        write(gen_io,*) 'velocity scale', istep,gammaji
        do iatom = 1,nnuc
           vel(:,iatom) = vel(:,iatom) - gammaji*dar(:,iatom)/amass(iatom)
        end do
        write(hop_io,*) 'velocity scale', istep,gammaji
        
     end if
     
     write(gen_io,*) 'bji,aji,-deltaE,gammaji'
     write(gen_io,*) bji,aji,-deltae,gammaji

  !*  +---------------------------------------------------------------+
  !*  |  add in the weighted nonadiabatic vector components to the    |
  !*  |  velocities                                                   |  
  !*  +---------------------------------------------------------------+

  !*  +---------------------------------------------------------------+
  !*  |  do the solvent                                               |
  !*  +---------------------------------------------------------------+

  !c        if(solv) then
  !c           do i = 1, molcs
  !c              vx(i,np) = vx(i,np) - gammaji*dsx(i)/wac
  !c              vy(i,np) = vy(i,np) - gammaji*dsy(i)/wac
  !c              vz(i,np) = vz(i,np) - gammaji*dsz(i)/wac
  !c              vx(i+molcs,np) = vx(i+molcs,np) - gammaji*dsx(i+molcs)
  !c    $              /wame
  !c              vy(i+molcs,np) = vy(i+molcs,np) - gammaji*dsy(i+molcs)
  !c    $              /wame
  !c              vz(i+molcs,np) = vz(i+molcs,np) - gammaji*dsz(i+molcs)
  !c    $              /wame
  !c              vx(i+mol2,np) = vx(i+mol2,np) - gammaji*dsx(i+mol2)/waN
  !c              vy(i+mol2,np) = vy(i+mol2,np) - gammaji*dsy(i+mol2)/waN
  !c              vz(i+mol2,np) = vz(i+mol2,np) - gammaji*dsz(i+mol2)/waN
  !c           enddo
  !c        endif

  !*  +---------------------------------------------------------------+
  !*  |  do the solute                                                |
  !*  +---------------------------------------------------------------+
  
end if

END SUBROUTINE TRANSITION
