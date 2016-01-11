!c------------------------------------------------------------
!c                                                           c
!c  subroutine computing first step in velocity verlet alg.  c
!c                                                           c
!c  1) r(t+dt)=r(t)+v(t)*dt+0.5*(f(t)/m)*dt^2                c
!c  2) v(t+dt/2) = v(t) + 0.5*(f(t)/m)*dt^2                  c
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
!* []  DECK VELVERLEI1 -- VELOCITY VERLET FOR BETAINE-30    
!* =====================================================================

SUBROUTINE verlet_s1(eki)
!SUBROUTINE VELVERLEI1_PPV(istep,nc_ppv,nap_ppv)!,vx_pi,vy_pi,vz_pi)

  USE commondata, ONLY: DP,Tmb,ttol,nnuc,log_confine,log_scalevel,iscale,&
       boltz,istep,deltat,navg
  USE qcfio, ONLY: gen_io,temp_io
  USE classical, ONLY: mass_tot,vel,vel_l,pos,pos_l,pos_m,amass,fra
  
  IMPLICIT NONE
  
  !input variables
  REAL(DP), INTENT(OUT) :: eki

  !local variables
  INTEGER :: i,j,k,iter
  REAL(DP), DIMENSION(3) :: com,ar,vtot,sub,vr,totr,ftemp
  REAL(DP) :: efactor,tfactor,eki_COM,vscale,vsq,tempi,eke,fact
  REAL :: ftemp2,gausvar

  !*  +---------------------------------------------------------------+
  !*  |  find kinetic energy of betaine-30 and save old positions     |
  !*  +---------------------------------------------------------------+

  eki=0.0
  !APW units, conversion factor for energy below 
  efactor=.01 !to kJ/mol
  !APW units, conversion factor for temperature below
  tfactor=2.0d0*120.27/(3.d0*dfloat(nnuc)-3.d0) !to Kelvin
  !APW units
  vscale=2.454058E+11 !(sqrt(J/amu)/(Angstrom/picosecond))
!!$na_loop=nat_sol
  
  do  i=1,nnuc

     vel_l(:,i) = vel(:,i)
     
     eki=eki+(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)*amass(i)
     !c        else
     !c           eki=eki+(vxa(i)*vxa(i)+vya(i)*vya(i)+vza(i)*vza(i))*winv(i)
     !c        endif
  end do
  
  eki=0.5*eki*efactor
  eke=eki
  tempi=eke*tfactor
  !eki=(eki)*few
  write(gen_io,*) 'solute temperature', eki*tfactor
  write(temp_io,'(i6,1x,f16.8)')istep,eki*tfactor
  
!8585 FORMAT(i6,1x,f16.8)


  !************************************************************************
  !************************************************************************
  !C ------------calculate COM velocity and KE----------------

  vr=0.0d0

  do i=1,nnuc
     vr(:)=vr(:)+vel(:,i)*amass(i)
  end do

  !*	 wfaciA = 1 / molecular mass
  
  vr(:)=vr(:)/mass_tot

  vsq = vr(1)**2+vr(2)**2+vr(3)**2
  eki_COM = 0.5*vsq*efactor*mass_tot
  
  write(gen_io,*) 'COM KE:',eki_COM,eki,eki-eki_COM,eki_COM*tfactor
  !************************************************************************
  !************************************************************************

  !*  +---------------------------------------------------------------+
  !*  |  set up primary and secondary arrays for adma forces          |
  !*  +---------------------------------------------------------------+

  !c-------------------------------------------------------------------
  !c For Meh-PPV Not use constraints. Application of v.v. alg. with
  !c internal force comulated in fxa,fya,fza
  !c-------------------------------------------------------------------

  !c        if(dump_vel) then
  !c           do i=1,nap_ppv
  !c              ii=nc_ppv(i)
  !c              vx_pi(i)=vxa(ii)
  !c              vy_pi(i)=vya(ii)
  !c              vz_pi(i)=vza(ii)
  !c           enddo
  !c        endif
  
  !c-------------------------------------------------------------------
  !c Here Scale Option for velocities
  !c--------------------------------------------------------------------
  
  if(log_scalevel.and.(mod(istep,iscale).eq.0)) then
     !APW already done above
     !eki=0.d0
     !do i=1,nat_sol
     !   eki=eki+(vxa(i)*vxa(i)+vya(i)*vya(i)+vza(i)*vza(i))*watom(i)
     !enddo
     !eki=(eki)*few
     !tempi=eki*ftemi3
     !fact=abs(tempi-TMPMWL)
     fact = abs(eki*tfactor-Tmb)
     if(fact.gt.ttol) then
        iter=0
        gausvar=dsqrt(boltz*Tmb)
        do while((tempi-tmb)**2.gt.9.0)
           iter=iter+1
           totr(:)=0
           eke=0.0d0
           vtot=0.0d0
           do i=1,nnuc
              ftemp(1)=dsqrt(amass(i))
              do k=1,3
                 call gasdev_s(ftemp2)
                 vel(k,i) =  vscale*ftemp2*gausvar/ftemp(1)
                 vtot(k)=vtot(k)+vel(k,i)*amass(i)
              end do
           end do
           
           sub=vtot/mass_tot
           
           do i=1,nnuc
              vel(:,i)=vel(:,i)-sub(:)
              eke=eke+amass(i)*(vel(1,i)**2+vel(2,i)**2+vel(3,i)**2)
           end do
           eke=0.5*eke*efactor
           tempi=eke*tfactor
           
           if (iter .gt. 5000) then
              !c write(*,*) ' sample velocities solute ' 
              write(gen_io,*) "***WARNING: velocity problem ***"
              !go to 1
              stop
           end if
           
           write(gen_io,*) '---> SCALING VELOCITIES FOR SOLUTE iter=',iter
           write(gen_io,*) '---> NEW TEMPERATURE T=',tempi
           
           !APW done above
           !eke=0.
           !do i=1,nnuc
           !   eke = eke+(vxa(i)*vxa(i)+vya(i)*vya(i)+vza(i)*vza(i))
           !+          *watom(i)
           ! enddo
        end do
     end if
     
     write(gen_io,*) '---> OLD AND NEW KE',eki,eke
     eki=eke
  end if
  
!c--- End Rescaling Procedure ----------------------------------
  
  !c FIRST TWO PARTS OF VELOCITY VERLET ALGORITHM:
  !c   -> calculate x(t+dt)   (new positions)
  !c   -> calculate v(t+dt/2) (half velocity)
  !c   -> requires x(t), v(t), f(t)
  do  i=1,nnuc
     
     !APW HERE pos_m and pos_l need to be added in
     ftemp(:)=pos_m(:,i)
     pos_m(:,i)=pos(:,i)
     pos(:,i)=pos_m(:,i)+vel(:,i)*deltat+0.5d0*100.0d0*deltat**2*fra(:,i)/amass(i)
     pos_l(:,i)=ftemp(:)

     vel(:,i)=vel(:,i) + 0.5d0*100.0d0*deltat*fra(:,i)/amass(i)

  end do

  !c----- Coordinates shifted to the center of the cell -----
  !c----- Periodic boundary conditions -----
  
  com=0.0
  do  i=1,nnuc
     com(:)=com(:)+pos(:,i)*amass(i)
  end do
  
  !APW units
  com(:) = com(:)/(150.*mass_tot)
  !xcmg=xxppv*wfaciA
  !ycmg=yyppv*wfaciA
  !zcmg=zzppv*wfaciA

  ar=2.0*dint(com)
  
  !ay=2.0*dint(ycmg)
  !az=2.0*dint(zcmg)
  
  !*     note the +2 here
  
  if(.not.log_confine) then
     do  i=1,nnuc
        pos(:,i)=pos(:,i)-ar(:)*150.
        pos_m(:,i)=pos_m(:,i)-ar(:)*150.
        pos_l(:,i)=pos_l(:,i)-ar(:)*150.
     end do
  end if
  
END SUBROUTINE verlet_s1
