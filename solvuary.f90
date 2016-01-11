!c------------------------------------------------------------
!c                                                           c
!c  subroutine computing the force in exc. state              
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
!c-----------------------------------------------------------c

!* =====================================================================
!* []  DECK SOLVUARY -- FIND FIRST DERIVATIVE OF UNITARY 
!*                      TRANSFORMATION MATRIX                    
!* =====================================================================

!*  +---------------------------------------------------------------+
!*  |  Pople, JA, et al,Int. J. of Quantum Chemistry: Quantum       |
!*  |  Chem. Symp. 13, 225 (1979)                                   |
!*  +---------------------------------------------------------------+

SUBROUTINE SOLVUARY(gamma_nuc,gamma_pi,ener)

  USE commondata, ONLY: DP,npi,a0,nnuc,norb,log_delre,log_ppv
  USE classical, ONLY: length,pos,gamma,pi_iden,atype,locate_atom,npispec,pispec_code,species
  USE quantum
  USE PPV, ONLY: junction,ntor_junc,torsion
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE

  !input variables
  REAL(DP), DIMENSION(npi), INTENT(IN) :: ener
  REAL(DP), DIMENSION(npi,npi), INTENT(IN) :: gamma_pi
  REAL(DP), DIMENSION(nnuc,nnuc), INTENT(IN) :: gamma_nuc
  
  !local variables
  LOGICAL :: ispi
  INTEGER :: i,j,k,a,i1,i2,i3,i4,iaa,ibb,ib1,ib2,ii1,ii2,ii3,ii4
  INTEGER :: id1,id2,idx1,idx2,idh,idx,info,ir,jj,p,q,u,v,rho
  INTEGER, DIMENSION(6) :: ind
  INTEGER, DIMENSION(npi*(npi-1)/2) :: ipvt
  REAL(DP), DIMENSION(npi,npi,npi) :: hx,hy,hz     ! AY
  REAL(DP), DIMENSION(npi,norb+1:npi,norb) :: xai  ! AY
  REAL(DP), DIMENSION(npi*(npi-1)/2,npi*(npi-1)/2) :: aea
  REAL(DP), DIMENSION(3) :: gamfac,th,th1,th2,th3,ftemp,for_temp,dr
  REAL(DP), DIMENSION(3,6,ntor_junc) :: ftemp2
  !REAL(DP), DIMENSION(npi*(npi-1)/2) :: work
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: work
  INTEGER :: lwork
  REAL(DP) :: angle,fcos,csuc,csuc1,csuc2,factor,gam2,rij,s,tmp1,tmp2,&
       tmp3,xaiuu,xaiuv,xairhorho
  REAL(DP) :: time_f,time_i
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: index1
  INTEGER :: ind1,n
!-------------------------------------------------------------------
  integer clock_max,clock_rate,counti,count1,count2,count3
  
  !APW for scalapack
  !INTEGER ictxt, nprow, npcol, myrow, mycol
  !External Subroutines
  !External BLACS_EXIT, BLACS_GRIDEXIT, BLACS_GRIDINFO, DESCINIT, MATINIT, pdgetrf, pdgetri, sl_init

  
  call system_clock ( counti, clock_rate, clock_max )
  call cpu_time(time_i)
  id2 = npi*(npi-1)/2
  
  lwork=id2*64
  ALLOCATE(work(lwork))
  
  aea=0.0d0
  
  !*  +---------------------------------------------------------------+
  !*  |  construct elements of the B vector:                          |
  !*  |                                                               |  
  !*  |  B'(qp) = H'(qp) - 0.5*sum(u .ne. v) C_u^q C_v^p P_uv Y'_uv   |  
  !*  |                  + 0.5*sum(u)        C_u^q C_u^p P_uu Y'_uu   |  
  !*  |                  +     sum(u .ne. k) C_u^q C_u^p P_kk Y'_uk   |  
  !*  |                                                               |  
  !*  |  H'(qp) = sum(u;v) C_u^q [ d H_uv / dR ] C_v^p                |  
  !*  +---------------------------------------------------------------+

  !c--------------------------------------------------------------------
  !c NC associate to index PI atom, real index in of atom in Meh-PPV
  !c--------------------------------------------------------------------
  !do i= 1,npi
!!$            rposx(i) = xna(nc_ppv(i))*unitl/0.529177
!!$            rposy(i) = yna(nc_ppv(i))*unitl/0.529177
!!$            rposz(i) = zna(nc_ppv(i))*unitl/0.529177
!!$
!!$CCC    MBH comment out solvent forces
!!$c           do j = 1,nap
!!$c              dvsdrx(nc_ppv(i),nc_ppv(j)) = dvsdrx(nc_ppv(i),nc_ppv(j))
!!$c    $              /forcefac
!!$c              dvsdry(nc_ppv(i),nc_ppv(j)) = dvsdry(nc_ppv(i),nc_ppv(j))
!!$c    $              /forcefac
!!$c              dvsdrz(nc_ppv(i),nc_ppv(j)) = dvsdrz(nc_ppv(i),nc_ppv(j))
!!$c    $              /forcefac
!!$c           enddo
!!$            
!!$      enddo

  !*  +---------------------------------------------------------------+
  !*  |  we don't zero derivatives of the unitary transformation      |
  !*  |  matrix elements out here because we want the zeroth order    |
  !*  |  values to be the values from the previous time steps         |
  !*  |  instead of using  A_0 matrix - we zero them out in setup.f   |
  !*  |  and the zeroth order guess will be the A_0 matrix on the     |
  !*  |  1st time step                                                |
  !*  +---------------------------------------------------------------+

  hx=0.0d0
  hy=0.0d0
  hz=0.0d0
  qr=0.0d0
  
  !APW computed but not used
  !if(log_delre) then
  !   hnopix=0.0d0
  !   hnopiy=0.0d0
  !   hnopiz=0.0d0
  !end if

  !APW computed but not used
  !hnopidihx=0.0d0
  !hnopidihy=0.0d0
  !hnopidihz=0.0d0
  
  !*  +---------------------------------------------------------------+
  !*  |  contribution from diagonal elements of the one-electron      |  
  !*  |  core matrix:                                                 |  
  !*  |                                                               |  
  !*  |  H_uu = alpha_u - sum(k .ne. u) Z_k Y_uk                      |  
  !*  |                 - sum(j) e*q_j / R_ju                         |  
  !*  +---------------------------------------------------------------+
  
  !*  +---------------------------------------------------------------+
  !*  |  solvent-electron charge interaction part:                    |
  !*  +---------------------------------------------------------------+
  

  !CCC   MBH comment out solvent contributions
  !c        do v = 1,loop
  !c           do u = 1,loop
  !c              do p = 1,loop
  !c                 do q = p,loop
  !c                    
  !c                    csuc = mos(u,q)*mos(u,p)
  !c                    Hx(q,p,v) = Hx(q,p,v) + csuc*dvsdrx(nc_ppv(u)
  !c    $                       ,nc_ppv(v))
  !c                    Hy(q,p,v) = Hy(q,p,v) + csuc*dvsdry(nc_ppv(u)
  !c    $                       ,nc_ppv(v))
  !c                    Hz(q,p,v) = Hz(q,p,v) + csuc*dvsdrz(nc_ppv(u)
  !c    $                       ,nc_ppv(v))
  !c                 enddo
  !c              enddo
  !c           enddo
  !c        endd

  !*  +---------------------------------------------------------------+
  !*  |  electron-core interaction part (depends on Y_uk)             |
  !*  +---------------------------------------------------------------+
  !c-------------------------------------------------------------------
  !c Now rposix run over just PI atom. See definition above.
  !c I use Gamma running over just PI atom.
  !c Be carefull with charge. Facion follow all atom sequence
  !c-------------------------------------------------------------------

  !c-------------------------------------------------------------------
  !c I need to add contribution. e NoPI atom. This because I have 
  !c delRe sigma charge over all atom.
  !c-------------------------------------------------------------------
  do u =1,npi
     do rho =u+1,npi

        dr = (pos(:,pi_iden(rho))-pos(:,pi_iden(u)))/a0
        rij=dsqrt(length(dr))
        
        gam2 = -  gamma_pi(u,rho)**2/rij
        gamfac = gam2*(pos(:,pi_iden(rho))-pos(:,pi_iden(u)))/a0
        
        th1 = z_atom(pi_iden(rho))*gamfac
        th2 = z_atom(pi_iden(u))*gamfac
        do p = 1,npi
           do q = p,npi
              csuc1 = mocoef(u,q)*mocoef(u,p)
              csuc2 = mocoef(rho,q)*mocoef(rho,p)
              hx(q,p,rho) = hx(q,p,rho) - csuc1*th1(1) - csuc2*th2(1) 
              hx(q,p,u) = hx(q,p,u) + csuc1*th1(1) + csuc2*th2(1) 

              hy(q,p,rho) = hy(q,p,rho) - csuc1*th1(2) - csuc2*th2(2)  
              hy(q,p,u) = hy(q,p,u) + csuc1*th1(2) + csuc2*th2(2) 
              
              hz(q,p,rho) = hz(q,p,rho) - csuc1*th1(3) - csuc2*th2(3)  
              hz(q,p,u) = hz(q,p,u) + csuc1*th1(3) + csuc2*th2(3) 
           end do
        end do
     end do
  end do
  !call system_clock ( count1, clock_rate, clock_max )
  !print *,'shit1',(count1-counti)/real(clock_rate)
  !c-------------------------------------------------------------------
  !c Here I add contribution from e and NoPI atoms
  !c
  !c NoPI atom position is defined in all atoms index. So I define
  !c rposx_noPI, rposy_noPI, rposz_noPI
  !c I use gamma_nuc that cover all atom index
  !c-------------------------------------------------------------------
  
  
  if(log_delre) then
     
     do u=1,npi
        do v=1,nnopi
           jj=ind_nopi(v)
           
           !rposx_noPI=xna(jj)*unitl/0.529177
           !rposy_noPI=yna(jj)*unitl/0.529177
           !rposz_noPI=zna(jj)*unitl/0.529177
           
           dr=(pos(:,jj)-pos(:,pi_iden(u)))/a0
           rij=dsqrt(length(dr))
           
           gam2=-gamma_nuc(pi_iden(u),jj)**2/rij
           gamfac=gam2*(pos(:,jj)-pos(:,pi_iden(u)))/a0
!!$                  gam2 = -  gamma_nuc(nc_ppv(u),jj)**2/rij
!!$                  gamfacx = gam2*(rposx_noPI-rposx(u))
!!$                  gamfacy = gam2*(rposy_noPI-rposy(u))
!!$                  gamfacz = gam2*(rposz_noPI-rposz(u))
           
           th3=z_atom(jj)*gamfac !APW vectors
           
!!$                  thx3 = facion(jj)*gamfacx
!!$                  thy3 = facion(jj)*gamfacy
!!$                  thz3 = facion(jj)*gamfacz
           
              do p=1,npi
                 do q=p,npi
                    csuc1=mocoef(u,q)*mocoef(u,p)
                    hx(q,p,u) = hx(q,p,u) + csuc1*th3(1)
                    hy(q,p,u) = hy(q,p,u) + csuc1*th3(2)
                    hz(q,p,u) = hz(q,p,u) + csuc1*th3(3)
                    !c----------------------------------------------------------
                    !c Here I consider the force over the noPI atoms
                    !c I will treat as I treat solvent atoms 
                    !c in addexstsoleps_ppv.f
                    !c---------------------------------------------------------- 
                    
                    !APW computed but not used
                    !hnopix(q,p,v)=hnopix(q,p,v) - csuc1*th3(1)
                    !hnopiy(q,p,v)=hnopiy(q,p,v) - csuc1*th3(2)
                    !hnopiz(q,p,v)=hnopiz(q,p,v) - csuc1*th3(3)
                    
                 end do
              end do

           end do
        end do
     end if

     !*  +---------------------------------------------------------------+
     !*  |  off-diagonal elements of the one-electron core matrix:       |  
     !*  |                                                               |  
     !*  |  H_uv = Beta_uv         (Lindberg approximation)              |  
     !*  +---------------------------------------------------------------+
     
     do i=1, ibeta_pairs
        !atomic index
        i1=ibeta(1,i)
        i2=ibeta(2,i)
        !pi index
        ii1=inv_pi(i1)
        ii2=inv_pi(i2)
        
!!$        angle=1.0d0
!!$        if(log_ppv) then
!!$           if(junc_flg(i1).and.junc_flg(i2)) then
!!$              idh=inv_dih(i1,i2)
!!$              if(idh.eq.0) then
!!$                 write(gen_io,*) '***Warning: error in subroutine solvuary.f90, idh=0***'
!!$                 stop
!!$              end if
!!$              angle = torsion(idh)
!!$           end if
!!$        end if
!!$           
!!$        dr=(pos(:,i1)-pos(:,i2))/a0
!!$        ftemp(1)=dsqrt(length(dr))
!!$        ftemp(2)=betalin(ftemp(1),species(atype(i1)),species(atype(i2)))
!!$        for_temp=fbetalin(dr,species(atype(i1)),species(atype(i2)),ftemp(2),angle)
        for_temp=st_fbetalin(i,:)
        
        do p=1,npi
           do q=p,npi
              csuc = mocoef(ii1,q)*mocoef(ii2,p)+mocoef(ii1,p)*mocoef(ii2,q)
              
              hx(q,p,ii1) = hx(q,p,ii1) - for_temp(1)*csuc
              hx(q,p,ii2) = hx(q,p,ii2) + for_temp(1)*csuc
              hy(q,p,ii1) = hy(q,p,ii1) - for_temp(2)*csuc
              hy(q,p,ii2) = hy(q,p,ii2) + for_temp(2)*csuc
              hz(q,p,ii1) = hz(q,p,ii1) - for_temp(3)*csuc
              hz(q,p,ii2) = hz(q,p,ii2) + for_temp(3)*csuc
              
           end do
        end do
     end do
     !call system_clock ( count2, clock_rate, clock_max )
     !print *,'shit2',(count2-count1)/real(clock_rate)
     !*  +---------------------------------------------------------------+
     !*  |  Beta_uv ==> Beta_uv * cos(theta)                             |
     !*  +---------------------------------------------------------------+
     !c------------------------------------------------------------------
     !c Some dihedral involves NoPI atoms.
     !c 1. Use of xna yna zna
     !c 2. Definition of Q_noPI_dih to be treated in solvent scheme
     !c------------------------------------------------------------------

     !exd=.false.
     
     if(log_ppv) then
        ftemp2=0.0d0
        do ir=1,ntor_junc
           i1=junction(3,ir)
           i2=junction(4,ir)
           i3=junction(5,ir)
           i4=junction(6,ir)
           
           !c--- Atom in Central Bond -----------------------------------------
           ib1=junction(1,ir)
           ib2=junction(2,ir)
           
           iaa=inv_pi(ib1)
           ibb=inv_pi(ib2)

           !c-------------------------------------------------------------------
           !c---  Compute force over all atoms involved in torsion -------------
           !c---  Here I pass logical exd and scale xna position inside died ---
           !c---  No any supplementar conversion for forces               ------
           !c-------------------------------------------------------------------
           
           ind=junction(:,ir)
           call died_ppv(ind,ftemp2(:,:,ir),fcos)
           
           if(fcos.lt.0) s=-1.d0               
           if(fcos.gt.0) s=1.d0
           
           !dr=(pos(:,ib1)-pos(:,ib2))/a0
           !ftemp(1)=dsqrt(length(dr))
           !factor = s*betalin(ftemp(1),species(atype(ib1)),species(atype(ib2)))
           !APW changed
           factor=s*st_betalin(inv_ibeta(ib1,ib2),1)

           !APW resolve this factor!
           !factor=s*factor/(unitl/0.529177)
           
           !APW hnopidihx below is computed but not used

           !c-------------------------------------------------------------
           !c BABY HERE BE CAREFULL. IT TOURNS OUT THAT SOME DIHEDRAL
           !c INVOLVE NO PI ATOMS (I.E. Hydrogen). So use same strategy
           !c     of solvent atoms.....define a noPI_dih list and generate
           !C     Q_noPI_dih(u,p)
           !c     treat this in solvent frame ciao Today BBKing so no stress!!!!
           !c----------------------------------------------------------------
           
           do p=1, npi
              do q=p,npi
                 !c---  Check if external atoms are PI or NoPI --------------
                 
                 !APW modify for pi sulfurs (or general)
                 ispi=.FALSE.
                 do i=1,npispec
                    if(atype(i1).eq.pispec_code(i)) ispi=.true.
                 end do
                 if(ispi) then
                    !if(atype(i1).eq.locate_atom('A')) then
                    ii1=inv_pi(i1)
                    
                    hx(q,p,ii1)=hx(q,p,ii1)&
                         + mocoef(iaa,q)*(factor*ftemp2(1,3,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(1,3,ir))*mocoef(iaa,p)
                    hy(q,p,ii1)=hy(q,p,ii1)&
                         + mocoef(iaa,q)*(factor*ftemp2(2,3,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(2,3,ir))*mocoef(iaa,p)
                    hz(q,p,ii1)=hz(q,p,ii1)&
                         + mocoef(iaa,q)*(factor*ftemp2(3,3,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(3,3,ir))*mocoef(iaa,p)
                    
                 else
                    ii1=inv_dih_nopi(i1)
                    
                    !APW computed but not used
                    !hnopidihx(q,p,ii1)=hnopidihx(q,p,ii1)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(1,3,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(1,3,ir))*mocoef(iaa,p)
                    !hnopidihy(q,p,ii1)=hnopidihy(q,p,ii1)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(2,3,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(2,3,ir))*mocoef(iaa,p)
                    !hnopidihz(q,p,ii1)=hnopidihz(q,p,ii1)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(3,3,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(3,3,ir))*mocoef(iaa,p)
                    
                 end if
                 
                 !APW modify for pi sulfur (or general)
                 ispi=.FALSE.
                 do i=1,npispec
                    if(atype(i2).eq.pispec_code(i)) ispi=.true.
                 end do
                 if(ispi) then
                    !if(atype(i2).eq.locate_atom('A')) then
                    ii2=inv_pi(i2)
                    
                    hx(q,p,ii2)=hx(q,p,ii2)&
                         + mocoef(iaa,q)*(factor*ftemp2(1,4,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(1,4,ir))*mocoef(iaa,p)
                    hy(q,p,ii2)=hy(q,p,ii2)&
                         + mocoef(iaa,q)*(factor*ftemp2(2,4,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(2,4,ir))*mocoef(iaa,p)
                    hz(q,p,ii2)=hz(q,p,ii2)&
                         + mocoef(iaa,q)*(factor*ftemp2(3,4,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(3,4,ir))*mocoef(iaa,p)
                 else
                    
                    ii2=inv_dih_nopi(i2)
                    
                    !APW computed but not used
                    !hnopidihx(q,p,ii2)=hnopidihx(q,p,ii2)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(1,4,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(1,4,ir))*mocoef(iaa,p)
                    !hnopidihy(q,p,ii2)=hnopidihy(q,p,ii2)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(2,4,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(2,4,ir))*mocoef(iaa,p)
                    !hnopidihz(q,p,ii2)=hnopidihz(q,p,ii2)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(3,4,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(3,4,ir))*mocoef(iaa,p)
                    
                 end if
                 
                 hx(q,p,ibb) = hx(q,p,ibb)&
                      + mocoef(iaa,q)*(factor*ftemp2(1,2,ir))*mocoef(ibb,p)&
                      + mocoef(ibb,q)*(factor*ftemp2(1,2,ir))*mocoef(iaa,p)
                 hy(q,p,ibb) = hy(q,p,ibb)&
                      + mocoef(iaa,q)*(factor*ftemp2(2,2,ir))*mocoef(ibb,p)&
                      + mocoef(ibb,q)*(factor*ftemp2(2,2,ir))*mocoef(iaa,p)
                 hz(q,p,ibb) = hz(q,p,ibb)&
                      + mocoef(iaa,q)*(factor*ftemp2(3,2,ir))*mocoef(ibb,p)&
                      + mocoef(ibb,q)*(factor*ftemp2(3,2,ir))*mocoef(iaa,p)
                 
                 hx(q,p,iaa) = hx(q,p,iaa)&
                      + mocoef(iaa,q)*(factor*ftemp2(1,1,ir))*mocoef(ibb,p)&
                      + mocoef(ibb,q)*(factor*ftemp2(1,1,ir))*mocoef(iaa,p)
                 hy(q,p,iaa) = hy(q,p,iaa)&
                      + mocoef(iaa,q)*(factor*ftemp2(2,1,ir))*mocoef(ibb,p)&
                      + mocoef(ibb,q)*(factor*ftemp2(2,1,ir))*mocoef(iaa,p)
                 hz(q,p,iaa) = hz(q,p,iaa)&
                      + mocoef(iaa,q)*(factor*ftemp2(3,1,ir))*mocoef(ibb,p)&
                      + mocoef(ibb,q)*(factor*ftemp2(3,1,ir))*mocoef(iaa,p)
                 
                 !APW modify for pi sulfur or (general)
                 ispi=.FALSE.
                 do i=1,npispec
                    if(atype(i3).eq.pispec_code(i)) ispi=.true.
                 end do
                 if(ispi) then
                    !if(atype(i3).eq.locate_atom('A')) then
                    ii3=inv_pi(i3)
                    
                    hx(q,p,ii3)=hx(q,p,ii3)&
                         + mocoef(iaa,q)*(factor*ftemp2(1,5,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(1,5,ir))*mocoef(iaa,p)
                    hy(q,p,ii3)=hy(q,p,ii3)&
                         + mocoef(iaa,q)*(factor*ftemp2(2,5,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(2,5,ir))*mocoef(iaa,p)
                    hz(q,p,ii3)=hz(q,p,ii3)&
                         + mocoef(iaa,q)*(factor*ftemp2(3,5,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(3,5,ir))*mocoef(iaa,p)
                 else
                    
                    ii3=inv_dih_nopi(i3)
                    
                    !APW computed but not used
                    !hnopidihx(q,p,ii3)=hnopidihx(q,p,ii3)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(1,5,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(1,5,ir))*mocoef(iaa,p)
                    !hnopidihy(q,p,ii3)=hnopidihy(q,p,ii3)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(2,5,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(2,5,ir))*mocoef(iaa,p)
                    !hnopidihz(q,p,ii3)=hnopidihz(q,p,ii3)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(3,5,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(3,5,ir))*mocoef(iaa,p)
                    
                 end if
                 
                 !APW modify for pi sulfur (or general)
                 ispi=.FALSE.
                 do i=1,npispec
                    if(atype(i4).eq.pispec_code(i)) ispi=.true.
                 end do
                 if(ispi) then
                    !if(atype(i4).eq.locate_atom('A')) then
                    ii4=inv_pi(i4)
                    
                    hx(q,p,ii4)=hx(q,p,ii4)&
                         + mocoef(iaa,q)*(factor*ftemp2(1,6,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(1,6,ir))*mocoef(iaa,p)
                    hy(q,p,ii4)=hy(q,p,ii4)&
                         + mocoef(iaa,q)*(factor*ftemp2(2,6,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(2,6,ir))*mocoef(iaa,p)
                    hz(q,p,ii4)=hz(q,p,ii4)&
                         + mocoef(iaa,q)*(factor*ftemp2(3,6,ir))*mocoef(ibb,p)&
                         + mocoef(ibb,q)*(factor*ftemp2(3,6,ir))*mocoef(iaa,p)
                 else
                    
                    ii4=inv_dih_nopi(i4)
                    
                    !APW computed but not used
                    !hnopidihx(q,p,ii4)=hnopidihx(q,p,ii4)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(1,6,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(1,6,ir))*mocoef(iaa,p)
                    !hnopidihy(q,p,ii4)=hnopidihy(q,p,ii4)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(2,6,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(2,6,ir))*mocoef(iaa,p)
                    !hnopidihz(q,p,ii4)=hnopidihz(q,p,ii4)&
                    !     + mocoef(iaa,q)*(factor*ftemp2(3,6,ir))*mocoef(ibb,p)&
                    !     + mocoef(ibb,q)*(factor*ftemp2(3,6,ir))*mocoef(iaa,p)
                    
                 end if
              end do
           end do
        end do
     end if
     !call system_clock ( count3, clock_rate, clock_max )
     !print *,'shit3',(count3-count2)/real(clock_rate)
     
     !*  +---------------------------------------------------------------+
     !*  |  - 0.5*sum(u .ne. v) C_u^q C_v^p P_uv Y'_uv                   |
     !*  +---------------------------------------------------------------+
     do u =1,npi
        do v  = u+1,npi

           dr=(pos(:,pi_iden(u))-pos(:,pi_iden(v)))/a0
           rij=dsqrt(length(dr))
           
           gam2=-bomat(u,v)*gamma_pi(u,v)**2/rij
           th=gam2*dr!(pos(:,pi_iden(u))-pos(:,pi_iden(v)))/a0

           do p=1,npi
              do q=p,npi
                 csuc=0.5d0*mocoef(u,q)*mocoef(v,p) + 0.5d0*mocoef(v,q)*mocoef(u,p)
                 qr(:,q,p,u)=qr(:,q,p,u) - csuc*th(:)
                 qr(:,q,p,v)=qr(:,q,p,v) + csuc*th(:)
              end do
           end do
        end do
     end do
     
     !*  +---------------------------------------------------------------+
     !*  |  sum(u .ne. k) C_u^q C_u^p P_kk Y'_uk                         |
     !*  +---------------------------------------------------------------+
     
     do u =1,npi
        do rho  = u+1,npi
           dr=(pos(:,pi_iden(u))-pos(:,pi_iden(rho)))/a0
           rij=dsqrt(length(dr))
           gam2 = -gamma_pi(u,rho)**2/rij
           th1=bomat(rho,rho)*gam2*dr!(pos(:,pi_iden(u))-pos(:,pi_iden(rho)))/a0
           th2=bomat(u,u)*gam2*dr!(pos(:,pi_iden(u))-pos(:,pi_iden(rho)))/a0
           do p=1,npi
              do q=p,npi
                 csuc1 = mocoef(u,q)*mocoef(u,p)
                 csuc2 = mocoef(rho,q)*mocoef(rho,p)
                 qr(:,q,p,u)=qr(:,q,p,u) + csuc1*th1(:) + csuc2*th2(:)
                 qr(:,q,p,rho)=qr(:,q,p,rho) - csuc1*th1(:) - csuc2*th2(:)
                 
              end do
           end do
        end do
     end do
     
     !*  +---------------------------------------------------------------+
     !*  |  populate the B vector                                        |
     !*  +---------------------------------------------------------------+
     
     do u=1,npi
        do p=1,npi
           do q=p,npi
              qr(1,q,p,u) = qr(1,q,p,u) + hx(q,p,u)
              qr(2,q,p,u) = qr(2,q,p,u) + hy(q,p,u)
              qr(3,q,p,u) = qr(3,q,p,u) + hz(q,p,u)
           end do
        end do
        
        do p=1,npi
           do q=p+1,npi

              qr(:,p,q,u) = qr(:,q,p,u)
              
           end do
        end do
     end do
     !*  +---------------------------------------------------------------+
     !*  |  construct elements of the G matrix:                          |  
     !*  |                                                               |  
     !*  |  G(ai,qp) = (E_p - E_q) - sum(a,virt.) sum(i,real) A(ai,qp)   |  
     !*  |                                                               |  
     !*  |  A(ai,qp) = -0.5 sum(u .ne. v) X_ai^uv C_u^q C_v^p Y_uv       |  
     !*  |             +0.5 sum(u)        X_ai^uu C_u^q C_u^p Y_uu       |  
     !*  |             +    sum(u .ne. k) X_ai^kk C_u^q C_u^p Y_uk       |  
     !*  |                                                               |  
     !*  |  X_ai^uv = 2 * (C_u^i C_v^a + C_u^a C_v^i)                    |  
     !*  +---------------------------------------------------------------+
     
     !print *,'G matrix nine step'
     do u = 1,npi
        do i = 1, norb
           do a = norb+1,npi
              xai(u,a,i) = 2.0d0*(mocoef(u,a)*mocoef(u,i))
           end do
        end do
     end do
     !call system_clock ( count1, clock_rate, clock_max )
     !print *,'shit4',(count1-count3)/real(clock_rate)
!!$         norbit=norb
!!$         do i = 1,norbit
!!$            do a = norbit+1,idim
!!$               xai(u,a,i) = 2.0d0*(mos(u,a)*mos(u,i))
!!$            enddo
!!$         enddo
!!$      enddo

     ind1=0
     do q=1,npi
        do p=q,npi
           do i=1,norb
              do a=norb+1,npi
                 ind1=ind1+1
              end do
           end do
        end do
     end do
     ALLOCATE(index1(4,ind1))
     
     ind1=0
     do q=1,npi
        do p=q,npi
           do i=1,norb
              do a=norb+1,npi
                 ind1=ind1+1
                 index1(1,ind1)=q
                 index1(2,ind1)=p
                 index1(3,ind1)=i
                 index1(4,ind1)=a
              end do
           end do
        end do
     end do
     
     !$OMP parallel do private(q,p,i,a,csuc1,csuc2,tmp1,tmp2,tmp3,xaiuv,xairhorho)
     do n=1,ind1
        q=index1(1,n)
        p=index1(2,n)
        i=index1(3,n)
        a=index1(4,n)
              
        tmp1=0.0d0
        tmp2=0.0d0
        tmp3=0.0d0
        do u=1,npi
           do v=u+1,npi
              csuc1=mocoef(u,q)*mocoef(v,p)
              csuc2=mocoef(u,p)*mocoef(v,q)
              xaiuv = (mocoef(u,a)*mocoef(v,i)+mocoef(v,a)*mocoef(u,i))
              tmp1=tmp1-(csuc1+csuc2)*xaiuv*gamma_pi(u,v)
           end do
        end do
        
        do u=1,npi
           csuc1=mocoef(u,q)*mocoef(u,p)
           xaiuv = xai(u,a,i)
           tmp2=tmp2 + csuc1*xaiuv*gamma_pi(u,u)
        end do
        utmp(a,i,q,p) = utmp(a,i,q,p) + tmp1
        
        do u=1,npi
           csuc1 = 2.0d0*mocoef(u,q)*mocoef(u,p)
           xaiuv = 2.0d0*xai(u,a,i)
           do rho=u+1,npi
              xairhorho = xai(rho,a,i)
              csuc2=mocoef(rho,q)*mocoef(rho,p)
              tmp3 = tmp3 + (csuc1*xairhorho+csuc2*xaiuv)*gamma_pi(rho,u)
           end do
        end do
        utmp(a,i,q,p) = tmp1+tmp2+tmp3
        utmp(a,i,p,q) = utmp(a,i,q,p)
     end do
     !call system_clock ( count2, clock_rate, clock_max )
     !print *,'shit5',(count2-count1)/real(clock_rate)
     !*  +---------------------------------------------------------------+
     !*  |  populate the G matrix                                        |
     !*  +---------------------------------------------------------------+
     !write(*,*)"populating G matrix"
     !$OMP parallel do private(idx1,idx2)
     do q=1,npi
        do p=1,q-1
           do i=1,norb
              do a=norb+1,npi
                 idx1=(a-2)*(a-1)/2 + i
                 idx2 = (q-2)*(q-1)/2+p
                 aea(idx1,idx2) = -utmp(a,i,p,q)
              end do
           end do
        end do
     end do
     
     do idx=1,id2
        aea(idx,idx) = aea(idx,idx) + ener(idxzvecp(idx)) - ener(idxzvecq(idx))
     end do
     !*  +---------------------------------------------------------------+
     !*  |  find inverse of the G matrix                                 |
     !*  +---------------------------------------------------------------+
     
     !*     do i = 1,id2
     !*         do j = 1,id2
     !*             aeainv(i,j) = aea(i,j) 
     !*         enddo
     !*     enddo

     !print *,'matrix inversion 11 step'
     !call system_clock ( count3, clock_rate, clock_max )
     !print *,'shit6',(count3-count2)/real(clock_rate)
     
     call dcopy(id2*id2,aea,1,aeainv,1)  
     call dgetrf(id2,id2,aeainv,id2,ipvt,info)
      
     !call system_clock ( count1, clock_rate, clock_max )
     !print *,'shit7',(count1-count3)/real(clock_rate)
     
     if (info.ne.0) then
        write(*,*) 'matrix inversion problem in dgetrf solvuary.f90',info
         stop
      end if

      call dgetri(id2,aeainv,id2,ipvt,work,lwork,info)
      !call system_clock ( count2, clock_rate, clock_max )
      !print *,'shit8',(count2-count1)/real(clock_rate)
      
      if (info.ne.0) then
         write(*,*) 'matrix inversion problem in dgetri solvuary.f90'
         stop
      end if
      
      
      !call cpu_time(time_f)
      !print *,'solvuary cpu',time_f-time_i
      DEALLOCATE(index1,work)
    END SUBROUTINE SOLVUARY
 
