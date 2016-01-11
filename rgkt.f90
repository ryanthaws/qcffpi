!c------------------------------------------------------------
!c                                                           c
!c  subroutine integration quantum degree of freedom         c
!c  For MEH-PPV internal degree are explicitely treat:       c
!c  We don't neet the extra ter for vibronic motion add      c
!c  a posterio                                               c
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
!* []  DECK RGKT -- RUNGE-KUTTA INTEGRATION OF TULLY COEFFICIENTS
!* =====================================================================

SUBROUTINE rgkt
  
  USE commondata, ONLY: navg,deltat,nnuc,DP
  USE classical, ONLY: fra
  USE quantum, ONLY: iexstate,iexstate_l,ciener,ciener_l,fdv,fdv_l,fexa,ntully,ifinish,dar,deltae,ctully,fdv_m,&
       cieig,cieig_old,ciind,ciind_old,nci,ciener,fexa_l
  USE qcfio, ONLY: gen_io,tul_io,fssh_io,swap_io
  USE nr, ONLY: ran2
  
  
  !* x     initial value of x on input final value of x on output
  !* fxn   1st derivative of x (function)
  !* t     initial time on input -> final time on output   
  !* n     number of time steps
  !* h     time step
  
  !* note that integration is done in program units 
  !* forces,velocities, energies are expressed in program units

  IMPLICIT NONE

  !local variables
  LOGICAL :: tie
  INTEGER :: i,j,k,n,iatom,kk
  REAL(DP) :: thing,t,totprob,tprob,h,efactor,ddnorm,bjk1,bjk2,akk,fdotv,ftemp,max
  REAL :: random
  REAL(DP), DIMENSION(0:ntully) :: addit
  COMPLEX*16, DIMENSION(0:ntully) :: cdot,ctully_l,f1,f2,f3,f4,c1,c2,c3
  REAL(DP), DIMENSION(3,nnuc,0:ntully) :: tfexa
  
  
  !*  +---------------------------------------------------------------+
!*  |  converts energy in au (i.e. hartress) to program units       |
!*  +---------------------------------------------------------------+
  
  iexstate_l=iexstate
  
  !enerfac = 4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !enerfac = enerfac/fkjmol
  efactor = 2625.5 !Hartree to kJ/mol

  fdv=0.0d0
  do i = 0,ntully
     do j = i+1,ntully
        call fdvsub(i,j)
        fdv(j,i) = -fdv(i,j) 
        fdv_l(j,i) = -fdv_l(i,j)
     end do
  end do
  
  fdv_m=0.0d0
  do i=1,ntully
     do j = 1,ntully
        if(fdv(i,j).eq.0 .and. fdv_l(i,j).eq.0) then
           call fdvsub2(i,j)
        end if
     end do
  end do
  
    !APW update old cicoefs
  do i=1,nci
     ciind_old(i,1)=ciind(i,1)
     ciind_old(i,2)=ciind(i,2)
     do j=1,nci
        cieig_old(j,i)=cieig(j,i)
     end do
  end do

  !C JJBL
  !c     call dgemm('T','N',ieveg,ieveg,nci,1.0d0,cieig,nci
  !c    x          ,cieigold,nci,0.0d0,dotm2,ieveg)
  !c
  !c     do i = 1,ntully
  !c      do j = i,ntully
  !c         d_ijvold(i,j)=d_ijv(i,j)
  !c         d_ijvold(j,i)=-d_ijv(i,j)
  !c         d_ijv(i,j)=0.5d0*(dotm2(i,j)-dotm2(j,i))
  !c         d_ijv(j,i)=-d_ijv(i,j)
  !c         write(*,*)' ',i,j,d_ijv(i,j),dotm2(i,j)
  !c      enddo
  !c     enddo
  
  !c     do i = 1,ntully
  !c      do j = i+1,ntully
  !c        fdv(i,j) = d_ijv(i,j)
  !c        fdv(i,j) = -d_ijv(j,i)
  !c        fdvl(i,j) = d_ijvold(i,j)
  !c        fdvl(i,j) = -d_ijvold(j,i)
  !c      enddo
  !c     enddo
  !c     pause
  
  write(gen_io,*) 'rgkt: non adiabatic coupling'

  !c----------------------------------------------------------------------
  !c For MEH-PPV no normal modes terms
  !c----------------------------------------------------------------------
  do i=0,ntully
     addit(i)=0.0d0
  enddo

  !do i = 0,ntully
  !   do j = i,ntully
  !      if (i.ne.j) then
  !         call fdotvout(navg,j,i)
  !      end if
  !   end do
  !end do
      
  write(gen_io,*) 'interpolation done'
  
  !*  +---------------------------------------------------------------+
  !*  |  we are integrating from t to t+dt where dt = 1.0d0 in        |
  !*  |  program units; the initial time origin t is set to zero      |
  !*  |  and the Runge-Kutta time step is dt/n                        |
  !*  +---------------------------------------------------------------+
  !*     RK time step 
  n = 1000
  !h = 1.0d0 / dfloat(n)
  h=1/dfloat(n)
  t = 0.0d0
  !ctully=(0.0d0,0.0d0)
  do k = 1,n
     
     call findcdot(t,cdot,ctully)
     
     do i = 0,ntully
        f1(i) = deltat*h*cdot(i)
        c1(i) = ctully(i) + f1(i)/2.0d0
     end do
     
     call findcdot(t+h/2.0d0,cdot,c1)
     
     do i = 0,ntully
        f2(i) = deltat*h*cdot(i)
        c2(i) = ctully(i) + f2(i)/2.0d0
     end do
     
     call findcdot(t+h/2.0d0,cdot,c2)
     
     do i = 0,ntully
        f3(i) = deltat*h*cdot(i)
        c3(i) = ctully(i) + f3(i)
     end do
     
     call findcdot(t+h,cdot,c3)
     
     do i = 0,ntully
        f4(i) = deltat*h*cdot(i)
     end do
     
     !*  +---------------------------------------------------------------+
     !*  |  f1,f2,f3,f4() are the values of ctully() at intermediate     |
     !*  |  points between t+nh and t + (n+1)h                           |
     !*  +---------------------------------------------------------------+
  
     do i = 0,ntully
        ctully_l(i) = ctully(i) 
        ctully(i) = ctully(i) &
             + (f1(i) + 2.0d0*f2(i) + 2.0d0*f3(i) + f4(i))/6.0d0
     end do
     
     t = dble(k)*h
     
     !write out tully coefficients
     if (mod(k,(n/10)).eq.0) then             
        do i = 0, ntully
           write(tul_io,'(F10.7,2x,I2,2x,4(F22.16,2x))') dble(navg-1)*deltat+deltat*t,i,dreal(ctully(i)),dimag(ctully(i)),dreal(ctully(i)*dconjg(ctully(i))),ciener_l(i)+t*(ciener(i)-ciener_l(i))
        end do
        !write(tul_io,1000) (dble(navg-1)*deltat+t),(ctully(j),j=0,ntully)
        
     end if
1000 format(40(f16.8,1x))
     
     !*  +---------------------------------------------------------------+
     !*  |  integrate the transition probabilities (iexstat - > i)       |
     !*  |  along the way as well using simple mid-point integrator      |
     !*  +---------------------------------------------------------------+
     do i = 0,ntully
        if (i.ne.iexstate) then
           fdotv=(fdv(iexstate,i)-fdv_l(iexstate,i))*(t-h)+fdv_l(iexstate,i)
           bjk1= -2.0d0*dreal(dconjg(ctully_l(i))*ctully_l(iexstate)*fdotv)
           fdotv=(fdv(iexstate,i)-fdv_l(iexstate,i))*t+fdv_l(iexstate,i)
           bjk2= -2.0d0*dreal(dconjg(ctully(i))*ctully(iexstate)*fdotv)
           
           addit(i) = addit(i) + (0.5d0*bjk1 + 0.5d0*bjk2)*h*deltat
           !addit(i) = addit(i) + (0.5d0*bjk(t-h,i,iexstat,ctullyold) &
           !     + 0.5d0*bjk(t,i,iexstat,ctully))*h
           if(i.eq.9) then
              print '(3f15.5)', addit(i),bjk1,bjk2
           end if
        end if
     end do
      
     !**************************************************************************************!
     !APW for non-coupled crossings
     if(k.eq.n/2) then
        i=1
        do i=1,ntully-1
           max=0.0d0
           do j=i+1,ntully
              ftemp=1.0d-09*(dabs(fdv_m(i,j))+dabs(fdv_m(j,i)))*(1-fdv_m(i,i)**2)*(1-fdv_m(j,j)**2)/((ciener(i)-ciener(j)+ciener_l(i)-ciener_l(j))**2) 
              if(ftemp>max) then
                 kk=j
                 max=ftemp
              end if
              
              !write(swap_io,'(3i5,6f18.10)') navg,i,j,ftemp,fdv_m(i,j),fdv_m(j,i),fdv_m(i,i),fdv_m(j,j),(ciener(i)-ciener(j)+ciener_l(i)-ciener_l(j))
           end do
           if(max .gt. 1.0d0) then
              write(swap_io, '(4i5)') navg,i,kk,iexstate
              do j = 1,ntully-1
                 write(swap_io, '(1i5,4f18.10)') j,fdv_l(i,j),fdv(i,j),fdv_l(kk,j),fdv(kk,j)
              end do
              ctully_l(i)=ctully(i)
              ctully(i)=ctully(kk)
              ctully(kk)=ctully_l(i)
              ftemp=addit(i)
              addit(i)=addit(kk)
              addit(kk)=ftemp
              
              !APW need to move fexa too!!
              !tfexa(:,:,:)=fexa(:,:,:,i)
              !fexa(:,:,:,i)=fexa(:,:,:,kk)
              !fexa(:,:,:,kk)=tfexa(:,:,:)
              !tfexa(:,:,:)=fexa(:,:,i,:)
              !fexa(:,:,i,:)=fexa(:,:,kk,:)
              !fexa(:,:,kk,:)=tfexa(:,:,:)
              
              !tfexa(:,:,:)=fexa_l(:,:,:,i)
              !fexa_l(:,:,:,i)=fexa_l(:,:,:,kk)
              !fexa_l(:,:,:,kk)=tfexa(:,:,:)
              !tfexa(:,:,:)=fexa_l(:,:,i,:)
              !fexa_l(:,:,i,:)=fexa_l(:,:,kk,:)
              !fexa_l(:,:,kk,:)=tfexa(:,:,:)
              
              if(iexstate.eq.i) then 
                 iexstate=kk
                 do iatom = 1,nnuc 
                    fra(:,iatom) = fra(:,iatom) + fexa(:,iatom,kk,kk) - fexa(:,iatom,i,i)
                 end do
              else
                 if(iexstate.eq.kk) then
                    iexstate=i
                    do iatom = 1,nnuc 
                       fra(:,iatom) = fra(:,iatom) + fexa(:,iatom,i,i) - fexa(:,iatom,kk,kk)
                    end do
                 end if
              end if
           end if
        end do
     end if
     !**************************************************************************************!
     
  end do !do k=1,n  (end RGKT steps)

  
  
  
  
  !***** make sure the trace of the tully matrix = 1 *****
  thing=0.d0
  do i=0,ntully
     thing=thing+dreal(ctully(i)*dconjg(ctully(i)))
  end do
  write(gen_io,*)' '
  write(gen_io,*)' ',thing,'  dat ting should be 1!!' !MJBH
  write(gen_io,*)' '
    
  write(gen_io,*) 'RK integration complete'
  
  !*  +---------------------------------------------------------------+
  !*  |  set negative transition probabilities  to zero               |
  !*  +---------------------------------------------------------------+
  
  !*     random = ran1(iseed)
  !random = ran2(iseed) 
  !call ran2_s(random)
  call ran2(random)
  !APW change
  !random=-1.0d0
  
  totprob  = 0.0d0
  ifinish = .false.
  tie = .false.
  akk = ctully(iexstate)*dconjg(ctully(iexstate))

  write(gen_io,*) 'currently occupied state',iexstate_l 
  
  
  !c     goto 9981

  !CCCCCCCC MJBH Added on 13.Oct.09 -- renormalization of FSSH
  !CCCCCCCC probabilities for when the transition probabilites
  !CCCCCCCC sum to something greater than one
  do i = 0,ntully
     tprob = addit(i)
     if (tprob.lt.0) tprob = 0.d0
     if (iexstate .ne. i) totprob = totprob + tprob/akk
  end do

  if (totprob .gt. 1) then
     ddnorm = totprob
     write(gen_io,*) "Going to re-normalize tully coefs"
     do i = 0,ntully
        if (iexstate .ne. i) write(gen_io,*) "    ",i,addit(i)/ddnorm/akk
     end do
  else
     ddnorm = 1.d0
  end if
  
  totprob = 0.d0
  
  do i = 0,ntully
     tprob = addit(i)/ddnorm
     if (tprob.lt.0) tprob = 0.0d0
     if (.not.ifinish) then
        if (iexstate.ne.i) then
           if ((random.ge.totprob).and.(random.lt.(totprob+(tprob/(akk))))) then
              iexstate_l = iexstate
              iexstate = i
              ifinish = .true.
              write(gen_io,*) 'my random number was',random
              write(gen_io,*) 'I was between',totprob
              write(gen_io,*) 'and ' , totprob+(tprob/(akk))
              write(gen_io,*) 'for state #',iexstate,iexstate_l
           end if
        end if
        totprob = totprob + tprob/akk
     end if
  end do
  
  do i=0,ntully
     write(fssh_io,3100)navg,iexstate,i,addit(i)/(ddnorm*akk),addit(i)/akk,akk,random
  end do
3100 format(3(i5,1x),4(f16.8,1x))
  
  write(gen_io,*) 'Finished FSSH check'
  !c     stop

  !c9981 continue
  !c     write(*,*)'Skipping FSSH check... code - 9981'


  !*  +---------------------------------------------------------------+
!*  |  if we make a transition then adjust the velocities along     |  
!*  |  the nonadiabatic counpling vector to conserve energy         |  
!*  +---------------------------------------------------------------+
  if (ifinish) then
     if (iexstate.eq.iexstate_l) then 
        write(gen_io,*) 'DABOO You ask the impossible master'
        write(gen_io,*) ' I quvit'
        stop
     end if
     
     !*  +---------------------------------------------------------------+
     !*  |  convert energy to program units                              |
     !*  +---------------------------------------------------------------+ 
     if (iexstate_l.eq.0) then
        deltaE = (-ciener(iexstate))
     else
        if (iexstate.eq.0) then
           deltaE = ciener(iexstate_l)
        else
           deltaE = (-ciener(iexstate)+ciener(iexstate_l))
        end if
     end if
          
     
     !*  +---------------------------------------------------------------+
     !*  |  get the nonadiabatic coupling vectors                        |
     !*  +---------------------------------------------------------------+ 
     !c        if (solv) then
     !c        do iatom = 1,natoms
     !c           dsx(iatom) = fexx(iatom,iexstat,iexstatold)/deltaE
     !c           dsy(iatom) = fexy(iatom,iexstat,iexstatold)/deltaE
     !c           dsz(iatom) = fexz(iatom,iexstat,iexstatold)/deltaE
     !c        enddo
     !c         
     !c        if (deltaE .eq. 0.d0) then
     !c         write(*,*)'0 s  ',dsx(iatom),fexx(iatom,iexstat,iexstatold)
     !c         write(*,*)'0 s  ',dsy(iatom),fexy(iatom,iexstat,iexstatold)
     !c         write(*,*)'0 s  ',dsz(iatom),fexz(iatom,iexstat,iexstatold)
     !cc        pause
     !c        endif
     !c        endif
     
     
     do iatom = 1,nnuc
        dar(:,iatom) = fexa(:,iatom,iexstate,iexstate_l)/(deltaE*efactor)
        !write(gen_io,*) "shit",iatom,iexstate,iexstate_l,dar(1,iatom)
     end do
     
     
     
     if (deltaE .eq. 0.d0) then
        write(gen_io,*)'0 a  ',dar(1,iatom),fexa(1,iatom,iexstate,iexstate_l)
        write(gen_io,*)'0 a  ',dar(2,iatom),fexa(2,iatom,iexstate,iexstate_l)
        write(gen_io,*)'0 a  ',dar(3,iatom),fexa(3,iatom,iexstate,iexstate_l)
        !c         pause
     end if
     !*  +---------------------------------------------------------------+
     !*  |  apply constraint conditions so that the nonadiabatic         |  
     !*  |  coupling vectors have no components along the constraints    |  
     !*  +---------------------------------------------------------------+
     !cccccc   if(solv) call fixsolvent(dsx,dsy,dsz)
         
  end if
      
      
END SUBROUTINE rgkt
