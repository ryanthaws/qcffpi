!c------------------------------------------------------------
!c                                                           c
!c  subroutine adding the forces in exc. state               c
!c             for the solute                                c
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


SUBROUTINE ADD_FORCE_EX(cit,cib,gamma_pi)
  USE commondata, ONLY: npi,DP,norb,a0,nnuc,navg
  USE classical, ONLY : length,pos,pi_iden
  USE quantum, ONLY: fexa,ntully,nci,neig,mocoef,ciener,qr,aeainv,utmp,cieig
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  !input variables
  INTEGER, DIMENSION(norb*(npi-norb)), INTENT(in) :: cit,cib
  !REAL(DP), DIMENSION(nci,nci), INTENT(in) :: cieig
  REAL(DP), DIMENSION(npi,npi), INTENT(in) :: gamma_pi
  
  !local variables
  LOGICAL :: test1,test2
  LOGICAL, DIMENSION(nci) :: topbotlog
  LOGICAL, DIMENSION(nci*(nci-1)/2) :: tbmnlog
  INTEGER :: i,j,k,idxsize,b,iatom,idx,idxf,idxpnt,idz,ifake,ifind,isaw,itully,jtully,&
       jdz,l,m,n,norbit,p,q,u,v,idxsize2,idxsize3
  INTEGER, ALLOCATABLE, DIMENSION(:) :: idxpntf,idxpoint
  REAL(DP), DIMENSION(3) :: dr
  REAL(DP), DIMENSION(3,npi) :: fr
  REAL(DP), DIMENSION(npi,npi) :: rij
  REAL(DP), DIMENSION(npi,nci) :: bot,top
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: botn,topn,botm,topm
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: xea,xes
  REAL(DP), DIMENSION(npi,npi,0:ntully,0:ntully) :: epsz
  REAL(DP) :: sqrt2,ftmp,cnm,forcefac,th,th1,th2,tmp,yth,yth1,yth2,zth,zth1,zth2
  REAL(DP) :: time_f,time_i,ftemp1a,ftemp1b,ftemp2a,ftemp2b
  INTEGER, DIMENSION(4,nci*nci*(ntully+1)*(ntully+1)/4) :: ciind
  INTEGER, DIMENSION(2,(ntully+1)*(ntully+1)/2) :: ciind1
  INTEGER, DIMENSION(2,nci*nci/2,(ntully+1)*(ntully+1)/2) :: ciind2
  INTEGER, DIMENSION((ntully+1)*(ntully+1)/2) :: cinum2
  INTEGER :: cinum1
  INTEGER :: cinum,tind,t2ind
  INTEGER :: count1,count2,count3,clock_rate,clock_max
  
  call system_clock ( count1, clock_rate, clock_max )
  
  forcefac=2625.5/a0 !from Hartree/bohr radius to kJ/(mol*angstrom)
  
  sqrt2=dsqrt(2.0d0)
  idxsize=nci*nci*ntully*ntully
  idxsize2=npi*(npi-1)/2
  idxsize3=nci*(nci-1)/2
  ALLOCATE(idxpntf(idxsize),idxpoint(idxsize))
  ALLOCATE(botn(npi,idxsize3),topn(npi,idxsize3),botm(npi,idxsize3),topm(npi,idxsize3))
  ALLOCATE(xea(idxsize2,0:ntully,0:ntully),xes(idxsize2,0:ntully,0:ntully))
  
  topbotlog = .true.
  bot=0.0d0
  top=0.0d0
  !do i = 1,nci
  !   topbotlog(i) = .FALSE.
  !end do
  
  tbmnlog = .true.
  botn=0.0d0
  topn=0.0d0
  botm=0.0d0
  topm=0.0d0
  !do i = 1,nci*(nci-1)/2 
  !   tbmnlog(i) = .FALSE.
  !end do
  
  call cpu_time(time_i)
  
  !APW determine distances between pi-atoms
  do i = 1,npi
     do j = i+1,npi
        dr = (pos(:,pi_iden(i)) - pos(:,pi_iden(j)))/a0
        rij(i,j) = dsqrt(length(dr))
        rij(j,i) = rij(i,j)
     end do
  end do
  
  idxpnt = 0                        
  idxf = 0 
  
  
  epsz = 0.0d0;
  
  !*  +---------------------------------------------------------------+
  !*  |  construct elements of the Langrangian matrix:                |
  !*  |                                                               |
  !*  |  X_ij = sum(m,MO) Q_jm H_im                                   |
  !*  |       + 2 * sum(m;k;l,MO) G_jmkl (im|kl)                      |  
  !*  +---------------------------------------------------------------+
  !*  |  note that matrix elements are not symmetric anymore          |  
  !*  |  cieig(M,itully)*cieig(N,jtully)                              |  
  !*  |                   != cieig(M,jtully)*cieig(N,itully)          |  
  !*  +---------------------------------------------------------------+
  cinum=0
  cinum1=0
  cinum2=0
  do itully = 1, ntully
     do jtully = itully,ntully
        cinum1=cinum1+1
        ciind1(1,cinum1)=itully
        ciind1(2,cinum1)=jtully
        do n=1,nci
           test1=.false.
           if(dabs(cieig(n,itully)*cieig(n,jtully)) .gt. 1.0d-04) then
              test1=.true.
              cinum=cinum+1
              ciind(1,cinum)=itully
              ciind(2,cinum)=jtully
              ciind(3,cinum)=n
              cinum2(cinum1)=cinum2(cinum1)+1
              ciind2(1,cinum2(cinum1),cinum1)=n
           end if
           if(test1 .and. topbotlog(n)) then
              topbotlog=.false.
              bot(:,n) = 0.0d0
              top(:,n) = 0.0d0
              
              !$OMP parallel do private(ftemp1a,ftemp2a)
              do l = 1,npi
                 ftemp1a=0.0d0
                 ftemp2a=0.0d0
                 do u =1,npi
                    do v = 1,npi
                       !bot(l,n) = bot(l,n)&
                       ftemp1a = ftemp1a&
                            + gamma_pi(u,v)*(mocoef(u,l)*mocoef(v,cit(n))&
                            *(2.0d0*mocoef(u,cit(n))*mocoef(v,cib(n))&
                            - mocoef(u,cib(n))*mocoef(v,cit(n)))&
                            + mocoef(u,cib(n))*mocoef(v,cit(n))&
                            *(2.0d0*mocoef(u,cit(n))*mocoef(v,l)&
                            - mocoef(u,l)*mocoef(v,cit(n))))
                       
                       !top(l,n) = top(l,n)&
                       ftemp2a = ftemp2a&    
                            + gamma_pi(u,v)*(mocoef(u,cib(n))*mocoef(v,l)&
                            *(2.0d0*mocoef(u,cit(n))*mocoef(v,cib(n))&
                            - mocoef(u,cib(n))*mocoef(v,cit(n)))&
                            + mocoef(u,cib(n))*mocoef(v,cit(n))&
                            * (2.0d0*mocoef(u,l)*mocoef(v,cib(n))&
                            - mocoef(u,cib(n))*mocoef(v,l)))
                    end do
                 end do
                 bot(l,n) = bot(l,n) + ftemp1a
                 top(l,n) = top(l,n) + ftemp2a
              end do
           end if
        end do
     end do
  end do
  
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'ex1',(count2-count1)/real(clock_rate)
    
  !$OMP parallel do private(itully,jtully,n,cnm)
  do l = 1,npi
     do tind=1,cinum
        itully = ciind(1,tind)
        jtully = ciind(2,tind)
        n = ciind(3,tind)
        
        !*  +-----------------------------------------------------------+
        !*  |  calculating the weird gamma stuff for diagonal elements  |        
        !*  |  of the CI matrix                                         | 
        !*  +-----------------------------------------------------------+
        cnm = cieig(n,jtully)*cieig(n,itully)
        
        epsz(l,cib(n),itully,jtully) = epsz(l,cib(n),itully,jtully)&
             - cnm*bot(l,n)
        epsz(l,cit(n),itully,jtully) = epsz(l,cit(n),itully,jtully)&
             - cnm*top(l,n)
     end do
  end do
  
  !*  +---------------------------------------------------------------+
  !*  |  1. upper fock eigenvalue                                     |
  !*  |  for the derivatives of the bond order matrix elements        |
  !*  +---------------------------------------------------------------+
  !*  +---------------------------------------------------------------+
  !*  |  2. lower fock eigenvalue                                     |
  !*  |  for the derivatives of the bond order matrix elements        |
  !*  +---------------------------------------------------------------+
  
  !$OMP parallel do private(itully,jtully,n,cnm)
  do tind=1,cinum
     itully = ciind(1,tind)
     jtully = ciind(2,tind)
     n = ciind(3,tind)
     do j = 1,norb
        do b = norb+1,npi
           epsz(b,j,itully,jtully) = epsz(b,j,itully,jtully)&
                + cnm*(utmp(b,j,cib(n),cib(n)) - utmp(b,j,cit(n),cit(n)))
        end do
     end do
  end do
  
  !call system_clock ( count3, clock_rate, clock_max )
  !print *,'ex2',(count3-count2)/real(clock_rate)
  
     !*  +---------------------------------------------------------------+
     !*  |  CIS matrix: elements are A_nm                                |
     !*  |                                                               |
     !*  |  < Psi(i->k) | H_el - V_pi | Psi(j->l) >                      |
     !*  |              = (E_k - E_i) delta(i,j) delta(k,l)              |
     !*  |              + sum(u;v) C_u^j C_v^k ( 2 C_u^l C_v^i           |
     !*  |              - C_u^i C_v^l) Y_uv                              |
     !*  +---------------------------------------------------------------+
     
     !*  +---------------------------------------------------------------+
     !*  |  Gradient of A_nm:                                            |
     !*  |                                                               |
     !*  |  A_nm^(R) = sum(u;v) C_u^j C_v^k ( 2 C_u^l C_v^i              |  
     !*  |           - C_u^i C_v^l) (d Y_uv / dR)                        |  
     !*  +---------------------------------------------------------------+
  
  
  !do tind=1,cinum
  !   itully = ciind(1,tind)
  !   jtully = ciind(2,tind)
  !   n = ciind(3,tind)
  !$OMP parallel do private(itully,jtully,n,fr)
  do tind=1,cinum1
     itully=ciind1(1,tind)
     jtully=ciind1(2,tind)
     do t2ind=1,cinum2(tind)
        n=ciind2(1,t2ind,tind)
        fr=0.0d0
        
        !APW terms in here only depend on n
        do u =1,npi
           do v = 1,npi
              if (u.ne.v) then
                 ftmp = mocoef(u,cib(n))*mocoef(v,cit(n))& 
                      * (2.0d0*mocoef(u,cit(n))*mocoef(v,cib(n))&
                      - mocoef(u,cib(n))*mocoef(v,cit(n)))&
                      * gamma_pi(u,v)**2/rij(u,v)
                 
                 !Hartree/Bohr
                 fr(:,u) = fr(:,u) + ftmp*(pos(:,pi_iden(u))-pos(:,pi_iden(v)))/a0
                 fr(:,v) = fr(:,v) - ftmp*(pos(:,pi_iden(u))-pos(:,pi_iden(v)))/a0
                 
              end if
           end do
        end do
        
        !APW terms in here only depend on n
        do iatom =1,npi 
           fr(:,iatom) = fr(:,iatom)&
                -qr(:,cit(n),cit(n),iatom)&
                +qr(:,cib(n),cib(n),iatom)  
        end do
        
        !*  +---------------------------------------------------------------+
        !*  |  convert back to the other number ordering and add in         |  
        !*  |  the forces weighted by the CI coefficients:                  |
        !*  |                                                               |  
        !*  |  sum(m;n) C_m^I C_n^I A_nm^(R)                                |  
        !*  +---------------------------------------------------------------+
        
        do iatom = 1,npi
           fexa(:,pi_iden(iatom),itully,jtully)=fexa(:,pi_iden(iatom),itully,jtully)&
                + cieig(n,itully)*cieig(n,jtully) * forcefac*fr(:,iatom)
        end do
        
     end do
  end do
     
  !call system_clock ( count1, clock_rate, clock_max )
  !print *,'ex3',(count1-count3)/real(clock_rate)
  
  cinum=0
  !*********************************compute some matrix elements ******************************!
  !call cpu_time(time_f)
  !print *,'shit1a',time_f-time_i
  !time_i=time_f
  cinum1=0
  cinum2=0
  do itully = 1, ntully
     do jtully = itully,ntully
        cinum1=cinum1+1
        ciind1(1,cinum1)=itully
        ciind1(2,cinum1)=jtully
        do n=1,nci-1
           do m = n+1,nci
              idx = nci*(n-1) - n*(n-1)/2 + (m-n)
              test1=.false.
              if(dabs(cieig(n,itully)*cieig(m,jtully)) .gt. 1.0d-04) then
                 test1=.true.
                 cinum=cinum+1
                 ciind(1,cinum)=itully
                 ciind(2,cinum)=jtully
                 ciind(3,cinum)=n
                 ciind(4,cinum)=m
                 cinum2(cinum1)=cinum2(cinum1)+1
                 ciind2(1,cinum2(cinum1),cinum1)=n
                 ciind2(2,cinum2(cinum1),cinum1)=m
              else
                 if(dabs(cieig(n,jtully)*cieig(m,itully)) .gt. 1.0d-04) then
                    test1=.true.
                    cinum=cinum+1
                    ciind(1,cinum)=itully
                    ciind(2,cinum)=jtully
                    ciind(3,cinum)=n
                    ciind(4,cinum)=m
                    cinum2(cinum1)=cinum2(cinum1)+1
                    ciind2(1,cinum2(cinum1),cinum1)=n
                    ciind2(2,cinum2(cinum1),cinum1)=m
                 end if
              end if
              if(test1 .and. tbmnlog(idx)) then
                 
                 tbmnlog(idx)=.false.
                 botn(:,idx) = 0.0d0
                 topn(:,idx) = 0.0d0
                 botm(:,idx) = 0.0d0
                 topm(:,idx) = 0.0d0
                 
                 !$OMP parallel do private(ftemp1a,ftemp1b,ftemp2a,ftemp2b)
                 do l = 1,npi
                    ftemp1a=0.0d0
                    ftemp1b=0.0d0
                    ftemp2a=0.0d0
                    ftemp2b=0.0d0
                    do u =1,npi
                       do v = 1,npi
                          
                          !botm(l,idxf) = botm(l,idxf)&
                          ftemp1a=ftemp1a&
                               +gamma_pi(u,v)*(mocoef(u,l)*mocoef(v,cit(n))&
                               *(2.0d0*mocoef(u,cit(m))*mocoef(v,cib(n))&
                               -mocoef(u,cib(n))*mocoef(v,cit(m))))
                          
                          !botn(l,idxf) = botn(l,idxf)& 
                          ftemp1b=ftemp1b&
                               + gamma_pi(u,v)*(mocoef(u,cib(m))*mocoef(v,cit(n))&
                               *(2.0d0*mocoef(u,cit(m))*mocoef(v,l)&
                               - mocoef(u,l)*mocoef(v,cit(m))))
                          
                          !topm(l,idxf) = topm(l,idxf)& 
                          ftemp2a=ftemp2a&
                               + gamma_pi(u,v)*(mocoef(u,cib(m))*mocoef(v,cit(n))&
                               *(2.0d0*mocoef(u,l)*mocoef(v,cib(n))&
                               -mocoef(u,cib(n))*mocoef(v,l)))
                          
                          !topn(l,idxf) = topn(l,idxf)& 
                          ftemp2b=ftemp2b&
                               + gamma_pi(u,v)*(mocoef(u,cib(m))*mocoef(v,l)&
                               *(2.0d0*mocoef(u,cit(m))*mocoef(v,cib(n))&
                               -mocoef(u,cib(n))*mocoef(v,cit(m))))
                       end do
                    end do
                    botm(l,idx) = botm(l,idx)+ftemp1a
                    botn(l,idx) = botn(l,idx)+ftemp1b
                    topm(l,idx) = topm(l,idx)+ftemp2a
                    topn(l,idx) = topn(l,idx)+ftemp2b
                 end do
              end if
           end do
        end do
     end do
  end do
  
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'ex4',(count2-count1)/real(clock_rate)
  !call cpu_time(time_f)
  !print *,'shit1b',time_f-time_i,idxf
  !time_i=time_f
  
  !$OMP parallel do private(cnm,itully,jtully,n,m,idx)
  do tind=1,cinum1
     itully=ciind1(1,tind)
     jtully=ciind1(2,tind)
     do t2ind=1,cinum2(tind)
        n=ciind2(1,t2ind,tind)
        m=ciind2(2,t2ind,tind)
        
        idx = nci*(n-1) - n*(n-1)/2 + (m-n)
        !*  +---------------------------------------------------------------+
        !*  |  for the off-diagonal CI matrix elements                      |
        !*  |  calculating the weird gamma stuff for diagonal elements of   |  
        !*  |  the CI matrix                                                |  
        !*  +---------------------------------------------------------------+
        
        cnm = cieig(n,itully)*cieig(m,jtully)+cieig(n,jtully)*cieig(m,itully)
        
        do l = 1,npi
           epsz(l,cib(n),itully,jtully)=epsz(l,cib(n),itully,jtully)& 
                - cnm*botn(l,idx)
           
           epsz(l,cit(n),itully,jtully)=epsz(l,cit(n),itully,jtully)& 
                - cnm*topn(l,idx)
           
           epsz(l,cib(m),itully,jtully)=epsz(l,cib(m),itully,jtully)& 
                - cnm*botm(l,idx)
           
           epsz(l,cit(m),itully,jtully)=epsz(l,cit(m),itully,jtully)& 
                - cnm*topm(l,idx)
        end do
     end do
  end do

  !call system_clock ( count3, clock_rate, clock_max )
  !print *,'ex5',(count3-count2)/real(clock_rate)
  
  !call cpu_time(time_f)
  !print *,'shit1c',time_f-time_i,cinum,nci*nci*ntully*ntully/4
  !time_i=time_f
  !*  +---------------------------------------------------------------+
  !*  |  Gradient of A_nm:                                            |
  !*  |                                                               |
  !*  |  A_nm^(R) = sum(u;v) C_u^j C_v^k ( 2 C_u^l C_v^i              |
  !*  |           - C_u^i C_v^l) (d Y_uv / dR)                        |
  !*  +---------------------------------------------------------------+
  
  
  !
  !!$OMP parallel do private(itully,jtully,n,m,fr,ftmp)
  !do tind=1,cinum
  !   itully=ciind(1,tind)
  !   jtully=ciind(2,tind)
  !   n=ciind(3,tind)
  !   m=ciind(4,tind)
  !$OMP parallel do private(itully,jtully,n,m,fr,ftmp)
  do tind=1,cinum1
     itully=ciind1(1,tind)
     jtully=ciind1(2,tind)
     do t2ind=1,cinum2(tind)
        n=ciind2(1,t2ind,tind)
        m=ciind2(2,t2ind,tind)
        fr = 0.0d0
        
        do u =1,npi
           do v = 1,npi
              if (u.ne.v) then
                 
                 !APW only depends on n and m
                 ftmp = mocoef(u,cib(m))*mocoef(v,cit(n))& 
                      *(2.0d0*mocoef(u,cit(m))*mocoef(v,cib(n))&
                      -mocoef(u,cib(n))*mocoef(v,cit(m)))&
                      * gamma_pi(u,v)**2/rij(u,v)
                 
                 fr(:,u) = fr(:,u) + ftmp*(pos(:,pi_iden(u))-pos(:,pi_iden(v)))/a0
                 fr(:,v) = fr(:,v) - ftmp*(pos(:,pi_iden(u))-pos(:,pi_iden(v)))/a0
                 
              end if
           end do
        end do
        
        !*  +---------------------------------------------------------------+
        !*  |  convert back to the other number ordering and add in         |
        !*  |  the forces weighted by the CI coefficients:                  |
        !*  |                                                               |
        !*  |  sum(m;n) C_m^I C_n^I A_nm^(R)                                |
        !*  +---------------------------------------------------------------+
        do iatom = 1,npi
           fexa(:,pi_iden(iatom),itully,jtully)=fexa(:,pi_iden(iatom),itully,jtully)& 
                + cieig(n,itully)*cieig(m,jtully)*forcefac*fr(:,iatom)&
                + cieig(n,jtully)*cieig(m,itully)*forcefac*fr(:,iatom)
        end do
        
     end do
  end do
  
  !call system_clock ( count1, clock_rate, clock_max )
  !print *,'ex6',(count1-count3)/real(clock_rate)
  
  !*  +---------------------------------------------------------------+
  !*  |  enddos for the tully indices; now add in the unitary stuff   |
  !*  +---------------------------------------------------------------+
  
  !*  +---------------------------------------------------------------+
  !*  |  additional CIS NA coupling component  <Psi_I|dPsi_J/dR>      |  
  !*  |                                                               | 
  !*  |  forces are in program units; we multiply here by the         |  
  !*  |  transition energy in atomic units; we will convert the       |  
  !*  |  force back to atomic units and divide by an energy in        |  
  !*  |  atomic units to get the nonadibatic coupling vector          |  
  !*  +---------------------------------------------------------------+
  
  !call cpu_time(time_f)
  !print *,'shit1',time_f-time_i,idxf
  !time_i=time_f
  
  do itully = 1,ntully
     do m = 1,nci
        if (dabs(cieig(m,itully)).gt.1.0d-04) then
           th = ciener(itully)*cieig(m,itully)*sqrt2
           do u = 1,npi
              yth = th*mocoef(u,cib(m))
              do l = 1,npi
                 zth = yth*mocoef(u,l)
                 epsz(l,cit(m),0,itully)=epsz(l,cit(m),0,itully) + zth
              end do
           end do
        end if
     end do
  end do
  
  !call cpu_time(time_f)
  !print *,'shit2',time_f-time_i
  !time_i=time_f
  do itully = 1,ntully
     do jtully = itully+1,ntully
        
        do m = 1,nci
           do n = m+1,nci
              
              test1 = dabs(cieig(m,jtully)*cieig(n,itully)).gt.1.0d-04                 
              test2 = dabs(cieig(m,itully)*cieig(n,jtully)).gt.1.0d-04  
              
              if (test1 .or. test2) then 
                 
                 if ((cib(m).eq.cib(n)).and.(cit(m).eq.cit(n))) then
                    
                    continue
                    
                 else
                    
                    !*  +---------------------------------------------------------------+
                    !*  |  MISTAKE HERE OCTOBER 8 1997 ---> SHOULD HAVE USED - OF THE   |
                    !*  |  BELOW EXPRESSIONS when the virtual  spin - orbitals are      |
                    !*  |  equal; WHY? BECAUSE <PSI_i->a|d/dR PSI_j->a> = <j|d/dR i>    |
                    !*  |  makes no difference to to the <PSI_0 d/dR |PSI_excited>      |
                    !*  |  NA coupling; so probably a minor mistake                     |
                    !*  +---------------------------------------------------------------+
                    
                    if (cit(m).eq.cit(n)) then
                       
                       th1 = (ciener(jtully)-ciener(itully))&
                            *cieig(m,itully)*cieig(n,jtully)*sqrt2
                       
                       th2 = (ciener(jtully)-ciener(itully))&
                            *cieig(n,itully)*cieig(m,jtully)*sqrt2
                       
                       do u = 1,npi
                          
                          yth1 = th1*mocoef(u,cib(m))
                          yth2 = th2*mocoef(u,cib(n))
                          
                          do l = 1,npi
                             
                             zth1 = yth1*mocoef(u,l)
                             
                             epsz(l,cib(n),itully,jtully)=epsz(l,cib(n),itully,jtully) + zth1
                             
                             zth2 = yth2*mocoef(u,l)
                             
                             epsz(l,cib(m),itully,jtully)=epsz(l,cib(m),itully,jtully) + zth2
                             
                          end do
                       end do
                    end if
                    
                    if (cib(m).eq.cib(n)) then
                       
                       th1 = (ciener(jtully)-ciener(itully))&
                            *cieig(m,itully)*cieig(n,jtully)*sqrt2
                       
                       th2 = (ciener(jtully)-ciener(itully))&
                            *cieig(n,itully)*cieig(m,jtully)*sqrt2
                       
                       do u = 1,npi
                          
                          yth1 = th1*mocoef(u,cit(m))
                          yth2 = th2*mocoef(u,cit(n))
                          
                          do l = 1,npi
                             
                             zth1 = yth1*mocoef(u,l)
                             
                             epsz(l,cit(n),itully,jtully)=epsz(l,cit(n),itully,jtully) + zth1
                             
                             zth2 = yth2*mocoef(u,l)
                             
                             epsz(l,cit(m),itully,jtully)=epsz(l,cit(m),itully,jtully) + zth2
                          end do
                       end do
                    end if
                 end if
              end if
           end do
        end do
        
     end do
  end do
  
  !call system_clock ( count2, clock_rate, clock_max )
  !print *,'ex7',(count2-count1)/real(clock_rate)
  
  !*  +---------------------------------------------------------------+
  !*  |  Z vector method:  X G_inv = Z                                |  
  !*  |                                                               |  
  !*  |  total CI force:                                              |  
  !*  |      sum(m;n) C_m^I C_n^I d/dR A_nm =                         |  
  !*  |      sum(m;n) C_m^I C_n^I A_nm^(R) + 2 Z B^(R)                |  
  !*  +---------------------------------------------------------------+
  
  !call cpu_time(time_f)
  !print *,'shit3',time_f-time_i
  !time_i=time_f
  do itully = 0,ntully
     do jtully = itully,ntully
        
        do q=1,npi
           do p = 1,q-1
              
              idx = (q-2)*(q-1)/2 + p
              
              xea(idx,itully,jtully) = epsz(q,p,itully,jtully)-epsz(p,q,itully,jtully)
           end do
        end do
        
        !$OMP parallel do
        do idz = 1, npi*(npi-1)/2
           
           xes(idz,itully,jtully) = 0.0d0
           do jdz = 1, npi*(npi-1)/2
              
              xes(idz,itully,jtully) = xes(idz,itully,jtully)+aeainv(idz,jdz)*xea(jdz,itully,jtully)
              
           end do
        end do
        
        !$OMP parallel do private(idx)
        do iatom = 1,npi
           
           do q=1,npi
              do p = 1,q-1
                 
                 idx = (q-2)*(q-1)/2 + p
                 
                 fexa(:,pi_iden(iatom),itully,jtully)=fexa(:,pi_iden(iatom),itully,jtully)& 
                      + forcefac*xes(idx,itully,jtully)*qr(:,q,p,iatom)
                 
              end do
           end do
        end do
     end do
  end do
  
  !call system_clock ( count3, clock_rate, clock_max )
  !print *,'ex8',(count3-count2)/real(clock_rate)
  
  !APW?
  !*     +---------------------------------------------------------------+
  !*     |  should this be done here??? Need to check back on this!!     |  
  !*     +---------------------------------------------------------------+
  
  !call cpu_time(time_f)
  !print *,'shit4',time_f-time_i
  !time_i=time_f
  do itully = 0,ntully
     do jtully = itully,ntully
        
        do iatom = 1,npi
           
           fexa(:,pi_iden(iatom),jtully,itully)=fexa(:,pi_iden(iatom),itully,jtully)
        end do
     end do
  end do
  !stop
  !call cpu_time(time_f)
  !print *,'shit5',time_f-time_i
  !time_i=time_f
  
  !call system_clock ( count1, clock_rate, clock_max )
  !print *,'ex9',(count1-count3)/real(clock_rate)
  
  DEALLOCATE(idxpntf,idxpoint,botn,topn,botm,topm,xea,xes)  
  
  !call cpu_time(time_f)
  !print *,'add_force_ex cpu',time_f-time_i

  
END SUBROUTINE ADD_FORCE_EX
