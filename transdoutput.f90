!c------------------------------------------------------------
!c                                                           c
!c  subroutine computing the transition dipole               c
!c                                                           c
!c  For MEH-PPV RAB is zero because involves only C atoms    c
!c                                                           c
!c-----------------------------------------------------------c
!c                                                           c
!c   Modified from original subroutine of JL                 c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm/utexas.edu                                  c
!c   Modified further by M. Bedard-Hearn 2008.Jan.15         c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c

SUBROUTINE transdoutput(cit,cib,cimatdiag)
  
  USE commondata, ONLY: DP,norb,npi,a0
  USE classical, ONLY: pi_iden,length,pos
  USE quantum, ONLY: osc,osc_vec,inv_pi,nci,neig,ciener,mocoef,deltaab,cieig
  USE qcfio, ONLY: gen_io
  USE PPV, ONLY: junction,ntor_junc,torsion

  IMPLICIT NONE
  
  !c neig is the number of eigenvalues and eigen vectors to find.
  !c        it is defined in pppc.h, can be as high as NCI
  
  !input variables
  INTEGER, DIMENSION(norb*(npi-norb)), INTENT(IN) :: cit,cib
  REAL(DP), DIMENSION(norb*(npi-norb)), INTENT(IN) :: cimatdiag
  !REAL(DP), DIMENSION(nci,neig), INTENT(IN) :: cieig
  
  !local variables
  INTEGER :: i,j,k,ii,jj,kk,ibottom,itop,iconfig,ioutput,imain,u,v,ib1,ib2
  REAL(DP), DIMENSION(npi,npi) :: s12
  REAL(DP), DIMENSION(3,nci) :: transr,transd
  REAL(DP), DIMENSION(nci) :: f_0r,f_0del
  REAL(DP), DIMENSION(3) :: dr,mr,md
  REAL(DP) :: rdist,s12cc,sumd,sumr,totald,totalr,rab,rabcc,biggest,ftemp

  !do i=1,nat_sol
  !rposx(i) = xna(i)*unitl/0.529177
  !rposy(i) = yna(i)*unitl/0.529177
  !rposz(i) = zna(i)*unitl/0.529177
  !enddo
  
  !set up the transition AO matrix elements             
  !c---------------------------------------------------------------------
  !c RAB is <phi_n|r|phi_m>=d is the dipole moment in the reference system
  !c center in the bond n-m. 
  !c For same atom n and m Rab=0. 
  !c So for C-C RAC=0 and is different for C-N and C-O.  
  !c
  !c Now for MEH-PPV we don't have different atom contributing to PI system
  !c For MEH-PPV RAB is always Zero. We don't add any extra contribution
  !c to first loop
  !c
  !c S12 no used. Just don't add term S12 R_mn^o. This because must be zero
  !c for ortonormality
  !c 
  !c DELTA for MEH-PPV is defined before (check if done or not)
  !c---------------------------------------------------------------------
  s12=0.0e0
  rab=0.0e0

  !c--- Diagonal Elements ----------------------------------
  do i = 1,npi
     s12(i,i) = 1.0e0
     !deltaab(i,i) = 0.0e0
     !rab(i,i) = 0.0e0
  end do

  S12CC = 0.d0
  RABCC = 0.d0
  
  !c-----------------------------------------------------------------c
  !c    Note for MEH-PPV all are CC type                             c
  !c    I fill all the matrix element here                           c  
  !c    no need to extra copy                                        c
  !c-----------------------------------------------------------------c

  do i=1,npi
     do j=1,npi
        if(i.ne.j) S12(i,j)=S12CC
     enddo
  enddo

  do i=1,ntor_junc
     ii=inv_pi(junction(1,i))
     jj=inv_pi(junction(2,i))
     s12(ii,jj)=s12(ii,jj)*torsion(i)
     deltaab(ii,jj)=deltaab(ii,jj)*torsion(i)
  end do
  
  !c--- Delta Matrix ------------------------------------------------c
  
  !C nci === number of excited state configs being included.         
  !C neig === number of excited states that we bother to calculate
  !C          eigenvectors for 
  !C 1st get the matrix element <a|o|r> where a is the "bottom" MO
  !C                                          r is the "top"    MO
  !C we evaluate the transition dipole element with o = r operator
  
  
  do ii = 1,nci
     transr(:,ii) = 0.0
     !c------------------------------------------------------------
     !c cib and cit are defined in qcffsol, as well as in ppp_ppv.
     !c Here, we use the one defined in ppp_ppv after sorting the
     !c CI matrix in ascending order (using sort2)...they are sorted
     !c such that the first pair corresponds to the smallest element
     !c on the diagonal of the CIS excitation matrix.  cib stands for
     !c the BOTTOM end of a single CI excitation, and cit stands for
     !c the TOP end, where BOTTOM and TOP refer to the MO that you will
     !c be emptying and the TOP refers to the MO that the electron is
     !c being promoted to.  the benefit of doing the sorting of pairs
     !c this way is that the pairs are effectively sorted in terms of
     !c their relative "importance" to the lowest-lying excited states
     !c (eg., the first cib-cit pair describes the homo -> lumo
     !c single excitation, which has the smallest 1-e energy gap
     !c & often has the largest CI coefficient in the linear expansion)
     !c
     !c dimensions of cib and cit are defined in pppc.h
     !c------------------------------------------------------------
     
     do u = 1,npi          !cycles over pi-carbon atoms
        kk=pi_iden(u)
        transr(:,ii) = transr(:,ii)+mocoef(u,cib(ii))*mocoef(u,cit(ii))*pos(:,kk)/a0
     end do
     
     !C we evaluate the transition dipole element with o = d /dr operator
     
     transd(:,ii) = 0.0
     
     do u = 1,npi
        do v = u+1,npi
           kk=pi_iden(u)
           jj=pi_iden(v)
           ftemp=mocoef(u,cib(ii))*mocoef(v,cit(ii))-mocoef(u,cit(ii))*mocoef(v,cib(ii))
           dr=(pos(:,jj)-pos(:,kk))/a0
           rdist=sqrt(length(dr))
           transd(:,ii)=transd(:,ii)+ftemp*(pos(:,jj)-pos(:,kk))*deltaab(u,v)/(a0*rdist)
!!$                  transdx(ii) = transdx(ii) 
!!$     &                 + mos(u,ib(ii))*mos(v,it(ii))*(rposx(jj)
!!$     $                 -rposx(kk))/rdist*deltaab(u,v)- mos(u,it(ii))
!!$     $                 *mos(v,ib(ii))*(rposx(jj)-rposx(kk))/rdist
!!$     $                 *deltaab(u,v)
!!$                  transdy(ii) = transdy(ii) 
!!$     &                 + mos(u,ib(ii))*mos(v,it(ii))*(rposy(jj)
!!$     $                 -rposy(kk))/rdist*deltaab(u,v)- mos(u,it(ii))
!!$     $                 *mos(v,ib(ii))*(rposy(jj)-rposy(kk))/rdist
!!$     $                 *deltaab(u,v) 
!!$                  transdz(ii) = transdz(ii) 
!!$     &                 + mos(u,ib(ii))*mos(v,it(ii))*(rposz(jj)
!!$     $                 -rposz(kk))/rdist*deltaab(u,v)- mos(u,it(ii))
!!$     $                 *mos(v,ib(ii))*(rposz(jj)-rposz(kk))/rdist
!!$     $                 *deltaab(u,v)
        end do
     end do
  end do     !do ii = 1,nci

  !C choose the # of states we wish to output data on.

  !c        ioutput = 5
  ioutput = nci
  
  do jj=1,nci
     do ii=1,3
        osc_vec(ii,jj)=0.d0
     end do
  end do
  
  do iconfig = 1,neig
     !C     now get the transition moment between ground and iconfig excited CI state

     biggest = 0.0
     
     mr=0.d0
     md=0.d0
!!$            MxR = 0.0
!!$            MyR = 0.0
!!$            MzR = 0.0

!!$            Mxd = 0.0
!!$            Myd = 0.0
!!$            Mzd = 0.0

     if (iconfig.lt. 5) then
        write(gen_io,*) ' ******************************************'
        write(gen_io,*) 'transition dipole moment between',iconfig,'th CI state and ground state'
        write(gen_io,*) ' ******************************************'
     end if

     do ii = 1,nci
        
        mr(:)=mr(:)+dsqrt(2.0d0)*cieig(ii,iconfig)*transr(:,ii)
!!$               MxR = MxR 
!!$     &              + sqrt(2.0)*cieig(ii,iconfig)*transRx(ii)
!!$               MyR = MyR 
!!$     &              + sqrt(2.0)*cieig(ii,iconfig)*transRy(ii)
!!$               MzR = MzR 
!!$     &              + sqrt(2.0)*cieig(ii,iconfig)*transRz(ii)
        
        md(:) = md(:) + dsqrt(2.0d0)*cieig(ii,iconfig)*transd(:,ii)
!!$               Mxd = Mxd + sqrt(2.0)*cieig(ii,iconfig)*transdx(ii)
!!$               Myd = Myd + sqrt(2.0)*cieig(ii,iconfig)*transdy(ii)
!!$               Mzd = Mzd + sqrt(2.0)*cieig(ii,iconfig)*transdz(ii)
        
        if (abs(cieig(ii,iconfig)).gt.biggest) then

           ibottom = cib(ii)
           itop = cit(ii)
           imain = ii

           biggest = abs(cieig(ii,iconfig))
        end if

        if (iconfig.lt.5) then
           if (abs(cieig(ii,iconfig)).gt.0.5) then
              write(gen_io,*) ' configuration # ' ,ii  
              write(gen_io,*) 'a - > r ' ,cib(ii),cit(ii)
              write(gen_io,*) ' contribution ' , cieig(ii,iconfig)
           end if
        end if
     end do

     mr=mr*(ciener(iconfig))
     !APW HERE
!!$            MxR = MxR * (ciener(iconfig))
!!$            MyR = MyR * (ciener(iconfig))
!!$            MzR = MzR * (ciener(iconfig))
     
     !APW not sure what the initial values of the sr array is here
     !sr=sr*(ciener(iconfig))
!!$            SxR = SxR * (ciener(iconfig))
!!$            SyR = SyR * (ciener(iconfig))
!!$            SzR = SzR * (ciener(iconfig))

     totalr=mr(1)**2+mr(2)**2+mr(3)**2
     totald=md(1)**2+md(2)**2+md(3)**2
!!$            totalR = (MxR**2+MyR**2+MzR**2)
!!$            totald = (MxD**2+MyD**2+MzD**2)

     f_0r(iconfig) = 2.0/(3.0*ciener(iconfig))*totalR !what's this 2/3 business?
!!$           f_0R(iconfig) =  2.0/(3.0*ciener(iconfig)) !what's this 2/3 business?
!!$     &           *totalR  
     f_0del(iconfig) =2.0/(3.0*ciener(iconfig))*totald
!!$            f_0del(iconfig) =  2.0/(3.0*ciener(iconfig))
!!$     &           *totald  
     
     if(iconfig.lt.5) then
        write(gen_io,*) '<PSI_0|r|PSI_CI^',iconfig,'> X Y Z'
        write(gen_io,*) mr(1),mr(2),mr(3)
        write(gen_io,*) 'FROM the',imain,' config. ',ibottom,'->',itop
        
        write(gen_io,*) transr(1,imain)*sqrt(2.0)*cimatdiag(imain),&
             transr(2,imain)*sqrt(2.0)*cimatdiag(imain),&
             transr(3,imain)*sqrt(2.0)*cimatdiag(imain)
        !APW not sure what sr is?
        !write(gen_io,*) 'SX comp.'
        !write(gen_io,*) sr(1),sr(2),sr(3)
        write(gen_io,*) '<PSI_0|d /dr|PSI_CI^',iconfig,'> X Y Z'
        write(gen_io,*) md(1),md(2),md(3)
        write(gen_io,*) 'FROM the',imain,' config. ',ibottom,'->',itop
        write(gen_io,*) transd(1,imain)*sqrt(2.0),&
             transd(2,imain)*sqrt(2.0),&
             transd(3,imain)*sqrt(2.0)

        write(gen_io,*) 'ratio of f_r / f_del ' ,totalR/totalD
        write(gen_io,*) 'energy of ',iconfig,' state=',ciener(iconfig)
        write(gen_io,*) 'energy of ',imain,' config mat=',cimatdiag(imain)
        write(gen_io,*) 'f_0',iconfig,'(R) = ' , f_0R(iconfig)
        write(gen_io,*) 'f_0',iconfig,'(del) = ' ,f_0del(iconfig)
     end if
     
     osc_vec(:,iconfig)=mr !components of the transition dipole moment in Hartree*bohr
  end do      !do iconfig = 1,ieveg

  !*     now check the Thomas-Kuhn sum rule
  sumD = 0.0
  sumR = 0.0
  do iconfig = 1,nci
     sumD = sumD + f_0del(iconfig)
     sumR = sumR + f_0R(iconfig)
  end do
     
  write(gen_io,*) ' T-K sum rule for r operator ' , sumR
  write(gen_io,*) ' T-K sum rule for del operator ' , sumD

  do i=1,nci
     osc(i)=f_0R(i)   !oscillator strength units are
  enddo                   !Hartree*bohr^2

END SUBROUTINE transdoutput
