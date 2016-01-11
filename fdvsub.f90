!* =====================================================================
!* []  DECK FDVSUB  -- CALCULATE THE NONADIABATIC COUPLING VECTOR
!* =====================================================================

SUBROUTINE fdvsub(i,j)
  
  USE commondata, ONLY: DP,nnuc,navg,npi,deltat,norb
  USE classical, ONLY: vel,vel_l
  USE quantum, ONLY: ciener,fexa,fexa_l,ciener_l,fdv,fdv_l,deltae,cieig,cieig_old,&
       ciind,ciind_old,nci!,mocoef,mocoef_old
  USE qcfio, ONLY: gen_io,fdotv_io
  

  IMPLICIT NONE
  
  !input variables
  INTEGER, INTENT(IN) :: i,j

  !local variables
  INTEGER :: iatom,n,m,l,u
  REAL(DP), DIMENSION(3) :: fa,va
  REAL(DP) :: dot,efactor,sum1,sum2,sum3,ftemp1,ftemp2,sq2,dtn
  REAL(DP), DIMENSION(npi,npi) :: cid,cian,ciao
  
  sq2=sqrt(2.0d0)
  if (i.eq.j) then 
     write(gen_io,*) 'DABOO You ask the impossible master'
     write(gen_io,*) ' I quvit'
     stop
  end if
      
!************* FIND FDV FOR CURRENT TIME STEP ***************************
!*  +---------------------------------------------------------------+
!*  |  converts energy in au (i.e. hartress) to program units       |
!*  +---------------------------------------------------------------+
      !enerfac = 4.3598d-18*1.0d-3/dfloat(molcs)*avogad
      !enerfac = enerfac/fkjmol
  efactor = 2625.5   !APW keep energy in units of Hartree?
  
!*  +---------------------------------------------------------------+
!*  |  convert energy to program units                              |
!*  +---------------------------------------------------------------+
  if (i.eq.0) then
     deltae = -ciener(j)
  else
     if (j.eq.0) then
        deltae = ciener(i)
     else
        deltae = (ciener(i)-ciener(j))
     endif
  end if
  !deltaE = deltaE*efactor
      
  dot = 0.0d0

!*     do the solvent first
!c     if (solv) then
!c     do iatom = 1,natoms
!c         fsx = fexx(iatom,j,i)
!c         fsy = fexy(iatom,j,i)
!c         fsz = fexz(iatom,j,i)
!c         vsx = vx(iatom,np)
!c         vsy = vy(iatom,np)
!c         vsz = vz(iatom,np)
!c         dot = dot + fsx*vsx + fsy*vsy + fsz*vsz
!c     enddo
!c     endif

!*     now for ppv (solute)
      
  do iatom = 1,nnuc
     fa(:)=fexa(:,iatom,j,i)
     !fax = fexxa(iatom,j,i)
     !fay = fexya(iatom,j,i)
     !faz = fexza(iatom,j,i)
     va(:)=vel(:,iatom)
     !vax = vxa(iatom)
     !vay = vya(iatom)
     !vaz = vza(iatom)
     dot = dot + fa(1)*va(1) + fa(2)*va(2)+ fa(3)*va(3)
     !print *,'shit3a',navg,i,j,iatom
     !print *,'shit3b',navg,fa(1),fa(2),fa(3),va(1),va(2),va(3)
     if(j.eq.10 .and. i.eq.9) then
        print '(1i5,6f15.5)',iatom,fa(1),fa(2),fa(3),va(1),va(2),va(3)
     end if
  end do
  !print *,'shit3c',navg,dot,deltaE*efactor
  fdv(i,j) = dot / (deltaE*efactor)
  if(j.eq.10 .and. i.eq.9) then
     print *,fdv(i,j),dot
  end if
  if (dabs(deltaE) .le. 1.d-15) then
     write(gen_io,*)'yo! mtv sucks',i,j,dot,deltaE,fdv(i,j)
     fdv(i,j)=0.d0
  end if
  
  
  !************* FIND FDV FOR PREVIOUS TIME STEP **************************
  !*  +---------------------------------------------------------------+
  !*  |  convert energy to program units                              |
  !*  +---------------------------------------------------------------+
  if (i.eq.0) then
     deltaE = -ciener_l(j)
  else
     if (j.eq.0) then
        deltaE = ciener_l(i)
     else
        deltaE = (ciener_l(i)-ciener_l(j))
     end if
  end if
  !deltaE = deltaE*efactor
  
  dot = 0.0d0
  
  !*      do the solvent
  !c      if (solv) then
  !c      do iatom = 1,natoms
  !c         fsx = fexxl(iatom,j,i)
  !c         fsy = fexyl(iatom,j,i)
  !c         fsz = fexzl(iatom,j,i)
  !c         vsx = vxl(iatom,np)
  !c         vsy = vyl(iatom,np)
  !c         vsz = vzl(iatom,np)
  !c         dot = dot + fsx*vsx + fsy*vsy + fsz*vsz
  !c      enddo
  !c      endif
  
  !*     do the solute (ppv)
  do iatom = 1,nnuc
     fa(:) = fexa_l(:,iatom,j,i)
     !fax = fexxal(iatom,j,i)
     !fay = fexyal(iatom,j,i)
     !faz = fexzal(iatom,j,i)
     va(:) = vel_l(:,iatom)
     !vax = vxal(iatom)
     !vay = vyal(iatom)
     !vaz = vzal(iatom)
     dot = dot + fa(1)*va(1) + fa(2)*va(2)+ fa(3)*va(3)
     !if(j.eq.10 .and. i.eq.9) then
     !   print '(1i5,6f15.5)',iatom,fa(1),fa(2),fa(3),va(1),va(2),va(3)
     !end if
  end do
  
  fdv_l(i,j) = dot / (deltae*efactor)
  
  !c	fdvl(i,j) = 0.d0
  !c	fdv(i,j) = 0.d0
  
  if (dabs(deltaE) .le. 1.d-15) then
     write(*,*)'yo! mtvl sucks',i,j,deltaE,dot,fdv_l(i,j)
     fdv_l(i,j)=0.d0
     !c      pause
  end if
  
  !fdv_m(i,j)=(fdv(i,j)+fdv_l(i,j))/2.0d0
  !if(fdv(i,j)**2.eq.0 .and. fdv_l(i,j)**2.eq.0) then
  !   cid=0.0d0
  !   cian=0.0d0
  !   ciao=0.0d0
  !   if(i.ne.0 .and. j.ne.0) then
  !      do n=1,nci
  !         cid(ciind(n,1),ciind(n,2))=cieig(n,i)
  !         cian(ciind(n,1),ciind(n,2))=cieig(n,j)
  !         ciao(ciind_old(n,1),ciind_old(n,2))=cieig_old(n,j)
  !      end do
        
        !dtn=deltat/1000.0d0
        
  !      sum1=0.0d0
  !      sum2=0.0d0
  !      sum3=0.0d0
  !      do n=1,norb
  !         do m=norb+1,npi
  !            sum1=sum1+cid(n,m)*(cian(n,m)-ciao(n,m))/deltat
  !         end do
  !      end do
        
        !do n=1,norb
        !   do l=1,norb
        !      do m=norb+1,npi
        !         ftemp1=cid(n,m)*cian(l,m)
        !         do u=1, npi
        !            ftemp2=sq2*mocoef(u,l)*(mocoef(u,n)-mocoef_old(u,n))/deltat
        !            sum2=sum2+ftemp1*ftemp2;
        !         end do
        !      end do
        !   end do
        !end do
        
        !do n=1,norb
        !   do l=norb+1,npi
        !      do m=norb+1,npi
        !         ftemp1=cid(n,l)*cian(n,m)
        !         do u=1, npi
        !            ftemp2=sq2*mocoef(u,l)*(mocoef(u,m)-mocoef_old(u,m))/deltat
        !            sum3=sum3+ftemp1*ftemp2;
        !         end do
        !      end do
        !   end do
        !end do
        
  !      fdv_m(i,j)=sum1+sum2+sum3
  !      
  !   end if
  !end if
     
  write(fdotv_io,'(3i5,3f25.15)') navg,i,j,fdv(i,j),fdv_l(i,j)
  
END SUBROUTINE fdvsub
      
SUBROUTINE fdvsub2(i,j)
  
  USE commondata, ONLY: DP,norb,npi,navg,deltat
  USE quantum, ONLY: fdv,fdv_l,fdv_m,cieig,cieig_old,ciind,ciind_old,nci!,mocoef,mocoef_old
  USE qcfio, ONLY: gen_io,swap_io
  

  IMPLICIT NONE
  
  !input variables
  INTEGER, INTENT(IN) :: i,j

  !local variables
  INTEGER :: iatom,n,m,l,u
  REAL(DP), DIMENSION(3) :: fa,va
  REAL(DP) :: dot,efactor,sum1,sum2,sum3,sum4,ftemp1,ftemp2,sq2,dtn
  REAL(DP), DIMENSION(npi,npi) :: cid,cian,ciao
  
  fdv_m(i,j)=0.0d0
  sq2=dsqrt(2.0d0)
  
  cid=0.0d0
  cian=0.0d0
  ciao=0.0d0
  
  do n=1,nci
     cid(ciind(n,1),ciind(n,2))=cieig(n,i)
     cian(ciind(n,1),ciind(n,2))=cieig(n,j)
     ciao(ciind_old(n,1),ciind_old(n,2))=cieig_old(n,j)  
  end do
  
  
  sum1=0.0d0
  sum2=0.0d0
  sum3=0.0d0
  sum4=0.0d0
  do n=1,norb
     do m=norb+1,npi
        sum1=sum1+cid(n,m)*(cian(n,m)-ciao(n,m))/deltat
        sum2=sum2+cian(n,m)*ciao(n,m)
     end do
  end do
  
  !do n=1,norb
  !   do l=1,norb
  !      do m=norb+1,npi
  !         ftemp1=cid(n,m)*cian(l,m)
  !         do u=1, npi
  !            ftemp2=sq2*mocoef(u,l)*(mocoef(u,n)-mocoef_old(u,n))/deltat
  !            sum3=sum3+ftemp1*ftemp2;
  !         end do
  !      end do
  !   end do
  !end do
  
  !do n=1,norb
  !   do l=norb+1,npi
  !      do m=norb+1,npi
  !         ftemp1=cid(n,l)*cian(n,m)
  !         do u=1, npi
  !            ftemp2=sq2*mocoef(u,l)*(mocoef(u,m)-mocoef_old(u,m))/deltat
  !            sum4=sum4+ftemp1*ftemp2;
  !         end do
  !      end do
  !   end do
  !end do
  
  
  !fdv_m(i)=sum1*(1-sum2*sum2)
  fdv_m(i,j)=sum1
  if(i.eq.j) fdv_m(i,i)=sum2
  
  
  !write(swap_io,'(3i5,4f25.15)') navg,i,j,sum1,sum2,sum3,sum4
  !write(fdotv_io,'(3i5,3f25.15)') navg,i,j,fdv(i,j),fdv_m(i,j),fdv_l(i,j)

END SUBROUTINE fdvsub2
