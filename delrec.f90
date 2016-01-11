!APW this routine is not called within the version of the code I recieved. I'll try to code it up just in case but I am unlikely to error check.

!c--------------------------------------------------------------------
!c---
!c--- del re sigma charges                                   ---------
!c--- pocket of routines:
!c--- delrec , fordet , ludcmp , lubksb , chafo 
!c--- 
!c--------------------------------------------------------------------
!
!c---------------------------------------------------------------------
!c     
!c     this subroutine calculates the sigma charges by the del re method
!c     
!c---------------------------------------------------------------------
SUBROUTINE delrec()
!!$      subroutine delrec(natom,iaco,iac,cha,gadr,do,delta,dha,do1,eadr
!!$     $     ,ilink,iprint)
  
  USE commondata, ONLY: nnuc
  USE classical, ONLY: ilink,links,qatom,eadr,gadr,atype,new_do
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  !local variables
  INTEGER :: il,im,icod,jcod,iatom,jatom,imi,ikot,ii,ll
  INTEGER :: i,j,k,m,ip,imax,ikon,ifin
  INTEGER, DIMENSION(3) :: itemp
  INTEGER, DIMENSION(nnuc) :: indx
  DOUBLE PRECISION :: d,tiny,sum,r_arg,aamax
  DOUBLE PRECISION, DIMENSION(2) :: ftemp
  DOUBLE PRECISION, DIMENSION(nnuc) :: do1,vv
  DOUBLE PRECISION, DIMENSION(nnuc,nnuc) :: odelta,delta
  
  tiny=1.0d-20
  
!!$c     Form matrix delta = 1 - gamma
  !call fordet(delta,iaco,iac,natom,do1,gadr,do,ilink) 
  !c------------------------------------------------------------------
  !
  !     find the atoms involved in the sigma-charge determination of
  !     a given atom.
  !
  !c------------------------------------------------------------------
  itemp(1)=34345
  itemp(2)=231
  ikot=0
  
  delta=0.000000001d0 !matrix(nnuc,nnuc)

  do iatom=1,nnuc !1
     il=ilink(iatom)
     do m=1,il !2
        ikot=ikot+1
        im=4*(m-1)
     end do
     icod=atype(iatom)
     do1(iatom)=-new_do(icod)
     r_arg=iatom*77.
     !delta(iatom,iatom)=-1.d0+(rand(itemp(1),itemp(2))*0.001)
     !APW this seems to replace random number generator, doesn't seem kosher
     delta(iatom,iatom)=-1.d0+(dcos(r_arg)*0.001)
     itemp(1)=il*4 !ifin
     itemp(2)=0 !iammin
     itemp(3)=0 !icarb
     do imi=1,itemp(1) !3
        if(links(imi,ikot).ne.0) then
           jatom=links(imi,ikot)
           jcod=atype(jatom)
           if(jcod.eq.1)  itemp(2)=itemp(2)+1
           if(jcod.ne.1) itemp(3)=itemp(3)+1
           delta(iatom,jatom)=gadr(icod,jcod)
        end if
     end do
     
     if(icod.eq.3.and.itemp(2).eq.2) do1(iatom)=-new_do(15)
     if(icod.eq.5.and.itemp(3).eq.3) do1(iatom)=-0.07
  end do

  odelta=delta
  
  d=1.0
  do i=1,nnuc
     ftemp(1)=0.0 !aamax
     do j=1,nnuc
        if (abs(delta(i,j)).gt.ftemp(1)) ftemp(1)=abs(delta(i,j))
     end do
     if(ftemp(1).eq.0) then
        write(gen_io,*) "***warning: singular matrix ***"
        stop
     end if
     vv(i)=1./ftemp(1)
  end do
  
  do j=1,nnuc
     if(j.gt.1) then
        do i=1,j-1
           sum=delta(i,j)
           if(i.gt.1) then
              do k=1,i-1
                 sum=sum-delta(i,k)*delta(k,j)
              end do
              delta(i,j)=sum
           end if
        end do
     end if
     ftemp(1)=0.0
     do i=j,nnuc
        sum=delta(i,j)
        if(j.gt.1) then
           do k=1,j-1
              sum=sum-delta(i,k)*delta(k,j)
           end do
           delta(i,j)=sum
        end if
        ftemp(2)=vv(i)*abs(sum) !dum
        if(ftemp(2).ge.ftemp(1)) then
           imax=i
           aamax=ftemp(2)
        end if
     end do
     if(j.ne.imax) then
        do k=1,nnuc
           ftemp(2)=delta(imax,k)
           delta(imax,k)=delta(j,k)
           delta(j,k)=ftemp(2)
        end do
        d=-d
        vv(imax)=vv(j)
     end if
     indx(j)=imax
     if(j.ne.nnuc) then
        if(delta(j,j).eq.0.) delta(j,j)=tiny
        ftemp(2)=1./delta(j,j)
        do i=j+1,nnuc
           delta(i,j)=delta(i,j)*ftemp(2)
        end do
     end if
  end do
  if(delta(nnuc,nnuc).eq.0.) delta(nnuc,nnuc)=tiny
  
  !call lubksb(delta,natom,idim_q,indx,do1)
  !subroutine lubksb(a,n,np,indx,b)
  
  ii=0
  do i=1,nnuc
     ll=indx(i)
     sum=do1(ll)
     do1(ll)=do1(i)
     if(ii.ne.0) then
        do j=ii,i-1
           sum=sum-delta(i,j)*do1(j)
        end do
     else
        if(sum.ne.0) ii=1
     end if
     do1(i)=sum
  end do
  do i=nnuc,1,-1
     sum=do1(i)
     if(i.lt.nnuc) then
        do j=i+1,nnuc
           sum=sum-delta(i,j)*do1(j)
        end do
     end if
     do1(i)=sum/delta(i,i)
  end do

  !determine sigma charges from dha
  !call chafo(do1,cha,natom,iaco,iac,eadr,ilink)
  !c------------------------------------------------------------------
  !c
  !c      unpack the atoms involved in the sigma charge determination
  !c      of a given atom
  !c
  !c------------------------------------------------------------------
  !subroutine chafo(dha,cha,v_natom,v_wlin,v_iac,v_eadr,v_ilink)

  ikon=0
  do iatom=1,nnuc !1
     icod=atype(iatom)
     il=ilink(iatom)
     !do ip=1,il !2
     ikon=ikon+1
     ifin=il*4
     qatom(iatom)=0.
     if(species(atype(iatom)).eq.'Z') qatom(iatom)=2.
     !if(atype(iatom).eq.14) qatom(iatom)=2.
     do ip=1,ifin
        if(links(ip,ikon).ne.0) then
           jatom=links(ip,ikon)
           jcod=atype(jatom)
           ftemp(1)=(do1(jatom)-do1(iatom))/(2.*eadr(icod,jcod))
           qatom(iatom)=qatom(iatom)+ftemp(1)/dsqrt(1.0+ftemp(1)**2)
        end if
     end do
  end do
END SUBROUTINE delrec
