  !c----------------------------------------------------------------------
  !c
  !c     this section calculates the nonbond and electrostatic energy
  !c     contributions to e,d,dd.
  !c
  !c----------------------------------------------------------------------
SUBROUTINE het_nonbon(eng)
  !subroutine het_nonbon(d_,x_,e,iac2_ppv,iblo_ppv,wnb_ppv)
  
  USE commondata, ONLY:nnuc
  USE classical
  USE qcfio, ONLY:gen_io
  
  IMPLICIT NONE
  
  !c     input 
  DOUBLE PRECISION, INTENT(out) :: eng
!!$      real*8 d_(1),x_(1),e
!!$c      integer iac(1),iblo(1),wnb(1)
!!$      integer iac2_ppv(1),iblo_ppv(1)
!!$      integer wnb_ppv(mxbond,4)

  !c     local
  INTEGER, DIMENSION(2) :: ind
  INTEGER :: i,j,k,l,btype,ncex,hit,iatom,jatom,ibeg,i1,i2,i3,i4,iind1,iind2
  DOUBLE PRECISION :: rsq,r,f,df,ddf,charge,r1i,r2i,r6i
  DOUBLE PRECISION, DIMENSION(3) :: dx,ftemp,atemp
  DOUBLE PRECISION, DIMENSION(3,2) :: dr
  DOUBLE PRECISION, DIMENSION(3,3,3) :: ddr
  INTEGER, DIMENSION(2,nnuc*nnuc/2) :: nbind
  DOUBLE PRECISION, DIMENSION(nnuc*nnuc/2) :: find
  INTEGER :: nb1,nb2
  
  do i=1,3
     dr(i,1)=0.0
     dr(i,2)=0.0
     do j=1,3
        do k=1,3
           ddr(k,j,i)=0.0
        end do
     end do
  end do

  nb1=0
  do iind1=1,nnuc-1
     do iind2=iind1+1,nnuc
        if(nonbs(iind1,iind2).eq.1) then
           nb1=nb1+1
           nbind(1,nb1)=iind1
           nbind(2,nb1)=iind2
        end if
     end do
  end do

  !do iind1=1,nnuc-1
  !ind(1)=iind1
  !!call het_prnbf2(iatom) !APW old code?
  !do iind2=iind1+1,nnuc
  !ind(2)=iind2
  !if(nonbs(ind(1),ind(2)).eq.1) then

  !$OMP parallel do private(ind,dx,rsq,r,r1i,dr,ddr,btype,hit,i,ftemp,charge,r2i,r6i,atemp,f,df,ddf,ncex)
  do nb2=1,nb1
     ind(1)=nbind(1,nb2)
     ind(2)=nbind(2,nb2)
     
     
     do k=1,3
        dx(k)=pos(k,ind(1))-pos(k,ind(2))
     end do
     rsq=length(dx)
     !rsq=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
     r=dsqrt(rsq)
     r1i=1.d0/r
     call difbond(dx,dr,ddr,r)
     !c--- The electrostatic part is computed in PPV module ----------
     !c      call het_nbfun2(iatom,jatom,bi)
     !c      e=e+f
     !c      call drd(i3,j3,newton,df,dr,d_,ddr)
     !c----------------------------------------------------------------
     
     !call het_nbfun1(iatom,jatom,bi,iac2_ppv)
     !subroutine het_nbfun1(i1,j1,bi,iac2_ppv)
     
     ncex=0
     
     btype=pair_code(atype(ind(1)),atype(ind(2)))
     hit=0
     i=0
     do while(i.lt.nnonb_type .and. hit.eq.0)
        i=i+1
        if(btype.eq.nonb_code(i)) hit=i
     end do
     if(hit.eq.0) then
        !APW previosly set hit=1
        write(gen_io,*) '***warning: nonbond pair',atype(ind(1)),&
             atype(ind(2)),'is missing from the nonbonding parameters***'
        stop
     end if
     
     ftemp(1)=nonb_pms(1,hit)
     ftemp(2)=nonb_pms(2,hit)
     ftemp(3)=nonb_pms(3,hit)
     
     charge=qatom(atype(ind(1)))*qatom(atype(ind(2)))
     if(r.lt.1) r1i=1.0
     r2i=r1i*r1i
     r6i=r2i*r2i*r2i
     atemp(1)=r1i**ncex
     atemp(2)=ftemp(2)*dexp(-ftemp(3)*r)*atemp(1) ! A*exp(-\mu*R_{i,j})
     atemp(3)=-ftemp(3)*atemp(2)
     
     f=atemp(2)-ftemp(1)*r6i+charge*r1i
     find(nb2)=f
     df=atemp(3)+(6.d0*ftemp(1)*r6i-charge*r1i-ncex*atemp(2))*r1i
     ddf=-ftemp(3)*atemp(3)-2.d0*ncex*atemp(3)*r1i&
          -(42.d0*ftemp(1)*r6i-2.d0*charge*r1i-ncex*(ncex+1)*atemp(2))*r2i
     
     !eng=eng+f
     call drd(ind,df,ddf,dr,ddr)
     
     !end if
     !end do
  end do
  
  
  do nb2=1,nb1
     eng=eng+find(nb2);
  end do

END SUBROUTINE het_nonbon

