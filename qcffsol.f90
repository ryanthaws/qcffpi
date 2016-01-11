! QCFFSOL calculates the classical PE (e_ppv) and forces (d_ppv).
!------------------------------------------------------------
!                                                           !
!  subroutine calculates the classical PE(e_ppv) and        !
!  and forces (d_ppv).                                      !
!                                                           !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
!   Revised by Adam P. Willard 2010                         !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!  code is very bulk, breaking it into seperate subfiles    ! 
!     qcfblkdta.f90 (now all within qcf module)             ! 
!     npar.f90                                              ! 
!     molein.f90                                            ! 
!     read3d.f90                                            ! 
!                                                           ! 
!-----------------------------------------------------------!
!APW break this up so the first stuff is in a seperate subroutine
SUBROUTINE qcffsol(eqcff,first)
  !subroutine qcffsol(first,x_ppv,d_ppv,ib_ppv,jb_ppv,nap_ppv,nc_ppv
  !   $     ,nbond_ppv,nphi_ppv,wp_ppv,ntheta_ppv,wt_ppv,nat_sol
  !   $     ,sigma_ppv,iac_ppv,e_ppv,iac2_ppv,wb_ppv,wnb_ppv,wlin_ppv
  !   $     ,ilink_ppv,icg_ppv,iblo_ppv,use_delre,acode_ppv
  !   $     ,ineigh_ppv,iachem_ppv,istep)
  
  USE commondata, ONLY: nnuc,DP
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  !input variables
  REAL(DP), INTENT(OUT) :: eqcff
  LOGICAL, INTENT(IN) :: first
  
  !local variables
  INTEGER :: i,j,k,iprint,istep
  
  !APW moved to setup_vars.f90
  !if(first) then
  !call npar(1)           !0: no parameters printing
  !call molein()
  !first=.false.
  !end if

  call molecu(eqcff)
  
END SUBROUTINE qcffsol
  
!!$c----------------------------------------------------------------------
!!$c
!!$c     this subroutine builds the off diagonal matrix elements  of h
!!$c     which  are of the form f(b)*g(phi)
!!$c     the contributions from the s**2 correction for orthogonalization
!!$c     are  constructed in the r matrix and in the diagonal part of the
!!$c      h matrix .  all the corresponding energy contributions
!!$c     and their derivatives are also constructed
!!$c
!!$c----------------------------------------------------------------------
!!$      subroutine het_vphi(ap,i1,j1,ih,jh,lp,ic)
!!$      implicit none
!!$      real*8 length
!!$      include 'qcf.h'
!!$      include 'qcfio.h'
!!$c:::  input
!!$      integer i1,j1,ih,jh
!!$      integer lp,ic
!!$
!!$      real*8 ap
!!$c:::  local
!!$      integer i,i3
!!$      integer j,j3
!!$      integer ii(2)
!!$      integer iaci,iacj,ic1
!!$      integer hit
!!$      real*8 alphp2,ap1,apt1,s
!!$      real*8 pi,pi2,pi22,c_1
!!$      real*8 cp1,cp2,sp1,bo1,dr1,bdexp,bet1,vbb
!!$      real*8 dnex,pcon,efs,frac,cefs,cph,vpb,p3ii,sdg
!!$      real*8 aij1,aij2,pci,alph1,adexp,ba,cona,vpa,hs,rs
!!$      real*8 dvpbdp,ddpbdp,dvdadp,ddpadp,dvpadb,dvpadp,dvpbdb,
!!$     $     ddpadb,ddpbdb
!!$      real*8 dfp,ddfp,dbdp
!!$      real*8 dr(6), ddr(6,6), dx(6)
!!$      real*8 bi
!!$c------------------------------------------------------------------------
!!$      i3=(i1-1)*3
!!$      j3=(j1-1)*3
!!$c
!!$      do j=1,3
!!$         dx(j)=x(i3+j)-x(j3+j)
!!$      enddo 
!!$c
!!$      bi=length(dx)
!!$      iaci=iacod(iac(i1))
!!$      iacj=iacod(iac(j1))
!!$
!!$      ic1 = 0
!!$      hit = 0
!!$      do while(ic1.lt.nbeta.and.hit.eq.0)
!!$         ic1 = ic1 + 1
!!$         if(ic.eq.ibtcd(ic1)) hit = ic1
!!$      enddo
!!$
!!$      if(hit.eq.0) ic1=1
!!$
!!$      alphp2=alpha(2,iaci)+alpha(2,iacj)
!!$      alphp2=0.5*alphp2
!!$
!!$      ap1  = ap
!!$      pi   = dacos(-1.0d0)
!!$      pi2  = pi*0.5
!!$      pi22 = pi*1.5
!!$      apt1 = dabs(ap1)
!!$      cp2  = dcos(ap1)
!!$      s = 1.0
!!$
!!$      if(apt1.ge.pi2.and.apt1.le.pi22) then
!!$         s   = -1.0
!!$         ap1 = ap1-pi
!!$      endif
!!$
!!$      cp1=dcos(ap1)
!!$      sp1=dsin(ap1)
!!$
!!$      bo1=beta(2,ic1)
!!$      dr1=bi-bo1
!!$      bdexp=dexp(-beta(3,ic1)*dr1)
!!$      bet1=beta(1,ic1)*bdexp
!!$      vbb=bet1*(1.+beta(4,ic1)*dr1)
!!$      
!!$      if(nex.gt.0) then
!!$         dnex=1.0
!!$      else
!!$         dnex=0.0
!!$      endif
!!$
!!$c      pcon = 2.0*p3(ih,jh)+dnex*p1(ih,jh)
!!$      efs  = beta(5,ic1)
!!$
!!$c:::  modification of ref.1,where efs=beta(5,ic1)*p3(ih,jh)
!!$
!!$      frac=1.0/iatg(lp)
!!$      cefs=frac/(1.0 + efs)
!!$      cph=(1.0 + efs*cp1)*cefs
!!$      vpb=con*pcon*cph*cp1
!!$
!!$c:::  the coeficients of the contributions from the s**2 terms.
!!$c     for (nex.gt.0) excited state contributions are included
!!$
!!$c      p3ii = 0.5*(p3(ih,ih) + p3(jh,jh))
!!$c      sdg  = 0.125*gamma(4,ic1)*(p3ii + p3(ih,jh)**2)
!!$
!!$c      if(nex.gt.0) then
!!$c         aij1 = ps(ih)     + ps(jh)
!!$c         aij2 = ps(ih+nap) + ps(jh+nap)
!!$c         pci=0.25*(ps(ih)*p3(ih,ih) + ps(jh)*p3(jh,jh))
!!$c     $        + 0.5*aij2 - 0.5*pc(ih,jh)
!!$c         sdg = sdg + 0.5*(gamma(4,ic1)*pci+alphp2*aij1)
!!$c      endif
!!$
!!$c      sdg   = sdg-0.25*gamma(4,ic1)*(z(ih)-p3(ih,ih))*(z(jh)-p3(jh,jh))
!!$c      alph1 = alphp2*p3ii+sdg
!!$c      adexp = dexp(-2.*beta(3,ic1)*dr1)
!!$c      ba    = alph1*adexp
!!$c      cona  = frac*2*con
!!$c      vpa   = cona*cp1*cp1
!!$
!!$c:::  builds the off diagonal h matrix and adds the contributions off th
!!$c     s**2 terms to the r matrix and to the diagonal elements of h
!!$
!!$c      h(ih,jh)=h(ih,jh)+cph*vbb*cp1
!!$c      h(jh,ih)=h(ih,jh)
!!$c      hs=frac*alphp2*adexp*cp1*cp1
!!$c      h(ih,ih)=h(ih,ih)+hs
!!$c      h(jh,jh)=h(jh,jh)+hs
!!$c      rs=.5*hs*gamma(4,ic1)/alphp2
!!$c      r(ih,ih)=r(ih,ih)+rs
!!$c      r(jh,jh)=r(jh,jh)+rs
!!$c      r(ih,jh)=r(ih,jh)-rs
!!$c      r(jh,ih)=r(ih,jh)
!!$
!!$      f=f+vbb*vpb+ba*vpa
!!$
!!$      c_1=con*pcon*cefs
!!$      dvpbdp=2.0*c_1*efs*cp2+c_1*s
!!$      ddpbdp=2.0*c_1*efs
!!$
!!$      dvpadp=2.0*cona*cp2
!!$      ddpadp=2.0*cona
!!$
!!$      dvpadb=alph1*(-2.*beta(3,ic1))*adexp
!!$      ddpadb=alph1*4.*beta(3,ic1)**2*adexp
!!$
!!$      dvpbdb=bet1*beta(4,ic1)-beta(3,ic1)*vbb
!!$      ddpbdb=beta(3,ic1)*(beta(3,ic1)*vbb-2*bet1*beta(4,ic1))
!!$
!!$      df=df+dvpbdp*vbb+dvpadp*ba
!!$      ddf=ddf+ddpbdp*vbb+ddpadp*ba
!!$
!!$      dfp=df
!!$      ddfp=ddf
!!$
!!$      df=dvpadb*vpa+dvpbdb*vpb
!!$      ddf=ddpadb*vpa+ddpbdb*vpb
!!$
!!$      call difbon(dx,dr,ddr,bi,newton)
!!$
!!$      call drd(i3,j3,newton,df,dr,d,ddr)
!!$
!!$      df=dfp
!!$      ddf=ddfp
!!$
!!$c:::  cross derivatives of v with respect to phi and b
!!$
!!$      dbdp=dvpadb*dvpadp+dvpbdb*dvpbdp
!!$      cro(1,4)=dbdp
!!$      do i=1,12
!!$         dtp(1,i)=0.
!!$      enddo
!!$
!!$      do i=1,6
!!$         dtp(1,i+3)=dr(i)
!!$      enddo
!!$
!!$      return
!!$      end 
!!$
!!$
!!$
!!$c----------------------------------------------------------------------

!!$      end
