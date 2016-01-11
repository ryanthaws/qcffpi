!c------------------------------------------------------------c
!c    Compute the cosine of dihedral Angle and its gradient   c
!c------------------------------------------------------------c
!c                                                            c
!c   (OLD NOTATION)                                           c
!c    cos(teta)=a*c/|a||c|                                    c
!c                                                            c
!c    a   = v1-(v1*b12/|b12|)b12/|b12|                        c
!c    c   = v2-(v2*b12/|b12|)b12/|b12|                        c
!c                                                            c
!c    v1  = j - i                                             c
!c    v2  = l - k                                             c
!c    b12 = ib2 - ib1                                         c
!c                                                            c
!c      i          k                                          c
!c       \        /                                           c
!c        ib1--ib2                                            c
!c       /        \                                           c
!c      j          l                                          c
!c                                                            c
!c------------------------------------------------------------c
!c   (PRESENT NOTATION)                                       c
!c    cos(teta)=a1*a2/|a1||a2|                                c
!c                                                            c
!c    a1   = bnd1-(bnd1*bnd3/|bnd3|)bnd3/|bnd32|              c
!c    a2   = bnd2-(bnd2*bnd3/|bnd3|)bnd3/|bnd32|              c
!c                                                            c
!c    bnd1  = j - i                                           c
!c    bnd2  = l - k                                           c
!c    bnd3 = ib2 - ib1                                        c
!c                                                            c
!c      i          k                                          c
!c       \        /                                           c
!c        ib1--ib2                                            c
!c       /        \                                           c
!c      j          l                                          c
!c                                                            c
!-------------------------------------------------------------!
!APW this subroutine was origanally carried out in units of   !
! angstrom/150, why I don't know but I plan to correct it     !
!-------------------------------------------------------------!


SUBROUTINE died_ppv(ind,ftemp,fcos)
  !AY
  !  ind(6)     : six atoms defining a torsion angle 
  !  ftemp(3,6) : gradient (x,y,z directions) of six atoms  (= dcos(phi)/dx)
  !  fcos       : cosine of torsional angle (= cos(phi))

  
  USE commondata, ONLY:DP,a0
  USE classical, ONLY: pos
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE 

  !input variables
  INTEGER, DIMENSION(6), INTENT(in) :: ind
  REAL(DP), INTENT(out) :: fcos
  REAL(DP), DIMENSION(3,6), INTENT(out) :: ftemp
  
  !local variables
  REAL(DP) :: a1a2,a1a,a1a1i,a1a2i,a2a2,a2a,a2a1i,a2a2i,dprod,&
       b1b2,fac1,fac2,bb2,bb,bb1i,bb2i,b1b3,b2b3
  REAL(DP), DIMENSION(3) :: a1,a2,bnd1,bnd2,bnd3
  
  !APW apparently excited dynamics are done in atomic units of length
  !c--- unit length for exicted dynamic -------------------------
  !if(exd) then
  !   do mm=1,nat_sol
  !      xna(mm)=xna(mm)*unitl/0.529177
  !      yna(mm)=yna(mm)*unitl/0.529177
  !      zna(mm)=zna(mm)*unitl/0.529177
  !   enddo
  !endif
      
  !c--- set force zero ------------------------------------------
  
  fcos=0.d0
  ftemp=0.d0
  
  !fxib1=0.d0
  !fyib1=0.d0
  !fzib1=0.d0
  
  !fxib2=0.d0
  !fyib2=0.d0
  !fzib2=0.d0
  
  !c--- Bonds ---------------------------------------------------
  
  !APW change the length unit conversion (150) soon
  !v1xyz
  bnd1=(pos(:,ind(4))-pos(:,ind(3)))/a0
  !v2xyz
  bnd2=(pos(:,ind(6))-pos(:,ind(5)))/a0
  !bxyz
  bnd3=(pos(:,ind(2))-pos(:,ind(1)))/a0
  
  !v1v2
  b1b2=bnd1(1)*bnd2(1)+bnd1(2)*bnd2(2)+bnd1(3)*bnd2(3)
  
  bb2=bnd3(1)**2+bnd3(2)**2+bnd3(3)**2
  bb=sqrt(bb2)
  bb1i=1/bb
  bb2i=1/bb2
  
  !v1b
  b1b3=bnd1(1)*bnd3(1)+bnd1(2)*bnd3(2)+bnd1(3)*bnd3(3)
  !v2b
  b2b3=bnd2(1)*bnd3(1)+bnd2(2)*bnd3(2)+bnd2(3)*bnd3(3)
  
  !axyz
  a1=bnd1-bb2i*bnd3*b1b3
  
  a1a2=a1(1)**2+a1(2)**2+a1(3)**2
  a1a=dsqrt(a1a2)
  a1a1i=1/a1a
  a1a2i=1/a1a2
  
  !cxyz
  a2=bnd2-bb2i*bnd3*b2b3
  a2a2=a2(1)**2+a2(2)**2+a2(3)**2
  a2a=dsqrt(a2a2)
  a2a1i=1/a2a
  a2a2i=1/a2a2
  
  !c--- Cos(teta) --------------------------------------------
  
  dprod=(a1(1)*a2(1)+a1(2)*a2(2)+a1(3)*a2(3))
  fcos = dprod / (a1a*a2a)  
  !N=(ax*cx)+(ay*cy)+(az*cz)
  !D=aa*cc
  !D2=D*D
  !invD=1/D
  !invD2=1/D2
  
  !cos=N/D
  
  !c--- Forces on the 6 atoms --------------------------------
  
  !c--- Atoms i and j  -----------------------------------------------
  
  !dNx=(-cx)
  !dNy=(-cy)
  !dNz=(-cz)
  
  !dDx=(ax)
  !dDy=(ay)
  !dDz=(az)
  
  fac1=a1a2i*dprod
  
  ftemp(:,3)=(fac1*a1-a2)/(a1a*a2a)

  !fxi=invD*(dNx+fac1*dDx)
  !fyi=invD*(dNy+fac1*dDy)
  !fzi=invD*(dNz+fac1*dDz)
  ftemp(:,4)=-ftemp(:,3)
  !fxj=-fxi
  !fyj=-fyi
  !fzj=-fzi
      
  !c--- Atoms k and l  ------------------------------------------------

  !dNx=(-ax)
  !dNy=(-ay)
  !dNz=(-az)
  
  !dDx=(cx)
  !dDy=(cy)
  !dDz=(cz)
  
  fac2=a2a2i*dprod
  
  ftemp(:,5)=(fac2*a2-a1)/(a1a*a2a)
  !fxk=invD*(dNx+fac2*dDx)
  !fyk=invD*(dNy+fac2*dDy)
  !fzk=invD*(dNz+fac2*dDz)
  ftemp(:,6)=-ftemp(:,5)
  !fxl=-fxk
  !fyl=-fyk
  !fzl=-fzk
  
  !c--- Atoms ib1 and ib2 --------------------------------------------
  
  fac1=bb2i/(a1a*a2a)
  !facb1=invbb2*invD
  fac2=-bb2i*dprod/(a1a*a2a)**2
  !facb2=-invD2*N*invbb2
  
  ftemp(:,1)=fac1*(a2*b1b3+a1*b2b3)
  ftemp(:,1)=ftemp(:,1)+fac2*(a2a*a1a1i*a1*b1b3)
  ftemp(:,1)=ftemp(:,1)+fac2*(a1a*a2a1i*a2*b2b3)
  !fxib1=facb1*(cx*v1b+ax*v2b)
  !fxib1=fxib1+facb2*(cc*invaa*(ax*v1b))
  !fxib1=fxib1+facb2*(aa*invcc*(cx*v2b))

  ftemp(:,2)=-ftemp(:,1)
  
  !c--- Save in tmp array ---------------------------------------------
  
!!$  tmp1x(ir)=fxi
!!$  tmp1y(ir)=fyi
!!$  tmp1z(ir)=fzi
!!$  
!!$  tmp2x(ir)=fxj
!!$  tmp2y(ir)=fyj
!!$  tmp2z(ir)=fzj
!!$  
!!$  tmp3x(ir)=fxk
!!$  tmp3y(ir)=fyk
!!$  tmp3z(ir)=fzk
!!$  
!!$  tmp4x(ir)=fxl
!!$  tmp4y(ir)=fyl
!!$  tmp4z(ir)=fzl
!!$  
!!$  tmpb1x(ir)=fxib1
!!$  tmpb1y(ir)=fyib1
!!$  tmpb1z(ir)=fzib1
!!$  
!!$  tmpb2x(ir)=fxib2
!!$  tmpb2y(ir)=fyib2
!!$  tmpb2z(ir)=fzib2
  

  !if(exd) then
  !   do mm=1,nat_sol
  !      xna(mm)=xna(mm)/(unitl/0.529177)
  !      yna(mm)=yna(mm)/(unitl/0.529177)
  !      zna(mm)=zna(mm)/(unitl/0.529177)
  !   enddo
  !endif

END SUBROUTINE died_ppv
