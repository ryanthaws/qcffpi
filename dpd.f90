!c----------------------------------------------------------------------
!c
!c     pack dp(12) into  d(150)
!c
!c----------------------------------------------------------------------
SUBROUTINE dpd(ind,df,ddf,dc,ddc,tvec,cro,dtp)
  
  USE commondata, ONLY: log_dderiv
  USE classical, ONLY: deriv,dderiv,crosprod
  USE qcfio, ONLY: gen_io

  IMPLICIT NONE
!!$      subroutine dpd(ii,d,dp,df,ne)
  
  !input variables
  INTEGER, DIMENSION(4), INTENT(in) :: ind
  DOUBLE PRECISION, INTENT(in) :: df,ddf
  
  DOUBLE PRECISION, DIMENSION(3,2), INTENT(in) :: dc
  DOUBLE PRECISION, DIMENSION(3,3,3), INTENT(in) :: ddc
  DOUBLE PRECISION, DIMENSION(3,6), INTENT(in) :: tvec
  DOUBLE PRECISION, DIMENSION(4,4), INTENT(in) :: cro
  DOUBLE PRECISION, DIMENSION(3,4,4), INTENT(inout) :: dtp

  !DOUBLE PRECISION, DIMENSION(3,3) :: getmat !APW function

  
  !local variables
  INTEGER :: i,j,k,l,i2,j2,m1,is,in,ak1,ak2,aj1,aj2
  INTEGER :: ir,jr,iatom,jatom,i1,j1,ibeg
  INTEGER, DIMENSION(3,2) :: crd8,crd9
  DOUBLE PRECISION :: ft,ft2,cros
  DOUBLE PRECISION, DIMENSION(3) :: dftemp
  DOUBLE PRECISION, DIMENSION(3,3) :: avec,bvec,ddct4
  DOUBLE PRECISION, DIMENSION(3,3,4,4) :: ddp
  DOUBLE PRECISION, DIMENSION(3,4) :: dp
  
  dp(:,1)=crosprod(tvec(:,1),dc(:,1))
  dp(:,4)=crosprod(tvec(:,1),dc(:,2))
  dp(:,2)=crosprod(tvec(:,2),dc(:,1))
  dftemp=crosprod(tvec(:,3),dc(:,2))
  dp(:,2)=dp(:,2)+dftemp
  dp(:,3)=crosprod(tvec(:,4),dc(:,1))
  dftemp=crosprod(tvec(:,5),dc(:,2))
  dp(:,3)=dp(:,3)+dftemp
  dp(:,4)=-dp(:,4)
  
  do i=1,4
     deriv(:,ind(i))=deriv(:,ind(i))+dp(:,i)*df
     dtp(:,i,4)=dp(:,i)
  end do
  
  if(log_dderiv) then
     
     ddp=0.d0
     
     !c this section places the correct elements of ddc(6,6) in ddp(12,12)
     !c after multiplication by the correct elements of a set of vectors
     !c v(it,l)
     
     do i=1,3
        do j=1,3
           ddct4(i,j)=ddc(j,i,3)
        end do
     end do
     
     
     avec=getmat1(ddc(:,:,1),tvec(:,1))
     bvec=getmat2(avec,tvec(:,1))
     ddp(:,:,1,1)=ddp(:,:,1,1) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,2))
     ddp(:,:,1,2)=ddp(:,:,1,2) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,4))
     ddp(:,:,1,3)=ddp(:,:,1,3) + bvec(:,:)
     
     avec=getmat1(ddc(:,:,3),tvec(:,1))
     bvec=getmat2(avec,tvec(:,3))
     ddp(:,:,1,2)=ddp(:,:,1,2) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,5))
     ddp(:,:,1,3)=ddp(:,:,1,3) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,6))
     ddp(:,:,1,4)=ddp(:,:,1,4) + bvec(:,:)
     
     avec=getmat1(ddc(:,:,1),tvec(:,2))
     bvec=getmat2(avec,tvec(:,2))
     ddp(:,:,2,2)=ddp(:,:,2,2) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,4))
     ddp(:,:,2,3)=ddp(:,:,2,3) + bvec(:,:)

     avec=getmat1(ddc(:,:,3),tvec(:,2))
     bvec=getmat2(avec,tvec(:,3))
     ddp(:,:,2,2)=ddp(:,:,2,2) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,5))
     ddp(:,:,2,3)=ddp(:,:,2,3) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,6))
     ddp(:,:,2,4)=ddp(:,:,2,4) + bvec(:,:)
     
     avec=getmat1(ddct4,tvec(:,3))
     bvec=getmat2(avec,tvec(:,2))
     ddp(:,:,2,2)=ddp(:,:,2,2) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,4))
     ddp(:,:,2,3)=ddp(:,:,2,3) + bvec(:,:)

     avec=getmat1(ddc(:,:,2),tvec(:,3))
     bvec=getmat2(avec,tvec(:,3))
     ddp(:,:,2,2)=ddp(:,:,2,2) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,5))
     ddp(:,:,2,3)=ddp(:,:,2,3) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,6))
     ddp(:,:,2,4)=ddp(:,:,2,4) + bvec(:,:)
     
     avec=getmat1(ddc(:,:,1),tvec(:,4))
     bvec=getmat2(avec,tvec(:,4))
     ddp(:,:,3,3)=ddp(:,:,3,3) + bvec(:,:)
     
     avec=getmat1(ddc(:,:,3),tvec(:,4))
     bvec=getmat2(avec,tvec(:,5))
     ddp(:,:,3,3)=ddp(:,:,3,3) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,6))
     ddp(:,:,3,4)=ddp(:,:,3,4) + bvec(:,:)
     
     avec=getmat1(ddct4,tvec(:,5))
     bvec=getmat2(avec,tvec(:,4))
     ddp(:,:,3,3)=ddp(:,:,3,3) + bvec(:,:)

     avec=getmat1(ddc(:,:,2),tvec(:,5))
     bvec=getmat2(avec,tvec(:,5))
     ddp(:,:,3,3)=ddp(:,:,3,3) + bvec(:,:)
     bvec=getmat2(avec,tvec(:,6))
     ddp(:,:,3,4)=ddp(:,:,3,4) + bvec(:,:)
     
     avec=getmat1(ddc(:,:,2),tvec(:,6))
     bvec=getmat2(avec,tvec(:,6))
     ddp(:,:,4,4)=ddp(:,:,4,4) + bvec(:,:)

     !place in ddp(12,12) the contibution due to the 1st derivatives
     i2=3
     ft=1.0d0
     do i=1,2
        do j=i+1,3
           ddp(i,j,1,2)=ddp(i,j,1,2)-dc(i2,1)*ft
           ddp(j,i,1,2)=ddp(j,i,1,2)+dc(i2,1)*ft
           ddp(i,j,1,3)=ddp(i,j,1,3)+dc(i2,1)*ft
           ddp(j,i,1,3)=ddp(j,i,1,3)-dc(i2,1)*ft
           ddp(i,j,2,3)=ddp(i,j,2,3)-dc(i2,1)*ft
           ddp(j,i,2,3)=ddp(j,i,2,3)+dc(i2,1)*ft
           
           ddp(i,j,3,4)=ddp(i,j,3,4)+dc(i2,2)*ft
           ddp(j,i,3,4)=ddp(j,i,3,4)-dc(i2,2)*ft
           ddp(i,j,2,4)=ddp(i,j,2,4)-dc(i2,2)*ft
           ddp(j,i,2,4)=ddp(j,i,2,4)+dc(i2,2)*ft
           ddp(i,j,2,3)=ddp(i,j,2,3)+dc(i2,2)*ft
           ddp(j,i,2,3)=ddp(j,i,2,3)-dc(i2,2)*ft
        
           i2=i2-1
           ft=-ft
        end do
     end do
     
     do k=1,4
        do l=k,4
           do i=1,3
              do j=1,3
                 ddp(j,i,l,k)=ddp(i,j,k,l)
              end do
           end do
        end do
     end do
     
     !matpac
     !APW the routine below has been rewritten for a proper 
     !4-dimensional matrix...screw this 1-D crap, you know what I mean?
     !c-------------------------------------------------------------------
     !c
     !c     this routine packs the upper half of the matrix of energy second
     !c     derivatives w.r.t. cartesian in a vector called dderiv .
     !c     diagonal terms(diagonal derivatives with respect to internal
     !c     coordinates , s ) are constructed using df (df=df/ds) , ddf
     !c     (ddf=d2f/ds*ds) , dp (dp=ds/dx) and ddp (ddp=d2s/(dx(i)*dx(j)) .
     !c     cross terms (non diagonal derivatives with respect to internal
     !c     coordinates ) are constructed using cro(i,j) (where cro(i,j) is
     !c     d2f/(ds(i)*ds(j)) ) and dtp (where dtp(ir,l)=ds(ir)/dx(l))
     !c     an element (i,j) , of the second deriv matrix  is the m-th
     !c     element in the dd vector,where
     !c     m=(nat3+nat3-(i-2))*(i-1)/2+j-i+1=nat3*(i-1)-i*(i-1)/2+j
     !c     ne=4 for phi matrix , 3 for theta and 2 for bond matrices
     !c
     !c     intp=1 for diagonal terms ,intp=2,3 and 4 for cross terms
     !c     of the corresponding order (e.g. intp=2 for f(bond,theta) ).
     !c
     !c-------------------------------------------------------------------
     
     do iatom=1,4
        i1=ind(iatom)
        do i2=1,3
           !k=(iatom-1)*3+i2
           ibeg=i2
           do jatom=iatom,4
              j1=ind(jatom)
              if(jatom.gt.iatom) ibeg=1
              do j2=ibeg,3
                 !l=(jatom-1)*3+j2
                 dderiv(j2,i2,j1,i1) = dderiv(j2,i2,j1,i1)&
                      +ddp(j2,i2,jatom,iatom)*df+dp(i2,iatom)*dp(j2,jatom)*ddf
                 dderiv(i2,j2,i1,j1) = dderiv(j2,i2,j1,i1)
                 do ir=1,3
                    do jr=ir+1,4
                       cros=cro(ir,jr)*(dtp(i2,iatom,ir)*dtp(j2,jatom,jr)+dtp(j2,jatom,ir)*dtp(i2,iatom,jr))
                       dderiv(j2,i2,j1,i1)=dderiv(j2,i2,j1,i1)+cros
                       dderiv(i2,j2,i1,j1)=dderiv(j2,i2,j1,i1)
                    end do
                 end do

              end do
           end do
        end do
    end do
    
 end if
 
  CONTAINS
  
  FUNCTION getmat1(mat1,v1)
    
    USE classical, ONLY: crosprod
    
  IMPLICIT NONE
  
  DOUBLE PRECISION, DIMENSION(3,3) :: getmat1
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(in) :: mat1
  DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: v1
  
  INTEGER :: i
  DOUBLE PRECISION, DIMENSION(3,3) :: tmat
  
  do i=1,3
     getmat1(:,i)=crosprod(mat1(:,i),v1)
  end do
  
end  FUNCTION getmat1
  
FUNCTION getmat2(mat1,v1)
  
  USE classical, ONLY: crosprod
  
  IMPLICIT NONE
  
  DOUBLE PRECISION, DIMENSION(3,3) :: getmat2
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(in) :: mat1
  DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: v1
  
  INTEGER :: i
  
  do i=1,3
     getmat2(i,:)=crosprod(mat1(i,:),v1)
  end do
  
end FUNCTION getmat2

end SUBROUTINE dpd
