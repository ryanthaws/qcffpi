  SUBROUTINE drd(ind,df,ddf,dr,ddr)
    
    USE commondata, ONLY: log_dderiv
    USE classical, ONLY:deriv,dderiv
    USE qcfio, ONLY:gen_io
    
    IMPLICIT NONE
    
    !input variables
    INTEGER, DIMENSION(2), INTENT(IN) :: ind
    DOUBLE PRECISION, INTENT(IN) :: df,ddf
    DOUBLE PRECISION, DIMENSION(3,2), INTENT(IN) :: dr
    DOUBLE PRECISION, DIMENSION(3,3,3), INTENT(IN) :: ddr
    
    !local variables
    INTEGER :: i,j,k,l,iatom,jatom,i1,i2,j1,j2,ibeg
    DOUBLE PRECISION, DIMENSION(3,3,2,2) :: ddp
    
    ddp=0.d0
    
    !:::  pack 1st derivatives vector deriv
    do j=1,3
       deriv(j,ind(1))=deriv(j,ind(1))+df*dr(j,1)
       deriv(j,ind(2))=deriv(j,ind(2))-df*dr(j,1)
    end do

    if(log_dderiv) then    
       do i=1,3
          do j=1,3
             ddp(j,i,1,1)=ddr(i,j,1)
             ddp(j,i,2,2)=ddr(i,j,2)
             ddp(j,i,2,1)=ddr(i,j,3)
             ddp(j,i,1,2)=ddr(j,i,3)
          end do
       end do
       
       !c:::  pack 2nd derivatives matrix dd in matpac
       
       !call matpac(ind,2,dr,df,ddf)
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
              
       do iatom=1,2
          i1=ind(iatom)
          do i2=1,3
             !k=(iatom-1)*3+i2
             ibeg=i2
             do jatom=iatom,2
                j1=ind(jatom)
                if(jatom.gt.iatom) ibeg=1
                !do j2=1,3
                do j2=ibeg,3
                   !l=(jatom-1)*3+i4
                   dderiv(j2,i2,j1,i1) = dderiv(j2,i2,j1,i1)&
                        +ddp(j2,i2,jatom,iatom)*df+dr(i2,iatom)*dr(j2,jatom)*ddf
                   dderiv(i2,j2,i1,j1) = dderiv(j2,i2,j1,i1)
                end do
             end do
          end do
       end do
    end if
    
  END SUBROUTINE drd
