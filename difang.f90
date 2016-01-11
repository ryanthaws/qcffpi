  SUBROUTINE difang(dx,dy,ra,rb,cthet,dvcos,ddvcos,dsign,iarc)
    !subroutine difang  ( x,bi,bj,ct,dc,ddc,dsign,newton,iarc)
    !c----------------------------------------------------------------------
    !c     subroutine difang computes:
    !c     the derivatives of cos(theta),if iarc=.false.
    !c     the derivatives of theta,     if iarc=.true.
    !c     w.r.t. the components of the 2 vectors forming the theta angle
    !c----------------------------------------------------------------------    
    
    USE commondata, ONLY: log_dderiv
    USE qcfio, ONLY:gen_io
    
    IMPLICIT NONE
    
    !input variables
    LOGICAL, INTENT(in) :: iarc
    DOUBLE PRECISION, INTENT(in) :: ra,rb,cthet
    DOUBLE PRECISION, DIMENSION(3), INTENT(in) :: dx,dy
    DOUBLE PRECISION, DIMENSION(3,2), INTENT(out) :: dvcos
    DOUBLE PRECISION, DIMENSION(3,3,3), INTENT(out) :: ddvcos
    
    !local variables
    INTEGER :: i,j,k
    DOUBLE PRECISION :: ra1i,rb1i,&
         raorb,rbora,rb2i,ra2i,ra3i,rb3i,ra3ib1i,rb3ia1i,ra2ib2i
    DOUBLE PRECISION :: sthet,sthet2,dsign,ftemp
    
    dvcos=0.d0
    ddvcos=0.d0
    
        !c:::  ist derivatives
    sthet2=1.0-cthet**2
    if(dabs(sthet2) .lt. 1.0d-10) then
       !sthet2=1.0d-10   !APW old version why??
       !write(gen_io,*) '***sin^2(theta) is small...',sthet2,'***'
    end if
    sthet=dsqrt(sthet2)
    sthet=sthet*dsign !APW changed from if(dsign.lt.0) st=-st
    ra1i=1.d0/ra
    rb1i=1.d0/rb
    ra2i=ra1i**2
    rb2i=rb1i**2
    raorb=ra*rb1i
    rbora=rb*ra1i
    dvcos(:,1)=ra2i*(raorb*dy-cthet*dx)
    dvcos(:,2)=rb2i*(rbora*dx-cthet*dy)
    !APW old
    !do i=1,3
    !   dvcos(1,i)=ra2i*(raorb*dy(i)-cthet*dx(i))
    !   dvcos(2,i)=rb2i*(rbora*dx(i)-cthet*dy(i))
    !end do
    if(iarc) then
       dvcos=-dvcos/sthet
    end if
    
    if(log_dderiv) then
       do i=1,3
          ddvcos(i,i,1)=-ra2i*(2.d0*dx(i)*dvcos(i,1)&
               +cthet*(1.d0-ra2i*dx(i)**2))
          ddvcos(i,i,2)=-rb2i*(2.d0*dy(i)*dvcos(i,2)&
               +cthet*(1.d0-rb2i*dy(i)**2))
       end do
       ra3i=ra2i*ra1i
       rb3i=rb2i*rb1i
       ra3ib1i=ra3i*rb1i
       rb3ia1i=rb3i*ra1i
       ra2ib2i=ra2i*rb2i
       do i=1,3
          ddvcos(i,i,3)=-ra1i*(rb3i*dy(i)**2-rb1i+dx(i)*dvcos(i,2)*ra1i)
       end do
       
       ftemp=-dy(1)*dy(2)*rb3ia1i-dx(1)*dx(2)*ra3ib1i
       ddvcos(1,2,3)=ftemp+dx(1)*dy(2)*cthet*ra2ib2i
       ddvcos(2,1,3)=ftemp+dx(2)*dy(1)*cthet*ra2ib2i
       
       ftemp=-dy(2)*dy(3)*rb3ia1i-dx(2)*dx(3)*ra3ib1i
       ddvcos(2,3,3)=ftemp+dx(2)*dy(3)*cthet*ra2ib2i
       ddvcos(3,2,3)=ftemp+dx(3)*dy(2)*cthet*ra2ib2i
       
       ftemp=-dy(3)*dy(1)*rb3ia1i-dx(3)*dx(1)*ra3ib1i
       ddvcos(3,1,3)=ftemp+dx(3)*dy(1)*cthet*ra2ib2i
       ddvcos(1,3,3)=ftemp+dx(1)*dy(3)*cthet*ra2ib2i
       
       !APW d^2/dxidyi
       ddvcos(1,2,1)=-ra2i*(dx(1)*dvcos(2,1)+dx(2)*dvcos(1,1)&
            -cthet*dx(1)*dx(2)*ra2i)
       ddvcos(1,2,2)=-rb2i*(dy(1)*dvcos(2,2)+dy(2)*dvcos(1,2)&
            -cthet*dy(1)*dy(2)*rb2i)
       
       ddvcos(1,3,1)=-ra2i*(dx(1)*dvcos(3,1)+dx(3)*dvcos(1,1)&
            -cthet*dx(1)*dx(3)*ra2i)
       ddvcos(1,3,2)=-rb2i*(dy(1)*dvcos(3,2)+dy(3)*dvcos(1,2)&
            -cthet*dy(1)*dy(3)*rb2i)
           
       ddvcos(2,3,1)=-ra2i*(dx(2)*dvcos(3,1)+dx(3)*dvcos(2,1)&
            -cthet*dx(2)*dx(3)*ra2i)
       ddvcos(2,3,2)=-rb2i*(dy(2)*dvcos(3,2)+dy(3)*dvcos(2,2)&
            -cthet*dy(2)*dy(3)*rb2i)
       
       do i=1,3
          do j=1,3
             ddvcos(j,i,1)=ddvcos(i,j,1)
             ddvcos(j,i,2)=ddvcos(i,j,2)
          end do
       end do
       
       if(iarc) then
          do i=1,3
             do j=1,3
                ddvcos(i,j,1)=-(ddvcos(i,j,1)+dvcos(i,1)*dvcos(j,1)&
                     *cthet/sthet2)/sthet
                ddvcos(i,j,2)=-(ddvcos(i,j,2)+dvcos(i,2)*dvcos(j,2)&
                     *cthet/sthet2)/sthet
                ddvcos(i,j,3)=-(ddvcos(i,j,3)+dvcos(i,1)*dvcos(j,2)&
                     *cthet/sthet2)/sthet
             end do
          end do
       endif
    end if
    
  end SUBROUTINE difang
