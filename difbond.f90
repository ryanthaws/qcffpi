!!$c----------------------------------------------------------------------
!!$c     this subroutine returns the 1st and 2nd derivatives of r(i,j)
!!$c----------------------------------------------------------------------
 SUBROUTINE difbond(x,dr,ddr,r)
   !subroutine difbon (x,dr,ddr,r,newton)   
   
   USE commondata, ONLY: log_dderiv
   USE qcfio, ONLY: gen_io
   
   IMPLICIT NONE
   
   !c:::  input vars
   DOUBLE PRECISION, INTENT(IN) :: r
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: x
   DOUBLE PRECISION, DIMENSION(3,2), INTENT(OUT) :: dr
   DOUBLE PRECISION, DIMENSION(3,3,3), INTENT(OUT) :: ddr 
   
   !c:::  local vars
   INTEGER :: i,j,k
   DOUBLE PRECISION :: r1i,a
   
   dr=0.d0
   ddr=0.d0
   
   !c:::  1st derivatives
   r1i=1.0/r
   dr(:,1)=x*r1i
   dr(:,2)=-x*r1i
   !APW old
   !do k=1,3
   !   dr(k,1)=x(k)*r1i
   !   dr(k,2)=-dr(1,k)
   !end do
   
   if(log_dderiv) then
      do k=1,3
         a=(1.0-dr(k,1)*dr(k,1))*r1i
         ddr(k,k,1)=a
         ddr(k,k,2)=a
         ddr(k,k,3)=-a
      end  do
      
      do i=1,2
         do j=i+1,3
            a=-dr(i,1)*dr(j,1)*r1i
            ddr(i,j,1)=a
            ddr(i,j,2)=a
            ddr(i,j,3)=-a
         end  do
      end do
      !ddr(3,5)=ddr(2,6)
      !ddr(2,4)=ddr(1,5)
      !ddr(3,4)=ddr(1,6)
      
      do i=1,3
         do j=1,3
            ddr(j,i,1)=ddr(i,j,1)
            ddr(j,i,2)=ddr(i,j,2)
            ddr(j,i,3)=ddr(i,j,3)
         end do
      end do
   end if
   
 END SUBROUTINE difbond
