SUBROUTINE sort2(n,arr,brr,crr)
  
  IMPLICIT NONE
  
  !input variables
  INTEGER, INTENT(in) :: n
  DOUBLE PRECISION, DIMENSION(n), INTENT(inout) :: arr
  INTEGER, DIMENSION(n), INTENT(inout) :: brr,crr
  
  !local variables
  LOGICAL :: ltemp
  INTEGER :: jmax,i,j,itemp
  DOUBLE PRECISION :: ftemp
  

  jmax=n-1
  
  do i=1,n-1
     ltemp=.true.
     do j=1,jmax
        if(arr(j).gt.arr(j+1)) then
           ltemp=.false.
           ftemp=arr(j)
           arr(j)=arr(j+1)
           arr(j+1)=ftemp
           itemp=brr(j)
           brr(j)=brr(j+1)
           brr(j+1)=itemp
           itemp=crr(j)
           crr(j)=crr(j+1)
           crr(j+1)=itemp
        end if
     end do
     if(ltemp) return
     jmax=jmax-1
  end do
        
end SUBROUTINE sort2
