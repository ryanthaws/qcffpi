  !* =====================================================================
  !* []  DECK FINDCDOT  -- DERIVATIVE OF TULLY EXPANSION COEFFICIENTS
  !* =====================================================================
  !c----------------------------------------------------------------------
  !c For MEH-PPV we don't need to add the normal mode term
  !c----------------------------------------------------------------------
SUBROUTINE findcdot(t,cdot,ctul)
  
  USE commondata, ONLY: DP,hbar,navg
  USE quantum, ONLY: fdv,fdv_l,ciener,ciener_l,ntully,deltae

  IMPLICIT NONE
  
  !input variables
  COMPLEX*16, DIMENSION(0:ntully), INTENT(out) :: cdot
  COMPLEX*16, DIMENSION(0:ntully), INTENT(in) :: ctul
  REAL(DP), INTENT(in) :: t
  
  !local variables
  INTEGER :: i,j
  REAL(DP) :: efactor,fdotv
  
  !*  +---------------------------------------------------------------+
  !*  |  converts energy in au (i.e. hartress) to program units       |
  !*  +---------------------------------------------------------------+
  efactor=2625.5
  !enerfac = 4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !enerfac = enerfac/fkjmol

  do j = 0,ntully
     deltae=0.0d0
     if(j.gt.0) then
        deltae=(ciener(j)-ciener_l(j))*t+ciener_l(j)
     end if
     cdot(j)=-dcmplx(0.0d0,1.0d0)*efactor*deltae*ctul(j)/hbar
     
     do i = 0,ntully
        if (i.ne.j) then
           fdotv=(fdv(i,j)-fdv_l(i,j))*t+fdv_l(i,j)
           cdot(j) = cdot(j) - fdotv*ctul(i)
        end if
     end do
  end do

END SUBROUTINE findcdot


