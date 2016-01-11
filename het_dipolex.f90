!APW it appears there are no forces associated with this subroutine ?
SUBROUTINE het_dipolex(eng)
  
  USE commondata, ONLY:istep,npi
  USE classical, ONLY: pi_iden,pos,dfield
  USE qcfio, ONLY: gen_io
  USE quantum, ONLY: z_atom
  
  IMPLICIT NONE
  
  
  !input variables
  DOUBLE PRECISION :: eng

  !local variables
  INTEGER :: i,ind
  DOUBLE PRECISION :: com,dfkcal,ftemp
  
  !APW kcal/(mol-A)
  !dfkcal=dfield*627.509391d0/0.529177d0
  dfkcal=dfield*627.509491d0
  eng = 0.d0
  com=0.d0
  
  do i=1,npi
     ind=pi_iden(i)
     !APW only x-component for some reason
     com=com+pos(1,ind)
  end do
  com=com/npi
  !ftemp=dfkcal**2*(1.d0-exp(-(dble(istep)/30.d0)**2))
  !ftemp=dfkcal**2
  
  do i=1,npi
     ind=pi_iden(i)
     com=0.d0
     !APW only x-component for some reason
     !eng=eng+ftemp*(pos(1,ind)-com) !*(-1.d0)**i-ftemp*com
     eng=eng+dfkcal*z_atom(ind)*(pos(1,ind)-com)
  end do
END SUBROUTINE het_dipolex
