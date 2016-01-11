!c------------------------------------------------------------
!c                                                           c
!c  subroutine computing second step in velocity verlet alg. c
!c                                                           c
!c  2) v(t+dt)=v(t)+0.5*(f(t+dt))/m)dt                       c
!c                                                           c
!c-----------------------------------------------------------c
!c                                                           c
!c   Modified from original subroutine of JL                 c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm/utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c

SUBROUTINE verlet_s2
  !SUBROUTINE VELVERLEI2_PPV
  
  USE commondata, ONLY: nnuc,deltat,navg
  USE classical, ONLY: vel,fra,amass
  
  IMPLICIT NONE
  
  INTEGER :: i
  
  do i=1,nnuc
     vel(:,i) = vel(:,i) + 100.0d0*0.5d0*deltat*fra(:,i)/amass(i)
     !print *,'shit2e',navg,i,fra(:,i)*0.5d0*100.0d0*deltat/amass(i)
     !print *,'shit2f',navg,i,vel(:,i)
  end do
      
  
END SUBROUTINE verlet_s2


