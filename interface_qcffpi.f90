!c------------------------------------------------------------
!c                                                           c
!c  subroutine interfacing QCFFPI and PPV                    c
!c                                                           c
!c  unpack derivative for PPV calculation                    c
!c-----------------------------------------------------------c
!c                                                           c
!c   Modified from original subroutine of JL                 c
!c   Author Fabio Sterpone 2004                              c
!c   sterpone@cm/utexas.edu                                  c
!c                                                           c
!c   NoCopyright                                             c
!c                                                           c
!c-----------------------------------------------------------c

SUBROUTINE interface_qcffpi(in,out)
  !subroutine interface_qcffpi(in,out,xna,yna,zna,x_ppv,fqcffx,fqcffy
  !$     ,fqcffz,d_ppv,unitl_ppv)
  USE classical, ONLY: deriv,pos
  
  IMPLICIT none
  
  !input variables
  LOGICAL :: in,out
  
  !local variables
  INTEGER :: i
  !integer i,i1
  !real*8 fqcffx(*),fqcffy(*),fqcffz(*)
  !real*8 d_ppv(*),x_ppv(*)
  
  !real*8 xna(*),yna(*),zna(*),unitl_ppv
  !logical in,out,print_pdb
  
  
  
  
  !c--- Executive Stantements ----------------------------------
  
  
  !if(out) then
  !   do i=1,nnuc
  !      !c--- In qcffsol we take gradient of the U not  -gradU -------
  !      fqcffx(:,i)=-deriv(:,i)
  !   enddo
  !endif
      
  !if(in) then 
  !   do i=1,nat_sol
  !      i1=(i-1)*3
  !      x_ppv(i1+1)=xna(i)*unitl_ppv  !puts positions into angstrom
  !      x_ppv(i1+2)=yna(i)*unitl_ppv
  !      x_ppv(i1+3)=zna(i)*unitl_ppv
  !   enddo
  !endif
 

  return
END SUBROUTINE interface_qcffpi
      
