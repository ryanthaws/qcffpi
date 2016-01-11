!* =====================================================================
!* []  DECK FDOTVOUT -- RETURN THE LINEAR INTERPOLATION OF         
!*                      NONADIABATIC COUPLING VECTOR            
!* =====================================================================

SUBROUTINE fdotvout(n,i,j)
  !SUBROUTINE FDOTVOUT_PPV(T,J,I)    
  
  USE commondata, ONLY: DP,nnuc
  USE classical, ONLY: vel
  USE quantum, ONLY: ciener,fexa
  USE qcfio, ONLY: gen_io,fdotv_io
  
  IMPLICIT NONE
  
  !input variables
  INTEGER, INTENT(IN) :: n,i,j
  
  !local variables
  INTEGER :: iatom
  REAL(DP), DIMENSION(3) :: fa,va
  REAL(DP) :: efactor,deltae,dot,sdot,bdot,ftemp

  !*  +---------------------------------------------------------------+
  !*  |  converts energy in au (i.e. hartress) to program units       |
  !*  +---------------------------------------------------------------+
  !enerfac = 4.3598d-18*1.0d-3/dfloat(molcs)*avogad
  !enerfac = enerfac/fkjmol
  
  efactor = 2625.5

  if (i.eq.j) then 
     write(gen_io,*) 'DABOO You ask the impossible master'
     write(gen_io,*) ' I quit'
     stop
  end if

  if (i.eq.0) then
     deltae = -ciener(j)
  else
     if (j.eq.0) then
        deltae = ciener(i)
     else
        deltae = ciener(i)-ciener(j)
     end if
  end if

  !deltae = deltae*efactor

  dot = 0.0d0
  sdot = 0.d0
  bdot = 0.0d0
  
  !*  +---------------------------------------------------------------+
  !*  |  calculate F \cdot v for solvent atoms                        |
  !*  +---------------------------------------------------------------+
  !c     IF (solv) THEN
  !c     do iatom = 1,natoms
  !c         ifr(iatom) = iatom
  !c         ifr2(iatom) = iatom
  !c
  !c         fsx = fexx(iatom,j,i)
  !c         fsy = fexy(iatom,j,i)
  !c         fsz = fexz(iatom,j,i)
  !c
  !c         vsx = vx(iatom,np)
  !c         vsy = vy(iatom,np)
  !c         vsz = vz(iatom,np)
  !c
  !c         fr(iatom) = dabs((fsx*vsx + fsy*vsy + fsz*vsz)/deltae)
  !c         dot = dot + fsx*vsx + fsy*vsy + fsz*vsz
  !c     enddo
  !c     ENDIF
  !c     sdot = dot
  
  !*  +---------------------------------------------------------------+
  !*  |  now for PPV atoms                                            |
  !*  +---------------------------------------------------------------+
  
  do iatom = 1,nnuc
     !ifra(iatom) = iatom
     !c        ifra2(iatom) = iatom
     
     fa(:) = fexa(:,iatom,j,i)
     !fax = fexxa(iatom,j,i)
     !fay = fexya(iatom,j,i)
     !faz = fexza(iatom,j,i)
      
     va(:) = vel(:,iatom)
     !vax = vxa(iatom)
     !vay = vya(iatom)
     !vaz = vza(iatom)
     ftemp = fa(1)*va(1) + fa(2)*va(2)+ fa(3)*va(3)
     dot = dot + ftemp
     !bdot = bdot + ftemp
     !fra(iatom) = dabs((fax*vax + fay*vay + faz*vaz)/deltaE)
  end do
  
  ftemp=dot
  dot = dot /(deltae*efactor)                      !total
  !sdot = sdot /deltaE                    !solvent only
  !bdot = bdot /deltaE                    !solute only
  
  !APW deltae in Hartree
  !write(fdotv_io,1112)n,i,j,dot,deltae,fexa(1,1,i,j),&
  !     fexa(2,1,i,j),fexa(3,1,i,j),va(1),va(2),va(3)
  !write(fdotv_io,*)n,i,j,dot,deltae,fexa(1,1,i,j),&
  !     fexa(2,1,i,j),fexa(3,1,i,j),va(1),va(2),va(3)
  write(fdotv_io,'(3i5,3f25.15)') n,i,j,dot,ftemp,deltae*efactor
1112 format(f16.8,2i8,10f16.8)
!c     write(*,112) infnty,i,j,dot,deltaE/enerfac
!c112  format('this is coupling',3i10,4e16.4)

!*     call hpcisort(natoms,fr,ifr,ifr)
!*     call hpcisort(iatomsa,fra,ifra,ifra)


!c----------------------------------------------------------------------
!c Where the variables in sort2 are used?
!c sort2 sorts each atom's NA coupling strength from lowest
!c to highest, but for what purpose?
!c----------------------------------------------------------------------

!c     call sort2(natoms,fr,ifr,ifr)    !IF (SOLV) THEN use this line
  !APW is this necessary?
  !call sort2(nnuc,fra,ifra,ifra)


  !c     if (solv) then
  !c     do k = natoms,natoms-10,-1
  !c        if ((ifr(k).lt.0).or.(ifr(k).gt.natoms)) then
  !c           write(*,*) 'PROBLEM IN SOLVENT '
  !c           write(*,*) 'found one a problem',k,ifr(k),fr(k)
  !c        else
  !c           write(32,1111) t,fr(k),ifr(k),fexx(ifr(k),j,i),
  !c    x           fexy(ifr(k),j,i),fexz(ifr(k),j,i)
  !c        endif
  !c     enddo
  !c     endif
      
  
  !c--- vedere che significa iatoms-10
  !c--- anche qui c''e una parte di solvente che puo essere skippata
  !c--- con condizione!!!!!!vedi sopra!
  
  !cmbh - iatomsa is the total number of atoms in the solute, pi-carbons
  !cmbh - or not...  not really sure wtf this piece of code does, but
  !cmbh - i don't particularly like it, so i commented it out (Feb. 2008)
  
  !c     do k = iatomsa,iatomsa-10,-1
  !c        if ((ifra(k).lt.0).or.(ifra(k).gt.IMPA)) then
  !c           write(*,*) 'PROBLEM IN SOLUTE '
  !c           write(*,*) 'found one a problem',k,ifra(k),fra(k)
  !c        else
  !cmbh        write(32,1111) t,fra(k),ifra(k),fexxa(ifra(k),j,i),
  !cmbh x           fexya(ifra(k),j,i),fexza(ifra(k),j,i)
  !c        endif
  !c     enddo
      
  !1111 format(2e16.4,i5,10e16.4)
      
END SUBROUTINE fdotvout
