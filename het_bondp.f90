!c----------------------------------------------------------------------
!c     this subroutine computes the bond contributions to e , d and dd
!c----------------------------------------------------------------------
SUBROUTINE het_bondp(eng)
  !      subroutine het_bondp(v_d,v_x,e,v_b,wb_ppv,v_nc,v_iac,v_dsx)
  
  USE commondata, ONLY: npi,idsx,nexstate
  USE classical
  USE qcfio, ONLY: gen_io
  
  IMPLICIT NONE
  
  !c:::  input vars
  DOUBLE PRECISION, INTENT(out):: eng
  
  !c:::  local vars
  INTEGER :: i,j,k,l,ibond,btype,hit,hit1,hit2,typ1,typ2,idelt,ii
  INTEGER :: iatom,jatom,ibeg,i4,i3,i2,i1
  INTEGER, DIMENSION(2) :: itemp
  INTEGER, DIMENSION(2) :: ind,iind
  
  DOUBLE PRECISION :: rsq,r,bpms1,bpms2,f,df,ddf,ft,dft,ddft,sbond
  !DOUBLE PRECISION, DIMENSION(nbonds) :: bondl
  DOUBLE PRECISION, DIMENSION(3) :: dx
  DOUBLE PRECISION, DIMENSION(4) :: ftemp,atemp
  DOUBLE PRECISION, DIMENSION(3,2) :: dr
  DOUBLE PRECISION, DIMENSION(3,3,3) :: ddr
  
  do ibond=1,nbonds
     ind(1)=bonds(1,ibond)
     ind(2)=bonds(2,ibond)
     
     btype=pair_code(atype(ind(1)),atype(ind(2)))
     
     do k=1,3
        dx(k)=pos(k,ind(1))-pos(k,ind(2))
     end do
     
     rsq=length(dx)
     !rsq=dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)
     r=dsqrt(rsq)
     !bondl(ibond)=r
     
     call difbond(dx,dr,ddr,r)
     
     if(idsx.eq.0) then
        do j=1,3
           deriv_bond(1,j,ibond)=dr(j,1)
           deriv_bond(2,j,ibond)=-dr(j,1)
           !dsx(j,ibond)=dr(j)
           !dsx(3+j,ibond)=-dr(j)
        enddo
     else    !Pack ist derivative vector d and total energy e
        !call het_ebond(ic,ind1,ind2,r,v_nc,v_iac)
        !subroutine het_ebond(ic,i1,j1,bi,nc,iac,nap)
        hit1=0
        hit2=0
        
        do i=1,npi
           if(ind(1).eq.pi_iden(i)) hit1=i
           if(ind(2).eq.pi_iden(i)) hit2=i
        end do
        
        if(hit1.eq.0 .or. hit2.eq.0) then
           !quadratic bond energy function (f=cb*(b-b0)**2)
           typ1=atype(ind(1))
           typ2=atype(ind(2))
           !APW why?
           !if(typ2.eq.6 .and. (typ1.eq.5 .or. typ1.eq.3)) typ2=4
           if(species(typ1).eq.'B'.and.(species(typ2).eq.'A' .or. species(typ2).eq.'N')) typ1=locate_atom('C')
           hit=0
           i=0
           do while(i.lt.nbond_type .and. hit.eq.0)
              i=i+1
              if(btype.eq.bond_code(i)) hit=i
           end do
           if(hit.eq.0) then
              write(gen_io,*) '***warning: bond pair',typ1,typ2,'is missing from bonding parameters***'
              stop
           end if
           
           ftemp(1) = bond_pms(1,hit)
           ftemp(2) = bond_pms(2,hit)
           ftemp(3) = bond_pms(3,hit)
           
           f = ftemp(1)*(r-ftemp(2))**2-ftemp(3)
           df = 2*ftemp(1)*(r-ftemp(2))
           ddf = 2*ftemp(1)
        else
           !APW if bond is a pi-bond
           !Morse potential
           hit=0
           i=0
           do while(i.lt.nbond_type .and. hit.eq.0)
              i=i+1
              if(btype.eq.bond_code(i)) hit=i
           end do
           if(hit.ne.0) then
              ftemp(1) = bond_pms(1,hit)
              ftemp(2) = bond_pms(2,hit)
              ftemp(3) = bond_pms(3,hit)
              
              !APW trying this out
              do ii=1,nphib
                 iind(1)=phib_pairs(1,ii)
                 iind(2)=phib_pairs(2,ii)
                 if((ind(1).eq.iind(1).and.ind(2).eq.iind(2)).or.(ind(2).eq.iind(1).and.iind(1).eq.iind(2))) then
                    ftemp(2)=1.41
                 end if
              end do
              
              atemp(1)=dexp(-1.d0*ftemp(3)*(r-ftemp(2)))
              atemp(2)=atemp(1)**2
              f = ftemp(1)*(atemp(2)-2.d0*atemp(1))
              df = -ftemp(1)*ftemp(3)*(2.d0*atemp(2)-2d0*atemp(1))
              ddf = ftemp(1)*ftemp(3)**2*(4.d0*atemp(2)-2.d0*atemp(1))
           else
              write(gen_io,*) '***warning: bond pair',typ1,typ2,'is missing from bonding parameters***'
              stop
           end if
           !c:::  atype is used to check the special case
           !ccc iac=2 refers to an oxygen!!!
           itemp(1)=atype(ind(1))
           itemp(2)=atype(ind(2))
           if(species(itemp(1)).eq.'O' .or. species(itemp(2)).eq.'O') then
           !if(atype(ind(1)).eq.2 .or. atype(ind(2)).eq.2) then
              !c--------------------------------------------------------------
              !c
              !c     pi energy contributions for c=o bond .
              !c     the pi integrals and energy functions for c=o bond are calculated
              !c     here and not at vphi because there is no corresponding torsional
              !c     angle
              !c
              !c--------------------------------------------------------------
              hit=0
              i=0
              do while(i.lt.nbeta_type .and. hit.eq.0)
                 i=i+1
                 if(btype.eq.beta_code(i)) hit=i
              end do
              if(hit.eq.0) hit =1    !if not found use the default
              ftemp(1)=beta(1,hit)
              ftemp(2)=beta(2,hit)
              ftemp(3)=beta(3,hit)
              ftemp(4)=beta(4,hit)
              atemp(1) = ftemp(1)*dexp(-ftemp(3)*(r-ftemp(2))) !espo
              atemp(2) = r-ftemp(2) !dr
              atemp(3) = 1.d0+ftemp(4)*atemp(2) !trln
              ft = atemp(1)*atemp(3)
              dft = -ftemp(3)*ft+ftemp(4)*atemp(1)
              ddft = ftemp(3)*ftemp(3)*ft-2.d0*ftemp(3)*ftemp(4)*atemp(1)
              
              idelt=0
              if(nexstate.gt.0) then
                 idelt=1
              end if
              !sbond=con*(2.*p3(ihb,ihb1)+idelt*p1(ihb,ihb1))
              f=f+ft*sbond
              df=df+dft*sbond
              ddf=ddf+ddft*sbond
              !h(ihb,ihb1)=ft
              !h(ihb1,ihb)=ft
           end if
        end if
     end if
     eng=eng+f
     call drd(ind,df,ddf,dr,ddr)
     
  end do
  
END SUBROUTINE het_bondp
