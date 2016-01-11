!------------------------------------------------------------
!  This subroutine packs the arrays that contain the        !
!  different multibody interactions                         !
!  combines the subroutines linpak,thepak,phipak,nonbpk     !
!                                                           !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
!   Revised by Adam P. Willard 2010                         !
!   atom.willard@gmail.com                                  !
!                                                           !
!   NoCopyright                                             !
!                                                           !
!-----------------------------------------------------------!

SUBROUTINE allpak()
  
  USE commondata, ONLY:nnuc,log_delre
  USE classical
  USE qcfio, ONLY:gen_io
  
  IMPLICIT NONE

  !c:::  local vars
  CHARACTER*1 :: atom
  LOGICAL :: ltemp
  INTEGER :: i,j,k,icon,it,jt,kt,lt,l,m,n,i5,nb14,iat,jat,kat,lat,&
       exclude14,ii,lm2,lm1,icp2,iout,ntbond,icount,nrange,nrang1,&
       nrang2,lbn
  INTEGER, DIMENSION(4) :: itemp
  INTEGER, DIMENSION(4) :: tbon
  !INTEGER, DIMENSION(400) :: ilink
  INTEGER, DIMENSION(3,4000) :: tthetas
  INTEGER, DIMENSION(4,6000) :: tphis
  INTEGER, DIMENSION(2,8000) :: tbonds
  INTEGER, DIMENSION(4000) :: opln_atom
  INTEGER, DIMENSION(nnuc) :: nnbi
  INTEGER, DIMENSION(30,nnuc) :: ifij
  INTEGER, DIMENSION(nnuc) :: ifj,idj
  
  !      call linpak(nbond,natom,ib,jb,wlin_ppv,ilink,iprint)
  !c======================================================================
  !      subroutine linpak(nbond,natom,ib,jb,wlin_ppv,ilink,iprint)
  !c======================================================================
  
  icon=0
  do i=1,nnuc
     do j=1,4
        tbon(j)=0
     end do
     itemp(1)=0
     do j=1,nbonds
        if(bonds(1,j).eq.i) then
           itemp(1)=itemp(1)+1
           tbon(itemp(1))=bonds(2,j)
        end if
        if(bonds(2,j).eq.i) then
           itemp(1)=itemp(1)+1
           tbon(itemp(1))=bonds(1,j)
        end if
     end do
     
     !APW delre
     !if(use_delre.eq.1) then
     !   ilink(i)=(itemp(1)-1)/4+1
     !end if
     
     if(itemp(1).gt.4) then
        write(gen_io,*) '***too many bonds for atom ',i,'***'
        stop
     end if
     icon=icon+1
     !APW delre
     !do j=1,4
     !   links(j,icon)=tbon(j)
     !end do
  end do
  
  ntbond=nbonds
  do i=1,nbonds
     tbonds(1,i)=bonds(1,i)
     tbonds(2,i)=bonds(2,i)
  end do
     
  !      call thepak(nbond,nbo1,ib,jb,wt_ppv,ntheta)
  !c======================================================================
  !      subroutine thepak(nbond,nbo1,ib,jb,wt_ppv,ntheta)
  !c======================================================================
  !c     select atoms forming theta angles
  
  nthetas=0
  do i=1,(nbonds-1)
     itemp(2)=i+1
     do j=itemp(2),nbonds
        ltemp=.false.
        if(bonds(1,i).eq.bonds(1,j)) then
           it=bonds(2,i)
           jt=bonds(1,i)
           kt=bonds(2,j)
           ltemp=.true.
        end if
        if(bonds(1,i).eq.bonds(2,j)) then
           it=bonds(2,i)
           jt=bonds(1,i)
           kt=bonds(1,j)
           ltemp=.true.
        end if
        if(bonds(2,i).eq.bonds(1,j)) then
           it=bonds(1,i)
           jt=bonds(2,i)
           kt=bonds(2,j)
           ltemp=.true.
        end if
        if(bonds(2,i).eq.bonds(2,j)) then
           it=bonds(1,i)
           jt=bonds(2,i)
           kt=bonds(1,j)
           ltemp=.true.
        end if
        
        if(ltemp) then
           nthetas = nthetas+1
           ntbond=ntbond+1
           !APW not sure why the next two lines are necessary
           tbonds(1,ntbond)=it
           tbonds(2,ntbond)=kt
           tthetas(1,nthetas)=it
           tthetas(2,nthetas)=jt
           tthetas(3,nthetas)=kt
        end if
     end do
  end do
  
  ALLOCATE(thetas(3,nthetas))
  do i=1,nthetas
     thetas(1,i)=tthetas(1,i)
     thetas(2,i)=tthetas(2,i)
     thetas(3,i)=tthetas(3,i)
  end do
  
  !      call phipak(ib,jb,wp_ppv,iac,wt_ppv,ntheta,natom,nbond,nbo1,nphi,nof)
  !c---- select atom for phi angles ---------------------------
  !c======================================================================
  !      subroutine phipak(ib_,jb_,wp_ppv,iac_,wt_ppv,nth,
  !     $     natom_,nbond_,nbo1,nphi_,nof_)
  !c======================================================================
  exclude14=1
  nphis=0
  do i=1,nbonds
     !igen=0
     lm1=1
     it=atype(bonds(1,i))
     jt=atype(bonds(2,i))
     i5=pair_code(it,jt)
     do j=1,nbonds
        ltemp=.false.
        if(i.ne.j) then
           if(bonds(1,i).eq.bonds(2,j)) then
              it=bonds(1,j)
              jt=bonds(1,i)
              kt=bonds(2,i)
              ltemp=.true.
           end if
           if(bonds(1,i).eq.bonds(1,j)) then
              it=bonds(2,j)
              jt=bonds(1,j)
              kt=bonds(2,i)
              ltemp=.true.
           end if
           if(bonds(2,i).eq.bonds(2,j)) then
              it=bonds(1,j)
              jt=bonds(2,i)
              kt=bonds(1,i)
              ltemp=.true.
           end if
           if(bonds(2,i).eq.bonds(1,j)) then
              it=bonds(2,j)
              jt=bonds(2,i)
              kt=bonds(1,i)
              ltemp=.true.
           end if
           
           if(ltemp) then
              do k=j,nbonds
                 if(k.ne.i) then
                    lt=0
                    if(bonds(1,k).eq.kt) lt=bonds(2,k)
                    if(bonds(2,k).eq.kt) lt=bonds(1,k)
                    if(lt.ne.0) then
                       itemp(1)=atype(it)
                       itemp(2)=atype(jt)
                       itemp(3)=atype(kt)
                       itemp(4)=atype(lt)
                       
                       l=pair_code(itemp(2),itemp(3))
                       !c     prevent definition of torsional angle with b instead of c
                       if(species(itemp(1)).eq.'B') itemp(1)=locate_atom('H')
                       if(species(itemp(4)).eq.'B') itemp(4)=locate_atom('H')
                       if(species(itemp(1)).eq.'D') itemp(1)=locate_atom('H')
                       if(species(itemp(4)).eq.'D') itemp(4)=locate_atom('H')
                       !if(itemp(1).eq.6) itemp(1)=1
                       !if(itemp(4).eq.6) itemp(4)=1
                       !if(itemp(1).eq.7) itemp(1)=1
                       !if(itemp(4).eq.7) itemp(4)=1
                       lm2=itemp(1)+itemp(4)
                       !c     for 1-4 exclusion go to 110
                       if(exclude14.eq.0) then
                          ntbond=ntbond+1
                          tbonds(1,ntbond) = it
                          tbonds(2,ntbond) = lt
                       end if
                       ltemp=.true.
                       do ii=1,nphi_in
                          if(i5.eq.phi_in(ii)) then
                             ltemp=.false.
                             nphis=nphis+1
                             tphis(1,nphis)=it
                             tphis(2,nphis)=jt
                             tphis(3,nphis)=kt
                             tphis(4,nphis)=lt
                          end if
                       end do
                       !APW this code shouldn't get used unless possible phi-pairs are not defined in the input file.
                       !  if(lm2.gt.lm1 .and. ltemp) then
                       !     lm1=lm2
                       !     nphi_ppv=nphi_ppv+1
                       !     phis(1,nphi_ppv)=it
                       !     phis(2,nphi_ppv)=jt
                       !     phis(3,nphi_ppv)=kt
                       !     phis(4,nphi_ppv)=lt
                       !  end if
                    end if
                 end if
              end do
           end if
        end if
     end do
  end do
  
  if(nphis.eq.0) then 
     write(gen_io,*) "*** no phi angles ***"
     stop  !APW stop might not be necessary
  end if
  
  nphis_in=nphis
  !c
  !c     out of plane angles
  !c
  
  icp2=0
  iout=0
  do i=1,nthetas
     it=thetas(1,i)
     jt=thetas(2,i)
     kt=thetas(3,i)
     if(species(atype(jt)).eq.'N'.or. species(atype(jt)).eq.'A') then
        !APW generalizing atypes
        !if(atype(jt).eq.3 .or. atype(jt).eq.5) then
        ltemp=.true.
        if(iout.gt.0) then
           do j=1,iout
              if(jt.eq.opln_atom(j)) ltemp=.false.
           end do
        end if
        if(ltemp) then
           itemp(1)=atype(it)
           itemp(2)=atype(kt)
           if(species(itemp(1)).eq.'D') itemp(1)=locate_atom('H')
           !APW generalizing atypes
           !if(itemp(1).eq.7) itemp(1)=1
           do j=1,nbonds
              if(bonds(1,j).eq.jt .and. bonds(2,j).ne.it .and. bonds(2,j).ne.kt) then
                 lt=bonds(2,j)
                 icp2=icp2+1
                 nphis=nphis+1
                 tphis(1,nphis)=it
                 tphis(2,nphis)=jt
                 tphis(3,nphis)=kt
                 tphis(4,nphis)=lt
                 iout=iout+1
                 opln_atom(iout)=jt
              end if
              if(bonds(2,j).eq.jt .and. bonds(1,j).ne.it .and. bonds(1,j).ne.kt) then
                 lt=bonds(1,j)
                 icp2=icp2+1
                 nphis=nphis+1
                 tphis(1,nphis)=it
                 tphis(2,nphis)=jt
                 tphis(3,nphis)=kt
                 tphis(4,nphis)=lt
                 iout=iout+1
                 opln_atom(iout)=jt
              end if
           end do
        end if
     end if
  end do
  
  ALLOCATE(phis(4,nphis))
  
  do i=1,nphis
     phis(1,i)=tphis(1,i)
     phis(2,i)=tphis(2,i)
     phis(3,i)=tphis(3,i)
     phis(4,i)=tphis(4,i)
  end do

  ALLOCATE(deriv_bond(2,3,nbonds),deriv_theta(3,3,nthetas),deriv_phi(4,3,nphis))
  
  
  !APW totally modify this section, the nonbonded species will interact through nnucXnnuc matrix that is 1 if they interact or 0 if they don't
  
  do i=1,nnuc
     do j=1,nnuc
        nonbs(j,i)=1
     end do
  end do
  
  do i=1,ntbond
     nonbs(tbonds(1,i),tbonds(2,i))=0
     nonbs(tbonds(2,i),tbonds(1,i))=0
  end do
  

END SUBROUTINE allpak
